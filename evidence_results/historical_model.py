"""Dictionary learning with observation weights specified as a function of margin.

Sparse coding minimises 0.5 * ||X - A S||_F^2 + alpha * sum(A), with A >= 0.
The rows of S are unit atoms. After solving for A, the dictionary update uses

    M = A + (X - A S) S.T - alpha,       W = weight_fn(M).

At a coding optimum, A = M.clamp_min(0). Numerically this equality holds to the
coding tolerance; weights are always defined on M, including its negative part.
The dictionary step projects -W.T @ (X - A S) / X.numel() onto the tangent
spaces at the atoms. Adam uses one second moment per atom, and every dictionary
step is followed by row normalisation.

Example (NumPy data, Torch weighting function)::

    model = MarginDictionaryLearning(
        400, weight_fn=lambda m: m.clamp(0, 0.5), random_state=1)
    codes = model.fit_transform(X)
    new_codes = model.transform(X_new)

Other weighting functions can be supplied directly::

    ordinary = lambda m: m.clamp_min(0)
    sigmoid = lambda m: torch.sigmoid(m / 0.4)
    active_sigmoid = lambda m: (m > 0) * torch.sigmoid(m / 0.4)

Smooth clipped weights equal Phi_T(m) - Phi_T(m - kappa). This equivalent
expression avoids cancellation for large positive margins::

    def smooth_clipped(m):
        T, kappa = 0.05, 0.5
        softplus = torch.nn.functional.softplus
        return m.clamp(0, kappa) + T * (
            softplus(-m.abs() / T) - softplus(-(m - kappa).abs() / T))

A quantile cap is also a function of the margin matrix::

    def quantile_weights(m):
        positive = m[m > 0]
        if positive.numel() == 0:
            return torch.zeros_like(m)
        return m.clamp_min(0).clamp_max(torch.quantile(positive, 0.8))

Such a cap is computed once per dictionary update and held fixed during it.
Matrix-dependent weights do not in general define one fixed separable gain.
Numerical work uses float32 arrays with float64 coding-certificate reductions
and TF32 disabled during fitting. The optional CPU Gram backend uses Torch
sparse products and Numba-compiled FISTA, without approximate decompositions.
CUDA replays the same FISTA calculations in blocks of ten steps.
"""

from contextlib import contextmanager, nullcontext
from functools import lru_cache
import math
import time

import numpy as np
import torch
from tqdm.auto import tqdm

__all__ = ["MarginDictionaryLearning"]


@contextmanager
def _full_precision():
    previous = torch.get_float32_matmul_precision()
    torch.set_float32_matmul_precision("highest")
    try:
        yield
    finally:
        torch.set_float32_matmul_precision(previous)


def _lipschitz(s):
    # Both Gram matrices have the same largest eigenvalue. For a symmetric B,
    # lambda_max(B)^2 <= ||B^2||_infinity. Use the smaller Gram matrix.
    gram = s.T @ s if s.shape[1] <= s.shape[0] else s @ s.T
    return 1.01 * (gram @ gram).abs().sum(1).max().sqrt().clamp_min(1e-12)


def _certificate(x, s, a, alpha, row_scale):
    r = x - a @ s
    correlation = r @ s.T
    mismatch = (a - (a + correlation - alpha).clamp_min(0)).abs()
    relative_error = (mismatch.amax(1) / a.abs().amax(1).clamp_min(1)).max()

    # U = q R is dual feasible because U S.T <= alpha (no absolute value).
    # This gap avoids subtracting two nearly equal objective values.
    ad, cd, rd = a.double(), correlation.double(), r.double()
    q = alpha / cd.amax(1).clamp_min(alpha)
    gap = (0.5 * (1 - q).square() * rd.square().sum(1)
           + (ad * (alpha - q[:, None] * cd)).sum(1))
    values = torch.stack((mismatch.max().double(), relative_error.double(),
                          (gap.clamp_min(0) / row_scale).max(),
                          (gap / row_scale).min()))
    return values, r, correlation


class _FISTA:
    """Warm-started nonnegative lasso; CPU and CUDA share the same steps."""

    def __init__(self, x, s, initial, alpha, tol, max_iter):
        self.x, self.alpha, self.tol = x, alpha, tol
        self.max_iter, self.block = max_iter, min(10, max_iter)
        self.s, self.a, self.y = s.clone(), initial.clone(), initial.clone()
        self.lip = _lipschitz(self.s)
        self.row_scale = (0.5 * x.double().square().sum(1)).clamp_min(1e-12)
        momentum, schedule = 1.0, []
        for _ in range(max_iter):
            next_momentum = (1 + math.sqrt(1 + 4 * momentum * momentum)) / 2
            schedule.append((momentum - 1) / next_momentum)
            momentum = next_momentum
        self.schedule = torch.tensor(schedule, dtype=x.dtype, device=x.device)
        self.beta = self.schedule[:self.block].clone()
        self.graph = None
        if x.is_cuda:
            with torch.cuda.device(x.device):
                current = torch.cuda.current_stream(x.device)
                stream = torch.cuda.Stream(device=x.device)
                stream.wait_stream(current)
                with torch.cuda.stream(stream):
                    for _ in range(2):
                        self._chunk(self.block)
                current.wait_stream(stream)
                self._reset(s, initial)
                stream.wait_stream(current)
                self.graph = torch.cuda.CUDAGraph()
                with torch.cuda.graph(self.graph, stream=stream):
                    self.graph_result = self._chunk(self.block)
                current.wait_stream(stream)
                # Keep the supplied starting values intact after graph setup.
                self._reset(s, initial)

    def _reset(self, s, initial):
        self.s.copy_(s)
        self.a.copy_(initial)
        self.y.copy_(self.a)
        self.lip.copy_(_lipschitz(self.s))

    def _chunk(self, count):
        for j in range(count):
            a_next = torch.relu(self.y + (
                (self.x - self.y @ self.s) @ self.s.T - self.alpha) / self.lip)
            y_next = a_next + self.beta[j] * (a_next - self.a)
            self.a.copy_(a_next)
            self.y.copy_(y_next)
        return _certificate(self.x, self.s, self.a, self.alpha, self.row_scale)

    def solve(self, s, initial):
        context = torch.cuda.device(self.x.device) if self.x.is_cuda else nullcontext()
        with context:
            self._reset(s, initial)
            result = _certificate(self.x, self.s, self.a, self.alpha, self.row_scale)
            steps = 0
            while True:
                absolute, relative, gap, minimum_gap = result[0].cpu().tolist()
                metrics = dict(inner_steps=steps, margin_error=absolute,
                               relative_margin_error=relative, relative_gap=gap,
                               minimum_relative_gap=minimum_gap)
                finite = all(math.isfinite(v) for v in metrics.values())
                if (finite and relative <= self.tol and gap <= self.tol
                        and minimum_gap >= -self.tol):
                    return self.a, result[1], result[2], metrics
                if not finite or steps == self.max_iter:
                    raise RuntimeError(
                        f"Sparse coding did not converge to {self.tol:g}: {metrics}")
                count = min(self.block, self.max_iter - steps)
                self.beta[:count].copy_(self.schedule[steps:steps + count])
                if self.graph is not None and count == self.block:
                    self.graph.replay()
                    result = self.graph_result
                else:
                    result = self._chunk(count)
                steps += count


@lru_cache(maxsize=1)
def _compiled_gram_solver():
    """Compile the same FISTA steps; no fast-math or reduced precision."""
    from numba import njit

    @njit(nogil=True)
    def certificate(a, corr, cross, norm2, alpha):
        # Same coefficient--margin error and feasible-dual gap as _certificate.
        penalty = np.float32(alpha)
        absolute = 0.0
        relative = 0.0
        max_gap = 0.0
        min_gap = np.inf
        resnorm = 0.0
        for i in range(len(a)):
            maxa = 1.0
            maxc = np.float64(alpha)
            mismatch = 0.0
            rr = norm2[i]
            for j in range(a.shape[1]):
                aa = a[i, j]
                cc = corr[i, j]
                if not (np.isfinite(aa) and np.isfinite(cc)):
                    return np.nan, np.nan, np.nan, np.nan, np.nan
                margin = np.float32(aa + cc) - penalty
                mm = abs(aa - max(np.float32(0), margin))
                mismatch = max(mismatch, np.float64(mm))
                maxa = max(maxa, np.float64(abs(aa)))
                maxc = max(maxc, np.float64(cc))
                rr -= np.float64(aa) * (np.float64(cross[i, j]) + np.float64(cc))
            rr = max(rr, 0.0)
            resnorm += rr
            q = np.float64(alpha) / maxc
            gap = 0.5 * (1 - q) ** 2 * rr
            for j in range(a.shape[1]):
                gap += np.float64(a[i, j]) * (np.float64(alpha) - q * np.float64(corr[i, j]))
            gap /= max(0.5 * norm2[i], 1e-12)
            absolute = max(absolute, mismatch)
            relative = max(relative, mismatch / maxa)
            max_gap = max(max_gap, gap)
            min_gap = min(min_gap, gap)
        return (absolute, relative, max_gap, min_gap, resnorm)

    @njit(nogil=True)
    def solve(a, cross, gram, norm2, alpha, lip, schedule, tol, maxiter):
        a = a.copy()
        y = a.copy()
        steps = 0
        penalty = np.float32(alpha)
        while True:
            corr = cross - a @ gram
            vals = certificate(a, corr, cross, norm2, alpha)
            absolute, relative, gap, min_gap, resnorm = vals
            finite = np.isfinite(absolute + relative + gap + min_gap)
            if finite and relative <= tol and (gap <= tol) and (min_gap >= -tol):
                return (a, corr, vals, steps)
            if not finite or steps == maxiter:
                raise RuntimeError('Gram sparse coding did not converge')
            count = min(10, maxiter - steps)
            for j in range(count):
                yg = y @ gram
                beta = schedule[steps + j]
                for i in range(len(a)):
                    for k in range(a.shape[1]):
                        gradient = np.float32(cross[i, k] - yg[i, k]) - penalty
                        anew = max(np.float32(0), y[i, k] + gradient / lip)
                        y[i, k] = anew + beta * (anew - a[i, k])
                        a[i, k] = anew
            steps += count
    return solve


class _GramFISTA:
    """Exact residual algebra for a small dictionary and sparse CPU data.

    C = X S.T, G = S S.T. Thus R S.T = C - A G and
    ||R_i||^2 = ||X_i||^2 - A_i dot (C_i + (R S.T)_i).
    No entries, components, or singular values are discarded.
    """

    def __init__(self, x, s, initial, alpha, tol, max_iter):
        from scipy.sparse import csr_matrix
        import warnings

        self.x, self.alpha, self.tol, self.max_iter = x, alpha, tol, max_iter
        sparse = csr_matrix(x.numpy())
        def tensor(matrix):
            # PyTorch emits this API-status notice once, not a numerical warning.
            with warnings.catch_warnings():
                warnings.filterwarnings("ignore", message="Sparse CSR tensor support is in beta state")
                return torch.sparse_csr_tensor(
                    torch.from_numpy(matrix.indptr), torch.from_numpy(matrix.indices),
                    torch.from_numpy(matrix.data), size=matrix.shape)
        self.x_csr, self.xt_csr = tensor(sparse), tensor(sparse.T.tocsr())
        self.norm2 = x.double().square().sum(1).numpy()
        momentum, schedule = 1.0, []
        for _ in range(max_iter):
            next_momentum = (1 + math.sqrt(1 + 4 * momentum * momentum)) / 2
            schedule.append((momentum - 1) / next_momentum)
            momentum = next_momentum
        self.schedule = np.asarray(schedule, dtype=np.float32)
        self.compiled_solve = _compiled_gram_solver()

    def solve(self, s, initial):
        cross, gram = self.x_csr @ s.T, s @ s.T
        a, correlation, values, steps = self.compiled_solve(
            initial.numpy(), cross.numpy(), gram.numpy(), self.norm2,
            self.alpha, np.float32(_lipschitz(s).item()), self.schedule,
            self.tol, self.max_iter)
        absolute, relative, gap, minimum_gap, residual_norm2 = values
        metrics = dict(inner_steps=steps, margin_error=absolute,
                       relative_margin_error=relative, relative_gap=gap,
                       minimum_relative_gap=minimum_gap)
        return (torch.from_numpy(a), torch.tensor(residual_norm2, dtype=torch.float64),
                torch.from_numpy(correlation), metrics)

    def gradient(self, s, a, weights):
        # -W.T (X - A S) / (n p), without constructing an n-by-p residual.
        return ((weights.T @ a) @ s - (self.xt_csr @ weights).T) / self.x.numel()

    def audit(self, s, a, used_steps):
        # Independently check the original residual-based certificate. If rounding
        # changes a borderline decision, refine the final codes with that solver.
        remaining = self.max_iter - used_steps
        a, residual, correlation, metrics = _FISTA(
            self.x, s, a, self.alpha, self.tol, remaining).solve(s, a)
        return a, residual.double().square().sum(), correlation, metrics


def _dictionary_step(s, weights, residual, first, second, iteration, lr, gradient=None):
    if gradient is None:
        gradient = -(weights.T @ residual) / residual.numel()
    gradient -= s * (gradient * s).sum(1, keepdim=True)
    first.mul_(0.9).add_(gradient, alpha=0.1)
    first -= s * (first * s).sum(1, keepdim=True)
    second.mul_(0.999).add_(gradient.square().mean(1, keepdim=True), alpha=0.001)
    updated = s - lr * (first / (1 - 0.9 ** iteration)) / (
        (second / (1 - 0.999 ** iteration)).sqrt() + 1e-8)
    return updated / updated.norm(dim=1, keepdim=True).clamp_min(1e-12)


class MarginDictionaryLearning:
    """Unit-norm dictionary learning with a callable observation-weight rule.

    Parameters
    ----------
    n_components : int
        Number of atoms.
    weight_fn : callable
        Maps the entire (n_observations, n_components) margin tensor to finite,
        nonnegative weights with the same shape, dtype and device. It is called
        once per dictionary step, after sparse coding, without autograd.
    alpha : float, default=0.2
        Positive coefficient penalty in 0.5 * ||X - A S||^2 + alpha * sum(A).
    lr : float, default=0.01
        Dictionary Adam learning rate (betas 0.9, 0.999; epsilon 1e-8).
    n_iter : int, default=5000
        Dictionary updates, not a claim of dictionary convergence.
    code_tol, code_max_iter : default=1e-5, 5000
        Sparse-coding tolerance and maximum FISTA steps per solve. Both the
        relative coefficient--margin error and relative dual gap must pass.
    dict_init, code_init : array-like, optional
        Starting dictionary (K, p), normalised on input, and optional codes
        (n, K). Otherwise sample training rows and start codes at zero.
    random_state : int, optional
        Seed for sampling the initial atoms.
    device : str, default='auto'
        'auto' selects CUDA when available, otherwise CPU. An explicit device
        is used as requested; GPU errors are not silently redirected to CPU.
    solver : {"direct", "gram"}, default="direct"
        The default computes residuals explicitly on CPU or CUDA. The optional
        CPU Gram backend uses exact sparse products and compiled FISTA loops
        (SciPy and Numba required), with the same objective, updates and stopping
        thresholds. Floating-point operation order differs. Returned codes also
        pass the original residual-based certificate. Most useful when K << p.
    record_every : int, default=100
        Record iteration 0, every record_every updates, and the final state.
    verbose : bool, default=True
        Show one updating progress bar with elapsed time and an ETA.
    callback : callable, optional
        callback(iteration, dictionary) receives a NumPy copy at each recorded
        state. It may return a dict of additional metrics, e.g. a recovery
        score, with names distinct from the built-in history fields.

    Attributes
    ----------
    components_, coefficients_ : numpy.ndarray
        Final dictionary and converged training codes, both float32.
    history_ : list of dict
        Each row describes the dictionary AFTER 'iteration' updates, with its
        codes refitted. Includes objective per observation, coding diagnostics,
        mean inner steps, elapsed seconds, and maximum atom displacement in
        the most recent dictionary step (Euclidean distance between unit rows).
    final_coding_diagnostics_ : dict
        Certificate for the returned dictionary and training coefficients.
    n_iter_, training_time_, device_ : fitted metadata
    """

    def __init__(self, n_components, weight_fn, *, alpha=0.2, lr=0.01,
                 n_iter=5000, code_tol=1e-5, code_max_iter=5000,
                 dict_init=None, code_init=None, random_state=None,
                 device="auto", record_every=100, verbose=True, callback=None,
                 solver="direct"):
        self.n_components, self.weight_fn = n_components, weight_fn
        self.alpha, self.lr, self.n_iter = alpha, lr, n_iter
        self.code_tol, self.code_max_iter = code_tol, code_max_iter
        self.dict_init, self.code_init = dict_init, code_init
        self.random_state, self.device = random_state, device
        if solver not in ("direct", "gram"):
            raise ValueError("solver must be direct or gram")
        self.solver = solver
        self.record_every, self.verbose, self.callback = record_every, verbose, callback
        for name in ("n_components", "code_max_iter", "record_every", "n_iter"):
            value = getattr(self, name)
            if not isinstance(value, (int, np.integer)) or value < (0 if name == "n_iter" else 1):
                raise ValueError(f"{name} must be an integer >= {0 if name == 'n_iter' else 1}")
        if any(not math.isfinite(v) or v <= 0 for v in (alpha, lr, code_tol)):
            raise ValueError("alpha, lr and code_tol must be finite and positive")
        if not callable(weight_fn):
            raise TypeError("weight_fn must be callable")

    def _device(self):
        name = "cuda" if torch.cuda.is_available() else "cpu"
        return torch.device(name if self.device == "auto" else self.device)

    @staticmethod
    def _data(X, device):
        x = torch.as_tensor(X, dtype=torch.float32, device=device).detach().clone()
        if x.ndim != 2 or min(x.shape) == 0 or not torch.isfinite(x).all():
            raise ValueError("X must be a finite, nonempty two-dimensional array")
        return x

    @torch.no_grad()
    def fit(self, X):
        """Fit the dictionary and return self. Failed coding solves raise an error."""
        with _full_precision():
            device = self._device()
            x = self._data(X, device)
            if self.dict_init is None:
                rng = np.random.default_rng(self.random_state)
                indices = rng.choice(len(x), self.n_components, replace=self.n_components > len(x))
                initial = x[torch.as_tensor(indices, device=device)].double()
            else:
                initial = torch.as_tensor(self.dict_init, dtype=torch.float64, device=device)
            if (initial.shape != (self.n_components, x.shape[1])
                    or not torch.isfinite(initial).all() or (initial.norm(dim=1) == 0).any()):
                raise ValueError("Initial atoms must be finite, nonzero rows of shape (K, p)")
            s = (initial / initial.norm(dim=1, keepdim=True)).float()
            s /= s.norm(dim=1, keepdim=True)
            a = (x.new_zeros((len(x), self.n_components))
                 if self.code_init is None else torch.as_tensor(
                     self.code_init, dtype=torch.float32, device=device).clone())
            if a.shape != (len(x), self.n_components) or not (torch.isfinite(a) & (a >= 0)).all():
                raise ValueError("Initial codes must be finite, nonnegative and have shape (n, K)")

            start = time.perf_counter()
            if self.solver == "gram" and device.type != "cpu":
                raise ValueError("The Gram solver is a CPU backend; use direct for CUDA")
            solver_class = _GramFISTA if self.solver == "gram" else _FISTA
            solver = solver_class(x, s, a, self.alpha, self.code_tol, self.code_max_iter)
            first, second = torch.zeros_like(s), torch.zeros_like(s[:, :1])
            movement = s.new_zeros(())
            self.history_ = []
            inner_total, last_recorded = 0, 0
            with tqdm(total=self.n_iter, unit="update", desc="Dictionary learning",
                      disable=not self.verbose) as progress:
                for iteration in range(self.n_iter + 1):
                    a, residual, correlation, metrics = solver.solve(s, a)
                    if self.solver == "gram" and iteration == self.n_iter:
                        inner_steps = metrics["inner_steps"]
                        a, residual, correlation, metrics = solver.audit(s, a, inner_steps)
                        metrics["inner_steps"] += inner_steps
                    inner_total += metrics["inner_steps"]
                    if iteration % self.record_every == 0 or iteration == self.n_iter:
                        residual_norm2 = (residual if self.solver == "gram"
                                          else residual.double().square().sum())
                        objective = (0.5 * residual_norm2 + self.alpha * a.double().sum()) / len(x)
                        row = dict(iteration=iteration, objective=objective.item(),
                                   max_atom_change=movement.item(),
                                   mean_inner_steps=inner_total / (iteration + 1), **metrics)
                        extras = (self.callback(iteration, s.cpu().numpy().copy())
                                  if self.callback is not None else None)
                        if extras is not None:
                            if row.keys() & extras.keys() or "seconds" in extras:
                                raise ValueError("Callback metrics must have distinct names from history fields")
                            row.update(extras)
                        row["seconds"] = time.perf_counter() - start
                        self.history_.append(row)
                        display = dict(inner=metrics["inner_steps"],
                                       gap=f'{metrics["relative_gap"]:.1e}',
                                       move=f'{row["max_atom_change"]:.1e}')
                        display.update(extras or {})
                        progress.set_postfix(display, refresh=False)
                        progress.update(iteration - last_recorded)
                        last_recorded = iteration
                    if iteration == self.n_iter:
                        break
                    margins = a + correlation - self.alpha
                    weights = self.weight_fn(margins)
                    if (not isinstance(weights, torch.Tensor) or weights.shape != a.shape
                            or weights.device != a.device or weights.dtype != a.dtype):
                        raise ValueError("weight_fn must return a tensor with the margin shape, dtype and device")
                    if not (torch.isfinite(weights) & (weights >= 0)).all():
                        raise ValueError("weight_fn must return finite, nonnegative weights")
                    gradient = solver.gradient(s, a, weights) if self.solver == "gram" else None
                    updated = _dictionary_step(s, weights, residual, first, second,
                                               iteration + 1, self.lr, gradient=gradient)
                    movement = (updated - s).norm(dim=1).max()
                    s = updated

            self.components_ = s.cpu().numpy().copy()
            self.coefficients_ = a.cpu().numpy().copy()
            self.final_coding_diagnostics_ = metrics.copy()
            self.n_iter_, self.n_features_in_ = self.n_iter, x.shape[1]
            self.training_time_, self.device_ = time.perf_counter() - start, str(x.device)
        return self

    @torch.no_grad()
    def transform(self, X):
        """Solve nonnegative lasso for new data using the fitted dictionary."""
        if not hasattr(self, "components_"):
            raise RuntimeError("Call fit before transform")
        with _full_precision():
            device = self._device()
            x = self._data(X, device)
            if x.shape[1] != self.n_features_in_:
                raise ValueError("X must have the same number of features as the training data")
            s = torch.as_tensor(self.components_, device=device)
            initial = x.new_zeros((len(x), self.n_components))
            if self.solver == "gram" and device.type != "cpu":
                raise ValueError("The Gram solver is a CPU backend; use direct for CUDA")
            solver_class = _GramFISTA if self.solver == "gram" else _FISTA
            solver = solver_class(x, s, initial, self.alpha, self.code_tol, self.code_max_iter)
            a, _, _, metrics = solver.solve(s, initial)
            if self.solver == "gram":
                inner_steps = metrics["inner_steps"]
                a, _, _, metrics = solver.audit(s, a, inner_steps)
                metrics["inner_steps"] += inner_steps
            self.transform_diagnostics_ = metrics
            return a.cpu().numpy().copy()

    def fit_transform(self, X):
        """Return converged training codes for the final fitted dictionary."""
        return self.fit(X).coefficients_

    def inverse_transform(self, coefficients):
        """Reconstruct observations as coefficients @ components_."""
        return np.asarray(coefficients) @ self.components_
