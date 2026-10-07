#!/usr/bin/env python3
"""Prepare the Allen ageing atlas as one portable, compressed expression file.

Usage:
    python prepare_atlas.py
    python prepare_atlas.py --verify-only

Requires Python, numpy, scipy, pandas, h5py, tqdm. CPU only; no GPU required.
Downloads the pinned 20241130 raw-count H5AD (13.85 GB) and metadata (0.42 GB),
verifies publisher MD5 checksums, and processes the complete released atlas.

Frozen representation:
  - library total over all released genes;
  - natural log1p(10000 * counts / library_total);
  - positive counts in >=ceil(0.001 * atlas_cells) cells across >=3 mice;
  - positive pooled log-expression variance; NO per-gene SD scaling;
  - one global divisor: pooled median retained-gene L2 cell norm;
  - float64 transformations/statistics and float32 stored expression.

Output senescence_atlas_processed.h5 is already losslessly compressed and
contains CSR row data, identifiers, original library totals, full source cell
metadata, preprocessing parameters and content fingerprints. It can be read
in bounded slices by h5py without unpacking another complete atlas.

Interrupted downloads and preprocessing resume with the same command. Keep
_atlas_downloads and _atlas_work until successful completion. This script
never silently overwrites an incompatible processed atlas.
"""

PROCESSOR_SOURCE_SHA256 = '4b860685d2288d8b4ae3de1564a5eebfd2773e070fcbed84a8e29ebc11d670c6'

"""Bounded, resumable downloads with system-CA TLS and publisher checksums."""
import concurrent.futures
import hashlib
import http.client
import json
import math
import os
from pathlib import Path
import re
import shutil
import ssl
import threading
import time
import urllib.request


def _download_ssl_context():
    """Use OS-trusted CAs; retain certificate-chain and hostname verification.

    Python 3.13's extra strict-X.509 checks reject some otherwise trusted Windows
    enterprise roots. Disabling only that flag restores earlier Python behavior;
    CERT_REQUIRED and check_hostname remain enabled.
    """
    context = ssl.create_default_context()
    if os.name == "nt" and hasattr(ssl, "VERIFY_X509_STRICT"):
        context.verify_flags &= ~ssl.VERIFY_X509_STRICT
    assert context.verify_mode == ssl.CERT_REQUIRED and context.check_hostname
    return context


def _download_md5(path, bar=None):
    if bar is not None:
        bar.unit, bar.unit_scale, bar.unit_divisor = "B", True, 1024
        bar.reset(total=path.stat().st_size)
        bar.set_description("Verify MD5 " + path.name)
    digest = hashlib.md5()
    with path.open("rb") as handle:
        while True:
            block = handle.read(8 * 1024**2)
            if not block:
                break
            digest.update(block)
            if bar is not None:
                bar.update(len(block))
    return digest.hexdigest()


def _download_allocated_bytes(path):
    """Physical allocation, not logical size (a .part may be a sparse file)."""
    if not path.exists():
        return 0
    stat = path.stat()
    if hasattr(stat, "st_blocks"):
        return int(stat.st_blocks) * 512
    if os.name == "nt":
        import ctypes
        from ctypes import wintypes
        kernel = ctypes.WinDLL("kernel32", use_last_error=True)
        query = kernel.GetCompressedFileSizeW
        query.argtypes = [wintypes.LPCWSTR, ctypes.POINTER(wintypes.DWORD)]
        query.restype = wintypes.DWORD
        high = wintypes.DWORD()
        ctypes.set_last_error(0)
        low = query(str(path.resolve()), ctypes.byref(high))
        if low != 0xFFFFFFFF or ctypes.get_last_error() == 0:
            return (high.value << 32) | low
    return 0  # Conservative if the platform cannot expose physical allocation.


def download_cached(url, path, bar=None, expected_bytes=None, workers=4,
                    chunk_bytes=32 * 1024**2, expected_md5=None,
                    max_range_seconds=180.0):
    """Return (Path, info); reset/update a supplied, reusable tqdm bar on the main thread.

    A completed file plus its .download.json receipt is reused without network access
    when its size, modification time and requested publisher checksum still match.
    Interrupted transfers resume .part chunks recorded in .part.json. Validation checks
    HTTP range, object identity where an ETag is available, and exact response lengths;
    an optional publisher MD5 is verified from disk before the final rename.
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    workers, chunk_bytes = int(workers), int(chunk_bytes)
    if workers < 1 or chunk_bytes < 1:
        raise ValueError("workers and chunk_bytes must be positive")
    if not math.isfinite(max_range_seconds) or max_range_seconds <= 0:
        raise ValueError("max_range_seconds must be finite and positive")
    expected_bytes = None if expected_bytes is None else int(expected_bytes)
    if expected_md5 is not None:
        expected_md5 = str(expected_md5).lower()
        if not re.fullmatch(r"[0-9a-f]{32}", expected_md5):
            raise ValueError("expected_md5 must contain exactly 32 hexadecimal digits")
    started = time.perf_counter()
    receipt = path.with_name(path.name + ".download.json")
    part = path.with_name(path.name + ".part")
    state_path = path.with_name(path.name + ".part.json")

    def read_json(filename):
        try:
            value = json.loads(filename.read_text(encoding="utf-8"))
            return value if isinstance(value, dict) else {}
        except (OSError, ValueError):
            return {}

    def atomic_json(filename, value):
        temporary = filename.with_name(filename.name + ".tmp")
        with temporary.open("w", encoding="utf-8") as handle:
            json.dump(value, handle, sort_keys=True)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, filename)

    def prepare_bar(size, completed):
        if bar is not None:
            bar.unit, bar.unit_scale, bar.unit_divisor = "B", True, 1024
            bar.reset(total=size)
            bar.set_description("Download " + path.name)
            if completed:
                bar.update(completed)

    cached = read_json(receipt)
    if (path.is_file() and cached.get("complete") is True
            and cached.get("url") == url and cached.get("size") == path.stat().st_size
            and (expected_bytes is None or cached["size"] == expected_bytes)):
        # A timestamp change invalidates the checksum receipt, not the whole file.
        current_stat = path.stat()
        trusted = (cached.get("mtime_ns") == current_stat.st_mtime_ns
                   and (expected_md5 is None or cached.get("verified_md5") == expected_md5))
        if not trusted:
            if expected_md5 is None:
                raise RuntimeError(f"Cached file has changed and no checksum is supplied: {path}")
            observed = _download_md5(path, bar)
            if observed != expected_md5:
                raise RuntimeError(f"Publisher MD5 mismatch in existing cache: {path}")
            cached.update(verified_md5=observed, mtime_ns=current_stat.st_mtime_ns)
            atomic_json(receipt, cached)
        prepare_bar(cached["size"], cached["size"])
        return path, {**cached, "cache_hit": True, "downloaded_bytes": 0,
                      "seconds": time.perf_counter() - started}
    if path.exists():
        # Recover a verified completed file if interruption occurred between final
        # rename and receipt commit, or import a separately downloaded source.
        if (expected_md5 is not None and path.is_file()
                and (expected_bytes is None or path.stat().st_size == expected_bytes)
                and _download_md5(path, bar) == expected_md5):
            cached = dict(url=url, size=path.stat().st_size, complete=True,
                          verified_md5=expected_md5, mtime_ns=path.stat().st_mtime_ns)
            atomic_json(receipt, cached)
            return path, {**cached, "cache_hit": True, "downloaded_bytes": 0,
                          "seconds": time.perf_counter() - started}
        raise RuntimeError(f"Existing cache lacks a matching verified download receipt: {path}")

    tls_context = _download_ssl_context()

    def open_range(range_headers):
        request = urllib.request.Request(url, headers=range_headers)
        response = urllib.request.urlopen(request, context=tls_context, timeout=30)
        # Initial connection timeout is 30 seconds; permit slower large reads.
        raw = getattr(getattr(response, "fp", None), "raw", None)
        sock = getattr(raw, "_sock", None)
        if sock is not None:
            sock.settimeout(120)
        return response

    headers = {"Accept-Encoding": "identity", "Range": "bytes=0-0"}
    for attempt in range(3):
        try:
            with open_range(headers) as response:
                if response.status != 206:
                    raise RuntimeError(f"Server must support HTTP ranges (got {response.status})")
                match = re.fullmatch(r"bytes 0-0/(\d+)", response.headers.get("Content-Range", ""))
                if not match or response.headers.get("Content-Encoding", "identity") != "identity":
                    raise RuntimeError("Invalid range probe or encoded HTTP response")
                size, etag = int(match.group(1)), response.headers.get("ETag")
                modified = response.headers.get("Last-Modified")
                probe = response.read(2)
                if len(probe) != 1:
                    raise RuntimeError("Range probe did not contain exactly one byte")
            break
        except (RuntimeError, OSError, http.client.HTTPException):
            if attempt == 2:
                raise
            time.sleep(0.5 * 2**attempt)
    if expected_bytes is not None and size != expected_bytes:
        raise RuntimeError(f"Object size is {size:,}, expected {expected_bytes:,} bytes")
    count = math.ceil(size / chunk_bytes)
    identity = {"url": url, "size": size, "etag": etag,
                "last_modified": modified, "chunk_bytes": chunk_bytes}
    state = read_json(state_path)
    reusable = (part.is_file() and part.stat().st_size == size
                and all(state.get(key) == value for key, value in identity.items())
                and bool(etag or modified))
    raw_done = state.get("done", []) if reusable else []
    done = set(raw_done) if isinstance(raw_done, list) and all(isinstance(i, int) for i in raw_done) else set()
    if not all(isinstance(index, int) and 0 <= index < count for index in done):
        done = set()
    lengths = lambda index: min(chunk_bytes, size - index * chunk_bytes)
    completed = sum(lengths(index) for index in done)
    # Existing allocation is reusable, including on Windows where truncate() may
    # reserve the whole file. Sparse logical length alone is never counted.
    allocation_needed = max(0, size - _download_allocated_bytes(part))
    if shutil.disk_usage(path.parent).free < allocation_needed + 256 * 1024**2:
        raise RuntimeError(f"Insufficient disk space for {allocation_needed:,} additional bytes")
    if not reusable:
        with part.open("wb") as handle:
            handle.truncate(size)
    state = {**identity, "done": sorted(done)}
    atomic_json(state_path, state)
    prepare_bar(size, completed)
    stop = threading.Event()

    def fetch(index):
        begin = index * chunk_bytes
        end = min(size, begin + chunk_bytes) - 1
        range_headers = {"Accept-Encoding": "identity", "Range": f"bytes={begin}-{end}"}
        if etag:
            range_headers["If-Match"] = etag
        for attempt in range(3):
            if stop.is_set():
                raise RuntimeError("Download interrupted")
            try:
                range_started = time.monotonic()
                with open_range(range_headers) as response:
                    if response.status != 206:
                        raise RuntimeError(f"Expected HTTP 206, got {response.status}")
                    if response.headers.get("Content-Range") != f"bytes {begin}-{end}/{size}":
                        raise RuntimeError("Server returned a mismatched byte range")
                    if etag and response.headers.get("ETag") != etag:
                        raise RuntimeError("Object ETag changed during download")
                    if modified and response.headers.get("Last-Modified") != modified:
                        raise RuntimeError("Object modification time changed during download")
                    if response.headers.get("Content-Encoding", "identity") != "identity":
                        raise RuntimeError("Encoded range responses are unsupported")
                    if int(response.headers.get("Content-Length", end - begin + 1)) != end - begin + 1:
                        raise RuntimeError("Mismatched range Content-Length")
                    written = 0
                    with part.open("r+b", buffering=0) as handle:
                        handle.seek(begin)
                        while True:
                            block = response.read(1024**2)
                            if time.monotonic() - range_started > max_range_seconds:
                                raise TimeoutError("Range transfer exceeded its elapsed-time limit")
                            if not block:
                                break
                            if stop.is_set():
                                raise RuntimeError("Download interrupted")
                            if written + len(block) > end - begin + 1:
                                raise RuntimeError("Range response exceeds requested length")
                            if handle.write(block) != len(block):
                                raise OSError("Short file write")
                            written += len(block)
                        if written != end - begin + 1:
                            raise RuntimeError("Incomplete range response")
                        os.fsync(handle.fileno())
                    return index, written
            except (RuntimeError, OSError, ValueError, http.client.HTTPException):
                if attempt == 2 or stop.is_set():
                    raise
                time.sleep(0.5 * 2**attempt)

    downloaded = 0
    executor = concurrent.futures.ThreadPoolExecutor(max_workers=workers)
    try:
        futures = [executor.submit(fetch, index) for index in range(count) if index not in done]
        for future in concurrent.futures.as_completed(futures):
            index, written = future.result()
            done.add(index)
            state["done"] = sorted(done)
            atomic_json(state_path, state)
            downloaded += written
            if bar is not None:
                bar.update(written)
    except BaseException:
        stop.set()
        executor.shutdown(wait=True, cancel_futures=True)
        raise
    finally:
        executor.shutdown(wait=True)
    if len(done) != count or part.stat().st_size != size:
        raise RuntimeError("Final download validation failed")
    observed_md5 = _download_md5(part, bar) if expected_md5 is not None else None
    if expected_md5 is not None and observed_md5 != expected_md5:
        # Keep the suspect .part for inspection; invalidate its chunk bitmap so a
        # rerun fetches all bytes rather than repeatedly trusting a corrupt chunk.
        atomic_json(state_path, {**identity, "done": [], "checksum_mismatch": observed_md5})
        raise RuntimeError(f"Publisher MD5 mismatch: {observed_md5} != {expected_md5}")
    os.replace(part, path)
    completed_receipt = {**identity, "complete": True, "verified_md5": observed_md5,
                         "mtime_ns": path.stat().st_mtime_ns}
    atomic_json(receipt, completed_receipt)
    state_path.unlink(missing_ok=True)
    return path, {**completed_receipt, "cache_hit": False,
                  "downloaded_bytes": downloaded,
                  "seconds": time.perf_counter() - started}




"""Bounded, resumable preprocessing into a portable compressed CSR HDF5 file.

The input is the published, QC-filtered Allen raw-count H5AD. There is no
centering, gene-SD scaling, cell exclusion, dense atlas, or network access here.
"""

from pathlib import Path
import hashlib
import json
import math
import os
import time

import h5py
import numpy as np
import pandas as pd
from scipy import sparse


_PROCESS_FORMAT = "senescence-csr-v1"


def _process_json(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False)


def _process_atomic_json(path, value):
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(_process_json(value), encoding="utf-8")
    os.replace(temporary, path)


def _process_atomic_npz(path, **arrays):
    temporary = path.with_name(path.name + ".tmp")
    with temporary.open("wb") as handle:
        np.savez(handle, **arrays)
    os.replace(temporary, path)


def _process_strings(node):
    if isinstance(node, h5py.Group):
        codes = np.asarray(node["codes"], dtype=np.int64)
        if np.any(codes < 0):
            raise ValueError("Missing H5AD identifier")
        return _process_strings(node["categories"])[codes]
    values = node.asstr(encoding="utf-8")[:] if h5py.check_string_dtype(node.dtype) else node[:]
    return np.asarray(values, dtype=str)


def _process_text_hash(values):
    digest = hashlib.sha256()
    for value in values:
        digest.update(str(value).encode("utf-8") + b"\0")
    return digest.hexdigest()


def _process_array_hash(values):
    return hashlib.sha256(np.ascontiguousarray(values).view(np.uint8)).hexdigest()


def _process_md5(path):
    digest = hashlib.md5()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024**2), b""):
            digest.update(block)
    return digest.hexdigest()


def _process_progress(bar, phase, done, total):
    if bar is None:
        return
    if getattr(bar, "_atlas_phase", None) != phase:
        bar._atlas_phase = phase
        bar.unit, bar.unit_scale = "block", False
        bar.reset(total=total)
    bar.set_description(phase, refresh=False)
    bar.set_postfix_str(f"{done:,}/{total:,}", refresh=False)
    bar.update(done - bar.n)
    bar.refresh()


def _process_row_sums(values, pointers):
    answer = np.zeros(len(pointers) - 1, dtype=np.float64)
    occupied = np.diff(pointers) > 0
    if occupied.any():
        answer[occupied] = np.add.reduceat(values, pointers[:-1][occupied], dtype=np.float64)
    return answer


def _process_read_counts(source, pointers, start, stop, genes):
    lo, hi = int(pointers[start]), int(pointers[stop])
    if hi - lo >= 2**31:
        raise ValueError("A row block exceeds int32 capacity; reduce block_rows")
    values = np.asarray(source["data"][lo:hi], dtype=np.float64)
    columns = np.asarray(source["indices"][lo:hi], dtype=np.int64)
    if not np.isfinite(values).all() or np.any(values < 0):
        raise ValueError(f"Nonfinite/negative source counts in rows {start}:{stop}")
    if np.any(columns < 0) or np.any(columns >= genes):
        raise ValueError("Invalid source gene index")
    counts = sparse.csr_matrix((values, columns.astype(np.int32),
        (pointers[start:stop + 1] - lo).astype(np.int32)), shape=(stop - start, genes))
    counts.sum_duplicates()
    counts.eliminate_zeros()
    return counts


def _process_log_values(counts, libraries, target_sum):
    scales = np.divide(target_sum, libraries, out=np.zeros_like(libraries), where=libraries > 0)
    values = counts.data.astype(np.float64, copy=True)
    values *= np.repeat(scales, np.diff(counts.indptr))
    np.log1p(values, out=values)
    return values


def _process_numeric(parent, name, data=None, *, shape=None, maxshape=None, dtype=None):
    if data is not None:
        data = np.asarray(data)
        dtype, shape = data.dtype, data.shape
    if len(shape) != 1:
        raise ValueError("Portable arrays are one-dimensional")
    options = dict(compression="gzip", compression_opts=1, shuffle=True, fletcher32=True)
    chunk = 1048576 if maxshape and maxshape[0] is None else max(1, min(1048576, max(shape[0], 1)))
    return parent.create_dataset(name, data=data, shape=shape, maxshape=maxshape,
                                 dtype=dtype, chunks=(chunk,), **options)


def _process_text_dataset(parent, name, values):
    encoded = [str(value).encode("utf-8") for value in values]
    width = max((len(value) for value in encoded), default=1)
    return _process_numeric(parent, name, np.asarray(encoded, dtype=f"S{max(width, 1)}"))


def _process_metadata(path, cell_ids, cfg):
    cell_column = cfg.get("cell_column", "cell_label")
    animal_column = cfg.get("animal_column", "donor_label")
    frame = pd.read_csv(path, usecols=[cell_column, animal_column], dtype=str)
    if frame[cell_column].isna().any() or frame[cell_column].duplicated().any():
        raise ValueError("Metadata cell identifiers must be unique and nonmissing")
    aligned = frame.set_index(cell_column).reindex(cell_ids)
    if aligned[animal_column].isna().any() or aligned[animal_column].str.strip().eq("").any():
        raise ValueError("Each expression row must map to one mouse in the metadata")
    codes, labels = pd.factorize(aligned[animal_column], sort=True)
    return codes.astype(np.int32), np.asarray(labels, dtype=str), {
        "cell_column": cell_column, "animal_column": animal_column,
        "cells": len(cell_ids), "animals": len(labels),
        "metadata_rows": len(frame), "metadata_extra_rows": len(frame) - len(cell_ids),
    }


def process_atlas(raw_path, metadata_path, output_path, config=None, bar=None):
    """Create one directly readable, compressed HDF5 CSR atlas.

    Three passes over the compressed raw CSR input calculate gene support and
    statistics, then the exact median retained-gene log-expression norm, then X.
    Global arrays are small; expression values are only decoded a block at a time.
    Completed work blocks and output row commits are reused on rerun. Sources and
    configuration are checked by content identity, never absolute file paths.
    Published source_md5 values must have been verified by the calling downloader.
    """
    cfg = dict(cell_fraction=0.001, min_animals=3, block_rows=4096, target_sum=10000.0)
    cfg.update(config or {})
    fraction, target = float(cfg["cell_fraction"]), float(cfg["target_sum"])
    if not 0 < fraction <= 1 or not np.isfinite(target) or target <= 0:
        raise ValueError("cell_fraction must lie in (0,1] and target_sum must be positive")
    for key in ("min_animals", "block_rows"):
        if not isinstance(cfg[key], (int, np.integer)) or cfg[key] <= 0:
            raise ValueError(f"{key} must be a positive integer")
    raw_path, metadata_path, output_path = map(Path, (raw_path, metadata_path, output_path))
    output_path.parent.mkdir(parents=True, exist_ok=True)
    if output_path.exists():
        manifest = verify_processed(output_path, config=cfg, bar=bar)
        return manifest
    if not raw_path.is_file() or not metadata_path.is_file():
        raise FileNotFoundError("Raw H5AD and cell metadata must exist before processing")
    work = output_path.parent / "_atlas_work"
    work.mkdir(exist_ok=True)
    part_path = output_path.with_name(output_path.name + ".part")
    sources = {}
    for key, path in (("raw", raw_path), ("metadata", metadata_path)):
        supplied = dict(cfg.get("sources", {}).get(key, {}))
        supplied["url"] = cfg.get("source_urls", {}).get(key, supplied.get("url", ""))
        supplied["md5"] = cfg.get("source_md5", {}).get(key, supplied.get("md5")) or _process_md5(path)
        supplied["bytes"] = path.stat().st_size
        sources[key] = supplied
    with h5py.File(raw_path, "r") as raw:
        source = raw["X"]
        encoding = source.attrs.get("encoding-type", "")
        if isinstance(encoding, bytes):
            encoding = encoding.decode()
        if encoding != "csr_matrix":
            raise ValueError("Expected a CSR-encoded raw-count H5AD")
        n, p = map(int, source.attrs["shape"])
        if n < 2 or p < 1:
            raise ValueError("Atlas needs at least two rows and one gene")
        cell_ids = _process_strings(raw["obs"][cfg.get("h5ad_cell_key", "cell_label")])
        gene_ids = _process_strings(raw["var"][cfg.get("h5ad_gene_key", "gene_identifier")])
        if len(cell_ids) != n or len(gene_ids) != p:
            raise ValueError("Identifiers disagree with matrix dimensions")
        if len(np.unique(cell_ids)) != n or len(np.unique(gene_ids)) != p:
            raise ValueError("Cell and gene identifiers must be unique")
        donor_codes, donor_labels, metadata_summary = _process_metadata(metadata_path, cell_ids, cfg)
        min_cells = math.ceil(fraction * n)
        parameters = dict(cell_fraction=fraction, min_cell_fraction=fraction,
            min_cells=min_cells, min_animals=int(cfg["min_animals"]), target_sum=target,
            library_target=target, log="natural_log1p", gene_sd_scaling=False,
            centering=False, positive_variance=True,
            global_scale="median_cell_l2_after_gene_filter", output_dtype="float32")
        identity = dict(format=_PROCESS_FORMAT, sources=sources, parameters=parameters,
            raw_shape=[n, p], cell_order_sha256=_process_text_hash(cell_ids),
            raw_gene_order_sha256=_process_text_hash(gene_ids))
        fingerprint = hashlib.sha256(_process_json(identity).encode()).hexdigest()
        work_identity = dict(fingerprint=fingerprint, block_rows=int(cfg["block_rows"]))
        work_identity_path = work / "identity.json"
        if work_identity_path.exists():
            if json.loads(work_identity_path.read_text(encoding="utf-8")) != work_identity:
                raise ValueError("Existing _atlas_work belongs to other source/settings/block size")
        else:
            if any(work.iterdir()) or part_path.exists():
                raise ValueError("Unidentified partial work exists; use a fresh output directory")
            _process_atomic_json(work_identity_path, work_identity)
        pointers = np.asarray(source["indptr"], dtype=np.int64)
        if (len(pointers) != n + 1 or pointers[0] != 0 or np.any(np.diff(pointers) < 0)
            or pointers[-1] != len(source["data"]) or pointers[-1] != len(source["indices"])):
            raise ValueError("Source CSR pointers disagree with its data")
        blocks = [(start, min(n, start + cfg["block_rows"])) for start in range(0, n, cfg["block_rows"])]
        library_sizes = np.empty(n, dtype=np.float64)
        sums, squares = np.zeros(p), np.zeros(p)
        detected = np.zeros(p, dtype=np.int64)
        presence = np.zeros((len(donor_labels), p), dtype=bool)
        _process_progress(bar, "1/4 Gene support and normalization statistics", 0, len(blocks))
        for index, (start, stop) in enumerate(blocks):
            stat_path = work / f"stats_{index:05d}.npz"
            if not stat_path.exists():
                counts = _process_read_counts(source, pointers, start, stop, p)
                libraries = _process_row_sums(counts.data, counts.indptr)
                logged = _process_log_values(counts, libraries, target)
                local_presence = np.zeros(len(donor_labels) * p, dtype=bool)
                per_value_donor = np.repeat(donor_codes[start:stop], np.diff(counts.indptr)).astype(np.int64)
                local_presence[per_value_donor * p + counts.indices] = True
                _process_atomic_npz(stat_path, fingerprint=np.asarray(fingerprint),
                    libraries=libraries, sums=np.bincount(counts.indices, weights=logged, minlength=p),
                    squares=np.bincount(counts.indices, weights=logged * logged, minlength=p),
                    detected=np.bincount(counts.indices, minlength=p).astype(np.int64),
                    presence=np.packbits(local_presence))
                del counts, logged, local_presence, per_value_donor
            with np.load(stat_path, allow_pickle=False) as saved:
                if str(saved["fingerprint"]) != fingerprint:
                    raise ValueError("Statistics checkpoint identity differs")
                library_sizes[start:stop] = saved["libraries"]
                sums += saved["sums"]
                squares += saved["squares"]
                detected += saved["detected"]
                presence |= np.unpackbits(saved["presence"], count=presence.size).reshape(presence.shape).astype(bool)
            _process_progress(bar, "1/4 Gene support and normalization statistics", index + 1, len(blocks))
        variance = np.maximum(0.0, (squares - sums * sums / n) / (n - 1))
        animal_support = presence.sum(axis=0).astype(np.int32)
        retained = np.flatnonzero((detected >= min_cells) &
            (animal_support >= cfg["min_animals"]) & (variance > 0)).astype(np.int32)
        if not len(retained):
            raise ValueError("No genes pass the specified support and variance thresholds")
        row_norms = np.empty(n, dtype=np.float64)
        _process_progress(bar, "2/4 Global median cell norm", 0, len(blocks))
        for index, (start, stop) in enumerate(blocks):
            norm_path = work / f"norms_{index:05d}.npz"
            if not norm_path.exists():
                counts = _process_read_counts(source, pointers, start, stop, p)[:, retained]
                logged = _process_log_values(counts, library_sizes[start:stop], target)
                norms = np.sqrt(_process_row_sums(logged * logged, counts.indptr))
                _process_atomic_npz(norm_path, fingerprint=np.asarray(fingerprint), norms=norms)
                del counts, logged
            with np.load(norm_path, allow_pickle=False) as saved:
                if str(saved["fingerprint"]) != fingerprint:
                    raise ValueError("Norm checkpoint identity differs")
                row_norms[start:stop] = saved["norms"]
            _process_progress(bar, "2/4 Global median cell norm", index + 1, len(blocks))
        median_norm = float(np.median(row_norms))
        if not np.isfinite(median_norm) or median_norm <= 0:
            raise ValueError("Retained expression has no finite positive median cell norm")
        manifest = dict(format=_PROCESS_FORMAT, fingerprint=fingerprint, identity=identity,
            sources=sources, parameters=parameters, raw_shape=[n, p], shape=[n, len(retained)],
            processor_source_sha256=cfg.get("processor_source_sha256"),
            raw_nnz=int(pointers[-1]), median_cell_norm=median_norm,
            animals=len(donor_labels), metadata_summary=metadata_summary,
            processing_block_rows=int(cfg["block_rows"]),
            retained_gene_order_sha256=_process_text_hash(gene_ids[retained]),
            library_sizes_sha256=_process_array_hash(library_sizes),
            donor_codes_sha256=_process_array_hash(donor_codes),
            donor_labels_sha256=_process_text_hash(donor_labels),
            cells_with_zero_retained_expression=int(np.count_nonzero(row_norms == 0)),
            zero_library_cells=int(np.count_nonzero(library_sizes == 0)))
        if not part_path.exists():
            initialization_path = part_path.with_name(part_path.name + ".initializing")
            with h5py.File(initialization_path, "w") as out:
                out.attrs["format"] = _PROCESS_FORMAT
                out.attrs["complete"] = False
                out.attrs["fingerprint"] = fingerprint
                out.attrs["manifest_json"] = _process_json(manifest)
                out.attrs["committed_rows"] = 0
                out.attrs["committed_nnz"] = 0
                x = out.create_group("X")
                x.attrs["encoding-type"], x.attrs["shape"] = "csr_matrix", [n, len(retained)]
                for name, dtype in (("data", "<f4"), ("indices", "<i4")):
                    _process_numeric(x, name, shape=(0,), maxshape=(None,), dtype=dtype)
                _process_numeric(x, "indptr", np.zeros(n + 1, dtype=np.int64))
                obs = out.create_group("obs")
                _process_text_dataset(obs, "cell_ids", cell_ids)
                _process_numeric(obs, "donor_codes", donor_codes)
                _process_numeric(obs, "library_sizes", library_sizes)
                var = out.create_group("var")
                _process_text_dataset(var, "gene_ids", gene_ids[retained])
                for name, values in (("raw_gene_indices", retained), ("detected_cells", detected[retained]),
                    ("detected_animals", animal_support[retained]), ("log_variance", variance[retained])):
                    _process_numeric(var, name, values)
                for symbol_key in ("gene_symbol", "gene_name", "feature_name"):
                    if symbol_key in raw["var"]:
                        _process_text_dataset(var, "gene_symbols", _process_strings(raw["var"][symbol_key])[retained])
                        break
                _process_text_dataset(out.create_group("donors"), "labels", donor_labels)
                all_genes = out.create_group("gene_filter")
                _process_text_dataset(all_genes, "raw_gene_ids", gene_ids)
                _process_numeric(all_genes, "detected_cells", detected)
                _process_numeric(all_genes, "detected_animals", animal_support)
                _process_numeric(all_genes, "log_variance", variance)
                meta = out.create_group("source_metadata")
                csv = _process_numeric(meta, "cell_metadata_csv", shape=(metadata_path.stat().st_size,), dtype=np.uint8)
                with metadata_path.open("rb") as stream:
                    offset = 0
                    for chunk in iter(lambda: stream.read(1048576), b""):
                        csv[offset:offset + len(chunk)] = np.frombuffer(chunk, dtype=np.uint8)
                        offset += len(chunk)
                checks = out.create_group("verification")
                _process_numeric(checks, "block_start", np.asarray([a for a, b in blocks], dtype=np.int64))
                _process_numeric(checks, "block_stop", np.asarray([b for a, b in blocks], dtype=np.int64))
                for name in ("data", "indices", "indptr"):
                    _process_numeric(checks, name + "_sha256", np.zeros(len(blocks), dtype="S64"))
                out.flush()
                if cfg.get("_interrupt_during_initialization"):
                    raise InterruptedError("Requested test interruption before atomic initialization")
            os.replace(initialization_path, part_path)
        _process_progress(bar, "3/4 Writing compressed float32 atlas", 0, len(blocks))
        with h5py.File(part_path, "r+") as out:
            if out.attrs["fingerprint"] != fingerprint:
                raise ValueError("Partial output fingerprint differs")
            committed_rows = int(out.attrs["committed_rows"])
            committed_nnz = int(out.attrs["committed_nnz"])
            if committed_rows not in {a for a, b in blocks} | {n}:
                raise ValueError("Partial output is not committed at a row-block boundary")
            if int(out["X/indptr"][committed_rows]) != committed_nnz:
                raise ValueError("Partial output commit pointer is inconsistent")
            for name in ("data", "indices"):
                if len(out["X"][name]) < committed_nnz:
                    raise ValueError("Partial output is shorter than its committed row data")
                out["X"][name].resize((committed_nnz,))
            for index, (start, stop) in enumerate(blocks):
                if stop <= committed_rows:
                    _process_progress(bar, "3/4 Writing compressed float32 atlas", index + 1, len(blocks))
                    continue
                counts = _process_read_counts(source, pointers, start, stop, p)[:, retained]
                logged = _process_log_values(counts, library_sizes[start:stop], target)
                values = (logged / median_norm).astype(np.float32)
                columns = counts.indices.astype(np.int32, copy=False)
                block_pointers = counts.indptr.astype(np.int64) + committed_nnz
                new_nnz = committed_nnz + len(values)
                if not np.isfinite(values).all() or np.any(values <= 0):
                    raise ValueError("Float32 conversion produced a nonpositive/nonfinite stored value")
                for name, array in (("data", values), ("indices", columns)):
                    dataset = out["X"][name]
                    dataset.resize((new_nnz,))
                    dataset[committed_nnz:new_nnz] = array
                    out["verification"][name + "_sha256"][index] = _process_array_hash(array).encode()
                out["X/indptr"][start:stop + 1] = block_pointers
                out["verification/indptr_sha256"][index] = _process_array_hash(block_pointers).encode()
                out.flush()
                out.attrs.modify("committed_rows", stop)
                out.attrs.modify("committed_nnz", new_nnz)
                out.flush()
                committed_rows, committed_nnz = stop, new_nnz
                del counts, logged, values, columns, block_pointers
                _process_progress(bar, "3/4 Writing compressed float32 atlas", index + 1, len(blocks))
                # This private hook supports a deterministic tiny interruption test.
                if cfg.get("_interrupt_after_write_blocks") == index + 1:
                    raise InterruptedError("Requested test interruption after a committed block")
            manifest["nnz"] = committed_nnz
            out.attrs.modify("manifest_json", _process_json(manifest))
            out.flush()
    verified = verify_processed(part_path, config=cfg, bar=bar, require_complete=False)
    manifest.update(array_sha256=verified["array_sha256"],
                    processed_median_cell_norm=verified["processed_median_cell_norm"],
                    verified_utc=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()), complete=True)
    with h5py.File(part_path, "r+") as out:
        out.attrs.modify("manifest_json", _process_json(manifest))
        out.attrs.modify("complete", True)
        out.flush()
    os.replace(part_path, output_path)
    return manifest


def verify_processed(path, config=None, bar=None, require_complete=True):
    """Read every compressed expression block, verify hashes and scale, no raw input.

    Per-block SHA256 was computed from the arrays before writing. Numeric-array
    hashes are also accumulated over the complete readback for later verification.
    The original metadata CSV is checked against its source MD5 while streaming.
    """
    with h5py.File(path, "r") as handle:
        if handle.attrs.get("format") != _PROCESS_FORMAT:
            raise ValueError("Unrecognized processed atlas format")
        if require_complete and not bool(handle.attrs.get("complete", False)):
            raise ValueError("Atlas is incomplete")
        manifest = json.loads(handle.attrs["manifest_json"])
        expected_fingerprint = hashlib.sha256(_process_json(manifest["identity"]).encode()).hexdigest()
        if manifest["fingerprint"] != expected_fingerprint or handle.attrs["fingerprint"] != expected_fingerprint:
            raise ValueError("Portable identity fingerprint differs")
        params = manifest["parameters"]
        if params != manifest["identity"]["parameters"] or params["gene_sd_scaling"] or params["centering"]:
            raise ValueError("Unexpected preprocessing parameters")
        for key in ("cell_fraction", "min_animals", "target_sum"):
            if config and key in config and params[key] != config[key]:
                raise ValueError(f"Prepared atlas has a different {key}")
        n, p = map(int, manifest["shape"])
        if params["min_cells"] != math.ceil(params["cell_fraction"] * n):
            raise ValueError("Relative cell support threshold is inconsistent")
        if config:
            for key in ("raw", "metadata"):
                for config_name, source_field in (("source_urls", "url"), ("source_md5", "md5")):
                    expected = config.get(config_name, {}).get(key)
                    if expected is not None and expected != manifest["sources"][key][source_field]:
                        raise ValueError(f"Prepared atlas has a different source {key} {source_field}")
                for source_field, expected in config.get("sources", {}).get(key, {}).items():
                    if expected != manifest["sources"][key].get(source_field):
                        raise ValueError(f"Prepared atlas has a different source {key} {source_field}")
        if list(handle["X"].attrs["shape"]) != [n, p]:
            raise ValueError("Matrix shape disagrees with manifest")
        for name, dtype in (("data", "float32"), ("indices", "int32"), ("indptr", "int64")):
            if handle["X"][name].dtype != np.dtype(dtype):
                raise ValueError(f"Unexpected CSR {name} dtype")
        pointers = np.asarray(handle["X/indptr"], dtype=np.int64)
        nnz = len(handle["X/data"])
        if (len(pointers) != n + 1 or pointers[0] != 0 or np.any(np.diff(pointers) < 0)
            or pointers[-1] != nnz or nnz != len(handle["X/indices"]) or nnz != manifest["nnz"]
            or int(handle.attrs["committed_rows"]) != n or int(handle.attrs["committed_nnz"]) != nnz):
            raise ValueError("Output CSR pointers/commit state disagree")
        cell_ids, gene_ids = _process_strings(handle["obs/cell_ids"]), _process_strings(handle["var/gene_ids"])
        donors = _process_strings(handle["donors/labels"])
        libraries = np.asarray(handle["obs/library_sizes"], dtype=np.float64)
        codes = np.asarray(handle["obs/donor_codes"], dtype=np.int32)
        if len(cell_ids) != n or len(gene_ids) != p or len(libraries) != n or len(codes) != n:
            raise ValueError("Identifiers/covariates have incorrect lengths")
        if (_process_text_hash(cell_ids) != manifest["identity"]["cell_order_sha256"]
            or _process_text_hash(gene_ids) != manifest["retained_gene_order_sha256"]
            or _process_text_hash(donors) != manifest["donor_labels_sha256"]
            or _process_array_hash(libraries) != manifest["library_sizes_sha256"]
            or _process_array_hash(codes) != manifest["donor_codes_sha256"]):
            raise ValueError("Identifiers/covariates failed their hashes")
        if not np.isfinite(libraries).all() or np.any(libraries < 0) or np.any(codes < 0) or np.any(codes >= len(donors)):
            raise ValueError("Invalid library size or donor code")
        if (np.any(handle["var/detected_cells"][:] < params["min_cells"])
            or np.any(handle["var/detected_animals"][:] < params["min_animals"])
            or np.any(handle["var/log_variance"][:] <= 0)):
            raise ValueError("Retained genes do not satisfy the frozen support/variance rule")
        raw_genes = _process_strings(handle["gene_filter/raw_gene_ids"])
        raw_indices = np.asarray(handle["var/raw_gene_indices"], dtype=np.int32)
        expected_retained = np.flatnonzero(
            (handle["gene_filter/detected_cells"][:] >= params["min_cells"])
            & (handle["gene_filter/detected_animals"][:] >= params["min_animals"])
            & (handle["gene_filter/log_variance"][:] > 0))
        if (len(raw_genes) != manifest["raw_shape"][1]
            or _process_text_hash(raw_genes) != manifest["identity"]["raw_gene_order_sha256"]
            or np.any(raw_indices < 0) or np.any(raw_indices >= len(raw_genes))
            or not np.array_equal(raw_indices, expected_retained)
            or not np.array_equal(raw_genes[raw_indices], gene_ids)):
            raise ValueError("Retained gene mapping is inconsistent")
        starts = np.asarray(handle["verification/block_start"], dtype=np.int64)
        stops = np.asarray(handle["verification/block_stop"], dtype=np.int64)
        if not len(starts) or starts[0] != 0 or stops[-1] != n or np.any(stops <= starts) or not np.array_equal(starts[1:], stops[:-1]):
            raise ValueError("Verification blocks do not partition the cell axis")
        hashes = {key: hashlib.sha256() for key in ("data", "indices", "indptr")}
        hashes["indptr"].update(pointers.view(np.uint8))
        norms = np.empty(n, dtype=np.float64)
        _process_progress(bar, "4/4 Verifying compressed atlas", 0, len(starts))
        for index, (start, stop) in enumerate(zip(starts, stops)):
            lo, hi = int(pointers[start]), int(pointers[stop])
            data = np.asarray(handle["X/data"][lo:hi])
            columns = np.asarray(handle["X/indices"][lo:hi])
            block_ptr = pointers[start:stop + 1]
            for name, array in (("data", data), ("indices", columns), ("indptr", block_ptr)):
                expected = bytes(handle["verification"][name + "_sha256"][index]).decode()
                if _process_array_hash(array) != expected:
                    raise ValueError(f"Block {index} {name} failed SHA256 readback")
                if name != "indptr":
                    hashes[name].update(array.view(np.uint8))
            if not np.isfinite(data).all() or np.any(data <= 0) or np.any(columns < 0) or np.any(columns >= p):
                raise ValueError("Invalid processed expression value or gene index")
            values64 = data.astype(np.float64)
            norms[start:stop] = np.sqrt(_process_row_sums(values64 * values64, block_ptr - lo))
            _process_progress(bar, "4/4 Verifying compressed atlas", index + 1, len(starts))
        array_hashes = {key: value.hexdigest() for key, value in hashes.items()}
        if "array_sha256" in manifest and array_hashes != manifest["array_sha256"]:
            raise ValueError("Complete numeric-array hash mismatch")
        median_norm = float(np.median(norms))
        if not np.isfinite(median_norm) or abs(median_norm - 1.0) > 2e-6:
            raise ValueError(f"Processed median cell norm is not one: {median_norm}")
        metadata_digest = hashlib.md5()
        csv = handle["source_metadata/cell_metadata_csv"]
        if len(csv) != manifest["sources"]["metadata"]["bytes"]:
            raise ValueError("Embedded metadata size differs")
        for begin in range(0, len(csv), 8 * 1024**2):
            metadata_digest.update(np.asarray(csv[begin:begin + 8 * 1024**2], dtype=np.uint8).tobytes())
        if metadata_digest.hexdigest() != manifest["sources"]["metadata"]["md5"]:
            raise ValueError("Embedded original metadata failed MD5 verification")
        return {**manifest, "array_sha256": array_hashes, "processed_median_cell_norm": median_norm}


def main(argv=None):
    """Download the published counts, preprocess once, and verify a portable atlas."""
    import argparse
    import datetime
    import json
    import os
    from pathlib import Path
    import sys
    import time
    from tqdm import tqdm

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir', type=Path, default=Path(__file__).resolve().parent)
    parser.add_argument('--workers', type=int, default=6, help='Concurrent bounded HTTP range transfers')
    parser.add_argument('--block-rows', type=int, default=4096)
    parser.add_argument('--verify-only', action='store_true', help='Fully verify the existing processed file')
    args = parser.parse_args(argv)
    destination = args.output_dir.resolve()
    destination.mkdir(parents=True, exist_ok=True)
    output = destination / 'senescence_atlas_processed.h5'
    downloads = destination / '_atlas_downloads'
    status_path = destination / 'processing_status.json'
    started = time.perf_counter()
    utc = lambda: datetime.datetime.now(datetime.timezone.utc).isoformat()
    status = dict(pid=os.getpid(), started_utc=utc(), state='RUNNING', phase='Starting',
                  output=str(output), cell_fraction=0.001, min_animals=3,
                  gene_sd_scaling=False, elapsed_seconds=0.0)
    last_status_write = 0.0

    def write_status(force=False, **fields):
        nonlocal last_status_write
        status.update(fields)
        now = time.perf_counter()
        if not force and now - last_status_write < 3:
            return
        status.update(updated_utc=utc(), elapsed_seconds=now - started)
        temporary = status_path.with_suffix('.json.tmp')
        temporary.write_text(json.dumps(status, indent=2, default=str) + '\n', encoding='utf-8')
        os.replace(temporary, status_path)
        last_status_write = now

    class Progress:
        def __init__(self):
            self._bar = tqdm(total=1, desc='Starting', unit='phase', mininterval=2,
                             dynamic_ncols=True, file=sys.stderr)
        @property
        def n(self):
            return self._bar.n
        @n.setter
        def n(self, value):
            self._bar.n = value
        @property
        def total(self):
            return self._bar.total
        @total.setter
        def total(self, value):
            self._bar.total = value
        @property
        def unit(self):
            return self._bar.unit
        @unit.setter
        def unit(self, value):
            self._bar.unit = value
        @property
        def unit_scale(self):
            return self._bar.unit_scale
        @unit_scale.setter
        def unit_scale(self, value):
            self._bar.unit_scale = value
        @property
        def unit_divisor(self):
            return self._bar.unit_divisor
        @unit_divisor.setter
        def unit_divisor(self, value):
            self._bar.unit_divisor = value
        def set_description(self, description=None, refresh=True, **kwargs):
            self._bar.set_description(description, refresh=refresh)
            write_status(phase=description, completed=self.n, total=self.total, unit=self.unit)
        def set_postfix_str(self, value='', refresh=True):
            self._bar.set_postfix_str(value, refresh=refresh)
            write_status(detail=value, completed=self.n, total=self.total, unit=self.unit)
        def reset(self, total=None):
            self._bar.reset(total=total)
            write_status(force=True, completed=self.n, total=self.total, unit=self.unit)
        def update(self, value=1):
            self._bar.update(value)
            write_status(completed=self.n, total=self.total, unit=self.unit)
        def refresh(self):
            self._bar.refresh()
            write_status(completed=self.n, total=self.total, unit=self.unit)
        def close(self):
            self._bar.close()

    write_status(force=True)
    bar = Progress()
    try:
        if args.verify_only:
            write_status(force=True, phase='Verifying existing processed atlas')
            report = verify_processed(output, bar=bar)
        else:
            downloads.mkdir(exist_ok=True)
            base = 'https://allen-brain-cell-atlas.s3.us-west-2.amazonaws.com/'
            raw_url = (base + 'expression_matrices/Zeng-Aging-Mouse-10Xv3/20241130/'
                       'Zeng-Aging-Mouse-10Xv3-raw.h5ad')
            metadata_url = base + 'metadata/Zeng-Aging-Mouse-10Xv3/20241130/cell_metadata.csv'
            raw_md5 = 'ba3fc54629ab87207135f4322efb7a13'
            metadata_md5 = '20944f6e4c54f13aca3c9f471f7fcb04'
            config = dict(cell_fraction=0.001, min_animals=3, block_rows=args.block_rows,
                          target_sum=10000.0, animal_column='donor_label',
                          processor_source_sha256=PROCESSOR_SOURCE_SHA256,
                          source_urls=dict(raw=raw_url, metadata=metadata_url),
                          source_md5=dict(raw=raw_md5, metadata=metadata_md5))
            raw, raw_info = download_cached(
                raw_url, downloads / 'counts.h5ad', bar=bar, expected_bytes=13849954183,
                expected_md5=raw_md5, workers=args.workers)
            metadata, metadata_info = download_cached(
                metadata_url, downloads / 'cell_metadata.csv', bar=bar,
                expected_bytes=423063567, expected_md5=metadata_md5, workers=args.workers)
            write_status(force=True, phase='Preprocessing atlas',
                         download_seconds=raw_info['seconds'] + metadata_info['seconds'])
            report = process_atlas(raw, metadata, output, config, bar=bar)
        write_status(force=True, state='COMPLETE', phase='Complete and verified', report=report)
        print(json.dumps(dict(output=str(output), bytes=output.stat().st_size,
                              elapsed_seconds=time.perf_counter() - started, report=report),
                         indent=2, default=str))
    except BaseException as error:
        write_status(force=True, state='STOPPED' if isinstance(error, KeyboardInterrupt) else 'FAILED',
                     error=f'{type(error).__name__}: {error}')
        raise
    finally:
        bar.close()


if __name__ == '__main__':
    main()
