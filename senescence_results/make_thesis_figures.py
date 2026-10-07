"""Make thesis figures from the frozen, selected-program numerical CSVs.

Run from anywhere: python senescence_results/make_thesis_figures.py
No enrichment or SenMayo files are read. Exact zeros are retained as reported
probability masses; distribution curves are not renormalised after removing them.
"""
from pathlib import Path
import hashlib
import json

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.lines import Line2D
from matplotlib.ticker import LogLocator, MaxNLocator, NullLocator
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parent
SOURCE = ROOT / 'analysis' / 'top_7_4fa212641968'
OUTPUT = ROOT / 'thesis_figures'
OUTPUT.mkdir(exist_ok=True)
COLORS = {'identity': '#C47A44', 'activity': '#2678AC'}
AGE_COLORS = {'young': '#2678AC', 'aged': '#C47A44'}
INPUTS = {}
plt.rcParams.update({
    'font.family': 'STIXGeneral', 'mathtext.fontset': 'stix', 'font.size': 10,
    'axes.labelsize': 10, 'xtick.labelsize': 8, 'ytick.labelsize': 8,
    'axes.linewidth': .6, 'axes.edgecolor': '#A4ABB1',
    'xtick.color': '#39434A', 'ytick.color': '#39434A',
    'text.color': '#263238', 'axes.labelcolor': '#263238',
    'axes.spines.top': False, 'axes.spines.right': False,
    'pdf.fonttype': 42, 'ps.fonttype': 42,
})


def read_csv(path, **kwargs):
    INPUTS[str(path.relative_to(ROOT))] = hashlib.sha256(path.read_bytes()).hexdigest()
    return pd.read_csv(path, float_precision='round_trip', **kwargs)


def style(ax, grid='y'):
    ax.set_axisbelow(True)
    ax.grid(axis=grid, color='#E8EBEE', linewidth=.6)
    ax.tick_params(length=3, width=.6, pad=3)


def save(fig, name):
    assert fig.get_figwidth() <= 6.3 + 1e-10
    assert fig._suptitle is None and all(not ax.get_title() for ax in fig.axes)
    fig.savefig(OUTPUT / f'{name}.pdf', bbox_inches='tight', pad_inches=.04)
    fig.savefig(OUTPUT / f'{name}.png', dpi=220, bbox_inches='tight', pad_inches=.04)
    plt.close(fig)


scores = read_csv(ROOT / 'analysis' / 'screen' / 'program_scores.csv').set_index('program_id')
top = read_csv(SOURCE / 'selected_programs.csv').set_index('program_id')
assert top.index.tolist() == [179, 163, 122, 196, 188, 60, 124]
assert top['rank'].tolist() == list(range(1, 8))
np.testing.assert_allclose(scores.loc[top.index, 'score'], top['score'], rtol=1e-12)

# Full screen and a magnified view of the decision boundary; no genes inform it.
fig, axes = plt.subplots(1, 2, figsize=(6.3, 2.85), gridspec_kw={'width_ratios': [1.2, 1]})
for ax, data in zip(axes, [scores, scores.loc[scores['rank'] <= 20]]):
    style(ax)
    for label, color in COLORS.items():
        part = data.loc[data.label == label]
        ax.scatter(part['rank'], part.score, s=13, color=color, linewidths=.2,
                   edgecolors='white', label=label.capitalize(), zorder=3)
    ax.axvline(len(top) + .5, color='#929BA3', lw=.8, ls=(0, (3, 3)))
axes[0].axhline(0, color='#929BA3', lw=.7)
axes[0].set(xlabel='Rank among 240 programs', ylabel='Normalised age-tail score', xlim=(0, 244))
axes[0].set_xticks([1, 60, 120, 180, 240])
axes[0].legend(frameon=False, fontsize=8, loc='upper right', handletextpad=.3)
zoom = scores.loc[scores['rank'] <= 20]
axes[1].set(xlabel='Rank among the first 20', xlim=(.3, 20.8),
            ylim=(max(0, zoom.score.min() - .035), top.score.max() + .065))
axes[1].set_xticks([1, 5, 7, 10, 15, 20])
for program, row in top.iterrows():
    axes[1].annotate(str(program), (row['rank'], row.score), xytext=(4, 3),
                     textcoords='offset points', fontsize=7, va='bottom')
fig.tight_layout(pad=.6, w_pad=.8)
save(fig, 'fig_senescence_ranking')

# Two aligned class matrices: vertical stacking preserves all 24 readable labels.
usage = read_csv(SOURCE / 'candidate_class_usage.csv', index_col='program_id').loc[top.index]
contribution = read_csv(SOURCE / 'candidate_class_contributions.csv', index_col='program_id').loc[top.index]
assert usage.columns.tolist() == contribution.columns.tolist()
np.testing.assert_allclose(usage.sum(axis=1), 1, atol=1e-12)
np.testing.assert_allclose(contribution.sum(axis=1), top.score, rtol=1e-10, atol=1e-12)
abbreviations = ['IT/ET', 'NP/CT/L6b', 'OB-CR', 'DG-IMN', 'OB-IMN', 'CTX-CGE',
                 'CTX-MGE', 'CNU-LGE', 'CNU-HYa GABA', 'HY GABA', 'CNU-HYa Glut',
                 'HY Glut', 'HY-MM', 'MB Glut', 'MB GABA', 'MB Dopa', 'MB/HB Sero',
                 'P Glut', 'MY Glut', 'P GABA', 'Astro/Epen', 'OPC/Oligo', 'Vascular', 'Immune']
assert len(abbreviations) == len(usage.columns) == 24
class_labels = [f'{c.split()[0]} {short}' for c, short in zip(usage.columns, abbreviations)]
pd.DataFrame({'class_name': usage.columns, 'figure_label': class_labels}).to_csv(
    OUTPUT / 'class_label_key.csv', index=False)
diverging = LinearSegmentedColormap.from_list('age_contribution', ['#2678AC', '#FAFAFA', '#B6553C'])
fig = plt.figure(figsize=(6.3, 4.2))
grid = fig.add_gridspec(2, 2, width_ratios=[1, .027], hspace=.19, wspace=.045)
limit = np.abs(contribution.to_numpy()).max()
for row, (frame, cmap, vmin, vmax, barlabel) in enumerate([
    (usage, 'Blues', 0, 1, 'Class weight'),
    (contribution, diverging, -limit, limit, 'Score contribution'),
]):
    ax, bar = fig.add_subplot(grid[row, 0]), fig.add_subplot(grid[row, 1])
    artist = ax.imshow(frame, aspect='auto', interpolation='nearest', cmap=cmap, vmin=vmin, vmax=vmax)
    ax.set_yticks(np.arange(len(top)), top.index.astype(str))
    ax.set_ylabel('Program', labelpad=5)
    ax.set_xticks(np.arange(24), class_labels if row else [''] * 24,
                  rotation=90, fontsize=7.5)
    ax.tick_params(length=0, pad=3)
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.axvline(19.5, color='#A4ABB1', lw=.65)
    colourbar = fig.colorbar(artist, cax=bar)
    colourbar.set_label(barlabel, fontsize=9, labelpad=5)
    colourbar.ax.tick_params(labelsize=8, length=2, width=.5)
fig.subplots_adjust(left=.075, right=.875, top=.985, bottom=.285)
save(fig, 'fig_senescence_class_contributions')

# Full-data points and omission ranges; the ranges are not confidence intervals.
influence = read_csv(SOURCE / 'candidate_mouse_influence.csv', index_col='program_id').loc[top.index]
fig, axes = plt.subplots(1, 2, figsize=(6.3, 2.8), sharey=True,
                         gridspec_kw={'width_ratios': [1.3, 1]})
positions = np.arange(len(top))
for ax, column in zip(axes, ['score', 'rank']):
    style(ax, 'x')
    ax.hlines(positions, influence[f'minimum_{column}'], influence[f'maximum_{column}'],
              color=COLORS['identity'], lw=1.5, zorder=2)
    ax.scatter(influence[column], positions, color=COLORS['identity'], s=25,
                edgecolors='white', lw=.5, zorder=3)
    ax.set_yticks(positions, top.index.astype(str))
axes[0].invert_yaxis()
axes[0].set(xlabel='Normalised age-tail score', ylabel='Program', xlim=(0, None))
axes[0].axvline(0, color='#929BA3', ls=(0, (3, 3)), lw=.7)
axes[1].set_xlabel('Rank among all programs')
axes[1].xaxis.set_major_locator(MaxNLocator(integer=True, nbins=7))
axes[1].axvline(len(top), color='#929BA3', ls=(0, (3, 3)), lw=.7)
fig.tight_layout(pad=.6, w_pad=.8)
save(fig, 'fig_senescence_mouse_influence')

# The original 50-bin mixtures, with exact zero masses displayed separately.
# Removing zeros from the first bin changes neither bin edges nor positive mass.
histograms = read_csv(SOURCE / 'candidate_activity_histograms.csv')
zeros = read_csv(SOURCE / 'candidate_zero_fractions.csv')
zeros = zeros.loc[zeros.class_name.eq('mixture')].set_index(['program_id', 'age']).zero_fraction
fig, axes = plt.subplots(4, 2, figsize=(6.3, 7.1))
positive_histograms = []
for ax, program in zip(axes.flat, top.index):
    style(ax)
    all_positive = []
    zero_values = []
    for age, color in AGE_COLORS.items():
        data = histograms.loc[histograms.program_id.eq(program) & histograms.age.eq(age)].sort_values('bin_index')
        assert data.bin_index.tolist() == list(range(50))
        edges = np.r_[data.bin_left.to_numpy(), data.bin_right.iloc[-1]]
        mass = data.bin_fraction.to_numpy().copy()
        zero = float(zeros.loc[program, age])
        np.testing.assert_allclose(mass.sum(), 1, atol=1e-10)
        mass[0] -= zero
        assert mass.min() >= -1e-12
        mass = np.maximum(mass, 0)
        np.testing.assert_allclose(mass.sum() + zero, 1, atol=1e-10)
        positive_histograms.append(data.assign(positive_bin_fraction=mass, exact_zero_fraction=zero))
        # A zero-probability bin has no height on a logarithmic axis.
        ax.stairs(np.where(mass > 0, mass, np.nan), edges, baseline=None,
                  color=color, lw=1.0, label=age.capitalize())
        all_positive.extend(mass[mass > 0])
        zero_values.append(zero)
    lower = 10.0 ** np.floor(np.log10(min(all_positive)))
    upper = 10.0 ** np.ceil(np.log10(max(all_positive)))
    ax.set(yscale='log', xlim=(0, edges[-1]), ylim=(lower, upper * 1.35),
            xlabel='Program activity', ylabel=f'{program}\nCell fraction per bin')
    ax.yaxis.set_major_locator(LogLocator(base=10, numticks=4))
    ax.yaxis.set_minor_locator(NullLocator())
    ax.xaxis.set_major_locator(MaxNLocator(nbins=4, min_n_ticks=3))
    ax.text(.98, .97, f'Zero: {zero_values[0]:.1%} / {zero_values[1]:.1%}',
             transform=ax.transAxes, ha='right', va='top', fontsize=7)
axes[-1, -1].axis('off')
axes[-1, -1].legend(handles=[Line2D([0], [0], color=color, lw=1.3, label=age.capitalize())
                           for age, color in AGE_COLORS.items()],
                    loc='center', frameon=False, fontsize=10)
fig.tight_layout(pad=.6, h_pad=1.0, w_pad=1.0)
save(fig, 'fig_senescence_activity_distributions')
pd.concat(positive_histograms, ignore_index=True).to_csv(
    OUTPUT / 'positive_activity_histogram_data.csv', index=False)

(OUTPUT / 'figure_manifest.json').write_text(json.dumps({
    'source_selection': str(SOURCE.relative_to(ROOT)), 'program_ids': top.index.tolist(),
    'input_sha256': INPUTS, 'script_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    'matplotlib_version': matplotlib.__version__, 'numpy_version': np.__version__,
    'pandas_version': pd.__version__,
    'distribution_display': 'Exact zero mass annotated young / aged; unrenormalised positive-bin mass on log y; all 50 original bins and full x range retained.',
}, indent=2), encoding='utf-8')
print(f'Saved four PDF/PNG figure pairs to {OUTPUT}')
