#!/usr/bin/env python3
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

BASE = Path('/mnt/data/submission_work/Neurochemical_Research_#U76f4#U63a5#U6295#U7a3f#U5305_20260623_v5/08_Data_and_analysis_files/01_processed_result_tables')
OUT = Path('/mnt/data/optimized_submission_work')
OUT.mkdir(parents=True, exist_ok=True)
priority = pd.read_csv(BASE/'v3_integrated_8gene_priority_20260604.csv')
components = {
    'Human FDR': priority['human_FDR_comparisons'].astype(float),
    'Human nominal': priority['human_nominal_comparisons'].astype(float),
    'Acute brain effect': priority['GSE209552_acute_logFC'].clip(lower=0).astype(float),
    'Mouse recurrence': priority['mouse_nominal_units'].astype(float),
    'Peripheral |r|': priority['GSE223245_severity_r_all'].abs().astype(float),
}
# min-max normalize each component; zero if constant
norm = pd.DataFrame({'gene_symbol': priority['gene_symbol']})
for name, series in components.items():
    lo, hi = np.nanmin(series), np.nanmax(series)
    if not np.isfinite(lo) or not np.isfinite(hi) or hi == lo:
        norm[name] = 0.0
    else:
        norm[name] = (series - lo) / (hi - lo)
component_names = list(components.keys())
norm['mean_normalized_score'] = norm[component_names].mean(axis=1)
# ranks under leave-one-layer-out; rank 1 best, deterministic tie handling
rank_table = pd.DataFrame({'gene_symbol': norm['gene_symbol']})
rank_table['Full'] = norm['mean_normalized_score'].rank(method='min', ascending=False).astype(int)
for name in component_names:
    score = norm[[c for c in component_names if c != name]].mean(axis=1)
    rank_table['No ' + name] = score.rank(method='min', ascending=False).astype(int)

out_table = priority.merge(norm, on='gene_symbol').merge(rank_table, on='gene_symbol')
out_table.to_csv(OUT/'Table_S9_Robustness_sensitivity.csv', index=False)

# Sort by base score for display
order = norm.sort_values('mean_normalized_score', ascending=False)['gene_symbol'].tolist()
norm_plot = norm.set_index('gene_symbol').loc[order]
rank_plot = rank_table.set_index('gene_symbol').loc[order]

plt.rcParams.update({
    'font.size': 9,
    'axes.titlesize': 11,
    'axes.labelsize': 9,
    'xtick.labelsize': 8,
    'ytick.labelsize': 8,
})
fig = plt.figure(figsize=(8.5, 8.4), constrained_layout=False)
grid = fig.add_gridspec(2, 2, width_ratios=[1.12, 1.0], height_ratios=[1.0, 1.05], wspace=0.55, hspace=0.72)
fig.suptitle('Online Resource 3. Robustness audit for the eight-gene transporter-actin panel', fontsize=14, fontweight='bold', y=0.985)

# A normalized component heatmap
ax1 = fig.add_subplot(grid[0,0])
im1 = ax1.imshow(norm_plot[component_names].values, aspect='auto', vmin=0, vmax=1, cmap='YlGnBu')
ax1.set_xticks(np.arange(len(component_names)), labels=component_names, rotation=35, ha='right')
ax1.set_yticks(np.arange(len(order)), labels=order)
ax1.set_title('A. Normalized dry-evidence components', loc='left', fontweight='bold')
for i in range(len(order)):
    for j in range(len(component_names)):
        v = norm_plot.iloc[i][component_names[j]]
        ax1.text(j, i, f'{v:.2f}', ha='center', va='center', fontsize=6.5)
cb1 = fig.colorbar(im1, ax=ax1, fraction=0.046, pad=0.025)
cb1.set_label('normalized support')

# B bar chart
ax2 = fig.add_subplot(grid[0,1])
y = np.arange(len(order))
scores = norm_plot['mean_normalized_score'].values
ax2.barh(y, scores)
ax2.set_yticks(y, labels=order)
ax2.invert_yaxis()
ax2.set_xlim(0, 1)
ax2.set_xlabel('mean normalized evidence')
ax2.set_title('B. Integrated robustness score', loc='left', fontweight='bold')
for yi, s in zip(y, scores):
    ax2.text(min(s + 0.02, 0.98), yi, f'{s:.2f}', va='center', fontsize=8)

# C rank sensitivity heatmap (1 best)
ax3 = fig.add_subplot(grid[1,0])
rank_cols = list(rank_plot.columns)
# custom reversed colormap so low ranks are darker
im3 = ax3.imshow(rank_plot[rank_cols].values, aspect='auto', vmin=1, vmax=8, cmap='viridis_r')
ax3.set_xticks(np.arange(len(rank_cols)), labels=rank_cols, rotation=35, ha='right')
ax3.set_yticks(np.arange(len(order)), labels=order)
ax3.set_title('C. Leave-one-evidence-layer sensitivity', loc='left', fontweight='bold')
for i in range(len(order)):
    for j in range(len(rank_cols)):
        r = int(rank_plot.iloc[i][rank_cols[j]])
        ax3.text(j, i, str(r), ha='center', va='center', fontsize=7)
cb3 = fig.colorbar(im3, ax=ax3, fraction=0.046, pad=0.025)
cb3.set_label('rank (1 = highest)')

# D text boundary
ax4 = fig.add_subplot(grid[1,1])
ax4.axis('off')
ax4.set_title('D. Interpretation boundary', loc='left', fontweight='bold')
boundary = (
    'This audit reweights only the public-data results already generated for the manuscript.\n\n'
    'It does not add wet-lab evidence, diagnose TBI, or prove disulfidptotic cell death.\n\n'
    'Candidates that retain favorable ranks after removal of one evidence layer are better suited for validation planning. In this audit, transporter-regulatory candidates remain comparatively stable, whereas cytoskeletal endpoint genes remain important acute/mechanical-stress readouts.'
)
ax4.text(0.0, 0.95, boundary, ha='left', va='top', wrap=True, fontsize=9)

fig.subplots_adjust(bottom=0.12)
fig.savefig(OUT/'Online_Resource_3_FigS1_Robustness_Audit.png', dpi=300, bbox_inches='tight')
fig.savefig(OUT/'Online_Resource_3_FigS1_Robustness_Audit.tif', dpi=600, bbox_inches='tight', pil_kwargs={'compression':'tiff_lzw'})
