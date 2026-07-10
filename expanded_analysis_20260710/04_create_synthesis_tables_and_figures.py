#!/usr/bin/env python3
"""Create synthesis tables and a compact main/supplementary figure set."""

from __future__ import annotations

import argparse
import math
from pathlib import Path

import anndata as ad
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize, TwoSlopeNorm
from matplotlib.lines import Line2D
from matplotlib.patches import Circle, Ellipse, FancyArrowPatch, FancyBboxPatch
import numpy as np
import pandas as pd
from PIL import Image, ImageDraw, ImageFont, JpegImagePlugin  # noqa: F401


COLORS = {
    "blue": "#0072B2",
    "orange": "#D55E00",
    "green": "#009E73",
    "purple": "#7B61A8",
    "sky": "#56B4E9",
    "yellow": "#E69F00",
    "pink": "#CC79A7",
    "gray": "#6B7280",
    "light_gray": "#E5E7EB",
    "dark": "#202124",
}

MODULE_ORDER = [
    "Transporter_entry",
    "Regulatory_bridge",
    "Cytoskeletal_endpoint",
    "Disulfidptosis_core",
    "PPP_NADPH",
    "Glutathione_redox",
    "WRC_actin_branching",
    "Focal_adhesion_tension",
    "Actomyosin_endpoint",
    "Ferroptosis_comparator",
    "Apoptosis_comparator",
    "Pyroptosis_comparator",
    "Inflammation_context",
    "BBB_endothelial_context",
]

MODULE_LABELS = {
    "Transporter_entry": "Cystine\ntransport",
    "Regulatory_bridge": "Regulatory\nbridge",
    "Cytoskeletal_endpoint": "Cytoskeletal\nendpoint",
    "Disulfidptosis_core": "Disulfidptosis\ncore",
    "PPP_NADPH": "PPP /\nNADPH",
    "Glutathione_redox": "Glutathione\nredox",
    "WRC_actin_branching": "WRC actin\nbranching",
    "Focal_adhesion_tension": "Focal adhesion\ntension",
    "Actomyosin_endpoint": "Actomyosin\nendpoint",
    "Ferroptosis_comparator": "Ferroptosis",
    "Apoptosis_comparator": "Apoptosis",
    "Pyroptosis_comparator": "Pyroptosis",
    "Inflammation_context": "Inflammation",
    "BBB_endothelial_context": "BBB /\nendothelium",
}

HUMAN_CONTRASTS = [
    ("GSE209552", "acute_severe_TBI_vs_control", "Acute severe TBI"),
    ("GSE193407", "continuous_CTE_stage", "CTE stage trend"),
    ("GSE193407", "late_CTE_stage3_4_vs_stage0", "Late CTE vs stage 0"),
    ("GSE319253", "CTE_vs_control_donor_collapsed", "Independent CTE"),
    ("GSE223245", "ordinal_peripheral_severity", "Peripheral severity"),
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--bulk-dir", required=True, type=Path)
    parser.add_argument("--snrna-dir", required=True, type=Path)
    parser.add_argument("--h5ad", required=True, type=Path)
    parser.add_argument("--base-fig1", required=True, type=Path)
    parser.add_argument("--figure-dir", required=True, type=Path)
    parser.add_argument("--table-dir", required=True, type=Path)
    return parser.parse_args()


def configure_style() -> None:
    mpl.rcParams.update(
        {
            "font.family": "Arial",
            "font.size": 8,
            "axes.titlesize": 9,
            "axes.labelsize": 8,
            "xtick.labelsize": 7,
            "ytick.labelsize": 7,
            "axes.linewidth": 0.8,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "legend.fontsize": 7,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "figure.facecolor": "white",
            "axes.facecolor": "white",
            "savefig.facecolor": "white",
        }
    )


def save_figure(fig: plt.Figure, outdir: Path, stem: str) -> None:
    outdir.mkdir(parents=True, exist_ok=True)
    fig.savefig(outdir / f"{stem}.pdf", bbox_inches="tight")
    fig.savefig(outdir / f"{stem}.png", dpi=300, bbox_inches="tight")
    tiff_path = outdir / f"{stem}.tif"
    fig.savefig(
        tiff_path,
        dpi=600,
        bbox_inches="tight",
        pil_kwargs={"compression": "tiff_lzw"},
    )
    with Image.open(tiff_path) as image:
        rgba = image.convert("RGBA")
        white = Image.new("RGBA", rgba.size, "white")
        rgb = Image.alpha_composite(white, rgba).convert("RGB")
        rgb.save(tiff_path, compression="tiff_lzw", dpi=(600, 600))
    plt.close(fig)


def add_panel_label(ax: plt.Axes, label: str) -> None:
    ax.text(-0.13, 1.08, label, transform=ax.transAxes, fontsize=11, fontweight="bold", va="bottom")


def find_font(size: int, bold: bool = False) -> ImageFont.FreeTypeFont | ImageFont.ImageFont:
    name = "arialbd.ttf" if bold else "arial.ttf"
    path = Path("C:/Windows/Fonts") / name
    if path.exists():
        return ImageFont.truetype(str(path), size=size)
    return ImageFont.load_default()


def draw_rodent(draw: ImageDraw.ImageDraw, center: tuple[int, int], scale: float, color: str, impact: bool) -> None:
    x, y = center
    body = [x - int(55 * scale), y - int(22 * scale), x + int(40 * scale), y + int(25 * scale)]
    draw.ellipse(body, fill="#F7F7F7", outline=color, width=max(2, int(3 * scale)))
    draw.ellipse(
        [x + int(25 * scale), y - int(26 * scale), x + int(63 * scale), y + int(10 * scale)],
        fill="#F7F7F7",
        outline=color,
        width=max(2, int(3 * scale)),
    )
    draw.ellipse(
        [x + int(42 * scale), y - int(32 * scale), x + int(55 * scale), y - int(19 * scale)],
        fill="#F7F7F7",
        outline=color,
        width=max(1, int(2 * scale)),
    )
    draw.ellipse([x + int(51 * scale), y - int(14 * scale), x + int(56 * scale), y - int(9 * scale)], fill="#202124")
    tail = [(x - int(52 * scale), y + int(10 * scale)), (x - int(82 * scale), y + int(18 * scale)), (x - int(90 * scale), y - int(3 * scale))]
    draw.line(tail, fill=color, width=max(2, int(3 * scale)), joint="curve")
    if impact:
        draw.line(
            [(x + int(5 * scale), y - int(70 * scale)), (x + int(5 * scale), y - int(31 * scale))],
            fill=COLORS["orange"],
            width=max(3, int(5 * scale)),
        )
        draw.polygon(
            [
                (x - int(4 * scale), y - int(38 * scale)),
                (x + int(14 * scale), y - int(38 * scale)),
                (x + int(5 * scale), y - int(25 * scale)),
            ],
            fill=COLORS["orange"],
        )


def make_figure1(old_path: Path, outdir: Path) -> None:
    old = Image.open(old_path).convert("RGB")
    width, height = old.size
    strip_h = int(height * 0.22)
    canvas = Image.new("RGB", (width, height + strip_h), "white")
    canvas.paste(old, (0, 0))
    draw = ImageDraw.Draw(canvas)
    title_font = find_font(max(22, width // 145), bold=True)
    body_font = find_font(max(17, width // 185), bold=False)
    small_font = find_font(max(15, width // 215), bold=False)
    draw.text((30, height + 18), "C", font=find_font(max(28, width // 115), bold=True), fill=COLORS["dark"])
    draw.text((85, height + 22), "Expanded temporal and heterogeneity validation", font=title_font, fill=COLORS["dark"])

    boxes = [
        (int(width * 0.34), COLORS["purple"], "7", "Repeated closed-head injury", "GSE298240  |  sex and p38-inhibitor strata", False),
        (int(width * 0.72), COLORS["green"], "8", "Rat fluid-percussion injury", "GSE111452  |  24 h to 12 mo", True),
    ]
    box_w, box_h = int(width * 0.34), int(strip_h * 0.66)
    top = height + int(strip_h * 0.25)
    for cx, color, number, title, subtitle, impact in boxes:
        left = cx - box_w // 2
        draw.rounded_rectangle(
            [left, top, left + box_w, top + box_h],
            radius=18,
            fill="#FBFBFD",
            outline=color,
            width=3,
        )
        draw.ellipse([left + 16, top + 14, left + 56, top + 54], fill=color)
        number_font = find_font(24, bold=True)
        number_width = draw.textbbox((0, 0), number, font=number_font)[2]
        draw.text((left + 36 - number_width / 2, top + 20), number, font=number_font, fill="white")
        draw_rodent(draw, (left + 165, top + box_h // 2 + 8), 0.8, color, impact)
        draw.text((left + 285, top + 18), title, font=body_font, fill=COLORS["dark"])
        draw.text((left + 285, top + 54), subtitle, font=small_font, fill=color)
        line_y = top + box_h - 30
        draw.line([(left + 290, line_y), (left + box_w - 32, line_y)], fill="#B9BEC8", width=3)
        for fraction in (0.0, 0.33, 0.66, 1.0):
            px = int(left + 290 + fraction * (box_w - 322))
            draw.ellipse([px - 6, line_y - 6, px + 6, line_y + 6], fill=color)

    canvas.save(outdir / "Fig1_Study_design_and_hypothesis.tif", compression="tiff_lzw", dpi=(600, 600))
    canvas.save(outdir / "Fig1_Study_design_and_hypothesis.png", dpi=(300, 300))
    canvas.save(outdir / "Fig1_Study_design_and_hypothesis.pdf", resolution=600)


def signed_evidence(frame: pd.DataFrame) -> pd.Series:
    q = frame["FDR"].clip(lower=1e-6).fillna(1.0)
    return np.sign(frame["effect"].fillna(0)) * np.minimum(-np.log10(q), 4)


def human_figure(module_effects: pd.DataFrame, core: pd.DataFrame, scores: pd.DataFrame, coupling: pd.DataFrame, outdir: Path) -> None:
    selected_modules = [
        "Transporter_entry",
        "Regulatory_bridge",
        "Cytoskeletal_endpoint",
        "Disulfidptosis_core",
        "PPP_NADPH",
        "Glutathione_redox",
        "Ferroptosis_comparator",
        "Inflammation_context",
    ]
    fig, axes = plt.subplots(2, 2, figsize=(7.2, 6.5), gridspec_kw={"height_ratios": [1.05, 0.95]})
    ax = axes[0, 0]
    matrix = np.full((len(HUMAN_CONTRASTS), len(selected_modules)), np.nan)
    for i, (dataset, contrast, _) in enumerate(HUMAN_CONTRASTS):
        sub = module_effects[(module_effects.dataset == dataset) & (module_effects.contrast == contrast)].set_index("module")
        for j, module in enumerate(selected_modules):
            if module in sub.index:
                row = sub.loc[[module]].iloc[0]
                matrix[i, j] = np.sign(row.effect) * min(-math.log10(max(row.FDR, 1e-6)), 4)
    image = ax.imshow(matrix, cmap="RdBu_r", norm=TwoSlopeNorm(vmin=-4, vcenter=0, vmax=4), aspect="auto")
    ax.set_xticks(range(len(selected_modules)), [MODULE_LABELS[x] for x in selected_modules], rotation=45, ha="right")
    ax.set_yticks(range(len(HUMAN_CONTRASTS)), [x[2] for x in HUMAN_CONTRASTS])
    ax.set_title("Human module evidence across complementary designs")
    for i, (dataset, contrast, _) in enumerate(HUMAN_CONTRASTS):
        sub = module_effects[(module_effects.dataset == dataset) & (module_effects.contrast == contrast)].set_index("module")
        for j, module in enumerate(selected_modules):
            if module in sub.index and sub.loc[[module]].iloc[0].FDR < 0.05:
                ax.text(j, i, "*", ha="center", va="center", color="white" if abs(matrix[i, j]) > 2 else "black", fontsize=9)
    cbar = fig.colorbar(image, ax=ax, fraction=0.035, pad=0.02)
    cbar.set_label("Signed -log10(FDR)")
    add_panel_label(ax, "A")

    ax = axes[0, 1]
    genes = ["SLC3A2", "SLC7A11", "WASF2", "TLN1", "ACTB", "MYH9", "MYL6", "FLNA"]
    gene_contrasts = HUMAN_CONTRASTS[:4]
    plot_rows = []
    for yi, (dataset, contrast, label) in enumerate(gene_contrasts):
        sub = core[(core.dataset == dataset) & (core.contrast == contrast)].set_index("gene_symbol")
        for xi, gene in enumerate(genes):
            if gene not in sub.index:
                continue
            row = sub.loc[[gene]].iloc[0]
            plot_rows.append((xi, yi, row.logFC, row.FDR, gene, label))
    color_limit = max(abs(x[2]) for x in plot_rows if np.isfinite(x[2]))
    for xi, yi, effect, q, _, _ in plot_rows:
        strength = min(-math.log10(max(float(q), 1e-6)), 4)
        ax.scatter(xi, yi, s=20 + 35 * strength, c=[effect], cmap="RdBu_r", norm=TwoSlopeNorm(vmin=-color_limit, vcenter=0, vmax=color_limit), edgecolor="white", linewidth=0.5)
    ax.set_xticks(range(len(genes)), genes, rotation=45, ha="right")
    ax.set_yticks(range(len(gene_contrasts)), [x[2] for x in gene_contrasts])
    ax.set_xlim(-0.6, len(genes) - 0.4)
    ax.set_ylim(len(gene_contrasts) - 0.4, -0.6)
    ax.grid(color="#ECEFF3", linewidth=0.6)
    ax.set_axisbelow(True)
    ax.set_title("Eight-gene direction and inferential strength")
    scalar = mpl.cm.ScalarMappable(norm=TwoSlopeNorm(vmin=-color_limit, vcenter=0, vmax=color_limit), cmap="RdBu_r")
    cbar = fig.colorbar(scalar, ax=ax, fraction=0.035, pad=0.02)
    cbar.set_label("log2FC or stage coefficient")
    add_panel_label(ax, "B")

    ax = axes[1, 0]
    stage_scores = scores[(scores.dataset == "GSE193407") & (scores.analysis_contrast == "continuous_CTE_stage")]
    modules = ["Transporter_entry", "Regulatory_bridge", "Cytoskeletal_endpoint"]
    module_colors = [COLORS["blue"], COLORS["orange"], COLORS["green"]]
    for module, color in zip(modules, module_colors):
        sub = stage_scores[stage_scores.module == module].copy()
        summary = sub.groupby("stage").score.agg(["mean", "sem"]).reset_index()
        ax.errorbar(summary.stage, summary["mean"], yerr=1.96 * summary["sem"], marker="o", color=color, linewidth=1.5, capsize=2, label=MODULE_LABELS[module].replace("\n", " "))
    ax.axhline(0, color="#B8BEC8", linewidth=0.7)
    ax.set_xlabel("CTE stage")
    ax.set_ylabel("Mean standardized module score (95% CI)")
    ax.set_xticks([0, 1, 2, 3, 4])
    ax.legend(frameon=False, loc="upper left")
    ax.set_title("Stage-resolved chronic human trajectories")
    add_panel_label(ax, "C")

    ax = axes[1, 1]
    coupling = coupling.sort_values("adjusted_beta")
    y = np.arange(len(coupling))
    ax.barh(y, coupling.adjusted_beta, color=COLORS["blue"], alpha=0.82)
    ax.scatter(coupling.stage_interaction_beta, y, marker="D", s=23, color=COLORS["orange"], label="Transport x stage")
    for yi, row in enumerate(coupling.itertuples()):
        if row.FDR < 0.05:
            ax.text(row.adjusted_beta + (0.025 if row.adjusted_beta >= 0 else -0.025), yi, "*", ha="left" if row.adjusted_beta >= 0 else "right", va="center")
    ax.axvline(0, color="#60646C", linewidth=0.8)
    ax.set_yticks(y, [MODULE_LABELS.get(x, x).replace("\n", " ") for x in coupling.outcome_module])
    ax.set_xlabel("Adjusted standardized coefficient")
    ax.legend(frameon=False, loc="lower right")
    ax.set_title("Transporter-module coupling, adjusted for stage, age and RIN")
    add_panel_label(ax, "D")
    fig.tight_layout(w_pad=1.2, h_pad=2.0)
    save_figure(fig, outdir, "Fig2_Human_cross_cohort_evidence")


def matrix_heatmap(ax: plt.Axes, frame: pd.DataFrame, row_key: str, col_key: str, value: str, q_key: str, row_order: list[str], col_order: list[str], title: str, vlim: float | None = None) -> mpl.image.AxesImage:
    pivot = frame.pivot_table(index=row_key, columns=col_key, values=value, aggfunc="first").reindex(index=row_order, columns=col_order)
    q = frame.pivot_table(index=row_key, columns=col_key, values=q_key, aggfunc="first").reindex(index=row_order, columns=col_order)
    values = pivot.to_numpy(dtype=float)
    if vlim is None:
        finite = np.abs(values[np.isfinite(values)])
        vlim = np.quantile(finite, 0.95) if finite.size else 1
    vlim = max(float(vlim), 0.1)
    image = ax.imshow(values, cmap="RdBu_r", norm=TwoSlopeNorm(vmin=-vlim, vcenter=0, vmax=vlim), aspect="auto")
    ax.set_xticks(range(len(col_order)), [MODULE_LABELS.get(x, x) for x in col_order], rotation=45, ha="right")
    ax.set_yticks(range(len(row_order)), row_order)
    for i in range(len(row_order)):
        for j in range(len(col_order)):
            if i < q.shape[0] and j < q.shape[1] and np.isfinite(q.iloc[i, j]) and q.iloc[i, j] < 0.05:
                ax.text(j, i, "*", ha="center", va="center", color="black", fontsize=8)
    ax.set_title(title)
    return image


def snrna_figure(snrna_dir: Path, h5ad_path: Path, outdir: Path) -> None:
    localization = pd.read_csv(snrna_dir / "Table_snRNA_8gene_localization.csv")
    effects = pd.read_csv(snrna_dir / "Table_snRNA_donor_level_condition_effects.csv")
    coupling = pd.read_csv(snrna_dir / "Table_snRNA_module_coupling.csv")
    adata = ad.read_h5ad(h5ad_path, backed="r")
    umap = np.asarray(adata.obsm["X_umap"])
    cell_types = adata.obs["cell_type"].astype(str).to_numpy()
    palette = {
        "Astrocyte": COLORS["purple"],
        "Endothelial": COLORS["sky"],
        "Immune": COLORS["pink"],
        "Microglia": COLORS["yellow"],
        "Neuron": COLORS["green"],
        "OPC": "#F0C95A",
        "Oligodendrocyte": "#8BCF5B",
    }

    fig, axes = plt.subplots(2, 2, figsize=(7.2, 6.4))
    ax = axes[0, 0]
    for cell_type in palette:
        mask = cell_types == cell_type
        ax.scatter(umap[mask, 0], umap[mask, 1], s=1.1, color=palette[cell_type], alpha=0.7, rasterized=True, label=cell_type)
    ax.set_xlabel("UMAP1")
    ax.set_ylabel("UMAP2")
    ax.set_title("Annotated acute TBI and control nuclei")
    ax.legend(frameon=False, ncol=2, markerscale=4, loc="best")
    add_panel_label(ax, "A")

    ax = axes[0, 1]
    genes = ["SLC3A2", "SLC7A11", "WASF2", "TLN1", "ACTB", "MYH9", "MYL6", "FLNA"]
    types = ["Astrocyte", "Endothelial", "Microglia", "Neuron", "OPC", "Oligodendrocyte"]
    loc = localization[localization.cell_type.isin(types)].copy()
    norm = Normalize(vmin=0, vmax=max(loc.median_donor_mean_expression.max(), 0.1))
    for row in loc.itertuples():
        xi, yi = genes.index(row.gene_symbol), types.index(row.cell_type)
        ax.scatter(xi, yi, s=15 + 240 * row.mean_detection_fraction, c=[row.median_donor_mean_expression], cmap="Blues", norm=norm, edgecolor="#D5D9E0", linewidth=0.4)
    ax.set_xticks(range(len(genes)), genes, rotation=45, ha="right")
    ax.set_yticks(range(len(types)), types)
    ax.set_ylim(len(types) - 0.5, -0.5)
    ax.set_title("Donor-aggregated localization of the eight genes")
    cbar = fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap="Blues"), ax=ax, fraction=0.035, pad=0.02)
    cbar.set_label("Median donor mean expression")
    ax.legend(
        handles=[plt.scatter([], [], s=15 + 240 * f, facecolor="white", edgecolor="#60646C", label=f"{int(f*100)}%") for f in (0.1, 0.3, 0.5)],
        title="Detected nuclei",
        frameon=False,
        loc="lower right",
    )
    add_panel_label(ax, "B")

    ax = axes[1, 0]
    modules = ["Transporter_entry", "Regulatory_bridge", "Cytoskeletal_endpoint", "Disulfidptosis_core", "PPP_NADPH", "Glutathione_redox"]
    types_effect = [x for x in types if x in effects.cell_type.unique()]
    image = matrix_heatmap(ax, effects, "cell_type", "module", "hedges_g", "FDR", types_effect, modules, "Exploratory donor-level TBI-control effects", vlim=2)
    cbar = fig.colorbar(image, ax=ax, fraction=0.035, pad=0.02)
    cbar.set_label("Hedges g")
    for row in effects.itertuples():
        if row.cell_type in types_effect and row.module in modules and row.permutation_p < 0.05 and row.FDR >= 0.05:
            ax.text(modules.index(row.module), types_effect.index(row.cell_type), "+", ha="center", va="center", fontsize=9)
    add_panel_label(ax, "C")

    ax = axes[1, 1]
    pairs = [
        ("Transporter_entry", "Regulatory_bridge", "Transport -> bridge"),
        ("Transporter_entry", "Cytoskeletal_endpoint", "Transport -> cytoskeleton"),
        ("PPP_NADPH", "Disulfidptosis_core", "PPP -> disulfidptosis"),
        ("Glutathione_redox", "Disulfidptosis_core", "GSH -> disulfidptosis"),
        ("WRC_actin_branching", "Actomyosin_endpoint", "WRC -> actomyosin"),
        ("Focal_adhesion_tension", "Actomyosin_endpoint", "Adhesion -> actomyosin"),
    ]
    coupling = coupling.copy()
    coupling["pair"] = coupling.apply(lambda r: next((label for source, target, label in pairs if r.source_module == source and r.target_module == target), None), axis=1)
    coupling = coupling[coupling.pair.notna()]
    types_coupling = [x for x in types if x in coupling.cell_type.unique()]
    pair_order = [x[2] for x in pairs]
    image = matrix_heatmap(ax, coupling, "cell_type", "pair", "partial_spearman_rho", "FDR", types_coupling, pair_order, "Condition-, abundance- and capture-adjusted coupling", vlim=1)
    ax.set_xticklabels(pair_order, rotation=45, ha="right")
    cbar = fig.colorbar(image, ax=ax, fraction=0.035, pad=0.02)
    cbar.set_label("Partial Spearman rho")
    add_panel_label(ax, "D")
    fig.tight_layout(w_pad=1.4, h_pad=2.1)
    save_figure(fig, outdir, "Fig3_Donor_level_single_nucleus_context")


def animal_figure(module_effects: pd.DataFrame, interactions: pd.DataFrame, outdir: Path) -> None:
    modules = ["Transporter_entry", "Regulatory_bridge", "Cytoskeletal_endpoint", "Disulfidptosis_core", "PPP_NADPH", "Glutathione_redox", "Inflammation_context", "Pyroptosis_comparator"]
    fig, axes = plt.subplots(2, 2, figsize=(7.2, 7.0), gridspec_kw={"height_ratios": [0.9, 1.1]})

    ax = axes[0, 0]
    sub = module_effects[module_effects.dataset == "GSE163415"].copy()
    order = [
        "3DPI_Hipp_TBI_vs_NoTBI_adjusted_for_treatment",
        "3DPI_Hypo_TBI_vs_NoTBI_adjusted_for_treatment",
        "3DPI_Thal_TBI_vs_NoTBI_adjusted_for_treatment",
        "29DPI_Hipp_TBI_vs_NoTBI_adjusted_for_treatment",
        "29DPI_Hypo_TBI_vs_NoTBI_adjusted_for_treatment",
        "29DPI_Thal_TBI_vs_NoTBI_adjusted_for_treatment",
    ]
    labels = ["3 d Hipp", "3 d Hypo", "3 d Thal", "29 d Hipp", "29 d Hypo", "29 d Thal"]
    sub["row_label"] = sub.contrast.map(dict(zip(order, labels)))
    image = matrix_heatmap(ax, sub, "row_label", "module", "effect", "FDR", labels, modules, "Mouse CCI time-region module effects", vlim=4)
    cbar = fig.colorbar(image, ax=ax, fraction=0.035, pad=0.02)
    cbar.set_label("Hedges g")
    add_panel_label(ax, "A")

    ax = axes[0, 1]
    sub = module_effects[module_effects.dataset == "GSE111452"].copy()
    rat_order = [
        "GPL15084_Hippocampus_24h_TBI_vs_sham",
        "GPL22740_Cortex_24h_TBI_vs_sham",
        "GPL22740_Cortex_2wk_TBI_vs_sham",
        "GPL22740_Hippocampus_2wk_TBI_vs_sham",
        "GPL15084_Hippocampus_3mo_TBI_vs_sham",
        "GPL22740_Cortex_3mo_TBI_vs_sham",
        "GPL22740_Cortex_6mo_TBI_vs_sham",
        "GPL22740_Hippocampus_6mo_TBI_vs_sham",
        "GPL22740_Cortex_12mo_TBI_vs_sham",
        "GPL22740_Hippocampus_12mo_TBI_vs_sham",
    ]
    rat_labels = ["24 h Hipp", "24 h Cortex", "2 wk Cortex", "2 wk Hipp", "3 mo Hipp", "3 mo Cortex", "6 mo Cortex", "6 mo Hipp", "12 mo Cortex", "12 mo Hipp"]
    sub["row_label"] = sub.contrast.map(dict(zip(rat_order, rat_labels)))
    image = matrix_heatmap(ax, sub, "row_label", "module", "effect", "FDR", rat_labels, modules, "Rat FPI trajectory from 24 h to 12 mo", vlim=3)
    cbar = fig.colorbar(image, ax=ax, fraction=0.035, pad=0.02)
    cbar.set_label("Hedges g")
    add_panel_label(ax, "B")

    ax = axes[1, 0]
    sub = module_effects[(module_effects.dataset == "GSE298240") & (module_effects.contrast == "chronic_5xCHI_vs_sham_adjusted")].set_index("module").reindex(modules)
    colors = [COLORS["orange"] if value >= 0 else COLORS["blue"] for value in sub.effect]
    y = np.arange(len(sub))
    ax.barh(y, sub.effect, color=colors, alpha=0.86)
    for yi, row in enumerate(sub.itertuples()):
        if row.FDR < 0.05:
            ax.text(row.effect + (0.025 if row.effect >= 0 else -0.025), yi, "*", ha="left" if row.effect >= 0 else "right", va="center")
    ax.axvline(0, color="#666B73", linewidth=0.8)
    ax.set_yticks(y, [MODULE_LABELS[x].replace("\n", " ") for x in sub.index])
    ax.invert_yaxis()
    ax.set_xlabel("Adjusted injury coefficient")
    ax.set_title("Chronic repetitive injury adjusted for sex and p38 inhibitor")
    add_panel_label(ax, "C")

    ax = axes[1, 1]
    interaction_modules = interactions[interactions.module.notna()].copy()
    display_rows = []
    for row in interaction_modules.itertuples():
        if row.module not in modules:
            continue
        if row.dataset == "GSE163415":
            label = row.contrast.replace("_treatment_by_injury", "").replace("3DPI_", "3 d ").replace("29DPI_", "29 d ").replace("Hipp", "Hipp").replace("Hypo", "Hypo").replace("Thal", "Thal")
        elif row.dataset == "GSE298240":
            label = "5xCHI " + row.contrast.replace("_by_injury", " x injury")
        else:
            continue
        display_rows.append({"row_label": label, "module": row.module, "effect": row.effect, "FDR": row.FDR})
    interaction_plot = pd.DataFrame(display_rows)
    row_order = list(dict.fromkeys(interaction_plot.row_label.tolist()))
    image = matrix_heatmap(ax, interaction_plot, "row_label", "module", "effect", "FDR", row_order, modules, "Treatment and sex interaction sensitivity", vlim=0.8)
    cbar = fig.colorbar(image, ax=ax, fraction=0.035, pad=0.02)
    cbar.set_label("Interaction coefficient")
    add_panel_label(ax, "D")
    fig.tight_layout(w_pad=1.4, h_pad=2.0)
    save_figure(fig, outdir, "Fig4_Rodent_spatiotemporal_and_heterogeneity_analysis")


def supplementary_figures(module_effects: pd.DataFrame, enrichment: pd.DataFrame, core: pd.DataFrame, outdir: Path) -> None:
    selected = enrichment[enrichment.pathway.isin(MODULE_ORDER)].copy()
    selected["row_label"] = selected.dataset + " | " + selected.contrast.str.replace("_", " ", regex=False)
    row_order = list(dict.fromkeys(selected.row_label.tolist()))
    fig, ax = plt.subplots(figsize=(7.2, max(5.5, len(row_order) * 0.22)))
    image = matrix_heatmap(ax, selected, "row_label", "pathway", "NES", "padj", row_order, MODULE_ORDER, "Ranked gene-set enrichment across all bulk contrasts", vlim=2.5)
    fig.colorbar(image, ax=ax, fraction=0.018, pad=0.01, label="Normalized enrichment score")
    fig.tight_layout()
    save_figure(fig, outdir, "FigS1_Full_ranked_module_enrichment")

    core = core.copy()
    core["row_label"] = core.dataset + " | " + core.contrast.str.replace("_", " ", regex=False)
    row_order = list(dict.fromkeys(core.row_label.tolist()))
    genes = ["SLC3A2", "SLC7A11", "WASF2", "TLN1", "ACTB", "MYH9", "MYL6", "FLNA"]
    fig, ax = plt.subplots(figsize=(7.2, max(5.5, len(row_order) * 0.22)))
    image = matrix_heatmap(ax, core, "row_label", "gene_symbol", "logFC", "FDR", row_order, genes, "Complete eight-gene effect map", vlim=1.5)
    fig.colorbar(image, ax=ax, fraction=0.018, pad=0.01, label="log2FC or model coefficient")
    fig.tight_layout()
    save_figure(fig, outdir, "FigS2_Full_eight_gene_effect_map")

    boundary = module_effects[
        ((module_effects.dataset == "GSE223245") & module_effects.module.isin(MODULE_ORDER))
    ].copy()
    remote = core[core.dataset == "GSE104687"].copy()
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.0))
    ax = axes[0]
    boundary = boundary.set_index("module").reindex(MODULE_ORDER[:9])
    ax.barh(np.arange(len(boundary)), boundary.effect, color=[COLORS["orange"] if x >= 0 else COLORS["blue"] for x in boundary.effect])
    ax.axvline(0, color="#666B73", linewidth=0.8)
    ax.set_yticks(np.arange(len(boundary)), [MODULE_LABELS[x].replace("\n", " ") for x in boundary.index])
    ax.invert_yaxis()
    ax.set_xlabel("Spearman r with peripheral severity")
    ax.set_title("Peripheral severity boundary")
    add_panel_label(ax, "A")
    ax = axes[1]
    remote["region"] = remote.contrast.str.extract(r"sample_(.*?)_TBI")
    pivot = remote.pivot_table(index="region", columns="gene_symbol", values="logFC", aggfunc="first").reindex(columns=["SLC3A2", "SLC7A11", "WASF2", "TLN1", "ACTB", "MYH9", "MYL6", "FLNA"])
    image = ax.imshow(pivot, cmap="RdBu_r", norm=TwoSlopeNorm(vmin=-0.03, vcenter=0, vmax=0.03), aspect="auto")
    ax.set_xticks(range(pivot.shape[1]), pivot.columns, rotation=45, ha="right")
    ax.set_yticks(range(pivot.shape[0]), pivot.index)
    ax.set_title("Remote-region human boundary")
    fig.colorbar(image, ax=ax, fraction=0.04, pad=0.02, label="log2FC")
    add_panel_label(ax, "B")
    fig.tight_layout(w_pad=1.5)
    save_figure(fig, outdir, "FigS3_Peripheral_and_remote_negative_boundaries")


def synthesis_tables(module_effects: pd.DataFrame, enrichment: pd.DataFrame, outdir: Path) -> None:
    outdir.mkdir(parents=True, exist_ok=True)
    evidence = module_effects.merge(
        enrichment[["dataset", "contrast", "pathway", "NES", "padj", "leadingEdge"]],
        left_on=["dataset", "contrast", "module"],
        right_on=["dataset", "contrast", "pathway"],
        how="left",
    )
    evidence["signed_log10_FDR"] = signed_evidence(evidence)
    evidence.to_csv(outdir / "Table_cross_dataset_module_evidence.csv", index=False)

    recurrence = []
    for module, frame in evidence.groupby("module", sort=False):
        recurrence.append(
            {
                "module": module,
                "module_class": frame.module_class.iloc[0],
                "n_bulk_contrasts": len(frame),
                "n_datasets": frame.dataset.nunique(),
                "n_positive_effects": int((frame.effect > 0).sum()),
                "n_effect_FDR_lt_0_05": int((frame.FDR < 0.05).sum()),
                "n_positive_enrichment_FDR_lt_0_05": int(((frame.NES > 0) & (frame.padj < 0.05)).sum()),
                "n_negative_enrichment_FDR_lt_0_05": int(((frame.NES < 0) & (frame.padj < 0.05)).sum()),
            }
        )
    pd.DataFrame(recurrence).to_csv(outdir / "Table_module_recurrence_descriptive.csv", index=False)

    claims = pd.DataFrame(
        [
            ["Chronic human transporter signal", "GSE193407 adjusted stage models and donor-collapsed GSE319253", "SLC3A2/SLC7A11 and transporter-module RNA associations recur", "No protein transport or cystine-flux measurement"],
            ["Acute cytoskeletal response", "GSE209552 ranked enrichment and GSE163415 3-d hippocampal replication", "WRC, focal-adhesion and actomyosin programs are prominent early", "Acute human cohort is small and tissue sources differ"],
            ["Long-horizon regional heterogeneity", "GSE111452 24-h to 12-mo rat FPI and GSE163415 3/29-d mouse CCI", "Support is strongest in early/subacute hippocampal windows and is not spatially uniform", "Cross-species models and platforms are not interchangeable"],
            ["Cellular context", "GSE209552 donor-by-cell-type aggregation", "The genes localize mainly to glial/vascular contexts and show exploratory module coupling", "Only three control donors; no condition effect survives FDR correction"],
            ["Mechanistic specificity", "Ferroptosis, apoptosis, pyroptosis and inflammatory comparator modules", "The candidate axis co-occurs with broader injury and cell-death responses", "Transcriptomics cannot establish disulfidptosis as the operative death mechanism"],
        ],
        columns=["claim", "supporting_analysis", "bounded_interpretation", "remaining_gap"],
    )
    claims.to_csv(outdir / "Table_key_claims_and_inferential_limits.csv", index=False)


def main() -> None:
    args = parse_args()
    configure_style()
    args.figure_dir.mkdir(parents=True, exist_ok=True)
    main_dir = args.figure_dir / "main"
    supp_dir = args.figure_dir / "supplementary"
    main_dir.mkdir(parents=True, exist_ok=True)
    supp_dir.mkdir(parents=True, exist_ok=True)

    module_effects = pd.read_csv(args.bulk_dir / "Table_module_effects.csv")
    enrichment = pd.read_csv(args.bulk_dir / "Table_ranked_module_enrichment.csv")
    core = pd.read_csv(args.bulk_dir / "Table_core_8gene_effects.csv")
    scores = pd.read_csv(args.bulk_dir / "Table_sample_module_scores.csv.gz")
    coupling = pd.read_csv(args.bulk_dir / "GSE193407_transport_actin_coupling.csv")
    interactions = pd.read_csv(args.bulk_dir / "Table_sex_treatment_interactions.csv")

    make_figure1(args.base_fig1, main_dir)
    human_figure(module_effects, core, scores, coupling, main_dir)
    snrna_figure(args.snrna_dir, args.h5ad, main_dir)
    animal_figure(module_effects, interactions, main_dir)
    supplementary_figures(module_effects, enrichment, core, supp_dir)
    synthesis_tables(module_effects, enrichment, args.table_dir)


if __name__ == "__main__":
    main()
