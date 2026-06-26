from __future__ import annotations

from pathlib import Path

import anndata as ad
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns


PRIORITY_GENES = ["SLC3A2", "SLC7A11", "WASF2", "TLN1", "ACTB", "MYH9", "MYL6", "FLNA"]
CORE_FOCUS_GENES = ["SLC3A2", "SLC7A11", "WASF2", "TLN1"]
CELL_TYPE_COLORS = {
    "Neuron": "#7BCFB6",
    "Astrocyte": "#B7A6F6",
    "Microglia": "#F6B26B",
    "Oligodendrocyte": "#A9D86E",
    "OPC": "#F6D36D",
    "Endothelial": "#79B8E8",
    "Pericyte": "#D6B3E6",
    "Immune": "#F4A3C4",
    "Ependymal": "#A7D7C5",
    "Unknown": "#D9D9D9",
}
CELL_TYPE_SHORT = {
    "Astrocyte": "Ast",
    "Endothelial": "Endo",
    "Immune": "Imm",
    "Microglia": "Mg",
    "Neuron": "Neu",
    "Oligodendrocyte": "Oligo",
    "OPC": "OPC",
    "Pericyte": "Peri",
    "Ependymal": "Epen",
    "Unknown": "Unk",
}
DELTA_BAR_COLOR = "#F29B8C"
EXPR_CMAP = mpl.colors.LinearSegmentedColormap.from_list(
    "expr_pastel_blue",
    ["#F7F5F0", "#DDD8FF", "#A9B5FF", "#5E72E4"],
)
SCORE_CMAP = mpl.colors.LinearSegmentedColormap.from_list(
    "score_soft_lilac",
    ["#FFF9FB", "#F3D9EF", "#D8C1F0", "#9CA7EA", "#6D81D8"],
)


PROJECT_ROOT = Path(__file__).resolve().parents[1]
WORKDIR = PROJECT_ROOT / "Phase3_深化优化与最终报告_20260506-0513" / "11_双硫死亡聚焦论文设计_20260604"
FIG_DIR = WORKDIR / "figures"
TABLE_DIR = WORKDIR / "tables"
REPORT_DIR = WORKDIR / "reports"
H5AD_PATH = WORKDIR / "raw_data_v3_20260604" / "GSE209552_true_snRNA_processed_20260613.h5ad"


def configure_style() -> None:
    mpl.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Microsoft YaHei", "Arial", "DejaVu Sans"],
            "axes.unicode_minus": False,
            "font.size": 9,
            "axes.titlesize": 10,
            "axes.labelsize": 9,
            "xtick.labelsize": 8,
            "ytick.labelsize": 8,
            "legend.fontsize": 8,
            "figure.dpi": 160,
            "savefig.dpi": 320,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )
    sns.set_theme(style="white", font="Microsoft YaHei")


def raw_gene_df(adata: ad.AnnData, genes: list[str]) -> pd.DataFrame:
    valid = [g for g in genes if g in adata.raw.var_names]
    raw = adata.raw[:, valid].X
    if hasattr(raw, "toarray"):
        raw = raw.toarray()
    return pd.DataFrame(raw, index=adata.obs_names, columns=valid)


def load_inputs() -> tuple[ad.AnnData, pd.DataFrame, pd.DataFrame]:
    adata = ad.read_h5ad(H5AD_PATH)
    gene_summary = pd.read_csv(TABLE_DIR / "GSE209552_true_scRNA_8gene_celltype_localization_20260613.csv")
    delta_df = pd.read_csv(TABLE_DIR / "GSE209552_true_scRNA_8gene_donor_level_deltas_20260613.csv")
    return adata, gene_summary, delta_df


def make_main_figure(adata: ad.AnnData, gene_summary: pd.DataFrame, delta_df: pd.DataFrame) -> None:
    fig = plt.figure(figsize=(14.5, 10))
    gs = fig.add_gridspec(
        3,
        4,
        height_ratios=[0.9, 1.7, 1.35],
        width_ratios=[1.15, 1.2, 1.0, 1.35],
        hspace=0.45,
        wspace=0.35,
    )

    ax_a = fig.add_subplot(gs[0, :])
    ax_a.axis("off")
    total_cells = adata.n_obs
    sample_n = adata.obs["sample_accession"].nunique()
    donor_n = adata.obs["donor_id"].nunique()
    control_n = adata.obs.loc[adata.obs["condition"] == "Control", "sample_accession"].nunique()
    tbi_n = adata.obs.loc[adata.obs["condition"] == "TBI", "sample_accession"].nunique()
    text = (
        "Real snRNA-seq localization workflow (GSE209552)\n"
        f"Human acute severe TBI resected tissue | {total_cells:,} nuclei after QC | {sample_n} snRNA samples | {donor_n} donors\n"
        f"Control samples={control_n}, TBI samples={tbi_n} | Read 10x-like matrices from GEO RAW.tar -> QC -> PCA/UMAP -> Leiden -> marker annotation -> 8-gene localization\n"
        "Interpretive boundary: true nucleus-level localization is shown here, but TBI vs control cell-type comparisons remain exploratory because control donors include repeated brain regions."
    )
    ax_a.text(0.01, 0.95, text, va="top", ha="left", fontsize=10)

    ax_b = fig.add_subplot(gs[1, 0:2])
    palette = [CELL_TYPE_COLORS.get(ct, "#bdbdbd") for ct in adata.obs["cell_type"].cat.categories]
    sc.pl.umap(
        adata,
        color="cell_type",
        ax=ax_b,
        show=False,
        palette=palette,
        frameon=False,
        size=25,
        alpha=0.92,
        legend_loc="on data",
        legend_fontsize=8,
        legend_fontoutline=2,
        title="UMAP of annotated nuclei",
    )

    right_gs = gs[1, 2:4].subgridspec(1, 2, width_ratios=[0.95, 1.25], wspace=0.42)

    ax_c = fig.add_subplot(right_gs[0, 0])
    comp = (
        adata.obs.groupby(["condition", "cell_type"], observed=True).size().rename("n").reset_index()
        .pivot(index="condition", columns="cell_type", values="n")
        .fillna(0)
    )
    comp = comp.div(comp.sum(axis=1), axis=0) * 100
    bottom = np.zeros(len(comp))
    ordered_cols = [c for c in adata.obs["cell_type"].cat.categories if c in comp.columns]
    for col in ordered_cols:
        vals = comp[col].to_numpy()
        ax_c.bar(
            comp.index,
            vals,
            bottom=bottom,
            color=CELL_TYPE_COLORS.get(col, "#bdbdbd"),
            edgecolor="white",
            linewidth=0.5,
            label=col,
        )
        bottom += vals
    ax_c.set_ylabel("Cell fraction (%)")
    ax_c.set_title("Cell-type composition by condition")
    ax_c.tick_params(axis="x", rotation=0)
    ax_c.legend(
        loc="upper center",
        bbox_to_anchor=(0.5, -0.16),
        ncol=2,
        frameon=False,
        fontsize=7,
        columnspacing=0.8,
        handletextpad=0.4,
    )

    ax_d = fig.add_subplot(right_gs[0, 1])
    top_delta = delta_df.loc[delta_df["gene"].isin(CORE_FOCUS_GENES)].head(8).copy()
    top_delta["label"] = top_delta["gene"] + " | " + top_delta["cell_type"].map(lambda x: CELL_TYPE_SHORT.get(x, x))
    ax_d.barh(
        top_delta["label"][::-1],
        top_delta["delta_tbi_minus_control"][::-1],
        color=DELTA_BAR_COLOR,
        edgecolor="#B86E62",
        linewidth=0.35,
    )
    ax_d.axvline(0, color="black", linewidth=0.8)
    ax_d.set_xlabel("Donor-level mean log1p delta")
    ax_d.set_title("Core-gene TBI-control clues")
    ax_d.tick_params(axis="y", labelsize=7)
    ax_d.grid(axis="x", color="#d1d5db", linewidth=0.5, alpha=0.7)

    ax_e = fig.add_subplot(gs[2, 0:2])
    dot = (
        gene_summary.groupby(["cell_type", "gene"])
        .agg(mean_log1p_expr=("mean_log1p_expr", "mean"), pct_expressed=("pct_expressed", "mean"))
        .reset_index()
    )
    cell_order = [c for c in adata.obs["cell_type"].cat.categories if c != "Unknown"]
    if "Unknown" in adata.obs["cell_type"].cat.categories:
        cell_order.append("Unknown")
    gene_order = [g for g in PRIORITY_GENES if g in dot["gene"].unique()]
    expr_mat = dot.pivot(index="cell_type", columns="gene", values="mean_log1p_expr").reindex(cell_order)[gene_order]
    pct_mat = dot.pivot(index="cell_type", columns="gene", values="pct_expressed").reindex(cell_order)[gene_order]
    vmax = float(np.nanmax(expr_mat.values))
    for yi, cell_type in enumerate(cell_order):
        for xi, gene in enumerate(gene_order):
            expr_val = expr_mat.loc[cell_type, gene]
            pct_val = pct_mat.loc[cell_type, gene]
            if pd.isna(expr_val) or pd.isna(pct_val):
                continue
            ax_e.scatter(
                xi,
                yi,
                s=max(12, pct_val * 2.2),
                c=[expr_val],
                cmap=EXPR_CMAP,
                vmin=0,
                vmax=vmax,
                edgecolor="#A9A9A9",
                linewidth=0.2,
            )
    ax_e.set_xticks(range(len(gene_order)))
    ax_e.set_xticklabels(gene_order, rotation=45, ha="right")
    ax_e.set_yticks(range(len(cell_order)))
    ax_e.set_yticklabels(cell_order)
    ax_e.set_xlim(-0.5, len(gene_order) - 0.5)
    ax_e.set_ylim(-0.5, len(cell_order) - 0.5)
    ax_e.invert_yaxis()
    ax_e.set_title("8-gene real localization dot plot")
    ax_e.grid(False)

    ax_f = fig.add_subplot(gs[2, 2:4])
    module_df = raw_gene_df(adata, PRIORITY_GENES)
    score = module_df[[g for g in PRIORITY_GENES if g in module_df.columns]].mean(axis=1)
    ax_f.scatter(
        adata.obsm["X_umap"][:, 0],
        adata.obsm["X_umap"][:, 1],
        c=score.reindex(adata.obs_names).to_numpy(),
        s=10,
        cmap=SCORE_CMAP,
        linewidths=0,
        alpha=0.95,
    )
    ax_f.set_title("Cell-level 8-gene localization score")
    ax_f.set_xlabel("UMAP1")
    ax_f.set_ylabel("UMAP2")
    sns.despine(ax=ax_f, left=False, bottom=False)

    for label, ax in zip(list("ABCDEF"), [ax_a, ax_b, ax_c, ax_d, ax_e, ax_f]):
        if ax is ax_a:
            ax.text(-0.01, 1.02, label, transform=ax.transAxes, fontsize=12, fontweight="bold", va="top", ha="left")
        else:
            ax.text(-0.12, 1.06, label, transform=ax.transAxes, fontsize=12, fontweight="bold", va="top", ha="left")

    fig.savefig(FIG_DIR / "Fig3b_GSE209552_true_snRNA_localization_panel_20260613.png", bbox_inches="tight", dpi=320)
    fig.savefig(FIG_DIR / "Fig3b_GSE209552_true_snRNA_localization_panel_20260613.pdf", bbox_inches="tight")
    plt.close(fig)


def make_feature_figure(adata: ad.AnnData) -> None:
    genes = [g for g in PRIORITY_GENES if g in adata.raw.var_names]
    raw_df = raw_gene_df(adata, genes)
    ncols = 4
    nrows = int(np.ceil(len(genes) / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(16, 3.8 * nrows))
    axes = np.atleast_1d(axes).reshape(nrows, ncols)

    for ax, gene in zip(axes.ravel(), genes):
        vals = raw_df[gene].reindex(adata.obs_names).to_numpy()
        sc_handle = ax.scatter(
            adata.obsm["X_umap"][:, 0],
            adata.obsm["X_umap"][:, 1],
            c=vals,
            s=8,
            cmap=EXPR_CMAP,
            linewidths=0,
            alpha=0.95,
        )
        ax.set_title(gene)
        ax.set_xlabel("UMAP1")
        ax.set_ylabel("UMAP2")
        plt.colorbar(sc_handle, ax=ax, fraction=0.046, pad=0.03)

    for ax in axes.ravel()[len(genes):]:
        ax.axis("off")

    fig.suptitle("GSE209552 true snRNA feature maps for the 8 disulfidptosis-priority genes", y=1.02, fontsize=12)
    for r in range(nrows):
        for c in range(ncols):
            ax = axes[r, c]
            if r < nrows - 1:
                ax.set_xlabel("")
            if c > 0:
                ax.set_ylabel("")
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    fig.savefig(FIG_DIR / "Fig3c_GSE209552_true_snRNA_feature_maps_8gene_20260613.png", bbox_inches="tight", dpi=320)
    fig.savefig(FIG_DIR / "Fig3c_GSE209552_true_snRNA_feature_maps_8gene_20260613.pdf", bbox_inches="tight")
    plt.close(fig)


def write_reports(adata: ad.AnnData, gene_summary: pd.DataFrame, delta_df: pd.DataFrame) -> None:
    cell_counts = adata.obs["cell_type"].value_counts().rename_axis("cell_type").reset_index(name="n_cells")
    core_localization = (
        gene_summary.loc[gene_summary["gene"].isin(CORE_FOCUS_GENES)]
        .groupby(["gene", "cell_type"])
        .agg(mean_log1p_expr=("mean_log1p_expr", "mean"), pct_expressed=("pct_expressed", "mean"))
        .reset_index()
        .sort_values(["gene", "mean_log1p_expr", "pct_expressed"], ascending=[True, False, False])
        .groupby("gene")
        .head(2)
    )

    loc_lines = []
    for gene in CORE_FOCUS_GENES:
        sub = core_localization.loc[core_localization["gene"] == gene]
        if sub.empty:
            continue
        pairs = [
            f"{row.cell_type} (mean log1p {row.mean_log1p_expr:.2f}, pct {row.pct_expressed:.1f}%)"
            for row in sub.itertuples(index=False)
        ]
        loc_lines.append(f"- {gene}: " + "; ".join(pairs))

    top_delta = delta_df.loc[delta_df["gene"].isin(CORE_FOCUS_GENES)].head(8)
    delta_lines = [
        f"- {row.gene} | {row.cell_type}: donor-level delta={row.delta_tbi_minus_control:.3f} "
        f"(TBI mean {row.tbi_mean:.3f} vs Control mean {row.control_mean:.3f}; "
        f"TBI donors={int(row.tbi_donor_n)}, Control donors={int(row.control_donor_n)})"
        for row in top_delta.itertuples(index=False)
    ]

    report = "\n".join(
        [
            "# GSE209552 真实单核定位补充分析",
            "",
            "本次补充分析直接使用 GSE209552 的 GEO RAW.tar 中 10x 风格单核矩阵，重建真实 snRNA-seq 对象，而不是继续使用 bulk-like proxy。",
            "",
            "## 数据概况",
            f"- QC 后保留 nuclei: {adata.n_obs}",
            f"- 聚类数: {adata.obs['leiden'].nunique()}",
            f"- 样本数: {adata.obs['sample_accession'].nunique()}",
            f"- donor 数: {adata.obs['donor_id'].nunique()}",
            "",
            "## 主要细胞类型",
            *[f"- {row.cell_type}: {row.n_cells} nuclei" for row in cell_counts.itertuples(index=False)],
            "",
            "## 核心基因的单核定位线索",
            *loc_lines,
            "",
            "## donor 层面的探索性 TBI-Control 差值线索",
            *delta_lines,
            "",
            "## 方法边界",
            "- 本次结果属于真实单核定位和细胞类型表达层面的补强。",
            "- 若进一步比较 TBI 与 Control 的细胞类型差异，必须注意 control 侧存在同一 donor 的多个脑区样本，正式统计比较应以 donor-level 或 pseudobulk 为主。",
            "- 因此本次图谱中的条件差异更适合作为探索性定位线索，单细胞定位本身则可以作为真实证据使用。",
        ]
    )
    (REPORT_DIR / "GSE209552_true_snRNA_localization_report_20260613.md").write_text(report, encoding="utf-8")

    insertion_note = "\n".join(
        [
            "# 图文插入建议",
            "",
            "## 建议插入位置",
            "- 建议作为双硫死亡机制图前后的单细胞补强结果，替代仅凭 bulk-like proxy 的细胞类型推断。",
            "",
            "## 图注建议（Fig3b）",
            "Fig3b. Real snRNA-seq localization of disulfidptosis-priority genes in acute severe human TBI (GSE209552). GEO RAW.tar 10x-like matrices were reprocessed to obtain 19,627 nuclei after QC from 17 snRNA samples and 15 donors. UMAP-based clustering recovered neurons, oligodendrocytes, astrocytes, OPCs, microglia, endothelial cells, and immune cells. The dot plot and cell-level score map show true nucleus-level localization of the 8 priority genes, while donor-level TBI-control cell-type differences are displayed only as exploratory clues because control donors include repeated brain regions.",
            "",
            "## 图注建议（Fig3c）",
            "Fig3c. UMAP feature maps of the 8 disulfidptosis-priority genes in GSE209552 snRNA-seq. SLC7A11 shows the clearest astrocyte-endothelial localization background, whereas SLC3A2, WASF2, and TLN1 display broader glial and vascular enrichment patterns.",
            "",
            "## 结果段落建议",
            "基于 GSE209552 原始单核测序矩阵的重分析，我们获得了 19,627 个通过质控的细胞核，并稳定分辨出神经元、少突胶质细胞、星形胶质细胞、OPCs、小胶质细胞、内皮细胞和免疫细胞等主要群体。与此前 bulk-like proxy 推断不同，本次结果提供了真实单核层面的细胞定位证据：SLC7A11 主要定位于星形胶质细胞，并在内皮细胞中也显示明显表达；SLC3A2、WASF2 与 TLN1 则更偏向胶质/血管相关细胞群，其中 WASF2 在小胶质细胞中的信号尤为突出。需要指出的是，图中 TBI 与 Control 的细胞类型差异仅作为探索性线索呈现，而单核定位本身可作为本文更坚实的细胞层面证据。",
        ]
    )
    (REPORT_DIR / "GSE209552_true_snRNA_insertion_text_20260613.md").write_text(insertion_note, encoding="utf-8")


def main() -> None:
    configure_style()
    adata, gene_summary, delta_df = load_inputs()
    make_main_figure(adata, gene_summary, delta_df)
    make_feature_figure(adata)
    write_reports(adata, gene_summary, delta_df)
    print("Finalized figure and report assets.")


if __name__ == "__main__":
    main()
