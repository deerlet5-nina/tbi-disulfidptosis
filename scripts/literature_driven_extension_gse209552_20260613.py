from __future__ import annotations

from pathlib import Path

import anndata as ad
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


PROJECT_ROOT = Path(__file__).resolve().parents[1]
WORKDIR = PROJECT_ROOT / "Phase3_深化优化与最终报告_20260506-0513" / "11_双硫死亡聚焦论文设计_20260604"
FIG_DIR = WORKDIR / "figures"
TABLE_DIR = WORKDIR / "tables"
REPORT_DIR = WORKDIR / "reports"
H5AD_PATH = WORKDIR / "raw_data_v3_20260604" / "GSE209552_true_snRNA_processed_20260613.h5ad"

EXTENSION_GENES = ["ACTN4", "MYH10", "DSTN", "NCKAP1", "TLN1"]
INFLAMMATORY_GENES = ["NFKBIA", "IL1B", "TNF", "CCL2", "HMOX1", "FOS", "JUN"]
CORE_GENES = ["SLC7A11", "SLC3A2", "WASF2", "TLN1"]
ALL_PLOT_GENES = CORE_GENES + [g for g in EXTENSION_GENES if g not in CORE_GENES] + INFLAMMATORY_GENES
CELL_ORDER = ["Astrocyte", "Endothelial", "Immune", "Microglia", "Neuron", "OPC", "Oligodendrocyte"]
CONDITION_COLORS = {"Control": "#D7DCEA", "TBI": "#F2B7AA"}
EXPR_CMAP = mpl.colors.LinearSegmentedColormap.from_list(
    "expr_pastel_blue",
    ["#F7F5F0", "#DDD8FF", "#A9B5FF", "#5E72E4"],
)


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


def summarize_gene_localization(adata: ad.AnnData, gene_df: pd.DataFrame, genes: list[str]) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for cell_type in CELL_ORDER:
        idx = adata.obs["cell_type"] == cell_type
        if idx.sum() == 0:
            continue
        for gene in genes:
            if gene not in gene_df.columns:
                continue
            vals = gene_df.loc[idx, gene]
            rows.append(
                {
                    "cell_type": cell_type,
                    "gene": gene,
                    "n_cells": int(len(vals)),
                    "mean_expr": float(vals.mean()),
                    "pct_expressed": float((vals > 0).mean() * 100),
                }
            )
    return pd.DataFrame(rows)


def donor_level_deltas(adata: ad.AnnData, gene_df: pd.DataFrame, genes: list[str]) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    obs = adata.obs[["donor_id", "condition", "cell_type"]].copy()
    for gene in genes:
        if gene not in gene_df.columns:
            continue
        tmp = pd.concat([obs, gene_df[[gene]]], axis=1)
        donor_means = (
            tmp.groupby(["donor_id", "condition", "cell_type"], observed=True)[gene]
            .mean()
            .reset_index()
        )
        for cell_type, sub in donor_means.groupby("cell_type", observed=True):
            tbi = sub.loc[sub["condition"] == "TBI", gene]
            ctrl = sub.loc[sub["condition"] == "Control", gene]
            if len(tbi) == 0 or len(ctrl) == 0:
                continue
            rows.append(
                {
                    "gene": gene,
                    "cell_type": cell_type,
                    "tbi_donor_n": int(len(tbi)),
                    "control_donor_n": int(len(ctrl)),
                    "tbi_mean": float(tbi.mean()),
                    "control_mean": float(ctrl.mean()),
                    "delta_tbi_minus_control": float(tbi.mean() - ctrl.mean()),
                }
            )
    return pd.DataFrame(rows).sort_values("delta_tbi_minus_control", ascending=False)


def inflammation_module_summary(adata: ad.AnnData, gene_df: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    valid = [g for g in INFLAMMATORY_GENES if g in gene_df.columns]
    module_score = gene_df[valid].mean(axis=1)
    obs = adata.obs[["donor_id", "condition", "cell_type"]].copy()
    obs["inflammation_score"] = module_score.reindex(adata.obs_names).to_numpy()
    obs["SLC7A11"] = gene_df["SLC7A11"].reindex(adata.obs_names).to_numpy()

    celltype_condition = (
        obs.groupby(["cell_type", "condition"], observed=True)
        .agg(
            inflammation_score=("inflammation_score", "mean"),
            slc7a11_mean=("SLC7A11", "mean"),
            n_cells=("inflammation_score", "size"),
        )
        .reset_index()
    )

    donor_level = (
        obs.groupby(["donor_id", "condition", "cell_type"], observed=True)
        .agg(
            inflammation_score=("inflammation_score", "mean"),
            slc7a11_mean=("SLC7A11", "mean"),
            n_cells=("inflammation_score", "size"),
        )
        .reset_index()
    )

    corr_rows: list[dict[str, object]] = []
    for cell_type, sub in donor_level.groupby("cell_type", observed=True):
        corr_rows.append(
            {
                "cell_type": cell_type,
                "donor_n": int(len(sub)),
                "slc7a11_inflammation_corr": float(sub["slc7a11_mean"].corr(sub["inflammation_score"])),
                "tbi_donor_n": int((sub["condition"] == "TBI").sum()),
                "control_donor_n": int((sub["condition"] == "Control").sum()),
                "mean_inflammation_score": float(sub["inflammation_score"].mean()),
                "mean_slc7a11": float(sub["slc7a11_mean"].mean()),
            }
        )
    corr_df = pd.DataFrame(corr_rows).sort_values("mean_inflammation_score", ascending=False)
    return celltype_condition, corr_df


def make_dotplot(ax: plt.Axes, df: pd.DataFrame, genes: list[str], title: str) -> None:
    plot_df = df.loc[df["gene"].isin(genes)].copy()
    plot_df["cell_type"] = pd.Categorical(plot_df["cell_type"], categories=CELL_ORDER, ordered=True)
    plot_df["gene"] = pd.Categorical(plot_df["gene"], categories=genes, ordered=True)
    plot_df = plot_df.sort_values(["cell_type", "gene"])
    vmax = float(plot_df["mean_expr"].max()) if len(plot_df) else 1.0
    for yi, cell_type in enumerate(CELL_ORDER):
        for xi, gene in enumerate(genes):
            sub = plot_df.loc[(plot_df["cell_type"] == cell_type) & (plot_df["gene"] == gene)]
            if sub.empty:
                continue
            row = sub.iloc[0]
            ax.scatter(
                xi,
                yi,
                s=max(12, float(row["pct_expressed"]) * 2.0),
                c=[float(row["mean_expr"])],
                cmap=EXPR_CMAP,
                vmin=0,
                vmax=vmax,
                edgecolor="#A9A9A9",
                linewidth=0.2,
            )
    ax.set_xticks(range(len(genes)))
    ax.set_xticklabels(genes, rotation=45, ha="right")
    ax.set_yticks(range(len(CELL_ORDER)))
    ax.set_yticklabels(CELL_ORDER)
    ax.set_xlim(-0.5, len(genes) - 0.5)
    ax.set_ylim(-0.5, len(CELL_ORDER) - 0.5)
    ax.invert_yaxis()
    ax.set_title(title)
    ax.grid(False)


def make_figure(
    ext_df: pd.DataFrame,
    inflam_df: pd.DataFrame,
    celltype_condition_df: pd.DataFrame,
    delta_df: pd.DataFrame,
) -> None:
    fig = plt.figure(figsize=(14, 10))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.1, 1.0], wspace=0.32, hspace=0.38)

    ax_a = fig.add_subplot(gs[0, 0])
    make_dotplot(ax_a, ext_df, EXTENSION_GENES, "Literature-expanded actin/remodeling genes")

    ax_b = fig.add_subplot(gs[0, 1])
    make_dotplot(ax_b, inflam_df, INFLAMMATORY_GENES, "Inflammatory-response genes linked to the TBI paper")

    ax_c = fig.add_subplot(gs[1, 0])
    bar_df = celltype_condition_df.loc[celltype_condition_df["cell_type"].isin(CELL_ORDER)].copy()
    bar_df["cell_type"] = pd.Categorical(bar_df["cell_type"], categories=CELL_ORDER, ordered=True)
    sns.barplot(
        data=bar_df,
        x="cell_type",
        y="inflammation_score",
        hue="condition",
        palette=CONDITION_COLORS,
        ax=ax_c,
    )
    ax_c.set_title("Inflammation module score by cell type")
    ax_c.set_xlabel("")
    ax_c.set_ylabel("Mean module score")
    ax_c.tick_params(axis="x", rotation=35)
    ax_c.legend(title="", frameon=False, loc="upper right")

    ax_d = fig.add_subplot(gs[1, 1])
    plot_delta = delta_df.head(10).copy()
    plot_delta["label"] = plot_delta["gene"] + " | " + plot_delta["cell_type"]
    ax_d.barh(
        plot_delta["label"][::-1],
        plot_delta["delta_tbi_minus_control"][::-1],
        color="#F2B7AA",
        edgecolor="#B86E62",
        linewidth=0.35,
    )
    ax_d.axvline(0, color="black", linewidth=0.8)
    ax_d.set_title("Top exploratory donor-level TBI-control clues")
    ax_d.set_xlabel("Mean log1p delta")
    ax_d.tick_params(axis="y", labelsize=7)
    ax_d.grid(axis="x", color="#d1d5db", linewidth=0.5, alpha=0.7)

    for label, ax in zip(list("ABCD"), [ax_a, ax_b, ax_c, ax_d]):
        ax.text(-0.12, 1.05, label, transform=ax.transAxes, fontsize=12, fontweight="bold", va="top", ha="left")

    fig.suptitle(
        "GSE209552 literature-driven extension: actin-remodeling expansion and SLC7A11-neuroinflammation context",
        y=0.995,
        fontsize=12,
    )
    fig.savefig(FIG_DIR / "Fig4_literature_driven_singlecell_extension_20260613.png", bbox_inches="tight", dpi=320)
    fig.savefig(FIG_DIR / "Fig4_literature_driven_singlecell_extension_20260613.pdf", bbox_inches="tight")
    plt.close(fig)


def write_report(
    ext_df: pd.DataFrame,
    inflam_df: pd.DataFrame,
    celltype_condition_df: pd.DataFrame,
    corr_df: pd.DataFrame,
    delta_df: pd.DataFrame,
) -> None:
    ext_top = (
        ext_df.sort_values(["gene", "mean_expr", "pct_expressed"], ascending=[True, False, False])
        .groupby("gene", observed=True)
        .head(2)
    )
    ext_lines = []
    for gene in EXTENSION_GENES:
        sub = ext_top.loc[ext_top["gene"] == gene]
        pairs = [f"{row.cell_type} (mean {row.mean_expr:.2f}, pct {row.pct_expressed:.1f}%)" for row in sub.itertuples(index=False)]
        ext_lines.append(f"- {gene}: " + "; ".join(pairs))

    inflam_top = (
        inflam_df.sort_values(["gene", "mean_expr", "pct_expressed"], ascending=[True, False, False])
        .groupby("gene", observed=True)
        .head(2)
    )
    inflam_lines = []
    for gene in INFLAMMATORY_GENES:
        sub = inflam_top.loc[inflam_top["gene"] == gene]
        pairs = [f"{row.cell_type} (mean {row.mean_expr:.2f}, pct {row.pct_expressed:.1f}%)" for row in sub.itertuples(index=False)]
        inflam_lines.append(f"- {gene}: " + "; ".join(pairs))

    module_lines = []
    module_df = celltype_condition_df.sort_values("inflammation_score", ascending=False)
    for row in module_df.itertuples(index=False):
        module_lines.append(
            f"- {row.cell_type} | {row.condition}: inflammation score={row.inflammation_score:.3f}, "
            f"SLC7A11 mean={row.slc7a11_mean:.3f}, n_cells={int(row.n_cells)}"
        )

    corr_lines = []
    for row in corr_df.itertuples(index=False):
        corr_lines.append(
            f"- {row.cell_type}: donor-level corr(SLC7A11, inflammation)={row.slc7a11_inflammation_corr:.3f}, "
            f"mean inflammation={row.mean_inflammation_score:.3f}, mean SLC7A11={row.mean_slc7a11:.3f}"
        )

    delta_lines = []
    for row in delta_df.head(10).itertuples(index=False):
        delta_lines.append(
            f"- {row.gene} | {row.cell_type}: delta={row.delta_tbi_minus_control:.3f} "
            f"(TBI {row.tbi_mean:.3f} vs Control {row.control_mean:.3f})"
        )

    report = "\n".join(
        [
            "# 文献驱动扩展分析：GSE209552 真实单核层面的机制升级",
            "",
            "## 参考文献",
            "- You H, Zou Q, Zhu T, et al. *Experimental Neurology* (2026): tDCS 通过抑制 SLC7A11 介导的 disulfidptosis 和 neuroinflammation 促进 TBI 大鼠功能恢复。DOI: 10.1016/j.expneurol.2026.115849. PubMed: https://pubmed.ncbi.nlm.nih.gov/42176875/",
            "- Bo X, Fan J, Xu S, et al. *Journal of Inflammation Research* (2025): 在 DCM 中用机器学习和单细胞分析提出 ACTN4/MYH10/TLN1/DSTN/NCKAP1 及免疫-纤维化耦合。DOI: 10.2147/JIR.S525114. URL: https://www.dovepress.com/disulfidptosis-related-genes-as-novel-biomarkers-and-therapeutic-targe-peer-reviewed-fulltext-article-JIR",
            "",
            "## 对当前 TBI 研究最重要的文献启发",
            "- TBI 文章把 SLC7A11、NQO1、NADPH/NADP+、氧化还原失衡、神经炎症和干预可逆性放在同一条机制链上，因此我们的讨论不应只盯住 8 基因表达，还应显式加入 NQO1/NADPH/炎症放大环。",
            "- DCM 文章虽然不是脑损伤研究，但它提供了一个可借鉴的“扩展骨架模块”：ACTN4、MYH10、DSTN、NCKAP1、TLN1，并且强调这些基因和免疫细胞浸润、组织重塑并不是分离现象。",
            "- 对我们来说，这两篇文献共同支持把当前模型从“转运入口-骨架终点”升级成“转运入口-氧化还原前提-炎症放大-骨架重塑”四层模型。",
            "",
            "## 扩展骨架基因在 GSE209552 单核图谱中的定位",
            *ext_lines,
            "",
            "## 炎症基因在 GSE209552 单核图谱中的定位",
            *inflam_lines,
            "",
            "## 炎症模块与 SLC7A11 的细胞类型关系",
            *module_lines,
            "",
            "## donor 层面的 SLC7A11-炎症耦合线索",
            *corr_lines,
            "",
            "## exploratory TBI-Control 差值线索",
            *delta_lines,
            "",
            "## 对当前文章叙事的直接升级建议",
            "- 机制主轴建议改写为：SLC7A11/SLC3A2 转运入口主要落在星形胶质细胞和内皮细胞；WASF2/TLN1 与扩展骨架基因把微胶质、内皮和神经元/OPC 的骨架重塑层连接起来；炎症放大则主要由微胶质和免疫细胞承担。",
            "- 因此，TBI 中更像是“跨细胞室耦合”的 disulfidptosis-like stress：astrocyte/endothelial 提供 transporter-redox 背景，microglia/immune 放大神经炎症，actin-remodeling 基因在血管、胶质和神经元不同程度分布。",
            "- 这也解释了为什么当前单核数据里 SLC7A11 并不在微胶质中最强，但微胶质却是 HMOX1/IL1B/TNF/WASF2/TLN1 更活跃的细胞群。",
            "",
            "## 对湿实验设计的升级建议",
            "- 在原有 SLC3A2/SLC7A11、WASF2/TLN1、ACTB/MYH9/MYL6/FLNA、F-actin、NADPH/GSH 之外，补入 NQO1、HMOX1、IL-1β、TNF-α、IBA1 作为炎症放大层读数。",
            "- 共定位优先顺序建议从单一 NeuN/GFAP 扩展为 CD31 + GFAP + IBA1 三轴优先，因为 transporter 背景更偏 astrocyte/endothelial，而炎症放大更偏 microglia/immune。",
            "- 如需做干预分支，tDCS 文章支持把 “SLC7A11/NQO1/NADPH + neuroinflammation” 作为一条可逆 readout 链；DCM 文章中的 quercetin/resveratrol 只能作为跨疾病启发，不应在本文中写成 TBI 已验证治疗。",
            "",
            "## 证据边界",
            "- 当前新增结果仍然属于真实单核定位和探索性 donor-level 差值线索，不能直接写作 TBI 已经证明存在和 DCM 相同的诊断基因模型。",
            "- DCM 文献中的机器学习诊断框架不能直接迁移到 TBI，但其中 ACTN4/MYH10/DSTN/NCKAP1/TLN1 这一“扩展骨架模块”可以作为我们的机制扩展层。",
        ]
    )
    (REPORT_DIR / "literature_driven_extension_report_20260613.md").write_text(report, encoding="utf-8")

    insertion = "\n".join(
        [
            "# 文献驱动插入段落",
            "",
            "## 结果段落建议",
            "参考近期 TBI 和双硫死亡相关疾病研究，我们进一步将 GSE209552 的真实单核图谱扩展到炎症放大和附加骨架重塑层。与 DCM 文献中 ACTN4/MYH10/TLN1/DSTN/NCKAP1 在纤维化和免疫浸润中的作用类似，这些扩展基因在 TBI 单核图谱中同样显示出明确的细胞偏好，但定位背景并不完全相同：ACTN4、DSTN 和 NCKAP1 更偏内皮/神经元相关区室，MYH10 主要见于神经元，而 TLN1 更偏微胶质和内皮细胞。另一方面，参考 Experimental Neurology 的 TBI 研究所强调的 SLC7A11-disulfidptosis-neuroinflammation 轴，我们观察到 SLC7A11 的表达背景主要位于星形胶质细胞和内皮细胞，而炎症基因模块则主要集中在微胶质和免疫细胞。这提示 TBI 中更可能存在一种跨细胞室耦合的 disulfidptosis-like stress，而非单一细胞类型内完成的简单同位表达事件。",
            "",
            "## 讨论段落建议",
            "这一文献驱动扩展使当前模型从“转运入口-肌动蛋白骨架终点”进一步升级为“转运入口-氧化还原前提-炎症放大-骨架重塑”四层结构。也就是说，astrocyte/endothelial 区室可能提供 SLC7A11/SLC3A2 相关的转运与还原力背景，microglia/immune 区室承担炎症放大，ACTN4/MYH10/DSTN/NCKAP1/TLN1/WASF2 等骨架重塑基因则把这些过程投射到血管、胶质和神经元的结构响应中。对后续验证而言，这意味着仅检测 F-actin 或仅检测 transporter 蛋白都不足以完成机制升级，更合理的方案应将 SLC7A11/SLC3A2、NQO1、NADPH/GSH、HMOX1/IL-1β/TNF-α 以及 CD31/GFAP/IBA1 共定位纳入同一验证框架。",
        ]
    )
    (REPORT_DIR / "literature_driven_extension_insertion_text_20260613.md").write_text(insertion, encoding="utf-8")


def main() -> None:
    configure_style()
    adata = ad.read_h5ad(H5AD_PATH)
    gene_df = raw_gene_df(adata, ALL_PLOT_GENES)

    ext_df = summarize_gene_localization(adata, gene_df, EXTENSION_GENES)
    inflam_df = summarize_gene_localization(adata, gene_df, INFLAMMATORY_GENES)
    delta_df = donor_level_deltas(adata, gene_df, EXTENSION_GENES + INFLAMMATORY_GENES + CORE_GENES)
    celltype_condition_df, corr_df = inflammation_module_summary(adata, gene_df)

    ext_df.to_csv(TABLE_DIR / "GSE209552_literature_extension_celltype_localization_20260613.csv", index=False, encoding="utf-8-sig")
    inflam_df.to_csv(TABLE_DIR / "GSE209552_literature_inflammation_celltype_localization_20260613.csv", index=False, encoding="utf-8-sig")
    delta_df.to_csv(TABLE_DIR / "GSE209552_literature_extension_donor_deltas_20260613.csv", index=False, encoding="utf-8-sig")
    celltype_condition_df.to_csv(TABLE_DIR / "GSE209552_literature_inflammation_module_summary_20260613.csv", index=False, encoding="utf-8-sig")
    corr_df.to_csv(TABLE_DIR / "GSE209552_literature_slc7a11_inflammation_correlations_20260613.csv", index=False, encoding="utf-8-sig")

    make_figure(ext_df, inflam_df, celltype_condition_df, delta_df)
    write_report(ext_df, inflam_df, celltype_condition_df, corr_df, delta_df)
    print("Literature-driven extension assets generated.")


if __name__ == "__main__":
    main()
