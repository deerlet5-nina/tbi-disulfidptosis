from __future__ import annotations

import importlib.util
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


DATE = "2026-06-05"
PRIORITY_GENES = ["SLC3A2", "SLC7A11", "WASF2", "TLN1", "ACTB", "MYH9", "MYL6", "FLNA"]
TRANSPORTER_GENES = ["SLC3A2", "SLC7A11"]
ACTIN_GENES = ["WASF2", "TLN1", "ACTB", "MYH9", "MYL6", "FLNA"]

OKABE = {
    "blue": "#0072B2",
    "orange": "#E69F00",
    "green": "#009E73",
    "red": "#D55E00",
    "purple": "#CC79A7",
    "sky": "#56B4E9",
    "yellow": "#F0E442",
    "black": "#111827",
    "gray": "#6B7280",
    "lightgray": "#E5E7EB",
}


ROOT = Path.cwd()
V3_PATH = next(ROOT.rglob("build_disulfidptosis_manuscript_v3_20260604.py"))
spec = importlib.util.spec_from_file_location("disulfidptosis_v3", V3_PATH)
v3 = importlib.util.module_from_spec(spec)
assert spec.loader is not None
spec.loader.exec_module(v3)

WORKBOOK = [
    p
    for p in ROOT.rglob("TBI_disulfidptosis_focused_design_tables_20260604.xlsx")
    if "Phase3_深化优化与最终报告_20260506-0513" in str(p)
    and "11_双硫死亡聚焦论文设计_20260604" in str(p)
][0]
WORKDIR = WORKBOOK.parents[1]
TABLEDIR = WORKDIR / "tables"
FIGDIR = WORKDIR / "figures"
REPORTDIR = WORKDIR / "reports"
for d in [TABLEDIR, FIGDIR, REPORTDIR]:
    d.mkdir(parents=True, exist_ok=True)

v3.WORKDIR = WORKDIR
v3.TABLEDIR = TABLEDIR
v3.FIGDIR = FIGDIR
v3.REPORTDIR = REPORTDIR


def configure_style() -> None:
    v3.configure_style()
    mpl.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Microsoft YaHei", "SimHei", "Arial", "DejaVu Sans"],
            "axes.unicode_minus": False,
            "font.size": 8.5,
            "axes.labelsize": 8.5,
            "axes.titlesize": 9.5,
            "xtick.labelsize": 7.5,
            "ytick.labelsize": 7.5,
            "legend.fontsize": 7.5,
            "figure.dpi": 160,
            "savefig.dpi": 320,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )


def fmt(x: object) -> str:
    if pd.isna(x):
        return ""
    if isinstance(x, (float, np.floating)):
        ax = abs(float(x))
        if ax != 0 and ax < 0.001:
            return f"{float(x):.2e}"
        if ax < 10:
            return f"{float(x):.3g}"
        return f"{float(x):.3f}"
    return str(x)


def markdown_table(df: pd.DataFrame, cols: list[str], max_rows: int | None = None) -> str:
    work = df.copy()
    if max_rows is not None:
        work = work.head(max_rows)
    lines = ["| " + " | ".join(cols) + " |", "| " + " | ".join(["---"] * len(cols)) + " |"]
    for _, row in work[cols].iterrows():
        lines.append("| " + " | ".join(fmt(row[c]).replace("\n", " ") for c in cols) + " |")
    return "\n".join(lines)


def figure_block(n: int, filename: str, caption: str) -> str:
    return f'<a id="fig-{n}"></a>\n\n**Fig. {n}. {caption}**\n\n![Fig. {n}](../figures/{filename}.png)'


def savefig(fig: plt.Figure, name: str) -> tuple[Path, Path]:
    png = FIGDIR / f"{name}.png"
    pdf = FIGDIR / f"{name}.pdf"
    fig.savefig(png, bbox_inches="tight", dpi=320)
    fig.savefig(pdf, bbox_inches="tight")
    plt.close(fig)
    return png, pdf


def load_inputs() -> dict[str, pd.DataFrame]:
    return {
        "evidence": pd.read_csv(TABLEDIR / "time_severity_extension_matrix_20260604.csv"),
        "human": pd.read_csv(TABLEDIR / "v3_human_brain_8gene_focused_results_20260604.csv"),
        "human_summary": pd.read_csv(TABLEDIR / "v3_human_brain_8gene_evidence_summary_20260604.csv"),
        "gse209": pd.read_csv(TABLEDIR / "v3_GSE209552_8gene_acute_severe_results_20260604.csv"),
        "marker_scores": pd.read_csv(TABLEDIR / "v3_GSE209552_bulk_marker_proxy_scores_20260604.csv"),
        "marker_corr": pd.read_csv(TABLEDIR / "v3_GSE209552_bulk_marker_proxy_correlations_20260604.csv"),
        "mouse": pd.read_csv(TABLEDIR / "v3_GSE163415_8gene_focused_DE_results_20260604.csv"),
        "sev_gene": pd.read_csv(TABLEDIR / "v3_GSE223245_8gene_severity_results_20260604.csv"),
        "sev_module": pd.read_csv(TABLEDIR / "v3_GSE223245_module_severity_results_20260604.csv"),
        "sev_scores": pd.read_csv(TABLEDIR / "v3_GSE223245_module_scores_long_20260604.csv"),
        "priority": pd.read_csv(TABLEDIR / "v3_integrated_8gene_priority_20260604.csv"),
        "wet": pd.read_csv(TABLEDIR / "zhang_mingyang_wet_validation_20260604.csv"),
    }


def build_extended_tables(data: dict[str, pd.DataFrame]) -> dict[str, pd.DataFrame]:
    human = data["human"].copy()
    human_comp = (
        human.groupby(["dataset", "comparison_short"])
        .agg(
            mean_effect=("effect", "mean"),
            positive_genes=("effect", lambda s: int((s > 0).sum())),
            nominal_genes=("significant_nominal_0_05", "sum"),
            fdr_genes=("significant_FDR_0_05", "sum"),
            n_genes=("gene_symbol", "nunique"),
        )
        .reset_index()
        .sort_values(["fdr_genes", "nominal_genes", "mean_effect"], ascending=[False, False, False])
    )

    mouse = data["mouse"].copy()
    mouse_units = (
        mouse.groupby(["analysis_unit", "time", "scope", "region", "treatment"])
        .agg(
            mean_logFC=("logFC_TBI_vs_NoTBI", "mean"),
            median_logFC=("logFC_TBI_vs_NoTBI", "median"),
            positive_genes=("logFC_TBI_vs_NoTBI", lambda s: int((s > 0).sum())),
            nominal_genes=("significant_nominal_0_05", "sum"),
            panel_FDR_genes=("significant_FDR_disulfidptosis_0_05", "sum"),
            min_p=("p_value", "min"),
            n_genes=("gene_symbol", "nunique"),
        )
        .reset_index()
        .sort_values(["panel_FDR_genes", "nominal_genes", "mean_logFC"], ascending=[False, False, False])
    )

    paired = mouse[mouse["treatment"].isin(["Drug", "Veh"])].pivot_table(
        index=["time", "scope", "region", "gene_symbol", "module"],
        columns="treatment",
        values="logFC_TBI_vs_NoTBI",
        aggfunc="mean",
    )
    paired = paired.reset_index()
    paired["drug_minus_vehicle_logFC"] = paired["Drug"] - paired["Veh"]
    treatment_module = (
        paired.dropna(subset=["drug_minus_vehicle_logFC"])
        .groupby(["time", "scope", "region", "module"])
        .agg(
            mean_drug_minus_vehicle=("drug_minus_vehicle_logFC", "mean"),
            median_drug_minus_vehicle=("drug_minus_vehicle_logFC", "median"),
            genes=("gene_symbol", "nunique"),
        )
        .reset_index()
        .sort_values(["mean_drug_minus_vehicle"], ascending=False)
    )

    corr = data["marker_corr"].copy()
    corr["abs_r"] = corr["spearman_r"].abs()
    marker_priority = corr.sort_values("abs_r", ascending=False).copy()
    marker_priority["interpretation"] = np.where(
        marker_priority["FDR"] < 0.05,
        "FDR-supported proxy clue",
        np.where(marker_priority["p_value"] < 0.05, "nominal proxy clue", "directional proxy clue"),
    )

    validation_matrix = pd.DataFrame(
        [
            {
                "evidence_gap": "Transporter entry",
                "readout": "SLC3A2/SLC7A11 mRNA and protein",
                "best_current_support": "CTE stage FDR support plus 3DPI hippocampus mouse support",
                "required_upgrade": "qPCR/WB/IHC in same region and time window",
                "priority": "High",
            },
            {
                "evidence_gap": "Actin cytoskeletal endpoint",
                "readout": "WASF2/TLN1 plus ACTB/MYH9/MYL6/FLNA and F-actin morphology",
                "best_current_support": "Acute severe TBI actin module and 3DPI hippocampus CCI",
                "required_upgrade": "Phalloidin staining, non-reducing gel, cytoskeletal fraction WB",
                "priority": "High",
            },
            {
                "evidence_gap": "Redox prerequisite",
                "readout": "NADPH/NADP+, GSH/GSSG, cystine/cysteine",
                "best_current_support": "Mechanistic plausibility from disulfidptosis literature",
                "required_upgrade": "Biochemical assays matched to transcriptomic window",
                "priority": "High",
            },
            {
                "evidence_gap": "Cell-type localization",
                "readout": "NeuN, GFAP, IBA1, OLIG2 and CD31 co-localization",
                "best_current_support": "Bulk marker-proxy suggests endothelial and neuronal prioritization",
                "required_upgrade": "IF/IHC or annotated snRNA-seq extraction",
                "priority": "High",
            },
            {
                "evidence_gap": "Severity relation",
                "readout": "Brain tissue severity strata and peripheral blood replication",
                "best_current_support": "GSE223245 blood trend is exploratory and sample-limited",
                "required_upgrade": "Independent severe/moderate/mild cohort or animal graded-injury design",
                "priority": "Medium",
            },
        ]
    )

    outputs = {
        "human_comp": human_comp,
        "mouse_units": mouse_units,
        "treatment_module": treatment_module,
        "marker_priority": marker_priority,
        "validation_matrix": validation_matrix,
    }
    for name, df in outputs.items():
        df.to_csv(TABLEDIR / f"v4_{name}_20260605.csv", index=False, encoding="utf-8-sig")
    return outputs


def draw_graphical_abstract(data: dict[str, pd.DataFrame]) -> None:
    fig, ax = plt.subplots(figsize=(12.0, 5.8))
    ax.axis("off")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    steps = [
        ("Public TBI/CTE\ntranscriptomes", 0.12, 0.58, OKABE["blue"]),
        ("Fixed 8-gene\ntransporter-actin panel", 0.37, 0.58, OKABE["green"]),
        ("Time, region,\nseverity and cell proxy", 0.62, 0.58, OKABE["orange"]),
        ("Validation window:\n3DPI hippocampus/cortex", 0.87, 0.58, OKABE["red"]),
    ]
    for text, x, y, color in steps:
        rect = mpl.patches.FancyBboxPatch(
            (x - 0.09, y - 0.12),
            0.18,
            0.24,
            boxstyle="round,pad=0.02,rounding_size=0.025",
            fc=color,
            ec=OKABE["black"],
            lw=1.0,
            alpha=0.14,
        )
        ax.add_patch(rect)
        ax.text(x, y, text, ha="center", va="center", fontsize=11, fontweight="bold")
    centers = [step[1] for step in steps]
    for left, right in zip(centers[:-1], centers[1:]):
        ax.annotate("", xy=(right - 0.11, 0.58), xytext=(left + 0.11, 0.58), arrowprops=dict(arrowstyle="->", lw=1.6, color=OKABE["gray"]))

    ax.text(0.5, 0.88, "Disulfidptosis-like transporter-actin cytoskeletal stress after TBI", ha="center", fontsize=15, fontweight="bold")
    ax.text(0.5, 0.80, "The manuscript asks when, where, at what injury severity, and in which cellular context the signal is most plausible.", ha="center", fontsize=10)

    findings = [
        ("Human chronic CTE", "FDR-supported\nstage-associated axis"),
        ("Acute severe TBI", "directional actin\nmodule increase"),
        ("Mouse CCI", "3DPI hippocampus\nvalidation window"),
        ("Peripheral blood", "severity clue only;\nnot brain mechanism"),
    ]
    y0 = 0.28
    for i, (label, text) in enumerate(findings):
        x = 0.11 + i * 0.255
        ax.plot([x - 0.075, x + 0.075], [y0 + 0.06, y0 + 0.06], color=list(OKABE.values())[i], lw=3)
        ax.text(x, y0, label, ha="center", va="center", fontsize=10, fontweight="bold")
        ax.text(x, y0 - 0.078, text, ha="center", va="center", fontsize=8.2, linespacing=1.2)
    ax.text(0.5, 0.055, "Mechanistic proof still requires protein, F-actin, redox and co-localization readouts in the same tissue window.", ha="center", fontsize=9.5, color=OKABE["black"])
    savefig(fig, "GraphicalAbstract_v4_TBI_disulfidptosis_20260605")


def draw_fig7(data: dict[str, pd.DataFrame], ext: dict[str, pd.DataFrame]) -> None:
    fig = plt.figure(figsize=(11.8, 8.0))
    gs = fig.add_gridspec(2, 2, hspace=0.48, wspace=0.55)
    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[0, 1])
    ax_c = fig.add_subplot(gs[1, 0])
    ax_d = fig.add_subplot(gs[1, 1])

    human = ext["human_comp"].copy()
    human["label"] = human["comparison_short"].str.replace("GSE209552 human acute brain severe", "209552 acute", regex=False)
    human = human.sort_values("mean_effect", ascending=False)
    colors = [OKABE["red"] if v > 0 else OKABE["gray"] for v in human["mean_effect"]]
    ax_a.barh(np.arange(len(human)), human["mean_effect"], color=colors, alpha=0.85)
    ax_a.set_yticks(np.arange(len(human)))
    ax_a.set_yticklabels(human["label"], fontsize=7)
    ax_a.invert_yaxis()
    ax_a.set_xlabel("Mean effect across 8 genes")
    ax_a.set_title("Human-brain comparison-level directionality")
    for i, r in human.reset_index(drop=True).iterrows():
        ax_a.text(r["mean_effect"] + 0.02, i, f"FDR {int(r['fdr_genes'])}; nom {int(r['nominal_genes'])}", va="center", fontsize=6.5)
    v3.add_panel_label(ax_a, "A")

    mouse = ext["mouse_units"].head(12).iloc[::-1]
    y = np.arange(len(mouse))
    ax_b.barh(y, mouse["nominal_genes"], color=OKABE["sky"], label="Nominal genes")
    ax_b.barh(y, mouse["panel_FDR_genes"], color=OKABE["blue"], label="Panel-FDR genes")
    labels = [f"{r.time} {r.region} {r.treatment}" for r in mouse.itertuples()]
    ax_b.set_yticks(y)
    ax_b.set_yticklabels(labels, fontsize=7)
    ax_b.set_xlabel("Number of 8-gene panel hits")
    ax_b.set_title("Mouse CCI validation-window ranking")
    ax_b.legend(frameon=False, loc="lower right")
    v3.add_panel_label(ax_b, "B")

    tm = ext["treatment_module"].copy()
    tm = tm[tm["scope"].isin(["region_treatment", "treatment"])]
    tm["unit"] = tm["time"] + " " + tm["region"]
    heat = tm.pivot_table(index="unit", columns="module", values="mean_drug_minus_vehicle", aggfunc="mean")
    sns.heatmap(heat, cmap="RdBu_r", center=0, annot=True, fmt=".2f", linewidths=0.4, cbar_kws={"label": "Drug - vehicle logFC"}, ax=ax_c)
    ax_c.set_title("Treatment-stratified modulation in GSE163415")
    ax_c.set_xlabel("")
    ax_c.set_ylabel("")
    v3.add_panel_label(ax_c, "C")

    val = ext["validation_matrix"].copy()
    priority_map = {"High": 3, "Medium": 2, "Low": 1}
    val["priority_score"] = val["priority"].map(priority_map)
    ax_d.barh(np.arange(len(val)), val["priority_score"], color=[OKABE["red"], OKABE["red"], OKABE["orange"], OKABE["orange"], OKABE["green"]])
    ax_d.set_xlim(0, 3.35)
    ax_d.set_xticks([1, 2, 3])
    ax_d.set_xticklabels(["Low", "Medium", "High"])
    ax_d.set_yticks(np.arange(len(val)))
    ax_d.set_yticklabels(val["evidence_gap"], fontsize=8)
    ax_d.invert_yaxis()
    ax_d.set_title("Mechanistic upgrade priorities")
    for i, r in val.reset_index(drop=True).iterrows():
        ax_d.text(3.08, i, r["priority"], va="center", fontsize=7, fontweight="bold")
    v3.add_panel_label(ax_d, "D")

    sns.despine(fig=fig)
    savefig(fig, "Fig7_v4_evidence_audit_and_validation_priorities_20260605")


def draw_fig8(data: dict[str, pd.DataFrame], ext: dict[str, pd.DataFrame]) -> None:
    fig = plt.figure(figsize=(12.0, 8.2))
    gs = fig.add_gridspec(2, 2, hspace=0.48, wspace=0.45)
    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[0, 1])
    ax_c = fig.add_subplot(gs[1, 0])
    ax_d = fig.add_subplot(gs[1, 1])

    ax_a.axis("off")
    tiers = [
        ("Mechanistic proof", "protein disulfidation; F-actin collapse;\nredox depletion; cell death phenotype", 0.80, OKABE["red"]),
        ("Tissue validation", "qPCR/WB/IF in matched\n3DPI brain region", 0.62, OKABE["orange"]),
        ("Public brain transcriptome", "human CTE stage; acute severe TBI;\nmouse CCI support", 0.44, OKABE["blue"]),
        ("Peripheral severity clue", "blood/PBMC trend only;\ncontext, not brain mechanism", 0.26, OKABE["green"]),
    ]
    for label, text, y, color in tiers:
        ax_a.add_patch(mpl.patches.Rectangle((0.05, y - 0.075), 0.88, 0.135, fc=color, alpha=0.14, ec=color, lw=1.0))
        ax_a.text(0.08, y + 0.048, label, va="top", ha="left", fontsize=8.0, fontweight="bold")
        ax_a.text(0.08, y + 0.014, text, va="top", ha="left", fontsize=6.5, linespacing=1.03)
    ax_a.set_title("Claim ladder used in the revised manuscript")
    v3.add_panel_label(ax_a, "A")

    corr = data["marker_corr"].copy()
    heat = corr.pivot_table(index="module", columns="celltype_proxy", values="spearman_r")
    sns.heatmap(heat, cmap="RdBu_r", center=0, vmin=-1, vmax=1, annot=True, fmt=".2f", linewidths=0.4, cbar_kws={"label": "Spearman r"}, ax=ax_b)
    ax_b.set_title("Bulk marker-proxy correlation map")
    ax_b.set_xlabel("")
    ax_b.set_ylabel("")
    ax_b.tick_params(axis="x", rotation=35)
    v3.add_panel_label(ax_b, "B")

    timepoints = ["1h", "6h", "12h", "1D", "2D", "3D", "7D", "29D"]
    rows = ["transporter", "actin/F-actin", "redox", "co-localization"]
    mat = pd.DataFrame(0.15, index=rows, columns=timepoints)
    for tp in ["3D", "7D"]:
        mat.loc[:, tp] = [0.95, 0.95, 0.90, 0.85]
    mat.loc[:, "29D"] = [0.55, 0.35, 0.35, 0.40]
    mat.loc[:, ["1h", "6h", "12h", "1D", "2D"]] = [0.35, 0.55, 0.60, 0.70, 0.80]
    sns.heatmap(mat, cmap="YlGnBu", vmin=0, vmax=1, annot=False, linewidths=0.4, cbar_kws={"label": "Suggested validation priority"}, ax=ax_c)
    ax_c.set_title("Wet-lab validation time-course design")
    ax_c.set_xlabel("Post-injury time")
    ax_c.set_ylabel("Readout class")
    v3.add_panel_label(ax_c, "C")

    sev = data["sev_gene"].copy().sort_values("spearman_r_severity_all_groups")
    colors = [OKABE["blue"] if p < 0.05 else OKABE["gray"] for p in sev["p_value_severity_all_groups"]]
    ax_d.barh(sev["gene_symbol"], sev["spearman_r_severity_all_groups"], color=colors, alpha=0.85)
    ax_d.axvline(0, color="black", lw=0.8)
    ax_d.set_xlim(-0.9, 0.55)
    ax_d.set_xlabel("Spearman r across severity 0-3")
    ax_d.set_title("Peripheral severity gene trends")
    for i, r in enumerate(sev.itertuples()):
        ax_d.text(0.51, i, f"P={r.p_value_severity_all_groups:.2g}", ha="right", va="center", fontsize=6.5)
    v3.add_panel_label(ax_d, "D")

    sns.despine(fig=fig)
    savefig(fig, "Fig8_v4_claim_ladder_cell_proxy_and_wetlab_design_20260605")


def build_key_values(data: dict[str, pd.DataFrame], ext: dict[str, pd.DataFrame]) -> dict[str, str]:
    human = data["human"]
    stage = human[human["comparison"].eq("GSE193407_human_prefrontal_BA9_CTE_stage_trend")]
    stage_fdr = stage[stage["significant_FDR_0_05"] == True].sort_values("FDR")
    late = human[human["comparison"].eq("GSE193407_human_prefrontal_BA9_late_CTE_stage3_4_vs_stage0")]
    late_fdr = late[late["significant_FDR_0_05"] == True].sort_values("FDR")
    g209 = data["gse209"].set_index("gene_symbol")
    mouse_top = ext["mouse_units"].iloc[0]
    mouse_hipp = ext["mouse_units"][ext["mouse_units"]["analysis_unit"].eq("GSE163415_3DPI_Hipp_all_treatments")].iloc[0]
    sev_mod = data["sev_module"].sort_values("p_value_severity_all_groups").iloc[0]
    sev_gene = data["sev_gene"].sort_values("p_value_severity_all_groups").iloc[0]
    marker_top = ext["marker_priority"].iloc[0]
    treat_top = ext["treatment_module"].iloc[0]
    priority_top = data["priority"].iloc[0]
    return {
        "stage_details": "; ".join(f"{r.gene_symbol}: r={r.stage_correlation:.3f}, FDR={r.FDR:.3g}" for r in stage_fdr.itertuples()),
        "late_details": "; ".join(f"{r.gene_symbol}: logFC={r.logFC:.3f}, FDR={r.FDR:.3g}" for r in late_fdr.itertuples()),
        "acute_genes": "; ".join(
            f"{gene} logFC={g209.loc[gene, 'GSE209552_logFC_TBI_vs_Control']:.3f}, P={g209.loc[gene, 'GSE209552_p_value_ttest']:.3g}"
            for gene in ["MYH9", "MYL6", "ACTB", "TLN1"]
        ),
        "mouse_top": f"{mouse_top.time}/{mouse_top.region}/{mouse_top.treatment}: mean logFC={mouse_top.mean_logFC:.3f}, positive genes={int(mouse_top.positive_genes)}/8, nominal={int(mouse_top.nominal_genes)}, panel-FDR={int(mouse_top.panel_FDR_genes)}",
        "mouse_hipp": f"3DPI hippocampus all treatments: mean logFC={mouse_hipp.mean_logFC:.3f}, positive genes={int(mouse_hipp.positive_genes)}/8, nominal={int(mouse_hipp.nominal_genes)}, panel-FDR={int(mouse_hipp.panel_FDR_genes)}",
        "severity_module": f"{sev_mod.module}: r={sev_mod.spearman_r_severity_all_groups:.3f}, P={sev_mod.p_value_severity_all_groups:.3g}, FDR={sev_mod.FDR_severity_all_groups:.3g}",
        "severity_gene": f"{sev_gene.gene_symbol}: r={sev_gene.spearman_r_severity_all_groups:.3f}, P={sev_gene.p_value_severity_all_groups:.3g}, FDR={sev_gene.FDR_severity_all_groups:.3g}",
        "marker_top": f"{marker_top.module} vs {marker_top.celltype_proxy}: r={marker_top.spearman_r:.3f}, P={marker_top.p_value:.3g}, FDR={marker_top.FDR:.3g}",
        "treatment_top": f"{treat_top.time}/{treat_top.region}/{treat_top.module}: mean Drug-Veh logFC difference={treat_top.mean_drug_minus_vehicle:.3f}",
        "priority_top": f"{priority_top.gene_symbol}, integrated evidence score={priority_top.integrated_evidence_score:.1f}",
    }


def build_table_strings(data: dict[str, pd.DataFrame], ext: dict[str, pd.DataFrame]) -> dict[str, str]:
    dataset_cols = ["dataset_or_source", "species", "tissue_or_model", "time_or_course", "severity_or_stage", "current_role_in_manuscript", "guardrail"]
    human_cols = ["dataset", "comparison_short", "mean_effect", "positive_genes", "nominal_genes", "fdr_genes", "n_genes"]
    mouse_cols = ["time", "region", "treatment", "mean_logFC", "positive_genes", "nominal_genes", "panel_FDR_genes", "min_p"]
    priority_cols = ["gene_symbol", "human_FDR_comparisons", "human_nominal_comparisons", "GSE209552_acute_logFC", "mouse_nominal_units", "GSE223245_severity_r_all", "integrated_evidence_score"]
    severity_cols = ["gene_symbol", "spearman_r_severity_all_groups", "p_value_severity_all_groups", "FDR_severity_all_groups", "spearman_r_TBI_only", "p_value_TBI_only"]
    marker_cols = ["module", "celltype_proxy", "spearman_r", "p_value", "FDR", "interpretation"]
    validation_cols = ["evidence_gap", "readout", "best_current_support", "required_upgrade", "priority"]
    return {
        "datasets": markdown_table(data["evidence"], dataset_cols),
        "human": markdown_table(ext["human_comp"], human_cols),
        "mouse": markdown_table(ext["mouse_units"], mouse_cols, max_rows=12),
        "priority": markdown_table(data["priority"], priority_cols),
        "severity": markdown_table(data["sev_gene"].sort_values("p_value_severity_all_groups"), severity_cols),
        "marker": markdown_table(ext["marker_priority"], marker_cols, max_rows=12),
        "validation": markdown_table(ext["validation_matrix"], validation_cols),
    }


def zh_manuscript(vals: dict[str, str], tables: dict[str, str]) -> str:
    ga = "GraphicalAbstract_v4_TBI_disulfidptosis_20260605"
    return f"""# TBI 后双硫死亡样转运-肌动蛋白骨架应激的时空定位、损伤程度关联与细胞类型优先级：基于人类 TBI/CTE 和小鼠 CCI 公共转录组的整合生物信息学研究

## 摘要

背景：创伤性脑损伤（traumatic brain injury, TBI）后的继发性损伤具有显著的时间、脑区、损伤程度和细胞背景异质性。双硫死亡是一种由 SLC7A11/SLC3A2 介导的胱氨酸输入、NADPH 消耗和肌动蛋白骨架二硫键压力驱动的新型细胞死亡形式，但公共疾病转录组尚不能直接证明其在 TBI 后发生。为了避免把广义氧化应激或细胞死亡通路误写为双硫死亡，本研究预先固定 8 个基因，聚焦“胱氨酸转运入口-肌动蛋白骨架终点”这一可检验轴线。

方法：本研究整合 GSE209552、GSE193407、GSE319253、GSE104687、GSE163415 和 GSE223245。SLC3A2 与 SLC7A11 定义转运入口，WASF2、TLN1、ACTB、MYH9、MYL6 和 FLNA 定义肌动蛋白骨架调控和结构终点。分析按证据层分开解释：GSE209552 用于急性 severe TBI 人脑方向性证据，GSE193407 和 GSE319253 用于慢性 CTE 病程证据，GSE104687 用于远期多脑区探索，小鼠 GSE163415 用于 CCI 时间-脑区旁证，GSE223245 仅用于外周血损伤程度线索。统计报告区分 FDR 支持和 nominal 探索性信号；同一 donor 多脑区样本不作为独立主分析；动物和外周血结果不被提升为人脑因果结论。

结果：人脑慢性病程证据最稳健。GSE193407 中 CTE stage trend 的 FDR 支持结果为 {vals['stage_details']}；late CTE stage 3-4 vs stage 0 的 FDR 支持结果为 {vals['late_details']}。急性 severe TBI 脑组织 GSE209552 显示所有 8 个基因方向均为正，较强的单基因方向包括 {vals['acute_genes']}，但样本量小，不能写作诊断标志物。小鼠 CCI 中，验证窗口排序显示 {vals['mouse_top']}，其中 {vals['mouse_hipp']}，提示 3DPI hippocampus 是当前最集中的动物验证窗口。外周血 GSE223245 提供损伤程度补充线索，模块层面为 {vals['severity_module']}，单基因层面最强为 {vals['severity_gene']}。GSE209552 bulk marker-proxy 中最强线索为 {vals['marker_top']}，只能指导后续共定位优先级。综合证据排序首位为 {vals['priority_top']}，但该排序用于实验优先级，不用于临床预测。

结论：现有公共转录组支持 TBI/CTE 后存在“候选双硫死亡样转运-肌动蛋白骨架应激”，其中慢性 CTE 的 SLC7A11/SLC3A2/WASF2/TLN1 轴、急性 severe TBI 的骨架方向性改变，以及小鼠 CCI 的 3DPI hippocampus 窗口构成当前最一致的证据链。下一步机制升级必须在同一时间窗和脑区中同时检测 SLC3A2/SLC7A11 蛋白、WASF2/TLN1 和 ACTB/MYH9/MYL6/FLNA 终点、F-actin 形态、NADPH/GSH 氧化还原读数，以及 NeuN/GFAP/IBA1/OLIG2/CD31 细胞类型共定位。

关键词：创伤性脑损伤；慢性创伤性脑病；双硫死亡；SLC7A11；SLC3A2；肌动蛋白骨架；损伤时间；损伤严重程度；公共转录组

## 图形摘要

<a id="fig-0"></a>

**Fig. 0. Graphical abstract of the revised manuscript.**

![Fig. 0](../figures/{ga}.png)

## 引言

TBI 的病理生理过程并不是一次性完成的事件。机械性原发损伤之后，局部脑组织会迅速进入能量代谢危机、离子稳态破坏、兴奋性毒性、线粒体损伤和氧化还原压力状态；随后数天至数周内，血脑屏障破坏、胶质细胞反应、细胞骨架重塑、轴突和突触结构改变以及炎症持续化共同塑造继发性损伤。慢性或反复轻度损伤背景下，CTE 又表现为长期神经炎症、细胞组成变化和神经退行性病理。正因为这些过程具有明显的时空分层，单纯比较 TBI 与 control 容易把不同阶段的信号混在一起。一个更适合机制推进的论文设计，应当把问题拆成“何时激活、何处更明显、与何种损伤程度相关、在哪类细胞背景中更值得验证”。

双硫死亡为这一问题提供了新的机制入口。原始机制研究指出，在 SLC7A11 高表达且葡萄糖不足的条件下，细胞持续摄入胱氨酸会消耗 NADPH，还原力不足使二硫键压力积累，最终导致肌动蛋白骨架蛋白异常二硫键化和 F-actin 骨架塌陷。TBI 后脑组织存在能量危机、线粒体功能障碍和氧化还原失衡，因此从病理生理背景看，双硫死亡样应激具有合理性。然而，公共转录组只能告诉我们 mRNA 或模块层面的变化，不能直接证明蛋白二硫键化、F-actin 塌陷或某一种细胞死亡形式。因此，本研究采用更谨慎的术语：“双硫死亡样转运-肌动蛋白骨架应激”。

本研究没有把所有含硫代谢、氧化应激和细胞死亡基因都纳入主线，而是将问题收束到 8 个预先固定的关键基因。SLC3A2 和 SLC7A11 代表胱氨酸转运入口；WASF2 和 TLN1 代表 WRC 介导的肌动蛋白调控、黏附和张力传递；ACTB、MYH9、MYL6 与 FLNA 代表更接近终点的骨架结构和收缩/张力节点。这种设计牺牲了通路覆盖面，但换来了清晰的可检验性：如果 TBI/CTE 确实存在与双硫死亡机制相邻的应激过程，应当能在特定时间窗、脑区或病程阶段看到转运入口与骨架终点的协同或阶段性改变。

本研究的目标是基于现有公共数据形成一篇完整、可继续投向湿实验验证的生物信息学论文。文章主线并不是证明 TBI 后已经发生双硫死亡，而是建立一个严格分层的证据模型：人脑慢性病程证据用于回答 CTE stage 关联，急性 severe TBI 脑组织用于回答早期重型损伤方向，小鼠 CCI 用于定位可验证时间窗和脑区，外周血严重程度数据用于补充 severity-aware 叙事，bulk marker-proxy 用于提出细胞类型验证优先级。这样的写法既保留了论文的机制想象力，也避免超出公共转录组能够支持的证据边界。

## 材料与方法

### 数据集纳入与证据分层

本研究纳入六个 GEO 公共转录组数据层。GSE209552 代表人类急性 severe TBI 外科切除脑组织，样本来自损伤后约 4 小时至 8 天，适合提供急性重型损伤方向性证据。GSE193407 代表人类 BA9 脑区 CTE stage 0-4，适合分析慢性病程梯度。GSE319253 作为 superior frontal cortex CTE vs control 的外部慢性验证层。GSE104687 是远期 TBI 多脑区死后脑组织数据，由于同一 donor 的多个脑区样本不应视为独立观测，因此只作为空间探索层。GSE163415 是小鼠 controlled cortical impact（CCI）数据，覆盖 3DPI/29DPI、hippocampus/thalamus/hypothalamus 和 vehicle/drug 分层，用于跨物种时间-脑区旁证。GSE223245 是人类 whole blood/PBMC mild、moderate 和 severe TBI 微阵列数据，仅作为外周严重程度线索。

{tables['datasets']}

### 候选基因、模块与统计原则

候选基因在分析前固定为 SLC3A2、SLC7A11、WASF2、TLN1、ACTB、MYH9、MYL6 和 FLNA。模块层面定义三个读数：8 基因总模块、SLC3A2/SLC7A11 转运入口模块，以及 WASF2/TLN1/ACTB/MYH9/MYL6/FLNA 肌动蛋白骨架模块。样本级模块评分采用每个基因在表达矩阵中按样本 z-score 标准化后求均值。单基因比较使用现有差异分析结果或 Welch t 检验；严重程度趋势使用 Spearman 相关；多重检验使用 Benjamini-Hochberg FDR。正文中凡 FDR<0.05 称为 FDR 支持，P<0.05 但 FDR 未达标称为 nominal 探索性信号。

### GSE223245 外周严重程度 focused 再分析

GSE223245 通过 series matrix 和平台 SOFT 注释重建表达矩阵。探针折叠到 gene symbol 时，保留同一基因平均表达最高的探针。样本按标题和 metadata 分为 Control、Mild、Moderate 和 Severe，并编码为 0、1、2 和 3。对 8 个基因分别计算全组 severity 0-3 Spearman 相关和 TBI-only severity 1-3 Spearman 相关，同时进行 Mild/Moderate/Severe vs Control 的 Welch t 检验。因为每组仅 n=4 且样本来自外周血，本分析只用于补充损伤程度相关外周线索，不用于推断脑内机制。

### GSE209552 bulk marker-proxy 细胞类型优先级

GSE209552 当前可用处理结果为 bulk-like 脑组织基因矩阵，尚不足以直接进行细胞类型归属。为了给后续 IF/IHC 共定位排序，本研究计算神经元、星形胶质细胞、小胶质细胞、少突胶质细胞、OPC 和内皮细胞 marker-proxy 分数，并与 8 基因模块、转运入口模块和骨架模块进行 Spearman 相关。该分析只能表示“在哪些 marker 背景下模块更相关”，不能表示某一细胞类型已经发生双硫死亡。

### 综合优先级与论文式证据审计

综合优先级评分由人脑 FDR 支持次数、人脑 nominal 支持次数、小鼠 CCI 支持单元、GSE209552 急性方向和 GSE223245 severity 线索构成。该评分只用于选择下一轮实验读数，而不是统计推断或临床预测模型。为提高论文可读性，v4 版新增 human comparison-level 审计、GSE163415 time-region-treatment 验证窗口排序、treatment 分层差异以及机制升级矩阵，并将其写入图表和结果段落。

## 结果

### 多证据层设计将问题限定为可检验的转运-骨架轴

图形摘要和 [Fig. 1](#fig-1) 展示了本研究的整体逻辑。研究没有把“双硫死亡”作为已经发生的事实，而是将其拆解为转运入口、骨架终点、氧化还原前提和细胞类型定位四类可验证读数。人类急性 severe TBI、慢性 CTE stage、远期多脑区 TBI、小鼠 CCI 和外周血 severe-aware 数据分别承担不同证据角色。这样的层级化设计避免了两个常见问题：其一，把外周血或动物结果直接写成人脑机制；其二，把 nominal 候选直接写成因果结论。

{figure_block(1, 'Fig1_v3_study_design_evidence_layers_20260604', 'Study design, evidence layers, pre-fixed 8-gene panel and interpretation guardrails.')}

### 人类脑组织显示慢性 CTE 病程中 SLC7A11/SLC3A2/WASF2/TLN1 轴最稳定

在人类脑组织中，8 基因面板的信号并不均一。[Fig. 2](#fig-2) 和 Table 2 显示，GSE193407 的 CTE stage trend 是当前最强的人脑证据层，达到 FDR 支持的基因为 {vals['stage_details']}；late stage 3-4 vs stage 0 中达到 FDR 支持的基因为 {vals['late_details']}。相比之下，GSE104687 多脑区远期 TBI 数据仅提供少量 nominal 空间线索，且受到同 donor 多脑区非独立性的限制。GSE319253 可作为慢性外部验证层，但当前更适合作为辅助证据而非主结论。

{figure_block(2, 'Fig2_v3_human_brain_multidataset_8gene_20260604', 'Human brain multidataset evidence for the 8-gene transporter-actin panel.')}

Table 2. Human comparison-level evidence audit.

{tables['human']}

### 急性 severe TBI 中骨架模块方向性增强，但仍处于探索性层级

GSE209552 代表急性 severe TBI 脑组织窗口。[Fig. 3](#fig-3) 显示，TBI 相对 control 的 8 个基因方向均为正，较强的单基因方向包括 {vals['acute_genes']}。这提示急性重型损伤脑组织中已经存在转运-骨架轴的方向性变化，尤其是骨架终点更明显。不过，TBI n=4、control n=3 的样本量决定了它不能承担“正式阳性发现”或“诊断标志物”的角色。更准确的写法是：急性 severe TBI 为后续实验提供方向性证据，其中骨架模块比转运模块更容易在早期窗口被观察到。

{figure_block(3, 'Fig3_v3_GSE209552_acute_severe_and_marker_proxy_20260604', 'Acute severe TBI evidence in GSE209552, including module scores, gene-level directionality and bulk marker-proxy prioritization.')}

### 小鼠 CCI 将优先验证窗口定位到 3DPI hippocampus

GSE163415 是本研究回答“何时、何处”的关键动物旁证。[Fig. 4](#fig-4) 和 [Fig. 7](#fig-7) 显示，time-region-treatment 单元排序中最强窗口为 {vals['mouse_top']}。按脑区和时间汇总时，3DPI hippocampus 的支持最集中，{vals['mouse_hipp']}。该结果说明，早期亚急性阶段比 29DPI 更适合作为第一轮湿实验验证窗口，hippocampus 比 thalamus 和 hypothalamus 更容易呈现转运-骨架模块协同改变。29DPI 仍有部分 SLC7A11 和转运/代谢轴残留，但不适合作为第一轮机制验证的唯一时间点。

{figure_block(4, 'Fig4_v3_GSE163415_mouse_CCI_spatiotemporal_suite_20260604', 'Mouse CCI spatiotemporal evidence in GSE163415.')}

Table 3. GSE163415 top time-region-treatment validation windows.

{tables['mouse']}

### Treatment 分层提示 3DPI 与 29DPI 的模块性质不同

GSE163415 的 treatment 分层不是本研究的主要疗效结论，但可用于理解模块的阶段性。[Fig. 7](#fig-7)C 显示 treatment-stratified Drug-Veh 差异在 29DPI core transporter/stress 模块更明显，当前最强分层差异为 {vals['treatment_top']}。这不能被写作药物机制结论，因为本研究没有围绕该治疗进行完整的因果设计；但它提示 29DPI 可能更偏向转运/代谢残留，而 3DPI hippocampus 更偏向转运入口与骨架终点共同增强。湿实验设计上，3DPI 应承担机制发现，7D 或 29DPI 可承担持续性观察。

{figure_block(7, 'Fig7_v4_evidence_audit_and_validation_priorities_20260605', 'Evidence audit, mouse validation-window ranking, treatment-stratified modulation and mechanistic upgrade priorities.')}

### 外周血 GSE223245 补充损伤程度叙事，但不能替代脑组织证据

为补足损伤程度维度，v3/v4 重新分析了 GSE223245。结果显示模块层面最强趋势为 {vals['severity_module']}，单基因层面最强为 {vals['severity_gene']}。[Fig. 5](#fig-5) 和 [Fig. 8](#fig-8)D 显示，SLC3A2 和 FLNA 在全组 severity 0-3 中呈负相关，MYL6 在 TBI-only 1-3 中呈正相关。由于样本来自 whole blood/PBMC，且每组仅 n=4，这些结果只能被写作 peripheral severity-aware clue。它们有助于提醒后续临床或动物设计纳入 mild/moderate/severe 或 graded-injury 维度，但不能被写作脑内双硫死亡机制证据。

{figure_block(5, 'Fig5_v3_GSE223245_peripheral_severity_suite_20260604', 'Peripheral severity-focused analysis of GSE223245 whole blood/PBMC data.')}

Table 4. GSE223245 8-gene severity trends.

{tables['severity']}

### Bulk marker-proxy 指向内皮和神经元优先共定位，但不是细胞归属

GSE209552 marker-proxy 分析显示，当前最强相关为 {vals['marker_top']}。此外，8 基因模块与 neuron proxy 呈负相关，actin 模块也与 endothelial proxy 呈较强正相关。[Fig. 8](#fig-8)B 将这些相关性作为热图展示。需要强调的是，bulk marker-proxy 受细胞比例、组织采样位置、坏死/血管成分和细胞状态共同影响，不能直接说明内皮细胞或神经元发生双硫死亡。它真正的价值是帮助决定后续共定位优先级：SLC3A2/SLC7A11 和 F-actin 读数应优先与 CD31、NeuN、GFAP、IBA1 和 OLIG2 同时检测，而不是只做单一细胞类型染色。

Table 5. GSE209552 top bulk marker-proxy prioritization clues.

{tables['marker']}

### 综合证据排序将表达读数和终点读数分开处理

综合排序显示，当前第一优先级为 {vals['priority_top']}。[Fig. 6](#fig-6) 和 Table 6 提示，SLC3A2/SLC7A11/WASF2/TLN1 更适合作为表达层面主读数，因为它们在人脑 CTE stage 中有更强 FDR 支持；ACTB/MYH9/MYL6/FLNA 更接近骨架终点，mRNA 改变应与蛋白、F-actin 形态和非还原检测联合解释。换言之，SLC3A2 和 SLC7A11 的意义在于“入口是否打开”，WASF2 和 TLN1 的意义在于“骨架调控和张力传递是否参与”，ACTB/MYH9/MYL6/FLNA 的意义在于“终点是否出现结构或机械异常”。

{figure_block(6, 'Fig6_v3_integrated_priority_validation_suite_20260604', 'Integrated 8-gene priority matrix, validation ranking, mechanistic upgrade workflow and marker-proxy summary.')}

Table 6. Integrated 8-gene evidence priority.

{tables['priority']}

### 机制升级需要同时满足转运、骨架、氧化还原和细胞定位四类证据

[Fig. 8](#fig-8) 将本文的结论强度分为四级：外周严重程度线索、公共脑组织转录组、组织层验证和机制证明。当前论文最强只能达到公共脑组织转录组与动物旁证层级，不能越级宣称双硫死亡已经发生。要把候选机制升级为机制支持，必须在同一时间窗和脑区中同时观察到：SLC3A2/SLC7A11 蛋白改变，WASF2/TLN1 与 ACTB/MYH9/MYL6/FLNA 终点异常，F-actin 形态改变，NADPH/NADP+ 与 GSH/GSSG 改变，以及细胞类型共定位。Table 7 给出建议的实验升级矩阵。

{figure_block(8, 'Fig8_v4_claim_ladder_cell_proxy_and_wetlab_design_20260605', 'Claim ladder, cell-proxy heatmap, wet-lab time-course design and peripheral severity gene trends.')}

Table 7. Mechanistic upgrade matrix for the next experimental stage.

{tables['validation']}

## 讨论

本研究的核心贡献不是发现某一个单独基因在 TBI 后上调，而是把一个新型细胞死亡机制转化为可检验的 TBI 转录组问题。双硫死亡的关键并不只是 SLC7A11 或 SLC3A2 表达变化，而是胱氨酸输入、还原力消耗、二硫键压力和肌动蛋白骨架崩塌之间的连续链条。公共转录组最多只能捕捉这条链的入口和部分终点，因此本文使用“双硫死亡样转运-肌动蛋白骨架应激”而不是“发生双硫死亡”。这种表述看起来更保守，但更适合真正的论文写作，因为它把可证明内容和待证明内容分开。

人类脑组织结果支持一个阶段性模型。在慢性 CTE 病程中，SLC7A11/SLC3A2/WASF2/TLN1 轴具有较稳定的 FDR 支持，说明反复或慢性损伤背景下，胱氨酸转运入口和 WRC/张力调控节点可能与病程推进相关。急性 severe TBI 中，骨架模块方向性增强更明显，而转运入口模块统计较弱，这可能反映急性期组织中骨架破坏、血管和胶质反应、细胞组成改变以及机械应力共同作用。两者并不矛盾：慢性病程更容易积累转运入口和调控轴信号，急性重型损伤更容易表现为骨架终点应激。

小鼠 CCI 分析为湿实验设计提供了最直接的信息。3DPI hippocampus 的 8 基因面板正向性和显著性最集中，因此第一轮实验不应平均铺开所有时间点和脑区。更经济、也更贴近证据链的设计是：先在 3DPI cortex/hippocampus 做 qPCR、WB、F-actin、NADPH/GSH 和细胞共定位，再用 7D 或 29DPI 判断信号是否持续或转向。若第一轮结果不支持转运入口与骨架终点共同改变，则不应继续用“双硫死亡”作为主机制，而应退回到更广义的骨架重塑或氧化还原应激。

损伤程度是另一个需要谨慎处理的维度。GSE209552 的 severe TBI、GSE193407 的 CTE stage 和 GSE223245 的 mild/moderate/severe 血液分组并不是同一种 severity。把它们合并为单一严重程度评分会制造统计幻觉。更合理的论文结构是分层叙述：GSE209552 表示急性重型脑组织窗口，GSE193407 表示慢性病程梯度，GSE223245 表示外周血严重程度线索。这样既回应了“损伤程度”的问题，也避免把不同生物学尺度混为一谈。

细胞类型解释仍是当前证据链中最大的缺口。内皮 marker-proxy 与 8 基因模块的相关性提示血管/屏障背景值得优先验证；神经元 proxy 的负相关提示神经元成分或神经元状态变化也不能忽视。但 bulk marker-proxy 无法区分细胞比例变化和同一细胞类型内表达状态变化。后续若能从 GSE209552 的 snRNA-seq 原始或处理对象中提取细胞注释，或者在动物组织中完成 SLC3A2/SLC7A11 与 NeuN、GFAP、IBA1、OLIG2、CD31 的共定位，论文的机制强度会显著提高。

本研究仍有明确局限。第一，各数据集的平台、物种、组织来源、疾病定义和协变量不同，不能进行简单合并。第二，GSE104687 存在同 donor 多脑区非独立性，不能把区域样本当作独立主分析。第三，GSE209552 样本量小且代表 severe TBI 外科切除脑组织，外推性有限。第四，GSE223245 是外周血数据，不能证明脑内机制。第五，8 基因综合评分是实验优先级工具，不是统计预测模型。第六，转录组证据不能替代蛋白二硫键化、F-actin 形态和细胞死亡检测。

## 结论

本研究将 TBI 后双硫死亡问题整理为一个可检验的转运-肌动蛋白骨架应激模型。当前公共转录组证据支持慢性 CTE stage 中 SLC7A11/SLC3A2/WASF2/TLN1 轴增强，急性 severe TBI 中骨架模块出现方向性改变，小鼠 CCI 中 3DPI hippocampus 是最集中的验证窗口，外周血严重程度数据可作为补充线索但不能替代脑组织证据。下一步应围绕 3DPI cortex/hippocampus 开展 SLC3A2/SLC7A11、WASF2/TLN1、ACTB/MYH9/MYL6/FLNA、F-actin、NADPH/GSH 和细胞共定位实验，以判断该候选转录组信号能否升级为真正的双硫死亡机制证据。

## 数据可用性

本研究使用 GEO 公共数据 GSE104687、GSE209552、GSE193407、GSE319253、GSE163415 和 GSE223245。所有重分析表格、PNG/PDF 图件、Markdown 和 Word 稿件均保存在 `Phase3_深化优化与最终报告_20260506-0513/11_双硫死亡聚焦论文设计_20260604/`。

## 伦理声明

本阶段仅使用公开去标识化数据和既有文献资料，不涉及新增人体样本。后续动物实验需在开展前获得所在单位动物伦理委员会批准。

## 参考文献

1. Liu X, Nie L, Zhang Y, et al. Actin cytoskeleton vulnerability to disulfide stress mediates disulfidptosis. Nature Cell Biology. 2023;25:404-414. https://doi.org/10.1038/s41556-023-01091-2  
2. Machesky LM. Deadly actin collapse by disulfidptosis. Nature Cell Biology. 2023.  
3. Zhao G, Zhao J, Lang J, Sun G. Nrf2 functions as a pyroptosis-related mediator in traumatic brain injury and is correlated with cytokines and disease severity: a bioinformatics analysis and retrospective clinical study. Frontiers in Neurology. 2024;15:1341342. https://doi.org/10.3389/fneur.2024.1341342  
4. Thomas I, Dickens AM, Posti JP, et al. Serum metabolome associated with severity of acute traumatic brain injury. Nature Communications. 2022;13:2545. https://doi.org/10.1038/s41467-022-30227-5  
5. Zhang M, Shan H, Chang P, et al. Hydrogen Sulfide Offers Neuroprotection on Traumatic Brain Injury in Parallel with Reduced Apoptosis and Autophagy in Mice. PLoS ONE. 2014;9:e87241. https://doi.org/10.1371/journal.pone.0087241  
6. Labadorf A, Agus F, Aytan N, et al. Inflammation and neuronal gene expression changes differ in early versus late chronic traumatic encephalopathy brain. BMC Medical Genomics. 2023;16:49. https://doi.org/10.1186/s12920-023-01471-5  
7. Garza R, Sharma Y, Atacho DAM, et al. Single-cell transcriptomics of human traumatic brain injury reveals activation of endogenous retroviruses in oligodendroglia. Cell Reports. 2023;42:113395. https://doi.org/10.1016/j.celrep.2023.113395  
8. Attilio PJ, Snapper DM, Rusnak M, et al. Transcriptomic Analysis of Mouse Brain After Traumatic Brain Injury Reveals That the Angiotensin Receptor Blocker Candesartan Acts Through Novel Pathways. Frontiers in Neuroscience. 2021;15:636259. https://doi.org/10.3389/fnins.2021.636259
"""


def en_manuscript(vals: dict[str, str], tables: dict[str, str]) -> str:
    ga = "GraphicalAbstract_v4_TBI_disulfidptosis_20260605"
    return f"""# Spatiotemporal Localization, Severity Context and Cell-Type Prioritization of Disulfidptosis-Like Transporter-Actin Cytoskeletal Stress After Traumatic Brain Injury: An Integrative Bioinformatic Study of Human TBI/CTE and Mouse CCI Transcriptomes

## Abstract

Background: Secondary injury after traumatic brain injury (TBI) is heterogeneous across time, brain region, injury severity and cellular context. Disulfidptosis is a recently defined cell-death mechanism driven by SLC7A11/SLC3A2-mediated cystine uptake, NADPH depletion and disulfide stress on the actin cytoskeleton. Public disease transcriptomes cannot directly prove disulfidptosis, but they can test whether a transporter-actin cytoskeletal stress program is recurrently perturbed after TBI.

Methods: We pre-specified an 8-gene panel. SLC3A2 and SLC7A11 represented the cystine-import entry point, whereas WASF2, TLN1, ACTB, MYH9, MYL6 and FLNA represented actin remodeling, adhesion/tension transfer and cytoskeletal endpoint biology. We integrated human acute severe TBI brain tissue (GSE209552), human CTE stage data (GSE193407), external CTE data (GSE319253), remote multi-region human TBI data (GSE104687), mouse controlled cortical impact (CCI) time-region data (GSE163415), and human peripheral blood mild/moderate/severe TBI data (GSE223245). FDR-supported findings were separated from nominal exploratory signals; donor/sample independence was respected; animal and peripheral blood data were interpreted only as supportive or contextual evidence.

Results: The strongest human-brain evidence came from chronic CTE stage analysis. In GSE193407, FDR-supported stage-trend signals were {vals['stage_details']}, and FDR-supported late CTE stage 3-4 versus stage 0 signals were {vals['late_details']}. In acute severe TBI brain tissue, all eight genes changed in the positive direction, with notable gene-level effects including {vals['acute_genes']}, but the small sample size limits inference. Mouse CCI validation-window ranking identified {vals['mouse_top']}; specifically, {vals['mouse_hipp']}, supporting 3DPI hippocampus as the most concentrated validation window. Peripheral GSE223245 added a severity-aware clue ({vals['severity_module']}; top gene-level trend: {vals['severity_gene']}), but this blood-based result cannot substitute for brain-mechanism evidence. Bulk marker-proxy analysis in GSE209552 highlighted {vals['marker_top']}, which should be used only to prioritize co-localization experiments. The top integrated validation candidate was {vals['priority_top']}.

Conclusions: Public transcriptomic evidence supports a candidate disulfidptosis-like transporter-actin cytoskeletal stress program after TBI/CTE. The most coherent dry-lab model is a chronic CTE-stage SLC7A11/SLC3A2/WASF2/TLN1 axis, directional acute severe TBI actin-module change, and a mouse CCI 3DPI hippocampal validation window. Mechanistic upgrading requires matched protein, F-actin, redox and cell-type co-localization evidence in the same time window and brain region.

Keywords: traumatic brain injury; chronic traumatic encephalopathy; disulfidptosis; SLC7A11; SLC3A2; actin cytoskeleton; injury time; injury severity; transcriptomics

## Graphical Abstract

<a id="fig-0"></a>

**Fig. 0. Graphical abstract of the revised manuscript.**

![Fig. 0](../figures/{ga}.png)

## Introduction

TBI pathophysiology is not a single molecular event. After the primary mechanical insult, local brain tissue rapidly enters metabolic crisis, ionic imbalance, excitotoxic stress, mitochondrial dysfunction and redox disequilibrium. Over subsequent days to weeks, blood-brain barrier disruption, glial activation, cytoskeletal remodeling, axonal and synaptic injury and persistent inflammation jointly shape secondary injury. In chronic or repetitive injury contexts, chronic traumatic encephalopathy (CTE) adds long-term neuroinflammation, cellular-composition shifts and neurodegenerative pathology. A mechanism-oriented transcriptomic study must therefore ask not only whether a pathway is perturbed, but also when, where, at what injury severity and in which cellular context that perturbation is plausible.

Disulfidptosis provides a mechanistic entry point linking amino-acid transport, redox failure and cytoskeletal vulnerability. In the original model, high SLC7A11 expression under glucose limitation promotes cystine uptake, consumes NADPH, produces abnormal disulfide stress and collapses the actin cytoskeleton. TBI brain tissue is characterized by energy failure, mitochondrial injury and redox imbalance, making a disulfidptosis-like state biologically plausible. However, transcriptomic evidence alone cannot demonstrate protein disulfidation, F-actin collapse or a specific cell-death modality. The present manuscript therefore uses the more precise term "disulfidptosis-like transporter-actin cytoskeletal stress."

Rather than mixing broad sulfur metabolism, oxidative stress and cell-death gene sets, this study narrows the hypothesis to eight pre-specified genes. SLC3A2 and SLC7A11 represent the cystine-import entry point; WASF2 and TLN1 represent WRC-linked actin regulation, adhesion and tension transfer; ACTB, MYH9, MYL6 and FLNA represent structural and contractile cytoskeletal endpoints. This focus sacrifices pathway breadth but improves testability: if a disulfidptosis-adjacent stress process is present after TBI/CTE, specific windows, regions or disease stages should show coordinated or stage-dependent perturbation across transporter-entry and actin-endpoint genes.

The goal of this study was to transform available public transcriptomes into a complete, experimentally actionable bioinformatic manuscript. The central claim is not that TBI has already been proven to cause disulfidptosis. Instead, we construct a tiered evidence model in which human chronic CTE data address disease-course association, acute severe TBI brain tissue addresses early severe injury directionality, mouse CCI identifies a testable time-region window, peripheral blood data add severity-aware context, and bulk marker-proxy analysis prioritizes cell-type localization experiments.

## Materials and Methods

### Dataset Inclusion and Evidence Tiers

Six GEO transcriptomic evidence layers were used. GSE209552 represented human acute severe TBI surgically sampled brain tissue obtained approximately 4 h to 8 d after injury. GSE193407 represented human BA9 CTE stage 0-4 and was used for chronic disease-course analysis. GSE319253 served as an external superior frontal cortex CTE versus control layer. GSE104687 represented remote multi-region postmortem TBI brain tissue and was treated only as regional exploratory evidence because multiple brain regions from the same donor are not independent observations. GSE163415 represented mouse controlled cortical impact across 3DPI/29DPI, hippocampus/thalamus/hypothalamus and vehicle/drug strata. GSE223245 represented whole blood/PBMC mild, moderate and severe TBI and was used only as a peripheral severity-aware clue.

{tables['datasets']}

### Candidate Genes, Modules and Statistical Principles

The candidate panel was fixed before analysis as SLC3A2, SLC7A11, WASF2, TLN1, ACTB, MYH9, MYL6 and FLNA. Three module readouts were calculated: the 8-gene total module, the SLC3A2/SLC7A11 transporter-entry module, and the WASF2/TLN1/ACTB/MYH9/MYL6/FLNA actin module. Module scores were calculated as the mean of sample-wise z-scored expression values. Gene-level comparisons used existing differential-expression outputs or Welch tests. Severity trends used Spearman correlation. Multiple testing was controlled using the Benjamini-Hochberg FDR. In the manuscript, FDR<0.05 is described as FDR-supported, whereas P<0.05 without FDR support is described as nominal exploratory evidence.

### Focused GSE223245 Severity Reanalysis

For GSE223245, the series matrix and platform SOFT annotation were used to reconstruct gene-level expression. Probes were collapsed to gene symbols by retaining the probe with the highest mean expression for each gene. Samples were classified as Control, Mild, Moderate and Severe and encoded as 0, 1, 2 and 3. For each gene, Spearman correlation was calculated across all groups and across TBI-only severity scores. Welch tests compared Mild, Moderate and Severe groups with controls. Because each group contained only four samples and the tissue source was peripheral blood, this analysis was used only as a severity-aware external context.

### GSE209552 Bulk Marker-Proxy Prioritization

The currently usable GSE209552 matrix is bulk-like brain tissue rather than an annotated single-nucleus object. To prioritize later tissue co-localization, marker-proxy scores were calculated for neurons, astrocytes, microglia, oligodendrocytes, OPCs and endothelial cells and correlated with the 8-gene, transporter and actin modules. This analysis cannot assign disulfidptosis to a cell type. It only identifies marker contexts that should be prioritized in IF/IHC or annotated snRNA-seq follow-up.

### Integrated Prioritization and Evidence Audit

The integrated priority score combined human FDR support, human nominal support, mouse CCI support, acute GSE209552 directionality and GSE223245 severity clues. It is intended for experimental prioritization, not clinical prediction or formal statistical inference. In this v4 revision, we also added comparison-level human evidence auditing, GSE163415 time-region-treatment validation-window ranking, treatment-stratified modulation and a mechanistic upgrade matrix.

## Results

### A Tiered Study Design Converted Disulfidptosis Into a Testable Transporter-Actin Axis

The graphical abstract and [Fig. 1](#fig-1) summarize the logic of the study. Disulfidptosis was not treated as an established fact in TBI; instead, it was decomposed into testable transcriptomic and experimental readouts: transporter entry, actin endpoint, redox prerequisite and cell-type localization. Each dataset had a distinct evidentiary role. This structure prevents animal or peripheral blood findings from being overinterpreted as human brain mechanism and prevents nominal signals from being upgraded into causal conclusions.

{figure_block(1, 'Fig1_v3_study_design_evidence_layers_20260604', 'Study design, evidence layers, pre-fixed 8-gene panel and interpretation guardrails.')}

### Human Brain Data Identified a Chronic CTE-Stage SLC7A11/SLC3A2/WASF2/TLN1 Axis

Across human brain datasets, the 8-gene panel was recurrent but heterogeneous. [Fig. 2](#fig-2) and Table 2 show that GSE193407 provided the strongest evidence tier. In the stage-trend analysis, FDR-supported genes were {vals['stage_details']}; in late CTE stage 3-4 versus stage 0, FDR-supported genes were {vals['late_details']}. GSE104687 provided only regional nominal clues and is limited by non-independent multi-region sampling from the same donors. GSE319253 can support chronic external validation but should remain secondary until a refreshed formal gene-level table is generated for submission.

{figure_block(2, 'Fig2_v3_human_brain_multidataset_8gene_20260604', 'Human brain multidataset evidence for the 8-gene transporter-actin panel.')}

Table 2. Human comparison-level evidence audit.

{tables['human']}

### Acute Severe TBI Showed Directional Actin-Axis Perturbation With Small-Sample Constraints

GSE209552 represents the acute severe TBI brain window. [Fig. 3](#fig-3) shows that all eight panel genes changed in the positive direction in TBI versus control. Notable gene-level effects included {vals['acute_genes']}. This supports the idea that acute severe TBI brain tissue already contains a transporter-actin directional perturbation, with cytoskeletal endpoints being particularly visible. However, the small sample size prevents biomarker or causal claims. The correct inference is directional support for experimental follow-up, not definitive molecular diagnosis.

{figure_block(3, 'Fig3_v3_GSE209552_acute_severe_and_marker_proxy_20260604', 'Acute severe TBI evidence in GSE209552, including module scores, gene-level directionality and bulk marker-proxy prioritization.')}

### Mouse CCI Localized the Highest-Priority Validation Window to 3DPI Hippocampus

GSE163415 was the key dataset for the "when and where" question. [Fig. 4](#fig-4) and [Fig. 7](#fig-7) show that the strongest time-region-treatment unit was {vals['mouse_top']}. When summarized by time and region, 3DPI hippocampus was the most concentrated window: {vals['mouse_hipp']}. This argues for a focused first wet-lab validation design centered on early subacute hippocampus and cortex rather than an evenly spread time-course. Later 29DPI signals can be used to assess persistence or transition, but they are less suitable as the sole discovery window.

{figure_block(4, 'Fig4_v3_GSE163415_mouse_CCI_spatiotemporal_suite_20260604', 'Mouse CCI spatiotemporal evidence in GSE163415.')}

Table 3. GSE163415 top time-region-treatment validation windows.

{tables['mouse']}

### Treatment-Stratified Mouse Data Suggested Different Early and Late Module Behavior

The treatment strata in GSE163415 were not used to make a treatment-mechanism claim, but they helped characterize stage-specific module behavior. [Fig. 7](#fig-7)C shows that Drug-Veh differences were most visible in the 29DPI core transporter/stress module, with the strongest stratified effect being {vals['treatment_top']}. This result should not be interpreted as drug efficacy. Instead, it suggests that 29DPI may retain transport/metabolic-axis modulation, whereas 3DPI hippocampus shows more concentrated transporter-actin co-perturbation. A practical validation strategy is to use 3DPI for mechanism discovery and 7D/29DPI for persistence checks.

{figure_block(7, 'Fig7_v4_evidence_audit_and_validation_priorities_20260605', 'Evidence audit, mouse validation-window ranking, treatment-stratified modulation and mechanistic upgrade priorities.')}

### Peripheral GSE223245 Added Severity Context Without Replacing Brain Evidence

Focused GSE223245 reanalysis added a severity-aware peripheral layer. The strongest module-level trend was {vals['severity_module']}, and the strongest gene-level trend was {vals['severity_gene']}. [Fig. 5](#fig-5) and [Fig. 8](#fig-8)D show that SLC3A2 and FLNA decreased across all-group severity, whereas MYL6 showed a positive TBI-only trend. Because the samples are whole blood/PBMC and each group contains only four samples, these results are best interpreted as peripheral severity clues. They justify including graded injury or severity strata in future designs but cannot demonstrate brain disulfidptosis.

{figure_block(5, 'Fig5_v3_GSE223245_peripheral_severity_suite_20260604', 'Peripheral severity-focused analysis of GSE223245 whole blood/PBMC data.')}

Table 4. GSE223245 8-gene severity trends.

{tables['severity']}

### Bulk Marker-Proxy Analysis Prioritized Endothelial and Neuronal Co-Localization

The strongest GSE209552 marker-proxy result was {vals['marker_top']}. The 8-gene module also showed a negative association with the neuronal proxy. [Fig. 8](#fig-8)B visualizes the full correlation matrix. Because bulk marker-proxy correlations are affected by cell proportion, sampling, vascular content and cell-state shifts, they cannot assign the pathway to endothelial cells or neurons. Their value is experimental prioritization: SLC3A2/SLC7A11 and F-actin readouts should be co-stained with CD31, NeuN, GFAP, IBA1 and OLIG2 rather than tested in only one cell class.

Table 5. Top GSE209552 bulk marker-proxy prioritization clues.

{tables['marker']}

### Integrated Gene Prioritization Separated Expression Readouts From Cytoskeletal Endpoints

The integrated score ranked {vals['priority_top']} as the leading validation candidate. [Fig. 6](#fig-6) and Table 6 show that SLC3A2/SLC7A11/WASF2/TLN1 are best suited as primary expression readouts because they receive stronger human CTE-stage support. ACTB/MYH9/MYL6/FLNA are closer to cytoskeletal endpoints and should be interpreted together with protein abundance, non-reducing assays, cytoskeletal fractionation and F-actin morphology. This distinction is important: transporter-entry genes ask whether the entry point is activated, whereas cytoskeletal endpoint genes ask whether the structural consequence is present.

{figure_block(6, 'Fig6_v3_integrated_priority_validation_suite_20260604', 'Integrated 8-gene priority matrix, validation ranking, mechanistic upgrade workflow and marker-proxy summary.')}

Table 6. Integrated 8-gene evidence priority.

{tables['priority']}

### Mechanistic Upgrading Requires Transporter, Cytoskeletal, Redox and Localization Evidence

[Fig. 8](#fig-8) formalizes the claim ladder used in the manuscript. Current evidence reaches public brain transcriptome plus animal supportive evidence, but not mechanistic proof. To upgrade the model, the same time window and brain region must show transporter protein changes, cytoskeletal endpoint abnormalities, F-actin morphology changes, NADPH/NADP+ and GSH/GSSG alterations, and cell-type co-localization. Table 7 summarizes the experimental upgrade matrix.

{figure_block(8, 'Fig8_v4_claim_ladder_cell_proxy_and_wetlab_design_20260605', 'Claim ladder, cell-proxy heatmap, wet-lab time-course design and peripheral severity gene trends.')}

Table 7. Mechanistic upgrade matrix for the next experimental stage.

{tables['validation']}

## Discussion

The main contribution of this study is not the nomination of a single TBI gene, but the translation of a new cell-death mechanism into a testable TBI transcriptomic framework. Disulfidptosis is not defined by SLC7A11 or SLC3A2 expression alone. It requires a chain linking cystine uptake, reducing-power depletion, disulfide stress and actin cytoskeletal collapse. Public transcriptomes can capture only the entry and part of the endpoint, so the manuscript deliberately uses "disulfidptosis-like transporter-actin cytoskeletal stress" rather than claiming that disulfidptosis has already occurred.

Human brain results support a stage-dependent model. In chronic CTE, the SLC7A11/SLC3A2/WASF2/TLN1 axis shows relatively stable FDR-supported association with disease stage, suggesting that cystine transport and actin regulatory/tension nodes may be linked to chronic disease progression. In acute severe TBI, cytoskeletal endpoint directionality is more visible, which may reflect mechanical stress, vascular/glial components, cytoskeletal disruption and cell-composition shifts. These two observations are complementary rather than contradictory: chronic disease may accumulate transporter and regulatory-axis changes, whereas acute severe injury may expose cytoskeletal endpoint stress more directly.

The mouse CCI analysis gives the clearest experimental guidance. Because 3DPI hippocampus concentrates the strongest panel-wide signal, the first wet-lab round should not diffuse effort across too many time points and regions. A coherent design would prioritize 3DPI cortex/hippocampus for qPCR, WB, F-actin staining, NADPH/GSH assays and co-localization, with 7D or 29DPI used to assess persistence or transition. If this first round does not show concordant transporter-entry and cytoskeletal-endpoint changes, the model should be reframed as broader cytoskeletal remodeling or redox stress rather than disulfidptosis-like stress.

Severity requires careful interpretation. Severe TBI in GSE209552, CTE stage in GSE193407 and mild/moderate/severe blood groups in GSE223245 are not the same biological scale. Combining them into a single severity score would create a misleading statistical construct. A better manuscript structure is layered: GSE209552 represents acute severe brain injury, GSE193407 represents chronic disease course, and GSE223245 represents peripheral severity context. This organization directly addresses injury severity while preserving biological meaning.

Cell-type localization remains the largest evidence gap. The endothelial marker-proxy association suggests that vascular or barrier-related compartments deserve priority, while the neuronal negative association indicates that neuronal abundance or state shifts should not be ignored. Nevertheless, bulk marker-proxy analysis cannot distinguish cell proportion from within-cell expression change. Annotated snRNA-seq extraction from GSE209552, if available, or tissue co-localization of SLC3A2/SLC7A11 with NeuN, GFAP, IBA1, OLIG2 and CD31 would materially strengthen the manuscript.

This study has limitations. Public datasets differ in platform, species, tissue source, injury definition and covariate availability. GSE104687 multi-region samples are not independent donor-level observations. GSE209552 is small and represents severe surgically sampled brain tissue. GSE223245 is peripheral blood and cannot replace brain evidence. The integrated score is an experimental prioritization tool rather than a formal predictive model. Finally, transcriptomics cannot substitute for protein disulfidation, F-actin morphology or cell-death assays.

## Conclusion

This integrative analysis organizes post-TBI disulfidptosis into a testable transporter-actin cytoskeletal stress model. Current public transcriptomic evidence supports a chronic CTE-stage SLC7A11/SLC3A2/WASF2/TLN1 axis, directional acute severe TBI cytoskeletal perturbation, and a mouse CCI 3DPI hippocampal validation window. Peripheral severity data are useful context but not brain-mechanism evidence. The next experimental step should test SLC3A2/SLC7A11, WASF2/TLN1, ACTB/MYH9/MYL6/FLNA, F-actin, NADPH/GSH and cell-type co-localization in matched 3DPI cortex/hippocampus tissue.

## Data Availability

This study used public GEO datasets GSE104687, GSE209552, GSE193407, GSE319253, GSE163415 and GSE223245. Reanalysis tables, PNG/PDF figures, Markdown manuscripts and Word files are stored under `Phase3_深化优化与最终报告_20260506-0513/11_双硫死亡聚焦论文设计_20260604/`.

## Ethics Statement

This stage used only public de-identified datasets and literature materials. Any future animal experiment should obtain institutional animal ethics approval before initiation.

## References

1. Liu X, Nie L, Zhang Y, et al. Actin cytoskeleton vulnerability to disulfide stress mediates disulfidptosis. Nature Cell Biology. 2023;25:404-414. https://doi.org/10.1038/s41556-023-01091-2  
2. Machesky LM. Deadly actin collapse by disulfidptosis. Nature Cell Biology. 2023.  
3. Zhao G, Zhao J, Lang J, Sun G. Nrf2 functions as a pyroptosis-related mediator in traumatic brain injury and is correlated with cytokines and disease severity. Frontiers in Neurology. 2024;15:1341342. https://doi.org/10.3389/fneur.2024.1341342  
4. Thomas I, Dickens AM, Posti JP, et al. Serum metabolome associated with severity of acute traumatic brain injury. Nature Communications. 2022;13:2545. https://doi.org/10.1038/s41467-022-30227-5  
5. Zhang M, Shan H, Chang P, et al. Hydrogen Sulfide Offers Neuroprotection on Traumatic Brain Injury in Parallel with Reduced Apoptosis and Autophagy in Mice. PLoS ONE. 2014;9:e87241. https://doi.org/10.1371/journal.pone.0087241  
6. Labadorf A, Agus F, Aytan N, et al. Inflammation and neuronal gene expression changes differ in early versus late chronic traumatic encephalopathy brain. BMC Medical Genomics. 2023;16:49. https://doi.org/10.1186/s12920-023-01471-5  
7. Garza R, Sharma Y, Atacho DAM, et al. Single-cell transcriptomics of human traumatic brain injury reveals activation of endogenous retroviruses in oligodendroglia. Cell Reports. 2023;42:113395. https://doi.org/10.1016/j.celrep.2023.113395  
8. Attilio PJ, Snapper DM, Rusnak M, et al. Transcriptomic Analysis of Mouse Brain After Traumatic Brain Injury Reveals That the Angiotensin Receptor Blocker Candesartan Acts Through Novel Pathways. Frontiers in Neuroscience. 2021;15:636259. https://doi.org/10.3389/fnins.2021.636259
"""


def write_reference_and_orchestra_notes(vals: dict[str, str]) -> None:
    refs = [
        ("Zhao 2024 TBI bioinformatics-validation template", "09_Zhao2024_TBI_bioinfo_clinical_validation_template.pdf", "borrowed candidate-screening, module ranking, clinical severity and validation narrative structure"),
        ("Thomas 2022 acute TBI severity omics template", "10_Thomas2022_acute_TBI_metabolomics_template.pdf", "borrowed severity/time-window caution and high-density results prose"),
        ("Liu 2023 disulfidptosis mechanism", "01_Liu2023_original_disulfidptosis.pdf", "borrowed SLC7A11/SLC3A2-NADPH-actin mechanistic logic"),
        ("Machesky 2023 disulfidptosis commentary", "05_Machesky2023_NCB_deadly_actin_collapse_disulfidptosis.pdf", "borrowed conceptual emphasis on actin collapse"),
        ("Zhang 2014 mouse TBI model", "04_Zhang2013_H2S_neuroprotection_TBI.pdf", "borrowed mouse TBI time-course and cortex/hippocampus validation style"),
    ]
    lines = ["# v4 reference-format papers and writing-tool audit", ""]
    for title, fname, note in refs:
        pdf = WORKDIR / "references" / fname
        lines.append(f"- [{title}]({pdf.resolve()}): {note}.")
    lines.extend(
        [
            "",
            "PaperOrchestra-style organization used in v4:",
            "- Idea: a fixed 8-gene disulfidptosis-like transporter-actin hypothesis after TBI.",
            "- Experimental log: human brain, mouse CCI, peripheral severity, marker-proxy and integrated ranking tables.",
            "- Plotting plan: graphical abstract plus Fig. 1-8, including two new v4 evidence-audit figures.",
            "- Refinement rule: every claim is assigned to one of four levels: peripheral clue, public brain transcriptome, tissue validation, or mechanistic proof.",
            "",
            f"Most important numerical anchors: {vals['stage_details']}; {vals['mouse_top']}; {vals['severity_module']}; {vals['marker_top']}.",
        ]
    )
    (REPORTDIR / "reference_format_papers_for_review_v4_20260605.md").write_text("\n".join(lines), encoding="utf-8")


def write_paper_orchestra_inputs(data: dict[str, pd.DataFrame], vals: dict[str, str]) -> None:
    ws = WORKDIR / "paper_orchestra_workspace_v4_20260605"
    inputs = ws / "inputs"
    inputs.mkdir(parents=True, exist_ok=True)
    (inputs / "idea.md").write_text(
        f"""# Dense Idea

We study whether TBI/CTE public transcriptomes support a disulfidptosis-like transporter-actin cytoskeletal stress program. The hypothesis is deliberately constrained to eight genes: SLC3A2, SLC7A11, WASF2, TLN1, ACTB, MYH9, MYL6 and FLNA. The manuscript asks when, where, at what injury severity and in which cellular context the signal should be validated.

Key claim boundary: transcriptomics can support a candidate disulfidptosis-like state but cannot prove disulfidptosis without protein, F-actin, redox and cell-type co-localization evidence.
""",
        encoding="utf-8",
    )
    (inputs / "experimental_log.md").write_text(
        f"""# Experimental Log

Human chronic CTE stage evidence: {vals['stage_details']}; late CTE comparison: {vals['late_details']}.

Acute severe TBI evidence: {vals['acute_genes']}.

Mouse CCI validation window: {vals['mouse_top']}; {vals['mouse_hipp']}.

Peripheral severity clue: {vals['severity_module']}; top gene trend {vals['severity_gene']}.

Cell-type prioritization clue: {vals['marker_top']}.
""",
        encoding="utf-8",
    )
    (inputs / "conference_guidelines.md").write_text(
        "Use IMRAD structure, full paragraphs, conservative causal language, clear FDR-vs-nominal labeling, complete figure legends, and explicit limitations for bulk, animal and peripheral datasets.\n",
        encoding="utf-8",
    )
    (inputs / "template.tex").write_text(
        r"""\documentclass[11pt]{article}
\usepackage{graphicx}
\usepackage{booktabs}
\begin{document}
\section{Introduction}
\section{Materials and Methods}
\section{Results}
\section{Discussion}
\section{Conclusion}
\end{document}
""",
        encoding="utf-8",
    )


def main() -> None:
    configure_style()
    data = load_inputs()
    ext = build_extended_tables(data)
    draw_graphical_abstract(data)
    draw_fig7(data, ext)
    draw_fig8(data, ext)
    vals = build_key_values(data, ext)
    tables = build_table_strings(data, ext)
    zh = zh_manuscript(vals, tables)
    en = en_manuscript(vals, tables)

    zh_md = REPORTDIR / "TBI_disulfidptosis_optimized_manuscript_v4_ZH_20260605.md"
    en_md = REPORTDIR / "TBI_disulfidptosis_optimized_manuscript_v4_EN_20260605.md"
    zh_docx = REPORTDIR / "TBI_disulfidptosis_optimized_manuscript_v4_ZH_20260605.docx"
    en_docx = REPORTDIR / "TBI_disulfidptosis_optimized_manuscript_v4_EN_20260605.docx"
    zh_md.write_text(zh, encoding="utf-8")
    en_md.write_text(en, encoding="utf-8")
    v3.markdown_to_docx(zh, zh_docx)
    v3.markdown_to_docx(en, en_docx)

    write_reference_and_orchestra_notes(vals)
    write_paper_orchestra_inputs(data, vals)

    artifacts = [
        ("report", zh_md, "Chinese optimized v4 Markdown manuscript"),
        ("report", zh_docx, "Chinese optimized v4 Word manuscript"),
        ("report", en_md, "English optimized v4 Markdown manuscript"),
        ("report", en_docx, "English optimized v4 Word manuscript"),
        ("report", REPORTDIR / "reference_format_papers_for_review_v4_20260605.md", "Reference-format and writing-tool audit"),
        ("figure", FIGDIR / "GraphicalAbstract_v4_TBI_disulfidptosis_20260605.png", "Graphical abstract PNG"),
        ("figure", FIGDIR / "GraphicalAbstract_v4_TBI_disulfidptosis_20260605.pdf", "Graphical abstract PDF"),
        ("figure", FIGDIR / "Fig7_v4_evidence_audit_and_validation_priorities_20260605.png", "New v4 evidence-audit multi-panel figure PNG"),
        ("figure", FIGDIR / "Fig7_v4_evidence_audit_and_validation_priorities_20260605.pdf", "New v4 evidence-audit multi-panel figure PDF"),
        ("figure", FIGDIR / "Fig8_v4_claim_ladder_cell_proxy_and_wetlab_design_20260605.png", "New v4 claim-ladder and wet-lab design multi-panel figure PNG"),
        ("figure", FIGDIR / "Fig8_v4_claim_ladder_cell_proxy_and_wetlab_design_20260605.pdf", "New v4 claim-ladder and wet-lab design multi-panel figure PDF"),
    ]
    for name in ["human_comp", "mouse_units", "treatment_module", "marker_priority", "validation_matrix"]:
        artifacts.append(("table", TABLEDIR / f"v4_{name}_20260605.csv", f"v4 derived analysis table: {name}"))
    artifacts.append(("script", Path(__file__), "v4 manuscript generation script"))
    manifest = pd.DataFrame(
        [{"type": kind, "path": str(path), "description": desc, "exists": path.exists()} for kind, path, desc in artifacts]
    )
    manifest_path = REPORTDIR / "optimized_v4_artifact_manifest_20260605.csv"
    manifest.to_csv(manifest_path, index=False, encoding="utf-8-sig")
    print(f"Wrote {zh_docx}")
    print(f"Wrote {en_docx}")
    print(f"Wrote {manifest_path}")


if __name__ == "__main__":
    main()
