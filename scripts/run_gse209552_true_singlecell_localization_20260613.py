from __future__ import annotations

import json
import math
import tarfile
from dataclasses import dataclass
from pathlib import Path

import anndata as ad
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
from scipy import sparse

np.random.seed(0)


PRIORITY_GENES = ["SLC3A2", "SLC7A11", "WASF2", "TLN1", "ACTB", "MYH9", "MYL6", "FLNA"]
CORE_FOCUS_GENES = ["SLC3A2", "SLC7A11", "WASF2", "TLN1"]
CELL_TYPE_MARKERS: dict[str, list[str]] = {
    "Neuron": ["RBFOX3", "SNAP25", "SYT1", "MAP2", "STMN2", "GRIN1"],
    "Astrocyte": ["GFAP", "AQP4", "ALDH1L1", "SLC1A3", "GJA1"],
    "Microglia": ["P2RY12", "CX3CR1", "AIF1", "CSF1R", "C1QC", "TYROBP"],
    "Oligodendrocyte": ["MBP", "MOG", "PLP1", "MAG", "MOBP"],
    "OPC": ["PDGFRA", "CSPG4", "VCAN", "SOX10", "OLIG1"],
    "Endothelial": ["PECAM1", "VWF", "CLDN5", "FLT1", "KDR", "EMCN"],
    "Pericyte": ["PDGFRB", "RGS5", "MCAM", "CSPG4", "ACTA2"],
    "Immune": ["PTPRC", "LST1", "FCER1G", "TYROBP", "AIF1"],
    "Ependymal": ["FOXJ1", "PIFO", "TTR", "CFAP44", "S100B"],
}
OKABE = {
    "blue": "#0072B2",
    "orange": "#E69F00",
    "green": "#009E73",
    "red": "#D55E00",
    "purple": "#CC79A7",
    "sky": "#56B4E9",
    "yellow": "#F0E442",
    "black": "#000000",
    "gray": "#6B7280",
}
CELL_TYPE_COLORS = {
    "Neuron": "#1b9e77",
    "Astrocyte": "#7570b3",
    "Microglia": "#d95f02",
    "Oligodendrocyte": "#66a61e",
    "OPC": "#e6ab02",
    "Endothelial": "#1f78b4",
    "Pericyte": "#a6761d",
    "Immune": "#e7298a",
    "Ependymal": "#666666",
    "Unknown": "#bdbdbd",
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


def find_file(name: str, contains: str | None = None) -> Path:
    matches = sorted(Path.cwd().rglob(name))
    if contains:
        matches = [p for p in matches if contains in str(p)]
    if not matches:
        raise FileNotFoundError(name)
    return matches[0]


WORKDIR = find_file(
    "TBI_disulfidptosis_focused_design_tables_20260604.xlsx",
    contains="11_双硫死亡聚焦论文设计_20260604",
).parents[1]
RAW_BASE = WORKDIR / "raw_data_v3_20260604" / "GSE209552_scRNA_raw"
SAMPLE_DIR = RAW_BASE / "samples"
TABLE_DIR = WORKDIR / "tables"
FIG_DIR = WORKDIR / "figures"
REPORT_DIR = WORKDIR / "reports"
for directory in (RAW_BASE, SAMPLE_DIR, TABLE_DIR, FIG_DIR, REPORT_DIR):
    directory.mkdir(parents=True, exist_ok=True)

RAW_TAR = RAW_BASE / "GSE209552_RAW.tar"
GEO_JSON = find_file("geo_candidate_records_raw.json", contains="01_TBI公共数据集目录")


# Override the path discovery above with the stable project structure so reruns do
# not depend on terminal-encoding-sensitive substring matches.
PROJECT_ROOT = Path(__file__).resolve().parents[1]
WORKDIR = PROJECT_ROOT / "Phase3_深化优化与最终报告_20260506-0513" / "11_双硫死亡聚焦论文设计_20260604"
RAW_BASE = WORKDIR / "raw_data_v3_20260604" / "GSE209552_scRNA_raw"
SAMPLE_DIR = RAW_BASE / "samples"
TABLE_DIR = WORKDIR / "tables"
FIG_DIR = WORKDIR / "figures"
REPORT_DIR = WORKDIR / "reports"
for directory in (RAW_BASE, SAMPLE_DIR, TABLE_DIR, FIG_DIR, REPORT_DIR):
    directory.mkdir(parents=True, exist_ok=True)

RAW_TAR = RAW_BASE / "GSE209552_RAW.tar"
GEO_JSON = PROJECT_ROOT / "Phase1_初始筛选与差异分析_20260420-0424" / "01_TBI公共数据集目录" / "geo_candidate_records_raw.json"


@dataclass
class SampleMeta:
    accession: str
    title: str
    donor_id: str
    condition: str
    region_label: str


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


def sample_metadata() -> dict[str, SampleMeta]:
    data = json.loads(GEO_JSON.read_text(encoding="utf-8"))
    target = next(rec for rec in data if rec.get("accession") == "GSE209552")
    out: dict[str, SampleMeta] = {}
    for sample in target["samples"]:
        accession = sample["accession"]
        title = sample["title"]
        condition = "Control" if title.startswith("Control") else "TBI"
        digits = "".join(ch for ch in title.split(",")[0] if ch.isdigit())
        donor_id = digits or accession.replace("GSM", "")
        title_lower = title.lower()
        if "frontal" in title_lower:
            region_label = "Frontal"
        elif "temporal" in title_lower:
            region_label = "Temporal"
        else:
            region_label = "TBI_resection"
        out[accession] = SampleMeta(
            accession=accession,
            title=title,
            donor_id=donor_id,
            condition=condition,
            region_label=region_label,
        )
    return out


def extract_sample_matrices() -> list[Path]:
    if not RAW_TAR.exists():
        raise FileNotFoundError(f"Missing raw tar: {RAW_TAR}")
    sample_dirs = sorted([p for p in SAMPLE_DIR.iterdir() if p.is_dir() and (p / "matrix.mtx.gz").exists()])
    if sample_dirs:
        return sample_dirs

    needed_suffixes = {
        "_barcodes.tsv.gz": "barcodes.tsv.gz",
        "_features.tsv.gz": "features.tsv.gz",
        "_matrix.mtx.gz": "matrix.mtx.gz",
    }
    with tarfile.open(RAW_TAR, "r") as tar:
        members = [m for m in tar.getmembers() if any(m.name.endswith(suffix) for suffix in needed_suffixes)]
        for member in members:
            accession = member.name.split("_", 1)[0]
            out_dir = SAMPLE_DIR / accession
            out_dir.mkdir(parents=True, exist_ok=True)
            extracted = tar.extractfile(member)
            if extracted is None:
                continue
            for suffix, target_name in needed_suffixes.items():
                if member.name.endswith(suffix):
                    target = out_dir / target_name
                    if not target.exists():
                        target.write_bytes(extracted.read())
                    break
    sample_dirs = sorted([p for p in SAMPLE_DIR.iterdir() if p.is_dir() and (p / "matrix.mtx.gz").exists()])
    if not sample_dirs:
        raise RuntimeError("No sample matrices extracted from GSE209552_RAW.tar")
    return sample_dirs


def infer_location_code(accession_dir: Path) -> str:
    original = next(accession_dir.glob("*"), None)
    if original is None:
        return "NA"
    name = original.name
    tokens = name.replace(".gz", "").split("_")
    for token in reversed(tokens):
        if token in {"F", "T"}:
            return token
        if token in {"RF", "LF", "RT", "LT", "LFP", "IFP"}:
            return token
    return "NA"


def load_samples() -> ad.AnnData:
    metas = sample_metadata()
    sample_dirs = extract_sample_matrices()
    adatas: list[ad.AnnData] = []
    for sample_dir in sample_dirs:
        accession = sample_dir.name
        meta = metas.get(accession)
        if meta is None:
            continue
        adata = sc.read_10x_mtx(sample_dir, var_names="gene_symbols", make_unique=True, cache=False)
        adata.var_names_make_unique()
        adata.obs["sample_accession"] = accession
        adata.obs["sample_title"] = meta.title
        adata.obs["condition"] = meta.condition
        adata.obs["donor_id"] = meta.donor_id
        adata.obs["region_label"] = meta.region_label
        adata.obs["location_code"] = infer_location_code(sample_dir)
        adata.obs_names = [f"{cell}_{accession}" for cell in adata.obs_names.astype(str)]
        adatas.append(adata)
    if not adatas:
        raise RuntimeError("No GSE209552 single-cell matrices loaded")
    combined = ad.concat(adatas, join="outer", label="concat_batch", keys=[a.obs["sample_accession"].iloc[0] for a in adatas], merge="same")
    combined.var_names_make_unique()
    combined.obs_names = [f"{cell}_{sample}" for cell, sample in zip(combined.obs_names.astype(str), combined.obs["sample_accession"].astype(str))]
    combined.layers["counts"] = combined.X.copy()
    return combined


def save_qc_summaries(adata_raw: ad.AnnData, adata_filtered: ad.AnnData) -> None:
    before = adata_raw.obs.groupby("sample_accession").agg(
        cells_before=("sample_accession", "size"),
        median_genes_before=("n_genes_by_counts", "median"),
        median_counts_before=("total_counts", "median"),
        mt_pct_before=("pct_counts_mt", "median"),
    )
    after = adata_filtered.obs.groupby("sample_accession").agg(
        cells_after=("sample_accession", "size"),
        median_genes_after=("n_genes_by_counts", "median"),
        median_counts_after=("total_counts", "median"),
        mt_pct_after=("pct_counts_mt", "median"),
    )
    qc = before.join(after, how="left").reset_index()
    qc["retained_fraction"] = qc["cells_after"] / qc["cells_before"]
    qc.to_csv(TABLE_DIR / "GSE209552_true_scRNA_qc_summary_20260613.csv", index=False, encoding="utf-8-sig")


def raw_gene_df(adata: ad.AnnData, genes: list[str]) -> pd.DataFrame:
    available = [g for g in genes if g in adata.raw.var_names]
    raw_subset = adata.raw[:, available]
    matrix = raw_subset.X.toarray() if sparse.issparse(raw_subset.X) else np.asarray(raw_subset.X)
    return pd.DataFrame(matrix, index=adata.obs_names, columns=available)


def preprocess(adata: ad.AnnData) -> ad.AnnData:
    adata.var["mt"] = adata.var_names.str.upper().str.startswith("MT-")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)
    raw_copy = adata.copy()
    sc.pp.filter_cells(adata, min_genes=150)
    sc.pp.filter_genes(adata, min_cells=5)
    adata = adata[adata.obs["pct_counts_mt"] < 5].copy()
    save_qc_summaries(raw_copy, adata)

    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    adata.raw = adata
    sc.pp.highly_variable_genes(adata, flavor="seurat", n_top_genes=min(3000, max(500, adata.n_vars - 1)))
    adata = adata[:, adata.var["highly_variable"]].copy()
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, svd_solver="arpack")
    n_pcs = min(30, adata.obsm["X_pca"].shape[1])
    sc.pp.neighbors(adata, n_neighbors=15, n_pcs=n_pcs)
    sc.tl.umap(adata, min_dist=0.35, random_state=0)
    sc.tl.leiden(adata, resolution=0.55, key_added="leiden", flavor="igraph", n_iterations=2, directed=False)
    sc.tl.rank_genes_groups(adata, "leiden", method="wilcoxon", key_added="rank_genes_leiden")
    return adata


def cluster_marker_table(adata: ad.AnnData, n_genes: int = 12) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for cluster in sorted(adata.obs["leiden"].unique(), key=lambda x: int(x)):
        df = sc.get.rank_genes_groups_df(adata, group=cluster, key="rank_genes_leiden").head(n_genes).copy()
        df["cluster"] = cluster
        rows.append(df)
    markers = pd.concat(rows, ignore_index=True)
    markers.to_csv(TABLE_DIR / "GSE209552_true_scRNA_cluster_markers_20260613.csv", index=False, encoding="utf-8-sig")
    return markers


def annotate_cell_types(adata: ad.AnnData) -> tuple[ad.AnnData, pd.DataFrame]:
    for cell_type, genes in CELL_TYPE_MARKERS.items():
        present = [g for g in genes if g in adata.raw.var_names]
        if present:
            sc.tl.score_genes(adata, gene_list=present, score_name=f"score_{cell_type}", use_raw=True)
        else:
            adata.obs[f"score_{cell_type}"] = 0.0

    score_cols = [f"score_{ct}" for ct in CELL_TYPE_MARKERS]
    cluster_scores = adata.obs.groupby("leiden", observed=True)[score_cols].mean()
    cluster_assignments = {}
    for cluster, row in cluster_scores.iterrows():
        best_col = row.astype(float).idxmax()
        best_score = float(row[best_col])
        second_score = float(row.astype(float).sort_values(ascending=False).iloc[1]) if len(row) > 1 else -999
        best_type = best_col.replace("score_", "")
        if best_score < 0.05 or (best_score - second_score) < 0.02:
            cluster_assignments[cluster] = "Unknown"
        else:
            cluster_assignments[cluster] = best_type

    adata.obs["cell_type"] = adata.obs["leiden"].map(cluster_assignments).astype("category")
    cluster_summary = cluster_scores.copy()
    cluster_summary["assigned_cell_type"] = cluster_summary.index.map(cluster_assignments)
    cluster_summary.to_csv(TABLE_DIR / "GSE209552_true_scRNA_cluster_celltype_scores_20260613.csv", encoding="utf-8-sig")
    return adata, cluster_summary.reset_index().rename(columns={"leiden": "cluster"})


def summarize_gene_localization(adata: ad.AnnData) -> tuple[pd.DataFrame, pd.DataFrame]:
    expr = raw_gene_df(adata, PRIORITY_GENES)
    obs = adata.obs[["cell_type", "condition", "donor_id", "sample_accession", "sample_title", "region_label"]].copy()
    valid_cell_types = sorted(pd.Series(obs["cell_type"].astype(str)).replace("nan", np.nan).dropna().unique())
    rows = []
    donor_rows = []
    for cell_type in valid_cell_types:
        idx_cell = obs.index[obs["cell_type"].astype(str) == cell_type]
        if len(idx_cell) == 0:
            continue
        for condition in sorted(obs.loc[idx_cell, "condition"].unique()):
            idx = obs.index[(obs["cell_type"].astype(str) == cell_type) & (obs["condition"] == condition)]
            if len(idx) == 0:
                continue
            sub = expr.loc[idx, [g for g in PRIORITY_GENES if g in expr.columns]]
            for gene in sub.columns:
                rows.append(
                    {
                        "cell_type": cell_type,
                        "condition": condition,
                        "gene": gene,
                        "n_cells": len(sub),
                        "mean_log1p_expr": float(sub[gene].mean()),
                        "pct_expressed": float((sub[gene] > 0).mean() * 100),
                    }
                )
        for donor_id, donor_df in obs.loc[idx_cell].groupby("donor_id", observed=True):
            if donor_df.empty:
                continue
            idx_donor = donor_df.index
            sub = expr.loc[idx_donor, [g for g in PRIORITY_GENES if g in expr.columns]]
            condition = donor_df["condition"].iloc[0]
            for gene in sub.columns:
                donor_rows.append(
                    {
                        "donor_id": donor_id,
                        "condition": condition,
                        "cell_type": cell_type,
                        "gene": gene,
                        "n_cells": len(sub),
                        "mean_log1p_expr": float(sub[gene].mean()),
                        "pct_expressed": float((sub[gene] > 0).mean() * 100),
                    }
                )
    gene_summary = pd.DataFrame(rows)
    donor_summary = pd.DataFrame(donor_rows)
    gene_summary.to_csv(TABLE_DIR / "GSE209552_true_scRNA_8gene_celltype_localization_20260613.csv", index=False, encoding="utf-8-sig")
    donor_summary.to_csv(TABLE_DIR / "GSE209552_true_scRNA_8gene_donor_level_localization_20260613.csv", index=False, encoding="utf-8-sig")
    return gene_summary, donor_summary


def donor_level_deltas(donor_summary: pd.DataFrame) -> pd.DataFrame:
    rows = []
    grouped = donor_summary.groupby(["cell_type", "gene"])
    for (cell_type, gene), df in grouped:
        tbi = df.loc[df["condition"] == "TBI", "mean_log1p_expr"]
        ctrl = df.loc[df["condition"] == "Control", "mean_log1p_expr"]
        if len(tbi) == 0 or len(ctrl) == 0:
            continue
        rows.append(
            {
                "cell_type": cell_type,
                "gene": gene,
                "tbi_donor_n": int(len(tbi)),
                "control_donor_n": int(len(ctrl)),
                "tbi_mean": float(tbi.mean()),
                "control_mean": float(ctrl.mean()),
                "delta_tbi_minus_control": float(tbi.mean() - ctrl.mean()),
                "tbi_median": float(tbi.median()),
                "control_median": float(ctrl.median()),
                "pct_tbi_mean": float(df.loc[df["condition"] == "TBI", "pct_expressed"].mean()),
                "pct_control_mean": float(df.loc[df["condition"] == "Control", "pct_expressed"].mean()),
            }
        )
    delta_df = pd.DataFrame(rows).sort_values(["delta_tbi_minus_control", "tbi_mean"], ascending=[False, False])
    delta_df.to_csv(TABLE_DIR / "GSE209552_true_scRNA_8gene_donor_level_deltas_20260613.csv", index=False, encoding="utf-8-sig")
    return delta_df


def make_main_figure(
    adata: ad.AnnData,
    gene_summary: pd.DataFrame,
    delta_df: pd.DataFrame,
    cluster_summary: pd.DataFrame,
) -> None:
    fig = plt.figure(figsize=(14, 10))
    gs = fig.add_gridspec(3, 4, height_ratios=[0.9, 1.6, 1.3], width_ratios=[1.1, 1.2, 1.2, 1.2], hspace=0.4, wspace=0.35)

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
        legend_loc="on data",
        legend_fontsize=8,
        legend_fontoutline=2,
        title="UMAP of annotated nuclei",
    )

    ax_c = fig.add_subplot(gs[1, 2])
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
        ax_c.bar(comp.index, vals, bottom=bottom, color=CELL_TYPE_COLORS.get(col, "#bdbdbd"), edgecolor="white", linewidth=0.5, label=col)
        bottom += vals
    ax_c.set_ylabel("Cell fraction (%)")
    ax_c.set_title("Cell-type composition by condition")
    ax_c.tick_params(axis="x", rotation=0)
    ax_c.legend(loc="upper left", bbox_to_anchor=(1.02, 1.0), frameon=False)

    ax_d = fig.add_subplot(gs[1, 3])
    top_delta = delta_df.loc[delta_df["gene"].isin(CORE_FOCUS_GENES)].head(10).copy()
    top_delta["label"] = top_delta["gene"] + "·" + top_delta["cell_type"]
    ax_d.barh(top_delta["label"][::-1], top_delta["delta_tbi_minus_control"][::-1], color=OKABE["red"], edgecolor="black", linewidth=0.4)
    ax_d.axvline(0, color="black", linewidth=0.8)
    ax_d.set_xlabel("Donor-level mean log1p delta")
    ax_d.set_title("Core-gene TBI-control clues")

    ax_e = fig.add_subplot(gs[2, 0:2])
    dot = gene_summary.groupby(["cell_type", "gene"]).agg(
        mean_log1p_expr=("mean_log1p_expr", "mean"),
        pct_expressed=("pct_expressed", "mean"),
    ).reset_index()
    cell_order = [c for c in adata.obs["cell_type"].cat.categories if c != "Unknown"] + (["Unknown"] if "Unknown" in adata.obs["cell_type"].cat.categories else [])
    gene_order = [g for g in PRIORITY_GENES if g in dot["gene"].unique()]
    expr_mat = dot.pivot(index="cell_type", columns="gene", values="mean_log1p_expr").reindex(cell_order)[gene_order]
    pct_mat = dot.pivot(index="cell_type", columns="gene", values="pct_expressed").reindex(cell_order)[gene_order]
    for yi, cell_type in enumerate(cell_order):
        for xi, gene in enumerate(gene_order):
            expr_val = expr_mat.loc[cell_type, gene]
            pct_val = pct_mat.loc[cell_type, gene]
            if pd.isna(expr_val) or pd.isna(pct_val):
                continue
            size = max(12, pct_val * 2.2)
            color = expr_val
            ax_e.scatter(xi, yi, s=size, c=[color], cmap="YlOrRd", vmin=0, vmax=float(np.nanmax(expr_mat.values)), edgecolor="black", linewidth=0.25)
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
    present = [g for g in PRIORITY_GENES if g in module_df.columns]
    score = module_df[present].mean(axis=1)
    plot_df = pd.DataFrame(
        {
            "umap1": adata.obsm["X_umap"][:, 0],
            "umap2": adata.obsm["X_umap"][:, 1],
            "score_8gene_mean": score.reindex(adata.obs_names).to_numpy(),
        },
        index=adata.obs_names,
    )
    ax_f.scatter(plot_df["umap1"], plot_df["umap2"], c=plot_df["score_8gene_mean"], s=10, cmap="magma", linewidths=0)
    ax_f.set_title("Cell-level 8-gene localization score")
    ax_f.set_xlabel("UMAP1")
    ax_f.set_ylabel("UMAP2")
    sns.despine(ax=ax_f, left=False, bottom=False)

    for label, ax in zip(list("ABCDEF"), [ax_a, ax_b, ax_c, ax_d, ax_e, ax_f]):
        if ax is ax_a:
            ax.text(-0.01, 1.02, label, transform=ax.transAxes, fontsize=12, fontweight="bold", va="top", ha="left")
        else:
            ax.text(-0.12, 1.06, label, transform=ax.transAxes, fontsize=12, fontweight="bold", va="top", ha="left")

    out_png = FIG_DIR / "Fig3b_GSE209552_true_snRNA_localization_panel_20260613.png"
    out_pdf = FIG_DIR / "Fig3b_GSE209552_true_snRNA_localization_panel_20260613.pdf"
    fig.savefig(out_png, bbox_inches="tight", dpi=320)
    fig.savefig(out_pdf, bbox_inches="tight")
    plt.close(fig)


def make_feature_figure(adata: ad.AnnData) -> None:
    genes = [g for g in PRIORITY_GENES if g in adata.raw.var_names]
    raw_df = raw_gene_df(adata, genes)
    ncols = 4
    nrows = math.ceil(len(genes) / ncols)
    fig, axes = plt.subplots(nrows, ncols, figsize=(16, 3.8 * nrows))
    axes = np.atleast_1d(axes).reshape(nrows, ncols)
    for ax, gene in zip(axes.ravel(), genes):
        vals = raw_df[gene].reindex(adata.obs_names).to_numpy()
        sc_handle = ax.scatter(adata.obsm["X_umap"][:, 0], adata.obsm["X_umap"][:, 1], c=vals, s=8, cmap="viridis", linewidths=0)
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
    out_png = FIG_DIR / "Fig3c_GSE209552_true_snRNA_feature_maps_8gene_20260613.png"
    out_pdf = FIG_DIR / "Fig3c_GSE209552_true_snRNA_feature_maps_8gene_20260613.pdf"
    fig.savefig(out_png, bbox_inches="tight", dpi=320)
    fig.savefig(out_pdf, bbox_inches="tight")
    plt.close(fig)


def write_report(adata: ad.AnnData, cluster_summary: pd.DataFrame, delta_df: pd.DataFrame) -> None:
    cell_counts = adata.obs["cell_type"].value_counts().rename_axis("cell_type").reset_index(name="n_cells")
    top_lines = []
    for _, row in delta_df.head(8).iterrows():
        top_lines.append(
            f"- {row['gene']} in {row['cell_type']}: donor-level delta={row['delta_tbi_minus_control']:.3f} "
            f"(TBI mean {row['tbi_mean']:.3f} vs Control mean {row['control_mean']:.3f}; "
            f"TBI donors={int(row['tbi_donor_n'])}, Control donors={int(row['control_donor_n'])})"
        )
    report = "\n".join(
        [
            "# GSE209552真实单核定位分析补充说明",
            "",
            "本次补充分析直接使用 GSE209552 GEO RAW.tar 中的 10x 风格单核矩阵，重建真实 snRNA-seq 对象，而不是继续使用 bulk-like proxy。",
            "",
            "## 数据概况",
            f"- QC后保留 nuclei: {adata.n_obs}",
            f"- 聚类数: {adata.obs['leiden'].nunique()}",
            f"- 样本数: {adata.obs['sample_accession'].nunique()}",
            f"- donor数: {adata.obs['donor_id'].nunique()}",
            "",
            "## 主要细胞类型",
            *[f"- {row.cell_type}: {row.n_cells} nuclei" for row in cell_counts.itertuples(index=False)],
            "",
            "## 8基因单核层面最强定位线索",
            *top_lines,
            "",
            "## 方法边界",
            "- 本次结果属于真实单核定位和细胞类型表达层面的补强。",
            "- 若进一步比较 TBI 与 Control 的细胞类型差异，必须注意 control 侧存在同一 donor 的多个脑区样本，正式统计比较应以 donor-level/pseudobulk 为主。",
            "- 因此本次图谱中的条件差异更适合作为探索性定位线索，单细胞定位本身则可以作为真实证据使用。",
        ]
    )
    (REPORT_DIR / "GSE209552_true_snRNA_localization_report_20260613.md").write_text(report, encoding="utf-8")
    cluster_summary.to_csv(TABLE_DIR / "GSE209552_true_snRNA_cluster_annotation_table_20260613.csv", index=False, encoding="utf-8-sig")


def main() -> None:
    configure_style()
    adata = load_samples()
    adata = preprocess(adata)
    markers = cluster_marker_table(adata)
    adata, cluster_summary = annotate_cell_types(adata)
    gene_summary, donor_summary = summarize_gene_localization(adata)
    delta_df = donor_level_deltas(donor_summary)
    adata.write(WORKDIR / "raw_data_v3_20260604" / "GSE209552_true_snRNA_processed_20260613.h5ad")
    make_main_figure(adata, gene_summary, delta_df, cluster_summary)
    make_feature_figure(adata)
    write_report(adata, cluster_summary, delta_df)
    print(f"Processed nuclei: {adata.n_obs}")
    print(f"Clusters: {adata.obs['leiden'].nunique()}")
    print(f"Cell types: {adata.obs['cell_type'].cat.categories.tolist()}")


if __name__ == "__main__":
    main()
