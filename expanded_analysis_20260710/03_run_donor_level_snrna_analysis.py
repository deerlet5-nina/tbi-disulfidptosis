#!/usr/bin/env python3
"""Donor-level analysis of the GSE209552 acute TBI snRNA-seq dataset.

The analysis deliberately uses donors, rather than nuclei, as independent
replicates. The processed AnnData object stores log-normalized expression in
``raw.X``. We aggregate these values within donor-by-cell-type strata, score
prespecified modules, estimate condition effects with exact label
permutations, and test transporter-to-cytoskeleton coupling after residualizing
condition and the number of nuclei.
"""

from __future__ import annotations

import argparse
import itertools
import math
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse, stats
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests


CORE_MODULES = [
    "Eight_gene_panel",
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
    "BBB_endothelial_context",
    "Inflammation_context",
]

HOUSEKEEPING_GENES = [
    "RPLP0",
    "RPL13A",
    "RPL32",
    "RPS18",
    "RPS27A",
    "TBP",
    "PPIA",
    "HPRT1",
    "GUSB",
    "YWHAZ",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--h5ad", required=True, type=Path)
    parser.add_argument("--modules", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--min-nuclei", type=int, default=20)
    return parser.parse_args()


def hedges_g(case: np.ndarray, control: np.ndarray) -> float:
    case = np.asarray(case, dtype=float)
    control = np.asarray(control, dtype=float)
    n1, n0 = len(case), len(control)
    if n1 < 2 or n0 < 2:
        return np.nan
    pooled_num = (n1 - 1) * np.var(case, ddof=1) + (n0 - 1) * np.var(control, ddof=1)
    pooled_den = n1 + n0 - 2
    if pooled_den <= 0 or pooled_num <= 0:
        return 0.0 if np.isclose(np.mean(case), np.mean(control)) else np.nan
    d = (np.mean(case) - np.mean(control)) / math.sqrt(pooled_num / pooled_den)
    correction = 1 - 3 / (4 * (n1 + n0) - 9)
    return float(correction * d)


def exact_permutation_p(values: np.ndarray, labels: np.ndarray) -> float:
    """Two-sided exact test for the difference in means.

    With 3 controls and 12 TBI donors, only 455 allocations are required.
    For larger allocation spaces, a deterministic 20,000-label sample is used.
    """

    values = np.asarray(values, dtype=float)
    labels = np.asarray(labels, dtype=bool)
    n = len(values)
    n_case = int(labels.sum())
    observed = abs(values[labels].mean() - values[~labels].mean())
    n_allocations = math.comb(n, n_case)
    if n_allocations <= 50_000:
        exceed = 0
        total = 0
        all_idx = np.arange(n)
        for case_idx in itertools.combinations(range(n), n_case):
            mask = np.zeros(n, dtype=bool)
            mask[list(case_idx)] = True
            statistic = abs(values[mask].mean() - values[~mask].mean())
            exceed += statistic >= observed - 1e-12
            total += 1
        return float(exceed / total)

    rng = np.random.default_rng(20260710)
    exceed = 1
    total = 20_001
    all_idx = np.arange(n)
    for _ in range(total - 1):
        mask = np.zeros(n, dtype=bool)
        mask[rng.choice(all_idx, size=n_case, replace=False)] = True
        statistic = abs(values[mask].mean() - values[~mask].mean())
        exceed += statistic >= observed - 1e-12
    return float(exceed / total)


def residual_spearman(frame: pd.DataFrame, x: str, y: str) -> tuple[float, float, int]:
    work = frame[[x, y, "condition", "n_nuclei", "technical_expression_index"]].dropna().copy()
    if len(work) < 6 or work[x].nunique() < 3 or work[y].nunique() < 3:
        return np.nan, np.nan, len(work)
    covariates = pd.DataFrame(
        {
            "TBI": (work["condition"] == "TBI").astype(float),
            "log_n_nuclei": np.log1p(work["n_nuclei"].astype(float)),
            "technical_expression_index": work["technical_expression_index"].astype(float),
        },
        index=work.index,
    )
    design = sm.add_constant(covariates, has_constant="add")
    residual_x = sm.OLS(work[x].astype(float), design).fit().resid
    residual_y = sm.OLS(work[y].astype(float), design).fit().resid
    result = stats.spearmanr(residual_x, residual_y)
    return float(result.statistic), float(result.pvalue), len(work)


def main() -> None:
    args = parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    module_table = pd.read_csv(args.modules, sep="\t")
    module_table = module_table[module_table["module"].isin(CORE_MODULES)].copy()

    adata = ad.read_h5ad(args.h5ad, backed="r")
    required_obs = {"donor_id", "condition", "cell_type", "sample_accession"}
    missing_obs = required_obs.difference(adata.obs.columns)
    if missing_obs:
        raise ValueError(f"Missing AnnData annotations: {sorted(missing_obs)}")
    if adata.raw is None:
        raise ValueError("Expected log-normalized full-gene expression in adata.raw")

    module_genes = module_table["gene_symbol"].drop_duplicates().tolist()
    requested_genes = list(dict.fromkeys(module_genes + HOUSEKEEPING_GENES))
    present_genes = [gene for gene in requested_genes if gene in adata.raw.var_names]
    missing_module_genes = sorted(set(module_genes).difference(present_genes))
    missing_housekeeping_genes = sorted(set(HOUSEKEEPING_GENES).difference(present_genes))
    gene_indices = adata.raw.var_names.get_indexer(present_genes)
    order = np.argsort(gene_indices)
    gene_indices = gene_indices[order]
    present_genes = [present_genes[i] for i in order]
    expression = adata.raw.X[:, gene_indices]
    if not sparse.issparse(expression):
        expression = sparse.csr_matrix(expression)
    else:
        expression = expression.tocsr()

    obs = adata.obs[["donor_id", "condition", "cell_type", "sample_accession"]].copy()
    obs["donor_id"] = obs["donor_id"].astype(str)
    obs["condition"] = obs["condition"].astype(str)
    obs["cell_type"] = obs["cell_type"].astype(str)
    obs["donor_key"] = obs["condition"] + "_" + obs["donor_id"]

    records: list[dict[str, object]] = []
    grouped = obs.groupby(["donor_key", "donor_id", "condition", "cell_type"], sort=True)
    for (donor_key, donor_id, condition, cell_type), indices in grouped.indices.items():
        row_indices = np.asarray(indices, dtype=int)
        n_nuclei = len(row_indices)
        means = np.asarray(expression[row_indices].mean(axis=0)).ravel()
        detected = np.asarray((expression[row_indices] > 0).mean(axis=0)).ravel()
        record: dict[str, object] = {
            "donor_key": donor_key,
            "donor_id": donor_id,
            "condition": condition,
            "cell_type": cell_type,
            "n_nuclei": n_nuclei,
        }
        record.update({f"expr_{gene}": value for gene, value in zip(present_genes, means)})
        record.update({f"detected_{gene}": value for gene, value in zip(present_genes, detected)})
        records.append(record)
    pseudobulk = pd.DataFrame.from_records(records)
    pseudobulk.to_csv(args.outdir / "Table_S_snRNA_donor_celltype_gene_expression.csv.gz", index=False)

    eligible = pseudobulk[pseudobulk["n_nuclei"] >= args.min_nuclei].copy()
    expression_columns = [f"expr_{gene}" for gene in present_genes]
    gene_matrix = eligible[expression_columns].copy()
    gene_matrix.columns = present_genes
    gene_sd = gene_matrix.std(axis=0, ddof=0).replace(0, np.nan)
    gene_z = (gene_matrix - gene_matrix.mean(axis=0)) / gene_sd

    scores = eligible[["donor_key", "donor_id", "condition", "cell_type", "n_nuclei"]].copy()
    housekeeping_present = [gene for gene in HOUSEKEEPING_GENES if gene in gene_matrix.columns]
    scores["technical_expression_index"] = gene_matrix[housekeeping_present].mean(axis=1)
    coverage_records = []
    for module_name, subtable in module_table.groupby("module", sort=False):
        genes = subtable["gene_symbol"].drop_duplicates().tolist()
        covered = [gene for gene in genes if gene in gene_z.columns]
        scores[module_name] = gene_z[covered].mean(axis=1) if covered else np.nan
        coverage_records.append(
            {
                "module": module_name,
                "n_defined": len(genes),
                "n_covered": len(covered),
                "coverage_fraction": len(covered) / len(genes) if genes else np.nan,
                "covered_genes": ";".join(covered),
                "missing_genes": ";".join(sorted(set(genes).difference(covered))),
            }
        )
    scores.to_csv(args.outdir / "Table_snRNA_donor_celltype_module_scores.csv", index=False)
    pd.DataFrame(coverage_records).to_csv(args.outdir / "Table_snRNA_module_coverage.csv", index=False)

    effect_records = []
    for cell_type, cell_frame in scores.groupby("cell_type", sort=True):
        for module_name in CORE_MODULES:
            if module_name not in cell_frame:
                continue
            work = cell_frame[["condition", module_name]].dropna()
            case = work.loc[work["condition"] == "TBI", module_name].to_numpy(float)
            control = work.loc[work["condition"] == "Control", module_name].to_numpy(float)
            if len(case) < 3 or len(control) < 3:
                continue
            labels = (work["condition"].to_numpy() == "TBI")
            effect_records.append(
                {
                    "cell_type": cell_type,
                    "module": module_name,
                    "n_TBI": len(case),
                    "n_control": len(control),
                    "mean_TBI": float(case.mean()),
                    "mean_control": float(control.mean()),
                    "mean_difference": float(case.mean() - control.mean()),
                    "hedges_g": hedges_g(case, control),
                    "permutation_p": exact_permutation_p(work[module_name].to_numpy(float), labels),
                }
            )
    effects = pd.DataFrame(effect_records)
    if not effects.empty:
        effects["FDR"] = multipletests(effects["permutation_p"], method="fdr_bh")[1]
    effects.to_csv(args.outdir / "Table_snRNA_donor_level_condition_effects.csv", index=False)

    coupling_records = []
    coupling_pairs = [
        ("Transporter_entry", "Regulatory_bridge"),
        ("Transporter_entry", "Cytoskeletal_endpoint"),
        ("PPP_NADPH", "Disulfidptosis_core"),
        ("Glutathione_redox", "Disulfidptosis_core"),
        ("WRC_actin_branching", "Actomyosin_endpoint"),
        ("Focal_adhesion_tension", "Actomyosin_endpoint"),
    ]
    coupling_scores = scores.copy()
    module_gene_lookup = {
        module: set(subtable["gene_symbol"]).intersection(gene_z.columns)
        for module, subtable in module_table.groupby("module", sort=False)
    }
    target_columns: dict[tuple[str, str], tuple[str, str]] = {}
    for source, target in coupling_pairs:
        overlap = module_gene_lookup.get(source, set()).intersection(module_gene_lookup.get(target, set()))
        if overlap:
            disjoint_target = sorted(module_gene_lookup[target].difference(overlap))
            target_column = f"{target}__disjoint_from__{source}"
            coupling_scores[target_column] = gene_z[disjoint_target].mean(axis=1)
        else:
            target_column = target
        target_columns[(source, target)] = (target_column, ";".join(sorted(overlap)))
    for cell_type, cell_frame in scores.groupby("cell_type", sort=True):
        for source, target in coupling_pairs:
            target_column, overlap_removed = target_columns[(source, target)]
            cell_frame = coupling_scores.loc[cell_frame.index]
            if source not in cell_frame or target_column not in cell_frame:
                continue
            rho, pvalue, n = residual_spearman(cell_frame, source, target_column)
            coupling_records.append(
                {
                    "cell_type": cell_type,
                    "source_module": source,
                    "target_module": target,
                    "partial_spearman_rho": rho,
                    "p": pvalue,
                    "n_donors": n,
                    "overlap_removed_from_target": overlap_removed,
                }
            )
    coupling = pd.DataFrame(coupling_records)
    if not coupling.empty:
        valid = coupling["p"].notna()
        coupling.loc[valid, "FDR"] = multipletests(coupling.loc[valid, "p"], method="fdr_bh")[1]
    coupling.to_csv(args.outdir / "Table_snRNA_module_coupling.csv", index=False)

    # Cell fractions are retained as an exploratory sensitivity analysis because
    # two control donors contributed both frontal and temporal specimens.
    counts = (
        obs.groupby(["donor_key", "donor_id", "condition", "cell_type"], observed=True)
        .size()
        .rename("n_nuclei")
        .reset_index()
    )
    totals = counts.groupby("donor_key")["n_nuclei"].transform("sum")
    counts["fraction"] = counts["n_nuclei"] / totals
    counts.to_csv(args.outdir / "Table_S_snRNA_donor_cell_fractions.csv", index=False)

    localization_records = []
    for cell_type, cell_frame in pseudobulk.groupby("cell_type", sort=True):
        for gene in ["SLC3A2", "SLC7A11", "WASF2", "TLN1", "ACTB", "MYH9", "MYL6", "FLNA"]:
            expr_col = f"expr_{gene}"
            detect_col = f"detected_{gene}"
            if expr_col not in cell_frame:
                continue
            localization_records.append(
                {
                    "cell_type": cell_type,
                    "gene_symbol": gene,
                    "n_donor_celltype_strata": len(cell_frame),
                    "median_donor_mean_expression": float(cell_frame[expr_col].median()),
                    "mean_detection_fraction": float(cell_frame[detect_col].mean()),
                }
            )
    pd.DataFrame(localization_records).to_csv(args.outdir / "Table_snRNA_8gene_localization.csv", index=False)

    inventory = pd.DataFrame(
        {
            "dataset": ["GSE209552"],
            "n_nuclei": [adata.n_obs],
            "n_donors": [obs["donor_key"].nunique()],
            "n_TBI_donors": [obs.loc[obs["condition"] == "TBI", "donor_key"].nunique()],
            "n_control_donors": [obs.loc[obs["condition"] == "Control", "donor_key"].nunique()],
            "n_cell_types": [obs["cell_type"].nunique()],
            "min_nuclei_per_tested_stratum": [args.min_nuclei],
            "missing_module_genes": [";".join(missing_module_genes)],
            "missing_housekeeping_genes": [";".join(missing_housekeeping_genes)],
            "note": [
                "Donors are statistical units; condition effects are exploratory because only three control donors are available."
            ],
        }
    )
    inventory.to_csv(args.outdir / "Table_snRNA_dataset_inventory.csv", index=False)


if __name__ == "__main__":
    main()
