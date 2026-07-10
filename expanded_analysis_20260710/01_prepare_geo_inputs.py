from __future__ import annotations

import argparse
import csv
import gzip
import io
import re
from pathlib import Path

import numpy as np
import pandas as pd


def parse_geo_series_matrix(path: Path) -> tuple[pd.DataFrame, dict[str, list[str]]]:
    metadata: dict[str, list[str]] = {}
    table_lines: list[str] = []
    reading_table = False
    with gzip.open(path, "rt", errors="replace") as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n")
            if line == "!series_matrix_table_begin":
                reading_table = True
                continue
            if line == "!series_matrix_table_end":
                break
            if reading_table:
                table_lines.append(line)
            elif line.startswith("!Sample_"):
                fields = [value.strip().strip('"') for value in line.split("\t")]
                metadata[fields[0]] = fields[1:]
    if not table_lines:
        raise ValueError(f"No series matrix table found in {path}")
    expression = pd.read_csv(io.StringIO("\n".join(table_lines)), sep="\t", low_memory=False)
    return expression, metadata


def extract_platform_gene_map(path: Path) -> pd.DataFrame:
    rows: list[dict[str, str]] = []
    header: list[str] | None = None
    reading_table = False
    with gzip.open(path, "rt", errors="replace", newline="") as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\r\n")
            if line == "!platform_table_begin":
                reading_table = True
                continue
            if line == "!platform_table_end":
                break
            if not reading_table:
                continue
            fields = next(csv.reader([line], delimiter="\t"))
            if header is None:
                header = fields
                continue
            if len(fields) < len(header):
                fields.extend([""] * (len(header) - len(fields)))
            record = dict(zip(header, fields))
            probe_id = record.get("ID", "").strip()
            symbol = record.get("GENE_SYMBOL", "").strip()
            if probe_id and symbol and symbol not in {"---", "NA"}:
                rows.append({"probe_id": probe_id, "gene_symbol": clean_symbol(symbol)})
    if not rows:
        raise ValueError(f"No probe-to-gene mapping found in {path}")
    return pd.DataFrame(rows).dropna().drop_duplicates()


def clean_symbol(value: object) -> str | None:
    if pd.isna(value):
        return None
    symbol = str(value).strip()
    if not symbol:
        return None
    symbol = re.split(r"///|//|;|,|\s+", symbol)[0].strip()
    return symbol.upper() or None


def collapse_probes_by_iqr(expression: pd.DataFrame, gene_map: pd.DataFrame) -> pd.DataFrame:
    expression = expression.rename(columns={expression.columns[0]: "probe_id"}).copy()
    sample_columns = [column for column in expression.columns if column != "probe_id"]
    expression[sample_columns] = expression[sample_columns].apply(pd.to_numeric, errors="coerce")
    merged = expression.merge(gene_map, on="probe_id", how="inner")
    merged["probe_iqr"] = merged[sample_columns].quantile(0.75, axis=1) - merged[sample_columns].quantile(0.25, axis=1)
    merged["probe_mean"] = merged[sample_columns].mean(axis=1)
    merged = merged.sort_values(["gene_symbol", "probe_iqr", "probe_mean"], ascending=[True, False, False])
    collapsed = merged.drop_duplicates("gene_symbol").set_index("gene_symbol")[sample_columns]
    collapsed.index = collapsed.index.astype(str).str.upper()
    return collapsed.sort_index()


def metadata_frame(metadata: dict[str, list[str]], sample_columns: list[str]) -> pd.DataFrame:
    accessions = metadata.get("!Sample_geo_accession", sample_columns)
    titles = metadata.get("!Sample_title", metadata.get("!Sample_source_name_ch1", sample_columns))
    if len(accessions) != len(sample_columns):
        accessions = sample_columns
    if len(titles) != len(sample_columns):
        titles = sample_columns
    return pd.DataFrame({"sample_id": sample_columns, "geo_accession": accessions, "title": titles})


def prepare_gse111452(matrix: Path, platform_soft: Path, output_dir: Path, platform: str) -> None:
    expression, metadata = parse_geo_series_matrix(matrix)
    expression = expression.rename(columns={expression.columns[0]: "probe_id"})
    sample_columns = [column for column in expression.columns if column != "probe_id"]
    gene_map = extract_platform_gene_map(platform_soft)
    collapsed = collapse_probes_by_iqr(expression, gene_map)
    sample_metadata = metadata_frame(metadata, sample_columns)

    def classify(title: str) -> pd.Series:
        region = "Hippocampus" if "hippocampus" in title.lower() else "Cortex"
        if "naive" in title.lower():
            condition = "Naive"
        elif "sham" in title.lower():
            condition = "Sham"
        else:
            condition = "TBI"
        lower = title.lower()
        if "24hr" in lower:
            time = "24h"
        elif "2wk" in lower:
            time = "2wk"
        elif "3mo" in lower:
            time = "3mo"
        elif "6mo" in lower:
            time = "6mo"
        elif "12mo" in lower or "1yr" in lower:
            time = "12mo"
        else:
            time = "unknown"
        return pd.Series({"region": region, "condition": condition, "time": time})

    sample_metadata = pd.concat([sample_metadata, sample_metadata["title"].apply(classify)], axis=1)
    sample_metadata["platform"] = platform
    gene_map.to_csv(output_dir / f"GSE111452_{platform}_probe_gene_map.csv", index=False)
    collapsed.to_csv(output_dir / f"GSE111452_{platform}_gene_expression.tsv.gz", sep="\t", compression="gzip")
    sample_metadata.to_csv(output_dir / f"GSE111452_{platform}_metadata.csv", index=False)


def prepare_gse223245(matrix: Path, output_dir: Path, family_soft: Path | None = None, gene_map_path: Path | None = None) -> None:
    expression, metadata = parse_geo_series_matrix(matrix)
    expression = expression.rename(columns={expression.columns[0]: "probe_id"})
    sample_columns = [column for column in expression.columns if column != "probe_id"]
    if gene_map_path is not None:
        gene_map = pd.read_csv(gene_map_path)
        if not {"probe_id", "gene_symbol"}.issubset(gene_map.columns):
            gene_map = gene_map.rename(columns={gene_map.columns[0]: "probe_id", gene_map.columns[1]: "gene_symbol"})
        gene_map = gene_map[["probe_id", "gene_symbol"]].copy()
        gene_map["gene_symbol"] = gene_map["gene_symbol"].map(clean_symbol)
        gene_map = gene_map.dropna().drop_duplicates()
    elif family_soft is not None:
        gene_map = extract_platform_gene_map(family_soft)
    else:
        raise ValueError("GSE223245 requires either --gse223245-map or --gse223245-soft")
    collapsed = collapse_probes_by_iqr(expression, gene_map)
    sample_metadata = metadata_frame(metadata, sample_columns)

    def severity(title: str) -> tuple[str, int]:
        lower = title.lower()
        for group, value in (("Severe", 3), ("Moderate", 2), ("Mild", 1)):
            if group.lower() in lower:
                return group, value
        return "Control", 0

    groups = sample_metadata["title"].map(severity)
    sample_metadata["group"] = [group for group, _ in groups]
    sample_metadata["severity"] = [value for _, value in groups]
    gene_map.to_csv(output_dir / "GSE223245_probe_gene_map.csv", index=False)
    collapsed.to_csv(output_dir / "GSE223245_gene_expression.tsv.gz", sep="\t", compression="gzip")
    sample_metadata.to_csv(output_dir / "GSE223245_metadata.csv", index=False)


def main() -> None:
    parser = argparse.ArgumentParser(description="Prepare GEO series matrices for the expanded TBI analysis.")
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--gse111452-gpl22740-matrix", type=Path)
    parser.add_argument("--gpl22740-soft", type=Path)
    parser.add_argument("--gse111452-gpl15084-matrix", type=Path)
    parser.add_argument("--gpl15084-soft", type=Path)
    parser.add_argument("--gse223245-matrix", type=Path)
    parser.add_argument("--gse223245-soft", type=Path)
    parser.add_argument("--gse223245-map", type=Path)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    if args.gse111452_gpl22740_matrix and args.gpl22740_soft:
        prepare_gse111452(args.gse111452_gpl22740_matrix, args.gpl22740_soft, args.output_dir, "GPL22740")
    if args.gse111452_gpl15084_matrix and args.gpl15084_soft:
        prepare_gse111452(args.gse111452_gpl15084_matrix, args.gpl15084_soft, args.output_dir, "GPL15084")
    if args.gse223245_matrix and (args.gse223245_soft or args.gse223245_map):
        prepare_gse223245(
            args.gse223245_matrix,
            args.output_dir,
            family_soft=args.gse223245_soft,
            gene_map_path=args.gse223245_map,
        )


if __name__ == "__main__":
    main()
