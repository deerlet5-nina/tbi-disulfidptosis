# Reproducible analysis package

## Scope

This package reproduces the expanded public-transcriptomic analysis for the manuscript titled "Transcriptomic Prioritization of a Disulfidptosis-Like Transporter-Actin Stress Axis After Traumatic Brain Injury". It fits each dataset within its native design and does not pool cross-platform, cross-species or cross-tissue differential effects.

## Included files

- `01_prepare_geo_inputs.py`: parses GEO series matrices and platform annotations for GSE111452 and GSE223245, maps probes to gene symbols and selects one probe per gene by across-sample interquartile range.
- `02_run_expanded_bulk_analysis.R`: runs gene-level, module-level, ranked-enrichment, coupling and interaction analyses for the bulk and microarray datasets.
- `03_run_donor_level_snrna_analysis.py`: aggregates GSE209552 nuclei within donor-by-cell-type strata and runs localization, exact-permutation condition effects and adjusted module-coupling analyses.
- `04_create_synthesis_tables_and_figures.py`: creates the four main figures, three supplementary figures and synthesis tables from saved analysis results.
- `mechanism_modules.tsv`: prespecified module-to-gene definitions.
- `R_sessionInfo.txt`: R version and package environment used for the submitted run.
- `assets/Fig1_base_design.tif`: base study-design and hypothesis artwork included in submitted Online Resource 3. The binary artwork is not duplicated in this GitHub folder; supply it locally when composing Fig. 1.

## Public inputs

Download the relevant count matrices, series matrices and platform annotations from NCBI GEO for GSE193407, GSE319253, GSE209552, GSE104687, GSE163415, GSE223245, GSE111452 and GSE298240. The donor-annotated GSE209552 AnnData object used for the submitted single-nucleus analysis is represented by the complete donor-level output tables in Online Resource 2; users re-creating the AnnData object from GEO should preserve donor, condition, sample-accession and cell-type annotations.

## Python input preparation

Run `01_prepare_geo_inputs.py` with GEO series-matrix and platform-SOFT paths. The script accepts separate GSE111452 GPL22740/GPL15084 inputs and either a GSE223245 SOFT annotation or a two-column probe-to-gene map. Output is written to a user-selected prepared-data directory.

Example:

```text
python 01_prepare_geo_inputs.py --output-dir prepared \
  --gse111452-gpl22740-matrix <series_matrix.txt.gz> --gpl22740-soft <GPL22740.soft.gz> \
  --gse111452-gpl15084-matrix <series_matrix.txt.gz> --gpl15084-soft <GPL15084.soft.gz> \
  --gse223245-matrix <series_matrix.txt.gz> --gse223245-soft <platform.soft.gz>
```

## Bulk and microarray analysis

Set the environment variables consumed by `02_run_expanded_bulk_analysis.R`. Required variables are `NCR_OUTPUT_DIR`, `NCR_MODULES`, `NCR_GSE193_COUNTS`, `NCR_GSE193_META`, `NCR_GSE319_COUNTS`, `NCR_GSE209_COUNTS`, `NCR_GSE163_3D_COUNTS`, `NCR_GSE163_29D_COUNTS`, `NCR_GSE298_COUNTS` and `NCR_GSE298_SERIES`. Optional variables provide the prepared GSE111452, GSE223245 and GSE104687 inputs. Then run:

```text
Rscript 02_run_expanded_bulk_analysis.R
```

The workflow uses edgeR filtering and TMM normalization, voom precision weights, robust limma empirical Bayes models, Hedges g module effects, Spearman severity tests and fgseaMultilevel ranked enrichment. Dataset-specific covariates are described in the manuscript and Supplementary Information.

## Donor-level single-nucleus analysis

```text
python 03_run_donor_level_snrna_analysis.py \
  --h5ad <GSE209552_processed.h5ad> \
  --modules mechanism_modules.tsv \
  --outdir results/snrna \
  --min-nuclei 20
```

Donors are the independent units. The script does not treat nuclei as biological replicates. Coupling tests residualize TBI condition, log nuclei count and a housekeeping-expression index and remove shared genes from target modules where necessary.

## Figure and synthesis-table generation

```text
python 04_create_synthesis_tables_and_figures.py \
  --bulk-dir results/bulk \
  --snrna-dir results/snrna \
  --h5ad <GSE209552_processed.h5ad> \
  --base-fig1 assets/Fig1_base_design.tif \
  --figure-dir figures \
  --table-dir results/synthesis
```

Figures are written as PDF, PNG and 600-dpi LZW-compressed TIFF. The manuscript uses the TIFF files as separate upload items and PNG files inside the review DOCX.

## Statistical interpretation

Benjamini-Hochberg FDR below 0.05 is considered supported within the relevant analysis family. Unadjusted P below 0.05 without FDR support is nominal. Cross-dataset heatmaps are visual summaries and are not pooled statistical estimates. Single-nucleus coupling is non-causal. None of the analyses measures cystine flux, NADPH debt, protein disulfidation, F-actin collapse or rescue of a specific cell-death phenotype.
