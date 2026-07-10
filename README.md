# tbi-disulfidptosis

This public repository supports the manuscript **"Transcriptomic Prioritization of a Disulfidptosis-Like Transporter-Actin Stress Axis After Traumatic Brain Injury"** submitted to *Neurochemical Research*.

The study reanalyzes public transcriptomic data and does not contain newly generated human-participant data, human biospecimens, animal experiments, or wet-laboratory measurements.

## Current expanded analysis

The submission-version workflow is available in [`expanded_analysis_20260710/`](expanded_analysis_20260710/). It includes:

- preparation of GSE111452 and GSE223245 GEO matrices;
- covariate-aware bulk and microarray models;
- prespecified 15-module scoring and ranked enrichment;
- donor-level GSE209552 single-nucleus localization, permutation testing, and adjusted module coupling;
- four compact main figures, three supplementary figures, and synthesis tables.

The complete processed source tables and the Figure 1 base artwork are distributed with the manuscript as Online Resources 2 and 3. The repository code is intended to be used with the public GEO inputs documented in the workflow README.

## Public data sources

GSE209552, GSE193407, GSE319253, GSE104687, GSE163415, GSE223245, GSE111452, and GSE298240.

## Interpretation boundary

Results are kept within their native dataset designs. Cross-dataset displays are descriptive and are not pooled statistical estimates. Donors, not nuclei, are the inferential units for single-nucleus comparisons. The analyses do not measure cystine flux, NADPH debt, protein disulfidation, F-actin collapse, or rescue of a specific cell-death phenotype; they prioritize experiments rather than demonstrate disulfidptotic cell death.

## Correspondence

Mingyang Zhang, mingyangzhang@suda.edu.cn
