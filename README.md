# tbi-disulfidptosis

This repository contains processed source tables and analysis scripts for the manuscript **"Transcriptomic Prioritization of a Disulfidptosis-Like Transporter-Actin Stress Axis After Traumatic Brain Injury"** submitted to *Neurochemical Research*.

The study is a public-transcriptome prioritization analysis. It does not contain newly generated human participant data, human biospecimens, animal experiments or wet-laboratory measurements. The repository is intended to make the processed tables, figure-supporting files and analysis scripts available for editorial and reviewer inspection.

## Repository contents

- `data/source_tables/`: CSV source tables corresponding to the manuscript online source-table archive.
- `scripts/`: analysis and manuscript/figure construction scripts supplied with the submission package.
- `figures/`: supplementary robustness-audit figure supplied as an online resource.
- `requirements.txt`: core Python packages used across the workflow.

## Data sources

The analyses use public Gene Expression Omnibus resources cited in the manuscript: GSE209552, GSE193407, GSE319253, GSE104687, GSE163415 and GSE223245. Some scripts may require users to download GEO data files separately before rerunning the complete workflow.

## Interpretation boundary

The integrated score and robustness audit are prioritization tools for heterogeneous public transcriptomic evidence. They are not pooled statistical effect estimates, diagnostic models or direct evidence of disulfidptotic cell death. Protein, redox, F-actin and cell-type validation remain required for mechanistic confirmation.

## Correspondence

Corresponding author: Mingyang Zhang, mingyangzhang@suda.edu.cn
