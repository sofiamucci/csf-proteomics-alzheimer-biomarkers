# Methods

## Pipeline Overview

This repository reanalyzes the cerebrospinal fluid (CSF) proteome of
Alzheimer's disease (AD) patients and controls from PRIDE PXD016278
(Bader et al., 2020, *Mol Syst Biol* 16:e9356; three cohorts, 197
individuals). Samples were acquired by data-independent acquisition (DIA) and
processed by the authors with Spectronaut; the analysis starts from their
protein group report rather than from raw files. Downstream steps (filtering,
normalization, missing-value handling, differential abundance with limma and
a diagnostic classifier) are documented in the sections below as they are
completed.

## Data Source and Acquisition

Three files are downloaded from PRIDE by `scripts/00_download_data.sh`: the
protein group report used for the published analysis, the authors' sample
annotation, and the community-annotated SDRF file. Raw files, spectral
libraries and Spectronaut sessions are not needed. Input data are not
committed to the repository; the script tests zip integrity and records the
date, size and SHA256 checksum of each file in `data/raw/DOWNLOAD_LOG.txt`,
so the exact inputs of every run are traceable.

Because the data come from Spectronaut rather than MaxQuant, there are no
`Reverse`, `Potential contaminant` or `Only identified by site` columns:
decoys are removed upstream by the q-value cutoff, and contaminants are
filtered explicitly in this pipeline. Values that did not pass the q-value
cutoff are reported as `"Filtered"` and are read as missing.

## Data Import

`scripts/01_import.R` converts the Spectronaut report into a protein × run
intensity matrix (1,542 protein groups × 210 runs, untransformed) and a
protein annotation table, saved together as `data/processed/01_import.rds`.
Protein groups are identified by `PG.ProteinGroups`, which is unique;
`PG.Genes` is kept as annotation only, since it is missing for three groups
and shared by others. Run names are the raw-file names without Spectronaut's
index prefix and quantity suffix. No filtering, transformation or imputation
is applied at this step; the script stops if dimensions, identifier
uniqueness or positivity of the intensities differ from expected. 21.2% of
values are missing.

## Software Environment

R 4.6.1 is installed in a dedicated conda environment (`environment.yml`,
exported with `--from-history`, containing R only); R packages are pinned
with renv (`renv.lock`). At the time of setup conda-forge had no `r-renv`
build for R 4.6, so renv is installed from CRAN through the project's
bootstrap script (`renv/activate.R`). A Docker image restoring `renv.lock`
on `rocker/r-ver:4.6.1` will be added once the package set is final.

