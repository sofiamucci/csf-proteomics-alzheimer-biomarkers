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
filtered explicitly in this pipeline. The export applied a decoy filter and a
q-value filter (≤ 0.01 per run); values that did not pass the cutoff are
reported as `"Filtered"` and are read as missing. No imputation was performed
in Spectronaut. Intensities were already normalized by the authors with
Spectronaut's "Local Normalization" (local regression; Callister et al.,
2006), applied once across all 210 runs.

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

## Sample Annotation

`scripts/02_annotate.R` joins the authors' sample table (Table S2) to the
intensity matrix by raw-file name; all 210 runs are annotated. The authors
excluded 13 runs (comorbidity, 6; excessive missing values, 5; blood–brain
barrier leakage, 1; duplicate injection, 1). These runs are retained in the
data object with their exclusion reason, for inspection in QC, and removed
before analysis, leaving the published 197 individuals.

The four collection sites were measured in three MS batches: Sweden (n = 60),
Magdeburg (38) and Kiel (16) together, and Berlin (83). Site therefore
determines batch. Kiel contributes only controls, younger than those of the
other sites (median age 32 vs. 68–74 years).

Diagnosis is the biochemical AD classification (88 AD, 109 controls), which is
complete for all retained samples; the clinical diagnosis is missing for 60.
Because the biochemical classification is derived from CSF tau and Aβ, the
association of MAPT with diagnosis is partly circular and is used as a
pipeline check rather than as a finding.

## Sample Quality Control

`scripts/03_sample_qc.R` computes per-run metrics, contamination scores and
flags on the 197 retained runs; `reports/qc_report.qmd` displays and models
them. Comparisons between runs use the 385 protein groups quantified in every
run, because medians of observed values are inflated in runs with more missing
values (r = 0.85 between missing fraction and observed median, vs. −0.12 on
complete protein groups).

**Normalization.** Spectronaut's Local Normalization was retained. Residual
batch offsets on complete protein groups (−0.18 and −0.23 log2 relative to the
Swedish batch) are modelled as a covariate in limma rather than removed by
renormalization, because per-run medians also showed a small association with
diagnosis (+0.03 log2, 95% CI −0.004 to 0.068) that run-level centering would
remove. A renormalized analysis is reported as a sensitivity check.

**Missing values.** The fraction of missing values per run is lower in AD than
in controls (−3.5 percentage points, p = 4 × 10⁻⁵) after adjusting for batch and
age. Missingness is therefore not independent of the contrast of interest,
which rules out both left-censored imputation and analysis of observed values
only.

**Low coverage.** Five runs quantify fewer protein groups than their batch
median minus 3 MAD. Three are young Kiel controls (20–47 years), measured in a
batch dominated by older Magdeburg samples; none show reduced intensity of
abundant proteins. All are retained, with sample quality weights estimated in
limma (`arrayWeights`).

**Blood contamination.** An erythrocyte score (mean z-score of nine
erythrocyte-specific proteins; Geyer et al., 2019) does not differ by
diagnosis (p = 0.89). Four runs exceed median + 3 MAD, all controls; they are
retained in the main analysis and removed in a sensitivity analysis. Protein
groups correlated with the score are flagged in the results.

**Global structure.** PCA on the complete protein groups (unit variance) shows
batch dominating PC2–PC4. PC1 (24% of variance) contrasts neuronal and neural
matrix proteins with liver- and immune-derived plasma proteins, consistent with
blood–CSF barrier permeability; it is partly associated with diagnosis
(R² = 0.11). The runs most distant in PCA are plasma-rich with normal
erythrocyte scores and are retained. PC1 is not used as a covariate, as it
carries disease signal; a plasma-protein score is used in a sensitivity
analysis instead.

## Software Environment

R 4.6.1 is installed in a dedicated conda environment (`environment.yml`,
exported with `--from-history`, containing R only); R packages are pinned
with renv (`renv.lock`). At the time of setup conda-forge had no `r-renv`
build for R 4.6, so renv is installed from CRAN through the project's
bootstrap script (`renv/activate.R`). A Docker image restoring `renv.lock`
on `rocker/r-ver:4.6.1` will be added once the package set is final.
