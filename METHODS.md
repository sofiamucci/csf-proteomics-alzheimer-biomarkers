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

## Protein Filtering

`scripts/04_filter.R` removes 14 laboratory-contaminant protein groups
(12 skin and hair keratin groups, and trypsin PRSS1/PRSS3, which shares
peptides with the digestion enzyme). Lysozyme, a common laboratory contaminant
but also a normal CSF protein, is retained, as are the blood-panel proteins
used as QC markers.

Protein groups are retained if quantified in at least 50% of the runs of at
least one diagnostic group, leaving 1,317 of 1,542. Requiring the threshold in
both groups would have removed 75 further protein groups, 74 of them quantified
mainly in AD; these include NEFL, UCHL1, PTPA and EFHD2, i.e. intracellular
neuronal proteins that rise above the detection limit in AD. The number of
retained protein groups changes gradually with the threshold (1,401 at 30%,
1,317 at 50%, 1,153 at 70%), so the choice of 50% is not at an unstable point.
All protein groups are kept in the output object with their exclusion reason.

## Missing Values

`scripts/05_missing_value_methods.R` compares three treatments of missing
values on the same model (`~ group + batch + age + sex`, AD vs. control):
limma on observed values only, limma after left-censored imputation
(down-shifted normal; shift 1.8 SD, width 0.3 SD, fixed seed), and proDA
(Ahlmann-Eltze & Anders, bioRxiv 2019), which models the probability of a missing
value as a function of intensity without imputation. The comparison is shown
in `reports/analysis_report.qmd`.

At FDR < 5%, proDA identifies 244 protein groups (178 higher, 66 lower in AD),
limma on observed values 259, and imputation 290. Imputation inflated
fold changes across the proteome, added 47 protein groups higher in AD, and
disagreed with proDA on 132 protein groups. Limma on observed values agreed
with proDA on 239 of 244; the 25 discordant protein groups lie near the
threshold with near-identical fold changes. proDA is used for the main
analysis. Because proDA does not accept sample weights, limma on observed
values with `arrayWeights` is reported as a sensitivity analysis.

## Differential Abundance

`scripts/06_differential_abundance.R` fits proDA with
`~ group + site + age + sex` on the 1,317 retained protein groups and the 197
retained runs, testing AD vs. control with Benjamini-Hochberg FDR < 5%. Site is
used instead of MS batch: each site was measured in a single batch, and a site
term separates Kiel (controls only) from Magdeburg, so Kiel controls inform
variance and age effects but not the AD contrast.

The main model identifies 259 protein groups (192 higher, 67 lower in AD). The
original report's positive controls are recovered: SOD1 (log2FC 0.34) and PARK7
(0.25), with MAPT the strongest effect (1.11; partly circular, as diagnosis is
defined from CSF tau). UCHL1 is higher in AD (0.58); NEFL is not significant,
consistent with its quantification near the detection limit by MS.

Protein groups are flagged as blood- or plasma-associated if they belong to the
corresponding panel or correlate with its score across runs (Spearman
r > 0.5).

Four sensitivity analyses were compared with the main model:

| Analysis | r (log2FC) | Significant | Shared | Lost | Gained |
|---|---|---|---|---|---|
| limma, observed values, `arrayWeights` | 0.959 | 313 | 246 | 13 | 67 |
| Without 4 blood-flagged runs | 0.999 | 261 | 252 | 7 | 9 |
| Plasma-protein score as covariate | 0.842 | 53 | 51 | 208 | 2 |
| Per-run median renormalization | 0.998 | 251 | 209 | 50 | 42 |

Results are robust to blood contamination, normalization and the
missing-value model. The plasma-protein score is lower in AD (−0.34, p = 6 × 10⁻⁴,
adjusted for site, age and sex), and adjusting for it leaves 51 of the 259
protein groups significant. Without a CSF/serum albumin quotient it cannot be
determined whether the score is a confounder (barrier permeability) or lies on
the disease pathway, so results are reported in two tiers: the main model, and
a barrier-robust core significant in both models (45 higher, 6 lower in AD).
Of the 67 protein groups lower in AD, 61 lose significance after the
adjustment, so decreases are largely explained by the plasma-protein balance.
The core is dominated by glycolytic enzymes, 14-3-3 proteins and established
CSF markers of AD.

`scripts/07_export_results.R` writes the full results table
(`results/differential_abundance.csv`) with a column dictionary.

## Diagnostic Classifier

`scripts/08_classifier.R` fits an elastic-net logistic regression (glmnet,
α = 0.5, `lambda.1se` from inner 5-fold cross-validation) to classify AD vs.
control. MAPT is excluded from the predictors because diagnosis is defined from
CSF tau; Kiel is excluded because its controls are age-mismatched and it has no
AD cases, leaving 181 samples (88 AD). Candidate predictors are the 981 protein
groups with ≤ 20% missing values. Imputation (half the training minimum),
scaling, protein selection and penalty are learned within each training split,
so no information from test samples enters model building. AUC is computed
with a fixed direction (controls lower).

Repeated nested cross-validation (10 × 5-fold, stratified) gives a mean AUC of
0.88 (range 0.84–0.90), against 0.60 for an age + sex baseline. In
leave-one-site-out validation the AUC is 0.86 for Sweden, 0.88 for Magdeburg
and 0.74 for Berlin; Berlin is the largest cohort, so the model is trained on
fewer samples when it is held out, and its biochemical controls are mostly
patients with other clinical diagnoses. The most stably selected predictors
(YWHAZ, ALDOA and SPP1 in all 50 outer-fold models; YWHAG, GDA, PKM and SMOC1 in
≥ 90%) coincide with the barrier-robust core of the differential abundance
analysis.

## Software Environment

R 4.6.1 is installed in a dedicated conda environment (`environment.yml`,
exported with `--from-history`, containing R only); R packages are pinned
with renv (`renv.lock`). At the time of setup conda-forge had no `r-renv`
build for R 4.6, so renv is installed from CRAN through the project's
bootstrap script (`renv/activate.R`). A Docker image restoring `renv.lock`
on `rocker/r-ver:4.6.1` will be added once the package set is final.
