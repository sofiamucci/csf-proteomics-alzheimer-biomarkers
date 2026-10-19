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
committed to the repository; the script tests zip integrity, records the date,
size and SHA256 checksum of each file in `data/raw/DOWNLOAD_LOG.txt`, and stops
if the protein report or the sample annotation differ from the checksums of
the analysed files (`data/PXD016278.sha256`). The SDRF is community-maintained
and may be revised, so its checksum is logged but not enforced.

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
protein groups significant. The score admits three readings that these data
cannot separate: a confounder (blood–CSF barrier permeability; no CSF/serum
albumin quotient is available), a mediator on the disease pathway, or a
compositional effect of normalization, in which a rise of many brain-derived
proteins lowers the relative share of plasma proteins. Consistent with the
last, 67% of all protein groups correlate negatively with the score, 143 of the
192 protein groups higher in AD have r < −0.3, the adjustment shrinks their
fold changes by a median of 42%, and the score correlates with the core score
defined below (ρ = −0.66). The per-run median renormalization does not test
this, as it is also composition-based. Results are therefore reported in two
tiers: the main model, and a barrier-robust core significant in both models,
read as a conservative lower bound rather than as the specific signal. Of the
67 protein groups lower in AD, 61 lose significance after the adjustment.

The barrier-robust tier contains 51 protein groups; MAPT is flagged as circular
and excluded from the core used downstream, leaving 50 (44 higher, 6 lower in
AD). The core is dominated by glycolytic enzymes, 14-3-3 proteins and
established CSF markers of AD, mostly cytosolic proteins whose presence in CSF
most likely reflects release from damaged cells.

`scripts/07_export_results.R` writes the full results table
(`results/differential_abundance.csv`) with a column dictionary.

## Diagnostic Classifier

`scripts/08_classifier.R` fits an elastic-net logistic regression (glmnet,
α = 0.5, `lambda.1se` from inner 5-fold cross-validation) to classify AD vs.
control. MAPT is excluded from the predictors because diagnosis is defined from
CSF tau; Kiel is excluded because its controls are age-mismatched and it has no
AD cases, leaving 181 samples (88 AD). The missing-value filter (≤ 20% missing
in the training split; median 978 of 1,316 protein groups retained),
imputation (half the training minimum), scaling, protein selection and penalty
are learned within each training split, so no information from test samples
enters model building. AUC is computed with a fixed direction (controls
lower).

Repeated nested cross-validation (10 × 5-fold, stratified) gives a mean AUC of
0.88 (range 0.83–0.90 across repetitions), against 0.60 for an age + sex
baseline. Averaging out-of-fold predictions over repetitions gives one
prediction per sample, with an AUC of 0.90 (DeLong 95% CI 0.85–0.95; paired
test against the baseline, p = 1 × 10⁻¹⁰); at a threshold of 0.5, sensitivity
is 0.81 and specificity 0.83. Predictions are under-confident (calibration
slope 2.67), as expected from the strong penalty of `lambda.1se`: the model
ranks samples well, but its probabilities are not calibrated risks.

In leave-one-site-out validation the AUC is 0.88 for Magdeburg, 0.74 for Berlin
and 0.73 for Sweden. An earlier version applied the missing-value filter to all
samples, test site included, and gave 0.86 for Sweden: of the 1,088 protein
groups that pass the filter without Sweden, 226 have more than 20% missing
values in Sweden (100 more than 50%), where they are imputed from the training
minimum. Cross-site performance therefore depends on how well the new site
quantifies the model's proteins. The most stably selected predictors (YWHAZ,
ALDOA and SPP1 in all 50 outer-fold models; YWHAG, GDA, PKM in ≥ 90%; SMOC1 in
88%) coincide with the core of the differential abundance analysis.

## Circularity of the Diagnostic Label

The biochemical classification requires high CSF tau and low amyloid; the
amyloid measure and thresholds differ by site (in Sweden, Aβ42 alone separates
the groups; in Magdeburg and Berlin it does not, consistent with use of the
Aβ42/40 ratio). Every AD case has t-tau above 400 pg/ml, so proteins that rise
with neurodegeneration are associated with the label by construction.
`scripts/09_circularity.R` quantifies this using t-tau and Aβ42, the only
biomarkers available at all three sites with AD cases (n = 181; p-tau is
missing in Berlin, Aβ42/40 in Sweden, and Kiel has neither). The
barrier-robust core without MAPT (50 protein groups) is summarised as a signed
score (mean z-score, proteins lower in AD inverted).

The score rises with log2 t-tau within AD (0.36 SD per doubling, p = 6 × 10⁻⁷)
and within controls (0.32, p = 1 × 10⁻⁷), and about half of the core proteins
do so individually in each group (FDR < 5%). Adjusting for log2 t-tau (with
site, age and sex) reduces the AD effect on the score from 0.67 to 0.12 SD
(p = 0.17), and no core protein retains an AD effect at FDR < 5%. Because the
label combines high tau with low amyloid, the tau-adjusted effect mainly
reflects amyloid status; there is no evidence that the core tracks amyloid
beyond tau, although the confidence interval (−0.05 to 0.28) does not exclude a
small effect. In Berlin, the 10 controls with t-tau in the AD range have an
intermediate score (60% of the AD effect, bootstrap 95% CI 27–94%) but less
tau than AD (median 466 vs 804 pg/ml); the tau-adjusted result is the same
within Berlin. Within AD in Berlin, the score correlates with MMSE
(ρ = −0.35, p = 0.046; n = 33, exploratory). A logistic model of t-tau, Aβ42
and site predicts the label with a cross-validated AUC of 0.98, above the
proteome classifier (0.88).

The core is therefore interpreted as a neuronal injury signature that tracks
CSF t-tau, not as an AD-specific signal. Three alternatives were considered and
not used: adjustment for p-tau (unavailable in Berlin), comparison of the
proteome with biomarkers for predicting the label (uninformative, since the
biomarkers define it), and validation against the clinical diagnosis (missing
in Sweden; clinical AD is a subset of biochemical AD and the two disagree in
only 10 individuals).

## Reproducibility

`run_pipeline.sh` runs the full analysis from data download to rendered
reports (about 11 minutes) and logs the git commit, R and Quarto versions, and
the duration of each step. A run from a clean state (processed data, results
and reports deleted) reproduced every tier assignment and all limma-based and
descriptive values exactly. proDA estimates differed between runs by less than
3 × 10⁻⁵ in log2 fold change and 3 × 10⁻⁴ in adjusted p-value, owing to its
iterative numerical optimization, without changing any result tier.

## Software Environment

R 4.6.1 is installed in a dedicated conda environment (`environment.yml`,
exported with `--from-history`, containing R only); R packages are pinned
with renv (`renv.lock`). At the time of setup conda-forge had no `r-renv`
build for R 4.6, so renv is installed from CRAN through the project's
bootstrap script (`renv/activate.R`).

The `Dockerfile` provides the same environment without conda:
`rocker/r-ver:4.6.1`, the packages in `renv.lock` (CRAN binaries from Posit
Package Manager, Bioconductor packages built from source) and Quarto 1.10.18.
The image holds the environment only; the repository is mounted at run time.
A full run in the container from a clean state reproduced every tier
assignment of the conda run, with log2 fold changes differing by less than
3 × 10⁻⁵.
