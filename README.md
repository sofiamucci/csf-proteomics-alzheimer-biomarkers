# What does Alzheimer's leave in the spinal fluid? A metabolic fingerprint

*A reanalysis of the cerebrospinal fluid proteome of 197 people, from the raw
protein report to a validated classifier, with every analytical decision
documented.*

## The Question

Cerebrospinal fluid (CSF) bathes the brain, and what neurons and glia release
ends up there. When Alzheimer's disease (AD) damages the brain, which proteins
change in CSF, and which of those changes are specific to the disease rather
than to how the sample was taken or how much blood plasma leaks into the fluid?

This project reanalyzes a public mass spectrometry dataset (PRIDE PXD016278;
Bader et al., 2020) of CSF from three cohorts: 88 people with biochemically
defined AD and 109 controls.

## The Finding

**AD leaves a metabolic fingerprint in the CSF.** Of 1,317 quantified protein
groups, 259 differ between AD and controls (FDR < 5%), most of them higher in
AD. The changes that survive every robustness check form a core of 51 proteins
dominated by **glycolytic enzymes** (ALDOA, ALDOC, PKM, LDHA, GPI, ENO2, PGK1)
and **14-3-3 proteins** (YWHAZ, YWHAG, YWHAE), together with established AD
markers such as SOD1, UCHL1, SMOC1 and SPP1. These are intracellular brain
proteins, released into the CSF as the disease progresses.

![Volcano plot, AD vs control](figures/volcano.png)

**The specific signal is what the brain releases, not what the blood brings.**
CSF in AD contains relatively less plasma-derived protein. Adjusting for that
balance removes most of the proteins that are *lower* in AD (61 of 67), while
most of those *higher* in AD remain.

![Barrier-robust core across all samples](figures/heatmap_core.png)

## Is It Robust?

The main result was challenged with four sensitivity analyses: an alternative
statistical model with sample quality weights, removal of blood-contaminated
samples, renormalization, and adjustment for the plasma-protein balance. Fold
changes were essentially unchanged in the first three (r ≥ 0.96). The plasma
adjustment reduces the significant set to the 51-protein core, so results are
reported in two tiers. A clean rerun of the full pipeline reproduced every
result tier.

The original study's positive controls are recovered (SOD1, PARK7, and tau as
the strongest effect), and the glycolytic signature it reported is replicated
independently.

## Can It Diagnose?

An elastic-net classifier built only on the proteome (tau excluded, since it
defines the diagnosis) separates AD from controls with an **AUC of 0.88** in
nested cross-validation, against 0.60 for age and sex alone. Trained on two
hospitals and tested on the third, it reaches 0.86 and 0.88 for two cohorts
and 0.74 for the third, whose controls are mostly patients with other
neurological diagnoses. The proteins it relies on most are the same 14-3-3 and
glycolytic proteins found by the differential analysis, although the
classifier never saw those results.

![ROC curves](figures/roc.png)

## How This Was Built

The analysis starts from the authors' Spectronaut protein report (DIA) and runs
in nine scripts plus two reports. Key decisions, each justified in
[METHODS.md](METHODS.md):

| Step | Decision |
|---|---|
| Sample QC | Batch offsets, blood contamination and plasma content measured per sample; outliers explained by covariates before any exclusion |
| Filtering | Lab contaminants removed; proteins kept if quantified in ≥ 50% of at least one group, which retains NEFL and UCHL1 |
| Missing values | Missingness depends on diagnosis, so it is modelled with proDA instead of imputed; imputation was shown to inflate effects |
| Model | `~ diagnosis + site + age + sex`; site separates an all-control cohort from the others |
| Classifier | Every data-dependent step learned inside training folds; validated across hospitals |

**Reports:** [QC report](https://sofiamucci.github.io/csf-proteomics-alzheimer-biomarkers/qc_report.html) ·
[Analysis report](https://sofiamucci.github.io/csf-proteomics-alzheimer-biomarkers/analysis_report.html)

**Results table:** [results/differential_abundance.csv](results/differential_abundance.csv),
with a [column dictionary](results/differential_abundance_dictionary.csv).

## Reproduce

Requires conda (or mamba) and Quarto.

```bash
mamba env create -f environment.yml
mamba activate csf-proteomics
R -e "renv::restore()"
bash run_pipeline.sh
```

`run_pipeline.sh` downloads the data from PRIDE, runs every step and renders
both reports (about 11 minutes), logging the git commit and software versions.

## Repository Structure

```
├── run_pipeline.sh                 # runs everything, in order
├── scripts/
│   ├── 00_download_data.sh         # PRIDE download with checksums
│   ├── 01_import.R                 # Spectronaut report -> matrix
│   ├── 02_annotate.R               # sample annotation, author exclusions
│   ├── 03_sample_qc.R              # per-sample metrics and contamination scores
│   ├── 04_filter.R                 # contaminants and valid-value filter
│   ├── 05_missing_value_methods.R  # observed vs imputed vs proDA
│   ├── 06_differential_abundance.R # main model and sensitivity analyses
│   ├── 07_export_results.R         # results table and dictionary
│   └── 08_classifier.R             # elastic net, nested and cross-site validation
├── reports/                        # QC and analysis reports (Quarto)
├── results/                        # differential abundance results
├── figures/                        # figures used in this README
├── METHODS.md
├── environment.yml                 # R version (conda)
└── renv.lock                       # R package versions
```

## Limitations

- **Diagnosis is biochemical**, defined from CSF tau and amyloid-β. The tau
  result is circular, and proteins correlated with tau may be partly so.
- **Confounder or mediator?** Without a CSF/serum albumin quotient, the
  plasma-protein balance cannot be classified as barrier permeability or as
  part of the disease process; hence the two-tier reporting.
- **Low-abundance proteins** such as NEFL are quantified near the detection
  limit by mass spectrometry; their absence from the results is not evidence
  of no change.
- **Site and measurement batch coincide**, so technical and collection effects
  cannot be separated.
- **No external cohort.** Cross-site validation uses cohorts from the same
  study and laboratory.

## AI-Assisted Development

This analysis was developed with Claude (Anthropic) as an assistant. I used it
to structure the pipeline, draft and debug R and Bash code, and discuss
methodological options, such as how to treat missing values, which covariates
to model and how to validate the classifier. Each analytical decision was
taken after examining the data and is recorded with its justification in
[METHODS.md](METHODS.md). Code was run, checked against expected values
(assertions in every script, equivalence checks after refactoring, a clean
end-to-end rerun) and reviewed before committing. Where an expected result did
not hold, the data prevailed and the discrepancy is documented.

## Data Source

Bader JM, Geyer PE, et al. Proteome profiling in cerebrospinal fluid reveals
novel biomarkers of Alzheimer's disease. *Mol Syst Biol* 16:e9356 (2020).
[doi:10.15252/msb.20199356](https://doi.org/10.15252/msb.20199356). Data:
[PRIDE PXD016278](https://www.ebi.ac.uk/pride/archive/projects/PXD016278).

## License

MIT
