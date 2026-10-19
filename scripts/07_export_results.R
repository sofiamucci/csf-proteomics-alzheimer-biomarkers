# 07_export_results.R
# Export differential abundance results as a CSV deliverable with a data
# dictionary, and assign each protein group to a result tier.
#
# Input:  data/processed/06_differential_abundance.rds
# Output: results/differential_abundance.csv
#         results/differential_abundance_dictionary.csv

library(dplyr)
library(readr)

in_rds  <- "data/processed/06_differential_abundance.rds"
out_dir <- "results"
FDR     <- 0.05
CIRCULAR_PATTERN <- "(^|;)MAPT(;|$)"   # diagnosis is defined from CSF tau

da <- readRDS(in_rds)
s  <- da$sensitivity

out <- da$results |>
  left_join(select(s$plasma_adjusted, protein,
                   log2FC_plasma_adjusted = logFC, p_adj_plasma_adjusted = p_adj),
            by = "protein") |>
  left_join(select(s$limma_arrayweights, protein,
                   log2FC_limma_arrayweights = logFC, p_adj_limma_arrayweights = p_adj),
            by = "protein") |>
  left_join(select(s$without_blood_flag, protein, p_adj_without_blood_flag = p_adj),
            by = "protein") |>
  left_join(select(s$renormalized, protein, p_adj_renormalized = p_adj),
            by = "protein") |>
  mutate(tier = case_when(
    significant & p_adj_plasma_adjusted < FDR ~ "barrier-robust",
    significant                               ~ "main model only",
    .default                                  = "not significant"
  ),
  circular = grepl(CIRCULAR_PATTERN, gene)) |>
  arrange(p_adj) |>
  transmute(
    protein_group = protein, gene, tier, circular, direction,
    log2FC = logFC, se, p_value = p, p_adj,
    log2FC_plasma_adjusted, p_adj_plasma_adjusted,
    log2FC_limma_arrayweights, p_adj_limma_arrayweights,
    p_adj_without_blood_flag, p_adj_renormalized,
    missing_control, missing_AD, avg_abundance, n_obs,
    r_erythrocyte, r_plasma, blood_associated, plasma_associated
  )

dictionary <- tribble(
  ~column,                     ~description,
  "protein_group",             "Spectronaut protein group (UniProt accessions)",
  "gene",                      "Gene name(s) of the protein group",
  "tier",                      "barrier-robust: FDR < 5% in the main and plasma-adjusted models; main model only: FDR < 5% in the main model only; not significant",
  "circular",                  "TRUE for MAPT: associated with diagnosis by definition (CSF tau); excluded from the barrier-robust core in downstream analyses",
  "direction",                 "Sign of log2FC (AD relative to control)",
  "log2FC",                    "Main model (proDA, ~ group + site + age + sex): log2 fold change, AD vs control",
  "se",                        "Standard error of log2FC",
  "p_value",                   "Main model p-value",
  "p_adj",                     "Main model Benjamini-Hochberg adjusted p-value",
  "log2FC_plasma_adjusted",    "log2FC with the plasma-protein score as an additional covariate",
  "p_adj_plasma_adjusted",     "Adjusted p-value, plasma-adjusted model",
  "log2FC_limma_arrayweights", "log2FC from limma on observed values with sample quality weights",
  "p_adj_limma_arrayweights",  "Adjusted p-value, limma with sample quality weights",
  "p_adj_without_blood_flag",  "Adjusted p-value, main model without the 4 runs flagged for blood contamination",
  "p_adj_renormalized",        "Adjusted p-value, main model after per-run median centring",
  "missing_control",           "Fraction of control runs with a missing value",
  "missing_AD",                "Fraction of AD runs with a missing value",
  "avg_abundance",             "Average log2 abundance estimated by proDA",
  "n_obs",                     "Number of observed values",
  "r_erythrocyte",             "Spearman correlation with the erythrocyte contamination score",
  "r_plasma",                  "Spearman correlation with the plasma-protein score",
  "blood_associated",          "Erythrocyte-panel protein or r_erythrocyte > 0.5",
  "plasma_associated",         "Plasma-panel protein or r_plasma > 0.5"
)

# ---- Checks -------------------------------------------------------------------
stopifnot(
  identical(sort(names(out)), sort(dictionary$column)),
  sum(out$tier != "not significant") == sum(da$results$significant),
  !anyDuplicated(out$protein_group),
  sum(out$circular, na.rm = TRUE) == 1L
)

print(count(out, tier, direction))

# ---- Save ---------------------------------------------------------------------
dir.create(out_dir, showWarnings = FALSE)
write_csv(out, file.path(out_dir, "differential_abundance.csv"))
write_csv(dictionary, file.path(out_dir, "differential_abundance_dictionary.csv"))
