# 05_missing_value_methods.R
# Compare three treatments of missing values on the same linear model:
#   observed  - limma on observed values only
#   imputed   - limma after left-censored imputation (down-shifted normal)
#   proDA     - probabilistic dropout model, no imputation
# This is a method comparison that justifies the choice of proDA; the final
# differential abundance model and sensitivity analyses are fitted in 06.
#
# Input:  data/processed/04_filter.rds
# Output: data/processed/05_missing_value_methods.rds
#           tibble, one row per protein group: AD vs control log2FC, p and
#           adjusted p for each method, plus missingness per group

library(dplyr)
library(limma)
library(proDA)

in_rds   <- "data/processed/04_filter.rds"
out_file <- "data/processed/05_missing_value_methods.rds"

SEED          <- 20261012
IMPUTE_SHIFT  <- 1.8   # down-shift, in SDs of each run's observed log2 values
IMPUTE_WIDTH  <- 0.3   # width of the imputation distribution, in SDs

dat  <- readRDS(in_rds)
keep <- !dat$sample$excluded
Y    <- log2(dat$intensity[, keep])
smp  <- as.data.frame(dat$sample[keep, ])

# Same model for all three methods; only the handling of missing values differs
model  <- ~ group + batch + age + sex
design <- model.matrix(model, data = smp)

limma_ad <- function(M) {
  fit <- eBayes(lmFit(M, design))
  tt  <- topTable(fit, coef = "groupAD", number = Inf, sort.by = "none")
  tibble(protein = rownames(tt), logFC = tt$logFC, p = tt$P.Value, p_adj = tt$adj.P.Val)
}

impute_downshift <- function(M, shift, width) {
  for (j in seq_len(ncol(M))) {
    obs <- M[!is.na(M[, j]), j]
    na  <- is.na(M[, j])
    M[na, j] <- rnorm(sum(na), mean(obs) - shift * sd(obs), width * sd(obs))
  }
  M
}

# ---- Fit ----------------------------------------------------------------------
res_observed <- limma_ad(Y)

set.seed(SEED)
res_imputed <- limma_ad(impute_downshift(Y, IMPUTE_SHIFT, IMPUTE_WIDTH))

fit_proda <- proDA(Y, design = model, col_data = smp)
tp        <- test_diff(fit_proda, contrast = "groupAD")
res_proda <- tibble(protein = tp$name, logFC = tp$diff, p = tp$pval, p_adj = tp$adj_pval)

# ---- Combine ------------------------------------------------------------------
grp <- smp$group
comparison <- tibble(
  protein          = rownames(Y),
  gene             = dat$protein$PG.Genes,
  missing_control  = rowMeans(is.na(Y[, grp == "control"])),
  missing_AD       = rowMeans(is.na(Y[, grp == "AD"]))
) |>
  left_join(rename_with(res_observed, ~ paste0(.x, "_observed"), -protein), by = "protein") |>
  left_join(rename_with(res_imputed,  ~ paste0(.x, "_imputed"),  -protein), by = "protein") |>
  left_join(rename_with(res_proda,    ~ paste0(.x, "_proda"),    -protein), by = "protein")

# ---- Checks -------------------------------------------------------------------
stopifnot(
  nrow(comparison) == nrow(Y),
  identical(comparison$protein, dat$protein$PG.ProteinGroups),
  !anyNA(comparison$logFC_proda)
)

summarise(comparison,
          across(starts_with("p_adj"), ~ sum(.x < 0.05, na.rm = TRUE), .names = "n_sig_{.col}")) |>
  print()

# ---- Save ---------------------------------------------------------------------
saveRDS(comparison, out_file)
