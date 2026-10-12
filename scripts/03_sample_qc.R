# 03_sample_qc.R
# Per-run QC metrics, contamination scores and flags used downstream.
#
# Input:  data/processed/02_annotate.rds
# Output: data/processed/03_sample_qc.rds
#           list(intensity, protein, sample)
#           sample  gains per-run QC metrics and flags (NA for runs excluded by
#                   the authors, which are not scored)
#           protein gains contamination-panel membership

library(dplyr)

in_rds   <- "data/processed/02_annotate.rds"
out_file <- "data/processed/03_sample_qc.rds"

# Erythrocyte panel: Geyer et al. (2019) EMBO Mol Med
ERY_PANEL    <- c("HBA1", "HBB", "HBD", "CA1", "CA2", "CAT", "PRDX2", "BLVRB", "SLC4A1")
PLASMA_PANEL <- c("ALB", "TF", "HPX", "SERPINA1", "AHSG", "APOA1", "GC", "PLG", "F2", "HP")
MAD_K        <- 3   # flagging threshold, in median absolute deviations

dat <- readRDS(in_rds)

# Scores are computed on included runs only, so excluded runs do not shift
# the z-scores of the others
keep     <- !dat$sample$excluded
X        <- log2(dat$intensity[, keep])
complete <- rowSums(is.na(X)) == 0

# ---- Contamination panels -----------------------------------------------------
gene_list <- strsplit(dat$protein$PG.Genes, ";")
in_panel  <- function(panel) vapply(gene_list, function(g) any(g %in% panel), logical(1))

dat$protein$ery_panel    <- in_panel(ERY_PANEL)
dat$protein$plasma_panel <- in_panel(PLASMA_PANEL)

# Mean z-score of panel proteins per run; 0 = typical level in this cohort
panel_score <- function(rows) colMeans(t(scale(t(X[rows, ]))), na.rm = TRUE)

# ---- Per-run metrics and flags ------------------------------------------------
qc <- tibble(
  run             = colnames(X),
  n_quantified    = colSums(!is.na(X)),
  frac_missing    = colMeans(is.na(X)),
  median_observed = apply(X, 2, median, na.rm = TRUE),
  median_complete = apply(X[complete, ], 2, median),
  ery_score       = panel_score(dat$protein$ery_panel),
  plasma_score    = panel_score(dat$protein$plasma_panel)
) |>
  left_join(select(dat$sample, run, batch), by = "run") |>
  # Coverage is judged within batch, since batches differ in depth
  group_by(batch) |>
  mutate(low_coverage = n_quantified < median(n_quantified) - MAD_K * mad(n_quantified)) |>
  ungroup() |>
  mutate(ery_flag = ery_score > median(ery_score) + MAD_K * mad(ery_score)) |>
  select(-batch)

sample <- left_join(dat$sample, qc, by = "run")

# ---- Checks -------------------------------------------------------------------
stopifnot(
  identical(sample$run, colnames(dat$intensity)),
  sum(complete) == 385L,
  sum(dat$protein$ery_panel) == length(ERY_PANEL),
  !anyNA(sample$ery_score[keep]),
  all(is.na(sample$ery_score[!keep]))
)

print(table(low_coverage = sample$low_coverage[keep],
            ery_flag     = sample$ery_flag[keep]))

# ---- Save ---------------------------------------------------------------------
saveRDS(list(intensity = dat$intensity, protein = dat$protein, sample = sample),
        out_file)
