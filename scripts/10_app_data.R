# 10_app_data.R
# Data for the interactive viewer (app/): main-model results, per-cohort
# results, log2 intensities and minimal sample annotation.
#
# Per-cohort results fit the main model within each site with AD cases
# (~ group + age + sex). They show whether an effect is consistent across
# cohorts; with 38-83 samples per site they have less power than the main
# model and are exploratory. Kiel (controls only) contributes no contrast.
#
# Input:  data/processed/04_filter.rds, results/differential_abundance.csv
# Output: app/app_data.rds

library(dplyr)
library(readr)
library(proDA)

in_rds   <- "data/processed/04_filter.rds"
in_res   <- "results/differential_abundance.csv"
out_file <- "app/app_data.rds"

COHORT_MODEL <- ~ group + age + sex
COHORTS      <- c("Sweden", "Magdeburg", "Berlin")

dat  <- readRDS(in_rds)
keep <- !dat$sample$excluded
Y    <- log2(dat$intensity[, keep])
smp  <- as.data.frame(dat$sample[keep, ])
res  <- read_csv(in_res, show_col_types = FALSE)

# ---- Per-cohort models --------------------------------------------------------
fit_cohort <- function(site) {
  s  <- smp[smp$site == site, ]
  tp <- test_diff(proDA(Y[, s$run], design = COHORT_MODEL, col_data = s),
                  contrast = "groupAD")
  tibble(protein_group = tp$name, cohort = site,
         n_control = sum(s$group == "control"), n_AD = sum(s$group == "AD"),
         log2FC = tp$diff, se = tp$se, p_value = tp$pval, p_adj = tp$adj_pval)
}
cohort <- bind_rows(lapply(COHORTS, fit_cohort))

main <- res |>
  transmute(protein_group, cohort = "All (main model)",
            n_control = sum(smp$group == "control"), n_AD = sum(smp$group == "AD"),
            log2FC, se, p_value, p_adj)

# ---- Assemble -----------------------------------------------------------------
proteins <- res |>
  transmute(protein_group, gene = coalesce(gene, protein_group), tier, circular,
            core = tier == "barrier-robust" & !circular,
            r_plasma, plasma_associated, blood_associated)

samples <- smp |>
  transmute(run, group = factor(group, levels = c("control", "AD")),
            site = droplevels(site), t_tau)

intensity <- round(Y[proteins$protein_group, samples$run], 3)

app_data <- list(
  proteins  = proteins,
  results   = bind_rows(main, cohort),
  samples   = samples,
  intensity = intensity,
  meta      = list(source = "PRIDE PXD016278 (Bader et al., 2020)",
                   main_model = "proDA, ~ group + site + age + sex",
                   cohort_model = paste(deparse(COHORT_MODEL), collapse = ""),
                   created = format(Sys.Date()))
)

# ---- Checks -------------------------------------------------------------------
stopifnot(
  nrow(proteins) == 1317L,
  sum(proteins$core) == 50L,
  identical(rownames(intensity), proteins$protein_group),
  ncol(intensity) == 197L,
  all(table(app_data$results$cohort) == 1317L),
  setequal(unique(app_data$results$cohort), c("All (main model)", COHORTS))
)

print(app_data$results |>
        group_by(cohort) |>
        summarise(n = first(n_control) + first(n_AD),
                  significant = sum(p_adj < 0.05, na.rm = TRUE),
                  missing_estimates = sum(is.na(log2FC))))

# Agreement of each cohort with the main model on the core
print(app_data$results |>
        filter(protein_group %in% proteins$protein_group[proteins$core]) |>
        select(protein_group, cohort, log2FC) |>
        tidyr::pivot_wider(names_from = cohort, values_from = log2FC) |>
        summarise(across(all_of(COHORTS),
                         list(r_with_main = ~ cor(.x, `All (main model)`, use = "complete.obs"),
                              same_sign   = ~ mean(sign(.x) == sign(`All (main model)`), na.rm = TRUE)))))

# ---- Save ---------------------------------------------------------------------
dir.create(dirname(out_file), showWarnings = FALSE)
saveRDS(app_data, out_file, compress = "xz")
cat("Saved", out_file, "(", round(file.size(out_file) / 1e6, 1), "MB )\n")
