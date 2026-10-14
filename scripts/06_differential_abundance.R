# 06_differential_abundance.R
# Differential abundance, AD vs control: main proDA model, contamination flags
# and sensitivity analyses.
#
# Input:  data/processed/04_filter.rds
# Output: data/processed/06_differential_abundance.rds
#           list(results             one row per protein group, main model + flags
#                sensitivity         list of per-analysis result tibbles
#                sensitivity_summary agreement of each analysis with the main model
#                model               main model formula)

library(dplyr)
library(limma)
library(proDA)

in_rds   <- "data/processed/04_filter.rds"
out_file <- "data/processed/06_differential_abundance.rds"

# Site (4 levels) rather than MS batch: each site was measured in one batch, and
# site separates Kiel (controls only) from Magdeburg, so Kiel does not inform
# the AD vs control contrast
MODEL    <- ~ group + site + age + sex
FDR      <- 0.05
CONTAM_R <- 0.5   # Spearman r with a contamination score that flags a protein group

dat  <- readRDS(in_rds)
keep <- !dat$sample$excluded
Y    <- log2(dat$intensity[, keep])
smp  <- as.data.frame(dat$sample[keep, ])

fit_proda <- function(M, s, model = MODEL) {
  tp <- test_diff(proDA(M, design = model, col_data = s), contrast = "groupAD")
  tibble(protein = tp$name, logFC = tp$diff, se = tp$se, p = tp$pval,
         p_adj = tp$adj_pval, avg_abundance = tp$avg_abundance, n_obs = tp$n_obs)
}

# ---- Main model ---------------------------------------------------------------
main <- fit_proda(Y, smp)

# ---- Contamination flags ------------------------------------------------------
# Protein groups that track a contamination score across samples
cor_with <- function(score) {
  apply(Y, 1, function(y) cor(y, score, method = "spearman", use = "pairwise.complete.obs"))
}
r_ery    <- cor_with(smp$ery_score)
r_plasma <- cor_with(smp$plasma_score)

results <- tibble(
  protein         = rownames(Y),
  gene            = dat$protein$PG.Genes,
  missing_control = rowMeans(is.na(Y[, smp$group == "control"])),
  missing_AD      = rowMeans(is.na(Y[, smp$group == "AD"]))
) |>
  left_join(main, by = "protein") |>
  mutate(
    significant       = p_adj < FDR,
    direction         = if_else(logFC > 0, "higher in AD", "lower in AD"),
    r_erythrocyte     = r_ery,
    r_plasma          = r_plasma,
    blood_associated  = r_ery > CONTAM_R | dat$protein$ery_panel,
    plasma_associated = r_plasma > CONTAM_R | dat$protein$plasma_panel
  )

# ---- Sensitivity analyses -----------------------------------------------------
# 1. limma on observed values with sample quality weights (proDA has no weights)
design   <- model.matrix(MODEL, data = smp)
aw       <- arrayWeights(Y, design)
tt       <- topTable(eBayes(lmFit(Y, design, weights = aw)),
                     coef = "groupAD", number = Inf, sort.by = "none")
s_limma  <- tibble(protein = rownames(tt), logFC = tt$logFC, p = tt$P.Value, p_adj = tt$adj.P.Val)

# 2. Without the runs flagged for blood contamination
clean    <- !smp$ery_flag
s_noblood <- fit_proda(Y[, clean], smp[clean, ])

# 3. Adjusting for the plasma-protein score (blood-CSF barrier)
s_plasma <- fit_proda(Y, smp, model = update(MODEL, ~ . + plasma_score))

# 4. Renormalized: each run centred on the median of complete protein groups
complete <- rowSums(is.na(Y)) == 0
offset   <- apply(Y[complete, ], 2, median)
s_renorm <- fit_proda(sweep(Y, 2, offset - median(offset)), smp)

sensitivity <- list(
  limma_arrayweights = s_limma,
  without_blood_flag = s_noblood,
  plasma_adjusted    = s_plasma,
  renormalized       = s_renorm
)

agreement <- function(alt, name) {
  m <- inner_join(select(main, protein, logFC_main = logFC, sig_main = p_adj),
                  select(alt,  protein, logFC_alt  = logFC, sig_alt  = p_adj),
                  by = "protein")
  tibble(
    analysis      = name,
    r_logFC       = cor(m$logFC_main, m$logFC_alt, use = "complete.obs"),
    n_significant = sum(m$sig_alt < FDR, na.rm = TRUE),
    shared        = sum(m$sig_main < FDR & m$sig_alt < FDR, na.rm = TRUE),
    lost          = sum(m$sig_main < FDR & !(m$sig_alt < FDR), na.rm = TRUE),
    gained        = sum(!(m$sig_main < FDR) & m$sig_alt < FDR, na.rm = TRUE)
  )
}
sensitivity_summary <- bind_rows(Map(agreement, sensitivity, names(sensitivity)))

# ---- Checks -------------------------------------------------------------------
stopifnot(
  identical(results$protein, rownames(Y)),
  identical(results$protein, dat$protein$PG.ProteinGroups),
  !anyNA(results$logFC),
  all(levels(smp$site) %in% unique(as.character(smp$site)))
)

count(results, significant, direction) |> print()
print(sensitivity_summary)
results |>
  filter(gene %in% c("MAPT", "SOD1", "PARK7", "NEFL", "UCHL1")) |>
  select(gene, logFC, p_adj, blood_associated, plasma_associated) |>
  print()

# ---- Save ---------------------------------------------------------------------
saveRDS(list(results = results, sensitivity = sensitivity,
             sensitivity_summary = sensitivity_summary,
             model = deparse(MODEL)),
        out_file)
