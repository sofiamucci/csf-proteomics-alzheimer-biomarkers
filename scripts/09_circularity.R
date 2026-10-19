# 09_circularity.R
# How much of the AD signal is explained by the CSF biomarkers that define the
# diagnosis? AD was classified biochemically (high tau and low amyloid; the
# amyloid measure and thresholds differ by site), so proteins that track
# neurodegeneration are associated with the label by construction.
#
#   A. Within-group association of core proteins with t-tau (AD and controls
#      separately), adjusted for site, age and sex
#   B. Berlin: core score in controls with normal tau, controls with AD-range
#      tau and AD. Do tau-high controls look like AD?
#   B2. AD effect on the core score and on each core protein after adjusting
#      for t-tau (all sites with t-tau, and Berlin alone)
#   C. Berlin AD: core score vs MMSE (exploratory)
#   D. Cross-validated AUC of t-tau + Abeta42 for the label (degree of circularity)
#
# t-tau and Abeta42 are the only biomarkers measured at every site with AD
# cases (p-tau is missing in Berlin, Abeta42/40 in Sweden); Kiel has none.
# The core is the barrier-robust tier without MAPT (circular by definition).
#
# Input:  data/processed/04_filter.rds, results/differential_abundance.csv
# Output: data/processed/09_circularity.rds

library(dplyr)
library(readr)
library(limma)
library(pROC)

in_rds   <- "data/processed/04_filter.rds"
in_res   <- "results/differential_abundance.csv"
out_file <- "data/processed/09_circularity.rds"

SEED          <- 20261019
TAU_HIGH      <- 400   # pg/ml; below the lowest t-tau of any AD case at every site
CIRCULAR      <- "MAPT"
COVARIATES    <- ~ site + age + sex
OUTER_FOLDS   <- 5
OUTER_REPEATS <- 10
N_BOOT        <- 2000

dat  <- readRDS(in_rds)
keep <- !dat$sample$excluded
Y    <- log2(dat$intensity[, keep])
smp  <- as.data.frame(dat$sample[keep, ]) |>
  mutate(log2_tau = log2(t_tau), log2_ab42 = log2(abeta42))

# ---- Core and core score ------------------------------------------------------
res  <- read_csv(in_res, show_col_types = FALSE)
core <- res |>
  filter(tier == "barrier-robust",
         !grepl(paste0("(^|;)(", paste(CIRCULAR, collapse = "|"), ")(;|$)"), gene))

# Signed so that a higher score is more AD-like for every protein
Z          <- t(scale(t(Y[core$protein_group, ])))
core_sign  <- sign(core$log2FC)
smp$core_score <- colMeans(Z * core_sign, na.rm = TRUE)

# ---- A. Within-group association with t-tau -----------------------------------
tau_assoc <- function(g) {
  s  <- filter(smp, group == g, !is.na(t_tau)) |> mutate(site = droplevels(site))
  d  <- model.matrix(update(COVARIATES, ~ log2_tau + .), data = s)
  ft <- eBayes(lmFit(Y[core$protein_group, s$run], d))
  tt <- topTable(ft, coef = "log2_tau", number = Inf, sort.by = "none")
  score_fit <- summary(lm(update(COVARIATES, core_score ~ log2_tau + .), data = s))
  list(
    proteins = tibble(group = g, protein_group = rownames(tt),
                      gene = core$gene[match(rownames(tt), core$protein_group)],
                      slope = tt$logFC, p = tt$P.Value, p_adj = tt$adj.P.Val),
    score = tibble(group = g, n = nrow(s),
                   slope = coef(score_fit)["log2_tau", "Estimate"],
                   se    = coef(score_fit)["log2_tau", "Std. Error"],
                   p     = coef(score_fit)["log2_tau", "Pr(>|t|)"])
  )
}
assoc       <- lapply(c("AD", "control"), tau_assoc)
tau_protein <- bind_rows(lapply(assoc, `[[`, "proteins"))
tau_score   <- bind_rows(lapply(assoc, `[[`, "score"))

# ---- B. Berlin: tau-high controls ---------------------------------------------
berlin <- smp |>
  filter(site == "Berlin") |>
  mutate(tau_group = case_when(
    group == "AD"      ~ "AD",
    t_tau > TAU_HIGH   ~ "control, high tau",
    .default           = "control, normal tau"
  ),
  tau_group = factor(tau_group, levels = c("control, normal tau", "control, high tau", "AD")))

fit_b    <- lm(core_score ~ tau_group + age + sex, data = berlin)
coef_b   <- coef(summary(fit_b))
contrasts_b <- tibble(
  contrast = c("control, high tau", "AD"),
  estimate = coef_b[paste0("tau_group", contrast), "Estimate"],
  se       = coef_b[paste0("tau_group", contrast), "Std. Error"],
  p        = coef_b[paste0("tau_group", contrast), "Pr(>|t|)"]
)

# Position of tau-high controls between normal controls (0) and AD (1)
fraction_of_ad <- function(d) {
  b <- coef(lm(core_score ~ tau_group + age + sex, data = d))
  unname(b["tau_groupcontrol, high tau"] / b["tau_groupAD"])
}
set.seed(SEED)
boot_frac <- replicate(N_BOOT, {
  idx <- unlist(lapply(split(seq_len(nrow(berlin)), berlin$tau_group),
                       function(i) sample(i, replace = TRUE)))   # stratified
  fraction_of_ad(berlin[idx, ])
})
position_b <- tibble(fraction_of_ad = fraction_of_ad(berlin),
                     ci_low  = quantile(boot_frac, 0.025, names = FALSE),
                     ci_high = quantile(boot_frac, 0.975, names = FALSE))

# Per protein, descriptive
d_b  <- model.matrix(~ tau_group + age + sex, data = berlin)
ft_b <- eBayes(lmFit(Y[core$protein_group, berlin$run], d_b))
protein_b <- tibble(
  protein_group = core$protein_group,
  gene          = core$gene,
  high_tau_vs_normal = ft_b$coefficients[, "tau_groupcontrol, high tau"],
  ad_vs_normal       = ft_b$coefficients[, "tau_groupAD"]
)

# ---- B2. AD effect beyond tau -------------------------------------------------
# Tau-high controls have less tau than AD (Berlin medians 466 vs 804 pg/ml), so
# their intermediate position may reflect tau dose alone. Diagnosis is tested
# with log2 t-tau as a covariate, in all sites with t-tau and in Berlin alone.
tau_adjusted <- function(d) {
  f  <- update(COVARIATES, core_score ~ group + log2_tau + .)
  if (n_distinct(d$site) == 1) f <- update(f, ~ . - site)
  cf <- coef(summary(lm(f, data = d)))
  unadj <- coef(summary(lm(update(f, ~ . - log2_tau), data = d)))
  tibble(n = nrow(d),
         ad_unadjusted = unadj["groupAD", "Estimate"],
         ad_tau_adjusted = cf["groupAD", "Estimate"],
         se = cf["groupAD", "Std. Error"],
         p  = cf["groupAD", "Pr(>|t|)"],
         tau_slope = cf["log2_tau", "Estimate"])
}
with_tau <- filter(smp, !is.na(t_tau)) |> mutate(site = droplevels(site))
beyond_tau <- bind_rows(
  mutate(tau_adjusted(with_tau), sites = "Sweden, Magdeburg, Berlin"),
  mutate(tau_adjusted(filter(with_tau, site == "Berlin")), sites = "Berlin")
)

# Per protein: which core proteins keep an AD effect after adjusting for tau
d_e  <- model.matrix(update(COVARIATES, ~ group + log2_tau + .), data = with_tau)
ft_e <- eBayes(lmFit(Y[core$protein_group, with_tau$run], d_e))
tt_e <- topTable(ft_e, coef = "groupAD", number = Inf, sort.by = "none")
protein_beyond_tau <- tibble(protein_group = rownames(tt_e),
                             gene = core$gene[match(rownames(tt_e), core$protein_group)],
                             log2FC_tau_adjusted = tt_e$logFC,
                             p = tt_e$P.Value, p_adj = tt_e$adj.P.Val)

# ---- C. Berlin AD: cognition --------------------------------------------------
berlin_ad <- filter(berlin, group == "AD", !is.na(mmse))
mmse_test <- cor.test(berlin_ad$core_score, berlin_ad$mmse, method = "spearman", exact = FALSE)
mmse_b    <- tibble(n = nrow(berlin_ad), rho = unname(mmse_test$estimate), p = mmse_test$p.value)

# ---- D. Circularity: CSF biomarkers predicting the label ----------------------
bio <- filter(smp, !is.na(t_tau), !is.na(abeta42)) |> mutate(site = droplevels(site))
y   <- as.integer(bio$group == "AD")

stratified_folds <- function(y, k) {
  f <- integer(length(y))
  for (cl in unique(y)) {
    idx    <- which(y == cl)
    f[idx] <- sample(rep_len(seq_len(k), length(idx)))
  }
  f
}
auc_of <- function(y, prob) as.numeric(auc(roc(y, prob, direction = "<", quiet = TRUE)))

set.seed(SEED)
auc_bio <- sapply(seq_len(OUTER_REPEATS), function(r) {
  folds <- stratified_folds(y, OUTER_FOLDS)
  prob  <- numeric(length(y))
  for (k in seq_len(OUTER_FOLDS)) {
    test  <- folds == k
    fit   <- glm(group ~ log2_tau + log2_ab42 + site, family = binomial, data = bio[!test, ])
    prob[test] <- predict(fit, newdata = bio[test, ], type = "response")
  }
  auc_of(y, prob)
})
biomarker_auc <- tibble(n = nrow(bio), n_AD = sum(y),
                        auc_mean = mean(auc_bio), auc_min = min(auc_bio), auc_max = max(auc_bio))

# ---- Checks -------------------------------------------------------------------
stopifnot(
  nrow(core) == 50L,
  !any(grepl("(^|;)MAPT(;|$)", core$gene)),
  all(core$protein_group %in% rownames(Y)),
  sum(!is.na(smp$t_tau)) == 181L,
  all(smp$t_tau[smp$group == "AD"] > TAU_HIGH),               # threshold is below every AD case
  identical(as.integer(table(berlin$tau_group)), c(40L, 10L, 33L)),
  nrow(bio) == 181L
)

cat("Core proteins:", nrow(core), "\n\nA. Core score vs log2 t-tau, within group\n")
print(tau_score)
print(count(tau_protein, group, significant = p_adj < 0.05))
cat("\nB. Berlin, core score relative to controls with normal tau\n")
print(contrasts_b)
print(position_b)
cat("\nB2. AD effect on the core score, with and without adjusting for t-tau\n")
print(beyond_tau)
print(filter(protein_beyond_tau, p_adj < 0.05) |> arrange(p_adj))
cat("\nC. Berlin AD, core score vs MMSE\n")
print(mmse_b)
cat("\nD. Cross-validated AUC of t-tau + Abeta42 + site for the AD label\n")
print(biomarker_auc)

# ---- Save ---------------------------------------------------------------------
saveRDS(list(core = core,
             core_score = select(smp, run, site, group, t_tau, abeta42, mmse, core_score),
             tau_protein = tau_protein, tau_score = tau_score,
             berlin = select(berlin, run, tau_group, t_tau, core_score),
             berlin_contrasts = contrasts_b, berlin_position = position_b,
             berlin_protein = protein_b,
             beyond_tau = beyond_tau, protein_beyond_tau = protein_beyond_tau,
             mmse = mmse_b,
             biomarker_auc = biomarker_auc, biomarker_auc_runs = auc_bio,
             params = list(tau_high = TAU_HIGH, circular = CIRCULAR,
                           covariates = COVARIATES, seed = SEED, n_boot = N_BOOT)),
        out_file)
