# 08_classifier.R
# AD vs control classifier from the CSF proteome: elastic-net logistic
# regression, evaluated by repeated nested cross-validation and by
# leave-one-site-out validation, against an age + sex baseline.
#
# Every data-dependent step (missing-value filter, protein selection,
# imputation, scaling, penalty) is learned on the training data of each split
# only. Out-of-fold predictions are averaged over repetitions to give one
# prediction per sample, used for a DeLong confidence interval, sensitivity and
# specificity at a fixed threshold, and calibration.
#
# Input:  data/processed/04_filter.rds
# Output: data/processed/08_classifier.rds

library(dplyr)
library(glmnet)
library(pROC)

in_rds   <- "data/processed/04_filter.rds"
out_file <- "data/processed/08_classifier.rds"

SEED          <- 20261014
ALPHA         <- 0.5     # elastic net mixing (0 = ridge, 1 = lasso)
OUTER_FOLDS   <- 5
OUTER_REPEATS <- 10
INNER_FOLDS   <- 5
MAX_MISSING   <- 0.2     # in the training split; protein groups above are not used
EXCLUDE_GENES <- "MAPT"  # diagnosis is defined from CSF tau
EXCLUDE_SITES <- "Kiel"  # controls only and age-mismatched
THRESHOLD     <- 0.5     # for sensitivity and specificity
CALIB_BINS    <- 5

dat <- readRDS(in_rds)
use <- !dat$sample$excluded & !(dat$sample$site %in% EXCLUDE_SITES)
smp <- dat$sample[use, ] |> mutate(site = droplevels(site))
y   <- as.integer(smp$group == "AD")

X     <- t(log2(dat$intensity[, use]))                 # runs x protein groups
genes <- dat$protein$PG.Genes
excl  <- grepl(paste0("(^|;)(", paste(EXCLUDE_GENES, collapse = "|"), ")(;|$)"), genes)
# Candidate set: no protein group quantified in fewer than 20% of runs overall
# can pass the per-split filter, so these are dropped here only to save time
feat  <- colMeans(is.na(X)) < 0.8 & !excl
X     <- X[, feat]
genes <- genes[feat]

# ---- Preprocessing learned on training data only ------------------------------
prep_fit <- function(Xtr) {
  keep <- colMeans(is.na(Xtr)) <= MAX_MISSING
  Xtr  <- Xtr[, keep, drop = FALSE]
  fill <- apply(Xtr, 2, function(v) min(v, na.rm = TRUE) - 1)   # half the minimum
  Xf   <- prep_apply(Xtr, list(keep = rep(TRUE, ncol(Xtr)), fill = fill, mu = 0, sd = 1))
  list(keep = keep, fill = fill, mu = colMeans(Xf), sd = apply(Xf, 2, sd))
}
prep_apply <- function(M, p) {
  M <- M[, p$keep, drop = FALSE]
  for (j in seq_len(ncol(M))) M[is.na(M[, j]), j] <- p$fill[j]
  sweep(sweep(M, 2, p$mu), 2, p$sd, "/")
}

fit_predict <- function(train, test) {
  p  <- prep_fit(X[train, , drop = FALSE])
  cv <- cv.glmnet(prep_apply(X[train, , drop = FALSE], p), y[train],
                  family = "binomial", alpha = ALPHA,
                  nfolds = INNER_FOLDS, type.measure = "auc")
  coef_full <- numeric(ncol(X))                      # 0 for filtered-out proteins
  coef_full[p$keep] <- as.numeric(coef(cv, s = "lambda.1se"))[-1]
  list(prob = predict(cv, prep_apply(X[test, , drop = FALSE], p),
                      s = "lambda.1se", type = "response")[, 1],
       coef = coef_full, n_used = sum(p$keep))
}

baseline_predict <- function(train, test) {
  fit <- glm(group ~ age + sex, family = binomial, data = smp[train, ])
  predict(fit, newdata = smp[test, ], type = "response")
}

# Direction fixed: controls are expected to have lower predicted probability.
# Letting pROC choose the direction would inflate AUCs near 0.5.
auc_of <- function(y, prob) as.numeric(auc(roc(y, prob, direction = "<", quiet = TRUE)))

stratified_folds <- function(y, k) {
  f <- integer(length(y))
  for (cl in unique(y)) {
    idx    <- which(y == cl)
    f[idx] <- sample(rep_len(seq_len(k), length(idx)))
  }
  f
}

# ---- Repeated nested cross-validation -----------------------------------------
set.seed(SEED)
cv_runs   <- vector("list", OUTER_REPEATS)
selected  <- matrix(0, nrow = OUTER_REPEATS * OUTER_FOLDS, ncol = ncol(X))
n_used    <- integer(OUTER_REPEATS * OUTER_FOLDS)
i <- 0
for (r in seq_len(OUTER_REPEATS)) {
  folds <- stratified_folds(y, OUTER_FOLDS)
  prob  <- prob_base <- numeric(length(y))
  for (k in seq_len(OUTER_FOLDS)) {
    test  <- which(folds == k)
    train <- which(folds != k)
    fp    <- fit_predict(train, test)
    prob[test]      <- fp$prob
    prob_base[test] <- baseline_predict(train, test)
    i <- i + 1
    selected[i, ] <- fp$coef
    n_used[i]     <- fp$n_used
  }
  cv_runs[[r]] <- tibble(repetition = r, run = smp$run, y = y,
                         prob = prob, prob_baseline = prob_base)
}
cv_pred <- bind_rows(cv_runs)

cv_auc <- cv_pred |>
  group_by(repetition) |>
  summarise(auc_proteome = auc_of(y, prob),
            auc_baseline = auc_of(y, prob_baseline), .groups = "drop")

# One out-of-fold prediction per sample: mean over repetitions
oof <- cv_pred |>
  group_by(run, y) |>
  summarise(prob = mean(prob), prob_baseline = mean(prob_baseline), .groups = "drop")

roc_p <- roc(oof$y, oof$prob, direction = "<", quiet = TRUE)
roc_b <- roc(oof$y, oof$prob_baseline, direction = "<", quiet = TRUE)
ci_p  <- as.numeric(ci.auc(roc_p, method = "delong"))
ci_b  <- as.numeric(ci.auc(roc_b, method = "delong"))
oof_auc <- tibble(
  model    = c("proteome", "age + sex"),
  auc      = c(ci_p[2], ci_b[2]),
  ci_low   = c(ci_p[1], ci_b[1]),
  ci_high  = c(ci_p[3], ci_b[3]),
  p_vs_baseline = c(roc.test(roc_p, roc_b, method = "delong", paired = TRUE)$p.value, NA)
)

pred_ad <- oof$prob >= THRESHOLD
operating <- tibble(
  threshold   = THRESHOLD,
  sensitivity = mean(pred_ad[oof$y == 1]),
  specificity = mean(!pred_ad[oof$y == 0]),
  ppv         = mean(oof$y[pred_ad] == 1),
  npv         = mean(oof$y[!pred_ad] == 0)
)

# Calibration: observed AD fraction by bin of predicted probability, and
# calibration slope and intercept (ideal: 1 and 0)
calibration <- oof |>
  mutate(bin = cut(prob, breaks = seq(0, 1, length.out = CALIB_BINS + 1),
                   include.lowest = TRUE)) |>
  group_by(bin) |>
  summarise(n = n(), predicted = mean(prob), observed = mean(y), .groups = "drop")
lp <- qlogis(pmin(pmax(oof$prob, 1e-6), 1 - 1e-6))
calib_fit <- tibble(
  slope     = unname(coef(glm(oof$y ~ lp, family = binomial))[2]),
  intercept = unname(coef(glm(oof$y ~ offset(lp), family = binomial))[1])
)

selection <- tibble(
  protein      = colnames(X),
  gene         = genes,
  freq         = colMeans(selected != 0),   # fraction of outer fits selecting it
  mean_coef    = colMeans(selected)
) |>
  arrange(desc(freq), desc(abs(mean_coef)))

# ---- Leave-one-site-out -------------------------------------------------------
loso <- bind_rows(lapply(levels(smp$site), function(s) {
  test  <- which(smp$site == s)
  train <- which(smp$site != s)
  tibble(site = s, run = smp$run[test], y = y[test],
         prob = fit_predict(train, test)$prob,
         prob_baseline = baseline_predict(train, test))
}))

loso_auc <- loso |>
  group_by(site) |>
  summarise(n = n(), n_AD = sum(y),
            auc_proteome = auc_of(y, prob),
            auc_baseline = auc_of(y, prob_baseline), .groups = "drop")

# ---- Checks -------------------------------------------------------------------
stopifnot(
  nrow(smp) == 181L,
  sum(y) == 88L,
  !any(grepl("(^|;)MAPT(;|$)", genes)),
  all(cv_pred$prob >= 0 & cv_pred$prob <= 1),
  nrow(oof) == nrow(smp),
  sum(calibration$n) == nrow(smp)
)

cat("Candidate protein groups:", ncol(X),
    "| used per training split: median", median(n_used),
    "(range", min(n_used), "-", max(n_used), ")\n")
print(summarise(cv_auc, across(starts_with("auc"),
                               list(mean = mean, min = min, max = max))))
print(oof_auc)
print(operating)
print(calibration)
print(calib_fit)
print(loso_auc)
print(head(selection, 15))

# ---- Save ---------------------------------------------------------------------
saveRDS(list(cv_pred = cv_pred, cv_auc = cv_auc, selection = selection,
             oof = oof, oof_auc = oof_auc, operating = operating,
             calibration = calibration, calib_fit = calib_fit,
             loso = loso, loso_auc = loso_auc,
             params = list(alpha = ALPHA, outer_folds = OUTER_FOLDS,
                           outer_repeats = OUTER_REPEATS, inner_folds = INNER_FOLDS,
                           max_missing = MAX_MISSING, excluded_genes = EXCLUDE_GENES,
                           excluded_sites = EXCLUDE_SITES, seed = SEED,
                           threshold = THRESHOLD,
                           n_features = ncol(X), n_used_median = median(n_used))),
        out_file)
