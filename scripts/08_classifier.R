# 08_classifier.R
# AD vs control classifier from the CSF proteome: elastic-net logistic
# regression, evaluated by repeated nested cross-validation and by
# leave-one-site-out validation, against an age + sex baseline.
#
# Every data-dependent step (protein selection, imputation, scaling, penalty)
# is learned on the training data of each split only.
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
MAX_MISSING   <- 0.2     # protein groups with more missing values are not used
EXCLUDE_GENES <- "MAPT"  # diagnosis is defined from CSF tau
EXCLUDE_SITES <- "Kiel"  # controls only and age-mismatched

dat <- readRDS(in_rds)
use <- !dat$sample$excluded & !(dat$sample$site %in% EXCLUDE_SITES)
smp <- dat$sample[use, ] |> mutate(site = droplevels(site))
y   <- as.integer(smp$group == "AD")

X     <- t(log2(dat$intensity[, use]))                 # runs x protein groups
genes <- dat$protein$PG.Genes
excl  <- grepl(paste0("(^|;)(", paste(EXCLUDE_GENES, collapse = "|"), ")(;|$)"), genes)
feat  <- colMeans(is.na(X)) <= MAX_MISSING & !excl
X     <- X[, feat]
genes <- genes[feat]

# ---- Preprocessing learned on training data only ------------------------------
prep_fit <- function(Xtr) {
  fill <- apply(Xtr, 2, function(v) min(v, na.rm = TRUE) - 1)   # half the minimum
  Xf   <- prep_apply(Xtr, list(fill = fill, mu = 0, sd = 1))
  list(fill = fill, mu = colMeans(Xf), sd = apply(Xf, 2, sd))
}
prep_apply <- function(M, p) {
  for (j in seq_len(ncol(M))) M[is.na(M[, j]), j] <- p$fill[j]
  sweep(sweep(M, 2, p$mu), 2, p$sd, "/")
}

fit_predict <- function(train, test) {
  p  <- prep_fit(X[train, , drop = FALSE])
  cv <- cv.glmnet(prep_apply(X[train, , drop = FALSE], p), y[train],
                  family = "binomial", alpha = ALPHA,
                  nfolds = INNER_FOLDS, type.measure = "auc")
  list(prob = predict(cv, prep_apply(X[test, , drop = FALSE], p),
                      s = "lambda.1se", type = "response")[, 1],
       coef = as.numeric(coef(cv, s = "lambda.1se"))[-1])
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
  }
  cv_runs[[r]] <- tibble(repetition = r, run = smp$run, y = y,
                         prob = prob, prob_baseline = prob_base)
}
cv_pred <- bind_rows(cv_runs)

cv_auc <- cv_pred |>
  group_by(repetition) |>
  summarise(auc_proteome = auc_of(y, prob),
            auc_baseline = auc_of(y, prob_baseline), .groups = "drop")

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
  all(cv_pred$prob >= 0 & cv_pred$prob <= 1)
)

cat("Proteins used as candidate predictors:", ncol(X), "\n")
print(summarise(cv_auc, across(starts_with("auc"),
                               list(mean = mean, min = min, max = max))))
print(loso_auc)
print(head(selection, 15))

# ---- Save ---------------------------------------------------------------------
saveRDS(list(cv_pred = cv_pred, cv_auc = cv_auc, selection = selection,
             loso = loso, loso_auc = loso_auc,
             params = list(alpha = ALPHA, outer_folds = OUTER_FOLDS,
                           outer_repeats = OUTER_REPEATS, inner_folds = INNER_FOLDS,
                           max_missing = MAX_MISSING, excluded_genes = EXCLUDE_GENES,
                           excluded_sites = EXCLUDE_SITES, seed = SEED,
                           n_features = ncol(X))),
        out_file)
