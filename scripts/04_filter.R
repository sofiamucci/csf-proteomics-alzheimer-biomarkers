# 04_filter.R
# Protein-level filtering: laboratory contaminants and minimum valid values.
#
# Input:  data/processed/03_sample_qc.rds
# Output: data/processed/04_filter.rds
#           list(intensity, protein, sample, protein_filter)
#           intensity / protein: retained protein groups only
#           protein_filter:      all protein groups, with the reason for exclusion

library(dplyr)

in_rds   <- "data/processed/03_sample_qc.rds"
out_file <- "data/processed/04_filter.rds"

# Skin/hair keratins and trypsin (PRSS3 shares peptides with the digestion enzyme)
CONTAMINANT_PATTERN <- "^(KRT[0-9]+[A-Z]?|PRSS1|PRSS3)$"
# A protein group is kept if quantified in at least this fraction of the runs
# of at least one diagnostic group. Requiring it in both groups would discard
# proteins that fall below detection in one group, which are strong candidates.
MIN_VALID_FRAC <- 0.5

dat <- readRDS(in_rds)

keep_run <- !dat$sample$excluded
X        <- dat$intensity[, keep_run]
grp      <- dat$sample$group[keep_run]

# ---- Contaminants -------------------------------------------------------------
gene_list   <- strsplit(dat$protein$PG.Genes, ";")
contaminant <- vapply(gene_list,
                      function(g) any(grepl(CONTAMINANT_PATTERN, g[!is.na(g)])),
                      logical(1))

# ---- Valid values per group ---------------------------------------------------
valid_frac <- sapply(levels(grp), function(l) rowMeans(!is.na(X[, grp == l])))
enough     <- apply(valid_frac >= MIN_VALID_FRAC, 1, any)

protein_filter <- dat$protein |>
  mutate(
    contaminant        = contaminant,
    valid_frac_control = valid_frac[, "control"],
    valid_frac_AD      = valid_frac[, "AD"],
    filter = case_when(
      contaminant ~ "contaminant",
      !enough     ~ "insufficient valid values",
      .default    = "retained"
    )
  )

retained <- protein_filter$filter == "retained"

# ---- Checks -------------------------------------------------------------------
stopifnot(
  sum(contaminant) == 14L,   # 12 keratin groups + PRSS1 + PRSS3
  !any(contaminant & (dat$protein$ery_panel | dat$protein$plasma_panel)),
  identical(rownames(dat$intensity), protein_filter$PG.ProteinGroups)
)

print(table(protein_filter$filter))

# ---- Save ---------------------------------------------------------------------
saveRDS(
  list(intensity      = dat$intensity[retained, ],
       protein        = protein_filter[retained, ],
       sample         = dat$sample,
       protein_filter = protein_filter),
  out_file
)
