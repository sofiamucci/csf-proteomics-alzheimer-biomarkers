# 02_annotate.R
# Join run-level sample annotation (authors' Table S2) to the imported data.
#
# Input:  data/processed/01_import.rds
#         data/raw/annotationofsamples_AM1.5.11/annotation of samples_AM1.5.11.xlsx
# Output: data/processed/02_annotate.rds
#           list(intensity, protein, sample); sample rows aligned to intensity columns

library(readxl)
library(dplyr)
library(stringr)

in_rds   <- "data/processed/01_import.rds"
ann_file <- "data/raw/annotationofsamples_AM1.5.11/annotation of samples_AM1.5.11.xlsx"
out_file <- "data/processed/02_annotate.rds"

dat <- readRDS(in_rds)

# ---- Sample annotation --------------------------------------------------------
# Runs excluded by the authors are kept and flagged, so they can be inspected in QC
sample <- read_excel(ann_file) |>
  filter(!is.na(`sample name`)) |>
  transmute(
    run              = str_remove(`sample name`, "\\.raw\\.PG\\.Quantity$"),
    site             = factor(`collection site`,
                              levels = c("Sweden", "Magdeburg", "Kiel", "Berlin")),
    batch            = factor(str_sub(run, 1, 8)),   # date of the MS run series
    group            = factor(case_match(`biochemical AD classification`,
                                         "biochemical control" ~ "control",
                                         "biochemical AD"      ~ "AD"),
                              levels = c("control", "AD")),
    clinical_dx      = `clinical AD diagnosis`,
    age              = age,
    sex              = factor(gender, levels = c("f", "m")),
    t_tau            = `t-tau [ng/L]`,
    p_tau            = `p-tau [ng/L]`,
    abeta42          = `Abeta-42 [ng/L]`,
    abeta40          = `Abeta-40 [ng/L]`,
    abeta_ratio      = `Abeta ratio`,
    mmse             = `MMSE score`,
    cohort_id        = `cohort sample ID`,
    excluded         = !is.na(comment),
    exclusion_reason = comment
  )

sample <- sample[match(colnames(dat$intensity), sample$run), ]
included <- !sample$excluded

# ---- Checks -------------------------------------------------------------------
site_batch <- table(sample$site, sample$batch)
stopifnot(
  nrow(sample) == 210L,
  identical(sample$run, colnames(dat$intensity)),
  sum(included) == 197L,
  !anyNA(sample$group[included]),
  !anyNA(sample$age[included]),
  !anyNA(sample$sex[included]),
  all(rowSums(site_batch > 0) == 1)      # each site was measured in a single batch
)

print(table(sample$site[included], sample$group[included]))

# ---- Save ---------------------------------------------------------------------
saveRDS(c(dat, list(sample = sample)), out_file)

