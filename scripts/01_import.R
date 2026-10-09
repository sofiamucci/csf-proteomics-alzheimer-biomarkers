# 01_import.R
# Import the Spectronaut protein group report (PXD016278).
#
# Input:  data/raw/proteingroupsreport_usedforanalysis/
#           protein groups report_used for analysis.csv
# Output: data/processed/01_import.rds
#           list(intensity = protein x run matrix (raw scale, NA = not quantified),
#                protein   = protein annotation tibble)

library(readr)
library(dplyr)
library(stringr)

in_file  <- "data/raw/proteingroupsreport_usedforanalysis/protein groups report_used for analysis.csv"
out_file <- "data/processed/01_import.rds"

# ---- Read ---------------------------------------------------------------------
# Spectronaut reports values below the q-value cutoff as "Filtered"
pg_raw <- read_csv(in_file, na = c("", "NA", "Filtered"), show_col_types = FALSE)

# ---- Split annotation and quantities -----------------------------------------
protein <- select(pg_raw, starts_with("PG."))
quant   <- select(pg_raw, ends_with("PG.Quantity")) # 1542 3
stopifnot(ncol(protein) + ncol(quant) == ncol(pg_raw)) # 1542 210

# ---- Run names ----------------------------------------------------------------
run_names <- names(quant) |>
  str_remove("^\\[\\d+\\]\\s*") |>
  str_remove("\\.raw\\.PG\\.Quantity$")

# ---- Intensity matrix ---------------------------------------------------------
# PG.ProteinGroups is the row identifier: PG.Genes has missing and repeated values
intensity <- as.matrix(quant)
dimnames(intensity) <- list(protein$PG.ProteinGroups, run_names)

# ---- Checks -------------------------------------------------------------------
stopifnot(
  identical(dim(intensity), c(1542L, 210L)),
  !anyDuplicated(rownames(intensity)),
  !anyDuplicated(colnames(intensity)),
  is.numeric(intensity),
  all(intensity > 0, na.rm = TRUE)
)

# ---- Save ---------------------------------------------------------------------
dir.create(dirname(out_file), showWarnings = FALSE, recursive = TRUE)
saveRDS(list(intensity = intensity, protein = protein), out_file)
