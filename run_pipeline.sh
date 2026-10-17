#!/usr/bin/env bash
# Run the full analysis, from data download to rendered reports.
#
# Usage: bash run_pipeline.sh
# Requires the csf-proteomics conda environment to be active. Quarto is taken
# from PATH, or from RStudio's bundled copy; override with QUARTO=/path/to/quarto.

set -euo pipefail
cd "$(dirname "$0")"

LOG_DIR="logs"
mkdir -p "$LOG_DIR"
LOG_FILE="${LOG_DIR}/pipeline_$(date '+%Y%m%d_%H%M%S').log"
exec > >(tee -a "$LOG_FILE") 2>&1

QUARTO="${QUARTO:-$(command -v quarto || echo /usr/lib/rstudio/resources/app/bin/quarto/bin/quarto)}"

log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"; }

run() {
  log ">>> $*"
  local t0=$SECONDS
  "$@"
  log "<<< done in $((SECONDS - t0)) s"
}

# ---- Provenance ---------------------------------------------------------------
[[ -x "$QUARTO" ]] || { log "ERROR: quarto not found; set QUARTO=/path/to/quarto"; exit 1; }

log "=== Pipeline started ==="
log "Git commit: $(git rev-parse --short HEAD 2>/dev/null || echo 'not a git repository')"
if [[ -n "$(git status --porcelain 2>/dev/null)" ]]; then
  log "WARNING: uncommitted changes; results may not correspond to the commit above"
fi
log "R:      $(Rscript -e 'cat(R.version.string)')"
log "Quarto: $("$QUARTO" --version)"

# ---- Pipeline -----------------------------------------------------------------
run bash    scripts/00_download_data.sh
run Rscript scripts/01_import.R
run Rscript scripts/02_annotate.R
run Rscript scripts/03_sample_qc.R
run Rscript scripts/04_filter.R
run Rscript scripts/05_missing_value_methods.R
run Rscript scripts/06_differential_abundance.R
run Rscript scripts/07_export_results.R
run Rscript scripts/08_classifier.R
run "$QUARTO" render reports/qc_report.qmd
run "$QUARTO" render reports/analysis_report.qmd

# ---- Publish ------------------------------------------------------------------
# docs/ is served by GitHub Pages
mkdir -p docs
run cp reports/qc_report.html reports/analysis_report.html docs/

log "=== Pipeline finished (log: ${LOG_FILE}) ==="
