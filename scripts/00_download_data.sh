#!/usr/bin/env bash
# =============================================================================
# 00_download_data.sh
# Download the input files for this analysis from PRIDE (PXD016278).
#
# Bader et al. (2020) Mol Syst Biol 16:e9356. CSF proteome in Alzheimer's
# disease, three cohorts, DIA acquisition processed with Spectronaut.
#
# Only the protein-level report and the sample annotation are downloaded.
# Raw files, spectral libraries and Spectronaut sessions are not needed,
# because the analysis starts from the protein group report.
#
# Usage (from the repository root):
#   bash scripts/00_download_data.sh
# =============================================================================

set -euo pipefail   # stop on errors, undefined variables and failed pipes

BASE_URL="https://ftp.pride.ebi.ac.uk/pride/data/archive/2020/05/PXD016278"
RAW_DIR="data/raw"
LOG_FILE="${RAW_DIR}/DOWNLOAD_LOG.txt"

FILES=(
  "proteingroupsreport_usedforanalysis.zip"   # protein-level quantification
  "annotationofsamples_AM1.5.11.zip"          # sample annotation from the authors
  "PXD016278_community_annotated.sdrf.tsv"    # SDRF: standardized sample metadata
)

mkdir -p "${RAW_DIR}"

log() {
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" | tee -a "${LOG_FILE}"
}

log "=== Download started (PXD016278) ==="

for f in "${FILES[@]}"; do
  dest="${RAW_DIR}/${f}"

  # Idempotent: if the file is already here, do not download it again
  if [[ -s "${dest}" ]]; then
    log "SKIP     ${f} (already present)"
  else
    log "DOWNLOAD ${f}"
    # -f: fail on HTTP errors, -L: follow redirects,
    # --retry: retry on network errors, -C -: resume partial downloads
    curl -fL --retry 3 -C - -o "${dest}" "${BASE_URL}/${f}"
  fi

  # Integrity check: a truncated zip fails this test
  if [[ "${f}" == *.zip ]]; then
    if unzip -tq "${dest}" > /dev/null; then
      log "OK       ${f} (zip integrity test passed)"
      unzip -oq "${dest}" -d "${RAW_DIR}/${f%.zip}"
      log "UNZIPPED ${f} -> ${RAW_DIR}/${f%.zip}/"
    else
      log "ERROR    ${f} failed the zip integrity test; delete it and rerun"
      exit 1
    fi
  fi

  # Record size and checksum so any later change to the file is detectable
  log "SHA256   $(sha256sum "${dest}" | cut -d' ' -f1)  ${f}  ($(du -h "${dest}" | cut -f1))"
done

log "=== Download finished ==="
