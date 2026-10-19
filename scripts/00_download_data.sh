#!/usr/bin/env bash
# Download input files for PXD016278 (Bader et al. 2020, Mol Syst Biol).
# Only the protein group report and sample annotation are needed; raw files,
# spectral libraries and Spectronaut sessions are not downloaded.
#
# Usage (from the repository root): bash scripts/00_download_data.sh

set -euo pipefail

BASE_URL="https://ftp.pride.ebi.ac.uk/pride/data/archive/2020/05/PXD016278"
RAW_DIR="data/raw"
LOG_FILE="${RAW_DIR}/DOWNLOAD_LOG.txt"

FILES=(
  "proteingroupsreport_usedforanalysis.zip"
  "annotationofsamples_AM1.5.11.zip"
  "PXD016278_community_annotated.sdrf.tsv"
)

mkdir -p "${RAW_DIR}"

log() {
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" | tee -a "${LOG_FILE}"
}

log "=== Download started (PXD016278) ==="

for f in "${FILES[@]}"; do
  dest="${RAW_DIR}/${f}"

  if [[ -s "${dest}" ]]; then
    log "SKIP     ${f} (already present)"
  else
    log "DOWNLOAD ${f}"
    curl -fL --retry 3 -C - -o "${dest}" "${BASE_URL}/${f}"
  fi

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

  log "SHA256   $(sha256sum "${dest}" | cut -d' ' -f1)  ${f}  ($(du -h "${dest}" | cut -f1))"
done

# Verify the analysis inputs against the checksums recorded at first download.
# The SDRF is community-maintained and may be revised, so it is only logged.
CHECKSUMS="data/PXD016278.sha256"
if (cd "${RAW_DIR}" && sha256sum --quiet -c "../$(basename "${CHECKSUMS}")"); then
  log "OK       checksums match ${CHECKSUMS}"
else
  log "ERROR    checksum mismatch against ${CHECKSUMS}; the PRIDE files differ from those analysed"
  exit 1
fi

log "=== Download finished ==="
