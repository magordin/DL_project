#!/usr/bin/env bash
set -euo pipefail

############################################
# run_qc.sh
# Compute gene-level QC summary metrics
############################################

PROJECT_ROOT="/work3/s252608/DL_project"
SCRIPT_DIR="${PROJECT_ROOT}/scripts"
VENV_PATH="${PROJECT_ROOT}/.venv"

RAW_DATA_PATH="${PROJECT_ROOT}/data/mock"
PROCESSED_PATH="${PROJECT_ROOT}/data/qc"

mkdir -p "${PROCESSED_PATH}"

DATASET_NAME="${DATASET_NAME:-bulk_mock}"

RAW_GENE_H5AD="${RAW_DATA_PATH}/${DATASET_NAME}_genes.h5ad"
RAW_TX_H5AD="${RAW_DATA_PATH}/${DATASET_NAME}_transcripts.h5ad"
MAPPING_JSON="${RAW_DATA_PATH}/${DATASET_NAME}_gene_to_transcripts.json"

QC_TABLE="${PROCESSED_PATH}/${DATASET_NAME}_gene_summary_metrics.csv"

if [[ ! -d "${VENV_PATH}" ]]; then
  echo "ERROR: .venv not found at ${VENV_PATH}"
  exit 1
fi

# shellcheck disable=SC1091
source "${VENV_PATH}/bin/activate"
PYTHON_BIN="${VENV_PATH}/bin/python"

for path in "${RAW_GENE_H5AD}" "${RAW_TX_H5AD}" "${_JSON}"; do
  if [[ ! -f "${path}" ]]; then
    echo "ERROR: missing input file -> ${path}"
    exit 1
  fi
done

echo "========================================"
echo "[QC] Computing gene-level summary table"
echo "========================================"

"${PYTHON_BIN}" "${SCRIPT_DIR}/00_compute_qc.py" \
  --gene-h5ad "${RAW_GENE_H5AD}" \
  --tx-h5ad "${RAW_TX_H5AD}" \
  --mapping-json "${MAPPING_JSON}" \
  --out-csv "${QC_TABLE}" \
  --row-block-size 256

echo "Saved QC table: ${QC_TABLE}"