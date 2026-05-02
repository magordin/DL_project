#!/usr/bin/env bash
#BSUB -J qc
#BSUB -n 1
#BSUB -R "rusage[mem=48GB]"
#BSUB -W 24:00
#BSUB -o /work3/s252608/DL_project/logs/qc_%J.out
#BSUB -e /work3/s252608/DL_project/logs/qc_%J.err

set -euo pipefail

PROJECT_ROOT="/work3/s252608/DL_project"
SCRIPT_DIR="${PROJECT_ROOT}/scripts"
VENV_PATH="${PROJECT_ROOT}/.venv"

RAW_DATA_PATH="${PROJECT_ROOT}/data/raw"
QC_PATH="${PROJECT_ROOT}/data/qc"
LOG_DIR="${PROJECT_ROOT}/logs"

mkdir -p "${QC_PATH}" "${LOG_DIR}"

DATASET_NAME="${DATASET_NAME:-bulk_processed}"

RAW_GENE_H5AD="${RAW_DATA_PATH}/${DATASET_NAME}_genes.h5ad"
RAW_TX_H5AD="${RAW_DATA_PATH}/${DATASET_NAME}_transcripts.h5ad"
MAPPING_JSON="${RAW_DATA_PATH}/${DATASET_NAME}_gene_to_transcripts.json"

QC_TABLE="${QC_PATH}/bulk_qc.csv"

source "${VENV_PATH}/bin/activate"
PYTHON_BIN="${VENV_PATH}/bin/python"

cd "${PROJECT_ROOT}"

PYTHONPATH="${PROJECT_ROOT}" "${PYTHON_BIN}" "${SCRIPT_DIR}/00_compute_qc.py" \
  --gene-h5ad "${RAW_GENE_H5AD}" \
  --tx-h5ad "${RAW_TX_H5AD}" \
  --mapping-json "${MAPPING_JSON}" \
  --out-csv "${QC_TABLE}" \
  --gene-batch-size 500 \
  --row-block-size 512 \
  --checkpoint-every 500