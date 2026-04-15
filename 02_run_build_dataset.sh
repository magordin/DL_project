#!/usr/bin/env bash
set -euo pipefail

############################################
# run_build_dataset.sh
# Build final dataset after QC-based filtering
############################################

PROJECT_ROOT="/work3/s252608/DL_project"
SCRIPT_DIR="${PROJECT_ROOT}/scripts"
VENV_PATH="${PROJECT_ROOT}/.venv"

RAW_DATA_PATH="${PROJECT_ROOT}/data/raw"
PROCESSED_PATH="${PROJECT_ROOT}/data/processed"

mkdir -p "${PROCESSED_PATH}"

DATASET_NAME="${DATASET_NAME:-bulk}"

RAW_GENE_H5AD="${RAW_DATA_PATH}/${DATASET_NAME}_processed_genes.h5ad"
RAW_TX_H5AD="${RAW_DATA_PATH}/${DATASET_NAME}_processed_transcripts.h5ad"
MAPPING_JSON="${PROJECT_ROOT}/data/${DATASET_NAME}_gene_to_transcripts.json"

QC_TABLE="${PROCESSED_PATH}/${DATASET_NAME}_gene_summary_metrics.csv"

N_HVG="${N_HVG:-3000}"

# Basic filtering
MIN_ISOFORMS="${MIN_ISOFORMS:-2}"
MIN_CORR="${MIN_CORR:-0.80}"
MAX_MEDIAN_REL_ERROR="${MAX_MEDIAN_REL_ERROR:-0.50}"

# Triviality filter
MIN_MEAN_ENTROPY="${MIN_MEAN_ENTROPY:-0.05}"
MIN_STD_ENTROPY="${MIN_STD_ENTROPY:-0.01}"

# Modeling target mode
# absolute | proportion
TARGET_MODE="${TARGET_MODE:-absolute}"

# Normalization parameter for gene input only
NORM_TARGET_SUM="${NORM_TARGET_SUM:-10000}"

FILTERED_QC_TABLE="${PROCESSED_PATH}/${DATASET_NAME}_filtered_gene_summary_metrics_hvg${N_HVG}_${TARGET_MODE}.csv"
BENCHMARK_GENES="${PROCESSED_PATH}/${DATASET_NAME}_benchmark_genes_hvg${N_HVG}_${TARGET_MODE}.txt"
X_HVG="${PROCESSED_PATH}/${DATASET_NAME}_x_hvg_${N_HVG}_${TARGET_MODE}.h5ad"
Y_TARGET="${PROCESSED_PATH}/${DATASET_NAME}_y_target_${N_HVG}_${TARGET_MODE}.h5ad"

if [[ ! -d "${VENV_PATH}" ]]; then
  echo "ERROR: .venv not found at ${VENV_PATH}"
  exit 1
fi

# shellcheck disable=SC1091
source "${VENV_PATH}/bin/activate"
PYTHON_BIN="${VENV_PATH}/bin/python"

for path in "${RAW_GENE_H5AD}" "${RAW_TX_H5AD}" "${MAPPING_JSON}" "${QC_TABLE}"; do
  if [[ ! -f "${path}" ]]; then
    echo "ERROR: missing input file -> ${path}"
    exit 1
  fi
done

echo "========================================"
echo "[BUILD DATASET] Building final dataset"
echo "========================================"
echo "Target mode: ${TARGET_MODE}"
echo "N_HVG: ${N_HVG}"

"${PYTHON_BIN}" "${SCRIPT_DIR}/01_build_dataset.py" \
  --gene-h5ad "${RAW_GENE_H5AD}" \
  --tx-h5ad "${RAW_TX_H5AD}" \
  --mapping-json "${MAPPING_JSON}" \
  --qc-csv "${QC_TABLE}" \
  --n-hvg "${N_HVG}" \
  --min-isoforms "${MIN_ISOFORMS}" \
  --min-corr "${MIN_CORR}" \
  --max-median-rel-error "${MAX_MEDIAN_REL_ERROR}" \
  --min-mean-entropy "${MIN_MEAN_ENTROPY}" \
  --min-std-entropy "${MIN_STD_ENTROPY}" \
  --target-mode "${TARGET_MODE}" \
  --norm-target-sum "${NORM_TARGET_SUM}" \
  --out-filtered-qc "${FILTERED_QC_TABLE}" \
  --out-benchmark-genes "${BENCHMARK_GENES}" \
  --out-x-hvg "${X_HVG}" \
  --out-y-target "${Y_TARGET}"

echo "Saved filtered QC   : ${FILTERED_QC_TABLE}"
echo "Saved benchmark set : ${BENCHMARK_GENES}"
echo "Saved X input       : ${X_HVG}"
echo "Saved Y target      : ${Y_TARGET}"