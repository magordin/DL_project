#!/usr/bin/env bash
#BSUB -J qc
#BSUB -n 1
#BSUB -R "rusage[mem=32GB]"
#BSUB -W 24:00
#BSUB -o /work3/s252608/DL_project/logs/qc_%J.out
#BSUB -e /work3/s252608/DL_project/logs/qc_%J.err

set -euo pipefail

PROJECT_ROOT="/work3/s252608/DL_project"
SCRIPT_DIR="${PROJECT_ROOT}/scripts"
VENV_PATH="${PROJECT_ROOT}/.venv"
QC_PATH="${PROJECT_ROOT}/data/qc"
RAW_DATA_PATH="${PROJECT_ROOT}/data/raw"
PROCESSED_PATH="${PROJECT_ROOT}/data/processed"

mkdir -p "${PROCESSED_PATH}"

DATASET_NAME="${DATASET_NAME:-bulk_processed}"

RAW_GENE_H5AD="${RAW_DATA_PATH}/${DATASET_NAME}_genes.h5ad"
RAW_TX_H5AD="${RAW_DATA_PATH}/${DATASET_NAME}_transcripts.h5ad"
MAPPING_JSON="${RAW_DATA_PATH}/${DATASET_NAME}_gene_to_transcripts.json"

QC_TABLE="${QC_PATH}/bulk_qc.csv"

source "${VENV_PATH}/bin/activate"
PYTHON_BIN="${VENV_PATH}/bin/python"

for file in "$RAW_GENE_H5AD" "$RAW_TX_H5AD" "$MAPPING_JSON" "$QC_TABLE"; do
    if [ ! -f "$file" ]; then
        echo "ERROR: File not found: $file"
        exit 1
    fi
done

echo "Target mode: ${TARGET_MODE}"
echo "N_HVG: ${N_HVG}"

"${PYTHON_BIN}" "${SCRIPT_DIR}/01_preprocess_data.py" \
  --gene-h5ad "${RAW_GENE_H5AD}" \
  --tx-h5ad "${RAW_TX_H5AD}" \
  --mapping-json "${MAPPING_JSON}" \
  --qc-csv "${QC_TABLE}" \
  --out-dir "${PROCESSED_PATH}" \
  --data-type "bulk" \
  --n-top 3000

echo "Processing complete. Files saved to ${PROCESSED_PATH}"