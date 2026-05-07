#!/usr/bin/env bash
#BSUB -J repr_vae
#BSUB -n 1
#BSUB -R "rusage[mem=32GB]"
#BSUB -W 08:00
#BSUB -o /work3/s252608/DL_project/logs/qc_%J.out
#BSUB -e /work3/s252608/DL_project/logs/qc_%J.err

set -euo pipefail

PROJECT_ROOT="/work3/s252608/DL_project"
SCRIPT_DIR="${PROJECT_ROOT}/scripts"
VENV_PATH="${PROJECT_ROOT}/.venv"

PROCESSED_PATH="${PROJECT_ROOT}/data/processed"
REPRESENTATIONS_PATH="${PROJECT_ROOT}/data/representations"

mkdir -p "${REPRESENTATIONS_PATH}"

DATASET_NAME="${DATASET_NAME:-bulk}"
LATENT_DIM="${LATENT_DIM:-128}"
SEED="${SEED:-1}"
REPRESENTATIONS="${REPRESENTATIONS:-vae}" 

INPUT_X="${PROCESSED_PATH}/${DATASET_NAME}_normalized_x_input.h5ad"

if [[ ! -d "${VENV_PATH}" ]]; then
  echo "ERROR: .venv not found at ${VENV_PATH}"
  exit 1
fi

source "${VENV_PATH}/bin/activate"
PYTHON_BIN="${VENV_PATH}/bin/python"

if [[ ! -f "${INPUT_X}" ]]; then
  echo "ERROR: missing X input file -> ${INPUT_X}"
  exit 1
fi

echo "========================================"
echo "[REPRESENTATIONS] Building ${REPRESENTATIONS}"
echo "========================================"
echo "Input: ${INPUT_X}"

cd "${PROJECT_ROOT}"

for REPR in ${REPRESENTATIONS}; do
  echo "---- Target Representation: ${REPR}"

  REPR_DIR="${REPRESENTATIONS_PATH}"
  mkdir -p "${REPR_DIR}"

  REPR_FILE="${REPR_DIR}/${DATASET_NAME}_${REPR}_dim${LATENT_DIM}.npy"

  "${PYTHON_BIN}" -m scripts.02_build_representations \
    --input-h5ad "${INPUT_X}" \
    --repr-type "${REPR}" \
    --latent-dim "${LATENT_DIM}" \
    --seed "${SEED}" \
    --out-file "${REPR_FILE}"

  echo "Saved: ${REPR_FILE}"
done

echo "Representation building complete."