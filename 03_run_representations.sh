#!/usr/bin/env bash
set -euo pipefail

############################################
# run_representations.sh
# Build representations from final X input
############################################

PROJECT_ROOT="/work3/s252608/DL_project"
SCRIPT_DIR="${PROJECT_ROOT}/scripts"
VENV_PATH="${PROJECT_ROOT}/.venv"

PROCESSED_PATH="${PROJECT_ROOT}/data/processed"
REPRESENTATIONS_PATH="${PROJECT_ROOT}/data/representations"

mkdir -p "${REPRESENTATIONS_PATH}"

DATASET_NAME="${DATASET_NAME:-bulk}"
N_HVG="${N_HVG:-3000}"
TARGET_MODE="${TARGET_MODE:-absolute}"

X_HVG="${PROCESSED_PATH}/${DATASET_NAME}_x_hvg_${N_HVG}_${TARGET_MODE}.h5ad"

REPRESENTATIONS="${REPRESENTATIONS:-hvg pca vae geneformer}"
LATENT_DIM="${LATENT_DIM:-128}"
SEED="${SEED:-1}"

if [[ ! -d "${VENV_PATH}" ]]; then
  echo "ERROR: .venv not found at ${VENV_PATH}"
  exit 1
fi

# shellcheck disable=SC1091
source "${VENV_PATH}/bin/activate"
PYTHON_BIN="${VENV_PATH}/bin/python"

if [[ ! -f "${X_HVG}" ]]; then
  echo "ERROR: missing X input file -> ${X_HVG}"
  exit 1
fi

echo "========================================"
echo "[REPRESENTATIONS] Building representations"
echo "========================================"
echo "Input: ${X_HVG}"

for REPR in ${REPRESENTATIONS}; do
  echo "---- Representation: ${REPR}"

  REPR_DIR="${REPRESENTATIONS_PATH}/${REPR}"
  mkdir -p "${REPR_DIR}"

  REPR_FILE="${REPR_DIR}/${DATASET_NAME}_${REPR}_hvg${N_HVG}_${TARGET_MODE}_dim${LATENT_DIM}.npy"

  "${PYTHON_BIN}" "${SCRIPT_DIR}/02_build_representations.py" \
    --input-h5ad "${X_HVG}" \
    --repr-type "${REPR}" \
    --latent-dim "${LATENT_DIM}" \
    --seed "${SEED}" \
    --out-file "${REPR_FILE}"

  echo "Saved: ${REPR_FILE}"
done