#!/usr/bin/env bash
#BSUB -J repr_vae_beta
#BSUB -n 1
#BSUB -R "rusage[mem=32GB]"
#BSUB -W 08:00
#BSUB -o /work3/s252608/DL_project/logs/vae_beta_%J.out
#BSUB -e /work3/s252608/DL_project/logs/vae_beta_%J.err

set -euo pipefail

PROJECT_ROOT="/work3/s252608/DL_project"
VENV_PATH="${PROJECT_ROOT}/.venv"

PROCESSED_PATH="${PROJECT_ROOT}/data/processed"
REPRESENTATIONS_PATH="${PROJECT_ROOT}/data/representations"
HISTORY_PATH="${PROJECT_ROOT}/results/vae_history"

mkdir -p "${REPRESENTATIONS_PATH}" "${HISTORY_PATH}" "${PROJECT_ROOT}/logs"

DATASET_NAME="${DATASET_NAME:-bulk}"
LATENT_DIM="${LATENT_DIM:-128}"
SEED="${SEED:-1}"
BETAS="${BETAS:-0.1 0.3 0.5 1.0}"

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
echo "[REPRESENTATIONS] Building VAE beta sweep"
echo "========================================"
echo "Input: ${INPUT_X}"
echo "Betas: ${BETAS}"

cd "${PROJECT_ROOT}"

for BETA in ${BETAS}; do
  echo "---- Target Representation: vae | beta=${BETA}"

  BETA_TAG="${BETA//./p}"

  REPR_FILE="${REPRESENTATIONS_PATH}/${DATASET_NAME}_vae_dim${LATENT_DIM}_beta${BETA_TAG}.npy"
  HISTORY_FILE="${HISTORY_PATH}/${DATASET_NAME}_vae_dim${LATENT_DIM}_beta${BETA_TAG}.csv"

  "${PYTHON_BIN}" -m scripts.02_build_representations \
    --input-h5ad "${INPUT_X}" \
    --repr-type vae \
    --latent-dim "${LATENT_DIM}" \
    --seed "${SEED}" \
    --beta "${BETA}" \
    --out-file "${REPR_FILE}" \
    --history-file "${HISTORY_FILE}"

  echo "Saved representation: ${REPR_FILE}"
  echo "Saved history: ${HISTORY_FILE}"
done

echo "VAE beta sweep complete."