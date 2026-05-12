#!/usr/bin/env bash
set -euo pipefail

PROJECT_ROOT="/work3/s252608/DL_project"
SCRIPT_DIR="${PROJECT_ROOT}/scripts"
VENV_PATH="${PROJECT_ROOT}/.venv"

PROCESSED_PATH="${PROJECT_ROOT}/data/processed"
REPRESENTATIONS_PATH="${PROJECT_ROOT}/data/representations"
MODEL_ROOT="${PROJECT_ROOT}/data/output/model"
PLOTS_PATH="${PROJECT_ROOT}/data/output/plots"

mkdir -p "${MODEL_ROOT}" "${PLOTS_PATH}"

DATASET_NAME="${DATASET_NAME:-bulk}"
TARGET_MODE="${TARGET_MODE:-absolute}"
LATENT_DIM="${LATENT_DIM:-128}"

Y_TARGET="${PROCESSED_PATH}/${DATASET_NAME}_normalized_y_target.h5ad"

REPRESENTATIONS="${REPRESENTATIONS:-raw pca vae}"
SEEDS="${SEEDS:-1 2 3}"

EPOCHS="${EPOCHS:-100}"
BATCH_SIZE="${BATCH_SIZE:-128}"
LEARNING_RATE="${LEARNING_RATE:-1e-3}"
HIDDEN_DIM="${HIDDEN_DIM:-512}"
DROPOUT="${DROPOUT:-0.2}"
WEIGHT_DECAY="${WEIGHT_DECAY:-1e-5}"

RUN_PLOTS="${RUN_PLOTS:-1}"

if [[ ! -d "${VENV_PATH}" ]]; then
  echo "ERROR: .venv not found at ${VENV_PATH}"
  exit 1
fi

source "${VENV_PATH}/bin/activate"
PYTHON_BIN="${VENV_PATH}/bin/python"

if [[ ! -f "${Y_TARGET}" ]]; then
  echo "ERROR: missing target file -> ${Y_TARGET}"
  exit 1
fi

echo "Target: ${Y_TARGET}"

for REPR in ${REPRESENTATIONS}; do
  case "${REPR}" in
    raw)
      REPR_FILE="${REPRESENTATIONS_PATH}/${DATASET_NAME}_raw_dim${LATENT_DIM}.npy"
      ;;
    pca)
      REPR_FILE="${REPRESENTATIONS_PATH}/${DATASET_NAME}_pca_dim${LATENT_DIM}.npy"
      ;;
    vae)
      REPR_FILE="${REPRESENTATIONS_PATH}/${DATASET_NAME}_vae_dim${LATENT_DIM}.npy"
      ;;
    *)
      echo "ERROR: unknown representation: ${REPR}"
      exit 1
      ;;
  esac

  if [[ ! -f "${REPR_FILE}" ]]; then
    echo "ERROR: missing representation file -> ${REPR_FILE}"
    exit 1
  fi

  for SEED in ${SEEDS}; do
    echo
    echo "---- Representation: ${REPR} | seed=${SEED}"

    REPR_MODEL_DIR="${MODEL_ROOT}/${REPR}/${TARGET_MODE}"
    mkdir -p "${REPR_MODEL_DIR}"

    MODEL_FILE="${REPR_MODEL_DIR}/model_seed${SEED}.pt"
    HISTORY_FILE="${REPR_MODEL_DIR}/history_seed${SEED}.csv"
    SPLIT_FILE="${REPR_MODEL_DIR}/split_seed${SEED}.json"
    METRICS_FILE="${REPR_MODEL_DIR}/metrics_seed${SEED}.csv"
    PRED_FILE="${REPR_MODEL_DIR}/preds_seed${SEED}.h5ad"

    "${PYTHON_BIN}" "${SCRIPT_DIR}/03_train_models.py" \
      --input-repr "${REPR_FILE}" \
      --target-h5ad "${Y_TARGET}" \
      --target-mode "${TARGET_MODE}" \
      --seed "${SEED}" \
      --epochs "${EPOCHS}" \
      --batch-size "${BATCH_SIZE}" \
      --learning-rate "${LEARNING_RATE}" \
      --hidden-dim "${HIDDEN_DIM}" \
      --dropout "${DROPOUT}" \
      --weight-decay "${WEIGHT_DECAY}" \
      --out-model "${MODEL_FILE}" \
      --out-history "${HISTORY_FILE}" \
      --out-split "${SPLIT_FILE}"

    "${PYTHON_BIN}" "${SCRIPT_DIR}/04_evaluate_models.py" \
      --input-repr "${REPR_FILE}" \
      --target-h5ad "${Y_TARGET}" \
      --target-mode "${TARGET_MODE}" \
      --model-path "${MODEL_FILE}" \
      --split-json "${SPLIT_FILE}" \
      --out-metrics "${METRICS_FILE}" \
      --out-preds "${PRED_FILE}"
  done
done

if [[ "${RUN_PLOTS}" == "1" ]]; then
  "${PYTHON_BIN}" "${SCRIPT_DIR}/05_plot_results.py" \
    --model-root "${MODEL_ROOT}" \
    --plots-dir "${PLOTS_PATH}" \
    --representations ${REPRESENTATIONS} \
    --target-mode "${TARGET_MODE}"
fi

echo
echo "Finished training/evaluation."