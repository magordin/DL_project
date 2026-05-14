#!/usr/bin/env bash
#BSUB -J dl_eval
#BSUB -q gpua100
#BSUB -gpu "num=1:mode=exclusive_process"
#BSUB -n 4
#BSUB -R "rusage[mem=64GB]"
#BSUB -W 12:00
#BSUB -o /work3/s252608/DL_project/logs/models_eval_%J.out
#BSUB -e /work3/s252608/DL_project/logs/models_eval_%J.err

set -euo pipefail

PROJECT_ROOT="/work3/s252608/DL_project"
SCRIPT_DIR="${PROJECT_ROOT}/scripts"
VENV_PATH="${PROJECT_ROOT}/.venv"

PROCESSED_PATH="${PROJECT_ROOT}/data/processed"
REPRESENTATIONS_PATH="${PROJECT_ROOT}/data/representations"
MODEL_ROOT="${PROJECT_ROOT}/data/output/model"

DATASET_NAME="${DATASET_NAME:-bulk}"
TARGET_MODE="${TARGET_MODE:-absolute}"
LATENT_DIM="${LATENT_DIM:-128}"

Y_TARGET="${PROCESSED_PATH}/${DATASET_NAME}_normalized_y_target_CPM.h5ad"

REPRESENTATIONS="${REPRESENTATIONS:-raw pca vae}"
MODEL_TYPES="${MODEL_TYPES:-mse gaussian nb}"
VAE_BETAS="${VAE_BETAS:-0p1 0p3 0p5 1p0}"

SEEDS="${SEEDS:-1}"
SPLIT_SEED="${SPLIT_SEED:-42}"
BATCH_SIZE="${BATCH_SIZE:-256}"

source "${VENV_PATH}/bin/activate"
PYTHON_BIN="${VENV_PATH}/bin/python"

echo "Target: ${Y_TARGET}"
echo "Split seed: ${SPLIT_SEED}"

for REPR in ${REPRESENTATIONS}; do

  if [[ "${REPR}" == "vae" ]]; then
    REPR_VARIANTS="${VAE_BETAS}"
  else
    REPR_VARIANTS="base"
  fi

  for VARIANT in ${REPR_VARIANTS}; do

    case "${REPR}" in
      raw)
        REPR_FILE="${REPRESENTATIONS_PATH}/${DATASET_NAME}_raw_dim${LATENT_DIM}.npy"
        REPR_TAG="raw"
        ;;
      pca)
        REPR_FILE="${REPRESENTATIONS_PATH}/${DATASET_NAME}_pca_dim${LATENT_DIM}.npy"
        REPR_TAG="pca"
        ;;
      vae)
        REPR_FILE="${REPRESENTATIONS_PATH}/${DATASET_NAME}_vae_dim${LATENT_DIM}_beta${VARIANT}.npy"
        REPR_TAG="vae_beta${VARIANT}"
        ;;
      *)
        echo "ERROR: unknown representation: ${REPR}"
        exit 1
        ;;
    esac

    for MODEL_TYPE in ${MODEL_TYPES}; do
      for SEED in ${SEEDS}; do

        REPR_MODEL_DIR="${MODEL_ROOT}/${REPR_TAG}/${TARGET_MODE}/${MODEL_TYPE}"

        MODEL_FILE="${REPR_MODEL_DIR}/model_seed${SEED}.pt"
        SPLIT_FILE="${REPR_MODEL_DIR}/split_seed${SPLIT_SEED}.json"
        METRICS_FILE="${REPR_MODEL_DIR}/metrics_seed${SEED}.csv"
        PRED_FILE="${REPR_MODEL_DIR}/preds_seed${SEED}.h5ad"

        echo
        echo "---- Evaluating: ${REPR_TAG} | model=${MODEL_TYPE} | seed=${SEED}"

        if [[ ! -f "${MODEL_FILE}" ]]; then
          echo "Skipping: missing model ${MODEL_FILE}"
          continue
        fi

        if [[ ! -f "${SPLIT_FILE}" ]]; then
          echo "Skipping: missing split ${SPLIT_FILE}"
          continue
        fi

        if [[ -f "${METRICS_FILE}" && -f "${PRED_FILE}" ]]; then
          echo "Skipping existing evaluation: ${METRICS_FILE}"
          continue
        fi

        "${PYTHON_BIN}" "${SCRIPT_DIR}/04_evaluate_models.py" \
          --input-repr "${REPR_FILE}" \
          --target-h5ad "${Y_TARGET}" \
          --target-mode "${TARGET_MODE}" \
          --model-path "${MODEL_FILE}" \
          --split-json "${SPLIT_FILE}" \
          --out-metrics "${METRICS_FILE}" \
          --out-preds "${PRED_FILE}" \
          --batch-size "${BATCH_SIZE}" \
          --eval-split test
          

      done
    done
  done
done

echo
echo "Finished evaluation."