#!/bin/bash

#BSUB -J BulkFormer_Inf
#BSUB -q gpuv100
#BSUB -gpu "num=1:mode=exclusive_process"
#BSUB -n 4
#BSUB -W 02:00
#BSUB -R "rusage[mem=32GB]"
#BSUB -o /work3/s252608/DL_project/logs/bulkformer_%J.out
#BSUB -e /work3/s252608/DL_project/logs/bulkformer_%J.err

module load cuda/12.4
source /work3/s252608/DL_project/.venv/bin/activate

export PROJECT_ROOT="/work3/s252608/DL_project"
export BULKFORMER_BASE="$PROJECT_ROOT/BulkFormer"
export PYTHONPATH="$PROJECT_ROOT:$BULKFORMER_BASE:$PYTHONPATH"
cd $PROJECT_ROOT

python scripts/02_build_representations.py \
    --input-h5ad "$PROJECT_ROOT/data/processed/bulk_normalized_x_input_CPM.h5ad" \
    --repr-type bulkformer \
    --latent-dim 128 \
    --out-file "$PROJECT_ROOT/data/representations/bulkformer_representations.npy"