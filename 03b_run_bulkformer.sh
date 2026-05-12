#!/bin/bash

#BSUB -J BulkFormer_Inf
#BSUB -q gpuv100
#BSUB -gpu "num=1:mode=exclusive_process"
#BSUB -n 4
#BSUB -W 02:00
#BSUB -R "rusage[mem=32GB]"
#BSUB -o /zhome/bf/7/219671/projects/DL_project/results/logs/bulkformer_%J.out
#BSUB -e /zhome/bf/7/219671/projects/DL_project/results/logs/bulkformer_%J.err

module load cuda/12.4
source /work3/s252608/DL_project/.venv/bin/activate

export PYTHONPATH="/work3/s252608/DL_project/BulkFormer:$PYTHONPATH"
cd /work3/s252608/DL_project/BulkFormer
python /work3/s252608/DL_project/BulkFormer/utils/inference_final.py