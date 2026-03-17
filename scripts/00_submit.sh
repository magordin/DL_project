#!/bin/sh

#BSUB -J Preprocess_AnnData
#BSUB -n 1
#BSUB -R "rusage[mem=32GB]"
#BSUB -W 02:00
#BSUB -gpu "num=1:mode=exclusive_process"

#BSUB -o results/logs/preprocess_%J.out
#BSUB -e results/logs/preprocess_%J.err

source /zhome/bf/7/219671/projects/DL_project/.venv/bin/activate
python3 scripts/01_preprocess_data.py