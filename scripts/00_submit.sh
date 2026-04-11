#!/bin/sh

#BSUB -J Preprocess_AnnData
#BSUB -n 1
#BSUB -R "rusage[mem=32GB]"
#BSUB -W 02:00
#BSUB -gpu "num=1:mode=exclusive_process"

#BSUB -o results/logs/preprocess_%J.out
#BSUB -e results/logs/preprocess_%J.err


#### DIRECTORIES
## raw data: /work3/s252608/DL_project/data/raw
## normalized HVGs: /work3/s252608/DL_project/data/processed
## representations: /work3/s252608/DL_project/data/representations
## model weights: /work3/s252608/DL_project/data/output/model
## loss plots: /work3/s252608/DL_project/data/output/plots


source /zhome/bf/7/219671/projects/DL_project/.venv/bin/activate
python3 scripts/01_preprocess_data.py