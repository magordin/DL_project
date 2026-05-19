# DL_project

This repository contains a bulk RNA-seq isoform representation and modeling pipeline. It supports gene/transcript QC, preprocessing, multiple representation extraction strategies (`raw`, `PCA`, `VAE`, `BulkFormer`), model training, evaluation, and plotting.

## Project structure

- `01_run_qc.sh`, `02_run_build_dataset.sh`, `03_run_representations.sh`, `03a_run_vae.sh`, `03b_run_bulkformer.sh`, `04_run_models.sh`, `05_predictions.sh`
  - Shell wrappers for core project stages.
- `scripts/`
  - `00_compute_qc.py` — compute gene-level QC summary metrics from gene and transcript `.h5ad` inputs.
  - `01_preprocess_data.py` — align gene/transcript data and save filtered, normalized AnnData objects.
  - `02_build_representations.py` — generate representations from filtered gene data: raw, PCA, VAE, or BulkFormer.
  - `03_train_models.py` — train prediction models on representation inputs.
  - `04_evaluate_models.py` — evaluate trained models and save metrics and predictions.
  - `05_plot_results.py` — aggregate metrics and generate plots across representation experiments.
- `scripts/src/` — helper modules used by the pipeline.
- `BulkFormer/` — BulkFormer model code, pretrained checkpoint, and related utilities.
- `data/` — raw inputs, processed outputs, embeddings, and helper files.
- `notebooks/` — exploration and evaluation notebooks.
- `results/` — generated results, figures, and summaries.
- `logs/` — execution logs from pipeline runs.
- `requirements.txt` — Python dependencies.

## Workflow overview

The typical pipeline is:

1. **QC**
   - Run `scripts/00_compute_qc.py` to compute gene quality metrics.
   - Input: raw gene `.h5ad`, raw transcript `.h5ad`, mapping JSON.
   - Output: QC summary CSV.

2. **Preprocessing**
   - Run `scripts/01_preprocess_data.py` to subset and normalize retained genes and transcripts.
   - Input: raw gene `.h5ad`, raw transcript `.h5ad`, QC CSV, mapping JSON.
   - Output: processed AnnData files in the specified output directory.

3. **Representation building**
   - Run `scripts/02_build_representations.py` on filtered gene input.
   - Supported `--repr-type`: `raw`, `pca`, `vae`, `bulkformer`.
   - Output: NumPy `.npy` representation file.

4. **Model training**
   - Run `scripts/03_train_models.py` on representation `.npy` data and target transcript `.h5ad` data.
   - Output: checkpoint, training history CSV, split JSON.

5. **Evaluation**
   - Run `scripts/04_evaluate_models.py` with the trained model and split files.
   - Output: metrics CSV and prediction `.h5ad` file.

6. **Plotting**
   - Run `scripts/05_plot_results.py` to collect metrics and generate comparison plots.

## Example commands

```bash
python scripts/00_compute_qc.py \
  --gene-h5ad data/raw/bulk_mock_genes.h5ad \
  --tx-h5ad data/raw/bulk_mock_transcripts.h5ad \
  --mapping-json data/mock/bulk_mock_gene_to_transcripts.json \
  --out-csv data/qc/gene_qc.csv

python scripts/01_preprocess_data.py \
  --gene-h5ad data/raw/bulk_mock_genes.h5ad \
  --tx-h5ad data/raw/bulk_mock_transcripts.h5ad \
  --mapping-json data/mock/bulk_mock_gene_to_transcripts.json \
  --qc-csv data/qc/gene_qc.csv \
  --out-dir data/processed

python scripts/02_build_representations.py \
  --input-h5ad data/processed/bulk_normalized_x_input_CPM.h5ad \
  --repr-type bulkformer \
  --latent-dim 128 \
  --out-file data/representations/bulkformer_embeddings.npy
```

## BulkFormer integration

The `BulkFormer/` folder contains the model implementation, a pretrained checkpoint at `BulkFormer/model/Bulkformer_ckpt_epoch_29.pt`, and utilities used by `scripts/src/bulkformer.py`.

## Dependencies

Install project dependencies with:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip setuptools wheel
pip install -r requirements.txt
```

## Notes

- The repository is designed for reproducible bulk RNA-seq representation and model evaluation.
- `notebooks/` contains interactive analyses and result exploration.
