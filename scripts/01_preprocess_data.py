from scripts.src.preprocessing import *
import argparse
import logging
from pathlib import Path
import pandas as pd
import anndata as ad

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')


def parse_args():
    parser = argparse.ArgumentParser(description="Preprocess Gene and Transcript data.")
    parser.add_argument("--gene-h5ad", type=str, required=True, help="Path to raw gene h5ad")
    parser.add_argument("--tx-h5ad", type=str, required=True, help="Path to raw transcript h5ad")
    parser.add_argument("--mapping-json", type=str, required=True, help="Path to gene-to-transcript JSON")
    parser.add_argument("--qc-csv", type=str, required=True, help="Path to the computed QC metrics CSV")
    parser.add_argument("--out-dir", type=str, required=True, help="Directory to save processed data")
    
    parser.add_argument("--data-type", type=str, default="bulk", choices=["bulk", "sc"])
    return parser.parse_args()


def subset_and_align_data(gene_adata, iso_adata, qc_csv_path, mapping):
    qc_results = pd.read_csv(qc_csv_path)
    retained_genes = qc_results[qc_results['qc_pass_relaxed']]['gene_id'].to_list()

    valid_genes = [g for g in retained_genes if g in gene_adata.var_names]
    logging.info(f"Retaining {len(valid_genes)} genes after QC filtering.")

    gene_sub = gene_adata[:, valid_genes].to_memory()
    valid_gene_set = set(valid_genes)
    target_iso_ids = [t for g in valid_genes if g in mapping for t in mapping[g]]
    final_iso_ids = [t for t in target_iso_ids if t in iso_adata.var_names]
    iso_sub = iso_adata[gene_sub.obs_names, final_iso_ids].to_memory()
    
    return gene_sub, iso_sub


def main():
    args = parse_args()
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    logging.info("Loading h5ad files...")
    raw_genes = ad.read_h5ad(args.gene_h5ad, backed='r')
    raw_iso = ad.read_h5ad(args.tx_h5ad, backed='r')
    mapping = get_gene_mapping(args.mapping_json, raw_genes)
    gene_sub, iso_sub = subset_and_align_data(raw_genes, raw_iso, args.qc_csv, mapping)

    gene_sub = normalize(gene_sub)
    iso_sub = get_isoform_proportions(iso_sub, mapping)

    save_adata(gene_sub, out_dir, args.data_type, "x_input")
    save_adata(iso_sub, out_dir, args.data_type, "y_target")
    
    logging.info(f"Preprocessing for {args.data_type} complete.")


if __name__ == "__main__":
    main()