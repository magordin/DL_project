from scripts.src.preprocessing import *
import argparse
from pathlib import Path
import pandas as pd
import anndata as ad


def parse_args():
    parser = argparse.ArgumentParser(description="Preprocess Gene and Transcript data.")
    parser.add_argument("--gene-h5ad", type=str, required=True, help="Path to raw gene h5ad")
    parser.add_argument("--tx-h5ad", type=str, required=True, help="Path to raw transcript h5ad")
    parser.add_argument("--mapping-json", type=str, required=True, help="Path to gene-to-transcript JSON")
    parser.add_argument("--qc-csv", type=str, required=True, help="Path to the computed QC metrics CSV")
    parser.add_argument("--out-dir", type=str, required=True, help="Directory to save processed data")
    
    parser.add_argument("--data-type", type=str, default="bulk", choices=["bulk", "sc"])
    parser.add_argument("--n-top", type=int, default=3000, help="Number of HVGs to select")
    return parser.parse_args()


def main():
    args = parse_args()
    out_path = Path(args.out_dir)
    out_path.mkdir(parents=True, exist_ok=True)

    raw_genes_bulk = ad.read_h5ad(args.gene_h5ad, backed='r')
    raw_iso_bulk = ad.read_h5ad(args.tx_h5ad, backed='r')

    mapping = get_gene_mapping(args.mapping_json)

    qc_results = pd.read_csv(args.qc_csv)
    retained_genes = qc_results['gene_id'].to_list()

    valid_retained_genes = [g for g in retained_genes if g in raw_genes_bulk.var_names]
    raw_genes_bulk_sub = raw_genes_bulk[:, valid_retained_genes].to_memory()

    adata_bulk_hvg = select_hvgs(raw_genes_bulk_sub, n_top=args.n_top)
    hvg_list = adata_bulk_hvg.var_names.tolist()

    adata_bulk_input = normalize(adata_bulk_hvg.copy())
    save_adata(adata_bulk_input, args.out_dir, args.data_type, "x_input")

    target_iso_ids = [t for g in hvg_list if g in mapping for t in mapping[g]]
    final_iso_ids = [t for t in target_iso_ids if t in raw_iso_bulk.var_names]
    
    adata_iso_target = raw_iso_bulk[adata_bulk_input.obs_names, final_iso_ids].to_memory()
    adata_iso_target.obs['library_size'] = adata_iso_target.X.sum(axis=1).A1

    save_adata(adata_iso_target, args.out_dir, args.data_type, "y_target")
    print(f"{args.data_type.capitalize()} processing complete.")


if __name__ == "__main__":
    main()