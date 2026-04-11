from scripts.src.preprocessing import *

def main():
    processed_data.mkdir(parents=True, exist_ok=True)

    raw_genes_bulk = ad.read_h5ad(raw_data / "bulk_processed_genes.h5ad", backed='r')
    raw_iso_bulk = ad.read_h5ad(raw_data / "bulk_processed_transcripts.h5ad", backed='r')
    mapping = get_gene_mapping("bulk")

    qc_results = pd.read_csv(processed_data / "gene_summary_metrics.csv")

    adata_bulk_hvg = select_hvgs(raw_genes_bulk.to_memory(), qc_results, n_top=3000)
    hvg_list = adata_bulk_hvg.var_names.tolist()

    adata_bulk_input = normalize(adata_bulk_hvg)
    save_adata(adata_bulk_input, "bulk", "x_input")

    target_iso_ids = [t for g in hvg_list if g in mapping for t in mapping[g]]
    final_iso_ids = [t for t in target_iso_ids if t in raw_iso_bulk.var_names]
    
    adata_iso_target = raw_iso_bulk[adata_bulk_input.obs_names, final_iso_ids].to_memory()
    # adata_iso_target = prep_raw_target(adata_iso_target)
    save_adata(adata_iso_target, "bulk", "y_target")

    print("Bulk processing complete.")