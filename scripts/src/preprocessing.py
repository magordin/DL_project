import scanpy as sc
import json
import numpy as np
from pathlib import Path
import anndata as ad
import pandas as pd

#### these will be given by the command, CHANGE!
root_dir = Path(__file__).resolve().parent.parent
data_dir = Path("/work3/s252608/DL_project/data")

raw_data = data_dir / "raw"
processed_data = data_dir / "processed"


class NumpyEncoder(json.JSONEncoder):
    """ json encoder for numpy types """
    def default(self, obj):
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        if isinstance(obj, (np.int64, np.int32, np.int8)):
            return int(obj)
        if isinstance(obj, (np.float64, np.float32)):
            return float(obj)
        return json.JSONEncoder.default(self, obj)


def normalize(adata, target_sum = 1e4):
    adata.layers["raw_counts"] = adata.X.copy() 
    sc.pp.filter_genes(adata, min_cells=10)
    sc.pp.normalize_total(adata, target_sum=target_sum)
    sc.pp.log1p(adata)
    return adata


def gene_qc_table(adata_gene, adata_iso, gene_to_iso):
    results = []
    for gene, isoforms in gene_to_iso.items():
        if gene not in adata_gene.var_names:
            continue
        g_expr = adata_gene[:, gene].X.toarray().flatten() if hasattr(adata_gene.X, "toarray") else adata_gene[:, gene].X.flatten()
        valid_iso = [i for i in isoforms if i in adata_iso.var_names]
        if len(valid_iso) < 2:
            continue

    iso_matrix = adata_iso[:, valid_iso].X
    sum_iso_expr = iso_matrix.sum(axis=1).A1 if hasattr(iso_matrix, "toarray") else iso_matrix.sum(axis=1)
    corr = np.corrcoef(g_expr, sum_iso_expr)[0, 1] if np.std(sum_iso_expr) > 0 else 0
    iso_fractions = iso_matrix.toarray() / (sum_iso_expr[:, None] + 1e-6)
    mean_fractions = iso_fractions.mean(axis=0)
    dom_iso_frac = np.max(mean_fractions)
    rel_error = np.abs(g_expr - sum_iso_expr) / (g_expr + 1e-6)
        
    results.append({
        'gene_id': gene,
        'n_isoforms': len(valid_iso),
        'mean_gene_expr': np.mean(g_expr),
        'corr_gene_vs_sum_iso': corr,
        'median_rel_error': np.median(rel_error),
        'p90_rel_error': np.percentile(rel_error, 90),
        'dom_iso_frac': dom_iso_frac,
        'fraction_detected': np.mean(g_expr > 0)
    })

    qc_df = pd.DataFrame(results)
    qc_save_path = processed_data / "gene_summary_metrics.csv"
    qc_df.to_csv(qc_save_path, index=False)
    print('Summary table saved successfully.')


def select_hvgs(adata, qc_df, n_top=3000):
    """
    Combines benchmark filtering with HVG selection.
    """
    mask = (
        (qc_df['n_isoforms'] > 1) & 
        (qc_df['dom_iso_frac'] < 0.99) & 
        (qc_df['corr_gene_vs_sum_iso'] > 0.8) & 
        (qc_df['median_rel_error'] < 0.5)
    )
    
    benchmark_gene_ids = qc_df.loc[mask, 'gene_id'].tolist()
    valid_benchmark_genes = [g for g in benchmark_gene_ids if g in adata.var_names]
    
    print(f"Benchmark filtering complete: {len(valid_benchmark_genes)} genes passed QC.")
    adata_benchmark = adata[:, valid_benchmark_genes].copy()

    sc.pp.highly_variable_genes(
        adata_benchmark, 
        n_top_genes=n_top, 
        flavor='seurat_v3'
    )
    final_hvg_indices = adata_benchmark.var.highly_variable
    adata_final = adata_benchmark[:, final_hvg_indices].to_memory()
    
    print(f"Selected the top {adata_final.n_vars} HVGs.")
    return adata_final


def select_hvg(adata, n_top = 3000):
    sc.pp.highly_variable_genes(adata, n_top_genes=n_top, flavor='seurat_v3', layer='raw_counts')
    ## flavor='seurat_v3' expects raw counts as it handles the NB-like dispersion internally
    return adata[:, adata.var.highly_variable].copy()


def load_and_normalize(filename, data_type):
    print(f"Loading {data_type}: {filename}")
    adata = ad.read_h5ad(raw_data / filename, backed='r')
    return normalize(adata.to_memory(), data_type=data_type)


def save_adata(adata, data_type, suffix):
    out_name = f"{data_type}_normalized_{suffix}.h5ad"
    adata.write(processed_data / out_name)
    print(f"Saved {data_type} {suffix} with {adata.n_vars} variables.\n")


def get_gene_mapping(data_type):
    json_path = data_dir / f"{data_type}_gene_to_transcripts.json"
    
    if json_path.exists():
        with open(json_path, 'r') as f:
            return json.load(f)
    
    print(f"{json_path} not found. Mapping from genes data...")
    genes_filename = f'{data_type}_processed_genes.h5ad'
    adata_genes = load_and_normalize(genes_filename, data_type)
    
    mapping = adata_genes.uns['gene_to_transcripts']

    with open(json_path, "w") as f:
        json.dump(mapping, f, cls=NumpyEncoder, indent=4)
    return mapping


def process_genes(filename, data_type, hvg_genes=None):
    adata = load_and_normalize(filename, data_type)

    if data_type == 'bulk':
        adata_subset = select_hvg(adata, n_top=3000)
    else:
        selected_genes = set(adata.var_names)
        common_genes = [g for g in hvg_genes if g in selected_genes]
        print(f"Overlap: {len(common_genes)}")
        adata_subset = adata[:, common_genes].copy()

    save_adata(adata_subset, data_type, "genes")
    return adata_subset.var_names.tolist()


def process_transcripts(filename, data_type, hvg_genes=None):
    mapping = get_gene_mapping(data_type)
    target_transcripts = [t for g in hvg_genes if g in mapping for t in mapping[g]]
    adata = load_and_normalize(filename, data_type)

    present_vars = set(adata.var_names)
    final_selection = [t for t in target_transcripts if t in present_vars]
    adata_subset = adata[:, final_selection].copy()

    save_adata(adata_subset, data_type, "transcripts")
    return final_selection

