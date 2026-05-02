import scanpy as sc
import json
import numpy as np
from pathlib import Path
import anndata as ad
import pandas as pd


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


def normalize(adata, target_sum=1e4):
    adata.layers["raw_counts"] = adata.X.copy() 
    sc.pp.normalize_total(adata, target_sum=target_sum)
    sc.pp.log1p(adata)
    return adata


def select_hvgs(adata, n_top=3000):
    sc.pp.highly_variable_genes(
        adata, 
        n_top_genes=n_top, 
        flavor='seurat_v3'
    )
    final_hvg_indices = adata.var.highly_variable
    adata_final = adata[:, final_hvg_indices].copy()
    
    print(f"Selected the top {adata_final.n_vars} HVGs.")
    return adata_final


def save_adata(adata, out_dir, data_type, suffix):
    out_path = Path(out_dir)
    out_name = f"{data_type}_normalized_{suffix}.h5ad"
    adata.write(out_path / out_name)
    print(f"Saved {data_type} {suffix} with {adata.n_vars} variables.\n")


def get_gene_mapping(mapping_json_path):
    json_path = Path(mapping_json_path)
    if json_path.exists():
        with open(json_path, 'r') as f:
            return json.load(f)
    else:
        raise FileNotFoundError(f"Mapping file not found at {json_path}")