import scanpy as sc
import logging
import json
import numpy as np
from pathlib import Path
import anndata as ad
import pandas as pd
from scipy.sparse import issparse


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


def save_adata(adata, out_dir, data_type, suffix):
    out_path = Path(out_dir)
    out_name = f"{data_type}_normalized_{suffix}.h5ad"
    adata.write(out_path / out_name)
    print(f"Saved {data_type} {suffix} with {adata.n_vars} variables.\n")


def get_gene_mapping(mapping_json_path, adata):
    json_path = Path(mapping_json_path)
    if json_path.exists():
        with open(json_path, 'r') as f:
            return json.load(f)
    if adata is not None:
        logging.warning(f"Mapping {json_path} not found. Generating from adata...")
        
        if 'gene_id' not in adata.var.columns:
            raise KeyError("adata.var must have a 'gene_id' column to generate mapping.")

        mapping = {}
        for tx_id, gene_id in zip(adata.var_names, adata.var['gene_id']):
            if gene_id not in mapping:
                mapping[gene_id] = []
            mapping[gene_id].append(tx_id)

        json_path.parent.mkdir(parents=True, exist_ok=True)
        with open(json_path, 'w') as f:
            json.dump(mapping, f, indent=4, cls=NumpyEncoder)
            
        logging.info(f"Successfully generated and saved mapping for {len(mapping)} genes.")
        return mapping
    
    else:
        raise FileNotFoundError(
            f"Mapping file {json_path} not found; no adata provided to generate it.")
    

def get_isoform_proportions(adata, mapping):
    X = adata.X.toarray() if issparse(adata.X) else adata.X.copy()
    var_names = adata.var_names.values
    var_idx = {name: i for i, name in enumerate(var_names)}

    for gene, isoforms in mapping.items():
        indices = [var_idx[iso] for iso in isoforms if iso in var_idx]
        
        if not indices:
            continue
        gene_sum = X[:, indices].sum(axis=1, keepdims=True)

        X[:, indices] = X[:, indices] / (gene_sum + 1e-9)

    adata.layers["isoform_proportions"] = X.astype(np.float32)
    return adata