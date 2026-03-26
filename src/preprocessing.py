import scanpy as sc
import json
import numpy as np


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


def normalize(adata, data_type = None, target_sum = 1e6):
    adata.layers["raw_counts"] = adata.X.copy()
    sc.pp.filter_genes(adata, min_cells=10)
    sc.pp.filter_cells(adata, min_genes=10)

    if data_type == 'bulk':
        print('Model raw counts as a Negative Binomial')
    else:
        print('Model raw counts as a Poisson')

    return adata


def select_hvg(adata, n_top = 5000):
    sc.pp.highly_variable_genes(adata, n_top_genes=n_top, flavor='seurat_v3', layer='raw_counts') 
    ## flavor='seurat_v3' expects raw counts as it handles the NB-like dispersion internally
    
    return adata[:, adata.var.highly_variable].copy()
