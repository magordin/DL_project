import numpy as np
import anndata as ad
from sklearn.decomposition import PCA


def generate_pca_representation(adata, data_type="bulk", n_comps=50, out_dir = None):
    x_data = adata.X.toarray() if hasattr(adata.X, "toarray") else adata.X
    pca = PCA(n_components=n_comps)
    pca_coords = pca.fit_transform(x_data)
    
    return pca_coords