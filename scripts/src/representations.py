import numpy as np
import anndata as ad
from sklearn.decomposition import PCA, TruncatedSVD
from sklearn.preprocessing import StandardScaler
import logging 

from scipy.sparse import issparse
from scripts.src.vae import train_vae_representation

def generate_pca_representation(adata, n_comps=128):
    sparse_input = issparse(adata.X)
    
    if sparse_input:
        logging.info("Sparse matrix, performing TruncatedSVD.")
        scaler = StandardScaler(with_mean=False)
        model = TruncatedSVD(n_components=n_comps)
    else:
        logging.info("Dense matrix, performing PCA.")
        scaler = StandardScaler(with_mean=True)
        model = PCA(n_components=n_comps)

    scaled_data = scaler.fit_transform(adata.X)

    if isinstance(model, PCA) and issparse(scaled_data):
        scaled_data = scaled_data.toarray()

    pca_coords = model.fit_transform(scaled_data)

    return pca_coords.astype(np.float32)

    
def generate_vae_representation(adata, latent_dim=128, seed=1):
    return train_vae_representation(
        adata=adata,
        latent_dim=latent_dim,
        seed=seed,
        epochs=100,
        batch_size=128,
        hidden_dim=512,
        lr=1e-3,
        beta=1.0,
    )