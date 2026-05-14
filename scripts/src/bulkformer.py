import torch
import sys
import os
import scanpy as sc
import numpy as np
from tqdm import tqdm
import anndata as ad
import rich

from BulkFormer import BulkFormer
from BulkFormer_block import BulkFormer_block

def extract_bulkformer_features(adata, latent_dim=128, checkpoint_path=None, DEVICE=None):
    if DEVICE is None:
        DEVICE = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print(f"Running on: {DEVICE}")

    if checkpoint_path is None:
        checkpoint_path = '/work3/s252608/DL_project/BulkFormer/model/Bulkformer_ckpt_epoch_29.pt'

    X = torch.tensor(adata.X.todense() if hasattr(adata.X, "todense") else adata.X).float()

    N_samples, N_genes = X.shape
    print(f"Data shape: {N_samples} samples, {N_genes} genes.")

    graph = torch.zeros((2,0), dtype=torch.long).to(DEVICE)
    gene_emb = torch.nn.Parameter(torch.randn(N_genes, latent_dim).to(DEVICE))
    
    model = BulkFormer(
        dim=latent_dim, 
        graph=graph,
        gene_emb=gene_emb,
        gene_length=N_genes,
        gb_repeat=1
    ).to(DEVICE)

    if not os.path.exists(checkpoint_path):
        raise FileNotFoundError(f"Missing checkpoint at {checkpoint_path}")
        
    checkpoint = torch.load(checkpoint_path, map_location=DEVICE)
    state_dict = {k: v for k, v in checkpoint.items() if "gene_emb" not in k}
    model.load_state_dict(state_dict, strict=False)
    
    model.eval()
    batch_size = 64
    all_embeddings = []
    
    with torch.no_grad():
        for i in tqdm(range(0, N_samples, batch_size)):
            batch_x = X[i : i + batch_size].to(DEVICE)
            out = model(batch_x)
            pooled_out = out.max(dim=1).values
            all_embeddings.append(pooled_out.cpu().numpy())

    return np.concatenate(all_embeddings, axis=0)