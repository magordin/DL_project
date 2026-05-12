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

def run_inference():
    DEVICE = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print(f"Running on: {DEVICE}")

    adata_path = '/work3/s252608/DL_project/data/processed/bulk_normalized_x_input.h5ad'
    adata = ad.read_h5ad(adata_path)
    
    X = torch.tensor(adata.X.todense() if hasattr(adata.X, "todense") else adata.X).float()
    N_samples, N_genes = X.shape
    print(f"Data shape: {N_samples} samples, {N_genes} genes.")

    D = 128
    graph = torch.zeros((2,0), dtype=torch.long).to(DEVICE)
    gene_emb = torch.nn.Parameter(torch.randn(N_genes, D).to(DEVICE))
    
    model = BulkFormer(
        dim=D, 
        graph=graph,
        gene_emb=gene_emb,
        gene_length=N_genes,
        gb_repeat=1
    ).to(DEVICE)

    ckpt_path = '/work3/s252608/DL_project/BulkFormer/model/Bulkformer_ckpt_epoch_29.pt'
    checkpoint = torch.load(ckpt_path, map_location=DEVICE)

    state_dict = {k: v for k, v in checkpoint.items() if "gene_emb" not in k}
    
    msg = model.load_state_dict(state_dict, strict=False)
    print(f"Transformer weights loaded. Missing keys: {msg.missing_keys}")
    
    model.eval()

    batch_size = 64
    all_embeddings = []
    
    with torch.no_grad():
        for i in tqdm(range(0, N_samples, batch_size)):
            batch_x = X[i : i + batch_size].to(DEVICE)
            out = model(batch_x)
            pooled_out = out.mean(dim=1) 
            all_embeddings.append(pooled_out.cpu().numpy())

    final_embeddings = np.concatenate(all_embeddings, axis=0)
    save_path = "/work3/s252608/DL_project/data/bulkformer_embeddings.npy"
    np.save(save_path, final_embeddings)
    print(f"Final embeddings saved: {final_embeddings.shape}")

if __name__ == "__main__":
    run_inference()