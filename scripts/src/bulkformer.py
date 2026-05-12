import torch
import sys
import os
import scanpy as sc
import numpy as np
from tqdm import tqdm
import anndata as ad

bulkformer_path = '/work3/s252608/DL_project/BulkFormer'
sys.path.appen(bulkformer_path)
from utils.BulkFormer import BulkFormer

def run_inference():
    device = 'cuda' if torch.cuda.is_available() else 'cpu'
    print(f'Using device {device}')

    adata_path = '/work3/s252608/DL_project/data/processed/bulk_normalized_x_input.h5ad'
    adata = ad.read_h5ad(adata_path)
    X = torch.tensor(adata.X.todense() if hasattr(adata.X, "todense") else adata.X).float()
    
    N_samples, N_genes = X.shape
    D = 128

    dummy_graph = torch.zeros((2,0), dtype = torch.long).to(device)
    dummy_gene_emb = torch.randn(N_genes, D).to(device)
    model = BulkFormer(
        dim = D, 
        graph = dummy_graph,
        gene_emb = dummy_gene_emb,
        gene_length = N_genes,
        gb_repeat = 1
    ).to(device)