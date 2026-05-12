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
    bulkformer_root = '/work3/s252608/DL_project/BulkFormer'

    DEVICE = torch.device("cuda") if torch.cuda.is_available() else torch.device("cpu")
    rich.print(f"Device: [red]{DEVICE}")

    adata_path = '/work3/s252608/DL_project/data/processed/bulk_normalized_x_input.h5ad'
    adata = ad.read_h5ad(adata_path)
    X = torch.tensor(adata.X.todense() if hasattr(adata.X, "todense") else adata.X).float()
    
    N_samples, N_genes = X.shape
    D = 128

    dummy_graph = torch.zeros((2,0), dtype = torch.long).to(DEVICE)
    dummy_gene_emb = torch.randn(N_genes, D).to(DEVICE)
    model = BulkFormer(
        dim = D, 
        graph = dummy_graph,
        gene_emb = dummy_gene_emb,
        gene_length = N_genes,
        gb_repeat = 1
    ).to(DEVICE)

    ckpt_path = os.path.join(bulkformer_root, 'model_Bulkformer_ckpt_epoch_29.pt')
    if os.path.exists(ckpt_path):
        checkpoint = torch.load(ckpt_path, map_location = DEVICE)
        model.load_state_dict(checkpoint, strict = False)
        print('Weights loaded')
    model.eval()

    batch_size = 64
    all_embeddings = []
    
    with torch.no_grad():
        for i in tqdm(range(0, N_samples, batch_size)):
            batch_x = X[i : i + batch_size].to(DEVICE)
            out = model(batch_x)
            pooled_out = out.max(dim = 1)[0]
            all_embeddings.append(pooled_out.cpu().numpy())

    final_embeddings = np.concatenate(all_embeddings, axis = 0)
    save_path = "/work3/s252608/DL_project/data/bulkformer_embeddings_max.npy"
    np.save(save_path, final_embeddings)

if __name__ == "__main__":
    run_inference()