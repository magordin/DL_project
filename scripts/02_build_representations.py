import argparse
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

from scripts.src.representations import (
    generate_pca_representation,
    generate_vae_representation,
)

def parse_args():
    parser = argparse.ArgumentParser(description="Build representations from filtered data.")
    parser.add_argument("--input-h5ad", type=str, required=True, help="Path to input X h5ad")
    parser.add_argument("--repr-type", type=str, required=True, choices=["raw", "pca", "vae", "geneformer"])
    parser.add_argument("--latent-dim", type=int, default=128)
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument("--out-file", type=str, required=True, help="Specific path for .npy output")
    parser.add_argument("--beta", type=float, default=1.0)
    parser.add_argument("--history-file", type=str, default=None)
    return parser.parse_args()

def main():
    args = parse_args()
    np.random.seed(args.seed)
    adata = ad.read_h5ad(args.input_h5ad)

    print(f"Generating {args.repr_type} representation...")
  
    if args.repr_type == "raw":
        repr_matrix = adata.X.toarray() if hasattr(adata.X, "toarray") else adata.X
        
    elif args.repr_type == "pca":
        repr_matrix = generate_pca_representation(adata, n_comps=args.latent_dim)
        
    elif args.repr_type == "vae":
        repr_matrix, history = generate_vae_representation(
            adata,
            latent_dim=args.latent_dim,
            seed=args.seed,
            beta=args.beta,
        )
        
    else:
        raise ValueError(f"Representation {args.repr_type} not supported yet.")

    out_path = Path(args.out_file)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    
    np.save(out_path, repr_matrix)
    print(f"Successfully saved {args.repr_type} matrix of shape {repr_matrix.shape}")

    if args.repr_type == "vae" and args.history_file is not None:

        history_path = Path(args.history_file)
        history_path.parent.mkdir(parents=True, exist_ok=True)

        pd.DataFrame(history).to_csv(history_path, index=False)
        print(f"Saved VAE training history to {history_path}")

if __name__ == "__main__":
    main()