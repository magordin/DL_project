import sys
from pathlib import Path
import anndata as ad

sys.path.append(str(Path(__file__).resolve().parent.parent))
from src.config import raw_data, processed_data
from src.preprocessing import normalize, select_hvg

def main():
    adata = ad.read_h5ad(raw_data / "bulk_processed_genes.h5ad", backed='r')
    print('Bulk data loaded')

    processed = normalize(adata.to_memory())
    processed_hvg = select_hvg(processed, 5000)

    processed_data.mkdir(parents=True, exist_ok=True) 
    print(f"Saving to {processed_data}...")

    processed_hvg.write(processed_data / "bulk_hvg_normalized.h5ad")
    print('Preprocessing complete')

if __name__ == "__main__":
    main()