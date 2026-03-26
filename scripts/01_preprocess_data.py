import sys
from pathlib import Path
import anndata as ad
import json

sys.path.append(str(Path(__file__).resolve().parent.parent))
from src.config import *
from src.preprocessing import normalize, select_hvg, NumpyEncoder


def load_and_normalize(filename, data_type):
    print(f"Loading {data_type}: {filename}")
    adata = ad.read_h5ad(raw_data / filename, backed='r')
    return normalize(adata.to_memory(), data_type=data_type)


def save_adata(adata, data_type, suffix):
    out_name = f"{data_type}_normalized_{suffix}.h5ad"
    adata.write(processed_data / out_name)
    print(f"Saved {data_type} {suffix} with {adata.n_vars} variables.\n")


def get_gene_mapping(data_type):
    json_path = data_dir / f"{data_type}_gene_to_transcripts.json"
    
    if json_path.exists():
        with open(json_path, 'r') as f:
            return json.load(f)
    
    print(f"{json_path} not found. Mapping from genes data...")
    genes_filename = f'{data_type}_processed_genes.h5ad'
    adata_genes = load_and_normalize(genes_filename, data_type)
    
    mapping = adata_genes.uns['gene_to_transcripts']

    with open(json_path, "w") as f:
        json.dump(mapping, f, cls=NumpyEncoder, indent=4)
    
    return mapping


def process_genes(filename, data_type, hvg_genes=None):
    adata = load_and_normalize(filename, data_type)

    if data_type == 'bulk':
        adata_subset = select_hvg(adata, n_top=5000)
    else:
        selected_genes = set(adata.var_names)
        common_genes = [g for g in hvg_genes if g in selected_genes]
        print(f"Overlap: {len(common_genes)}")
        adata_subset = adata[:, common_genes].copy()

    save_adata(adata_subset, data_type, "genes")
    return adata_subset.var_names.tolist()


def process_transcripts(filename, data_type, hvg_genes=None):
    mapping = get_gene_mapping(data_type)
    target_transcripts = [t for g in hvg_genes if g in mapping for t in mapping[g]]
    adata = load_and_normalize(filename, data_type)

    present_vars = set(adata.var_names)
    final_selection = [t for t in target_transcripts if t in present_vars]
    adata_subset = adata[:, final_selection].copy()

    save_adata(adata_subset, data_type, "transcripts")
    return final_selection
    

def main():
    processed_data.mkdir(parents=True, exist_ok=True)
    print(f"Saving to {processed_data}...")

    hvg_list = process_genes("bulk_processed_genes.h5ad", "bulk")
    # _ = process_genes("sc_processed_genes.h5ad", "sc", hvg_genes=hvg_list)

    _ = process_transcripts("bulk_processed_transcripts.h5ad", "bulk", hvg_list)
    # _ = process_transcripts("sc_processed_transcripts.h5ad", "sc", hvg_list)

    print('Preprocessing complete')


if __name__ == "__main__":
    main()