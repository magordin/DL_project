import scanpy as sc

def normalize(adata, target_sum = 1e6):
    if 'log1p' in adata.uns:
        del adata.uns['log1p']
    sc.pp.filter_genes(adata, min_cells=10)
    sc.pp.filter_cells(adata, min_genes=10)
    sc.pp.log1p(adata, target_sum)
    return adata

def select_hvg(adata, n_top = 5000):
    sc.pp.highly_variable_genes(adata, n_top_genes=n_top)
    return adata[:, adata.var.highly_variable].to_memory()