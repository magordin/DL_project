import json
from pathlib import Path
from typing import Dict, Tuple
import numpy as np
import anndata as ad
from scipy.sparse import issparse

import anndata as ad
import numpy as np
import scipy.sparse as sp


def to_dense_float32(x) -> np.ndarray:
    if sp.issparse(x):
        x = x.toarray()
    return np.asarray(x, dtype=np.float32)


def load_representation(path: Path) -> np.ndarray:
    x = np.load(path)
    if x.ndim != 2:
        raise ValueError(f"Representation must be 2D. Got shape: {x.shape}")
    return x.astype(np.float32)


def load_target(
    h5ad_path,
    model_type: str = "mse",
):
    adata = ad.read_h5ad(h5ad_path)

    if model_type in ["mse", "gaussian"]:
        y = adata.X

    elif model_type == "nb":
        if "raw_counts" not in adata.layers:
            raise ValueError("Expected adata.layers['raw_counts'] for NB targets.")

        y = adata.layers["raw_counts"]

    else:
        raise ValueError(f"Unknown model_type: {model_type}")

    if issparse(y):
        y = y.toarray()

    y = np.asarray(y, dtype=np.float32)

    return y, adata


def make_split(
    n_samples: int,
    seed: int,
    train_frac: float = 0.70,
    val_frac: float = 0.15,
) -> Dict[str, list]:
    rng = np.random.default_rng(seed)
    idx = rng.permutation(n_samples)

    n_train = int(n_samples * train_frac)
    n_val = int(n_samples * val_frac)

    return {
        "train": idx[:n_train].tolist(),
        "val": idx[n_train:n_train + n_val].tolist(),
        "test": idx[n_train + n_val:].tolist(),
    }


def save_split(split: Dict[str, list], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as f:
        json.dump(split, f, indent=2)


def load_split(path: Path) -> Dict[str, list]:
    with open(path) as f:
        return json.load(f)


def standardize_train_only(
    x: np.ndarray,
    train_idx: np.ndarray,
    eps: float = 1e-8,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    mean = x[train_idx].mean(axis=0, keepdims=True).astype(np.float32)
    std = x[train_idx].std(axis=0, keepdims=True).astype(np.float32)
    std = np.where(std < eps, 1.0, std).astype(np.float32)
    x_scaled = ((x - mean) / std).astype(np.float32)
    return x_scaled, mean, std


def apply_standardization(
    x: np.ndarray,
    mean: np.ndarray,
    std: np.ndarray,
) -> np.ndarray:
    return ((x - mean) / std).astype(np.float32)