from pathlib import Path
from typing import Dict

import anndata as ad
import numpy as np
import pandas as pd
import torch
from scipy.stats import pearsonr, spearmanr
from torch.utils.data import DataLoader, TensorDataset

from src.models import MLPRegressor


def load_model_from_checkpoint(path: Path, device: torch.device) -> tuple:
    checkpoint = torch.load(path, map_location="cpu")

    model = MLPRegressor(
        input_dim=checkpoint["input_dim"],
        output_dim=checkpoint["output_dim"],
        hidden_dim=checkpoint["hidden_dim"],
        dropout=checkpoint["dropout"],
        n_layers=checkpoint["n_layers"],
    ).to(device)

    model.load_state_dict(checkpoint["model_state"])
    model.eval()

    return model, checkpoint


def predict_array(
    model: torch.nn.Module,
    x: np.ndarray,
    batch_size: int,
    device: torch.device,
) -> np.ndarray:
    dataset = TensorDataset(torch.from_numpy(x))
    loader = DataLoader(dataset, batch_size=batch_size, shuffle=False)

    preds = []

    with torch.no_grad():
        for (xb,) in loader:
            xb = xb.to(device)
            pred = model(xb).cpu().numpy()
            preds.append(pred)

    return np.vstack(preds).astype(np.float32)


def safe_corr(a: np.ndarray, b: np.ndarray, method: str) -> float:
    if np.std(a) == 0 or np.std(b) == 0:
        return np.nan

    if method == "pearson":
        return float(pearsonr(a, b)[0])

    if method == "spearman":
        return float(spearmanr(a, b).correlation)

    raise ValueError(f"Unknown correlation method: {method}")


def compute_metrics_for_split(
    y_true: np.ndarray,
    y_pred: np.ndarray,
    split_name: str,
) -> dict:
    residual = y_pred - y_true

    mse = float(np.mean(residual ** 2))
    mae = float(np.mean(np.abs(residual)))
    rmse = float(np.sqrt(mse))

    ss_res = float(np.sum(residual ** 2))
    ss_tot = float(np.sum((y_true - y_true.mean()) ** 2))
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else np.nan

    yt = y_true.ravel()
    yp = y_pred.ravel()

    return {
        "split": split_name,
        "mse": mse,
        "mae": mae,
        "rmse": rmse,
        "r2": r2,
        "pearson_global": safe_corr(yt, yp, "pearson"),
        "spearman_global": safe_corr(yt, yp, "spearman"),
    }


def evaluate_predictions(
    y_true: np.ndarray,
    y_pred: np.ndarray,
    split: Dict[str, list],
) -> pd.DataFrame:
    rows = []

    for split_name, idx in split.items():
        idx = np.asarray(idx)
        rows.append(
            compute_metrics_for_split(
                y_true=y_true[idx],
                y_pred=y_pred[idx],
                split_name=split_name,
            )
        )

    return pd.DataFrame(rows)


def save_predictions_h5ad(
    path: Path,
    y_pred: np.ndarray,
    y_true: np.ndarray,
    target_adata: ad.AnnData,
    metadata: dict,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)

    pred_adata = ad.AnnData(
        X=y_pred,
        obs=target_adata.obs.copy(),
        var=target_adata.var.copy(),
    )

    pred_adata.layers["true"] = y_true

    for key, value in metadata.items():
        pred_adata.uns[key] = value

    pred_adata.write_h5ad(path)