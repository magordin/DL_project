from typing import Dict, Tuple

import numpy as np
import pandas as pd
import torch
from torch import nn
from torch.utils.data import DataLoader, TensorDataset

from src.models import (
    MLPRegressor,
    MLPGaussianRegressor,
    MLPNegativeBinomial,
    gaussian_nll_loss,
    negative_binomial_loss,
)


def make_loaders(
    x: np.ndarray,
    y: np.ndarray,
    split: Dict[str, list],
    batch_size: int,
) -> Dict[str, DataLoader]:
    loaders = {}

    for split_name, idx in split.items():
        idx = np.asarray(idx)

        dataset = TensorDataset(
            torch.from_numpy(x[idx]).float(),
            torch.from_numpy(y[idx]).float(),
        )

        loaders[split_name] = DataLoader(
            dataset,
            batch_size=batch_size,
            shuffle=(split_name == "train"),
        )

    return loaders


def run_epoch(
    model: nn.Module,
    loader: DataLoader,
    optimizer,
    loss_fn,
    device: torch.device,
    train: bool,
    model_type: str,
) -> float:
    model.train(train)

    total_loss = 0.0
    total_n = 0

    for xb, yb in loader:
        xb = xb.to(device)
        yb = yb.to(device)

        if train:
            optimizer.zero_grad(set_to_none=True)

        output = model(xb)

        if model_type == "mse":
            loss = loss_fn(output, yb)

        elif model_type == "gaussian":
            mean, std = output
            loss = gaussian_nll_loss(yb, mean, std)

        elif model_type == "nb":
            mu, theta = output
            loss = negative_binomial_loss(yb, mu, theta)

        else:
            raise ValueError(f"Unknown model_type: {model_type}")

        if train:
            loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), max_norm=5.0)
            optimizer.step()

        total_loss += loss.item() * xb.size(0)
        total_n += xb.size(0)

    return total_loss / max(total_n, 1)


def train_mlp(
    x: np.ndarray,
    y: np.ndarray,
    split: Dict[str, list],
    hidden_dim: int,
    dropout: float,
    n_layers: int,
    learning_rate: float,
    weight_decay: float,
    batch_size: int,
    epochs: int,
    device: torch.device,
    model_type: str = "mse",
) -> Tuple[nn.Module, pd.DataFrame, float]:
    loaders = make_loaders(x, y, split, batch_size)

    if model_type == "mse":
        model = MLPRegressor(
            input_dim=x.shape[1],
            output_dim=y.shape[1],
            hidden_dim=hidden_dim,
            dropout=dropout,
            n_layers=n_layers,
        ).to(device)

        loss_fn = nn.MSELoss()

    elif model_type == "gaussian":
        model = MLPGaussianRegressor(
            input_dim=x.shape[1],
            output_dim=y.shape[1],
            hidden_dim=hidden_dim,
            dropout=dropout,
            n_layers=n_layers,
        ).to(device)

        loss_fn = gaussian_nll_loss

    elif model_type == "nb":
        model = MLPNegativeBinomial(
            input_dim=x.shape[1],
            output_dim=y.shape[1],
            hidden_dim=hidden_dim,
            dropout=dropout,
            n_layers=n_layers,
        ).to(device)

        loss_fn = negative_binomial_loss

    else:
        raise ValueError(f"Unknown model_type: {model_type}")

    optimizer = torch.optim.AdamW(
        model.parameters(),
        lr=learning_rate,
        weight_decay=weight_decay,
    )

    best_val = float("inf")
    best_state = None
    history = []

    for epoch in range(1, epochs + 1):
        train_loss = run_epoch(
            model=model,
            loader=loaders["train"],
            optimizer=optimizer,
            loss_fn=loss_fn,
            device=device,
            train=True,
            model_type=model_type,
        )

        val_loss = run_epoch(
            model=model,
            loader=loaders["val"],
            optimizer=optimizer,
            loss_fn=loss_fn,
            device=device,
            train=False,
            model_type=model_type,
        )

        history.append({
            "epoch": epoch,
            "model_type": model_type,
            "train_loss": train_loss,
            "val_loss": val_loss,
        })

        if val_loss < best_val:
            best_val = val_loss
            best_state = {
                key: value.detach().cpu().clone()
                for key, value in model.state_dict().items()
            }

        if epoch == 1 or epoch % 10 == 0:
            print(
                f"epoch={epoch:04d} "
                f"model_type={model_type} "
                f"train_loss={train_loss:.6f} "
                f"val_loss={val_loss:.6f}"
            )

    if best_state is None:
        raise RuntimeError("No best model state was saved during training.")

    model.load_state_dict(best_state)

    return model, pd.DataFrame(history), best_val


def save_checkpoint(
    path,
    model: nn.Module,
    x_mean,
    x_std,
    config: dict,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)

    torch.save({
        "model_state": model.state_dict(),
        "model_type": config["model_type"],
        "input_dim": config["input_dim"],
        "output_dim": config["output_dim"],
        "hidden_dim": config["hidden_dim"],
        "dropout": config["dropout"],
        "n_layers": config["n_layers"],
        "x_mean": x_mean,
        "x_std": x_std,
        "config": config,
    }, path)