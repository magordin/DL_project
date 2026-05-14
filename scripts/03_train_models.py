#!/usr/bin/env python3

import argparse
from pathlib import Path

import numpy as np

from src.data import (
    load_representation,
    load_target,
    make_split,
    save_split,
    standardize_train_only,
)
from src.training import train_mlp, save_checkpoint
from src.utils import set_seed, get_device


def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument("--input-repr", required=True, type=Path)
    parser.add_argument("--target-h5ad", required=True, type=Path)
    parser.add_argument("--target-mode", default="absolute")

    parser.add_argument(
        "--model-type",
        type=str,
        default="mse",
        choices=["mse", "gaussian", "nb"],
    )

    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument("--split-seed", type=int, default=42)

    parser.add_argument("--epochs", type=int, default=100)
    parser.add_argument("--batch-size", type=int, default=128)

    parser.add_argument("--learning-rate", type=float, default=1e-3)
    parser.add_argument("--hidden-dim", type=int, default=512)
    parser.add_argument("--dropout", type=float, default=0.1)
    parser.add_argument("--weight-decay", type=float, default=1e-5)
    parser.add_argument("--n-layers", type=int, default=2)

    parser.add_argument("--out-model", required=True, type=Path)
    parser.add_argument("--out-history", required=True, type=Path)
    parser.add_argument("--out-split", required=True, type=Path)

    return parser.parse_args()


def main():
    args = parse_args()
    set_seed(args.seed)

    x = load_representation(args.input_repr)
    y, _ = load_target(
        args.target_h5ad,
        model_type=args.model_type,
    )

    if x.shape[0] != y.shape[0]:
        raise ValueError(f"Sample mismatch: X={x.shape[0]}, Y={y.shape[0]}")

    if args.model_type == "nb" and np.any(y < 0):
        raise ValueError("Negative Binomial model requires non-negative count targets.")

    split = make_split(n_samples=x.shape[0], seed=args.split_seed)
    train_idx = np.asarray(split["train"])

    x_scaled, x_mean, x_std = standardize_train_only(x, train_idx)

    device = get_device()

    print(f"Device: {device}")
    print(f"Model type: {args.model_type}")
    print(f"Input shape: {x_scaled.shape}")
    print(f"Target shape: {y.shape}")
    print(f"Training seed: {args.seed}")
    print(f"Split seed: {args.split_seed}")

    model, history, best_val = train_mlp(
        x=x_scaled,
        y=y,
        split=split,
        hidden_dim=args.hidden_dim,
        dropout=args.dropout,
        n_layers=args.n_layers,
        learning_rate=args.learning_rate,
        weight_decay=args.weight_decay,
        batch_size=args.batch_size,
        epochs=args.epochs,
        device=device,
        model_type=args.model_type,
        seed=args.seed,
    )

    config = {
        "input_repr": str(args.input_repr),
        "target_h5ad": str(args.target_h5ad),
        "target_mode": args.target_mode,
        "model_type": args.model_type,
        "seed": args.seed,
        "split_seed": args.split_seed,
        "epochs": args.epochs,
        "batch_size": args.batch_size,
        "learning_rate": args.learning_rate,
        "hidden_dim": args.hidden_dim,
        "dropout": args.dropout,
        "weight_decay": args.weight_decay,
        "n_layers": args.n_layers,
        "input_dim": x.shape[1],
        "output_dim": y.shape[1],
        "best_val_loss": best_val,
        "train_n": len(split["train"]),
        "val_n": len(split["val"]),
        "test_n": len(split["test"]),
    }

    args.out_model.parent.mkdir(parents=True, exist_ok=True)
    args.out_history.parent.mkdir(parents=True, exist_ok=True)
    args.out_split.parent.mkdir(parents=True, exist_ok=True)

    save_checkpoint(args.out_model, model, x_mean, x_std, config)
    history.to_csv(args.out_history, index=False)
    save_split(split, args.out_split)

    print(f"Best validation loss: {best_val:.6f}")
    print(f"Saved model: {args.out_model}")
    print(f"Saved history: {args.out_history}")
    print(f"Saved split: {args.out_split}")


if __name__ == "__main__":
    main()