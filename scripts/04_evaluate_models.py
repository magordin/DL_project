#!/usr/bin/env python3

import argparse
from pathlib import Path

from src.data import (
    load_representation,
    load_target,
    load_split,
    apply_standardization,
)
from src.evaluation import (
    load_model_from_checkpoint,
    predict_array,
    evaluate_predictions,
    save_predictions_h5ad,
)
from src.utils import get_device


def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument("--input-repr", required=True, type=Path)
    parser.add_argument("--target-h5ad", required=True, type=Path)
    parser.add_argument("--target-mode", default="absolute")

    parser.add_argument("--model-path", required=True, type=Path)
    parser.add_argument("--split-json", required=True, type=Path)

    parser.add_argument("--out-metrics", required=True, type=Path)
    parser.add_argument("--out-preds", required=True, type=Path)

    parser.add_argument("--batch-size", type=int, default=256)
    parser.add_argument("--eval-split", default="test", choices=["train", "val", "test"])

    return parser.parse_args()


def main():
    args = parse_args()

    device = get_device()
    print(f"Device: {device}")

    x = load_representation(args.input_repr)
    split = load_split(args.split_json)

    model, checkpoint = load_model_from_checkpoint(args.model_path, device)
    model_type = checkpoint.get("model_type", "mse")

    y, target_adata = load_target(
        args.target_h5ad,
        model_type=model_type,
    )

    if args.eval_split not in split:
        raise ValueError(
            f"Split '{args.eval_split}' not found in {args.split_json}. "
            f"Available splits: {list(split.keys())}"
        )

    eval_idx = split[args.eval_split]

    print(f"Model type: {model_type}")
    print(f"Evaluation split: {args.eval_split}")
    print(f"Number of samples: {len(eval_idx)}")

    x_scaled = apply_standardization(
        x=x,
        mean=checkpoint["x_mean"],
        std=checkpoint["x_std"],
    )

    x_eval = x_scaled[eval_idx]
    y_eval = y[eval_idx]
    target_adata_eval = target_adata[eval_idx].copy()

    y_pred = predict_array(
        model=model,
        x=x_eval,
        batch_size=args.batch_size,
        device=device,
        model_type=model_type,
    )

    metrics = evaluate_predictions(
        y_true=y_eval,
        y_pred=y_pred,
        split={args.eval_split: list(range(len(eval_idx)))},
    )

    args.out_metrics.parent.mkdir(parents=True, exist_ok=True)
    metrics.to_csv(args.out_metrics, index=False)

    save_predictions_h5ad(
        path=args.out_preds,
        y_pred=y_pred,
        y_true=y_eval,
        target_adata=target_adata_eval,
        metadata={
            "target_mode": args.target_mode,
            "model_type": model_type,
            "model_path": str(args.model_path),
            "input_repr": str(args.input_repr),
            "eval_split": args.eval_split,
        },
    )

    print(metrics)
    print(f"Saved metrics: {args.out_metrics}")
    print(f"Saved predictions: {args.out_preds}")


if __name__ == "__main__":
    main()