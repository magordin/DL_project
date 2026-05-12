#!/usr/bin/env python3

import argparse
from pathlib import Path

from src.plotting import (
    collect_metrics,
    plot_test_metric,
    plot_training_histories,
)


def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument("--model-root", required=True, type=Path)
    parser.add_argument("--plots-dir", required=True, type=Path)
    parser.add_argument("--representations", nargs="+", required=True)
    parser.add_argument("--target-mode", default="absolute")

    return parser.parse_args()


def main():
    args = parse_args()
    args.plots_dir.mkdir(parents=True, exist_ok=True)

    metrics = collect_metrics(
        model_root=args.model_root,
        representations=args.representations,
        target_mode=args.target_mode,
    )

    summary_file = args.plots_dir / f"summary_metrics_{args.target_mode}.csv"
    metrics.to_csv(summary_file, index=False)

    for metric in [
        "mse",
        "mae",
        "rmse",
        "r2",
        "pearson_global",
        "spearman_global",
    ]:
        plot_test_metric(
            metrics=metrics,
            metric=metric,
            out_path=args.plots_dir / f"test_{metric}_{args.target_mode}.png",
        )

    plot_training_histories(
        model_root=args.model_root,
        representations=args.representations,
        target_mode=args.target_mode,
        plots_dir=args.plots_dir,
    )

    print(f"Saved summary: {summary_file}")
    print(f"Saved plots in: {args.plots_dir}")


if __name__ == "__main__":
    main()