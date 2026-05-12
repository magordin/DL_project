from pathlib import Path
from typing import List

import matplotlib.pyplot as plt
import pandas as pd


def collect_metrics(
    model_root: Path,
    representations: List[str],
    target_mode: str,
) -> pd.DataFrame:
    rows = []

    for representation in representations:
        repr_dir = model_root / representation / target_mode

        if not repr_dir.exists():
            print(f"WARNING: missing directory: {repr_dir}")
            continue

        for metrics_file in sorted(repr_dir.glob("metrics_seed*.csv")):
            seed = int(metrics_file.stem.replace("metrics_seed", ""))
            df = pd.read_csv(metrics_file)
            df["representation"] = representation
            df["seed"] = seed
            rows.append(df)

    if not rows:
        raise FileNotFoundError("No metrics files found.")

    return pd.concat(rows, ignore_index=True)


def plot_test_metric(
    metrics: pd.DataFrame,
    metric: str,
    out_path: Path,
) -> None:
    test_df = metrics[metrics["split"] == "test"].copy()

    grouped = (
        test_df
        .groupby("representation")[metric]
        .agg(["mean", "std"])
        .reset_index()
        .sort_values("mean", ascending=True)
    )

    plt.figure(figsize=(8, 5))
    plt.bar(
        grouped["representation"],
        grouped["mean"],
        yerr=grouped["std"],
        capsize=4,
    )
    plt.xlabel("Representation")
    plt.ylabel(metric)
    plt.title(f"Test {metric} by representation")
    plt.tight_layout()
    plt.savefig(out_path, dpi=200)
    plt.close()


def plot_training_histories(
    model_root: Path,
    representations: List[str],
    target_mode: str,
    plots_dir: Path,
) -> None:
    for representation in representations:
        repr_dir = model_root / representation / target_mode

        if not repr_dir.exists():
            continue

        history_files = sorted(repr_dir.glob("history_seed*.csv"))

        if not history_files:
            continue

        plt.figure(figsize=(8, 5))

        for history_file in history_files:
            seed = history_file.stem.replace("history_seed", "")
            history = pd.read_csv(history_file)

            plt.plot(
                history["epoch"],
                history["train_mse"],
                alpha=0.5,
                label=f"train seed {seed}",
            )
            plt.plot(
                history["epoch"],
                history["val_mse"],
                alpha=0.8,
                linestyle="--",
                label=f"val seed {seed}",
            )

        plt.xlabel("Epoch")
        plt.ylabel("MSE")
        plt.title(f"Training history: {representation}")
        plt.legend(fontsize=7)
        plt.tight_layout()
        plt.savefig(plots_dir / f"history_{representation}_{target_mode}.png", dpi=200)
        plt.close()