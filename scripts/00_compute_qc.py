#!/usr/bin/env python3

import argparse
from pathlib import Path

import anndata as ad

from src.qc import (
    validate_obs_alignment,
    load_mapping_json,
    compute_gene_qc_table_backed,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Compute gene-level QC summary metrics from gene and transcript AnnData files."
    )
    parser.add_argument(
        "--gene-h5ad",
        type=Path,
        required=True,
        help="Path to gene-level .h5ad file",
    )
    parser.add_argument(
        "--tx-h5ad",
        type=Path,
        required=True,
        help="Path to transcript-level .h5ad file",
    )
    parser.add_argument(
        "--mapping-json",
        type=Path,
        required=True,
        help="Path to gene->transcripts mapping JSON",
    )
    parser.add_argument(
        "--out-csv",
        type=Path,
        required=True,
        help="Output CSV path for gene-level QC summary",
    )
    parser.add_argument(
        "--eps",
        type=float,
        default=1e-12,
        help="Small constant used to avoid division/log issues",
    )
    parser.add_argument(
        "--min-isoforms",
        type=int,
        default=2,
        help="Minimum number of mapped transcripts required to compute QC for a gene",
    )
    parser.add_argument(
    "--row-block-size",
    type=int,
    default=256,
    help="Number of samples to process per block",
    )
    parser.add_argument(
        "--max-genes",
        type=int,
        default=None,
        help="Optional cap for debugging/benchmarking",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    if not args.gene_h5ad.exists():
        raise FileNotFoundError(f"Gene h5ad not found: {args.gene_h5ad}")
    if not args.tx_h5ad.exists():
        raise FileNotFoundError(f"Transcript h5ad not found: {args.tx_h5ad}")
    if not args.mapping_json.exists():
        raise FileNotFoundError(f"Mapping JSON not found: {args.mapping_json}")

    args.out_csv.parent.mkdir(parents=True, exist_ok=True)

    print(f"Loading gene AnnData: {args.gene_h5ad}")
    adata_gene = ad.read_h5ad(args.gene_h5ad, backed="r")

    print(f"Loading transcript AnnData: {args.tx_h5ad}")
    adata_tx = ad.read_h5ad(args.tx_h5ad, backed="r")

    print(f"Loading mapping JSON: {args.mapping_json}")
    gene_to_tx = load_mapping_json(args.mapping_json)

    print("Validating sample alignment between gene and transcript matrices...")
    validate_obs_alignment(adata_gene, adata_tx)

    print("Computing gene-level QC summary table...")
    qc_df = compute_gene_qc_table_backed(
        adata_gene=adata_gene,
        adata_tx=adata_tx,
        gene_to_tx=gene_to_tx,
        min_isoforms=args.min_isoforms,
        eps=args.eps,
    )

    if qc_df is None:
        raise RuntimeError("compute_gene_qc_table returned None")

    sort_cols = [c for c in ["n_isoforms", "corr_gene_vs_sum_iso"] if c in qc_df.columns]
    if sort_cols:
        qc_df = qc_df.sort_values(
            by=sort_cols,
            ascending=[False] * len(sort_cols),
        ).reset_index(drop=True)

    qc_df.to_csv(args.out_csv, index=False)

    print(f"Saved QC table to: {args.out_csv}")
    print(f"Rows written: {len(qc_df)}")
    print("Done.")


if __name__ == "__main__":
    main()