#!/usr/bin/env python3
import argparse
from pathlib import Path

import anndata as ad

from src.qc import (
    validate_obs_alignment,
    load_mapping_json,
    compute_gene_qc_table_backed,
    QCThresholds,
)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Compute gene-level QC summary metrics.")
    p.add_argument("--gene-h5ad",    type=Path, required=True)
    p.add_argument("--tx-h5ad",      type=Path, required=True)
    p.add_argument("--mapping-json", type=Path, required=True)
    p.add_argument("--out-csv",      type=Path, required=True)
    p.add_argument("--eps",          type=float, default=1e-12)
    p.add_argument("--min-isoforms", type=int,   default=2)
    p.add_argument("--max-genes",    type=int,   default=None,
                   help="Cap for debugging (process only first N genes)")

    # I/O tuning
    p.add_argument("--row-block-size",  type=int, default=512)
    p.add_argument("--gene-batch-size", type=int, default=500)
    p.add_argument("--checkpoint-every",type=int, default=500)

    # Expression threshold
    p.add_argument("--expr-threshold",            type=float, default=None)
    p.add_argument("--expr-threshold-percentile", type=float, default=10.0)

    # QC filter thresholds (map to QCThresholds fields)
    p.add_argument("--min-valid-samples",     type=int,   default=50,
                   help="F2: min samples where sum(tx) > 0")
    p.add_argument("--max-median-ratio",      type=float, default=1.5,
                   help="F3: max median(sum_tx / gene) — catches missing gene data")
    p.add_argument("--min-corr",              type=float, default=0.8,
                   help="F4: min Pearson corr between gene and sum(tx)")
    p.add_argument("--max-dominant-fraction", type=float, default=0.9,
                   help="F5: max mean dominant isoform fraction")
    p.add_argument("--min-active-samples",    type=int,   default=20,
                   help="F5: min active samples needed to evaluate dominance")

    return p.parse_args()


def main() -> None:
    args = parse_args()

    for path in [args.gene_h5ad, args.tx_h5ad, args.mapping_json]:
        if not path.exists():
            raise FileNotFoundError(f"Not found: {path}")

    args.out_csv.parent.mkdir(parents=True, exist_ok=True)

    thresholds = QCThresholds(
        min_valid_samples=args.min_valid_samples,
        max_median_ratio=args.max_median_ratio,
        min_corr=args.min_corr,
        max_dominant_fraction=args.max_dominant_fraction,
        min_active_samples=args.min_active_samples,
    )

    print(f"Loading gene AnnData:        {args.gene_h5ad}")
    adata_gene = ad.read_h5ad(args.gene_h5ad, backed="r")

    print(f"Loading transcript AnnData:  {args.tx_h5ad}")
    adata_tx = ad.read_h5ad(args.tx_h5ad, backed="r")

    print(f"Loading mapping JSON:         {args.mapping_json}")
    gene_to_tx = load_mapping_json(args.mapping_json)

    print("Validating sample alignment...")
    validate_obs_alignment(adata_gene, adata_tx)

    print("Computing gene-level QC table...")
    qc_df = compute_gene_qc_table_backed(
        adata_gene=adata_gene,
        adata_tx=adata_tx,
        gene_to_tx=gene_to_tx,
        thresholds=thresholds,
        min_isoforms=args.min_isoforms,
        eps=args.eps,
        row_block_size=args.row_block_size,
        gene_batch_size=args.gene_batch_size,
        max_genes=args.max_genes,
        expr_threshold=args.expr_threshold,
        expr_threshold_percentile=args.expr_threshold_percentile,
        out_csv=args.out_csv,
        checkpoint_every=args.checkpoint_every,
    )

    sort_cols = [c for c in ["n_isoforms", "corr_gene_vs_sum_iso"] if c in qc_df.columns]
    if sort_cols:
        qc_df = qc_df.sort_values(by=sort_cols, ascending=[False]*len(sort_cols)).reset_index(drop=True)

    qc_df.to_csv(args.out_csv, index=False)
    print(f"\nSaved: {args.out_csv}  ({len(qc_df)} genes passed QC)")
    print("Done.")


if __name__ == "__main__":
    main()