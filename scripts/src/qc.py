import json
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional

import numpy as np
import pandas as pd


def load_mapping_json(mapping_json: Path) -> Dict[str, List[str]]:
    with open(mapping_json, "r") as f:
        mapping = json.load(f)

    if not isinstance(mapping, dict):
        raise ValueError("Mapping JSON must contain a dictionary: gene_id -> list[transcript_id]")

    return mapping


def validate_obs_alignment(adata_gene, adata_tx) -> None:
    if adata_gene.n_obs != adata_tx.n_obs:
        raise ValueError(
            f"Gene and transcript matrices have different n_obs: "
            f"{adata_gene.n_obs} vs {adata_tx.n_obs}"
        )

    if not np.array_equal(adata_gene.obs_names.to_numpy(), adata_tx.obs_names.to_numpy()):
        raise ValueError(
            "Gene and transcript AnnData objects do not have identical obs_names order. "
            "You must align samples before computing QC."
        )


def _to_numpy_1d(x) -> np.ndarray:
    if hasattr(x, "toarray"):
        x = x.toarray()
    return np.asarray(x).reshape(-1)


def _to_numpy_2d(x) -> np.ndarray:
    if hasattr(x, "toarray"):
        x = x.toarray()
    return np.asarray(x)


@dataclass
class GeneRecord:
    gene_id: str
    gene_idx: int
    tx_ids: List[str]
    tx_idx: List[int]
    n_isoforms: int


from typing import Dict
import numpy as np


class GeneAccumulator:
    def __init__(self, gene_id: str, n_isoforms: int):
        self.gene_id = gene_id
        self.n_isoforms = n_isoforms

        self.n_samples = 0
        self.n_valid_samples = 0

        self.n_gene_detected = 0
        self.n_sum_tx_detected = 0

        # Store per-sample vectors needed for robust final metrics
        self.g_values = []
        self.sum_tx_values = []
        self.rel_error_values = []
        self.ratio_values = []
        self.entropy_values = []
        self.dom_frac_values = []

        # For mean isoform proportion across valid samples
        self.sum_tx_props = np.zeros(n_isoforms, dtype=np.float64)

    def update(self, g_expr: np.ndarray, tx_matrix: np.ndarray, eps: float) -> None:
        g_expr = np.asarray(g_expr).reshape(-1)
        tx_matrix = np.asarray(tx_matrix)

        if tx_matrix.ndim != 2:
            raise ValueError(
                f"tx_matrix must be 2D for gene {self.gene_id}, got shape {tx_matrix.shape}"
            )
        if tx_matrix.shape[0] != g_expr.shape[0]:
            raise ValueError(
                f"Mismatched samples for gene {self.gene_id}: "
                f"g_expr has {g_expr.shape[0]}, tx_matrix has {tx_matrix.shape[0]}"
            )

        sum_tx_expr = tx_matrix.sum(axis=1)

        self.n_samples += g_expr.shape[0]
        self.n_valid_samples += int(np.sum(sum_tx_expr > 0))

        self.n_gene_detected += int(np.sum(g_expr > 0))
        self.n_sum_tx_detected += int(np.sum(sum_tx_expr > 0))

        # Keep original signal for robust final correlation and ratio summaries
        self.g_values.append(g_expr.astype(np.float32, copy=False))
        self.sum_tx_values.append(sum_tx_expr.astype(np.float32, copy=False))

        # Symmetric relative error:
        # 0 means perfect match, bounded near 1 when one is much smaller than the other
        rel_error = np.abs(g_expr - sum_tx_expr) / (np.abs(g_expr) + np.abs(sum_tx_expr) + eps)
        self.rel_error_values.append(rel_error.astype(np.float32, copy=False))

        # Identity ratio: sum_tx / gene
        # This is the key metric for checking gene ≈ sum(transcripts)
        ratio = sum_tx_expr / (g_expr + eps)
        self.ratio_values.append(ratio.astype(np.float32, copy=False))

        valid = sum_tx_expr > 0
        if np.any(valid):
            tx_props = tx_matrix[valid] / (sum_tx_expr[valid, None] + eps)
            self.sum_tx_props += tx_props.sum(axis=0)

            dom_frac = tx_props.max(axis=1)
            self.dom_frac_values.append(dom_frac.astype(np.float32, copy=False))

            if tx_props.shape[1] < 2:
                entropy = np.zeros(tx_props.shape[0], dtype=np.float32)
            else:
                entropy = (
                    -(tx_props * np.log(tx_props + eps)).sum(axis=1) / np.log(tx_props.shape[1])
                ).astype(np.float32, copy=False)

            self.entropy_values.append(entropy)

    def finalize(self) -> Dict[str, float]:
        if self.n_samples == 0:
            raise ValueError(f"No samples accumulated for gene {self.gene_id}")

        g_all = np.concatenate(self.g_values) if self.g_values else np.array([], dtype=np.float32)
        s_all = np.concatenate(self.sum_tx_values) if self.sum_tx_values else np.array([], dtype=np.float32)
        rel_error = (
            np.concatenate(self.rel_error_values)
            if self.rel_error_values
            else np.array([], dtype=np.float32)
        )
        ratio = (
            np.concatenate(self.ratio_values)
            if self.ratio_values
            else np.array([], dtype=np.float32)
        )
        entropy = (
            np.concatenate(self.entropy_values)
            if self.entropy_values
            else np.array([], dtype=np.float32)
        )
        dom_frac = (
            np.concatenate(self.dom_frac_values)
            if self.dom_frac_values
            else np.array([], dtype=np.float32)
        )

        mean_g = float(np.mean(g_all)) if g_all.size else np.nan
        mean_s = float(np.mean(s_all)) if s_all.size else np.nan
        mean_abs_error = float(np.mean(np.abs(g_all - s_all))) if g_all.size else np.nan

        if g_all.size == 0 or s_all.size == 0 or np.std(g_all) == 0 or np.std(s_all) == 0:
            corr = np.nan
        else:
            corr = float(np.corrcoef(g_all, s_all)[0, 1])

        if ratio.size:
            mean_ratio = float(np.mean(ratio))
            median_ratio = float(np.median(ratio))
            p10_ratio = float(np.percentile(ratio, 10))
            p90_ratio = float(np.percentile(ratio, 90))
        else:
            mean_ratio = np.nan
            median_ratio = np.nan
            p10_ratio = np.nan
            p90_ratio = np.nan

        if self.n_valid_samples > 0:
            mean_tx_props = self.sum_tx_props / self.n_valid_samples
            max_mean_isoform_proportion = float(np.max(mean_tx_props))

            mean_entropy = float(np.mean(entropy))
            std_entropy = float(np.std(entropy))
            p10_entropy = float(np.percentile(entropy, 10))
            p90_entropy = float(np.percentile(entropy, 90))

            mean_dom = float(np.mean(dom_frac))
            std_dom = float(np.std(dom_frac))
            p90_dom = float(np.percentile(dom_frac, 90))
        else:
            max_mean_isoform_proportion = np.nan
            mean_entropy = np.nan
            std_entropy = np.nan
            p10_entropy = np.nan
            p90_entropy = np.nan
            mean_dom = np.nan
            std_dom = np.nan
            p90_dom = np.nan

        return {
            "gene_id": self.gene_id,
            "n_isoforms": int(self.n_isoforms),
            "n_valid_samples": int(self.n_valid_samples),
            "mean_gene_expr": mean_g,
            "mean_sum_iso_expr": mean_s,
            "mean_abs_error": mean_abs_error,
            "corr_gene_vs_sum_iso": corr,
            "median_rel_error": float(np.median(rel_error)) if rel_error.size else np.nan,
            "p90_rel_error": float(np.percentile(rel_error, 90)) if rel_error.size else np.nan,
            "mean_ratio_sumTx_to_gene": mean_ratio,
            "median_ratio_sumTx_to_gene": median_ratio,
            "p10_ratio_sumTx_to_gene": p10_ratio,
            "p90_ratio_sumTx_to_gene": p90_ratio,
            "fraction_gene_detected": float(self.n_gene_detected / self.n_samples),
            "fraction_sum_iso_detected": float(self.n_sum_tx_detected / self.n_samples),
            "max_mean_isoform_proportion": max_mean_isoform_proportion,
            "mean_dominant_fraction": mean_dom,
            "std_dominant_fraction": std_dom,
            "p90_dominant_fraction": p90_dom,
            "mean_entropy": mean_entropy,
            "std_entropy": std_entropy,
            "p10_entropy": p10_entropy,
            "p90_entropy": p90_entropy,
        }


def build_gene_records(
    adata_gene,
    adata_tx,
    gene_to_tx: Dict[str, List[str]],
    min_isoforms: int = 2,
) -> List[GeneRecord]:
    gene_name_to_idx = {name: i for i, name in enumerate(adata_gene.var_names)}
    tx_name_to_idx = {name: i for i, name in enumerate(adata_tx.var_names)}

    records: List[GeneRecord] = []

    for gene_id, tx_ids in gene_to_tx.items():
        if gene_id not in gene_name_to_idx:
            continue

        valid_tx_ids = [tx for tx in tx_ids if tx in tx_name_to_idx]
        if len(valid_tx_ids) < min_isoforms:
            continue

        records.append(
            GeneRecord(
                gene_id=gene_id,
                gene_idx=gene_name_to_idx[gene_id],
                tx_ids=valid_tx_ids,
                tx_idx=[tx_name_to_idx[tx] for tx in valid_tx_ids],
                n_isoforms=len(valid_tx_ids),
            )
        )

    if not records:
        raise ValueError(
            "No valid genes after filtering. Check mapping overlap and min_isoforms."
        )

    return records


def compute_gene_qc_table_backed(
    adata_gene,
    adata_tx,
    gene_to_tx: Dict[str, List[str]],
    min_isoforms: int = 2,
    eps: float = 1e-12,
    row_block_size: int = 256,
    max_genes: Optional[int] = None,
) -> pd.DataFrame:
    records = build_gene_records(
        adata_gene=adata_gene,
        adata_tx=adata_tx,
        gene_to_tx=gene_to_tx,
        min_isoforms=min_isoforms,
    )

    if max_genes is not None:
        records = records[:max_genes]

    accs = {
        rec.gene_id: GeneAccumulator(rec.gene_id, rec.n_isoforms)
        for rec in records
    }

    gene_idx_all = np.array([rec.gene_idx for rec in records], dtype=int)

    n_obs = adata_gene.n_obs
    n_genes = len(records)

    gene_id_to_block_col = {
        rec.gene_id: j for j, rec in enumerate(records)
    }

    print(f"Valid genes to process: {n_genes}")
    print(f"Samples: {n_obs}")
    print(f"Row block size: {row_block_size}")

    for start in range(0, n_obs, row_block_size):
        end = min(start + row_block_size, n_obs)

        G_block = _to_numpy_2d(adata_gene[start:end, gene_idx_all].X)

        for rec in records:
            g_expr = G_block[:, gene_id_to_block_col[rec.gene_id]]

            tx_block = _to_numpy_2d(adata_tx[start:end, rec.tx_idx].X)

            accs[rec.gene_id].update(
                g_expr=g_expr,
                tx_matrix=tx_block,
                eps=eps,
            )

        print(f"Processed rows {start}:{end}")

    rows = [accs[rec.gene_id].finalize() for rec in records]
    qc_df = pd.DataFrame(rows)

    if qc_df.empty:
        raise ValueError("QC table is empty after processing.")

    return qc_df