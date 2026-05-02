"""
qc.py — gene-level QC for isoform prediction datasets
======================================================

PIPELINE OVERVIEW
-----------------
We have two AnnData files (backed/HDF5, too large for RAM):
  - adata_gene : (n_samples × n_genes)   — bulk gene expression
  - adata_tx   : (n_samples × n_transcripts) — bulk transcript expression

Both share the same rows (samples). The mapping JSON tells us which
transcripts belong to each gene:
  { "GENE_A": ["TX1", "TX2", "TX3"], ... }

Goal: for each gene, compute summary statistics comparing gene expression
vs. the sum of its transcripts, then decide if the gene is usable for
isoform proportion prediction.

FILTERS (applied in order — gene is discarded at first failure)
---------------------------------------------------------------
F1  min_isoforms      ≥ 2          (pre-scan, no data read needed)
F2  n_valid_samples   ≥ 50         (enough samples with sum_tx > 0)
F3  median_ratio      ≤ 1.5        (sum_tx / gene not suspiciously large)
                                    catches genes with missing/zero gene
                                    expression but present transcript data
F4  corr(gene, sum_tx) ≥ 0.8      (gene and isoforms track each other)
F5  mean_dominant_fraction ≤ 0.9  (not trivially single-isoform)

F2–F5 are evaluated in finalize() after reading all samples, but the
function returns None for discarded genes so we never store their data.

I/O STRATEGY (why this is fast)
--------------------------------
Naive approach: for each gene, for each row-block → 1 HDF5 read
  → O(n_genes × n_blocks) reads  ≈ 10 million reads for 19K genes

Here: group genes into batches, collect ALL column indices, do ONE read
per row-block that covers all genes in the batch:
  → O(n_batches × n_blocks) reads  ≈ 38 × 550 = 21K reads  (~500× faster)

CHECKPOINTING
-------------
Every `checkpoint_every` completed genes, results are flushed to a
.checkpoint.csv next to out_csv. On restart, genes already there are
skipped — the job resumes where it left off.
"""

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd


# ══════════════════════════════════════════════════════════════════════════════
# I/O helpers
# ══════════════════════════════════════════════════════════════════════════════

def load_mapping_json(mapping_json: Path) -> Dict[str, List[str]]:
    """Load gene → [transcript, ...] mapping from JSON."""
    with open(mapping_json, "r") as f:
        mapping = json.load(f)
    if not isinstance(mapping, dict):
        raise ValueError("Mapping JSON must be a dict: gene_id -> list[transcript_id]")
    return mapping


def validate_obs_alignment(adata_gene, adata_tx) -> None:
    """Ensure gene and transcript AnnData have the same samples in the same order."""
    if adata_gene.n_obs != adata_tx.n_obs:
        raise ValueError(
            f"n_obs mismatch: gene={adata_gene.n_obs}, tx={adata_tx.n_obs}"
        )
    if not np.array_equal(adata_gene.obs_names.to_numpy(), adata_tx.obs_names.to_numpy()):
        raise ValueError(
            "obs_names differ between gene and transcript AnnData. Align samples first."
        )


def _to_dense(x) -> np.ndarray:
    """Convert sparse or backed array to a dense numpy array."""
    if hasattr(x, "toarray"):
        return x.toarray()
    return np.asarray(x)


# ══════════════════════════════════════════════════════════════════════════════
# Data structures
# ══════════════════════════════════════════════════════════════════════════════

@dataclass
class GeneRecord:
    """Lightweight descriptor for one gene and its transcripts."""
    gene_id:    str
    gene_idx:   int        # column index in adata_gene
    tx_ids:     List[str]
    tx_idx:     List[int]  # column indices in adata_tx
    n_isoforms: int


# ──────────────────────────────────────────────────────────────────────────────
# QC thresholds dataclass — keeps all filter parameters in one place so they
# are easy to pass around and document.
# ──────────────────────────────────────────────────────────────────────────────
@dataclass
class QCThresholds:
    """
    All filter thresholds in one place.

    F2  min_valid_samples      : samples where sum(tx) > 0
    F3  max_median_ratio       : max allowed median(sum_tx / gene)
                                 > 1.5 means isoforms >> gene → bad mapping
    F4  min_corr               : Pearson r between gene and sum(tx)
    F5  max_dominant_fraction  : max allowed mean fraction of dominant isoform
    F5  min_active_samples     : min samples with sum_tx > expr_threshold
                                 needed to trust F5
    """
    min_valid_samples:     int   = 50
    max_median_ratio:      float = 1.5
    min_corr:              float = 0.8
    max_dominant_fraction: float = 0.9
    min_active_samples:    int   = 20


# ──────────────────────────────────────────────────────────────────────────────
# GeneAccumulator — collects statistics across row-blocks for one gene,
# then decides in finalize() whether the gene passes QC.
# ──────────────────────────────────────────────────────────────────────────────

# Sentinel returned by finalize() when a gene fails QC early.
# Using None would work too, but a named constant is clearer at the call site.
_GENE_FAILED = None


class GeneAccumulator:
    """
    Stateful accumulator for one gene.

    Lifecycle
    ---------
    1. __init__  : allocate accumulators
    2. update()  : called once per row-block (chunk of samples)
    3. finalize(): concatenate everything, compute final metrics,
                   apply QC filters → return dict or None (failed)

    Why accumulate instead of computing per-block?
    ----------------------------------------------
    Pearson correlation and percentiles need all samples at once.
    Keeping float32 arrays per block is memory-cheap (a few MB per gene
    at most) and avoids two passes over the data.
    """

    def __init__(self, gene_id: str, n_isoforms: int):
        self.gene_id    = gene_id
        self.n_isoforms = n_isoforms

        # Counters (scalars, updated per block)
        self.n_samples          = 0
        self.n_valid_samples    = 0   # samples where sum_tx > 0
        self.n_gene_detected    = 0   # samples where gene > 0
        self.n_sum_tx_detected  = 0   # (same as n_valid_samples, kept for symmetry)

        # Raw per-sample vectors (float32 lists, concatenated in finalize)
        self.g_values:        List[np.ndarray] = []
        self.sum_tx_values:   List[np.ndarray] = []
        self.rel_error_values:List[np.ndarray] = []
        self.ratio_values:    List[np.ndarray] = []

        # Computed only for "active" samples (sum_tx > expr_threshold)
        self.entropy_values:  List[np.ndarray] = []
        self.dom_frac_values: List[np.ndarray] = []
        self.sum_tx_props = np.zeros(n_isoforms, dtype=np.float64)  # for mean isoform props
        self.n_active_accumulated = 0

    def update(
        self,
        g_expr:   np.ndarray,   # shape (block_size,)
        tx_matrix: np.ndarray,  # shape (block_size, n_isoforms)
        eps:            float,
        expr_threshold: float = 0.0,
    ) -> None:
        """
        Ingest one row-block of data.

        g_expr    : gene expression for this block of samples
        tx_matrix : transcript expression matrix for this block,
                    columns = isoforms of this gene
        eps       : small constant to avoid log(0) and division by zero
        expr_threshold : sum_tx must exceed this to count as "active"
                         (used for entropy/dominance only)
        """
        sum_tx = tx_matrix.sum(axis=1)   # (block_size,) — total isoform expression

        # ── Counters ──────────────────────────────────────────────────────────
        self.n_samples         += g_expr.shape[0]
        self.n_valid_samples   += int(np.sum(sum_tx > 0))
        self.n_gene_detected   += int(np.sum(g_expr > 0))
        self.n_sum_tx_detected += int(np.sum(sum_tx > 0))

        # ── Store raw signals (needed for correlation and ratio percentiles) ──
        self.g_values.append(g_expr.astype(np.float32, copy=False))
        self.sum_tx_values.append(sum_tx.astype(np.float32, copy=False))

        # Symmetric relative error ∈ [0, 1]:
        #   0 = perfect match, ~1 = one is much larger than the other
        rel = np.abs(g_expr - sum_tx) / (np.abs(g_expr) + np.abs(sum_tx) + eps)
        self.rel_error_values.append(rel.astype(np.float32, copy=False))

        # Ratio = sum_tx / gene.
        # Healthy range: ~1.  >> 1 means tx >> gene (missing gene data).
        # We use eps on the denominator only (not symmetric) because we
        # specifically want to catch large ratios when gene ≈ 0.
        ratio = sum_tx / (g_expr + eps)
        self.ratio_values.append(ratio.astype(np.float32, copy=False))

        # ── Isoform proportion stats — only for "active" samples ──────────────
        # Using expr_threshold > 0 instead of > 0 prevents near-zero noise
        # samples from dominating the entropy/dominance averages.
        active = sum_tx > expr_threshold
        if np.any(active):
            # Normalize: each row sums to 1 → isoform proportions
            tx_props = tx_matrix[active] / (sum_tx[active, None] + eps)

            self.sum_tx_props         += tx_props.sum(axis=0)   # for mean props
            self.n_active_accumulated += int(active.sum())

            # Dominant fraction = proportion of the most expressed isoform
            self.dom_frac_values.append(tx_props.max(axis=1).astype(np.float32, copy=False))

            if tx_props.shape[1] >= 2:
                # Normalized Shannon entropy ∈ [0,1]:
                #   0 = one isoform takes everything
                #   1 = all isoforms equally expressed
                p = tx_props + eps
                h = -(p * np.log(p)).sum(axis=1) / np.log(tx_props.shape[1])
                self.entropy_values.append(h.astype(np.float32, copy=False))

    def finalize(
        self,
        thresholds: QCThresholds,
        eps: float = 1e-12,
    ) -> Optional[dict]:
        """
        Concatenate accumulated data, compute final metrics, apply QC filters.

        Returns
        -------
        dict   : metrics for this gene (passes all filters)
        None   : gene failed at least one filter → discard, save nothing

        Early-exit logic
        ----------------
        Filters are checked in cheapest-first order.  As soon as one fails
        we return None immediately — no point computing correlation or
        entropy for a gene we are going to throw away anyway.
        """
        g_all  = np.concatenate(self.g_values)
        s_all  = np.concatenate(self.sum_tx_values)
        ratio  = np.concatenate(self.ratio_values)
        rel    = np.concatenate(self.rel_error_values)

        # ── F2: enough valid samples ──────────────────────────────────────────
        # A gene with very few expressed samples gives unreliable stats.
        if self.n_valid_samples < thresholds.min_valid_samples:
            return _GENE_FAILED

        # ── F3: median ratio ≤ max_median_ratio ──────────────────────────────
        # If sum(isoforms) is consistently >> gene, the gene column is likely
        # missing/near-zero while transcript data exists.  This produces huge
        # ratios and is biologically uninterpretable for our task.
        # We use the MEDIAN (not mean) because a few outlier samples with very
        # low gene expression would inflate the mean unfairly.
        median_ratio = float(np.median(ratio))
        if median_ratio > thresholds.max_median_ratio:
            return _GENE_FAILED

        # ── F4: correlation gene vs sum(tx) ──────────────────────────────────
        # A high correlation means gene and isoforms track each other across
        # samples — good signal for learning isoform proportions.
        # NaN correlation (constant gene expression) also fails.
        corr = np.nan
        if g_all.size >= 3 and np.std(g_all) > 1e-8 and np.std(s_all) > 1e-8:
            corr = float(np.corrcoef(g_all, s_all)[0, 1])
        if np.isnan(corr) or corr < thresholds.min_corr:
            return _GENE_FAILED

        # ── F5: isoform diversity ─────────────────────────────────────────────
        # Only checked if we have enough "active" samples for a reliable estimate.
        # If not enough active samples: we give the gene the benefit of the
        # doubt and keep it (conservative), but flag it.
        entropy  = np.concatenate(self.entropy_values)  if self.entropy_values  else np.array([], np.float32)
        dom      = np.concatenate(self.dom_frac_values) if self.dom_frac_values else np.array([], np.float32)
        n_act    = self.n_active_accumulated
        f5_assessed = n_act >= thresholds.min_active_samples

        if f5_assessed and dom.size > 0:
            mean_dom = float(np.mean(dom))
            if mean_dom > thresholds.max_dominant_fraction:
                return _GENE_FAILED
        else:
            mean_dom = float(np.mean(dom)) if dom.size else np.nan

        # ── All filters passed — build output dict ────────────────────────────
        mtp = self.sum_tx_props / n_act if n_act > 0 else None

        def _p(a, q): return float(np.percentile(a, q)) if a.size else np.nan

        return {
            "gene_id":    self.gene_id,
            "n_isoforms": self.n_isoforms,
            "n_samples":  self.n_samples,
            "n_valid_samples":  self.n_valid_samples,
            "n_active_samples": n_act,
            "f5_assessed": f5_assessed,

            # Gene vs sum-of-isoforms agreement
            "mean_gene_expr":      float(np.mean(g_all)),
            "mean_sum_iso_expr":   float(np.mean(s_all)),
            "mean_abs_error":      float(np.mean(np.abs(g_all - s_all))),
            "corr_gene_vs_sum_iso": corr,
            "median_rel_error":    _p(rel, 50),
            "p90_rel_error":       _p(rel, 90),

            # Ratio metrics (F3 key metric)
            "mean_ratio_sumTx_to_gene":   float(np.mean(ratio)),
            "median_ratio_sumTx_to_gene": median_ratio,
            "p10_ratio_sumTx_to_gene":    _p(ratio, 10),
            "p90_ratio_sumTx_to_gene":    _p(ratio, 90),

            # Detection rates
            "fraction_gene_detected":    float(self.n_gene_detected   / self.n_samples),
            "fraction_sum_iso_detected": float(self.n_sum_tx_detected / self.n_samples),

            # Isoform proportion stats (only meaningful when f5_assessed=True)
            "max_mean_isoform_proportion": float(np.max(mtp)) if mtp is not None else np.nan,
            "mean_dominant_fraction":      mean_dom,
            "std_dominant_fraction":       float(np.std(dom))  if dom.size else np.nan,
            "p90_dominant_fraction":       _p(dom, 90),
            "mean_entropy":                float(np.mean(entropy)) if entropy.size else np.nan,
            "std_entropy":                 float(np.std(entropy))  if entropy.size else np.nan,
            "p10_entropy":                 _p(entropy, 10),
            "p90_entropy":                 _p(entropy, 90),
        }


# ══════════════════════════════════════════════════════════════════════════════
# Gene record building
# ══════════════════════════════════════════════════════════════════════════════

def build_gene_records(
    adata_gene,
    adata_tx,
    gene_to_tx:   Dict[str, List[str]],
    min_isoforms: int = 2,
) -> List[GeneRecord]:
    """
    F1 (pre-scan): build GeneRecord list, discarding genes with < min_isoforms
    transcripts present in adata_tx.  No data is read from the matrices here.
    """
    gene_name_to_idx = {name: i for i, name in enumerate(adata_gene.var_names)}
    tx_name_to_idx   = {name: i for i, name in enumerate(adata_tx.var_names)}
    records = []
    skipped_missing = 0
    skipped_few_tx  = 0

    for gene_id, tx_ids in gene_to_tx.items():
        if gene_id not in gene_name_to_idx:
            skipped_missing += 1
            continue
        valid_tx = [tx for tx in tx_ids if tx in tx_name_to_idx]
        if len(valid_tx) < min_isoforms:
            skipped_few_tx += 1
            continue
        records.append(GeneRecord(
            gene_id=gene_id,
            gene_idx=gene_name_to_idx[gene_id],
            tx_ids=valid_tx,
            tx_idx=[tx_name_to_idx[tx] for tx in valid_tx],
            n_isoforms=len(valid_tx),
        ))

    print(
        f"build_gene_records: {len(records)} valid genes "
        f"(skipped: {skipped_missing} not in gene matrix, "
        f"{skipped_few_tx} with <{min_isoforms} transcripts)"
    )
    if not records:
        raise ValueError("No valid genes after F1 filter. Check mapping overlap and min_isoforms.")
    return records


# ══════════════════════════════════════════════════════════════════════════════
# Expression threshold for active-sample detection
# ══════════════════════════════════════════════════════════════════════════════

def compute_global_expr_threshold(
    adata_tx,
    gene_records:    List[GeneRecord],
    row_block_size:  int   = 512,
    percentile:      float = 10.0,
    max_genes_sample:int   = 200,
) -> float:
    """
    Estimate a global expression threshold = p-th percentile of non-zero
    sum_tx values, computed over a random subset of genes.

    Why not use > 0?
    ----------------
    Sequencing data has many near-zero values from noise.  A sample with
    sum_tx = 0.001 is effectively not expressing the gene, but would pass
    the > 0 filter and contribute a trivially dominated (entropy = 0)
    observation to the dominance/entropy averages.  Using the p10 of
    non-zero values removes this noise floor.

    We sample max_genes_sample genes (default 200) for speed — the
    percentile estimate is stable with a few thousand values.
    """
    rng = np.random.default_rng(42)
    subset = gene_records
    if len(gene_records) > max_genes_sample:
        idx    = rng.choice(len(gene_records), size=max_genes_sample, replace=False)
        subset = [gene_records[i] for i in sorted(idx)]

    # One combined HDF5 slice per row-block (same batched-I/O trick)
    all_tx_idx = sorted(set(i for rec in subset for i in rec.tx_idx))
    tx_local   = {g: i for i, g in enumerate(all_tx_idx)}
    sums = []

    for start in range(0, adata_tx.n_obs, row_block_size):
        end   = min(start + row_block_size, adata_tx.n_obs)
        block = _to_dense(adata_tx[start:end, all_tx_idx].X)
        for rec in subset:
            lc = [tx_local[i] for i in rec.tx_idx]
            s  = block[:, lc].sum(axis=1)
            nz = s[s > 0]
            if nz.size:
                sums.append(nz.astype(np.float32))

    if not sums:
        print("Warning: no non-zero sum_tx found in threshold sample — using 0.0")
        return 0.0

    threshold = float(np.percentile(np.concatenate(sums), percentile))
    print(f"Global expr threshold (p{percentile:.0f} of non-zero sum_tx): {threshold:.4f}")
    return threshold


# ══════════════════════════════════════════════════════════════════════════════
# Checkpointing
# ══════════════════════════════════════════════════════════════════════════════

def _ckpt_path(out_csv: Path) -> Path:
    return out_csv.with_suffix(".checkpoint.csv")


def _load_checkpoint(ckpt: Path) -> Optional[pd.DataFrame]:
    """Load checkpoint CSV if it exists, return None otherwise."""
    if ckpt.exists():
        df = pd.read_csv(ckpt)
        print(f"Checkpoint: {len(df)} genes already done — skipping them.")
        return df
    return None


def _save_checkpoint(df: pd.DataFrame, ckpt: Path) -> None:
    """
    Atomic checkpoint write: write to .tmp then rename.
    If the job is killed mid-write, the old checkpoint survives intact.
    """
    tmp = ckpt.with_suffix(".tmp")
    df.to_csv(tmp, index=False)
    tmp.replace(ckpt)


# ══════════════════════════════════════════════════════════════════════════════
# Main QC computation
# ══════════════════════════════════════════════════════════════════════════════

def compute_gene_qc_table_backed(
    adata_gene,
    adata_tx,
    gene_to_tx:   Dict[str, List[str]],
    thresholds:   Optional[QCThresholds] = None,
    min_isoforms: int   = 2,
    eps:          float = 1e-12,
    row_block_size:  int = 512,
    gene_batch_size: int = 500,
    max_genes:    Optional[int]  = None,
    expr_threshold: Optional[float] = None,
    expr_threshold_percentile: float = 10.0,
    out_csv:      Optional[Path] = None,
    checkpoint_every: int = 500,
) -> pd.DataFrame:
    """
    Main entry point. See module docstring for full explanation.

    Parameters
    ----------
    thresholds       : QCThresholds instance (uses defaults if None)
    gene_batch_size  : genes per I/O batch — tune to available RAM
                       (500 genes × 512 samples × ~100 tx/gene ≈ 25 MB)
    checkpoint_every : flush checkpoint after this many completed genes
    """
    if thresholds is None:
        thresholds = QCThresholds()

    # ── F1: pre-scan (no matrix I/O) ──────────────────────────────────────────
    all_records = build_gene_records(adata_gene, adata_tx, gene_to_tx, min_isoforms)
    if max_genes is not None:
        all_records = all_records[:max_genes]

    # ── Resume from checkpoint ─────────────────────────────────────────────────
    ckpt     = _ckpt_path(out_csv) if out_csv else None
    done_df  = _load_checkpoint(ckpt) if ckpt else None
    done_ids = set(done_df["gene_id"]) if done_df is not None else set()
    records  = [r for r in all_records if r.gene_id not in done_ids]
    print(f"Genes to process : {len(records)}  (checkpoint skip: {len(done_ids)})")

    if not records:
        print("All genes already in checkpoint — done.")
        return done_df if done_df is not None else pd.DataFrame()

    # ── Expression threshold for active-sample detection ───────────────────────
    if expr_threshold is None:
        expr_threshold = compute_global_expr_threshold(
            adata_tx, all_records,
            row_block_size=row_block_size,
            percentile=expr_threshold_percentile,
        )

    n_obs     = adata_gene.n_obs
    n_batches = (len(records) + gene_batch_size - 1) // gene_batch_size
    completed: List[dict] = []
    n_discarded = 0

    print(f"Gene-batches     : {n_batches} × ≤{gene_batch_size} genes")
    print(f"Row-block size   : {row_block_size}  |  samples: {n_obs}")
    print(f"Expr threshold   : {expr_threshold:.4f}")
    print(f"QC thresholds    : {thresholds}")

    for bi, bs in enumerate(range(0, len(records), gene_batch_size)):
        batch = records[bs : bs + gene_batch_size]

        # ── Build unified column-index maps for this batch ─────────────────────
        # Collecting all gene/tx indices lets us do ONE HDF5 slice per row-block
        # instead of one per gene.
        gene_cols  = sorted(set(rec.gene_idx for rec in batch))
        gene_local = {g: i for i, g in enumerate(gene_cols)}
        tx_cols    = sorted(set(i for rec in batch for i in rec.tx_idx))
        tx_local   = {t: i for i, t in enumerate(tx_cols)}

        accs = {rec.gene_id: GeneAccumulator(rec.gene_id, rec.n_isoforms) for rec in batch}

        # ── Row-block loop — ONE gene read + ONE tx read per block ─────────────
        for rs in range(0, n_obs, row_block_size):
            re = min(rs + row_block_size, n_obs)

            # Single HDF5 read covering all genes in this batch
            G = _to_dense(adata_gene[rs:re, gene_cols].X)  # (block, n_gene_cols)
            # Single HDF5 read covering all transcripts in this batch
            T = _to_dense(adata_tx[rs:re,   tx_cols].X)    # (block, n_tx_cols)

            for rec in batch:
                g_expr = G[:, gene_local[rec.gene_idx]]
                tx_mat = T[:, [tx_local[i] for i in rec.tx_idx]]
                accs[rec.gene_id].update(g_expr, tx_mat, eps, expr_threshold)

        # ── Finalize batch — apply filters, discard failures ───────────────────
        # finalize() returns None if the gene fails any filter (early exit).
        # We never store None entries — this saves memory and checkpoint space.
        for rec in batch:
            result = accs[rec.gene_id].finalize(thresholds=thresholds, eps=eps)
            if result is _GENE_FAILED:
                n_discarded += 1
            else:
                completed.append(result)

        # ── Checkpoint flush ───────────────────────────────────────────────────
        if ckpt and (len(completed) % checkpoint_every < gene_batch_size
                     or bi == n_batches - 1):
            new_df   = pd.DataFrame(completed)
            combined = pd.concat([done_df, new_df], ignore_index=True) if done_df is not None else new_df
            _save_checkpoint(combined, ckpt)
            total_done = len(done_ids) + len(completed)
            print(
                f"  [batch {bi+1}/{n_batches}] checkpoint: "
                f"{total_done} passed / {n_discarded} discarded so far"
            )
        else:
            print(
                f"  Gene-batch {bi+1}/{n_batches} — "
                f"passed: {len(completed)}, discarded: {n_discarded}"
            )

    new_df = pd.DataFrame(completed)
    final  = pd.concat([done_df, new_df], ignore_index=True) if done_df is not None else new_df

    print(f"\nFinal: {len(final)} genes passed QC out of {len(all_records)} candidates")
    if final.empty:
        raise ValueError("No genes passed QC. Check thresholds and data.")
    return final


# ══════════════════════════════════════════════════════════════════════════════
# Post-hoc filter summary (optional — for exploration in notebooks)
# ══════════════════════════════════════════════════════════════════════════════

def summarize_qc(qc_df: pd.DataFrame, thresholds: Optional[QCThresholds] = None) -> pd.DataFrame:
    """
    Add per-filter boolean columns to an already-computed QC table.
    Useful in notebooks to explore which threshold catches how many genes.
    Note: genes that failed during compute() are already absent from qc_df.
    This function is for post-hoc exploration only.
    """
    if thresholds is None:
        thresholds = QCThresholds()
    df = qc_df.copy()
    df["pass_f2"] = df["n_valid_samples"] >= thresholds.min_valid_samples
    df["pass_f3"] = df["median_ratio_sumTx_to_gene"] <= thresholds.max_median_ratio
    df["pass_f4"] = df["corr_gene_vs_sum_iso"].notna() & (df["corr_gene_vs_sum_iso"] >= thresholds.min_corr)
    has_enough    = df["n_active_samples"] >= thresholds.min_active_samples
    df["pass_f5"] = ~has_enough | (df["mean_dominant_fraction"] <= thresholds.max_dominant_fraction)
    df["pass_all"] = df["pass_f2"] & df["pass_f3"] & df["pass_f4"] & df["pass_f5"]
    n = len(df)
    for col, label in [
        ("pass_f2", f"F2 n_valid_samples ≥ {thresholds.min_valid_samples}"),
        ("pass_f3", f"F3 median_ratio ≤ {thresholds.max_median_ratio}"),
        ("pass_f4", f"F4 corr ≥ {thresholds.min_corr}"),
        ("pass_f5", f"F5 dom_frac ≤ {thresholds.max_dominant_fraction}"),
        ("pass_all","ALL filters"),
    ]:
        print(f"  {label}: {df[col].sum()} / {n} ({100*df[col].sum()/n:.1f}%)")
    return df