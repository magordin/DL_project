"""
qc.py — gene-level QC for isoform prediction datasets
======================================================

PIPELINE OVERVIEW
-----------------
We have two AnnData files (backed/HDF5, too large for RAM):
  - adata_gene : (n_samples × n_genes)      — bulk gene expression
  - adata_tx   : (n_samples × n_transcripts) — bulk transcript expression

Both share the same rows (samples). The mapping JSON tells us which
transcripts belong to each gene:
  { "GENE_A": ["TX1", "TX2", "TX3"], ... }

Goal: for each gene, compute summary statistics comparing gene expression
vs. the sum of its transcripts, then decide if the gene is usable for
isoform proportion prediction.

FILTERS (applied in order inside finalize())
--------------------------------------------
F1  min_isoforms         ≥ 2     pre-scan, no data read needed
F2  n_both_detected      ≥ 50    samples where BOTH gene > 0 AND sum_tx > 0
                                  (replaces the old sum_tx > 0 only check)
F3  median_ratio         ≤ 1.5   median(sum_tx / gene) over both-detected
                                  samples; catches genes with missing gene
                                  data but present transcript signal
F4  corr(gene, sum_tx)   ≥ 0.8   Pearson over all samples; gene and isoforms
                                  must co-vary
F5  mean_dominant_frac   ≤ 0.9   checked only when n_active ≥ min_active;
                                  if too few active samples → gene FAILS
                                  (not silently passed — see note below)

NOTE on F5 "benefit of the doubt"
----------------------------------
Previous versions kept genes with too few active samples. That is permissive,
not conservative: if we cannot assess isoform diversity, the gene is not safe
for training. We now FAIL such genes explicitly with reason "F5_low_active".

CHECKPOINTING (fixed)
---------------------
Both PASSED and FAILED genes are written to the checkpoint, with columns:
  qc_pass     : bool
  fail_reason : str  (e.g. "F3_median_ratio") or "" for passing genes

This means restarted jobs skip ALL previously seen genes, not just passing ones.
The final output CSV contains all genes with qc_pass column so the caller
can decide what to do with borderline cases.
"""

import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional

import numpy as np
import pandas as pd


# ══════════════════════════════════════════════════════════════════════════════
# I/O helpers
# ══════════════════════════════════════════════════════════════════════════════

def load_mapping_json(mapping_json: Path) -> Dict[str, List[str]]:
    with open(mapping_json, "r") as f:
        mapping = json.load(f)
    if not isinstance(mapping, dict):
        raise ValueError("Mapping JSON must be a dict: gene_id -> list[transcript_id]")
    return mapping


def validate_obs_alignment(adata_gene, adata_tx) -> None:
    if adata_gene.n_obs != adata_tx.n_obs:
        raise ValueError(f"n_obs mismatch: gene={adata_gene.n_obs}, tx={adata_tx.n_obs}")
    if not np.array_equal(adata_gene.obs_names.to_numpy(), adata_tx.obs_names.to_numpy()):
        raise ValueError("obs_names differ between gene and transcript AnnData.")


def _to_dense(x) -> np.ndarray:
    if hasattr(x, "toarray"):
        return x.toarray()
    return np.asarray(x)


# ══════════════════════════════════════════════════════════════════════════════
# Data structures
# ══════════════════════════════════════════════════════════════════════════════

@dataclass
class GeneRecord:
    gene_id:    str
    gene_idx:   int
    tx_ids:     List[str]
    tx_idx:     List[int]
    n_isoforms: int


@dataclass
class QCThresholds:
    """
    All filter thresholds in one place.

    F2  min_both_detected    : min samples where BOTH gene > 0 AND sum_tx > 0
    F3  max_median_ratio     : max median(sum_tx / gene) over both-detected samples
    F4  min_corr             : min Pearson r(gene, sum_tx) over ALL samples
    F5  max_dominant_frac    : max mean dominant isoform fraction
        min_active_samples   : min samples with sum_tx > expr_threshold to
                               evaluate F5; genes below this threshold FAIL F5
    """
    min_both_detected:  int   = 50
    max_median_ratio:   float = 1.5
    min_corr:           float = 0.8
    max_dominant_frac:  float = 0.9
    min_active_samples: int   = 20


# ══════════════════════════════════════════════════════════════════════════════
# GeneAccumulator
# ══════════════════════════════════════════════════════════════════════════════

class GeneAccumulator:
    """
    Accumulates per-block statistics for one gene, then applies QC filters
    in finalize().

    finalize() ALWAYS returns a dict — never None.
    The dict always contains 'qc_pass' (bool) and 'fail_reason' (str).
    This way EVERY processed gene is recorded in the checkpoint, so restarts
    skip both passing and failing genes.

    Early-exit order (cheapest first):
        F2 → F3 → F4 → F5
    Once a filter fails, we return immediately with partial metrics so the
    caller still knows why the gene was rejected.
    """

    def __init__(self, gene_id: str, n_isoforms: int):
        self.gene_id    = gene_id
        self.n_isoforms = n_isoforms

        self.n_samples         = 0
        self.n_valid_samples   = 0   # sum_tx > 0
        self.n_gene_detected   = 0   # gene > 0
        self.n_both_detected   = 0   # gene > 0 AND sum_tx > 0  ← key for F2/F3
        self.n_double_zero     = 0   # gene = 0 AND sum_tx = 0  ← valid, excluded from metrics

        self.g_values:         List[np.ndarray] = []
        self.sum_tx_values:    List[np.ndarray] = []
        self.rel_error_values: List[np.ndarray] = []
        # ratio only accumulated over both-detected samples (avoids 0/eps blowup)
        self.ratio_values:     List[np.ndarray] = []

        self.entropy_values:   List[np.ndarray] = []
        self.dom_frac_values:  List[np.ndarray] = []
        self.sum_tx_props  = np.zeros(n_isoforms, dtype=np.float64)
        self.n_active_accumulated = 0

    def update(
        self,
        g_expr:         np.ndarray,
        tx_matrix:      np.ndarray,
        eps:            float,
        expr_threshold: float = 0.0,
    ) -> None:
        sum_tx = tx_matrix.sum(axis=1)

        both = (g_expr > 0) & (sum_tx > 0)

        self.n_samples       += g_expr.shape[0]
        self.n_valid_samples += int(np.sum(sum_tx > 0))
        self.n_gene_detected += int(np.sum(g_expr > 0))
        self.n_both_detected += int(np.sum(both))
        self.n_double_zero   += int(np.sum((g_expr == 0) & (sum_tx == 0)))

        # Double-zero mask: gene=0 AND sum_tx=0 — biologically valid (gene not
        # expressed in this sample) but uninformative for every metric:
        #   rel_error : abs(0-0)/(0+0+eps) = 0 trivially → biases median toward 0
        #   ratio     : 0/eps = 0 → misleading (no real signal)
        #   entropy   : tx_matrix all-zero → meaningless proportions
        # They ARE counted in n_samples and excluded from correlation vectors.
        # Samples where exactly one of gene/sum_tx is zero ARE kept for
        # correlation — the discordance is real signal.
        double_zero  = (g_expr == 0) & (sum_tx == 0)
        informative  = ~double_zero   # gene>0 OR sum_tx>0

        # Correlation vectors: informative samples only
        if np.any(informative):
            self.g_values.append(g_expr[informative].astype(np.float32, copy=False))
            self.sum_tx_values.append(sum_tx[informative].astype(np.float32, copy=False))

        # Symmetric relative error over informative samples only.
        # Compute in float64 to avoid float32 precision loss in the denominator.
        # eps=1e-12 truncates to ~1e-12 in float32 but the real risk is
        # catastrophic cancellation when gi ≈ si and both are large.
        if np.any(informative):
            gi = g_expr[informative].astype(np.float64)
            si = sum_tx[informative].astype(np.float64)
            rel = np.abs(gi - si) / (np.abs(gi) + np.abs(si) + eps)
            # Clamp to [0, 1] — mathematically guaranteed but float64 is safe
            rel = np.clip(rel, 0.0, 1.0)
            self.rel_error_values.append(rel.astype(np.float32))

        # Ratio over both-detected samples only (gene>0 AND sum_tx>0).
        # float64 for the division — large sum_tx / small gene can lose bits in float32.
        if np.any(both):
            ratio = sum_tx[both].astype(np.float64) / g_expr[both].astype(np.float64)
            self.ratio_values.append(ratio.astype(np.float32))

        # Entropy / dominance — active samples (sum_tx > threshold).
        # This already excludes double-zeros since sum_tx=0 <= any threshold.
        active = sum_tx > expr_threshold
        if np.any(active):
            tx_props = tx_matrix[active] / (sum_tx[active, None] + eps)
            self.sum_tx_props         += tx_props.sum(axis=0)
            self.n_active_accumulated += int(active.sum())
            self.dom_frac_values.append(tx_props.max(axis=1).astype(np.float32, copy=False))
            if tx_props.shape[1] >= 2:
                p = tx_props + eps
                h = -(p * np.log(p)).sum(axis=1) / np.log(tx_props.shape[1])
                self.entropy_values.append(h.astype(np.float32, copy=False))

    # ── helpers ───────────────────────────────────────────────────────────────

    @staticmethod
    def _p(a: np.ndarray, q: float) -> float:
        return float(np.percentile(a, q)) if a.size else float("nan")

    def _base_metrics(self, g_all, s_all, rel) -> dict:
        """Metrics that are always computed regardless of QC outcome."""
        return {
            "gene_id":             self.gene_id,
            "n_isoforms":          self.n_isoforms,
            "n_samples":           self.n_samples,
            "n_valid_samples":     self.n_valid_samples,
            "n_gene_detected":     self.n_gene_detected,
            "n_both_detected":     self.n_both_detected,
            "mean_gene_expr":      float(np.mean(g_all)),
            "mean_sum_iso_expr":   float(np.mean(s_all)),
            "mean_abs_error":      float(np.mean(np.abs(g_all - s_all))),
            "median_rel_error":    self._p(rel, 50),
            "p90_rel_error":       self._p(rel, 90),
            "fraction_gene_detected":    float(self.n_gene_detected   / self.n_samples),
            "fraction_sum_iso_detected": float(self.n_valid_samples   / self.n_samples),
            "fraction_both_detected":    float(self.n_both_detected   / self.n_samples),
            "n_double_zero":             self.n_double_zero,
            "fraction_double_zero":      float(self.n_double_zero / self.n_samples),
        }

    @staticmethod
    def _failed(base: dict, reason: str) -> dict:
        """Stamp a gene as failed and return early."""
        base.update({
            "qc_pass":   False,
            "fail_reason": reason,
            # Fill remaining columns with NaN so CSV schema is uniform
            "corr_gene_vs_sum_iso":        float("nan"),
            "median_ratio_sumTx_to_gene":  float("nan"),
            "p10_ratio_sumTx_to_gene":     float("nan"),
            "p90_ratio_sumTx_to_gene":     float("nan"),
            "n_active_samples":            0,
            "max_mean_isoform_proportion": float("nan"),
            "mean_dominant_fraction":      float("nan"),
            "std_dominant_fraction":       float("nan"),
            "p90_dominant_fraction":       float("nan"),
            "mean_entropy":                float("nan"),
            "std_entropy":                 float("nan"),
            "p10_entropy":                 float("nan"),
            "p90_entropy":                 float("nan"),
        })
        return base

    # ── main method ───────────────────────────────────────────────────────────

    def finalize(self, thresholds: QCThresholds, eps: float = 1e-12) -> dict:
        """
        Compute ALL metrics first, then apply QC filters.

        Every gene always gets a fully populated row — no empty cells.
        'qc_pass' and 'fail_reason' tell you the outcome; the metrics tell
        you WHY, which is essential for exploratory analysis in notebooks.

        Filter order (cheapest first, checked AFTER computing everything):
            F2 → F3 → F4 → F5
        """
        # ── Concatenate accumulated arrays ────────────────────────────────────
        g_all  = np.concatenate(self.g_values)       if self.g_values       else np.array([], np.float32)
        s_all  = np.concatenate(self.sum_tx_values)  if self.sum_tx_values  else np.array([], np.float32)
        rel    = np.concatenate(self.rel_error_values) if self.rel_error_values else np.array([], np.float32)
        ratio  = np.concatenate(self.ratio_values)   if self.ratio_values   else np.array([], np.float32)
        dom    = np.concatenate(self.dom_frac_values) if self.dom_frac_values else np.array([], np.float32)
        entr   = np.concatenate(self.entropy_values) if self.entropy_values  else np.array([], np.float32)
        n_act  = self.n_active_accumulated

        # ── Compute every metric unconditionally ──────────────────────────────

        # Correlation — nan if not enough informative samples or zero variance
        corr = float("nan")
        if g_all.size >= 3 and np.std(g_all) > 1e-8 and np.std(s_all) > 1e-8:
            corr = float(np.corrcoef(g_all, s_all)[0, 1])

        # Ratio percentiles — nan if no both-detected samples exist
        median_ratio = self._p(ratio, 50)

        # Isoform proportion stats — nan if no active samples
        mtp      = self.sum_tx_props / n_act if n_act > 0 else None
        mean_dom = float(np.mean(dom)) if dom.size else float("nan")

        # Build complete row with a uniform schema.
        # Metrics that cannot be estimated are stored as NaN.
        row = self._base_metrics(g_all, s_all, rel)
        row.update({
            "corr_gene_vs_sum_iso":        corr,
            "median_ratio_sumTx_to_gene":  median_ratio,
            "p10_ratio_sumTx_to_gene":     self._p(ratio, 10),
            "p90_ratio_sumTx_to_gene":     self._p(ratio, 90),
            "n_active_samples":            n_act,
            "max_mean_isoform_proportion": float(np.max(mtp)) if mtp is not None else float("nan"),
            "mean_dominant_fraction":      mean_dom,
            "std_dominant_fraction":       float(np.std(dom))  if dom.size else float("nan"),
            "p90_dominant_fraction":       self._p(dom, 90),
            "mean_entropy":                float(np.mean(entr)) if entr.size else float("nan"),
            "std_entropy":                 float(np.std(entr))  if entr.size else float("nan"),
            "p10_entropy":                 self._p(entr, 10),
            "p90_entropy":                 self._p(entr, 90),
        })

        # ── Apply filters in order, stamp result, return full row ─────────────
        # We never return early anymore — every gene gets a complete row.
        # The fail_reason tells you at which filter the gene was rejected.

        if self.n_both_detected < thresholds.min_both_detected:
            row["qc_pass"]    = False
            row["fail_reason"] = "F2_n_both_detected"
            return row

        if np.isnan(median_ratio) or median_ratio > thresholds.max_median_ratio:
            row["qc_pass"]    = False
            row["fail_reason"] = "F3_median_ratio"
            return row

        if np.isnan(corr) or corr < thresholds.min_corr:
            row["qc_pass"]    = False
            row["fail_reason"] = "F4_corr"
            return row

        if n_act < thresholds.min_active_samples:
            row["qc_pass"]    = False
            row["fail_reason"] = "F5_low_active"
            return row

        if np.isnan(mean_dom) or mean_dom > thresholds.max_dominant_frac:
            row["qc_pass"]    = False
            row["fail_reason"] = "F5_dominant"
            return row

        row["qc_pass"]    = True
        row["fail_reason"] = ""
        return row


# ══════════════════════════════════════════════════════════════════════════════
# Gene record building  (F1)
# ══════════════════════════════════════════════════════════════════════════════

def build_gene_records(
    adata_gene,
    adata_tx,
    gene_to_tx:   Dict[str, List[str]],
    min_isoforms: int = 2,
) -> List[GeneRecord]:
    gene_name_to_idx = {name: i for i, name in enumerate(adata_gene.var_names)}
    tx_name_to_idx   = {name: i for i, name in enumerate(adata_tx.var_names)}
    records, n_missing, n_few = [], 0, 0

    for gene_id, tx_ids in gene_to_tx.items():
        if gene_id not in gene_name_to_idx:
            n_missing += 1
            continue
        valid_tx = [tx for tx in tx_ids if tx in tx_name_to_idx]
        if len(valid_tx) < min_isoforms:
            n_few += 1
            continue
        records.append(GeneRecord(
            gene_id=gene_id,
            gene_idx=gene_name_to_idx[gene_id],
            tx_ids=valid_tx,
            tx_idx=[tx_name_to_idx[tx] for tx in valid_tx],
            n_isoforms=len(valid_tx),
        ))

    print(
        f"F1: {len(records)} genes kept "
        f"({n_missing} not in gene matrix, {n_few} with <{min_isoforms} tx)"
    )
    if not records:
        raise ValueError("No valid genes after F1. Check mapping and min_isoforms.")
    return records


# ══════════════════════════════════════════════════════════════════════════════
# Expression threshold
# ══════════════════════════════════════════════════════════════════════════════

def compute_global_expr_threshold(
    adata_tx,
    gene_records:     List[GeneRecord],
    row_block_size:   int   = 512,
    percentile:       float = 10.0,
    max_genes_sample: int   = 200,
) -> float:
    rng = np.random.default_rng(42)
    subset = gene_records
    if len(gene_records) > max_genes_sample:
        idx    = rng.choice(len(gene_records), size=max_genes_sample, replace=False)
        subset = [gene_records[i] for i in sorted(idx)]

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
        print("Warning: no non-zero sum_tx found — using threshold=0.0")
        return 0.0

    thr = float(np.percentile(np.concatenate(sums), percentile))
    print(f"Expr threshold (p{percentile:.0f} of non-zero sum_tx): {thr:.4f}")
    return thr


# ══════════════════════════════════════════════════════════════════════════════
# Checkpointing
# ══════════════════════════════════════════════════════════════════════════════
# FIX: checkpoint now stores ALL processed genes (pass + fail) with a
# 'qc_pass' column.  Restarted jobs skip every gene already in the file,
# not just passing ones.  This prevents re-processing failed genes.

def _ckpt_path(out_csv: Path) -> Path:
    return out_csv.with_suffix(".checkpoint.csv")


def _load_checkpoint(ckpt: Path):
    """Returns (DataFrame, set_of_done_gene_ids) or (None, empty_set)."""
    if not ckpt.exists():
        return None, set()
    df = pd.read_csv(ckpt)
    done = set(df["gene_id"].tolist())
    n_pass = int(df["qc_pass"].sum()) if "qc_pass" in df.columns else "?"
    print(f"Checkpoint: {len(done)} genes seen ({n_pass} passed), skipping all.")
    return df, done


def _save_checkpoint(df: pd.DataFrame, ckpt: Path) -> None:
    tmp = ckpt.with_suffix(".tmp")
    df.to_csv(tmp, index=False)
    tmp.replace(ckpt)


# ══════════════════════════════════════════════════════════════════════════════
# Main QC computation
# ══════════════════════════════════════════════════════════════════════════════

def compute_gene_qc_table_backed(
    adata_gene,
    adata_tx,
    gene_to_tx:       Dict[str, List[str]],
    thresholds:       Optional[QCThresholds] = None,
    min_isoforms:     int   = 2,
    eps:              float = 1e-12,
    row_block_size:   int   = 512,
    gene_batch_size:  int   = 500,
    max_genes:        Optional[int]   = None,
    expr_threshold:   Optional[float] = None,
    expr_threshold_percentile: float  = 10.0,
    out_csv:          Optional[Path]  = None,
    checkpoint_every: int   = 500,
) -> pd.DataFrame:
    """
    Returns a DataFrame with ALL processed genes (pass + fail).
    Filter on qc_pass == True downstream.
    """
    if thresholds is None:
        thresholds = QCThresholds()

    all_records = build_gene_records(adata_gene, adata_tx, gene_to_tx, min_isoforms)
    if max_genes is not None:
        all_records = all_records[:max_genes]

    # ── Resume from checkpoint ─────────────────────────────────────────────────
    ckpt         = _ckpt_path(out_csv) if out_csv else None
    done_df, done_ids = _load_checkpoint(ckpt) if ckpt else (None, set())
    records      = [r for r in all_records if r.gene_id not in done_ids]
    print(f"Genes to process : {len(records)}  (skipped: {len(done_ids)})")

    if not records:
        print("All genes already processed — returning checkpoint.")
        return done_df if done_df is not None else pd.DataFrame()

    if expr_threshold is None:
        expr_threshold = compute_global_expr_threshold(
            adata_tx, all_records,
            row_block_size=row_block_size,
            percentile=expr_threshold_percentile,
        )

    n_obs     = adata_gene.n_obs
    n_batches = (len(records) + gene_batch_size - 1) // gene_batch_size
    completed: List[dict] = []   # ALL genes, pass + fail
    n_pass, n_fail = 0, 0

    print(f"Gene-batches     : {n_batches} × ≤{gene_batch_size}")
    print(f"Row-block size   : {row_block_size}  |  samples: {n_obs}")
    print(f"Expr threshold   : {expr_threshold:.4f}")
    print(f"QC thresholds    : {thresholds}")

    for bi, bs in enumerate(range(0, len(records), gene_batch_size)):
        batch = records[bs : bs + gene_batch_size]

        gene_cols  = sorted(set(rec.gene_idx for rec in batch))
        gene_local = {g: i for i, g in enumerate(gene_cols)}
        tx_cols    = sorted(set(i for rec in batch for i in rec.tx_idx))
        tx_local   = {t: i for i, t in enumerate(tx_cols)}

        accs = {rec.gene_id: GeneAccumulator(rec.gene_id, rec.n_isoforms) for rec in batch}

        for rs in range(0, n_obs, row_block_size):
            re = min(rs + row_block_size, n_obs)
            G  = _to_dense(adata_gene[rs:re, gene_cols].X)
            T  = _to_dense(adata_tx[rs:re,   tx_cols].X)
            for rec in batch:
                g_expr = G[:, gene_local[rec.gene_idx]]
                tx_mat = T[:, [tx_local[i] for i in rec.tx_idx]]
                accs[rec.gene_id].update(g_expr, tx_mat, eps, expr_threshold)

        # finalize() always returns a dict now — never None
        for rec in batch:
            result = accs[rec.gene_id].finalize(thresholds=thresholds, eps=eps)
            completed.append(result)
            if result["qc_pass"]:
                n_pass += 1
            else:
                n_fail += 1

        # Checkpoint: save ALL completed genes (pass + fail)
        should_ckpt = ckpt and (
            len(completed) % checkpoint_every == 0
            or bi == n_batches - 1
        )
        if should_ckpt:
            new_df   = pd.DataFrame(completed)
            combined = pd.concat([done_df, new_df], ignore_index=True) if done_df is not None else new_df
            _save_checkpoint(combined, ckpt)
            total_seen = len(done_ids) + len(completed)
            print(
                f"  [batch {bi+1}/{n_batches}] ckpt: {total_seen} seen, "
                f"{n_pass} passed, {n_fail} failed"
            )
        else:
            print(f"  [batch {bi+1}/{n_batches}] passed={n_pass} failed={n_fail}")

    new_df = pd.DataFrame(completed)
    final  = pd.concat([done_df, new_df], ignore_index=True) if done_df is not None else new_df

    # Print failure breakdown
    if "fail_reason" in final.columns:
        reasons = final[~final["qc_pass"]]["fail_reason"].value_counts()
        print("\nFailure breakdown:")
        for reason, count in reasons.items():
            print(f"  {reason}: {count}")

    total_pass = int(final["qc_pass"].sum())
    print(f"\nTotal: {total_pass}/{len(final)} genes passed QC")
    return final


# ══════════════════════════════════════════════════════════════════════════════
# Summary helper (honest version)
# ══════════════════════════════════════════════════════════════════════════════

def summarize_qc(qc_df: pd.DataFrame, thresholds: Optional[QCThresholds] = None) -> None:
    """
    Print a summary of QC results from the full output table.

    FIX: previous version re-applied filters only to survivors, giving the
    false impression that filters are soft.  Now we report on the full table
    including failed genes.
    """
    if thresholds is None:
        thresholds = QCThresholds()

    n_total = len(qc_df)
    n_pass  = int(qc_df["qc_pass"].sum())

    print(f"\n{'='*50}")
    print(f"QC SUMMARY  ({n_total} genes total)")
    print(f"{'='*50}")
    print(f"  Passed QC : {n_pass} ({100*n_pass/n_total:.1f}%)")
    print(f"  Failed QC : {n_total - n_pass} ({100*(n_total-n_pass)/n_total:.1f}%)")

    if "fail_reason" in qc_df.columns:
        print("\nFailure breakdown:")
        reasons = qc_df[~qc_df["qc_pass"]]["fail_reason"].value_counts()
        for reason, count in reasons.items():
            pct = 100 * count / n_total
            print(f"  {reason:<30} {count:>6}  ({pct:.1f}% of all genes)")

    print(f"\nThresholds used: {thresholds}")