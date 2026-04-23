#!/usr/bin/env python3
"""Score RT-primer candidates for self-complementarity and hairpin risk.

Stage 4.5 — runs after score_conservation.py, before score_accessibility.py
(or after accessibility if you prefer; it only needs the candidates TSV with
rt_primer_seq present).

Background
----------
A 25-nt RT primer can fold back on itself in two ways:

  1. **Hairpin** — the 3' end base-pairs with an internal region of the same
     strand, forming a stem-loop.  A stem of ≥ 4 bp lowers effective primer
     concentration and, if 3'-anchored, directly blocks polymerase extension.

  2. **Self-dimer** — two copies of the same primer hybridise to each other
     in antiparallel orientation.  This competes with target binding and can
     be extended by the polymerase into primer-dimer artefacts.

Both are detected by sliding the primer sequence against its own reverse
complement (RC) at every offset and finding the longest contiguous
IUPAC-aware match run.  The algorithm reuses the same approach already
present in select_panel.py for pairwise primer interaction.

Output columns added
--------------------
self_rc_max_run      : int    longest match run at any offset
self_rc_3p_run       : int    longest run involving the 3' terminal base
self_rc_max_fraction : float  highest match fraction at any offset
self_comp_penalty    : float  composite [0,1] penalty (higher = riskier)
self_comp_flag       : str    "ok" | "warn" | "fail"  — human-readable tier

Flag thresholds (configurable via CLI)
--------------------------------------
  --hard-max-self-run    : 3'-run length above which candidate is flagged FAIL
  --warn-self-run        : 3'-run length above which candidate is flagged WARN
  --hard-max-body-run    : whole-body run length above which → FAIL

The script does NOT drop candidates — it annotates them.  Hard exclusion is
left to select_panel.py (or the user inspecting the ranked TSV), so you can
always see what was flagged rather than having candidates silently disappear.
"""
from __future__ import annotations

import argparse
from pathlib import Path
from typing import Dict, Iterable, Tuple

import pandas as pd


# ---------------------------------------------------------------------------
# Sequence utilities (self-contained so this script has no imports from
# other pipeline scripts — avoids fragile relative-import issues)
# ---------------------------------------------------------------------------

DNA_COMP = str.maketrans({
    "A": "T", "C": "G", "G": "C", "T": "A", "U": "A",
    "R": "Y", "Y": "R", "S": "S", "W": "W", "K": "M", "M": "K",
    "B": "V", "D": "H", "H": "D", "V": "B", "N": "N",
})

IUPAC_MASK: Dict[str, int] = {
    "A": 0b1000, "C": 0b0100, "G": 0b0010, "T": 0b0001,
    "U": 0b0001,
    "R": 0b1010, "Y": 0b0101, "S": 0b0110, "W": 0b1001,
    "K": 0b0011, "M": 0b1100,
    "B": 0b0111, "D": 0b1011, "H": 0b1101, "V": 0b1110,
    "N": 0b1111, "-": 0b0000,
}


def clean_seq(seq: str) -> str:
    return "".join(
        ch if ch in IUPAC_MASK else "N"
        for ch in seq.upper().replace("U", "T")
    )


def revcomp(seq: str) -> str:
    return clean_seq(seq).translate(DNA_COMP)[::-1]


def bases_match(a: str, b: str) -> bool:
    """IUPAC-aware single-base compatibility check."""
    return bool(IUPAC_MASK.get(a, 0) & IUPAC_MASK.get(b, 0))


def clamp01(x: float) -> float:
    return max(0.0, min(1.0, float(x)))


# ---------------------------------------------------------------------------
# Core algorithm
# ---------------------------------------------------------------------------

def _iter_offsets(seq_a: str, seq_b: str) -> Iterable[Tuple[int, int, int, int, int]]:
    """Slide seq_b over seq_a at every offset; yield per-offset match stats.

    Yields
    ------
    offset, overlap_len, n_match, longest_run, longest_terminal_run
        where "terminal" means the run touches position 0 or -1 of either
        sequence (i.e. it is anchored at an end).
    """
    la, lb = len(seq_a), len(seq_b)
    for offset in range(-lb + 1, la):
        a_start = max(0, offset)
        b_start = max(0, -offset)
        overlap = min(la - a_start, lb - b_start)
        if overlap <= 0:
            continue

        # --- pass 1: count matches and find longest run ---
        n_match = 0
        cur_run = 0
        longest_run = 0
        for k in range(overlap):
            if bases_match(seq_a[a_start + k], seq_b[b_start + k]):
                n_match += 1
                cur_run += 1
                if cur_run > longest_run:
                    longest_run = cur_run
            else:
                cur_run = 0

        # --- pass 2: find longest terminal run ---
        # A run is "terminal" when it starts or ends at position 0 or
        # len-1 of either participant (covers both hairpin ends).
        cur_run = 0
        run_start_a = 0
        longest_terminal_run = 0
        for k in range(overlap):
            ia, ib = a_start + k, b_start + k
            if bases_match(seq_a[ia], seq_b[ib]):
                if cur_run == 0:
                    run_start_a = ia
                cur_run += 1
            else:
                if cur_run > 0:
                    run_end_a = ia - 1
                    if run_start_a == 0 or ib - cur_run == 0 \
                            or run_end_a == la - 1 \
                            or (b_start + k - 1) == lb - 1:
                        if cur_run > longest_terminal_run:
                            longest_terminal_run = cur_run
                cur_run = 0
        if cur_run > 0:
            run_end_a = a_start + overlap - 1
            if run_start_a == 0 or (b_start + overlap - cur_run) == 0 \
                    or run_end_a == la - 1 \
                    or (b_start + overlap - 1) == lb - 1:
                if cur_run > longest_terminal_run:
                    longest_terminal_run = cur_run

        yield offset, overlap, n_match, longest_run, longest_terminal_run


def self_comp_metrics(primer_seq: str) -> Dict[str, float]:
    """Compute self-complementarity metrics for a single primer.

    Folds the primer against its own reverse complement at every
    register offset; identifies the worst-case hairpin/self-dimer
    signal.

    Parameters
    ----------
    primer_seq : str
        The primer oligonucleotide sequence (5'→3'), typically
        rt_primer_seq (the actual oligo to be ordered).

    Returns
    -------
    dict with keys:
        self_rc_max_run      (int)   — longest match run at any offset
        self_rc_3p_run       (int)   — longest run touching the 3' end
        self_rc_max_fraction (float) — highest match/overlap fraction
    """
    seq = clean_seq(primer_seq)
    seq_rc = revcomp(primer_seq)

    best_run = 0
    best_3p_run = 0
    best_frac = 0.0
    n = len(seq)

    for offset, overlap, n_match, longest_run, longest_terminal_run in _iter_offsets(seq, seq_rc):
        frac = n_match / overlap if overlap else 0.0
        if longest_run > best_run:
            best_run = longest_run
        if frac > best_frac:
            best_frac = frac
        # "terminal" in the context of the original primer means the run
        # touches position n-1 of seq (the 3' end of the primer)
        a_start = max(0, offset)
        run_end_a = a_start + overlap - 1
        if run_end_a == n - 1 and longest_terminal_run > best_3p_run:
            best_3p_run = longest_terminal_run

    return {
        "self_rc_max_run":      best_run,
        "self_rc_3p_run":       best_3p_run,
        "self_rc_max_fraction": round(best_frac, 6),
    }


def _normalize_run(run_len: float, safe: int, risky: int) -> float:
    if risky <= safe:
        return 1.0 if run_len >= risky else 0.0
    return clamp01((run_len - safe) / (risky - safe))


def compute_penalty(metrics: Dict[str, float],
                    safe_body_run: int,
                    risky_body_run: int,
                    safe_3p_run: int,
                    risky_3p_run: int) -> float:
    """Map self-comp metrics to a [0, 1] penalty.

    The 3'-anchored run is weighted 60% because it directly competes
    with target binding at the initiation site of reverse transcription.
    The whole-body run is weighted 40% to catch severe internal hairpins
    that reduce effective primer availability.
    """
    body_pen = _normalize_run(metrics["self_rc_max_run"], safe_body_run, risky_body_run)
    term_pen = _normalize_run(metrics["self_rc_3p_run"],  safe_3p_run,   risky_3p_run)
    return round(clamp01(0.40 * body_pen + 0.60 * term_pen), 6)


def assign_flag(metrics: Dict[str, float],
                hard_max_body_run: int,
                hard_max_3p_run: int,
                warn_3p_run: int) -> str:
    """Return 'ok', 'warn', or 'fail' for human-readable tier output."""
    if (metrics["self_rc_max_run"] > hard_max_body_run
            or metrics["self_rc_3p_run"] > hard_max_3p_run):
        return "fail"
    if metrics["self_rc_3p_run"] >= warn_3p_run:
        return "warn"
    return "ok"


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(
        description="Score RT-primer candidates for self-complementarity and hairpin risk.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    ap.add_argument("--candidates-tsv", required=True,
                    help="Input candidates TSV (must contain rt_primer_seq or sequence_ref).")
    ap.add_argument("--output-tsv", required=True,
                    help="Output TSV with self-comp columns appended.")
    ap.add_argument("--assay-type", choices=["capture_long", "rt_primer_25"],
                    default="rt_primer_25")

    thr = ap.add_argument_group("flagging thresholds")
    thr.add_argument("--hard-max-body-run", type=int, default=7,
                     help="Whole-body self-RC run length → FAIL flag. "
                          "A 7-bp internal stem is generally problematic for a 25-mer.")
    thr.add_argument("--hard-max-3p-run", type=int, default=4,
                     help="3'-anchored self-RC run length → FAIL flag. "
                          "≥ 4 bp at the 3' end blocks polymerase initiation.")
    thr.add_argument("--warn-3p-run", type=int, default=3,
                     help="3'-anchored run length → WARN flag (below hard-max-3p-run).")
    thr.add_argument("--safe-body-run", type=int, default=4,
                     help="Run length below which body penalty is zero.")
    thr.add_argument("--risky-body-run", type=int, default=7,
                     help="Run length at/above which body penalty is maximum.")
    thr.add_argument("--safe-3p-run", type=int, default=3,
                     help="3'-run length below which 3' penalty is zero.")
    thr.add_argument("--risky-3p-run", type=int, default=5,
                     help="3'-run length at/above which 3' penalty is maximum.")
    return ap.parse_args()


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    args = parse_args()

    df = pd.read_csv(args.candidates_tsv, sep="\t")

    # Resolve which sequence column to use as the primer oligo.
    # For rt_primer_25 the oligo ordered is rt_primer_seq (revcomp of target).
    # For capture_long the oligo is sequence_ref directly.
    if "rt_primer_seq" in df.columns and args.assay_type == "rt_primer_25":
        primer_col = "rt_primer_seq"
    elif "sequence_ref" in df.columns:
        primer_col = "sequence_ref"
        if args.assay_type == "rt_primer_25":
            print(
                "WARNING: rt_primer_seq not found in input TSV — "
                "falling back to sequence_ref. "
                "Self-comp is being assessed on the TARGET sequence, not the primer oligo."
            )
    else:
        raise SystemExit(
            "Input TSV must contain either 'rt_primer_seq' or 'sequence_ref'."
        )

    print(f"Scoring {len(df)} candidates for self-complementarity "
          f"(primer column: {primer_col!r}) …")

    # Compute per-row metrics
    metrics_list = [
        self_comp_metrics(str(seq))
        for seq in df[primer_col]
    ]

    df["self_rc_max_run"]      = [m["self_rc_max_run"]      for m in metrics_list]
    df["self_rc_3p_run"]       = [m["self_rc_3p_run"]       for m in metrics_list]
    df["self_rc_max_fraction"] = [m["self_rc_max_fraction"] for m in metrics_list]

    df["self_comp_penalty"] = [
        compute_penalty(
            m,
            safe_body_run  = args.safe_body_run,
            risky_body_run = args.risky_body_run,
            safe_3p_run    = args.safe_3p_run,
            risky_3p_run   = args.risky_3p_run,
        )
        for m in metrics_list
    ]

    df["self_comp_flag"] = [
        assign_flag(
            m,
            hard_max_body_run = args.hard_max_body_run,
            hard_max_3p_run   = args.hard_max_3p_run,
            warn_3p_run       = args.warn_3p_run,
        )
        for m in metrics_list
    ]

    Path(args.output_tsv).parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.output_tsv, sep="\t", index=False)

    # Summary
    n_fail = int((df["self_comp_flag"] == "fail").sum())
    n_warn = int((df["self_comp_flag"] == "warn").sum())
    n_ok   = int((df["self_comp_flag"] == "ok").sum())
    print(f"  ok  : {n_ok:5d}")
    print(f"  warn: {n_warn:5d}  (3'-run ≥ {args.warn_3p_run})")
    print(f"  fail: {n_fail:5d}  (3'-run > {args.hard_max_3p_run} "
          f"or body-run > {args.hard_max_body_run})")
    print(f"Wrote → {args.output_tsv}")


if __name__ == "__main__":
    main()
