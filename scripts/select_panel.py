#!/usr/bin/env python3
"""Greedy panel selection for RSVB candidate binding sites.

Behavior:
- compute a composite site score from available columns
- greedily pick candidates that add genomic span while penalizing overlap
- penalize panel-level cross-dimer / cross-hybridization risk using sequence heuristics

Notes on interaction penalty:
- This script uses a lightweight heuristic based on reverse-complement matching.
- For each candidate vs each already-selected oligo, it evaluates:
  * longest contiguous reverse-complement run
  * longest terminal reverse-complement run (last N nt of either oligo)
  * best reverse-complement match fraction at any alignment offset
- These are converted into a normalized pairwise penalty in [0, 1].
- The panel-level penalty for a candidate is the max pairwise penalty against the
  already-selected panel, with an optional weaker sum component.
"""
from __future__ import annotations

import argparse
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Tuple

import pandas as pd


DNA_COMP = str.maketrans({
    "A": "T", "C": "G", "G": "C", "T": "A", "U": "A",
    "R": "Y", "Y": "R", "S": "S", "W": "W", "K": "M", "M": "K",
    "B": "V", "D": "H", "H": "D", "V": "B", "N": "N",
})

IUPAC_MATCH: Dict[str, set[str]] = {
    "A": {"A"}, "C": {"C"}, "G": {"G"}, "T": {"T"}, "U": {"T"},
    "R": {"A", "G"}, "Y": {"C", "T"}, "S": {"G", "C"}, "W": {"A", "T"},
    "K": {"G", "T"}, "M": {"A", "C"}, "B": {"C", "G", "T"},
    "D": {"A", "G", "T"}, "H": {"A", "C", "T"}, "V": {"A", "C", "G"},
    "N": {"A", "C", "G", "T"},
}


def clamp01(x: float) -> float:
    return max(0.0, min(1.0, float(x)))


def clean_seq(seq: str) -> str:
    seq = (seq or "").upper().replace(" ", "").replace("\n", "")
    return "".join(ch if ch in IUPAC_MATCH else "N" for ch in seq)


def revcomp(seq: str) -> str:
    return clean_seq(seq).translate(DNA_COMP)[::-1]


def bases_match(a: str, b: str) -> bool:
    sa = IUPAC_MATCH.get(a.upper(), {"A", "C", "G", "T"})
    sb = IUPAC_MATCH.get(b.upper(), {"A", "C", "G", "T"})
    return not sa.isdisjoint(sb)


def score_gc(gc: float, target: float = 0.48, tolerance: float = 0.20) -> float:
    return clamp01(1.0 - abs(gc - target) / tolerance)


def score_tm(tm: float, target: float = 240.0, tolerance: float = 80.0) -> float:
    # Wallace-like Tm for 80-mers is crude and large in absolute value; placeholder.
    return clamp01(1.0 - abs(tm - target) / tolerance)


def thermo_score(row: pd.Series) -> float:
    gc_part = score_gc(float(row.get("gc", 0.0)))
    tm_part = score_tm(float(row.get("tm_est", 0.0)))
    hp_penalty = min(float(row.get("max_homopolymer", 0)), 8) / 8
    return clamp01(0.45 * gc_part + 0.45 * tm_part + 0.10 * (1.0 - hp_penalty))


def synthesis_score(row: pd.Series) -> float:
    hp_penalty = min(float(row.get("max_homopolymer", 0)), 8) / 8
    lc_penalty = clamp01(float(row.get("low_complexity_fraction", 0.0)))
    return clamp01(1.0 - 0.5 * hp_penalty - 0.5 * lc_penalty)


def site_score(row: pd.Series) -> float:
    return clamp01(
        0.35 * float(row.get("robustness_score", 0.0))
        + 0.25 * float(row.get("accessibility_score", 0.0))
        + 0.15 * thermo_score(row)
        + 0.15 * float(row.get("specificity_score", 0.0))
        + 0.10 * synthesis_score(row)
    )


def interval_overlap(a: Tuple[int, int], b: Tuple[int, int]) -> int:
    start = max(a[0], b[0])
    end = min(a[1], b[1])
    return max(0, end - start + 1)


def added_coverage(candidate: Tuple[int, int], selected: List[Tuple[int, int]]) -> int:
    length = candidate[1] - candidate[0] + 1
    if not selected:
        return length
    covered = 0
    for pos in range(candidate[0], candidate[1] + 1):
        if any(s[0] <= pos <= s[1] for s in selected):
            covered += 1
    return length - covered


def iter_offset_matches(seq_a: str, seq_b: str) -> Iterable[Tuple[int, int, int, int, int]]:
    """Yield alignment summaries for seq_a against seq_b over all offsets.

    Returns tuples:
      (offset, n_overlap, n_match, longest_run, longest_terminal_run)

    Offset is relative to seq_a start against seq_b start. Positive offset means
    seq_b is shifted right relative to seq_a.
    """
    a = clean_seq(seq_a)
    b = clean_seq(seq_b)
    la = len(a)
    lb = len(b)
    for offset in range(-lb + 1, la):
        a_start = max(0, offset)
        b_start = max(0, -offset)
        overlap = min(la - a_start, lb - b_start)
        if overlap <= 0:
            continue

        n_match = 0
        longest_run = 0
        current_run = 0
        longest_terminal_run = 0

        for k in range(overlap):
            ia = a_start + k
            ib = b_start + k
            match = bases_match(a[ia], b[ib])
            if match:
                n_match += 1
                current_run += 1
                if current_run > longest_run:
                    longest_run = current_run
            else:
                current_run = 0

        # second pass for terminal runs so we can require end-contact on either oligo
        current_run = 0
        run_starts_a = 0
        run_starts_b = 0
        for k in range(overlap):
            ia = a_start + k
            ib = b_start + k
            match = bases_match(a[ia], b[ib])
            if match:
                if current_run == 0:
                    run_starts_a = ia
                    run_starts_b = ib
                current_run += 1
            else:
                if current_run > 0:
                    run_end_a = ia - 1
                    run_end_b = ib - 1
                    if run_starts_a == 0 or run_starts_b == 0 or run_end_a == la - 1 or run_end_b == lb - 1:
                        longest_terminal_run = max(longest_terminal_run, current_run)
                current_run = 0
        if current_run > 0:
            run_end_a = a_start + overlap - 1
            run_end_b = b_start + overlap - 1
            if run_starts_a == 0 or run_starts_b == 0 or run_end_a == la - 1 or run_end_b == lb - 1:
                longest_terminal_run = max(longest_terminal_run, current_run)

        yield offset, overlap, n_match, longest_run, longest_terminal_run


def pair_interaction_metrics(seq1: str, seq2: str) -> Dict[str, float]:
    """Estimate pairwise cross-dimer / cross-hybridization risk.

    We compare seq1 to reverse-complement(seq2), because hybridization occurs by
    reverse-complementary pairing.
    """
    s1 = clean_seq(seq1)
    s2_rc = revcomp(seq2)
    best_fraction = 0.0
    best_run = 0
    best_terminal_run = 0
    best_match_count = 0
    best_overlap = 0
    best_offset = 0

    for offset, overlap, n_match, longest_run, longest_terminal_run in iter_offset_matches(s1, s2_rc):
        frac = n_match / overlap if overlap else 0.0
        # rank by match fraction, then by run lengths.
        if (
            frac > best_fraction
            or (frac == best_fraction and longest_run > best_run)
            or (frac == best_fraction and longest_run == best_run and longest_terminal_run > best_terminal_run)
        ):
            best_fraction = frac
            best_run = longest_run
            best_terminal_run = longest_terminal_run
            best_match_count = n_match
            best_overlap = overlap
            best_offset = offset

    return {
        "best_match_fraction": float(best_fraction),
        "best_match_count": int(best_match_count),
        "best_overlap": int(best_overlap),
        "max_rc_run": int(best_run),
        "max_terminal_rc_run": int(best_terminal_run),
        "best_offset": int(best_offset),
    }


def normalize_run(run_len: float, safe: int, risky: int) -> float:
    if risky <= safe:
        return 1.0 if run_len >= risky else 0.0
    return clamp01((run_len - safe) / (risky - safe))


def pair_penalty_from_metrics(
    metrics: Dict[str, float],
    *,
    safe_run: int,
    risky_run: int,
    safe_terminal_run: int,
    risky_terminal_run: int,
    safe_fraction: float,
    risky_fraction: float,
) -> float:
    run_pen = normalize_run(metrics["max_rc_run"], safe_run, risky_run)
    terminal_pen = normalize_run(metrics["max_terminal_rc_run"], safe_terminal_run, risky_terminal_run)
    if risky_fraction <= safe_fraction:
        frac_pen = 1.0 if metrics["best_match_fraction"] >= risky_fraction else 0.0
    else:
        frac_pen = clamp01((metrics["best_match_fraction"] - safe_fraction) / (risky_fraction - safe_fraction))
    # Heavier weight on long contiguous runs and terminal runs; density is secondary.
    return clamp01(0.45 * run_pen + 0.35 * terminal_pen + 0.20 * frac_pen)


def candidate_panel_interaction(
    candidate_seq: str,
    selected_seqs: Sequence[str],
    *,
    safe_run: int,
    risky_run: int,
    safe_terminal_run: int,
    risky_terminal_run: int,
    safe_fraction: float,
    risky_fraction: float,
) -> Dict[str, float]:
    if not selected_seqs:
        return {
            "cross_penalty_max": 0.0,
            "cross_penalty_sum": 0.0,
            "max_rc_run": 0,
            "max_terminal_rc_run": 0,
            "max_rc_fraction": 0.0,
            "n_risky_pairs": 0,
        }

    penalties: List[float] = []
    max_run = 0
    max_terminal_run = 0
    max_fraction = 0.0
    n_risky_pairs = 0

    for seq in selected_seqs:
        metrics = pair_interaction_metrics(candidate_seq, seq)
        penalty = pair_penalty_from_metrics(
            metrics,
            safe_run=safe_run,
            risky_run=risky_run,
            safe_terminal_run=safe_terminal_run,
            risky_terminal_run=risky_terminal_run,
            safe_fraction=safe_fraction,
            risky_fraction=risky_fraction,
        )
        penalties.append(penalty)
        if penalty >= 0.5:
            n_risky_pairs += 1
        max_run = max(max_run, int(metrics["max_rc_run"]))
        max_terminal_run = max(max_terminal_run, int(metrics["max_terminal_rc_run"]))
        max_fraction = max(max_fraction, float(metrics["best_match_fraction"]))

    return {
        "cross_penalty_max": max(penalties),
        "cross_penalty_sum": sum(penalties),
        "max_rc_run": max_run,
        "max_terminal_rc_run": max_terminal_run,
        "max_rc_fraction": max_fraction,
        "n_risky_pairs": n_risky_pairs,
    }


def main() -> None:
    ap = argparse.ArgumentParser(description="Select RSVB binding-site panel")
    ap.add_argument("--scored-tsv", required=True, help="Candidate TSV after conservation/accessibility/specificity scoring")
    ap.add_argument("--ranked-output", required=True)
    ap.add_argument("--panel-output", required=True)
    ap.add_argument("--max-sites", type=int, default=12)
    ap.add_argument("--min-site-score", type=float, default=0.45)
    ap.add_argument("--min-gap", type=int, default=500, help="Preferred minimal spacing between selected sites")
    ap.add_argument("--hard-min-gap", type=int, default=120, help="Hard filter: skip candidates whose nearest gap to the current panel is below this many bases")
    ap.add_argument("--max-overlap-fraction", type=float, default=0.10, help="Hard filter: skip candidates overlapping any selected site by more than this fraction of candidate length")
    ap.add_argument("--min-new-bp", type=int, default=65, help="Hard filter: candidate must contribute at least this many newly covered bases")
    ap.add_argument("--cross-penalty-weight", type=float, default=0.15, help="Weight for pairwise cross-dimer/hybridization risk in greedy gain")
    ap.add_argument("--cross-penalty-sum-weight", type=float, default=0.05, help="Additional weaker weight for cumulative interaction burden across the panel")
    ap.add_argument("--safe-run", type=int, default=8, help="Longest reverse-complement run considered mostly safe")
    ap.add_argument("--risky-run", type=int, default=12, help="Longest reverse-complement run considered high-risk")
    ap.add_argument("--safe-terminal-run", type=int, default=5, help="Terminal reverse-complement run considered mostly safe")
    ap.add_argument("--risky-terminal-run", type=int, default=8, help="Terminal reverse-complement run considered high-risk")
    ap.add_argument("--safe-match-fraction", type=float, default=0.55, help="Best offset reverse-complement fraction considered mostly safe")
    ap.add_argument("--risky-match-fraction", type=float, default=0.75, help="Best offset reverse-complement fraction considered high-risk")
    ap.add_argument("--hard-max-terminal-run", type=int, default=10, help="Optional hard filter: skip candidates whose worst terminal reverse-complement run exceeds this against the current panel")
    ap.add_argument("--hard-max-run", type=int, default=14, help="Optional hard filter: skip candidates whose worst reverse-complement run exceeds this against the current panel")
    args = ap.parse_args()

    df = pd.read_csv(args.scored_tsv, sep="\t")
    if "candidate_keep" in df.columns:
        df = df[df["candidate_keep"].astype(bool)].copy()

    if "sequence_ref" not in df.columns and "seq" in df.columns:
        df["sequence_ref"] = df["seq"]
    if "sequence_ref" not in df.columns:
        raise ValueError("Input TSV must include 'sequence_ref' or 'seq' column for interaction scoring.")

    df["sequence_ref"] = df["sequence_ref"].astype(str).map(clean_seq)
    df["thermo_score"] = df.apply(thermo_score, axis=1)
    df["synthesis_score"] = df.apply(synthesis_score, axis=1)
    df["site_score"] = df.apply(site_score, axis=1)
    ranked = df.sort_values(["site_score", "robustness_score", "accessibility_score"], ascending=False).reset_index(drop=True)

    selected_rows = []
    selected_intervals: List[Tuple[int, int]] = []
    selected_seqs: List[str] = []

    for _, row in ranked.iterrows():
        if len(selected_rows) >= args.max_sites:
            break
        if float(row["site_score"]) < args.min_site_score:
            continue

        cand = (int(row["start"]), int(row["end"]))
        cand_seq = str(row["sequence_ref"])
        overlap_bp = sum(interval_overlap(cand, s) for s in selected_intervals)
        gain_bp = added_coverage(cand, selected_intervals)
        max_overlap_fraction = 0.0
        if selected_intervals:
            max_overlap_fraction = max(
                interval_overlap(cand, s) / max(int(row["length"]), 1)
                for s in selected_intervals
            )
        if max_overlap_fraction > args.max_overlap_fraction:
            continue
        if gain_bp < args.min_new_bp:
            continue

        min_distance_penalty = 0.0
        nearest = None
        if selected_intervals:
            nearest = min(min(abs(cand[0] - s[1]), abs(cand[1] - s[0])) for s in selected_intervals)
            if nearest < args.hard_min_gap:
                continue
            if nearest < args.min_gap:
                min_distance_penalty = (args.min_gap - nearest) / args.min_gap

        interaction = candidate_panel_interaction(
            cand_seq,
            selected_seqs,
            safe_run=args.safe_run,
            risky_run=args.risky_run,
            safe_terminal_run=args.safe_terminal_run,
            risky_terminal_run=args.risky_terminal_run,
            safe_fraction=args.safe_match_fraction,
            risky_fraction=args.risky_match_fraction,
        )

        if interaction["max_terminal_rc_run"] > args.hard_max_terminal_run:
            continue
        if interaction["max_rc_run"] > args.hard_max_run:
            continue

        cross_penalty = (
            args.cross_penalty_weight * float(interaction["cross_penalty_max"])
            + args.cross_penalty_sum_weight * float(interaction["cross_penalty_sum"])
        )

        greedy_gain = (
            0.60 * gain_bp / max(int(row["length"]), 1)
            + 0.22 * float(row["site_score"])
            - 0.10 * overlap_bp / max(int(row["length"]), 1)
            - 0.08 * min_distance_penalty
            - cross_penalty
        )
        if greedy_gain <= 0:
            continue

        selected_intervals.append(cand)
        selected_seqs.append(cand_seq)
        rr = row.to_dict()
        rr.update(interaction)
        rr["overlap_bp"] = int(overlap_bp)
        rr["gain_bp"] = int(gain_bp)
        rr["max_overlap_fraction"] = round(max_overlap_fraction, 6)
        rr["nearest_gap_bp"] = int(nearest) if nearest is not None else int(row["length"])
        rr["cross_penalty"] = round(cross_penalty, 6)
        rr["greedy_gain"] = round(greedy_gain, 6)
        rr["reason_selected"] = (
            f"site_score={row['site_score']:.3f}; gain_bp={gain_bp}; overlap_bp={overlap_bp}; max_overlap_fraction={max_overlap_fraction:.3f}; "
            f"nearest_gap_bp={nearest if nearest is not None else int(row['length'])}; spacing_penalty={min_distance_penalty:.3f}; cross_penalty={cross_penalty:.3f}; "
            f"max_rc_run={interaction['max_rc_run']}; max_terminal_rc_run={interaction['max_terminal_rc_run']}; "
            f"max_rc_fraction={interaction['max_rc_fraction']:.3f}"
        )
        selected_rows.append(rr)

    Path(args.ranked_output).parent.mkdir(parents=True, exist_ok=True)
    Path(args.panel_output).parent.mkdir(parents=True, exist_ok=True)

    ranked.to_csv(args.ranked_output, sep="\t", index=False)
    panel_df = pd.DataFrame(selected_rows)
    if not panel_df.empty:
        panel_df.insert(0, "panel_rank", range(1, len(panel_df) + 1))
    panel_df.to_csv(args.panel_output, sep="\t", index=False)

    print(f"Ranked candidates: {len(ranked)}")
    print(f"Selected panel size: {len(panel_df)}")
    print(f"Wrote ranked TSV -> {args.ranked_output}")
    print(f"Wrote panel TSV -> {args.panel_output}")


if __name__ == "__main__":
    main()
