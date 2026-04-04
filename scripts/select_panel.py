#!/usr/bin/env python3
"""Greedy panel selection for RSVB candidate binding sites.

Initial MVP behavior:
- compute a composite site score from available columns
- greedily pick candidates that add genomic span while penalizing overlap
- output ranked candidates and final selected panel
"""
from __future__ import annotations

import argparse
from pathlib import Path
from typing import List, Tuple

import pandas as pd


def clamp01(x: float) -> float:
    return max(0.0, min(1.0, float(x)))


def score_gc(gc: float, target: float = 0.48, tolerance: float = 0.20) -> float:
    return clamp01(1.0 - abs(gc - target) / tolerance)


def score_tm(tm: float, target: float = 240.0, tolerance: float = 80.0) -> float:
    # Wallace Tm for 80-mers is crude and large in absolute value; this is just a placeholder.
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


def main() -> None:
    ap = argparse.ArgumentParser(description="Select RSVB binding-site panel")
    ap.add_argument("--scored-tsv", required=True, help="Candidate TSV after conservation/accessibility/specificity scoring")
    ap.add_argument("--ranked-output", required=True)
    ap.add_argument("--panel-output", required=True)
    ap.add_argument("--max-sites", type=int, default=12)
    ap.add_argument("--min-site-score", type=float, default=0.45)
    ap.add_argument("--min-gap", type=int, default=500, help="Preferred minimal spacing between selected sites")
    args = ap.parse_args()

    df = pd.read_csv(args.scored_tsv, sep="\t")
    if "candidate_keep" in df.columns:
        df = df[df["candidate_keep"].astype(bool)].copy()

    df["thermo_score"] = df.apply(thermo_score, axis=1)
    df["synthesis_score"] = df.apply(synthesis_score, axis=1)
    df["site_score"] = df.apply(site_score, axis=1)
    ranked = df.sort_values(["site_score", "robustness_score", "accessibility_score"], ascending=False).reset_index(drop=True)

    selected_rows = []
    selected_intervals: List[Tuple[int, int]] = []

    for _, row in ranked.iterrows():
        if len(selected_rows) >= args.max_sites:
            break
        if float(row["site_score"]) < args.min_site_score:
            continue

        cand = (int(row["start"]), int(row["end"]))
        overlap_bp = sum(interval_overlap(cand, s) for s in selected_intervals)
        gain_bp = added_coverage(cand, selected_intervals)
        min_distance_penalty = 0.0
        if selected_intervals:
            nearest = min(
                min(abs(cand[0] - s[1]), abs(cand[1] - s[0]))
                for s in selected_intervals
            )
            if nearest < args.min_gap:
                min_distance_penalty = (args.min_gap - nearest) / args.min_gap

        greedy_gain = 0.45 * gain_bp / max(int(row["length"]), 1) + 0.30 * float(row["site_score"]) - 0.15 * overlap_bp / max(int(row["length"]), 1) - 0.10 * min_distance_penalty
        if greedy_gain <= 0:
            continue

        selected_intervals.append(cand)
        rr = row.to_dict()
        rr["greedy_gain"] = round(greedy_gain, 6)
        rr["reason_selected"] = f"site_score={row['site_score']:.3f}; gain_bp={gain_bp}; overlap_bp={overlap_bp}; spacing_penalty={min_distance_penalty:.3f}"
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
