#!/usr/bin/env python3
"""Attach accessibility scores to candidate sites.

Preferred mode:
- provide --unpaired-prob-tsv from RNAplfold-derived preprocessing
  with columns: position, p_unpaired

Fallback mode:
- if no external profile is supplied, use a simple heuristic based on GC.
  This is only a scaffold and should be replaced by real RNA structure output.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


def heuristic_accessibility(gc: float) -> float:
    # Very rough placeholder: moderate GC tends to be favored.
    target = 0.48
    score = 1.0 - abs(gc - target) / target
    return max(0.0, min(1.0, score))


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description="Score RSVB candidates for accessibility")
    ap.add_argument("--candidates-tsv", required=True)
    ap.add_argument("--output-tsv", required=True)
    ap.add_argument(
        "--unpaired-prob-tsv",
        help="Optional TSV with columns: position, p_unpaired",
        default=None,
    )
    return ap.parse_args()


def main() -> None:
    args = parse_args()
    df = pd.read_csv(args.candidates_tsv, sep="\t")

    if args.unpaired_prob_tsv:
        up = pd.read_csv(args.unpaired_prob_tsv, sep="\t")
        required = {"position", "p_unpaired"}
        if not required.issubset(up.columns):
            raise SystemExit(f"{args.unpaired_prob_tsv} must contain columns: {sorted(required)}")
        up_map = dict(zip(up["position"].astype(int), up["p_unpaired"].astype(float)))

        access_mean = []
        access_min = []
        for _, row in df.iterrows():
            positions = range(int(row["start"]), int(row["end"]) + 1)
            vals = [up_map.get(p, 0.0) for p in positions]
            access_mean.append(sum(vals) / len(vals) if vals else 0.0)
            access_min.append(min(vals) if vals else 0.0)

        df["access_mean"] = access_mean
        df["access_min"] = access_min
        df["accessibility_mode"] = "external_profile"
    else:
        df["access_mean"] = df["gc"].astype(float).map(heuristic_accessibility)
        df["access_min"] = (df["access_mean"] * 0.85).round(6)
        df["accessibility_mode"] = "heuristic_gc_placeholder"

    df["accessibility_score"] = (0.7 * df["access_mean"] + 0.3 * df["access_min"]).clip(0, 1)

    Path(args.output_tsv).parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.output_tsv, sep="\t", index=False)
    print(f"Added accessibility scores for {len(df)} candidates")
    print(f"Mode: {df['accessibility_mode'].iloc[0] if len(df) else 'n/a'}")
    print(f"Wrote -> {args.output_tsv}")


if __name__ == "__main__":
    main()
