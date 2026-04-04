#!/usr/bin/env python3
"""Attach off-target / specificity scores to candidate sites.

Preferred mode:
- parse external search results (BLAST/MFEprimer-preprocessed summaries)
  from TSV with columns: site_id, human_offtarget_hits, virus_bg_hits

Fallback mode:
- initialize specificity fields but mark them as not yet externally evaluated.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


def specificity_from_hits(human_hits: int, virus_hits: int) -> float:
    penalty = min(1.0, 0.7 * min(human_hits, 5) / 5 + 0.3 * min(virus_hits, 5) / 5)
    return max(0.0, 1.0 - penalty)


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description="Score RSVB candidates for specificity")
    ap.add_argument("--candidates-tsv", required=True)
    ap.add_argument("--output-tsv", required=True)
    ap.add_argument(
        "--hits-tsv",
        default=None,
        help="Optional TSV with columns: site_id, human_offtarget_hits, virus_bg_hits",
    )
    return ap.parse_args()


def main() -> None:
    args = parse_args()
    df = pd.read_csv(args.candidates_tsv, sep="\t")

    if args.hits_tsv:
        hits = pd.read_csv(args.hits_tsv, sep="\t")
        required = {"site_id", "human_offtarget_hits", "virus_bg_hits"}
        if not required.issubset(hits.columns):
            raise SystemExit(f"{args.hits_tsv} must contain columns: {sorted(required)}")
        df = df.merge(hits[list(required)], on="site_id", how="left")
        df["human_offtarget_hits"] = df["human_offtarget_hits"].fillna(0).astype(int)
        df["virus_bg_hits"] = df["virus_bg_hits"].fillna(0).astype(int)
        df["specificity_mode"] = "external_hits"
    else:
        df["human_offtarget_hits"] = 0
        df["virus_bg_hits"] = 0
        df["specificity_mode"] = "placeholder_no_external_hits"

    df["offtarget_penalty"] = (
        0.7 * df["human_offtarget_hits"].clip(upper=5) / 5
        + 0.3 * df["virus_bg_hits"].clip(upper=5) / 5
    ).clip(0, 1)
    df["specificity_score"] = [
        specificity_from_hits(int(h), int(v))
        for h, v in zip(df["human_offtarget_hits"], df["virus_bg_hits"])
    ]

    Path(args.output_tsv).parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.output_tsv, sep="\t", index=False)
    print(f"Added specificity scores for {len(df)} candidates")
    print(f"Wrote -> {args.output_tsv}")


if __name__ == "__main__":
    main()
