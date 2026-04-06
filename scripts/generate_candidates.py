#!/usr/bin/env python3
"""Generate sliding-window candidate binding sites from an RSVB reference genome.

Initial MVP behavior:
- read one reference FASTA sequence
- generate windows with fixed size and step
- compute simple sequence-level filters/features
- output candidate TSV
"""
from __future__ import annotations

import argparse
from pathlib import Path
from typing import List

import pandas as pd
from Bio import SeqIO


def gc_fraction(seq: str) -> float:
    seq = seq.upper()
    gc = seq.count("G") + seq.count("C")
    return gc / len(seq) if seq else 0.0


def max_homopolymer(seq: str) -> int:
    if not seq:
        return 0
    max_run = 1
    current = 1
    upper = seq.upper()
    for i in range(1, len(upper)):
        if upper[i] == upper[i - 1]:
            current += 1
            max_run = max(max_run, current)
        else:
            current = 1
    return max_run


def low_complexity_fraction(seq: str, k: int = 3) -> float:
    seq = seq.upper()
    if len(seq) < k:
        return 0.0
    kmers = [seq[i : i + k] for i in range(len(seq) - k + 1)]
    unique = len(set(kmers))
    return 1.0 - (unique / len(kmers))


def wallace_tm(seq: str) -> float:
    seq = seq.upper()
    a = seq.count("A")
    t = seq.count("T")
    g = seq.count("G")
    c = seq.count("C")
    return 2 * (a + t) + 4 * (g + c)


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description="Generate candidate RSVB binding sites")
    ap.add_argument("--ref-fasta", required=True)
    ap.add_argument("--output-tsv", required=True)
    ap.add_argument("--window", type=int, default=80)
    ap.add_argument("--step", type=int, default=10)
    ap.add_argument("--min-gc", type=float, default=0.22)
    ap.add_argument("--max-gc", type=float, default=0.68)
    ap.add_argument("--max-homopolymer", type=int, default=7)
    ap.add_argument("--max-low-complexity", type=float, default=0.88,
                    help="Upper bound on low-complexity fraction. RSVB defaults are intentionally relaxed and k=3-based.")
    return ap.parse_args()


def main() -> None:
    args = parse_args()
    ref_records = list(SeqIO.parse(args.ref_fasta, "fasta"))
    if len(ref_records) != 1:
        raise SystemExit("--ref-fasta must contain exactly one reference sequence")

    ref = ref_records[0]
    seq = str(ref.seq).upper()
    rows: List[dict] = []

    site_num = 0
    for start0 in range(0, len(seq) - args.window + 1, args.step):
        end0 = start0 + args.window
        window_seq = seq[start0:end0]
        gc = gc_fraction(window_seq)
        max_hp = max_homopolymer(window_seq)
        lc = low_complexity_fraction(window_seq, k=3)
        tm_est = wallace_tm(window_seq)

        keep = True
        fail_reasons = []
        if gc < args.min_gc:
            keep = False
            fail_reasons.append("low_gc")
        if gc > args.max_gc:
            keep = False
            fail_reasons.append("high_gc")
        if max_hp > args.max_homopolymer:
            keep = False
            fail_reasons.append("homopolymer")
        if lc > args.max_low_complexity:
            keep = False
            fail_reasons.append("low_complexity")

        site_num += 1
        rows.append(
            {
                "site_id": f"RSVB_SITE_{site_num:05d}",
                "ref_id": ref.id,
                "start": start0 + 1,
                "end": end0,
                "length": args.window,
                "sequence_ref": window_seq,
                "gc": round(gc, 6),
                "tm_est": round(tm_est, 2),
                "max_homopolymer": max_hp,
                "low_complexity_fraction": round(lc, 6),
                "candidate_keep": keep,
                "candidate_fail_reason": ";".join(fail_reasons),
                "candidate_passes_relaxed_rsvb_defaults": keep,
            }
        )

    df = pd.DataFrame(rows)
    Path(args.output_tsv).parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.output_tsv, sep="\t", index=False)
    print(f"Generated {len(df)} total candidates")
    print(f"Candidates passing hard filters: {int(df['candidate_keep'].sum())}")
    print(f"Wrote -> {args.output_tsv}")


if __name__ == "__main__":
    main()
