#!/usr/bin/env python3
"""Generate candidate binding sites from a reference genome.

This version supports two design modes:
- capture_long: long hybrid-capture baits (legacy behavior)
- rt_primer_25: short biotinylated RT primers (~25 nt)

For rt_primer_25, the generator adds short-primer-specific features such as
3' GC clamp, 3' terminal base, and stricter homopolymer / low-complexity
filters.
"""
from __future__ import annotations

import argparse
from pathlib import Path
from typing import Dict, List

import pandas as pd
from Bio import SeqIO

VALID_BASES = set("ACGT")


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


def basic_tm(seq: str) -> float:
    seq = seq.upper()
    n = len(seq)
    gc = seq.count("G") + seq.count("C")
    if n < 14:
        return wallace_tm(seq)
    return 64.9 + 41.0 * (gc - 16.4) / n


def count_valid_bases(seq: str) -> int:
    return sum(1 for b in seq.upper() if b in VALID_BASES)


def gc_count(seq: str) -> int:
    s = seq.upper()
    return s.count("G") + s.count("C")


def resolve_defaults(args: argparse.Namespace) -> Dict[str, float | int]:
    if args.assay_type == "rt_primer_25":
        return {
            "window": 25,
            "step": 1,
            "min_gc": 0.32,
            "max_gc": 0.64,
            "max_homopolymer": 4,
            "max_low_complexity": 0.65,
            "low_complexity_k": 3,
            "min_tm": 52.0,
            "max_tm": 68.0,
            "min_gc_3p5": 1,
            "max_gc_3p5": 4,
            "max_3p_homopolymer": 3,
        }
    return {
        "window": 80,
        "step": 10,
        "min_gc": 0.22,
        "max_gc": 0.68,
        "max_homopolymer": 7,
        "max_low_complexity": 0.88,
        "low_complexity_k": 3,
        "min_tm": 0.0,
        "max_tm": 999.0,
        "min_gc_3p5": 0,
        "max_gc_3p5": 5,
        "max_3p_homopolymer": 99,
    }


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description="Generate candidate binding sites")
    ap.add_argument("--ref-fasta", required=True)
    ap.add_argument("--output-tsv", required=True)
    ap.add_argument("--assay-type", choices=["capture_long", "rt_primer_25"], default="rt_primer_25")
    ap.add_argument("--window", type=int, default=None)
    ap.add_argument("--step", type=int, default=None)
    ap.add_argument("--min-gc", type=float, default=None)
    ap.add_argument("--max-gc", type=float, default=None)
    ap.add_argument("--max-homopolymer", type=int, default=None)
    ap.add_argument("--max-low-complexity", type=float, default=None)
    ap.add_argument("--low-complexity-k", type=int, default=None)
    ap.add_argument("--min-tm", type=float, default=None)
    ap.add_argument("--max-tm", type=float, default=None)
    ap.add_argument("--min-gc-3p5", type=int, default=None)
    ap.add_argument("--max-gc-3p5", type=int, default=None)
    ap.add_argument("--max-3p-homopolymer", type=int, default=None)
    ap.add_argument(
        "--require-terminal-matchable-base",
        action="store_true",
        help="For rt_primer_25, drop candidates with terminal N or non-ACGT at the 3' end",
    )
    return ap.parse_args()


def main() -> None:
    args = parse_args()
    defaults = resolve_defaults(args)
    for key, value in defaults.items():
        if getattr(args, key) is None:
            setattr(args, key, value)

    ref_records = list(SeqIO.parse(args.ref_fasta, "fasta"))
    if len(ref_records) != 1:
        raise SystemExit("--ref-fasta must contain exactly one reference sequence")

    ref = ref_records[0]
    seq = str(ref.seq).upper()
    rows: List[dict] = []

    site_num = 0
    for start0 in range(0, len(seq) - int(args.window) + 1, int(args.step)):
        end0 = start0 + int(args.window)
        window_seq = seq[start0:end0]
        gc = gc_fraction(window_seq)
        max_hp = max_homopolymer(window_seq)
        lc = low_complexity_fraction(window_seq, k=int(args.low_complexity_k))
        tm_est = basic_tm(window_seq)
        seq_3p5 = window_seq[-5:]
        seq_3p8 = window_seq[-8:] if len(window_seq) >= 8 else window_seq
        gc_3p5 = gc_count(seq_3p5)
        tm_3p8 = basic_tm(seq_3p8)
        hp_3p8 = max_homopolymer(seq_3p8)
        terminal_base = window_seq[-1]
        valid_bases = count_valid_bases(window_seq)

        keep = True
        fail_reasons = []
        if valid_bases != len(window_seq):
            keep = False
            fail_reasons.append("non_acgt")
        if gc < float(args.min_gc):
            keep = False
            fail_reasons.append("low_gc")
        if gc > float(args.max_gc):
            keep = False
            fail_reasons.append("high_gc")
        if max_hp > int(args.max_homopolymer):
            keep = False
            fail_reasons.append("homopolymer")
        if lc > float(args.max_low_complexity):
            keep = False
            fail_reasons.append("low_complexity")
        if tm_est < float(args.min_tm):
            keep = False
            fail_reasons.append("low_tm")
        if tm_est > float(args.max_tm):
            keep = False
            fail_reasons.append("high_tm")
        if gc_3p5 < int(args.min_gc_3p5):
            keep = False
            fail_reasons.append("weak_3p_gc_clamp")
        if gc_3p5 > int(args.max_gc_3p5):
            keep = False
            fail_reasons.append("too_strong_3p_gc_clamp")
        if hp_3p8 > int(args.max_3p_homopolymer):
            keep = False
            fail_reasons.append("three_prime_homopolymer")
        if args.require_terminal_matchable_base and terminal_base not in VALID_BASES:
            keep = False
            fail_reasons.append("bad_terminal_base")

        site_num += 1
        rows.append(
            {
                "site_id": f"RSVB_SITE_{site_num:05d}",
                "ref_id": ref.id,
                "assay_type": args.assay_type,
                "start": start0 + 1,
                "end": end0,
                "length": int(args.window),
                "sequence_ref": window_seq,
                "gc": round(gc, 6),
                "tm_est": round(tm_est, 2),
                "tm_3p8_est": round(tm_3p8, 2),
                "max_homopolymer": max_hp,
                "three_prime_max_homopolymer": hp_3p8,
                "low_complexity_fraction": round(lc, 6),
                "gc_3p5": gc_3p5,
                "terminal_base": terminal_base,
                "terminal_is_gc": terminal_base in {"G", "C"},
                "candidate_keep": keep,
                "candidate_fail_reason": ";".join(fail_reasons),
                "candidate_passes_rt25_defaults": keep if args.assay_type == "rt_primer_25" else False,
            }
        )

    df = pd.DataFrame(rows)
    Path(args.output_tsv).parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.output_tsv, sep="\t", index=False)
    print(f"Generated {len(df)} total candidates")
    print(f"Candidates passing hard filters: {int(df['candidate_keep'].sum())}")
    print(f"Assay type: {args.assay_type}")
    print(f"Wrote -> {args.output_tsv}")


if __name__ == "__main__":
    main()
