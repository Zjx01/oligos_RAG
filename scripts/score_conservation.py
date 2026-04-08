#!/usr/bin/env python3
"""Score candidate sites for conservation / mismatch robustness.

Supports two modes:
- capture_long: legacy long-bait scoring
- rt_primer_25: short RT-primer-aware scoring with emphasis on the 3' end
"""
from __future__ import annotations

import argparse
import math
from pathlib import Path
from typing import Dict, List

import pandas as pd
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment

VALID_BASES = set("ACGT")


def shannon_entropy(chars: List[str]) -> float:
    counts: Dict[str, int] = {}
    valid = [c for c in chars if c in VALID_BASES]
    n = len(valid)
    if n == 0:
        return 0.0
    for c in valid:
        counts[c] = counts.get(c, 0) + 1
    entropy = 0.0
    for count in counts.values():
        p = count / n
        entropy -= p * math.log2(p)
    return entropy


def build_ref_pos_to_aln_col(ref_aligned: str) -> Dict[int, int]:
    mapping: Dict[int, int] = {}
    ref_pos = 0
    for aln_col, base in enumerate(ref_aligned):
        if base != "-":
            ref_pos += 1
            mapping[ref_pos] = aln_col
    return mapping


def extract_window_from_alignment(seq: str, cols: List[int]) -> str:
    return "".join(seq[c] for c in cols)


def mismatch_count(window_aln_seq: str, ref_window: str) -> int | None:
    ungapped = window_aln_seq.replace("-", "")
    if len(ungapped) != len(ref_window):
        return None
    mism = 0
    for a, b in zip(ungapped.upper(), ref_window.upper()):
        if a not in VALID_BASES:
            return None
        if a != b:
            mism += 1
    return mism


def normalize_entropy(ent: float, max_ent: float = 2.0) -> float:
    return min(max(ent / max_ent, 0.0), 1.0)


def exact_suffix_len(seq: str, ref: str) -> int:
    n = min(len(seq), len(ref))
    count = 0
    for i in range(1, n + 1):
        a = seq[-i].upper()
        b = ref[-i].upper()
        if a not in VALID_BASES:
            break
        if a == b:
            count += 1
        else:
            break
    return count


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description="Score candidates for conservation")
    ap.add_argument("--alignment-fasta", required=True)
    ap.add_argument("--ref-id", required=True)
    ap.add_argument("--candidates-tsv", required=True)
    ap.add_argument("--output-tsv", required=True)
    ap.add_argument("--assay-type", choices=["capture_long", "rt_primer_25"], default="rt_primer_25")
    ap.add_argument("--three-prime-nt", type=int, default=8)
    return ap.parse_args()


def main() -> None:
    args = parse_args()
    aln: MultipleSeqAlignment = AlignIO.read(args.alignment_fasta, "fasta")
    candidates = pd.read_csv(args.candidates_tsv, sep="\t")

    id_to_seq = {rec.id: str(rec.seq).upper() for rec in aln}
    if args.ref_id not in id_to_seq:
        raise SystemExit(f"Reference ID {args.ref_id} not found in alignment")

    ref_aln = id_to_seq[args.ref_id]
    ref_map = build_ref_pos_to_aln_col(ref_aln)
    aligned_sequences = list(id_to_seq.values())

    out_rows = []
    for _, row in candidates.iterrows():
        start = int(row["start"])
        end = int(row["end"])
        ref_window = str(row["sequence_ref"]).upper()
        cols = [ref_map[pos] for pos in range(start, end + 1)]
        three_prime_nt = min(int(args.three_prime_nt), len(ref_window))

        col_entropies = []
        valid_windows = 0
        cov0 = cov1 = cov2 = 0
        cov_80pct = cov_90pct = 0
        cov_3p0 = cov_3p1 = cov_terminal = 0
        suffix_fracs: List[float] = []

        for aln_col in cols:
            chars = [seq[aln_col] for seq in aligned_sequences]
            col_entropies.append(shannon_entropy(chars))

        for seq in aligned_sequences:
            window_aln = extract_window_from_alignment(seq, cols)
            mism = mismatch_count(window_aln, ref_window)
            if mism is None:
                continue
            valid_windows += 1
            if mism == 0:
                cov0 += 1
            if mism <= 1:
                cov1 += 1
            if mism <= 2:
                cov2 += 1

            identity = 1.0 - (mism / len(ref_window))
            if identity >= 0.80:
                cov_80pct += 1
            if identity >= 0.90:
                cov_90pct += 1

            seq_3p = window_aln.replace("-", "")[-three_prime_nt:]
            ref_3p = ref_window[-three_prime_nt:]
            mism_3p = 0
            valid_3p = True
            for a, b in zip(seq_3p, ref_3p):
                if a not in VALID_BASES:
                    valid_3p = False
                    break
                if a != b:
                    mism_3p += 1
            if valid_3p:
                if mism_3p == 0:
                    cov_3p0 += 1
                if mism_3p <= 1:
                    cov_3p1 += 1
                if seq_3p[-1] == ref_3p[-1]:
                    cov_terminal += 1
                suffix_fracs.append(exact_suffix_len(seq_3p, ref_3p) / max(len(ref_3p), 1))

        denom = valid_windows if valid_windows else 1
        cov0mm = cov0 / denom
        cov1mm = cov1 / denom
        cov2mm = cov2 / denom
        cov80 = cov_80pct / denom
        cov90 = cov_90pct / denom
        cov3p0 = cov_3p0 / denom
        cov3p1 = cov_3p1 / denom
        covterm = cov_terminal / denom
        mean_suffix_frac = sum(suffix_fracs) / len(suffix_fracs) if suffix_fracs else 0.0
        ent_mean = sum(col_entropies) / len(col_entropies) if col_entropies else 0.0
        ent_max = max(col_entropies) if col_entropies else 0.0

        if args.assay_type == "rt_primer_25":
            robustness = max(
                0.0,
                min(
                    1.0,
                    0.22 * cov0mm
                    + 0.18 * cov1mm
                    + 0.20 * cov3p0
                    + 0.15 * cov3p1
                    + 0.10 * covterm
                    + 0.10 * mean_suffix_frac
                    + 0.05 * cov90
                    - 0.08 * normalize_entropy(ent_mean)
                    - 0.04 * normalize_entropy(ent_max),
                ),
            )
        else:
            robustness = max(
                0.0,
                min(
                    1.0,
                    (0.5 * cov1mm + 0.5 * cov2mm) - 0.15 * normalize_entropy(ent_mean) - 0.10 * normalize_entropy(ent_max),
                ),
            )

        out = row.to_dict()
        out.update(
            {
                "cov_0mm": round(cov0mm, 6),
                "cov_1mm": round(cov1mm, 6),
                "cov_2mm": round(cov2mm, 6),
                "cov_80pct": round(cov80, 6),
                "cov_90pct": round(cov90, 6),
                "cov_3p_0mm": round(cov3p0, 6),
                "cov_3p_1mm": round(cov3p1, 6),
                "cov_terminal_match": round(covterm, 6),
                "mean_exact_3p_suffix_frac": round(mean_suffix_frac, 6),
                "entropy_mean": round(ent_mean, 6),
                "entropy_max": round(ent_max, 6),
                "valid_alignment_windows": valid_windows,
                "robustness_score": round(robustness, 6),
            }
        )
        out_rows.append(out)

    out_df = pd.DataFrame(out_rows)
    Path(args.output_tsv).parent.mkdir(parents=True, exist_ok=True)
    out_df.to_csv(args.output_tsv, sep="\t", index=False)
    print(f"Scored {len(out_df)} candidate sites for conservation")
    print(f"Assay type: {args.assay_type}")
    print(f"Wrote -> {args.output_tsv}")


if __name__ == "__main__":
    main()
