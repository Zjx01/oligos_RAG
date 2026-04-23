#!/usr/bin/env python3
from __future__ import annotations
import argparse
from pathlib import Path
from typing import Dict, List, Tuple
import numpy as np
import pandas as pd
from Bio import SeqIO

DNA_COMP = str.maketrans({
    "A":"T","C":"G","G":"C","T":"A","U":"A",
    "R":"Y","Y":"R","S":"S","W":"W","K":"M","M":"K",
    "B":"V","D":"H","H":"D","V":"B","N":"N",
})
IUPAC = {
    "A":{"A"}, "C":{"C"}, "G":{"G"}, "T":{"T"}, "U":{"T"},
    "R":{"A","G"}, "Y":{"C","T"}, "S":{"G","C"}, "W":{"A","T"},
    "K":{"G","T"}, "M":{"A","C"}, "B":{"C","G","T"}, "D":{"A","G","T"},
    "H":{"A","C","T"}, "V":{"A","C","G"}, "N":{"A","C","G","T"}, "-":set()
}

def clean_seq(seq: str) -> str:
    return "".join(ch if ch in IUPAC else "N" for ch in seq.upper().replace("U","T"))

def revcomp(seq: str) -> str:
    return clean_seq(seq).translate(DNA_COMP)[::-1]

def bases_match(a: str, b: str) -> bool:
    return len(IUPAC.get(a, {"A","C","G","T"}).intersection(IUPAC.get(b, {"A","C","G","T"}))) > 0

def shannon_entropy(col: List[str]) -> float:
    from math import log2
    vals = [x for x in col if x in "ACGT"]
    if not vals:
        return 0.0
    from collections import Counter
    c = Counter(vals)
    total = len(vals)
    ent = 0.0
    for n in c.values():
        p = n / total
        ent -= p * log2(p)
    return ent

def parse_args():
    ap = argparse.ArgumentParser(description="Score conservation in primer orientation")
    ap.add_argument("--alignment-fasta", required=True)
    ap.add_argument("--ref-id", required=True)
    ap.add_argument("--candidates-tsv", required=True)
    ap.add_argument("--output-tsv", required=True)
    ap.add_argument("--assay-type", choices=["capture_long","rt_primer_25"], default="rt_primer_25")
    ap.add_argument("--three-prime-nt", type=int, default=8)
    return ap.parse_args()

def build_ref_map(ref_aln_seq: str) -> Dict[int, int]:
    m = {}
    pos = 0
    for idx, ch in enumerate(ref_aln_seq, start=1):
        if ch != "-":
            pos += 1
            m[pos] = idx
    return m

def mismatch_count(target_chars: List[str], seq_chars: List[str]) -> int:
    mm = 0
    for a, b in zip(target_chars, seq_chars):
        if not bases_match(a, b):
            mm += 1
    return mm

def main():
    args = parse_args()
    recs = list(SeqIO.parse(args.alignment_fasta, "fasta"))
    if not recs:
        raise SystemExit("Empty alignment")
    ref_rec = None
    for r in recs:
        if r.id == args.ref_id:
            ref_rec = r
            break
    if ref_rec is None:
        raise SystemExit(f"Reference ID {args.ref_id} not found in alignment")
    ref_aln = clean_seq(str(ref_rec.seq))
    ref_map = build_ref_map(ref_aln)
    aln_recs = [(r.id, clean_seq(str(r.seq))) for r in recs]
    if any(len(s) != len(ref_aln) for _, s in aln_recs):
        raise SystemExit("Alignment sequences have inconsistent lengths")

    cand = pd.read_csv(args.candidates_tsv, sep="\t")

    # ------------------------------------------------------------------
    # Build alignment matrix once: shape (n_seqs, aln_len), dtype uint8.
    # Each character is stored as its ASCII code so we can compare arrays
    # with NumPy without any Python-level character loops.
    # For IUPAC-aware matching we build a per-base compatibility table:
    # two positions are compatible when their IUPAC sets overlap.
    # We pre-expand every IUPAC code to a 4-bit mask over {A,C,G,T} and
    # compare masks with a bitwise AND — match iff result != 0.
    # ------------------------------------------------------------------
    _BASE_MASK: Dict[str, int] = {
        "A": 0b1000, "C": 0b0100, "G": 0b0010, "T": 0b0001,
        "U": 0b0001,
        "R": 0b1010, "Y": 0b0101, "S": 0b0110, "W": 0b1001,
        "K": 0b0011, "M": 0b1100,
        "B": 0b0111, "D": 0b1011, "H": 0b1101, "V": 0b1110,
        "N": 0b1111, "-": 0b0000,
    }

    aln_seqs = [s for _, s in aln_recs]
    aln_len = len(aln_seqs[0])
    n_seqs = len(aln_seqs)

    # mask_matrix[i, j] = IUPAC 4-bit mask for sequence i, column j
    mask_matrix = np.zeros((n_seqs, aln_len), dtype=np.uint8)
    for i, seq in enumerate(aln_seqs):
        for j, ch in enumerate(seq):
            mask_matrix[i, j] = _BASE_MASK.get(ch, 0b0000)

    rows = []
    for _, row in cand.iterrows():
        start = int(row["start"]); end = int(row["end"])
        # aln column indices (0-based) for this candidate window
        cols_0 = np.array([ref_map[p] - 1 for p in range(start, end + 1)], dtype=np.int32)
        # 3'-proximal columns (left side of target for RT primer)
        tp_len = len(range(start, min(end, start + args.three_prime_nt - 1) + 1))
        left_cols_0 = cols_0[:tp_len]

        # Reference masks for this window
        ref_masks_whole = mask_matrix[0, cols_0]        # shape (window,)
        ref_masks_left  = mask_matrix[0, left_cols_0]  # shape (tp_len,)

        # Alignment sub-matrix for all seqs: shape (n_seqs, window)
        aln_whole = mask_matrix[:, cols_0]              # (n_seqs, window)
        aln_left  = mask_matrix[:, left_cols_0]        # (n_seqs, tp_len)

        # IUPAC match: non-zero bitwise AND means compatible
        match_whole = (aln_whole & ref_masks_whole) != 0   # bool (n_seqs, window)
        match_left  = (aln_left  & ref_masks_left)  != 0  # bool (n_seqs, tp_len)

        # Mismatch counts per sequence
        mm_whole = (~match_whole).sum(axis=1)   # (n_seqs,)
        mm_left  = (~match_left).sum(axis=1)    # (n_seqs,)

        cov0  = float((mm_whole == 0).mean())
        cov1  = float((mm_whole <= 1).mean())
        cov2  = float((mm_whole <= 2).mean())
        cov3p0 = float((mm_left == 0).mean())
        cov3p1 = float((mm_left <= 1).mean())
        covterm = float(match_left[:, 0].mean())   # terminal base match

        # Mean per-sequence fraction of 3'-suffix positions that match
        mean_exact_3p_suffix_frac = float(match_left.mean(axis=1).mean())

        # Shannon entropy per alignment column (ACGT only)
        entropy_vals: List[float] = []
        for col_idx in cols_0:
            col_masks = mask_matrix[:, col_idx]
            # count unambiguous bases only (exactly one bit set)
            counts = {"A": 0, "C": 0, "G": 0, "T": 0}
            for mask in col_masks:
                if mask == 0b1000: counts["A"] += 1
                elif mask == 0b0100: counts["C"] += 1
                elif mask == 0b0010: counts["G"] += 1
                elif mask == 0b0001: counts["T"] += 1
            total = sum(counts.values())
            if total == 0:
                entropy_vals.append(0.0)
                continue
            ent = 0.0
            for cnt in counts.values():
                if cnt > 0:
                    p = cnt / total
                    ent -= p * np.log2(p)
            entropy_vals.append(ent)
        ent_mean = float(np.mean(entropy_vals)) if entropy_vals else 0.0
        ent_max  = float(np.max(entropy_vals))  if entropy_vals else 0.0

        if args.assay_type == "rt_primer_25":
            robustness = max(0.0, min(1.0,
                0.18 * cov0 + 0.22 * cov1 + 0.12 * cov2 +
                0.18 * cov3p0 + 0.18 * cov3p1 + 0.08 * covterm +
                0.04 * mean_exact_3p_suffix_frac
            ))
        else:
            robustness = max(0.0, min(1.0, 0.3*cov0 + 0.35*cov1 + 0.35*cov2))

        out = row.to_dict()
        out.update({
            "cov_0mm": round(cov0, 6),
            "cov_1mm": round(cov1, 6),
            "cov_2mm": round(cov2, 6),
            "cov_3p_0mm": round(cov3p0, 6),
            "cov_3p_1mm": round(cov3p1, 6),
            "cov_terminal_match": round(covterm, 6),
            "mean_exact_3p_suffix_frac": round(mean_exact_3p_suffix_frac, 6),
            "entropy_mean": round(ent_mean, 6),
            "entropy_max": round(ent_max, 6),
            "robustness_score": round(robustness, 6),
        })
        if "rt_primer_seq" not in out:
            out["rt_primer_seq"] = revcomp(out["sequence_ref"])
        rows.append(out)

    outdf = pd.DataFrame(rows)
    Path(args.output_tsv).parent.mkdir(parents=True, exist_ok=True)
    outdf.to_csv(args.output_tsv, sep="\t", index=False)
    print(f"Added conservation scores for {len(outdf)} candidates")
    print(f"Wrote -> {args.output_tsv}")

if __name__ == "__main__":
    main()
