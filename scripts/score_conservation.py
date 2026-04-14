#!/usr/bin/env python3
from __future__ import annotations
import argparse
from pathlib import Path
from typing import Dict, List, Tuple
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
    rows = []
    for _, row in cand.iterrows():
        start = int(row["start"]); end = int(row["end"])
        cols = [ref_map[p] for p in range(start, end + 1)]
        target_chars = list(str(row["sequence_ref"]))
        # For RT primer, primer 3' end lands on LEFT side of target window
        left_cols = [ref_map[p] for p in range(start, min(end, start + args.three_prime_nt - 1) + 1)]
        left_target = list(str(row["sequence_ref"])[:len(left_cols)])

        whole_mm = []
        left_mm = []
        terminal_ok = 0
        whole_cols_bypos: List[List[str]] = [[] for _ in cols]
        for sid, aln_seq in aln_recs:
            seq_chars = [aln_seq[c-1] for c in cols]
            seq_left = [aln_seq[c-1] for c in left_cols]
            whole_mm.append(mismatch_count(target_chars, seq_chars))
            left_mm.append(mismatch_count(left_target, seq_left))
            if bases_match(left_target[0], seq_left[0]):
                terminal_ok += 1
            for i, ch in enumerate(seq_chars):
                whole_cols_bypos[i].append(ch)

        n = len(aln_recs)
        cov0 = sum(mm == 0 for mm in whole_mm) / n
        cov1 = sum(mm <= 1 for mm in whole_mm) / n
        cov2 = sum(mm <= 2 for mm in whole_mm) / n
        cov3p0 = sum(mm == 0 for mm in left_mm) / n
        cov3p1 = sum(mm <= 1 for mm in left_mm) / n
        covterm = terminal_ok / n
        entropy_vals = [shannon_entropy(col) for col in whole_cols_bypos]
        ent_mean = sum(entropy_vals) / len(entropy_vals) if entropy_vals else 0.0
        ent_max = max(entropy_vals) if entropy_vals else 0.0
        mean_exact_3p_suffix_frac = sum(
            sum(1 for a, b in zip(left_target, [aln_seq[c-1] for c in left_cols]) if bases_match(a, b)) / len(left_cols)
            for _, aln_seq in aln_recs
        ) / n

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
