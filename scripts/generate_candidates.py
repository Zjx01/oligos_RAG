#!/usr/bin/env python3
from __future__ import annotations
import argparse
import math
from pathlib import Path
from typing import List, Tuple
from Bio import SeqIO
import pandas as pd

DNA_COMP = str.maketrans({
    "A":"T","C":"G","G":"C","T":"A","U":"A",
    "R":"Y","Y":"R","S":"S","W":"W","K":"M","M":"K",
    "B":"V","D":"H","H":"D","V":"B","N":"N",
})

def clean_seq(seq: str) -> str:
    return "".join(ch if ch in "ACGTRYSWKMBDHVN" else "N" for ch in seq.upper().replace("U","T"))

def revcomp(seq: str) -> str:
    return clean_seq(seq).translate(DNA_COMP)[::-1]

def read_single_fasta(fasta: str) -> Tuple[str, str]:
    recs = list(SeqIO.parse(fasta, "fasta"))
    if len(recs) != 1:
        raise SystemExit("--ref-fasta must contain exactly one sequence")
    return recs[0].id, clean_seq(str(recs[0].seq))

def gc_frac(seq: str) -> float:
    seq = clean_seq(seq)
    if not seq:
        return 0.0
    return (seq.count("G") + seq.count("C")) / len(seq)

# ---------------------------------------------------------------------------
# Nearest-neighbour Tm (SantaLucia 1998, unified parameters).
# Reference: SantaLucia J (1998) PNAS 95:1460-1465.
# Conditions: 50 mM Na+, 250 nM oligo (default for primer design).
# The Wallace rule (2*AT + 4*GC) is only valid for oligos ≤14 nt and
# produces errors of 5-10 °C for 25-mers; NN is used instead.
# ---------------------------------------------------------------------------
_NN_DH: dict[str, float] = {          # kcal/mol
    "AA": -7.9,  "TT": -7.9,
    "AT": -7.2,  "TA": -7.2,
    "CA": -8.5,  "TG": -8.5,
    "GT": -8.4,  "AC": -8.4,
    "CT": -7.8,  "AG": -7.8,
    "GA": -8.2,  "TC": -8.2,
    "CG": -10.6, "GC": -9.8,
    "GG": -8.0,  "CC": -8.0,
}
_NN_DS: dict[str, float] = {          # cal/mol/K
    "AA": -22.2, "TT": -22.2,
    "AT": -20.4, "TA": -21.3,
    "CA": -22.7, "TG": -22.7,
    "GT": -22.4, "AC": -22.4,
    "CT": -21.0, "AG": -21.0,
    "GA": -22.2, "TC": -22.2,
    "CG": -27.2, "GC": -24.4,
    "GG": -19.9, "CC": -19.9,
}
_R = 1.987          # cal/mol/K
_OLIGO_CONC = 250e-9   # 250 nM total strand
_SALT_CONC   = 0.05    # 50 mM Na+

def nn_tm(seq: str,
          oligo_conc: float = _OLIGO_CONC,
          salt_conc: float = _SALT_CONC) -> float:
    """Nearest-neighbour Tm (°C) using SantaLucia 1998 parameters.

    Falls back to Wallace rule only when the sequence contains ambiguity
    codes that prevent NN lookup (rare for a cleaned reference window).
    Salt correction applied via the empirical formula of Owczarzy 2004.
    """
    seq = clean_seq(seq).upper()
    # Strip any non-ACGT before NN; if >20% ambiguous, fall back.
    acgt = sum(seq.count(b) for b in "ACGT")
    if acgt < len(seq) * 0.80 or len(seq) < 2:
        at = seq.count("A") + seq.count("T")
        gc = seq.count("G") + seq.count("C")
        return float(2 * at + 4 * gc)  # Wallace fallback

    dH = 0.0  # kcal/mol
    dS = 0.0  # cal/mol/K

    # Initiation parameters (SantaLucia 1998 Table 2)
    for end_base in (seq[0], seq[-1]):
        if end_base in "GC":
            dH += 0.1;  dS += -2.8
        else:                           # A or T
            dH += 2.3;  dS += 4.1

    for i in range(len(seq) - 1):
        dinuc = seq[i:i+2]
        if dinuc in _NN_DH:
            dH += _NN_DH[dinuc]
            dS += _NN_DS[dinuc]
        else:
            # Ambiguous dinucleotide — skip (conservative, minor effect)
            pass

    dH_cal = dH * 1000.0   # convert to cal/mol
    # Tm in Kelvin (non-self-complementary duplex)
    tm_k = dH_cal / (dS + _R * math.log(oligo_conc / 4.0))

    # Salt correction (Owczarzy 2004, Eq. 22 — monovalent)
    gc_frac_val = (seq.count("G") + seq.count("C")) / max(len(seq), 1)
    ln_salt = math.log(salt_conc)
    inv_tm_corr = (1.0 / tm_k) + (
        (4.29 * gc_frac_val - 3.95) * 1e-5 * ln_salt
        + 9.40e-6 * ln_salt ** 2
    )
    tm_celsius = (1.0 / inv_tm_corr) - 273.15
    return round(tm_celsius, 2)

def max_homopolymer(seq: str) -> int:
    seq = clean_seq(seq)
    best = cur = 0
    prev = None
    for ch in seq:
        if ch == prev:
            cur += 1
        else:
            cur = 1
            prev = ch
        best = max(best, cur)
    return best

def low_complexity_fraction(seq: str, k: int = 3) -> float:
    seq = clean_seq(seq)
    if len(seq) < k or k <= 0:
        return 0.0
    kmers = [seq[i:i+k] for i in range(len(seq)-k+1)]
    if not kmers:
        return 0.0
    from collections import Counter
    c = Counter(kmers)
    return max(c.values()) / len(kmers)

def parse_args():
    ap = argparse.ArgumentParser(description="Generate candidate target binding sites and RT primers")
    ap.add_argument("--ref-fasta", required=True)
    ap.add_argument("--output-tsv", required=True)
    ap.add_argument("--assay-type", choices=["capture_long","rt_primer_25"], default="rt_primer_25")
    ap.add_argument("--window", type=int, default=None)
    ap.add_argument("--step", type=int, default=None)
    ap.add_argument("--min-gc", type=float, default=0.32)
    ap.add_argument("--max-gc", type=float, default=0.64)
    ap.add_argument("--min-tm", type=float, default=54.0)
    ap.add_argument("--max-tm", type=float, default=64.0)
    ap.add_argument("--max-homopolymer", type=int, default=5)
    ap.add_argument("--max-low-complexity", type=float, default=0.85)
    ap.add_argument("--low-complexity-k", type=int, default=3)
    ap.add_argument("--min-gc-3p5", type=int, default=1)
    ap.add_argument("--max-gc-3p5", type=int, default=4)
    ap.add_argument("--max-3p-homopolymer", type=int, default=3)
    ap.add_argument("--require-terminal-matchable-base", action="store_true")
    return ap.parse_args()

def main():
    args = parse_args()
    if args.window is None:
        args.window = 25 if args.assay_type == "rt_primer_25" else 80
    if args.step is None:
        args.step = 1 if args.assay_type == "rt_primer_25" else 10

    ref_id, ref_seq = read_single_fasta(args.ref_fasta)

    rows = []
    for start0 in range(0, len(ref_seq) - args.window + 1, args.step):
        end0 = start0 + args.window
        target_seq = ref_seq[start0:end0]
        primer_seq = revcomp(target_seq)
        term_base = primer_seq[-1]
        gc = gc_frac(primer_seq)
        tm = nn_tm(primer_seq)
        gc_3p5 = primer_seq[-5:].count("G") + primer_seq[-5:].count("C")
        three_prime_hp = max_homopolymer(primer_seq[-8:])
        hp = max_homopolymer(primer_seq)
        lc = low_complexity_fraction(primer_seq, args.low_complexity_k)

        keep = True
        reasons: List[str] = []
        if gc < args.min_gc:
            keep = False; reasons.append(f"low_gc<{args.min_gc}")
        if gc > args.max_gc:
            keep = False; reasons.append(f"high_gc>{args.max_gc}")
        if tm < args.min_tm:
            keep = False; reasons.append(f"low_tm<{args.min_tm}")
        if tm > args.max_tm:
            keep = False; reasons.append(f"high_tm>{args.max_tm}")
        if hp > args.max_homopolymer:
            keep = False; reasons.append(f"homopolymer>{args.max_homopolymer}")
        if lc > args.max_low_complexity:
            keep = False; reasons.append(f"low_complexity>{args.max_low_complexity}")
        if gc_3p5 < args.min_gc_3p5:
            keep = False; reasons.append(f"low_gc_3p5<{args.min_gc_3p5}")
        if gc_3p5 > args.max_gc_3p5:
            keep = False; reasons.append(f"high_gc_3p5>{args.max_gc_3p5}")
        if three_prime_hp > args.max_3p_homopolymer:
            keep = False; reasons.append(f"3p_homopolymer>{args.max_3p_homopolymer}")
        if args.require_terminal_matchable_base and term_base not in "ACGT":
            keep = False; reasons.append("noncanonical_terminal")

        rows.append({
            "ref_id": ref_id,
            "site_id": f"{ref_id}_SITE_{start0+1:05d}",
            "start": start0 + 1,
            "end": end0,
            "length": args.window,
            "sequence_ref": target_seq,
            "rt_primer_seq": primer_seq,
            "gc": round(gc, 6),
            "tm_est": round(tm, 3),
            "rt_primer_gc": round(gc, 6),
            "rt_primer_tm": round(tm, 3),
            "rt_primer_terminal_base": term_base,
            "rt_primer_terminal_is_gc": term_base in ("G","C"),
            "rt_primer_gc_3p5": gc_3p5,
            "max_homopolymer": hp,
            "three_prime_max_homopolymer": three_prime_hp,
            "low_complexity_fraction": round(lc, 6),
            "candidate_keep": keep,
            "candidate_fail_reason": ";".join(reasons),
        })

    df = pd.DataFrame(rows)
    Path(args.output_tsv).parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.output_tsv, sep="\t", index=False)
    print(f"Generated {len(df)} candidates")
    print(f"Kept {int(df['candidate_keep'].sum())} after hard filters")
    print(f"Wrote -> {args.output_tsv}")

if __name__ == "__main__":
    main()
