#!/usr/bin/env python3
"""QC and deduplication for RSVB whole-genome FASTA.

Initial MVP behavior:
- filter by minimum length
- filter by maximum N fraction
- deduplicate exact sequences
- write cleaned FASTA and a QC summary TSV
"""
from __future__ import annotations

import argparse
import hashlib
from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Iterable, List

import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


@dataclass
class QcRecord:
    seq_id: str
    length: int
    n_fraction: float
    kept: bool
    reason: str
    duplicate_of: str = ""


def n_fraction(seq: str) -> float:
    seq = seq.upper()
    if not seq:
        return 1.0
    return seq.count("N") / len(seq)


def hash_seq(seq: str) -> str:
    return hashlib.sha256(seq.upper().encode()).hexdigest()


def run_qc(
    records: Iterable[SeqRecord],
    min_length: int,
    max_n_fraction: float,
) -> tuple[List[SeqRecord], pd.DataFrame]:
    kept_records: List[SeqRecord] = []
    qc_rows: List[QcRecord] = []
    seen: dict[str, str] = {}

    for rec in records:
        seq = str(rec.seq).upper()
        length = len(seq)
        nf = n_fraction(seq)

        if length < min_length:
            qc_rows.append(QcRecord(rec.id, length, nf, False, f"length<{min_length}"))
            continue
        if nf > max_n_fraction:
            qc_rows.append(QcRecord(rec.id, length, nf, False, f"N_fraction>{max_n_fraction:.4f}"))
            continue

        digest = hash_seq(seq)
        if digest in seen:
            qc_rows.append(QcRecord(rec.id, length, nf, False, "exact_duplicate", duplicate_of=seen[digest]))
            continue

        seen[digest] = rec.id
        qc_rows.append(QcRecord(rec.id, length, nf, True, "kept"))
        kept_records.append(rec)

    qc_df = pd.DataFrame(asdict(r) for r in qc_rows)
    return kept_records, qc_df


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description="QC RSVB genomes FASTA")
    ap.add_argument("--input-fasta", required=True, help="Input RSVB genome FASTA")
    ap.add_argument("--output-fasta", required=True, help="Cleaned FASTA output")
    ap.add_argument("--qc-tsv", required=True, help="QC summary TSV")
    ap.add_argument("--min-length", type=int, default=14500)
    ap.add_argument("--max-n-fraction", type=float, default=0.02)
    return ap.parse_args()


def main() -> None:
    args = parse_args()
    input_path = Path(args.input_fasta)
    output_path = Path(args.output_fasta)
    qc_path = Path(args.qc_tsv)

    records = list(SeqIO.parse(str(input_path), "fasta"))
    if not records:
        raise SystemExit(f"No sequences found in {input_path}")

    kept_records, qc_df = run_qc(records, args.min_length, args.max_n_fraction)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    qc_path.parent.mkdir(parents=True, exist_ok=True)

    SeqIO.write(kept_records, str(output_path), "fasta")
    qc_df.to_csv(qc_path, sep="\t", index=False)

    kept_n = int(qc_df["kept"].sum()) if not qc_df.empty else 0
    print(f"Input sequences: {len(records)}")
    print(f"Kept sequences: {kept_n}")
    print(f"Removed sequences: {len(records) - kept_n}")
    print(f"Wrote cleaned FASTA -> {output_path}")
    print(f"Wrote QC table -> {qc_path}")


if __name__ == "__main__":
    main()
