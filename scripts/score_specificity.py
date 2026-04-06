#!/usr/bin/env python3
"""Attach off-target / specificity scores to candidate sites.

Real-tool modes:
1) Run BLASTn directly against one or two prebuilt databases.
2) Parse precomputed BLAST outfmt6 tables.
3) Optionally merge an external MFEprimer summary TSV.

Expected BLAST outfmt 6 columns used here:
qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen
"""
from __future__ import annotations

import argparse
import shutil
import subprocess
import tempfile
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


BLAST_OUTFMT = "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen"


def specificity_from_hits(human_hits: int, virus_hits: int, mfe_penalty: float = 0.0) -> float:
    penalty = min(1.0, 0.7 * min(human_hits, 5) / 5 + 0.3 * min(virus_hits, 5) / 5 + mfe_penalty)
    return max(0.0, 1.0 - penalty)


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description="Score RSVB candidates for specificity")
    ap.add_argument("--candidates-tsv", required=True)
    ap.add_argument("--output-tsv", required=True)

    ap.add_argument("--blastn-bin", default="blastn", help="blastn executable name/path")
    ap.add_argument("--task", default="auto", help="BLAST task: auto, blastn, or blastn-short")
    ap.add_argument("--word-size", type=int, default=7)
    ap.add_argument("--evalue", type=float, default=1000.0)
    ap.add_argument("--max-target-seqs", type=int, default=50)
    ap.add_argument("--num-threads", type=int, default=1)
    ap.add_argument("--min-pident", type=float, default=85.0, help="Minimum percent identity to count a hit")
    ap.add_argument(
        "--min-query-cover",
        type=float,
        default=0.80,
        help="Minimum aligned fraction of query length to count a hit (0-1)",
    )
    ap.add_argument(
        "--query-fasta",
        default=None,
        help="Optional existing query FASTA. If omitted, it is generated from sequence_ref in candidates TSV.",
    )

    ap.add_argument("--human-blast-db", default=None, help="Prebuilt BLAST database for human/background host")
    ap.add_argument("--virus-blast-db", default=None, help="Prebuilt BLAST database for viral background")
    ap.add_argument("--human-blast-tsv", default=None, help="Precomputed BLAST outfmt6 TSV for the human DB")
    ap.add_argument("--virus-blast-tsv", default=None, help="Precomputed BLAST outfmt6 TSV for the viral background DB")

    ap.add_argument(
        "--mfeprimer-tsv",
        default=None,
        help="Optional TSV with site_id and either mfeprimer_penalty or mfeprimer_hits columns",
    )
    return ap.parse_args()


def choose_blast_task(task: str, query_lengths: Iterable[int]) -> str:
    if task != "auto":
        return task
    max_len = max(query_lengths) if query_lengths else 0
    return "blastn-short" if max_len <= 60 else "blastn"


def write_query_fasta_from_candidates(df: pd.DataFrame, fasta_path: str) -> None:
    required = {"site_id", "sequence_ref"}
    if not required.issubset(df.columns):
        raise SystemExit(f"Candidates TSV must contain columns: {sorted(required)}")
    records = [SeqRecord(Seq(str(row.sequence_ref).upper()), id=str(row.site_id), description="") for row in df.itertuples(index=False)]
    Path(fasta_path).parent.mkdir(parents=True, exist_ok=True)
    SeqIO.write(records, fasta_path, "fasta")


def run_blast(
    query_fasta: str,
    db: str,
    out_path: str,
    args: argparse.Namespace,
    query_lengths: List[int],
) -> None:
    if shutil.which(args.blastn_bin) is None:
        raise SystemExit(f"blastn binary not found: {args.blastn_bin}")
    task = choose_blast_task(args.task, query_lengths)
    cmd = [
        args.blastn_bin,
        "-query",
        query_fasta,
        "-db",
        db,
        "-task",
        task,
        "-word_size",
        str(args.word_size),
        "-evalue",
        str(args.evalue),
        "-dust",
        "no",
        "-max_target_seqs",
        str(args.max_target_seqs),
        "-num_threads",
        str(args.num_threads),
        "-outfmt",
        BLAST_OUTFMT,
        "-out",
        out_path,
    ]
    result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False)
    if result.returncode != 0:
        raise SystemExit(
            "blastn failed\n"
            f"CMD: {' '.join(cmd)}\n"
            f"STDOUT:\n{result.stdout.decode('utf-8', errors='replace')}\n"
            f"STDERR:\n{result.stderr.decode('utf-8', errors='replace')}"
        )


BLAST_COLUMNS = [
    "qseqid",
    "sseqid",
    "pident",
    "length",
    "mismatch",
    "gapopen",
    "qstart",
    "qend",
    "sstart",
    "send",
    "evalue",
    "bitscore",
    "qlen",
    "slen",
]


def load_blast_hits(blast_tsv: str) -> pd.DataFrame:
    if Path(blast_tsv).stat().st_size == 0:
        return pd.DataFrame(columns=BLAST_COLUMNS)
    return pd.read_csv(blast_tsv, sep="\t", names=BLAST_COLUMNS, comment="#", header=None)


def summarize_blast_hits(blast_tsv: Optional[str], min_pident: float, min_query_cover: float, prefix: str) -> pd.DataFrame:
    if blast_tsv is None:
        return pd.DataFrame(columns=["site_id", f"{prefix}_hits", f"{prefix}_best_pident", f"{prefix}_best_qcov"])

    hits = load_blast_hits(blast_tsv)
    if hits.empty:
        return pd.DataFrame(columns=["site_id", f"{prefix}_hits", f"{prefix}_best_pident", f"{prefix}_best_qcov"])

    hits["qcov"] = hits["length"].astype(float) / hits["qlen"].replace(0, pd.NA).astype(float)
    filt = hits[(hits["pident"].astype(float) >= min_pident) & (hits["qcov"] >= min_query_cover)].copy()

    if filt.empty:
        return pd.DataFrame(columns=["site_id", f"{prefix}_hits", f"{prefix}_best_pident", f"{prefix}_best_qcov"])

    summary = (
        filt.groupby("qseqid", as_index=False)
        .agg(
            hits=("sseqid", "count"),
            best_pident=("pident", "max"),
            best_qcov=("qcov", "max"),
        )
        .rename(
            columns={
                "qseqid": "site_id",
                "hits": f"{prefix}_hits",
                "best_pident": f"{prefix}_best_pident",
                "best_qcov": f"{prefix}_best_qcov",
            }
        )
    )
    return summary


def load_mfeprimer_summary(path: Optional[str]) -> pd.DataFrame:
    if not path:
        return pd.DataFrame(columns=["site_id", "mfeprimer_penalty", "mfeprimer_hits"])
    df = pd.read_csv(path, sep="\t")
    if "site_id" not in df.columns:
        raise SystemExit("--mfeprimer-tsv must contain a site_id column")
    if "mfeprimer_penalty" not in df.columns:
        if "mfeprimer_hits" in df.columns:
            df["mfeprimer_penalty"] = (df["mfeprimer_hits"].clip(upper=5) / 5.0) * 0.2
        else:
            df["mfeprimer_penalty"] = 0.0
    if "mfeprimer_hits" not in df.columns:
        df["mfeprimer_hits"] = 0
    return df[["site_id", "mfeprimer_penalty", "mfeprimer_hits"]].copy()


def main() -> None:
    args = parse_args()
    df = pd.read_csv(args.candidates_tsv, sep="\t")
    query_lengths = [len(str(s)) for s in df.get("sequence_ref", [])]

    tempdir_obj = None
    query_fasta = args.query_fasta
    if not query_fasta:
        tempdir_obj = tempfile.TemporaryDirectory(prefix="specificity_")
        query_fasta = str(Path(tempdir_obj.name) / "candidates.fa")
        write_query_fasta_from_candidates(df, query_fasta)

    human_blast_tsv = args.human_blast_tsv
    virus_blast_tsv = args.virus_blast_tsv

    if args.human_blast_db:
        if human_blast_tsv is None:
            if tempdir_obj is None:
                tempdir_obj = tempfile.TemporaryDirectory(prefix="specificity_")
            human_blast_tsv = str(Path(tempdir_obj.name) / "human_hits.tsv")
        run_blast(query_fasta, args.human_blast_db, human_blast_tsv, args, query_lengths)

    if args.virus_blast_db:
        if virus_blast_tsv is None:
            if tempdir_obj is None:
                tempdir_obj = tempfile.TemporaryDirectory(prefix="specificity_")
            virus_blast_tsv = str(Path(tempdir_obj.name) / "virus_hits.tsv")
        run_blast(query_fasta, args.virus_blast_db, virus_blast_tsv, args, query_lengths)

    human_summary = summarize_blast_hits(human_blast_tsv, args.min_pident, args.min_query_cover, prefix="human_offtarget")
    virus_summary = summarize_blast_hits(virus_blast_tsv, args.min_pident, args.min_query_cover, prefix="virus_bg")
    mfe_df = load_mfeprimer_summary(args.mfeprimer_tsv)

    df = df.merge(human_summary, on="site_id", how="left")
    df = df.merge(virus_summary, on="site_id", how="left")
    df = df.merge(mfe_df, on="site_id", how="left")

    fill_zero_int = ["human_offtarget_hits", "virus_bg_hits", "mfeprimer_hits"]
    fill_zero_float = [
        "human_offtarget_best_pident",
        "human_offtarget_best_qcov",
        "virus_bg_best_pident",
        "virus_bg_best_qcov",
        "mfeprimer_penalty",
    ]
    for col in fill_zero_int:
        if col not in df.columns:
            df[col] = 0
        df[col] = pd.to_numeric(df[col], errors='coerce').fillna(0).astype(int)
    for col in fill_zero_float:
        if col not in df.columns:
            df[col] = 0.0
        df[col] = pd.to_numeric(df[col], errors='coerce').fillna(0.0).astype(float)

    df["specificity_mode"] = "blastn_direct_or_parsed"
    df["specificity_task"] = choose_blast_task(args.task, query_lengths)
    df["specificity_min_pident"] = args.min_pident
    df["specificity_min_query_cover"] = args.min_query_cover
    df["offtarget_penalty"] = (
        0.7 * df["human_offtarget_hits"].clip(upper=5) / 5
        + 0.3 * df["virus_bg_hits"].clip(upper=5) / 5
        + df["mfeprimer_penalty"]
    ).clip(0, 1)
    df["specificity_score"] = [
        specificity_from_hits(int(h), int(v), float(m))
        for h, v, m in zip(df["human_offtarget_hits"], df["virus_bg_hits"], df["mfeprimer_penalty"])
    ]

    Path(args.output_tsv).parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.output_tsv, sep="\t", index=False)
    print(f"Added specificity scores for {len(df)} candidates")
    print(f"Human BLAST source: {human_blast_tsv or 'none'}")
    print(f"Virus BLAST source: {virus_blast_tsv or 'none'}")
    print(f"Task used: {choose_blast_task(args.task, query_lengths)}")
    print(f"Wrote -> {args.output_tsv}")

    if tempdir_obj is not None:
        tempdir_obj.cleanup()


if __name__ == "__main__":
    main()
