#!/usr/bin/env python3
"""Attach off-target / specificity scores to candidate sites.

This version keeps the original BLAST-driven workflow but adds short-primer-aware
summary fields for rt_primer_25, especially 3'-anchored host hits.
"""
from __future__ import annotations

import argparse
import shutil
import subprocess
import tempfile
from pathlib import Path
from typing import Iterable, List, Optional

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


BLAST_OUTFMT = "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen"
BLAST_COLUMNS = [
    "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend",
    "sstart", "send", "evalue", "bitscore", "qlen", "slen",
]


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description="Score candidates for specificity")
    ap.add_argument("--candidates-tsv", required=True)
    ap.add_argument("--output-tsv", required=True)
    ap.add_argument("--assay-type", choices=["capture_long", "rt_primer_25"], default="rt_primer_25")

    ap.add_argument("--blastn-bin", default="blastn")
    ap.add_argument("--task", default="auto", help="auto, blastn, or blastn-short")
    ap.add_argument("--word-size", type=int, default=None)
    ap.add_argument("--evalue", type=float, default=1000.0)
    ap.add_argument("--max-target-seqs", type=int, default=50)
    ap.add_argument("--num-threads", type=int, default=1)
    ap.add_argument("--min-pident", type=float, default=None)
    ap.add_argument("--min-query-cover", type=float, default=None)
    ap.add_argument("--three-prime-anchor-nt", type=int, default=8)
    ap.add_argument("--query-fasta", default=None)

    ap.add_argument("--human-blast-db", default=None)
    ap.add_argument("--virus-blast-db", default=None)
    ap.add_argument("--human-blast-tsv", default=None)
    ap.add_argument("--virus-blast-tsv", default=None)
    ap.add_argument("--mfeprimer-tsv", default=None)
    return ap.parse_args()


def choose_blast_task(task: str, query_lengths: Iterable[int]) -> str:
    if task != "auto":
        return task
    max_len = max(query_lengths) if query_lengths else 0
    return "blastn-short" if max_len <= 60 else "blastn"


def resolve_defaults(args: argparse.Namespace, query_lengths: List[int]) -> None:
    max_len = max(query_lengths) if query_lengths else 0
    if args.word_size is None:
        args.word_size = 7 if max_len <= 60 else 11
    if args.min_pident is None:
        args.min_pident = 80.0 if args.assay_type == "rt_primer_25" else 85.0
    if args.min_query_cover is None:
        args.min_query_cover = 0.80 if args.assay_type == "rt_primer_25" else 0.80


def write_query_fasta_from_candidates(df: pd.DataFrame, fasta_path: str) -> None:
    required = {"site_id", "sequence_ref"}
    if not required.issubset(df.columns):
        raise SystemExit(f"Candidates TSV must contain columns: {sorted(required)}")
    records = [SeqRecord(Seq(str(row.sequence_ref).upper()), id=str(row.site_id), description="") for row in df.itertuples(index=False)]
    Path(fasta_path).parent.mkdir(parents=True, exist_ok=True)
    SeqIO.write(records, fasta_path, "fasta")


def run_blast(query_fasta: str, db: str, out_path: str, args: argparse.Namespace, query_lengths: List[int]) -> None:
    if shutil.which(args.blastn_bin) is None:
        raise SystemExit(f"blastn binary not found: {args.blastn_bin}")
    task = choose_blast_task(args.task, query_lengths)
    cmd = [
        args.blastn_bin,
        "-query", query_fasta,
        "-db", db,
        "-task", task,
        "-word_size", str(args.word_size),
        "-evalue", str(args.evalue),
        "-dust", "no",
        "-max_target_seqs", str(args.max_target_seqs),
        "-num_threads", str(args.num_threads),
        "-outfmt", BLAST_OUTFMT,
        "-out", out_path,
    ]
    result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False)
    if result.returncode != 0:
        raise SystemExit(
            "blastn failed\n"
            f"CMD: {' '.join(cmd)}\n"
            f"STDOUT:\n{result.stdout.decode('utf-8', errors='replace')}\n"
            f"STDERR:\n{result.stderr.decode('utf-8', errors='replace')}"
        )


def load_blast_hits(blast_tsv: str) -> pd.DataFrame:
    if Path(blast_tsv).stat().st_size == 0:
        return pd.DataFrame(columns=BLAST_COLUMNS)
    return pd.read_csv(blast_tsv, sep="\t", names=BLAST_COLUMNS, comment="#", header=None)


def summarize_blast_hits(blast_tsv: Optional[str], min_pident: float, min_query_cover: float, prefix: str, three_prime_anchor_nt: int) -> pd.DataFrame:
    if blast_tsv is None:
        return pd.DataFrame(columns=[
            "site_id", f"{prefix}_hits", f"{prefix}_best_pident", f"{prefix}_best_qcov",
            f"{prefix}_anchored_hits", f"{prefix}_best_anchored_pident", f"{prefix}_best_anchored_qcov",
            f"{prefix}_nearperfect_hits",
        ])

    hits = load_blast_hits(blast_tsv)
    if hits.empty:
        return pd.DataFrame(columns=[
            "site_id", f"{prefix}_hits", f"{prefix}_best_pident", f"{prefix}_best_qcov",
            f"{prefix}_anchored_hits", f"{prefix}_best_anchored_pident", f"{prefix}_best_anchored_qcov",
            f"{prefix}_nearperfect_hits",
        ])

    hits["qcov"] = hits["length"].astype(float) / hits["qlen"].replace(0, pd.NA).astype(float)
    filt = hits[(hits["pident"].astype(float) >= min_pident) & (hits["qcov"] >= min_query_cover)].copy()
    if filt.empty:
        return pd.DataFrame(columns=[
            "site_id", f"{prefix}_hits", f"{prefix}_best_pident", f"{prefix}_best_qcov",
            f"{prefix}_anchored_hits", f"{prefix}_best_anchored_pident", f"{prefix}_best_anchored_qcov",
            f"{prefix}_nearperfect_hits",
        ])

    filt["anchored"] = (filt["qend"] >= filt["qlen"]) & ((filt["qend"] - filt["qstart"] + 1) >= three_prime_anchor_nt)
    filt["nearperfect"] = (filt["pident"] >= 92.0) & (filt["qcov"] >= 0.92)

    summary = (
        filt.groupby("qseqid", as_index=False)
        .agg(
            hits=("sseqid", "count"),
            best_pident=("pident", "max"),
            best_qcov=("qcov", "max"),
            anchored_hits=("anchored", "sum"),
            best_anchored_pident=("pident", lambda s: float(filt.loc[s.index][filt.loc[s.index, "anchored"]]["pident"].max()) if filt.loc[s.index, "anchored"].any() else 0.0),
            best_anchored_qcov=("qcov", lambda s: float(filt.loc[s.index][filt.loc[s.index, "anchored"]]["qcov"].max()) if filt.loc[s.index, "anchored"].any() else 0.0),
            nearperfect_hits=("nearperfect", "sum"),
        )
        .rename(columns={
            "qseqid": "site_id",
            "hits": f"{prefix}_hits",
            "best_pident": f"{prefix}_best_pident",
            "best_qcov": f"{prefix}_best_qcov",
            "anchored_hits": f"{prefix}_anchored_hits",
            "best_anchored_pident": f"{prefix}_best_anchored_pident",
            "best_anchored_qcov": f"{prefix}_best_anchored_qcov",
            "nearperfect_hits": f"{prefix}_nearperfect_hits",
        })
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


def specificity_score_from_row(row: pd.Series, assay_type: str) -> float:
    human_hits = float(row.get("human_offtarget_hits", 0))
    human_anchored_hits = float(row.get("human_offtarget_anchored_hits", 0))
    human_nearperfect = float(row.get("human_offtarget_nearperfect_hits", 0))
    virus_hits = float(row.get("virus_bg_hits", 0))
    virus_anchored = float(row.get("virus_bg_anchored_hits", 0))
    mfe_penalty = float(row.get("mfeprimer_penalty", 0.0))

    if assay_type == "rt_primer_25":
        penalty = min(
            1.0,
            0.35 * min(human_hits, 5) / 5
            + 0.30 * min(human_anchored_hits, 3) / 3
            + 0.15 * min(human_nearperfect, 3) / 3
            + 0.10 * min(virus_hits, 5) / 5
            + 0.05 * min(virus_anchored, 3) / 3
            + 0.05 * mfe_penalty,
        )
    else:
        penalty = min(
            1.0,
            0.70 * min(human_hits, 5) / 5
            + 0.20 * min(virus_hits, 5) / 5
            + 0.10 * mfe_penalty,
        )
    return max(0.0, 1.0 - penalty)


def main() -> None:
    args = parse_args()
    df = pd.read_csv(args.candidates_tsv, sep="\t")
    query_lengths = [len(str(s)) for s in df.get("sequence_ref", [])]
    resolve_defaults(args, query_lengths)

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

    human_summary = summarize_blast_hits(human_blast_tsv, args.min_pident, args.min_query_cover, prefix="human_offtarget", three_prime_anchor_nt=args.three_prime_anchor_nt)
    virus_summary = summarize_blast_hits(virus_blast_tsv, args.min_pident, args.min_query_cover, prefix="virus_bg", three_prime_anchor_nt=args.three_prime_anchor_nt)
    mfe_summary = load_mfeprimer_summary(args.mfeprimer_tsv)

    out = df.copy()
    for summary in [human_summary, virus_summary, mfe_summary]:
        if not summary.empty:
            out = out.merge(summary, on="site_id", how="left")

    for col in [
        "human_offtarget_hits", "human_offtarget_best_pident", "human_offtarget_best_qcov",
        "human_offtarget_anchored_hits", "human_offtarget_best_anchored_pident", "human_offtarget_best_anchored_qcov",
        "human_offtarget_nearperfect_hits",
        "virus_bg_hits", "virus_bg_best_pident", "virus_bg_best_qcov",
        "virus_bg_anchored_hits", "virus_bg_best_anchored_pident", "virus_bg_best_anchored_qcov",
        "virus_bg_nearperfect_hits",
        "mfeprimer_penalty", "mfeprimer_hits",
    ]:
        if col not in out.columns:
            out[col] = 0.0
        out[col] = out[col].fillna(0.0)

    out["offtarget_penalty"] = out.apply(lambda r: round(1.0 - specificity_score_from_row(r, args.assay_type), 6), axis=1)
    out["specificity_score"] = out.apply(lambda r: round(specificity_score_from_row(r, args.assay_type), 6), axis=1)
    out["specificity_assay_type"] = args.assay_type

    Path(args.output_tsv).parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.output_tsv, sep="\t", index=False)
    print(f"Added specificity scores for {len(out)} candidates")
    print(f"Assay type: {args.assay_type}")
    print(f"Wrote -> {args.output_tsv}")


if __name__ == "__main__":
    main()
