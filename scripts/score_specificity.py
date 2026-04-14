#!/usr/bin/env python3
"""Attach off-target / specificity scores to candidate sites.

Speed-oriented update:
- deduplicates identical candidate sequences before BLAST
- optionally caches BLAST result TSVs
- uses faster BLAST defaults for rt_primer_25
- pushes identity filtering into BLAST when possible
- supports ungapped mode for short primer searches
"""
from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import subprocess
import tempfile
from pathlib import Path
from typing import Iterable, List, Optional, Tuple

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


BLAST_OUTFMT = "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen"
BLAST_COLUMNS = [
    "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend",
    "sstart", "send", "evalue", "bitscore", "qlen", "slen",
]
EMPTY_SUMMARY_COLUMNS = [
    "query_id", "hits", "best_pident", "best_qcov",
    "anchored_hits", "best_anchored_pident", "best_anchored_qcov",
    "nearperfect_hits",
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
    ap.add_argument("--max-target-seqs", type=int, default=None)
    ap.add_argument("--max-hsps", type=int, default=None)
    ap.add_argument("--num-threads", type=int, default=1)
    ap.add_argument("--min-pident", type=float, default=None)
    ap.add_argument("--min-query-cover", type=float, default=None)
    ap.add_argument("--three-prime-anchor-nt", type=int, default=8)
    ap.add_argument("--query-fasta", default=None)
    ap.add_argument("--ungapped", choices=["auto", "yes", "no"], default="auto",
                    help="Use ungapped BLAST alignments. Default: auto (yes for rt_primer_25)")
    ap.add_argument("--blast-perc-identity", choices=["auto", "yes", "no"], default="auto",
                    help="Pass -perc_identity to blastn. Default: auto (yes when min-pident is set)")
    ap.add_argument("--cache-dir", default=None,
                    help="Directory for caching BLAST TSVs. Default: <output dir>/.blast_cache")
    ap.add_argument("--no-deduplicate-queries", action="store_true",
                    help="Disable deduplication of identical candidate sequences before BLAST")

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
        args.min_query_cover = 0.80
    if args.max_target_seqs is None:
        args.max_target_seqs = 10 if args.assay_type == "rt_primer_25" else 50
    if args.max_hsps is None:
        args.max_hsps = 1 if args.assay_type == "rt_primer_25" else 0  # 0 => omit
    if args.cache_dir is None:
        args.cache_dir = str(Path(args.output_tsv).resolve().parent / ".blast_cache")


def validate_db_prefix(db: str) -> None:
    p = Path(db).expanduser()
    if p.is_dir():
        files = sorted(x.name for x in p.iterdir())[:10]
        raise SystemExit(
            f"--blast-db expects a BLAST database prefix, not a directory: {db}\n"
            f"Directory preview: {files}\n"
            f"Example of correct value: {p / 'human_bg'}"
        )


def build_query_table(df: pd.DataFrame, assay_type: str, deduplicate: bool = True) -> Tuple[pd.DataFrame, pd.DataFrame]:
    required = {"site_id", "sequence_ref"}
    if not required.issubset(df.columns):
        raise SystemExit(f"Candidates TSV must contain columns: {sorted(required)}")
    q = df[["site_id", "sequence_ref"]].copy()
    if assay_type == "rt_primer_25":
        if "rt_primer_seq" not in df.columns:
            q["rt_primer_seq"] = q["sequence_ref"].astype(str).str.upper().map(lambda s: str(Seq(s).reverse_complement()))
        else:
            q["rt_primer_seq"] = df["rt_primer_seq"].astype(str).str.upper()
        q["query_seq"] = q["rt_primer_seq"]
    else:
        q["query_seq"] = q["sequence_ref"].astype(str).str.upper()
    q["query_len"] = q["query_seq"].str.len()
    if deduplicate:
        uniq = q[["query_seq", "query_len"]].drop_duplicates().reset_index(drop=True)
        uniq["query_id"] = [f"Q{i:06d}" for i in range(1, len(uniq) + 1)]
        mapping = q.merge(uniq, on=["query_seq", "query_len"], how="left")
        return mapping, uniq[["query_id", "query_seq", "query_len"]]
    else:
        q = q.copy()
        q["query_id"] = q["site_id"].astype(str)
        return q[["site_id", "query_seq", "query_len", "query_id"]], q[["query_id", "query_seq", "query_len"]]


def write_query_fasta(query_df: pd.DataFrame, fasta_path: str) -> None:
    records = [
        SeqRecord(Seq(seq), id=str(qid), description="")
        for qid, seq in zip(query_df["query_id"], query_df["query_seq"])
    ]
    Path(fasta_path).parent.mkdir(parents=True, exist_ok=True)
    SeqIO.write(records, fasta_path, "fasta")


def blast_cache_key(db: str, args: argparse.Namespace, query_df: pd.DataFrame) -> str:
    payload = {
        "db": str(Path(db).expanduser()),
        "assay_type": args.assay_type,
        "task": choose_blast_task(args.task, query_df["query_len"].tolist()),
        "word_size": args.word_size,
        "evalue": args.evalue,
        "max_target_seqs": args.max_target_seqs,
        "max_hsps": args.max_hsps,
        "min_pident": args.min_pident,
        "min_query_cover": args.min_query_cover,
        "three_prime_anchor_nt": args.three_prime_anchor_nt,
        "ungapped": args.ungapped,
        "blast_perc_identity": args.blast_perc_identity,
        "queries": list(zip(query_df["query_id"], query_df["query_seq"])),
    }
    raw = json.dumps(payload, sort_keys=True).encode("utf-8")
    return hashlib.sha256(raw).hexdigest()[:20]


def should_use_ungapped(args: argparse.Namespace) -> bool:
    if args.ungapped == "yes":
        return True
    if args.ungapped == "no":
        return False
    return args.assay_type == "rt_primer_25"


def should_pass_perc_identity(args: argparse.Namespace) -> bool:
    if args.blast_perc_identity == "yes":
        return True
    if args.blast_perc_identity == "no":
        return False
    return args.min_pident is not None


def run_blast(
    query_fasta: str,
    query_df: pd.DataFrame,
    db: str,
    out_path: str,
    args: argparse.Namespace,
) -> None:
    if shutil.which(args.blastn_bin) is None:
        raise SystemExit(f"blastn binary not found: {args.blastn_bin}")
    validate_db_prefix(db)
    task = choose_blast_task(args.task, query_df["query_len"].tolist())
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
    if args.max_hsps and args.max_hsps > 0:
        cmd.extend(["-max_hsps", str(args.max_hsps)])
    if should_use_ungapped(args):
        cmd.append("-ungapped")
    if should_pass_perc_identity(args) and args.min_pident is not None:
        cmd.extend(["-perc_identity", str(args.min_pident)])

    result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False)
    if result.returncode != 0:
        raise SystemExit(
            "blastn failed\n"
            f"CMD: {' '.join(cmd)}\n"
            f"STDOUT:\n{result.stdout.decode('utf-8', errors='replace')}\n"
            f"STDERR:\n{result.stderr.decode('utf-8', errors='replace')}"
        )


def load_blast_hits(blast_tsv: str) -> pd.DataFrame:
    if not blast_tsv or not Path(blast_tsv).exists() or Path(blast_tsv).stat().st_size == 0:
        return pd.DataFrame(columns=BLAST_COLUMNS)
    return pd.read_csv(blast_tsv, sep="\t", names=BLAST_COLUMNS, comment="#", header=None)


def summarize_blast_hits_query_level(
    blast_tsv: Optional[str],
    min_pident: float,
    min_query_cover: float,
    three_prime_anchor_nt: int,
) -> pd.DataFrame:
    if blast_tsv is None:
        return pd.DataFrame(columns=EMPTY_SUMMARY_COLUMNS)

    hits = load_blast_hits(blast_tsv)
    if hits.empty:
        return pd.DataFrame(columns=EMPTY_SUMMARY_COLUMNS)

    hits["pident"] = hits["pident"].astype(float)
    hits["length"] = hits["length"].astype(float)
    hits["qstart"] = hits["qstart"].astype(int)
    hits["qend"] = hits["qend"].astype(int)
    hits["qlen"] = hits["qlen"].astype(float)
    hits["qcov"] = hits["length"] / hits["qlen"].replace(0, pd.NA).astype(float)

    filt = hits[(hits["pident"] >= min_pident) & (hits["qcov"] >= min_query_cover)].copy()
    if filt.empty:
        return pd.DataFrame(columns=EMPTY_SUMMARY_COLUMNS)

    # Query reaches the 3' end and includes at least three_prime_anchor_nt aligned bases there.
    filt["anchored"] = (filt["qend"] >= filt["qlen"]) & (((filt["qend"] - filt["qstart"]) + 1) >= three_prime_anchor_nt)
    filt["nearperfect"] = (filt["pident"] >= 92.0) & (filt["qcov"] >= 0.92)

    rows = []
    for qid, g in filt.groupby("qseqid", sort=False):
        anchored = g[g["anchored"]]
        rows.append({
            "query_id": qid,
            "hits": int(len(g)),
            "best_pident": float(g["pident"].max()),
            "best_qcov": float(g["qcov"].max()),
            "anchored_hits": int(g["anchored"].sum()),
            "best_anchored_pident": float(anchored["pident"].max()) if not anchored.empty else 0.0,
            "best_anchored_qcov": float(anchored["qcov"].max()) if not anchored.empty else 0.0,
            "nearperfect_hits": int(g["nearperfect"].sum()),
        })
    return pd.DataFrame(rows)


def rename_query_summary(summary: pd.DataFrame, prefix: str) -> pd.DataFrame:
    if summary.empty:
        return pd.DataFrame(columns=[
            "query_id", f"{prefix}_hits", f"{prefix}_best_pident", f"{prefix}_best_qcov",
            f"{prefix}_anchored_hits", f"{prefix}_best_anchored_pident", f"{prefix}_best_anchored_qcov",
            f"{prefix}_nearperfect_hits",
        ])
    return summary.rename(columns={
        "hits": f"{prefix}_hits",
        "best_pident": f"{prefix}_best_pident",
        "best_qcov": f"{prefix}_best_qcov",
        "anchored_hits": f"{prefix}_anchored_hits",
        "best_anchored_pident": f"{prefix}_best_anchored_pident",
        "best_anchored_qcov": f"{prefix}_best_anchored_qcov",
        "nearperfect_hits": f"{prefix}_nearperfect_hits",
    })


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


def prepare_blast_tsv(
    label: str,
    blast_db: Optional[str],
    blast_tsv: Optional[str],
    unique_query_df: pd.DataFrame,
    args: argparse.Namespace,
) -> Optional[str]:
    if blast_tsv:
        return blast_tsv
    if not blast_db:
        return None

    cache_dir = Path(args.cache_dir)
    cache_dir.mkdir(parents=True, exist_ok=True)
    cache_key = blast_cache_key(blast_db, args, unique_query_df)
    cached_tsv = cache_dir / f"{label}.{cache_key}.tsv"
    cached_fasta = cache_dir / f"{label}.{cache_key}.fa"

    if cached_tsv.exists():
        return str(cached_tsv)

    if args.query_fasta and not args.no_deduplicate_queries:
        # External query FASTA may not match deduplicated query IDs; only cache FASTA if we created it.
        pass

    with tempfile.TemporaryDirectory(prefix="specificity_") as td:
        td_path = Path(td)
        query_fasta = str(td_path / "queries.fa")
        write_query_fasta(unique_query_df, query_fasta)
        run_blast(query_fasta, unique_query_df, blast_db, str(cached_tsv), args)
        shutil.copyfile(query_fasta, cached_fasta)
    return str(cached_tsv)


def main() -> None:
    args = parse_args()
    df = pd.read_csv(args.candidates_tsv, sep="\t")
    if "site_id" not in df.columns or "sequence_ref" not in df.columns:
        raise SystemExit("Candidates TSV must contain site_id and sequence_ref")

    query_lengths = [len(str(s)) for s in df.get("sequence_ref", [])]
    resolve_defaults(args, query_lengths)

    if args.assay_type == "rt_primer_25" and "rt_primer_seq" not in df.columns:
        df["rt_primer_seq"] = df["sequence_ref"].astype(str).str.upper().map(lambda s: str(Seq(s).reverse_complement()))
    mapping_df, unique_query_df = build_query_table(df, args.assay_type, deduplicate=not args.no_deduplicate_queries)

    human_blast_tsv = prepare_blast_tsv("human", args.human_blast_db, args.human_blast_tsv, unique_query_df, args)
    virus_blast_tsv = prepare_blast_tsv("virus", args.virus_blast_db, args.virus_blast_tsv, unique_query_df, args)

    human_query_summary = summarize_blast_hits_query_level(
        human_blast_tsv, args.min_pident, args.min_query_cover, args.three_prime_anchor_nt
    )
    virus_query_summary = summarize_blast_hits_query_level(
        virus_blast_tsv, args.min_pident, args.min_query_cover, args.three_prime_anchor_nt
    )

    human_query_summary = rename_query_summary(human_query_summary, "human_offtarget")
    virus_query_summary = rename_query_summary(virus_query_summary, "virus_bg")

    out = df.copy()
    out = out.merge(mapping_df[["site_id", "query_id"]], on="site_id", how="left")

    for summary in [human_query_summary, virus_query_summary]:
        if not summary.empty:
            out = out.merge(summary, on="query_id", how="left")

    mfe_summary = load_mfeprimer_summary(args.mfeprimer_tsv)
    if not mfe_summary.empty:
        out = out.merge(mfe_summary, on="site_id", how="left")

    fill_cols = [
        "human_offtarget_hits", "human_offtarget_best_pident", "human_offtarget_best_qcov",
        "human_offtarget_anchored_hits", "human_offtarget_best_anchored_pident", "human_offtarget_best_anchored_qcov",
        "human_offtarget_nearperfect_hits",
        "virus_bg_hits", "virus_bg_best_pident", "virus_bg_best_qcov",
        "virus_bg_anchored_hits", "virus_bg_best_anchored_pident", "virus_bg_best_anchored_qcov",
        "virus_bg_nearperfect_hits",
        "mfeprimer_penalty", "mfeprimer_hits",
    ]
    for col in fill_cols:
        if col not in out.columns:
            out[col] = 0.0
        out[col] = out[col].fillna(0.0)

    out["offtarget_penalty"] = out.apply(lambda r: round(1.0 - specificity_score_from_row(r, args.assay_type), 6), axis=1)
    out["specificity_score"] = out.apply(lambda r: round(specificity_score_from_row(r, args.assay_type), 6), axis=1)
    if args.assay_type == "rt_primer_25" and "rt_primer_seq" not in out.columns:
        out["rt_primer_seq"] = out["sequence_ref"].astype(str).str.upper().map(lambda s: str(Seq(s).reverse_complement()))
    out["specificity_assay_type"] = args.assay_type
    out["specificity_query_count"] = len(unique_query_df)
    out["specificity_deduplicated"] = (not args.no_deduplicate_queries)
    out["specificity_max_target_seqs"] = args.max_target_seqs
    out["specificity_ungapped"] = should_use_ungapped(args)

    if "query_id" in out.columns:
        out = out.drop(columns=["query_id"])

    Path(args.output_tsv).parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.output_tsv, sep="\t", index=False)
    print(f"Added specificity scores for {len(out)} candidates")
    print(f"Assay type: {args.assay_type}")
    print(f"Unique query sequences searched: {len(unique_query_df)}")
    print(f"BLAST max_target_seqs: {args.max_target_seqs}")
    print(f"Ungapped mode: {should_use_ungapped(args)}")
    print(f"Wrote -> {args.output_tsv}")


if __name__ == "__main__":
    main()
