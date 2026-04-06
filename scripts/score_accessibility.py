#!/usr/bin/env python3
"""Attach accessibility scores to candidate sites.

Real-tool modes:
1) Parse an RNAplfold ``*_lunp`` file via ``--lunp-file``.
2) Run RNAplfold directly from a reference FASTA via ``--ref-fasta``.
3) Backward-compatible simple TSV input via ``--unpaired-prob-tsv``.

Recommended usage for capture / RT-extension oligos:
- choose ``--u-length`` to reflect the local accessible stretch you care about
  (for example 15-30 nt), not necessarily the full oligo length.
- compute site-level accessibility by averaging the start-position unpaired
  probabilities across all valid windows of length ``u`` inside each candidate.
"""
from __future__ import annotations

import argparse
import shutil
import subprocess
import tempfile
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

import pandas as pd
from Bio import SeqIO


def heuristic_accessibility(gc: float) -> float:
    target = 0.48
    score = 1.0 - abs(gc - target) / target
    return max(0.0, min(1.0, score))


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description="Score RSVB candidates for accessibility")
    ap.add_argument("--candidates-tsv", required=True)
    ap.add_argument("--output-tsv", required=True)

    group = ap.add_mutually_exclusive_group(required=False)
    group.add_argument(
        "--lunp-file",
        help="RNAplfold *_lunp output file to parse directly",
        default=None,
    )
    group.add_argument(
        "--ref-fasta",
        help="Reference FASTA used to run RNAplfold automatically",
        default=None,
    )
    group.add_argument(
        "--unpaired-prob-tsv",
        help="Backward-compatible TSV with columns: position, p_unpaired",
        default=None,
    )

    ap.add_argument(
        "--u-length",
        type=int,
        default=15,
        help="Unpaired stretch length to score from RNAplfold output (default: 15)",
    )
    ap.add_argument(
        "--rnaplfold-bin",
        default="RNAplfold",
        help="RNAplfold executable name/path when --ref-fasta is used",
    )
    ap.add_argument(
        "--plfold-window",
        type=int,
        default=150,
        help="RNAplfold -W window size when running directly (default: 150)",
    )
    ap.add_argument(
        "--plfold-max-span",
        type=int,
        default=100,
        help="RNAplfold -L max base-pair span when running directly (default: 100)",
    )
    ap.add_argument(
        "--plfold-no-lp",
        action="store_true",
        help="Pass --noLP to RNAplfold",
    )
    ap.add_argument(
        "--keep-temp",
        action="store_true",
        help="Keep temporary RNAplfold working directory for debugging",
    )
    return ap.parse_args()


def read_single_fasta_sequence(fasta_path: str) -> Tuple[str, str]:
    records = list(SeqIO.parse(fasta_path, "fasta"))
    if len(records) != 1:
        raise SystemExit("--ref-fasta must contain exactly one sequence")
    record = records[0]
    return record.id, str(record.seq).upper()


def _candidate_id_from_lunp_name(path: Path) -> str:
    stem = path.name
    if stem.endswith("_lunp"):
        stem = stem[: -len("_lunp")]
    return stem


def find_lunp_file(workdir: Path, expected_seq_id: str | None = None) -> Path:
    candidates = sorted(workdir.glob("*_lunp"))
    if not candidates:
        raise FileNotFoundError(f"No *_lunp file found in {workdir}")
    if expected_seq_id:
        for path in candidates:
            if _candidate_id_from_lunp_name(path) == expected_seq_id:
                return path
    if len(candidates) == 1:
        return candidates[0]
    raise FileNotFoundError(
        f"Multiple *_lunp files found in {workdir}; please inspect manually or use --lunp-file. Files: {[p.name for p in candidates]}"
    )


def run_rnaplfold(ref_fasta: str, args: argparse.Namespace) -> Tuple[Path, Path]:
    if shutil.which(args.rnaplfold_bin) is None:
        raise SystemExit(f"RNAplfold binary not found: {args.rnaplfold_bin}")

    seq_id, _seq = read_single_fasta_sequence(ref_fasta)
    tempdir_obj = tempfile.TemporaryDirectory(prefix="rnaplfold_")
    workdir = Path(tempdir_obj.name)

    cmd = [
        args.rnaplfold_bin,
        "-W",
        str(args.plfold_window),
        "-L",
        str(args.plfold_max_span),
        "-u",
        str(args.u_length),
    ]
    if args.plfold_no_lp:
        cmd.append("--noLP")

    with open(ref_fasta, "rb") as fin:
        result = subprocess.run(
            cmd,
            cwd=workdir,
            stdin=fin,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            check=False,
        )
    if result.returncode != 0:
        stderr = result.stderr.decode("utf-8", errors="replace")
        stdout = result.stdout.decode("utf-8", errors="replace")
        if not args.keep_temp:
            tempdir_obj.cleanup()
        raise SystemExit(
            "RNAplfold failed\n"
            f"CMD: {' '.join(cmd)}\n"
            f"STDOUT:\n{stdout}\nSTDERR:\n{stderr}"
        )

    lunp_path = find_lunp_file(workdir, expected_seq_id=seq_id)
    if args.keep_temp:
        return workdir, lunp_path

    # Persist a copy because TemporaryDirectory will be cleaned up.
    persisted_dir = Path(tempfile.mkdtemp(prefix="rnaplfold_persist_"))
    persisted_lunp = persisted_dir / lunp_path.name
    persisted_lunp.write_bytes(lunp_path.read_bytes())
    tempdir_obj.cleanup()
    return persisted_dir, persisted_lunp


def parse_lunp_file(lunp_file: str, u_length: int) -> Dict[int, float]:
    """Return mapping: 1-based start position -> p_unpaired(length=u_length)."""
    result: Dict[int, float] = {}
    target_col = u_length  # first column is position; remaining cols correspond to lengths 1..u

    with open(lunp_file, "r", encoding="utf-8", errors="replace") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue
            fields = line.split()
            if len(fields) < target_col + 1:
                raise SystemExit(
                    f"{lunp_file} does not contain enough probability columns for --u-length {u_length}. "
                    f"Observed line with {len(fields)-1} probability columns."
                )
            try:
                pos = int(fields[0])
                prob = float(fields[target_col])
            except ValueError:
                # Skip non-numeric rows if present.
                continue
            result[pos] = max(0.0, min(1.0, prob))

    if not result:
        raise SystemExit(f"No numeric RNAplfold accessibility rows parsed from {lunp_file}")
    return result


def parse_simple_probability_tsv(tsv_path: str) -> Dict[int, float]:
    up = pd.read_csv(tsv_path, sep="\t")
    required = {"position", "p_unpaired"}
    if not required.issubset(up.columns):
        raise SystemExit(f"{tsv_path} must contain columns: {sorted(required)}")
    return dict(zip(up["position"].astype(int), up["p_unpaired"].astype(float)))


def compute_site_accessibility(
    start: int,
    end: int,
    unpaired_map: Dict[int, float],
    u_length: int,
) -> Tuple[float, float, int]:
    valid_starts = list(range(start, end - u_length + 2))
    if not valid_starts:
        return 0.0, 0.0, 0
    vals = [float(unpaired_map.get(pos, 0.0)) for pos in valid_starts]
    access_mean = sum(vals) / len(vals)
    access_min = min(vals)
    return access_mean, access_min, len(vals)


def main() -> None:
    args = parse_args()
    df = pd.read_csv(args.candidates_tsv, sep="\t")
    if "sequence_ref" not in df.columns or "start" not in df.columns or "end" not in df.columns:
        raise SystemExit("Candidates TSV must contain at least: site_id, start, end, sequence_ref")

    source_mode = "heuristic_gc_placeholder"
    source_detail = "no_external_source"
    temp_dir_to_report: str | None = None

    if args.lunp_file:
        up_map = parse_lunp_file(args.lunp_file, args.u_length)
        source_mode = "rnaplfold_lunp"
        source_detail = str(args.lunp_file)
    elif args.ref_fasta:
        workdir, lunp_path = run_rnaplfold(args.ref_fasta, args)
        up_map = parse_lunp_file(str(lunp_path), args.u_length)
        source_mode = "rnaplfold_direct"
        source_detail = str(lunp_path)
        temp_dir_to_report = str(workdir)
    elif args.unpaired_prob_tsv:
        up_map = parse_simple_probability_tsv(args.unpaired_prob_tsv)
        source_mode = "external_profile_tsv"
        source_detail = str(args.unpaired_prob_tsv)
    else:
        up_map = {}

    if up_map:
        access_mean: List[float] = []
        access_min: List[float] = []
        access_window_count: List[int] = []
        for _, row in df.iterrows():
            mean_v, min_v, n_v = compute_site_accessibility(
                int(row["start"]),
                int(row["end"]),
                up_map,
                args.u_length,
            )
            access_mean.append(mean_v)
            access_min.append(min_v)
            access_window_count.append(n_v)
        df["access_mean"] = access_mean
        df["access_min"] = access_min
        df["access_window_count"] = access_window_count
        df["accessibility_mode"] = source_mode
        df["accessibility_source"] = source_detail
        df["accessibility_u_length"] = args.u_length
    else:
        df["access_mean"] = df["gc"].astype(float).map(heuristic_accessibility)
        df["access_min"] = (df["access_mean"] * 0.85).round(6)
        df["access_window_count"] = (df["end"] - df["start"] + 1).clip(lower=0)
        df["accessibility_mode"] = source_mode
        df["accessibility_source"] = source_detail
        df["accessibility_u_length"] = args.u_length

    df["accessibility_score"] = (0.7 * df["access_mean"] + 0.3 * df["access_min"]).clip(0, 1)

    Path(args.output_tsv).parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.output_tsv, sep="\t", index=False)
    print(f"Added accessibility scores for {len(df)} candidates")
    print(f"Mode: {df['accessibility_mode'].iloc[0] if len(df) else 'n/a'}")
    print(f"Source: {df['accessibility_source'].iloc[0] if len(df) else 'n/a'}")
    if temp_dir_to_report:
        print(f"RNAplfold temp/work dir: {temp_dir_to_report}")
    print(f"Wrote -> {args.output_tsv}")


if __name__ == "__main__":
    main()
