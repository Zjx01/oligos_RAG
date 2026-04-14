#!/usr/bin/env python3
"""Attach accessibility scores to candidate sites.

Supports:
1) Parse an RNAplfold *_lunp file via --lunp-file.
2) Run RNAplfold directly from a reference FASTA via --ref-fasta.
3) Backward-compatible TSV input via --unpaired-prob-tsv.

For rt_primer_25, the score gives extra weight to accessibility near the 3' end.
"""
from __future__ import annotations

import argparse
import shutil
import subprocess
import tempfile
from pathlib import Path
from typing import Dict, List, Tuple

import pandas as pd
from Bio import SeqIO


def heuristic_accessibility(gc: float) -> float:
    target = 0.48
    score = 1.0 - abs(gc - target) / target
    return max(0.0, min(1.0, score))


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description="Score candidates for accessibility")
    ap.add_argument("--candidates-tsv", required=True)
    ap.add_argument("--output-tsv", required=True)
    ap.add_argument("--assay-type", choices=["capture_long", "rt_primer_25"], default="rt_primer_25")

    group = ap.add_mutually_exclusive_group(required=False)
    group.add_argument("--lunp-file", default=None)
    group.add_argument("--ref-fasta", default=None)
    group.add_argument("--unpaired-prob-tsv", default=None)

    ap.add_argument("--u-length", type=int, default=None, help="Unpaired stretch length used for scoring")
    ap.add_argument("--three-prime-nt", type=int, default=8, help="Number of nucleotides at the 3' end to emphasize")
    ap.add_argument("--rnaplfold-bin", default="RNAplfold")
    ap.add_argument("--plfold-window", type=int, default=150)
    ap.add_argument("--plfold-max-span", type=int, default=100)
    ap.add_argument("--plfold-no-lp", action="store_true")
    ap.add_argument("--keep-temp", action="store_true")
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
    raise FileNotFoundError(f"Multiple *_lunp files found in {workdir}: {[p.name for p in candidates]}")


def run_rnaplfold(ref_fasta: str, args: argparse.Namespace) -> Tuple[Path, Path]:
    if shutil.which(args.rnaplfold_bin) is None:
        raise SystemExit(f"RNAplfold binary not found: {args.rnaplfold_bin}")

    seq_id, _seq = read_single_fasta_sequence(ref_fasta)
    tempdir_obj = tempfile.TemporaryDirectory(prefix="rnaplfold_")
    workdir = Path(tempdir_obj.name)

    cmd = [
        args.rnaplfold_bin,
        "-W", str(args.plfold_window),
        "-L", str(args.plfold_max_span),
        "-u", str(args.u_length),
    ]
    if args.plfold_no_lp:
        cmd.append("--noLP")

    with open(ref_fasta, "rb") as fin:
        result = subprocess.run(cmd, cwd=workdir, stdin=fin, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False)
    if result.returncode != 0:
        stderr = result.stderr.decode("utf-8", errors="replace")
        stdout = result.stdout.decode("utf-8", errors="replace")
        if not args.keep_temp:
            tempdir_obj.cleanup()
        raise SystemExit("RNAplfold failed\n" f"CMD: {' '.join(cmd)}\nSTDOUT:\n{stdout}\nSTDERR:\n{stderr}")

    lunp_path = find_lunp_file(workdir, expected_seq_id=seq_id)
    if args.keep_temp:
        return workdir, lunp_path

    persisted_dir = Path(tempfile.mkdtemp(prefix="rnaplfold_persist_"))
    persisted_lunp = persisted_dir / lunp_path.name
    persisted_lunp.write_bytes(lunp_path.read_bytes())
    tempdir_obj.cleanup()
    return persisted_dir, persisted_lunp


def parse_lunp_file(lunp_file: str, u_length: int) -> Dict[int, float]:
    result: Dict[int, float] = {}
    target_col = u_length
    with open(lunp_file, "r", encoding="utf-8", errors="replace") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue
            fields = line.split()
            if len(fields) < target_col + 1:
                raise SystemExit(
                    f"{lunp_file} does not contain enough probability columns for --u-length {u_length}."
                )
            try:
                pos = int(fields[0])
                prob = float(fields[target_col])
            except ValueError:
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


def compute_site_accessibility(start: int, end: int, unpaired_map: Dict[int, float], u_length: int, three_prime_nt: int) -> Tuple[float, float, int, float, float]:
    valid_starts = list(range(start, end - u_length + 2))
    if not valid_starts:
        return 0.0, 0.0, 0, 0.0, 0.0
    vals = [float(unpaired_map.get(pos, 0.0)) for pos in valid_starts]
    access_mean = sum(vals) / len(vals)
    access_min = min(vals)

    # RT primer 3' end binds to the LEFT side of the target binding region.
    terminal_start = start
    access_3p_terminal = float(unpaired_map.get(terminal_start, 0.0)) if terminal_start <= end else 0.0

    near_3p_limit = min(end - u_length + 1, start + max(three_prime_nt - u_length, 0))
    near_3p_starts = [pos for pos in valid_starts if pos <= max(start, near_3p_limit)]
    if near_3p_starts:
        access_3p_near_mean = sum(float(unpaired_map.get(pos, 0.0)) for pos in near_3p_starts) / len(near_3p_starts)
    else:
        access_3p_near_mean = access_3p_terminal
    return access_mean, access_min, len(vals), access_3p_terminal, access_3p_near_mean


def main() -> None:
    args = parse_args()
    if args.u_length is None:
        args.u_length = 8 if args.assay_type == "rt_primer_25" else 15

    df = pd.read_csv(args.candidates_tsv, sep="\t")
    if "sequence_ref" not in df.columns or "start" not in df.columns or "end" not in df.columns:
        raise SystemExit("Candidates TSV must contain at least: start, end, sequence_ref")
    if "rt_primer_seq" not in df.columns and args.assay_type == "rt_primer_25":
        df["rt_primer_seq"] = df["sequence_ref"].astype(str).map(lambda s: clean_seq(s).translate(DNA_COMP)[::-1])

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
        access_3p_terminal: List[float] = []
        access_3p_near_mean: List[float] = []
        for _, row in df.iterrows():
            mean_v, min_v, n_v, term_v, near_v = compute_site_accessibility(
                int(row["start"]),
                int(row["end"]),
                up_map,
                int(args.u_length),
                int(args.three_prime_nt),
            )
            access_mean.append(mean_v)
            access_min.append(min_v)
            access_window_count.append(n_v)
            access_3p_terminal.append(term_v)
            access_3p_near_mean.append(near_v)
        df["access_mean"] = access_mean
        df["access_min"] = access_min
        df["access_window_count"] = access_window_count
        df["access_3p_terminal"] = access_3p_terminal
        df["access_3p_near_mean"] = access_3p_near_mean
        df["accessibility_mode"] = source_mode
        df["accessibility_source"] = source_detail
        df["accessibility_u_length"] = args.u_length
    else:
        df["access_mean"] = df["gc"].astype(float).map(heuristic_accessibility)
        df["access_min"] = (df["access_mean"] * 0.85).round(6)
        df["access_window_count"] = (df["end"] - df["start"] + 1).clip(lower=0)
        df["access_3p_terminal"] = (df["access_mean"] * 0.95).round(6)
        df["access_3p_near_mean"] = (df["access_mean"] * 0.92).round(6)
        df["accessibility_mode"] = source_mode
        df["accessibility_source"] = source_detail
        df["accessibility_u_length"] = args.u_length

    if args.assay_type == "rt_primer_25":
        df["accessibility_score"] = (
            0.30 * df["access_mean"]
            + 0.15 * df["access_min"]
            + 0.35 * df["access_3p_terminal"]
            + 0.20 * df["access_3p_near_mean"]
        ).clip(0, 1)
    else:
        df["accessibility_score"] = (0.7 * df["access_mean"] + 0.3 * df["access_min"]).clip(0, 1)

    Path(args.output_tsv).parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.output_tsv, sep="\t", index=False)
    print(f"Added accessibility scores for {len(df)} candidates")
    print(f"Assay type: {args.assay_type}")
    print(f"Mode: {df['accessibility_mode'].iloc[0] if len(df) else 'n/a'}")
    if temp_dir_to_report:
        print(f"RNAplfold temp/work dir: {temp_dir_to_report}")
    print(f"Wrote -> {args.output_tsv}")


if __name__ == "__main__":
    main()
