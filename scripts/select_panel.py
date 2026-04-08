#!/usr/bin/env python3
"""Greedy panel selection for RSV candidate sites.

This version supports:
- capture_long: long capture baits
- rt_primer_25: short biotinylated RT primers (~25 nt)

For rt_primer_25, the selector now includes stricter prefilters to remove:
- very low-accessibility sites
- sites with 3'-anchored host off-targets
- sites with weak 3'-end conservation / terminal match support

This is meant to make the selector more conservative for RT-priming use cases.
"""
from __future__ import annotations

import argparse
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Tuple

import pandas as pd


DNA_COMP = str.maketrans({
    "A": "T", "C": "G", "G": "C", "T": "A", "U": "A",
    "R": "Y", "Y": "R", "S": "S", "W": "W", "K": "M", "M": "K",
    "B": "V", "D": "H", "H": "D", "V": "B", "N": "N",
})

IUPAC_MATCH: Dict[str, set[str]] = {
    "A": {"A"}, "C": {"C"}, "G": {"G"}, "T": {"T"}, "U": {"T"},
    "R": {"A", "G"}, "Y": {"C", "T"}, "S": {"G", "C"}, "W": {"A", "T"},
    "K": {"G", "T"}, "M": {"A", "C"}, "B": {"C", "G", "T"},
    "D": {"A", "G", "T"}, "H": {"A", "C", "T"}, "V": {"A", "C", "G"},
    "N": {"A", "C", "G", "T"},
}


def clamp01(x: float) -> float:
    return max(0.0, min(1.0, float(x)))


def clean_seq(seq: str) -> str:
    seq = (seq or "").upper().replace(" ", "").replace("\n", "")
    return "".join(ch if ch in IUPAC_MATCH else "N" for ch in seq)


def revcomp(seq: str) -> str:
    return clean_seq(seq).translate(DNA_COMP)[::-1]


def bases_match(a: str, b: str) -> bool:
    sa = IUPAC_MATCH.get(a.upper(), {"A", "C", "G", "T"})
    sb = IUPAC_MATCH.get(b.upper(), {"A", "C", "G", "T"})
    return not sa.isdisjoint(sb)


def score_gc(gc: float, target: float = 0.48, tolerance: float = 0.20) -> float:
    return clamp01(1.0 - abs(gc - target) / tolerance)


def score_tm(tm: float, target: float, tolerance: float) -> float:
    return clamp01(1.0 - abs(tm - target) / tolerance)


def thermo_score(row: pd.Series, assay_type: str) -> float:
    if assay_type == "rt_primer_25":
        gc_part = score_gc(float(row.get("gc", 0.0)), target=0.48, tolerance=0.18)
        tm_part = score_tm(float(row.get("tm_est", 0.0)), target=60.0, tolerance=8.0)
        clamp_part = clamp01(float(row.get("gc_3p5", 0)) / 3.0)
        terminal_part = 1.0 if bool(row.get("terminal_is_gc", False)) else 0.7
        hp_penalty = min(float(row.get("three_prime_max_homopolymer", row.get("max_homopolymer", 0))), 5) / 5
        return clamp01(
            0.28 * gc_part
            + 0.35 * tm_part
            + 0.22 * clamp_part
            + 0.10 * terminal_part
            + 0.05 * (1.0 - hp_penalty)
        )
    gc_part = score_gc(float(row.get("gc", 0.0)))
    tm_part = score_tm(float(row.get("tm_est", 0.0)), target=240.0, tolerance=80.0)
    hp_penalty = min(float(row.get("max_homopolymer", 0)), 8) / 8
    return clamp01(0.45 * gc_part + 0.45 * tm_part + 0.10 * (1.0 - hp_penalty))


def synthesis_score(row: pd.Series, assay_type: str) -> float:
    if assay_type == "rt_primer_25":
        hp_penalty = min(float(row.get("max_homopolymer", 0)), 5) / 5
        hp3_penalty = min(float(row.get("three_prime_max_homopolymer", 0)), 4) / 4
        lc_penalty = clamp01(float(row.get("low_complexity_fraction", 0.0)))
        return clamp01(1.0 - 0.35 * hp_penalty - 0.30 * hp3_penalty - 0.35 * lc_penalty)
    hp_penalty = min(float(row.get("max_homopolymer", 0)), 8) / 8
    lc_penalty = clamp01(float(row.get("low_complexity_fraction", 0.0)))
    return clamp01(1.0 - 0.5 * hp_penalty - 0.5 * lc_penalty)


def site_score(row: pd.Series, assay_type: str) -> float:
    if assay_type == "rt_primer_25":
        return clamp01(
            0.28 * float(row.get("robustness_score", 0.0))
            + 0.24 * float(row.get("accessibility_score", 0.0))
            + 0.18 * thermo_score(row, assay_type)
            + 0.22 * float(row.get("specificity_score", 0.0))
            + 0.08 * synthesis_score(row, assay_type)
        )
    return clamp01(
        0.35 * float(row.get("robustness_score", 0.0))
        + 0.25 * float(row.get("accessibility_score", 0.0))
        + 0.15 * thermo_score(row, assay_type)
        + 0.15 * float(row.get("specificity_score", 0.0))
        + 0.10 * synthesis_score(row, assay_type)
    )


def interval_overlap(a: Tuple[int, int], b: Tuple[int, int]) -> int:
    start = max(a[0], b[0])
    end = min(a[1], b[1])
    return max(0, end - start + 1)


def added_coverage(candidate: Tuple[int, int], selected: List[Tuple[int, int]]) -> int:
    length = candidate[1] - candidate[0] + 1
    if not selected:
        return length
    covered = 0
    for pos in range(candidate[0], candidate[1] + 1):
        if any(s[0] <= pos <= s[1] for s in selected):
            covered += 1
    return length - covered


def iter_offset_matches(seq_a: str, seq_b: str) -> Iterable[Tuple[int, int, int, int, int]]:
    a = clean_seq(seq_a)
    b = clean_seq(seq_b)
    la = len(a)
    lb = len(b)
    for offset in range(-lb + 1, la):
        a_start = max(0, offset)
        b_start = max(0, -offset)
        overlap = min(la - a_start, lb - b_start)
        if overlap <= 0:
            continue
        n_match = 0
        longest_run = 0
        current_run = 0
        longest_terminal_run = 0
        for k in range(overlap):
            ia = a_start + k
            ib = b_start + k
            match = bases_match(a[ia], b[ib])
            if match:
                n_match += 1
                current_run += 1
                longest_run = max(longest_run, current_run)
            else:
                current_run = 0
        current_run = 0
        run_starts_a = 0
        run_starts_b = 0
        for k in range(overlap):
            ia = a_start + k
            ib = b_start + k
            match = bases_match(a[ia], b[ib])
            if match:
                if current_run == 0:
                    run_starts_a = ia
                    run_starts_b = ib
                current_run += 1
            else:
                if current_run > 0:
                    run_end_a = ia - 1
                    run_end_b = ib - 1
                    if run_starts_a == 0 or run_starts_b == 0 or run_end_a == la - 1 or run_end_b == lb - 1:
                        longest_terminal_run = max(longest_terminal_run, current_run)
                current_run = 0
        if current_run > 0:
            run_end_a = a_start + overlap - 1
            run_end_b = b_start + overlap - 1
            if run_starts_a == 0 or run_starts_b == 0 or run_end_a == la - 1 or run_end_b == lb - 1:
                longest_terminal_run = max(longest_terminal_run, current_run)
        yield offset, overlap, n_match, longest_run, longest_terminal_run


def pair_interaction_metrics(seq1: str, seq2: str) -> Dict[str, float]:
    s1 = clean_seq(seq1)
    s2_rc = revcomp(seq2)
    best_fraction = 0.0
    best_run = 0
    best_terminal_run = 0
    best_match_count = 0
    best_overlap = 0
    best_offset = 0
    for offset, overlap, n_match, longest_run, longest_terminal_run in iter_offset_matches(s1, s2_rc):
        frac = n_match / overlap if overlap else 0.0
        if (
            frac > best_fraction
            or (frac == best_fraction and longest_run > best_run)
            or (frac == best_fraction and longest_run == best_run and longest_terminal_run > best_terminal_run)
        ):
            best_fraction = frac
            best_run = longest_run
            best_terminal_run = longest_terminal_run
            best_match_count = n_match
            best_overlap = overlap
            best_offset = offset
    return {
        "best_match_fraction": float(best_fraction),
        "best_match_count": int(best_match_count),
        "best_overlap": int(best_overlap),
        "max_rc_run": int(best_run),
        "max_terminal_rc_run": int(best_terminal_run),
        "best_offset": int(best_offset),
    }


def normalize_run(run_len: float, safe: int, risky: int) -> float:
    if risky <= safe:
        return 1.0 if run_len >= risky else 0.0
    return clamp01((run_len - safe) / (risky - safe))


def pair_penalty_from_metrics(metrics: Dict[str, float], *, safe_run: int, risky_run: int, safe_terminal_run: int, risky_terminal_run: int, safe_fraction: float, risky_fraction: float) -> float:
    run_pen = normalize_run(metrics["max_rc_run"], safe_run, risky_run)
    terminal_pen = normalize_run(metrics["max_terminal_rc_run"], safe_terminal_run, risky_terminal_run)
    if risky_fraction <= safe_fraction:
        frac_pen = 1.0 if metrics["best_match_fraction"] >= risky_fraction else 0.0
    else:
        frac_pen = clamp01((metrics["best_match_fraction"] - safe_fraction) / (risky_fraction - safe_fraction))
    return clamp01(0.45 * run_pen + 0.35 * terminal_pen + 0.20 * frac_pen)


def candidate_panel_interaction(candidate_seq: str, selected_seqs: Sequence[str], *, safe_run: int, risky_run: int, safe_terminal_run: int, risky_terminal_run: int, safe_fraction: float, risky_fraction: float) -> Dict[str, float]:
    if not selected_seqs:
        return {
            "cross_penalty_max": 0.0,
            "cross_penalty_sum": 0.0,
            "max_rc_run": 0,
            "max_terminal_rc_run": 0,
            "max_rc_fraction": 0.0,
            "n_risky_pairs": 0,
        }
    penalties: List[float] = []
    max_run = 0
    max_terminal_run = 0
    max_fraction = 0.0
    n_risky_pairs = 0
    for seq in selected_seqs:
        metrics = pair_interaction_metrics(candidate_seq, seq)
        penalty = pair_penalty_from_metrics(
            metrics,
            safe_run=safe_run,
            risky_run=risky_run,
            safe_terminal_run=safe_terminal_run,
            risky_terminal_run=risky_terminal_run,
            safe_fraction=safe_fraction,
            risky_fraction=risky_fraction,
        )
        penalties.append(penalty)
        if penalty >= 0.5:
            n_risky_pairs += 1
        max_run = max(max_run, int(metrics["max_rc_run"]))
        max_terminal_run = max(max_terminal_run, int(metrics["max_terminal_rc_run"]))
        max_fraction = max(max_fraction, float(metrics["best_match_fraction"]))
    return {
        "cross_penalty_max": max(penalties),
        "cross_penalty_sum": sum(penalties),
        "max_rc_run": max_run,
        "max_terminal_rc_run": max_terminal_run,
        "max_rc_fraction": max_fraction,
        "n_risky_pairs": n_risky_pairs,
    }


def resolve_assay_defaults(args: argparse.Namespace) -> None:
    if args.assay_type == "rt_primer_25":
        defaults = {
            "max_sites": 24,
            "min_site_score": 0.56,
            "min_gap": 250,
            "hard_min_gap": 40,
            "max_overlap_fraction": 0.20,
            "min_new_bp": 12,
            "cross_penalty_weight": 0.22,
            "cross_penalty_sum_weight": 0.08,
            "safe_run": 6,
            "risky_run": 9,
            "safe_terminal_run": 3,
            "risky_terminal_run": 5,
            "safe_match_fraction": 0.50,
            "risky_match_fraction": 0.75,
            "hard_max_terminal_run": 7,
            "hard_max_run": 10,
            # strict rt25 defaults
            "min_accessibility_score": 0.10,
            "min_access_3p_terminal": 0.15,
            "min_robustness_score": 0.90,
            "min_cov_3p_1mm": 0.95,
            "min_cov_terminal_match": 0.95,
            "max_human_hits": 10,
            "max_human_anchored_hits": 0,
            "max_human_nearperfect_hits": 0,
            "max_virus_anchored_hits": 0,
        }
    else:
        defaults = {
            "max_sites": 12,
            "min_site_score": 0.45,
            "min_gap": 500,
            "hard_min_gap": 120,
            "max_overlap_fraction": 0.10,
            "min_new_bp": 65,
            "cross_penalty_weight": 0.15,
            "cross_penalty_sum_weight": 0.05,
            "safe_run": 8,
            "risky_run": 12,
            "safe_terminal_run": 5,
            "risky_terminal_run": 8,
            "safe_match_fraction": 0.55,
            "risky_match_fraction": 0.75,
            "hard_max_terminal_run": 10,
            "hard_max_run": 14,
            "min_accessibility_score": None,
            "min_access_3p_terminal": None,
            "min_robustness_score": None,
            "min_cov_3p_1mm": None,
            "min_cov_terminal_match": None,
            "max_human_hits": None,
            "max_human_anchored_hits": None,
            "max_human_nearperfect_hits": None,
            "max_virus_anchored_hits": None,
        }
    for k, v in defaults.items():
        if getattr(args, k) is None:
            setattr(args, k, v)


def rt25_prefilter_reasons(row: pd.Series, args: argparse.Namespace) -> List[str]:
    reasons: List[str] = []
    if args.assay_type != "rt_primer_25":
        return reasons

    def f(name: str, default: float = 0.0) -> float:
        try:
            return float(row.get(name, default))
        except Exception:
            return float(default)

    if args.min_accessibility_score is not None and f("accessibility_score") < float(args.min_accessibility_score):
        reasons.append(f"low_accessibility<{float(args.min_accessibility_score):.2f}")
    if args.min_access_3p_terminal is not None and f("access_3p_terminal") < float(args.min_access_3p_terminal):
        reasons.append(f"low_3p_access<{float(args.min_access_3p_terminal):.2f}")
    if args.min_robustness_score is not None and f("robustness_score") < float(args.min_robustness_score):
        reasons.append(f"low_robustness<{float(args.min_robustness_score):.2f}")
    if args.min_cov_3p_1mm is not None and f("cov_3p_1mm") < float(args.min_cov_3p_1mm):
        reasons.append(f"low_cov_3p_1mm<{float(args.min_cov_3p_1mm):.2f}")
    if args.min_cov_terminal_match is not None and f("cov_terminal_match") < float(args.min_cov_terminal_match):
        reasons.append(f"low_terminal_match<{float(args.min_cov_terminal_match):.2f}")
    if args.max_human_hits is not None and f("human_offtarget_hits") > float(args.max_human_hits):
        reasons.append(f"human_hits>{int(args.max_human_hits)}")
    if args.max_human_anchored_hits is not None and f("human_offtarget_anchored_hits") > float(args.max_human_anchored_hits):
        reasons.append(f"human_anchored_hits>{int(args.max_human_anchored_hits)}")
    if args.max_human_nearperfect_hits is not None and f("human_offtarget_nearperfect_hits") > float(args.max_human_nearperfect_hits):
        reasons.append(f"human_nearperfect_hits>{int(args.max_human_nearperfect_hits)}")
    if args.max_virus_anchored_hits is not None and f("virus_bg_anchored_hits") > float(args.max_virus_anchored_hits):
        reasons.append(f"virus_anchored_hits>{int(args.max_virus_anchored_hits)}")
    return reasons


def main() -> None:
    ap = argparse.ArgumentParser(description="Select candidate panel")
    ap.add_argument("--scored-tsv", required=True)
    ap.add_argument("--ranked-output", required=True)
    ap.add_argument("--panel-output", required=True)
    ap.add_argument("--assay-type", choices=["capture_long", "rt_primer_25"], default="rt_primer_25")
    ap.add_argument("--max-sites", type=int, default=None)
    ap.add_argument("--min-site-score", type=float, default=None)
    ap.add_argument("--min-gap", type=int, default=None)
    ap.add_argument("--hard-min-gap", type=int, default=None)
    ap.add_argument("--max-overlap-fraction", type=float, default=None)
    ap.add_argument("--min-new-bp", type=int, default=None)
    ap.add_argument("--cross-penalty-weight", type=float, default=None)
    ap.add_argument("--cross-penalty-sum-weight", type=float, default=None)
    ap.add_argument("--safe-run", type=int, default=None)
    ap.add_argument("--risky-run", type=int, default=None)
    ap.add_argument("--safe-terminal-run", type=int, default=None)
    ap.add_argument("--risky-terminal-run", type=int, default=None)
    ap.add_argument("--safe-match-fraction", type=float, default=None)
    ap.add_argument("--risky-match-fraction", type=float, default=None)
    ap.add_argument("--hard-max-terminal-run", type=int, default=None)
    ap.add_argument("--hard-max-run", type=int, default=None)

    # new rt25 filter args
    ap.add_argument("--min-accessibility-score", type=float, default=None)
    ap.add_argument("--min-access-3p-terminal", type=float, default=None)
    ap.add_argument("--min-robustness-score", type=float, default=None)
    ap.add_argument("--min-cov-3p-1mm", type=float, default=None)
    ap.add_argument("--min-cov-terminal-match", type=float, default=None)
    ap.add_argument("--max-human-hits", type=int, default=None)
    ap.add_argument("--max-human-anchored-hits", type=int, default=None)
    ap.add_argument("--max-human-nearperfect-hits", type=int, default=None)
    ap.add_argument("--max-virus-anchored-hits", type=int, default=None)

    args = ap.parse_args()
    resolve_assay_defaults(args)

    df = pd.read_csv(args.scored_tsv, sep="\t")
    if "candidate_keep" in df.columns:
        df = df[df["candidate_keep"].astype(bool)].copy()
    if "sequence_ref" not in df.columns and "seq" in df.columns:
        df["sequence_ref"] = df["seq"]
    if "sequence_ref" not in df.columns:
        raise ValueError("Input TSV must include 'sequence_ref' or 'seq' column for interaction scoring.")

    df["sequence_ref"] = df["sequence_ref"].astype(str).map(clean_seq)
    df["thermo_score"] = df.apply(lambda r: thermo_score(r, args.assay_type), axis=1)
    df["synthesis_score"] = df.apply(lambda r: synthesis_score(r, args.assay_type), axis=1)
    df["site_score"] = df.apply(lambda r: site_score(r, args.assay_type), axis=1)

    # new prefilter diagnostics for rt25
    df["prefilter_fail_reason"] = df.apply(lambda r: ";".join(rt25_prefilter_reasons(r, args)), axis=1)
    df["prefilter_pass"] = df["prefilter_fail_reason"].eq("")

    sort_cols = ["prefilter_pass", "site_score", "robustness_score", "accessibility_score", "specificity_score"]
    ranked = df.sort_values(sort_cols, ascending=False).reset_index(drop=True)
    eligible = ranked[ranked["prefilter_pass"].astype(bool)].copy()

    selected_rows = []
    selected_intervals: List[Tuple[int, int]] = []
    selected_seqs: List[str] = []

    for _, row in eligible.iterrows():
        if len(selected_rows) >= int(args.max_sites):
            break
        if float(row["site_score"]) < float(args.min_site_score):
            continue

        cand = (int(row["start"]), int(row["end"]))
        cand_seq = str(row["sequence_ref"])
        overlap_bp = sum(interval_overlap(cand, s) for s in selected_intervals)
        gain_bp = added_coverage(cand, selected_intervals)
        max_overlap_fraction = 0.0
        if selected_intervals:
            max_overlap_fraction = max(interval_overlap(cand, s) / max(int(row["length"]), 1) for s in selected_intervals)
        if max_overlap_fraction > float(args.max_overlap_fraction):
            continue
        if gain_bp < int(args.min_new_bp):
            continue

        min_distance_penalty = 0.0
        nearest = None
        if selected_intervals:
            nearest = min(min(abs(cand[0] - s[1]), abs(cand[1] - s[0])) for s in selected_intervals)
            if nearest < int(args.hard_min_gap):
                continue
            if nearest < int(args.min_gap):
                min_distance_penalty = (int(args.min_gap) - nearest) / max(int(args.min_gap), 1)

        interaction = candidate_panel_interaction(
            cand_seq,
            selected_seqs,
            safe_run=int(args.safe_run),
            risky_run=int(args.risky_run),
            safe_terminal_run=int(args.safe_terminal_run),
            risky_terminal_run=int(args.risky_terminal_run),
            safe_fraction=float(args.safe_match_fraction),
            risky_fraction=float(args.risky_match_fraction),
        )
        if interaction["max_terminal_rc_run"] > int(args.hard_max_terminal_run):
            continue
        if interaction["max_rc_run"] > int(args.hard_max_run):
            continue

        cross_penalty = float(args.cross_penalty_weight) * float(interaction["cross_penalty_max"]) + float(args.cross_penalty_sum_weight) * float(interaction["cross_penalty_sum"])
        if args.assay_type == "rt_primer_25":
            greedy_gain = (
                0.40 * gain_bp / max(int(row["length"]), 1)
                + 0.32 * float(row["site_score"])
                - 0.08 * overlap_bp / max(int(row["length"]), 1)
                - 0.06 * min_distance_penalty
                - cross_penalty
            )
        else:
            greedy_gain = (
                0.60 * gain_bp / max(int(row["length"]), 1)
                + 0.22 * float(row["site_score"])
                - 0.10 * overlap_bp / max(int(row["length"]), 1)
                - 0.08 * min_distance_penalty
                - cross_penalty
            )
        if greedy_gain <= 0:
            continue

        selected_intervals.append(cand)
        selected_seqs.append(cand_seq)
        rr = row.to_dict()
        rr.update(interaction)
        rr["overlap_bp"] = int(overlap_bp)
        rr["gain_bp"] = int(gain_bp)
        rr["max_overlap_fraction"] = round(max_overlap_fraction, 6)
        rr["nearest_gap_bp"] = int(nearest) if nearest is not None else int(row["length"])
        rr["cross_penalty"] = round(cross_penalty, 6)
        rr["greedy_gain"] = round(greedy_gain, 6)
        rr["reason_selected"] = (
            f"site_score={row['site_score']:.3f}; gain_bp={gain_bp}; overlap_bp={overlap_bp}; max_overlap_fraction={max_overlap_fraction:.3f}; "
            f"nearest_gap_bp={nearest if nearest is not None else int(row['length'])}; spacing_penalty={min_distance_penalty:.3f}; cross_penalty={cross_penalty:.3f}; "
            f"max_rc_run={interaction['max_rc_run']}; max_terminal_rc_run={interaction['max_terminal_rc_run']}; "
            f"max_rc_fraction={interaction['max_rc_fraction']:.3f}; prefilter_pass={bool(row['prefilter_pass'])}"
        )
        selected_rows.append(rr)

    Path(args.ranked_output).parent.mkdir(parents=True, exist_ok=True)
    Path(args.panel_output).parent.mkdir(parents=True, exist_ok=True)
    ranked.to_csv(args.ranked_output, sep="\t", index=False)
    panel_df = pd.DataFrame(selected_rows)
    if not panel_df.empty:
        panel_df.insert(0, "panel_rank", range(1, len(panel_df) + 1))
    panel_df.to_csv(args.panel_output, sep="\t", index=False)

    print(f"Assay type: {args.assay_type}")
    print(f"Ranked candidates: {len(ranked)}")
    print(f"Prefilter-eligible candidates: {int(ranked['prefilter_pass'].sum())}")
    print(f"Selected panel size: {len(panel_df)}")
    print(f"Wrote ranked TSV -> {args.ranked_output}")
    print(f"Wrote panel TSV -> {args.panel_output}")


if __name__ == "__main__":
    main()
