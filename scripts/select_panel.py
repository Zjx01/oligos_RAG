#!/usr/bin/env python3
"""Even-coverage panel selection for RSV / HAV short RT primers.

This selector is designed for:
- rt_primer_25 panels where the goal is to distribute primers more evenly
  across the target genome instead of clustering on the highest-scoring local region.

Main ideas:
1. Apply hard rt25 prefilters first.
2. Compute a per-site score (robustness + accessibility + thermo + specificity + synthesis).
3. Create ideal slot positions across the genome based on max_sites.
4. Select at most one primer per slot in the first pass, favoring:
   - high site score
   - proximity to the slot center
   - low panel interaction
5. Optional second pass fills leftover slots / extra capacity using gap-filling heuristics.
6. For rt_primer_25, output `rt_primer_seq = revcomp(sequence_ref)` for ordering.
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


def base_site_score(row: pd.Series, assay_type: str) -> float:
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


def compute_centers(df: pd.DataFrame) -> pd.Series:
    return ((df["start"].astype(float) + df["end"].astype(float)) / 2.0).round(3)


def build_slot_centers(genome_start: int, genome_end: int, max_sites: int) -> List[float]:
    if max_sites <= 1 or genome_end <= genome_start:
        return [float((genome_start + genome_end) / 2.0)]
    step = (genome_end - genome_start) / float(max_sites - 1)
    return [genome_start + i * step for i in range(max_sites)]


def slot_proximity(center: float, slot_center: float, target_gap: float) -> float:
    tol = max(target_gap * 0.75, 100.0)
    return clamp01(1.0 - (abs(center - slot_center) / tol))


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
    return {
        "max_rc_run": int(best_run),
        "max_terminal_rc_run": int(best_terminal_run),
        "max_rc_fraction": float(best_fraction),
    }


def normalize_run(run_len: float, safe: int, risky: int) -> float:
    if risky <= safe:
        return 1.0 if run_len >= risky else 0.0
    return clamp01((run_len - safe) / (risky - safe))


def pair_penalty(metrics: Dict[str, float], safe_run: int, risky_run: int, safe_terminal_run: int, risky_terminal_run: int, safe_fraction: float, risky_fraction: float) -> float:
    run_pen = normalize_run(metrics["max_rc_run"], safe_run, risky_run)
    terminal_pen = normalize_run(metrics["max_terminal_rc_run"], safe_terminal_run, risky_terminal_run)
    frac_pen = clamp01((metrics["max_rc_fraction"] - safe_fraction) / max(risky_fraction - safe_fraction, 1e-6))
    return clamp01(0.45 * run_pen + 0.35 * terminal_pen + 0.20 * frac_pen)


def panel_interaction(candidate_seq: str, selected_seqs: Sequence[str], args: argparse.Namespace) -> Dict[str, float]:
    if not selected_seqs:
        return {
            "cross_penalty_max": 0.0,
            "cross_penalty_sum": 0.0,
            "max_rc_run": 0,
            "max_terminal_rc_run": 0,
            "max_rc_fraction": 0.0,
            "n_risky_pairs": 0,
        }
    penalties = []
    max_run = 0
    max_terminal_run = 0
    max_fraction = 0.0
    n_risky = 0
    for seq in selected_seqs:
        metrics = pair_interaction_metrics(candidate_seq, seq)
        pen = pair_penalty(
            metrics,
            args.safe_run,
            args.risky_run,
            args.safe_terminal_run,
            args.risky_terminal_run,
            args.safe_match_fraction,
            args.risky_match_fraction,
        )
        penalties.append(pen)
        if pen >= 0.5:
            n_risky += 1
        max_run = max(max_run, metrics["max_rc_run"])
        max_terminal_run = max(max_terminal_run, metrics["max_terminal_rc_run"])
        max_fraction = max(max_fraction, metrics["max_rc_fraction"])
    return {
        "cross_penalty_max": max(penalties),
        "cross_penalty_sum": sum(penalties),
        "max_rc_run": max_run,
        "max_terminal_rc_run": max_terminal_run,
        "max_rc_fraction": max_fraction,
        "n_risky_pairs": n_risky,
    }


def min_gap_to_selected(cand: Tuple[int, int], selected: List[Tuple[int, int]]) -> int | None:
    if not selected:
        return None
    gaps = [min(abs(cand[0] - s[1]), abs(cand[1] - s[0])) for s in selected]
    return min(gaps) if gaps else None


def added_bp(cand: Tuple[int, int], selected: List[Tuple[int, int]]) -> int:
    if not selected:
        return cand[1] - cand[0] + 1
    covered = 0
    for pos in range(cand[0], cand[1] + 1):
        if any(s[0] <= pos <= s[1] for s in selected):
            covered += 1
    return (cand[1] - cand[0] + 1) - covered


def largest_selected_gap(selected_centers: List[float], genome_start: float, genome_end: float) -> float:
    pts = [genome_start] + sorted(selected_centers) + [genome_end]
    return max((pts[i+1] - pts[i]) for i in range(len(pts) - 1))


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
            "min_accessibility_score": 0.10,
            "min_access_3p_terminal": 0.15,
            "min_robustness_score": 0.90,
            "min_cov_3p_1mm": 0.95,
            "min_cov_terminal_match": 0.95,
            "max_human_hits": 10,
            "max_human_anchored_hits": 0,
            "max_human_nearperfect_hits": 0,
            "max_virus_anchored_hits": 0,
            "slot_weight": 0.35,
            "gap_fill_weight": 0.20,
            "second_pass": True,
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
            "slot_weight": 0.20,
            "gap_fill_weight": 0.10,
            "second_pass": True,
        }
    for k, v in defaults.items():
        if not hasattr(args, k) or getattr(args, k) is None:
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


def candidate_is_feasible(row: pd.Series, selected_intervals: List[Tuple[int, int]], selected_seqs: List[str], args: argparse.Namespace) -> Tuple[bool, Dict[str, float]]:
    cand = (int(row["start"]), int(row["end"]))
    length = int(row["length"])
    overlap_bp = sum(interval_overlap(cand, s) for s in selected_intervals)
    max_overlap_fraction = 0.0
    if selected_intervals:
        max_overlap_fraction = max(interval_overlap(cand, s) / max(length, 1) for s in selected_intervals)
    if max_overlap_fraction > float(args.max_overlap_fraction):
        return False, {"overlap_bp": overlap_bp, "max_overlap_fraction": max_overlap_fraction}
    gain_bp = added_bp(cand, selected_intervals)
    if gain_bp < int(args.min_new_bp):
        return False, {"overlap_bp": overlap_bp, "max_overlap_fraction": max_overlap_fraction, "gain_bp": gain_bp}
    nearest_gap = min_gap_to_selected(cand, selected_intervals)
    if nearest_gap is not None and nearest_gap < int(args.hard_min_gap):
        return False, {"overlap_bp": overlap_bp, "max_overlap_fraction": max_overlap_fraction, "gain_bp": gain_bp, "nearest_gap_bp": nearest_gap}
    interaction = panel_interaction(str(row["sequence_ref"]), selected_seqs, args)
    if interaction["max_terminal_rc_run"] > int(args.hard_max_terminal_run):
        return False, {"overlap_bp": overlap_bp, "max_overlap_fraction": max_overlap_fraction, "gain_bp": gain_bp, "nearest_gap_bp": nearest_gap if nearest_gap is not None else length, **interaction}
    if interaction["max_rc_run"] > int(args.hard_max_run):
        return False, {"overlap_bp": overlap_bp, "max_overlap_fraction": max_overlap_fraction, "gain_bp": gain_bp, "nearest_gap_bp": nearest_gap if nearest_gap is not None else length, **interaction}
    return True, {
        "overlap_bp": overlap_bp,
        "max_overlap_fraction": max_overlap_fraction,
        "gain_bp": gain_bp,
        "nearest_gap_bp": nearest_gap if nearest_gap is not None else length,
        **interaction,
    }


def first_pass_slot_selection(eligible: pd.DataFrame, slot_centers: List[float], target_gap: float, genome_start: float, genome_end: float, args: argparse.Namespace):
    selected_rows: List[Dict] = []
    selected_site_ids = set()
    selected_intervals: List[Tuple[int, int]] = []
    selected_seqs: List[str] = []
    selected_centers: List[float] = []

    for slot_idx, slot_center in enumerate(slot_centers, start=1):
        remaining = eligible[~eligible["site_id"].isin(selected_site_ids)].copy()
        if remaining.empty:
            break

        remaining["slot_proximity"] = remaining["center"].map(lambda c: slot_proximity(float(c), float(slot_center), target_gap))
        # prefer candidates close to the slot center, but still require biological quality
        remaining["slot_select_score"] = (
            (1.0 - float(args.slot_weight)) * remaining["site_score"].astype(float)
            + float(args.slot_weight) * remaining["slot_proximity"].astype(float)
        )
        remaining = remaining.sort_values(
            ["slot_select_score", "site_score", "accessibility_score", "specificity_score", "robustness_score"],
            ascending=False,
        )

        chosen = None
        chosen_metrics = None
        for _, row in remaining.iterrows():
            if float(row["site_score"]) < float(args.min_site_score):
                continue
            feasible, metrics = candidate_is_feasible(row, selected_intervals, selected_seqs, args)
            if feasible:
                chosen = row
                chosen_metrics = metrics
                break

        if chosen is None:
            continue

        selected_site_ids.add(str(chosen["site_id"]))
        cand = (int(chosen["start"]), int(chosen["end"]))
        selected_intervals.append(cand)
        selected_seqs.append(str(chosen["sequence_ref"]))
        selected_centers.append(float(chosen["center"]))

        rowd = chosen.to_dict()
        rowd.update(chosen_metrics)
        cross_penalty = float(args.cross_penalty_weight) * float(chosen_metrics["cross_penalty_max"]) + float(args.cross_penalty_sum_weight) * float(chosen_metrics["cross_penalty_sum"])
        rowd["cross_penalty"] = round(cross_penalty, 6)
        rowd["slot_idx"] = slot_idx
        rowd["slot_center"] = round(slot_center, 3)
        rowd["slot_proximity"] = round(float(chosen["slot_proximity"]), 6)
        rowd["selection_phase"] = "slot_first_pass"
        rowd["reason_selected"] = (
            f"phase=slot_first_pass; slot_idx={slot_idx}; slot_center={slot_center:.1f}; "
            f"site_score={float(chosen['site_score']):.3f}; slot_proximity={float(chosen['slot_proximity']):.3f}; "
            f"gain_bp={chosen_metrics['gain_bp']}; nearest_gap_bp={chosen_metrics['nearest_gap_bp']}; "
            f"cross_penalty={cross_penalty:.3f}"
        )
        if args.assay_type == "rt_primer_25":
            rowd["rt_primer_seq"] = revcomp(str(chosen["sequence_ref"]))
        selected_rows.append(rowd)

    return selected_rows, selected_site_ids, selected_intervals, selected_seqs, selected_centers


def second_pass_gap_fill(eligible: pd.DataFrame, slot_centers: List[float], target_gap: float, genome_start: float, genome_end: float, selected_rows, selected_site_ids, selected_intervals, selected_seqs, selected_centers, args):
    while len(selected_rows) < int(args.max_sites):
        remaining = eligible[~eligible["site_id"].isin(selected_site_ids)].copy()
        if remaining.empty:
            break

        if selected_centers:
            remaining["nearest_selected_dist"] = remaining["center"].map(lambda c: min(abs(float(c) - s) for s in selected_centers))
        else:
            remaining["nearest_selected_dist"] = target_gap

        # favor candidates that help fill wide spaces; ideal nearest distance is around target_gap
        remaining["gap_fill_score"] = remaining["nearest_selected_dist"].map(
            lambda d: clamp01(1.0 - abs(float(d) - target_gap) / max(target_gap, 1.0))
        )
        # small bonus for helping extend current span early on
        current_largest_gap = largest_selected_gap(selected_centers, genome_start, genome_end) if selected_centers else (genome_end - genome_start)
        remaining["second_pass_score"] = (
            (1.0 - float(args.gap_fill_weight)) * remaining["site_score"].astype(float)
            + float(args.gap_fill_weight) * remaining["gap_fill_score"].astype(float)
        )
        remaining = remaining.sort_values(
            ["second_pass_score", "site_score", "gap_fill_score", "accessibility_score", "specificity_score"],
            ascending=False,
        )

        chosen = None
        chosen_metrics = None
        for _, row in remaining.iterrows():
            if float(row["site_score"]) < float(args.min_site_score):
                continue
            feasible, metrics = candidate_is_feasible(row, selected_intervals, selected_seqs, args)
            if feasible:
                chosen = row
                chosen_metrics = metrics
                break

        if chosen is None:
            break

        selected_site_ids.add(str(chosen["site_id"]))
        cand = (int(chosen["start"]), int(chosen["end"]))
        selected_intervals.append(cand)
        selected_seqs.append(str(chosen["sequence_ref"]))
        selected_centers.append(float(chosen["center"]))

        rowd = chosen.to_dict()
        rowd.update(chosen_metrics)
        cross_penalty = float(args.cross_penalty_weight) * float(chosen_metrics["cross_penalty_max"]) + float(args.cross_penalty_sum_weight) * float(chosen_metrics["cross_penalty_sum"])
        rowd["cross_penalty"] = round(cross_penalty, 6)
        rowd["slot_idx"] = int(chosen["nearest_slot_idx"])
        rowd["slot_center"] = round(float(chosen["nearest_slot_center"]), 3)
        rowd["slot_proximity"] = round(float(chosen["nearest_slot_proximity"]), 6)
        rowd["gap_fill_score"] = round(float(chosen["gap_fill_score"]), 6)
        rowd["selection_phase"] = "gap_fill_second_pass"
        rowd["reason_selected"] = (
            f"phase=gap_fill_second_pass; nearest_slot_idx={int(chosen['nearest_slot_idx'])}; "
            f"site_score={float(chosen['site_score']):.3f}; gap_fill_score={float(chosen['gap_fill_score']):.3f}; "
            f"gain_bp={chosen_metrics['gain_bp']}; nearest_gap_bp={chosen_metrics['nearest_gap_bp']}; "
            f"cross_penalty={cross_penalty:.3f}; largest_gap_before={current_largest_gap:.1f}"
        )
        if args.assay_type == "rt_primer_25":
            rowd["rt_primer_seq"] = revcomp(str(chosen["sequence_ref"]))
        selected_rows.append(rowd)

    return selected_rows


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description="Select evenly distributed panel")
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
    ap.add_argument("--min-accessibility-score", type=float, default=None)
    ap.add_argument("--min-access-3p-terminal", type=float, default=None)
    ap.add_argument("--min-robustness-score", type=float, default=None)
    ap.add_argument("--min-cov-3p-1mm", type=float, default=None)
    ap.add_argument("--min-cov-terminal-match", type=float, default=None)
    ap.add_argument("--max-human-hits", type=int, default=None)
    ap.add_argument("--max-human-anchored-hits", type=int, default=None)
    ap.add_argument("--max-human-nearperfect-hits", type=int, default=None)
    ap.add_argument("--max-virus-anchored-hits", type=int, default=None)
    ap.add_argument("--slot-weight", type=float, default=None, help="Weight on proximity to evenly spaced slots")
    ap.add_argument("--gap-fill-weight", type=float, default=None, help="Second-pass weight on filling larger spacing gaps")
    ap.add_argument("--target-gap", type=float, default=None, help="Override ideal gap between primer centers")
    ap.add_argument("--no-second-pass", action="store_true", help="Disable second-pass gap filling")
    return ap.parse_args()


def main() -> None:
    args = parse_args()
    # argparse defines --no-second-pass, not --second-pass
    # so initialize second_pass before filling defaults
    args.second_pass = False if bool(getattr(args, "no_second_pass", False)) else None
    resolve_assay_defaults(args)

    df = pd.read_csv(args.scored_tsv, sep="\t")
    if "candidate_keep" in df.columns:
        df = df[df["candidate_keep"].astype(bool)].copy()
    if "sequence_ref" not in df.columns and "seq" in df.columns:
        df["sequence_ref"] = df["seq"]
    if "sequence_ref" not in df.columns:
        raise SystemExit("Input TSV must include 'sequence_ref' or 'seq'.")
    for col in ["start", "end"]:
        if col not in df.columns:
            raise SystemExit("Input TSV must include start and end columns.")

    df["sequence_ref"] = df["sequence_ref"].astype(str).map(clean_seq)
    if "length" not in df.columns:
        df["length"] = df["end"].astype(int) - df["start"].astype(int) + 1
    df["thermo_score"] = df.apply(lambda r: thermo_score(r, args.assay_type), axis=1)
    df["synthesis_score"] = df.apply(lambda r: synthesis_score(r, args.assay_type), axis=1)
    df["site_score"] = df.apply(lambda r: base_site_score(r, args.assay_type), axis=1)
    df["prefilter_fail_reason"] = df.apply(lambda r: ";".join(rt25_prefilter_reasons(r, args)), axis=1)
    df["prefilter_pass"] = df["prefilter_fail_reason"].eq("")
    df["center"] = compute_centers(df)

    genome_start = int(df["start"].min())
    genome_end = int(df["end"].max())
    target_gap = float(args.target_gap) if args.target_gap is not None else ((genome_end - genome_start) / max(int(args.max_sites) - 1, 1))
    slot_centers = build_slot_centers(genome_start, genome_end, int(args.max_sites))

    # annotate nearest slot for diagnostics
    def nearest_slot(c):
        idx = min(range(len(slot_centers)), key=lambda i: abs(float(c) - slot_centers[i]))
        return idx + 1, slot_centers[idx], slot_proximity(float(c), slot_centers[idx], target_gap)

    nearest_info = df["center"].map(nearest_slot)
    df["nearest_slot_idx"] = nearest_info.map(lambda x: x[0])
    df["nearest_slot_center"] = nearest_info.map(lambda x: x[1])
    df["nearest_slot_proximity"] = nearest_info.map(lambda x: x[2])

    sort_cols = ["prefilter_pass", "site_score", "robustness_score", "accessibility_score", "specificity_score"]
    ranked = df.sort_values(sort_cols, ascending=False).reset_index(drop=True)
    eligible = ranked[ranked["prefilter_pass"].astype(bool)].copy()

    selected_rows, selected_site_ids, selected_intervals, selected_seqs, selected_centers = first_pass_slot_selection(
        eligible, slot_centers, target_gap, genome_start, genome_end, args
    )
    if bool(args.second_pass) and len(selected_rows) < int(args.max_sites):
        selected_rows = second_pass_gap_fill(
            eligible, slot_centers, target_gap, genome_start, genome_end,
            selected_rows, selected_site_ids, selected_intervals, selected_seqs, selected_centers, args
        )

    panel_df = pd.DataFrame(selected_rows)
    if not panel_df.empty:
        panel_df = panel_df.sort_values("center").reset_index(drop=True)
        panel_df.insert(0, "panel_rank", range(1, len(panel_df) + 1))
        panel_df["inter_primer_gap_from_prev"] = panel_df["start"].astype(int) - panel_df["end"].astype(int).shift(1).fillna(genome_start - 1).astype(int) - 1
        panel_df["genome_start"] = genome_start
        panel_df["genome_end"] = genome_end
        panel_df["target_gap"] = round(target_gap, 3)
    ranked = ranked.sort_values(["prefilter_pass", "nearest_slot_idx", "site_score"], ascending=[False, True, False]).reset_index(drop=True)

    Path(args.ranked_output).parent.mkdir(parents=True, exist_ok=True)
    Path(args.panel_output).parent.mkdir(parents=True, exist_ok=True)
    ranked.to_csv(args.ranked_output, sep="\t", index=False)
    panel_df.to_csv(args.panel_output, sep="\t", index=False)

    print(f"Assay type: {args.assay_type}")
    print(f"Genome span: {genome_start}-{genome_end}")
    print(f"Target gap: {target_gap:.2f}")
    print(f"Ranked candidates: {len(ranked)}")
    print(f"Prefilter-eligible candidates: {int(ranked['prefilter_pass'].sum())}")
    print(f"Selected panel size: {len(panel_df)}")
    print(f"Wrote ranked TSV -> {args.ranked_output}")
    print(f"Wrote panel TSV -> {args.panel_output}")


if __name__ == "__main__":
    main()
