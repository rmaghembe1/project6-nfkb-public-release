#!/usr/bin/env python3
"""Profile strict genotype-only frequency-variation candidates.

This is a descriptive QC/prioritisation-support step only.
It does not authorise population comparison, selection inference, functional
prioritisation or biological interpretation.
"""

from __future__ import annotations

import csv
import json
from collections import Counter, defaultdict
from pathlib import Path

GROUP_DIR = Path("11_results/phaseC/genotype_only_qc/strict_group_frequency_qc")
AUTO_BROAD = GROUP_DIR / "phaseC_strict_autosomal_broad_group_frequency_qc_v1.tsv"
CHRX_BROAD = GROUP_DIR / "phaseC_strict_IKBKG_chrX_broad_group_frequency_qc_v1.tsv"
AUTO_POP = GROUP_DIR / "phaseC_strict_autosomal_population_group_frequency_qc_v1.tsv"
CHRX_POP = GROUP_DIR / "phaseC_strict_IKBKG_chrX_population_group_frequency_qc_v1.tsv"

OUTDIR = Path("11_results/phaseC/genotype_only_qc/frequency_variation_candidates")
OUTJSON = Path("results/qc/phaseC_strict_frequency_variation_candidate_profile_v1.json")

EXPECTED_AUTO_BROAD_ROWS = 276066
EXPECTED_CHRX_BROAD_ROWS = 2681
EXPECTED_AUTO_POP_ROWS = 3155040
EXPECTED_CHRX_POP_ROWS = 30640

TOP_N = 250
MIN_AN_FOR_GROUP_AF = 20

OUTDIR.mkdir(parents=True, exist_ok=True)
OUTJSON.parent.mkdir(parents=True, exist_ok=True)

failures: list[str] = []


def read_rows(path: Path):
    with path.open(newline="", encoding="utf-8") as handle:
        yield from csv.DictReader(handle, delimiter="\t")


def safe_float(value: str) -> float:
    return float(value) if value not in {"", "NA", "nan"} else 0.0


def summarise_broad(path: Path, layer: str, expected_rows: int) -> tuple[list[dict[str, object]], dict[str, int]]:
    by_variant: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    n_rows = 0

    for row in read_rows(path):
        n_rows += 1
        by_variant[(row["gene_symbol"], row["allele_key"])].append(row)

    if n_rows != expected_rows:
        failures.append(f"{layer}: expected {expected_rows} broad rows; observed {n_rows}.")

    summaries = []
    for (gene, key), rows in by_variant.items():
        groups_with_ac = []
        groups_with_an = []
        af_values = []
        missing_by_group = {}

        for row in rows:
            group = row["group"]
            an = int(row["AN"])
            ac = int(row["AC"])
            af = safe_float(row["AF"])
            missing = int(row["missing_GT_cells"])

            if an >= MIN_AN_FOR_GROUP_AF:
                groups_with_an.append(group)
                af_values.append((group, af))
            if ac > 0:
                groups_with_ac.append(group)
            if missing > 0:
                missing_by_group[group] = missing

        if af_values:
            min_group, min_af = min(af_values, key=lambda x: x[1])
            max_group, max_af = max(af_values, key=lambda x: x[1])
            max_delta = max_af - min_af
        else:
            min_group = max_group = ""
            min_af = max_af = max_delta = 0.0

        summaries.append({
            "analysis_layer": layer,
            "gene_symbol": gene,
            "allele_key": key,
            "broad_groups_observed": len(rows),
            "broad_groups_with_AN_ge_20": len(groups_with_an),
            "broad_groups_with_AC_gt_0": len(set(groups_with_ac)),
            "observed_only_one_broad_group": "YES" if len(set(groups_with_ac)) == 1 else "NO",
            "max_AF_group": max_group,
            "max_AF": max_af,
            "min_AF_group": min_group,
            "min_AF": min_af,
            "max_minus_min_AF": max_delta,
            "groups_with_missing_GT": len(missing_by_group),
        })

    stats = {
        "input_rows": n_rows,
        "variant_records": len(summaries),
        "observed_only_one_broad_group": sum(r["observed_only_one_broad_group"] == "YES" for r in summaries),
    }
    return summaries, stats


auto_broad_summary, auto_stats = summarise_broad(AUTO_BROAD, "autosomal", EXPECTED_AUTO_BROAD_ROWS)
chrx_broad_summary, chrx_stats = summarise_broad(CHRX_BROAD, "chrX_sex_aware", EXPECTED_CHRX_BROAD_ROWS)

broad_summary = auto_broad_summary + chrx_broad_summary

broad_summary_path = OUTDIR / "phaseC_strict_broad_frequency_variation_summary_v1.tsv"
with broad_summary_path.open("w", newline="", encoding="utf-8") as handle:
    writer = csv.DictWriter(handle, fieldnames=list(broad_summary[0]), delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(broad_summary)

# Top descriptive candidates by broad max-minus-min AF.
top_delta = sorted(
    broad_summary,
    key=lambda r: (float(r["max_minus_min_AF"]), int(r["broad_groups_with_AC_gt_0"])),
    reverse=True,
)[:TOP_N]

top_delta_path = OUTDIR / "phaseC_strict_top_broad_AF_delta_candidates_v1.tsv"
with top_delta_path.open("w", newline="", encoding="utf-8") as handle:
    writer = csv.DictWriter(handle, fieldnames=list(top_delta[0]), delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(top_delta)

observed_one = [
    r for r in broad_summary
    if r["observed_only_one_broad_group"] == "YES"
]
observed_one = sorted(observed_one, key=lambda r: float(r["max_AF"]), reverse=True)

one_group_path = OUTDIR / "phaseC_strict_observed_only_one_broad_group_candidates_v1.tsv"
with one_group_path.open("w", newline="", encoding="utf-8") as handle:
    writer = csv.DictWriter(handle, fieldnames=list(broad_summary[0]), delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(observed_one[:TOP_N])

# Fine-population primary-reportable group count summary only; do not track full population table.
def summarise_population_table(path: Path, layer: str, expected_rows: int) -> dict[str, int]:
    n_rows = 0
    rows_meeting_min = 0
    rows_below_min = 0
    variants_seen = set()
    groups = set()
    groups_below_min = set()
    groups_meeting_min = set()

    for row in read_rows(path):
        n_rows += 1
        variants_seen.add((row["gene_symbol"], row["allele_key"]))
        group = row["group"]
        groups.add(group)
        if row["meets_min_n_20"] == "YES":
            rows_meeting_min += 1
            groups_meeting_min.add(group)
        else:
            rows_below_min += 1
            groups_below_min.add(group)

    if n_rows != expected_rows:
        failures.append(f"{layer}: expected {expected_rows} population rows; observed {n_rows}.")

    return {
        "input_rows": n_rows,
        "variant_records": len(variants_seen),
        "groups_total": len(groups),
        "groups_meeting_min_n_20": len(groups_meeting_min),
        "groups_below_min_n_20": len(groups_below_min),
        "rows_meeting_min_n_20": rows_meeting_min,
        "rows_below_min_n_20": rows_below_min,
    }


auto_pop_stats = summarise_population_table(AUTO_POP, "autosomal_population", EXPECTED_AUTO_POP_ROWS)
chrx_pop_stats = summarise_population_table(CHRX_POP, "chrX_population", EXPECTED_CHRX_POP_ROWS)

summary_rows = [
    {"analysis_layer": "autosomal_broad", **auto_stats},
    {"analysis_layer": "chrX_sex_aware_broad", **chrx_stats},
    {"analysis_layer": "autosomal_population", **auto_pop_stats},
    {"analysis_layer": "chrX_sex_aware_population", **chrx_pop_stats},
]

summary_path = OUTDIR / "phaseC_strict_frequency_variation_candidate_profile_summary_v1.tsv"
with summary_path.open("w", newline="", encoding="utf-8") as handle:
    fieldnames = sorted({k for row in summary_rows for k in row})
    writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(summary_rows)

payload = {
    "status": "PASS" if not failures else "FAIL",
    "operation": "Profile strict genotype-only descriptive frequency-variation candidates.",
    "top_n": TOP_N,
    "minimum_AN_for_group_AF": MIN_AN_FOR_GROUP_AF,
    "summary": summary_rows,
    "outputs": {
        "broad_summary": str(broad_summary_path),
        "top_broad_AF_delta_candidates": str(top_delta_path),
        "observed_only_one_broad_group_candidates": str(one_group_path),
        "profile_summary": str(summary_path),
    },
    "failures": failures,
    "boundary": (
        "This is a descriptive candidate-profiling/QC step only. Maximum AF differences and "
        "observed-only group labels are not interpreted as population differentiation, selection, "
        "function or biological importance."
    ),
}

OUTJSON.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
print(json.dumps(payload, indent=2))

if payload["status"] != "PASS":
    raise SystemExit(1)
