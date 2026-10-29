#!/usr/bin/env python3
"""Audit Phase C top candidate prioritisation table.

This checkpoint audits the applied technical prioritisation output for internal
consistency, interpretation boundaries, rank/score logic and follow-up labels.
It does not add biological interpretation.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path
from collections import Counter

APPLIED = Path(
    "11_results/phaseC/genotype_only_qc/candidate_prioritisation_applied/"
    "phaseC_candidate_prioritisation_applied_v1.tsv"
)
FOLLOW_UP = Path(
    "11_results/phaseC/genotype_only_qc/candidate_prioritisation_applied/"
    "phaseC_candidate_prioritisation_follow_up_candidates_v1.tsv"
)
CLAIM_COUNTS = Path(
    "11_results/phaseC/genotype_only_qc/candidate_prioritisation_applied/"
    "phaseC_candidate_prioritisation_claim_label_counts_v1.tsv"
)
SCORE_DIST = Path(
    "11_results/phaseC/genotype_only_qc/candidate_prioritisation_applied/"
    "phaseC_candidate_prioritisation_score_distribution_v1.tsv"
)

OUTDIR = Path("11_results/phaseC/genotype_only_qc/top_candidate_audit")
OUTJSON = Path("results/qc/phaseC_top_candidate_table_audit_v1.json")

OUTDIR.mkdir(parents=True, exist_ok=True)
OUTJSON.parent.mkdir(parents=True, exist_ok=True)

failures: list[str] = []


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_rows(path: Path, fieldnames: list[str], rows: list[dict[str, object]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row.get(key, ".") for key in fieldnames})


for path in [APPLIED, FOLLOW_UP, CLAIM_COUNTS, SCORE_DIST]:
    if not path.exists():
        failures.append(f"Required input missing: {path}")

applied = read_tsv(APPLIED) if APPLIED.exists() else []
follow_up = read_tsv(FOLLOW_UP) if FOLLOW_UP.exists() else []

if len(applied) != 500:
    failures.append(f"Expected 500 applied-prioritisation rows; observed {len(applied)}.")
if len(follow_up) != 173:
    failures.append(f"Expected 173 follow-up rows; observed {len(follow_up)}.")

# Rank checks.
ranks = [int(row["technical_priority_rank"]) for row in applied]
if ranks != list(range(1, len(applied) + 1)):
    failures.append("Technical priority ranks are not consecutive 1..500.")

# Score sort checks.
score_tuples = [
    (
        int(row["technical_priority_score"]),
        int(row["canonical_consequence_score"]),
        int(row["frequency_delta_score"]),
        row["gene"],
        row["allele_key"],
    )
    for row in applied
]
expected_sorted = sorted(score_tuples, key=lambda x: (-x[0], -x[1], -x[2], x[3], x[4]))
if score_tuples != expected_sorted:
    failures.append("Applied prioritisation table is not sorted by the declared score/rank logic.")

# Component sum checks and boundary checks.
component_mismatch_rows = []
boundary_failure_rows = []
label_logic_failure_rows = []

for row in applied:
    total = int(row["technical_priority_score"])
    expected_total = (
        int(row["frequency_delta_score"])
        + int(row["canonical_consequence_score"])
        + int(row["any_transcript_review_score"])
        + int(row["gene_context_score"])
        + int(row["robustness_penalty"])
    )
    if total != expected_total:
        component_mismatch_rows.append(
            {
                "allele_key": row["allele_key"],
                "technical_priority_score": total,
                "expected_component_sum": expected_total,
            }
        )

    if "not evidence of causality" not in row["interpretation_boundary"]:
        boundary_failure_rows.append(
            {
                "allele_key": row["allele_key"],
                "interpretation_boundary": row["interpretation_boundary"],
            }
        )

    label = row["claim_boundary_label"]
    canonical = int(row["canonical_consequence_score"])
    freq = int(row["frequency_delta_score"])
    any_score = int(row["any_transcript_review_score"])

    expected_label = "technical_candidate_only"
    if canonical >= 3:
        expected_label = "annotation_supported_candidate"
    elif freq >= 2:
        expected_label = "frequency_variation_candidate"
    elif any_score > 0:
        expected_label = "follow_up_candidate"

    if label != expected_label:
        label_logic_failure_rows.append(
            {
                "allele_key": row["allele_key"],
                "observed_label": label,
                "expected_label": expected_label,
                "canonical_consequence_score": canonical,
                "frequency_delta_score": freq,
                "any_transcript_review_score": any_score,
            }
        )

if component_mismatch_rows:
    failures.append(f"Component-sum mismatches observed: {len(component_mismatch_rows)}.")
if boundary_failure_rows:
    failures.append(f"Interpretation-boundary failures observed: {len(boundary_failure_rows)}.")
if label_logic_failure_rows:
    failures.append(f"Claim-label logic failures observed: {len(label_logic_failure_rows)}.")

# Follow-up table consistency.
follow_up_ids = {row["allele_key"] for row in follow_up}
expected_follow_up_ids = {
    row["allele_key"] for row in applied
    if row["claim_boundary_label"] != "technical_candidate_only"
}
if follow_up_ids != expected_follow_up_ids:
    failures.append(
        f"Follow-up table ID mismatch: expected {len(expected_follow_up_ids)}, observed {len(follow_up_ids)}."
    )

# Create compact audit extracts.
top_n = 25
top_rows = applied[:top_n]
top_path = OUTDIR / "phaseC_top25_candidate_audit_table_v1.tsv"
write_rows(top_path, list(applied[0].keys()), top_rows)

annotation_supported_rows = [
    row for row in applied if row["claim_boundary_label"] == "annotation_supported_candidate"
]
annotation_path = OUTDIR / "phaseC_annotation_supported_candidate_audit_table_v1.tsv"
write_rows(annotation_path, list(applied[0].keys()), annotation_supported_rows)

observed_only_rows = [
    row for row in applied if row["observed_only_or_group_restricted_flag"] == "YES"
]
observed_only_path = OUTDIR / "phaseC_observed_only_or_group_restricted_candidate_audit_table_v1.tsv"
write_rows(observed_only_path, list(applied[0].keys()), observed_only_rows)

# Audit summary counts.
label_counts = Counter(row["claim_boundary_label"] for row in applied)
score_counts = Counter(int(row["technical_priority_score"]) for row in applied)
observed_flag_counts = Counter(row["observed_only_or_group_restricted_flag"] for row in applied)
module_counts = Counter(row["module"] for row in applied)

summary_rows = [
    {"metric": "applied_candidate_rows", "value": len(applied)},
    {"metric": "follow_up_candidate_rows", "value": len(follow_up)},
    {"metric": "top_audit_rows", "value": len(top_rows)},
    {"metric": "annotation_supported_candidate_rows", "value": len(annotation_supported_rows)},
    {"metric": "observed_only_or_group_restricted_rows", "value": len(observed_only_rows)},
    {"metric": "component_sum_mismatch_rows", "value": len(component_mismatch_rows)},
    {"metric": "boundary_failure_rows", "value": len(boundary_failure_rows)},
    {"metric": "label_logic_failure_rows", "value": len(label_logic_failure_rows)},
    {"metric": "max_technical_priority_score", "value": max(score_counts) if score_counts else 0},
]
summary_path = OUTDIR / "phaseC_top_candidate_audit_summary_v1.tsv"
write_rows(summary_path, ["metric", "value"], summary_rows)

def counter_rows(counter: Counter, key_name: str, value_name: str) -> list[dict[str, object]]:
    return [
        {key_name: key, value_name: value}
        for key, value in sorted(counter.items(), key=lambda x: (str(x[0])))
    ]

label_path = OUTDIR / "phaseC_top_candidate_audit_claim_label_counts_v1.tsv"
write_rows(label_path, ["claim_boundary_label", "candidate_count"], counter_rows(label_counts, "claim_boundary_label", "candidate_count"))

score_path = OUTDIR / "phaseC_top_candidate_audit_score_counts_v1.tsv"
write_rows(score_path, ["technical_priority_score", "candidate_count"], counter_rows(score_counts, "technical_priority_score", "candidate_count"))

flag_path = OUTDIR / "phaseC_top_candidate_audit_observed_flag_counts_v1.tsv"
write_rows(flag_path, ["observed_only_or_group_restricted_flag", "candidate_count"], counter_rows(observed_flag_counts, "observed_only_or_group_restricted_flag", "candidate_count"))

module_path = OUTDIR / "phaseC_top_candidate_audit_module_counts_v1.tsv"
write_rows(module_path, ["module", "candidate_count"], counter_rows(module_counts, "module", "candidate_count"))

# Write empty failure tables with headers for audit completeness.
write_rows(
    OUTDIR / "phaseC_top_candidate_audit_component_sum_failures_v1.tsv",
    ["allele_key", "technical_priority_score", "expected_component_sum"],
    component_mismatch_rows,
)
write_rows(
    OUTDIR / "phaseC_top_candidate_audit_boundary_failures_v1.tsv",
    ["allele_key", "interpretation_boundary"],
    boundary_failure_rows,
)
write_rows(
    OUTDIR / "phaseC_top_candidate_audit_label_logic_failures_v1.tsv",
    [
        "allele_key",
        "observed_label",
        "expected_label",
        "canonical_consequence_score",
        "frequency_delta_score",
        "any_transcript_review_score",
    ],
    label_logic_failure_rows,
)

payload = {
    "status": "PASS" if not failures else "FAIL",
    "operation": "Audit Phase C applied top-candidate prioritisation table.",
    "applied_candidate_rows": len(applied),
    "follow_up_candidate_rows": len(follow_up),
    "top_audit_rows": len(top_rows),
    "annotation_supported_candidate_rows": len(annotation_supported_rows),
    "observed_only_or_group_restricted_rows": len(observed_only_rows),
    "component_sum_mismatch_rows": len(component_mismatch_rows),
    "boundary_failure_rows": len(boundary_failure_rows),
    "label_logic_failure_rows": len(label_logic_failure_rows),
    "claim_label_counts": dict(label_counts),
    "observed_flag_counts": dict(observed_flag_counts),
    "max_technical_priority_score": max(score_counts) if score_counts else 0,
    "outputs": {
        "audit_summary": str(summary_path),
        "top25_audit_table": str(top_path),
        "annotation_supported_candidate_audit_table": str(annotation_path),
        "observed_only_or_group_restricted_candidate_audit_table": str(observed_only_path),
        "claim_label_counts": str(label_path),
        "score_counts": str(score_path),
        "observed_flag_counts": str(flag_path),
        "module_counts": str(module_path),
    },
    "failures": failures,
    "boundary": (
        "This checkpoint audits technical prioritisation consistency only. "
        "It does not validate function, infer selection, make population-differentiation claims, "
        "or authorise biological interpretation."
    ),
}

OUTJSON.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
print(json.dumps(payload, indent=2))

if payload["status"] != "PASS":
    raise SystemExit(1)
