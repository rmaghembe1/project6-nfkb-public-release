#!/usr/bin/env python3
"""Apply predeclared Phase C candidate-prioritisation framework.

This checkpoint applies the frozen prioritisation framework to strict Phase C
candidate tables. It produces technical prioritisation tiers only. It does not
claim causality, selection, population differentiation, ancestry effects or
functional validation.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path

FREQ_TOP = Path(
    "11_results/phaseC/genotype_only_qc/frequency_variation_candidates/"
    "phaseC_strict_top_broad_AF_delta_candidates_v1.tsv"
)
FREQ_OBS_ONLY = Path(
    "11_results/phaseC/genotype_only_qc/frequency_variation_candidates/"
    "phaseC_strict_observed_only_one_broad_group_candidates_v1.tsv"
)
VEP_SUMMARY = Path(
    "11_results/phaseC/genotype_only_qc/vep_consequence_landscape/"
    "phaseC_vep_landscape_per_candidate_technical_summary_v1.tsv"
)
MODULE_CONTEXT = Path(
    "11_results/phaseC/genotype_only_qc/gene_module_context_map/"
    "phaseC_candidate_gene_module_context_v1.tsv"
)
FRAMEWORK_JSON = Path("results/qc/phaseC_candidate_prioritisation_framework_design_v1.json")
MODULE_JSON = Path("results/qc/phaseC_gene_module_context_map_freeze_v1.json")

OUTDIR = Path("11_results/phaseC/genotype_only_qc/candidate_prioritisation_applied")
OUTJSON = Path("results/qc/phaseC_candidate_prioritisation_applied_v1.json")

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


required = [FREQ_TOP, FREQ_OBS_ONLY, VEP_SUMMARY, MODULE_CONTEXT, FRAMEWORK_JSON, MODULE_JSON]
for path in required:
    if not path.exists():
        failures.append(f"Required input missing: {path}")

freq_rows = read_tsv(FREQ_TOP) if FREQ_TOP.exists() else []
obs_only_rows = read_tsv(FREQ_OBS_ONLY) if FREQ_OBS_ONLY.exists() else []
vep_rows = read_tsv(VEP_SUMMARY) if VEP_SUMMARY.exists() else []
module_rows = read_tsv(MODULE_CONTEXT) if MODULE_CONTEXT.exists() else []

# Robust allele-key detection across upstream tables.
def get_key(row: dict[str, str]) -> str:
    for key in ["allele_key", "variant_id", "variant_key", "Uploaded_variation", "uploaded_variation"]:
        if row.get(key):
            return row[key]
    # Fallback if separated VCF-like columns are present.
    for chrom_key in ["chrom", "CHROM", "#CHROM"]:
        if chrom_key in row and "pos" in {k.lower(): k for k in row}:
            lower_map = {k.lower(): k for k in row}
            pos = row[lower_map["pos"]]
            ref = row.get("ref", row.get("REF", ""))
            alt = row.get("alt", row.get("ALT", ""))
            return f"{row[chrom_key]}:{pos}:{ref}:{alt}"
    return ""


freq_by_key = {get_key(row): row for row in freq_rows if get_key(row)}
obs_only_keys = {get_key(row) for row in obs_only_rows if get_key(row)}
vep_by_key = {row["allele_key"]: row for row in vep_rows}
module_by_key = {row["allele_key"]: row for row in module_rows}

candidate_ids = sorted(set(vep_by_key) | set(module_by_key))
if len(candidate_ids) != 500:
    failures.append(f"Expected 500 candidate IDs from VEP/module context; observed {len(candidate_ids)}.")
if set(vep_by_key) != set(module_by_key):
    failures.append("VEP consequence summary and module context candidate IDs do not match exactly.")

# AF-delta field detection.
def detect_float(row: dict[str, str], candidates: list[str]) -> float | None:
    for col in candidates:
        if col in row and row[col] not in {"", "."}:
            try:
                return float(row[col])
            except ValueError:
                continue
    return None


delta_candidates = [
    "broad_AF_delta",
    "AF_delta",
    "max_AF_delta",
    "max_minus_min_AF",
    "af_delta",
    "delta_AF",
    "max_af_delta",
]

# Determine empirical tertile-like thresholds using detected AF deltas in available top table.
delta_values = []
for row in freq_by_key.values():
    value = detect_float(row, delta_candidates)
    if value is not None:
        delta_values.append(value)

delta_values_sorted = sorted(delta_values)
if not delta_values_sorted:
    failures.append("No AF-delta numeric field detected in top frequency-variation candidate table.")
    low_cut = high_cut = None
else:
    # Conservative ordinal bins over observed top-candidate distribution.
    low_cut = delta_values_sorted[len(delta_values_sorted) // 3]
    high_cut = delta_values_sorted[(2 * len(delta_values_sorted)) // 3]


def frequency_delta_score(allele_key: str) -> tuple[int, str, str]:
    row = freq_by_key.get(allele_key)
    if not row:
        return 0, ".", "not_in_top_broad_AF_delta_table"
    value = detect_float(row, delta_candidates)
    if value is None or low_cut is None or high_cut is None:
        return 1, ".", "in_top_table_delta_column_unresolved"
    if value >= high_cut:
        return 3, value, "top_AF_delta_bin"
    if value >= low_cut:
        return 2, value, "middle_AF_delta_bin"
    return 1, value, "lower_AF_delta_bin"


def canonical_consequence_score(row: dict[str, str]) -> tuple[int, str]:
    impact_union = row.get("canonical_like_impact_union", ".")
    impacts = set(impact_union.split(",")) if impact_union and impact_union != "." else set()
    if "HIGH" in impacts:
        return 4, "HIGH"
    if "MODERATE" in impacts:
        return 3, "MODERATE"
    if "LOW" in impacts:
        return 1, "LOW"
    return 0, "MODIFIER_or_unresolved"


def any_transcript_review_score(row: dict[str, str], canonical_score: int) -> tuple[int, str]:
    any_flag = row.get("has_moderate_or_high_impact_any_transcript", "NO")
    canonical_flag = row.get("has_moderate_or_high_impact_canonical_like", "NO")
    if any_flag == "YES" and canonical_flag != "YES":
        return 1, "noncanonical_moderate_or_high_review_flag"
    return 0, "no_extra_any_transcript_review_score"


def claim_label(freq_score: int, canon_score: int, any_score: int) -> str:
    if canon_score >= 3:
        return "annotation_supported_candidate"
    if freq_score >= 2:
        return "frequency_variation_candidate"
    if any_score > 0:
        return "follow_up_candidate"
    return "technical_candidate_only"


prioritised_rows = []
for allele_key in candidate_ids:
    vep = vep_by_key[allele_key]
    module = module_by_key[allele_key]

    freq_score, af_delta_value, freq_note = frequency_delta_score(allele_key)
    canon_score, canon_impact_label = canonical_consequence_score(vep)
    any_score, any_note = any_transcript_review_score(vep, canon_score)
    gene_context_score = 1 if module.get("context_available") == "YES" else 0
    robustness_penalty = 0

    total_score = freq_score + canon_score + any_score + gene_context_score + robustness_penalty
    label = claim_label(freq_score, canon_score, any_score)

    prioritised_rows.append(
        {
            "allele_key": allele_key,
            "gene": vep.get("gene", module.get("gene", ".")),
            "module": module.get("module", "."),
            "submodule": module.get("submodule", "."),
            "frequency_delta_score": freq_score,
            "frequency_delta_value": af_delta_value,
            "frequency_delta_note": freq_note,
            "observed_only_or_group_restricted_flag": "YES" if allele_key in obs_only_keys else "NO",
            "canonical_consequence_score": canon_score,
            "canonical_impact_label": canon_impact_label,
            "any_transcript_review_score": any_score,
            "any_transcript_review_note": any_note,
            "gene_context_score": gene_context_score,
            "robustness_penalty": robustness_penalty,
            "technical_priority_score": total_score,
            "claim_boundary_label": label,
            "has_moderate_or_high_impact_canonical_like": vep.get("has_moderate_or_high_impact_canonical_like", "."),
            "has_moderate_or_high_impact_any_transcript": vep.get("has_moderate_or_high_impact_any_transcript", "."),
            "canonical_like_consequence_union": vep.get("canonical_like_consequence_union", "."),
            "canonical_like_impact_union": vep.get("canonical_like_impact_union", "."),
            "interpretation_boundary": "Technical prioritisation only; not evidence of causality, selection or functional validation.",
        }
    )

prioritised_rows.sort(
    key=lambda row: (
        -int(row["technical_priority_score"]),
        -int(row["canonical_consequence_score"]),
        -int(row["frequency_delta_score"]),
        row["gene"],
        row["allele_key"],
    )
)

for idx, row in enumerate(prioritised_rows, start=1):
    row["technical_priority_rank"] = idx

fieldnames = [
    "technical_priority_rank",
    "allele_key",
    "gene",
    "module",
    "submodule",
    "frequency_delta_score",
    "frequency_delta_value",
    "frequency_delta_note",
    "observed_only_or_group_restricted_flag",
    "canonical_consequence_score",
    "canonical_impact_label",
    "any_transcript_review_score",
    "any_transcript_review_note",
    "gene_context_score",
    "robustness_penalty",
    "technical_priority_score",
    "claim_boundary_label",
    "has_moderate_or_high_impact_canonical_like",
    "has_moderate_or_high_impact_any_transcript",
    "canonical_like_consequence_union",
    "canonical_like_impact_union",
    "interpretation_boundary",
]

priority_path = OUTDIR / "phaseC_candidate_prioritisation_applied_v1.tsv"
write_rows(priority_path, fieldnames, prioritised_rows)

top_review_rows = [row for row in prioritised_rows if row["claim_boundary_label"] != "technical_candidate_only"]
top_review_path = OUTDIR / "phaseC_candidate_prioritisation_follow_up_candidates_v1.tsv"
write_rows(top_review_path, fieldnames, top_review_rows)

# Summary tables.
label_counts: dict[str, int] = {}
module_label_counts: dict[tuple[str, str], int] = {}
score_counts: dict[int, int] = {}

for row in prioritised_rows:
    label_counts[row["claim_boundary_label"]] = label_counts.get(row["claim_boundary_label"], 0) + 1
    module_label_counts[(row["module"], row["claim_boundary_label"])] = module_label_counts.get((row["module"], row["claim_boundary_label"]), 0) + 1
    score_counts[int(row["technical_priority_score"])] = score_counts.get(int(row["technical_priority_score"]), 0) + 1

label_rows = [
    {"claim_boundary_label": label, "candidate_count": count}
    for label, count in sorted(label_counts.items())
]
label_path = OUTDIR / "phaseC_candidate_prioritisation_claim_label_counts_v1.tsv"
write_rows(label_path, ["claim_boundary_label", "candidate_count"], label_rows)

module_label_rows = [
    {"module": module, "claim_boundary_label": label, "candidate_count": count}
    for (module, label), count in sorted(module_label_counts.items())
]
module_label_path = OUTDIR / "phaseC_candidate_prioritisation_module_by_claim_label_counts_v1.tsv"
write_rows(module_label_path, ["module", "claim_boundary_label", "candidate_count"], module_label_rows)

score_rows = [
    {"technical_priority_score": score, "candidate_count": count}
    for score, count in sorted(score_counts.items(), reverse=True)
]
score_path = OUTDIR / "phaseC_candidate_prioritisation_score_distribution_v1.tsv"
write_rows(score_path, ["technical_priority_score", "candidate_count"], score_rows)

payload = {
    "status": "PASS" if not failures else "FAIL",
    "operation": "Apply predeclared Phase C candidate-prioritisation framework.",
    "candidate_records_scored": len(prioritised_rows),
    "follow_up_candidate_rows": len(top_review_rows),
    "claim_label_counts": label_counts,
    "max_technical_priority_score": max((int(row["technical_priority_score"]) for row in prioritised_rows), default=0),
    "frequency_delta_numeric_values_detected": len(delta_values),
    "frequency_delta_low_cut": low_cut,
    "frequency_delta_high_cut": high_cut,
    "outputs": {
        "prioritised_candidates": str(priority_path),
        "follow_up_candidates": str(top_review_path),
        "claim_label_counts": str(label_path),
        "module_by_claim_label_counts": str(module_label_path),
        "score_distribution": str(score_path),
    },
    "failures": failures,
    "boundary": (
        "This checkpoint applies predeclared technical prioritisation scores only. "
        "It does not validate function, infer selection, make population-differentiation claims, "
        "or authorise biological interpretation."
    ),
}

OUTJSON.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
print(json.dumps(payload, indent=2))

if payload["status"] != "PASS":
    raise SystemExit(1)
