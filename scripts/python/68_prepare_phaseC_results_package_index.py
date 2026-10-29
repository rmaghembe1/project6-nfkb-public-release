#!/usr/bin/env python3
"""Prepare Phase C results package index.

This checkpoint maps manuscript-ready Phase C result statements to their
supporting files, QC JSONs and audit decisions. It does not perform new analysis,
write manuscript prose, infer selection, validate function or authorise broad
biological interpretation.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path

MANUSCRIPT_GATE_JSON = Path("results/qc/phaseC_manuscript_readiness_gate_v1.json")
MANUSCRIPT_ALLOWED = Path(
    "11_results/phaseC/genotype_only_qc/manuscript_readiness_gate/"
    "phaseC_allowed_manuscript_statements_v1.tsv"
)

OUTDIR = Path("11_results/phaseC/genotype_only_qc/results_package_index")
OUTJSON = Path("results/qc/phaseC_results_package_index_v1.json")

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


for path in [MANUSCRIPT_GATE_JSON, MANUSCRIPT_ALLOWED]:
    if not path.exists():
        failures.append(f"Required input missing: {path}")

gate_payload = json.loads(MANUSCRIPT_GATE_JSON.read_text(encoding="utf-8")) if MANUSCRIPT_GATE_JSON.exists() else {}
allowed_rows = read_tsv(MANUSCRIPT_ALLOWED) if MANUSCRIPT_ALLOWED.exists() else []

if gate_payload.get("status") != "PASS":
    failures.append("Manuscript-readiness gate does not have PASS status.")
if gate_payload.get("selection_or_population_claim_authorised") != "NO":
    failures.append("Manuscript-readiness gate does not block selection/population claims.")
if len(allowed_rows) != 6:
    failures.append(f"Expected 6 allowed manuscript statements; observed {len(allowed_rows)}.")

package_rows = [
    {
        "result_package_id": "PHASEC_RESULT_001",
        "result_topic": "strict_genotype_only_qc",
        "allowed_statement_id": "ALLOWED_TECH_001",
        "primary_result_file": "results/qc/phaseC_strict_genotype_only_group_frequency_qc_v1.json",
        "supporting_tables_or_dirs": "11_results/phaseC/genotype_only_qc/strict_group_frequency_qc",
        "audit_decision": "15_logs_and_audit/phaseC/genotype_only_qc/strict_group_frequency_qc/PhaseC_Strict_Genotype_Only_Group_Frequency_QC_Decision_2026-06-02.md",
        "manuscript_use": "Technical QC statement only.",
        "blocked_claims": "No biological interpretation, selection inference or population-specific claim.",
    },
    {
        "result_package_id": "PHASEC_RESULT_002",
        "result_topic": "candidate_prioritisation",
        "allowed_statement_id": "ALLOWED_TECH_002",
        "primary_result_file": "results/qc/phaseC_candidate_prioritisation_applied_v1.json",
        "supporting_tables_or_dirs": "11_results/phaseC/genotype_only_qc/candidate_prioritisation_applied",
        "audit_decision": "15_logs_and_audit/phaseC/genotype_only_qc/candidate_prioritisation_applied/PhaseC_Candidate_Prioritisation_Applied_Decision_2026-06-02.md",
        "manuscript_use": "Predeclared triage and follow-up prioritisation statement.",
        "blocked_claims": "No biological validation or causal ranking.",
    },
    {
        "result_package_id": "PHASEC_RESULT_003",
        "result_topic": "annotation_supported_combined_evidence",
        "allowed_statement_id": "ALLOWED_TECH_003",
        "primary_result_file": "results/qc/phaseC_annotation_supported_combined_evidence_summary_v1.json",
        "supporting_tables_or_dirs": "11_results/phaseC/genotype_only_qc/annotation_supported_evidence_review/combined_summary",
        "audit_decision": "15_logs_and_audit/phaseC/genotype_only_qc/annotation_supported_evidence_review/combined_summary/PhaseC_Annotation_Supported_Combined_Evidence_Summary_Decision_2026-06-02.md",
        "manuscript_use": "Evidence-grade summary for four annotation-supported candidates.",
        "blocked_claims": "No new evidence review or broad biological interpretation.",
    },
    {
        "result_package_id": "PHASEC_RESULT_004",
        "result_topic": "TRAF3_evidence",
        "allowed_statement_id": "ALLOWED_EVID_001",
        "primary_result_file": "results/qc/phaseC_TRAF3_chr14_102875712_T_C_evidence_review_v1.json",
        "supporting_tables_or_dirs": "11_results/phaseC/genotype_only_qc/annotation_supported_evidence_review/TRAF3_chr14_102875712_T_C",
        "audit_decision": "15_logs_and_audit/phaseC/genotype_only_qc/annotation_supported_evidence_review/TRAF3_chr14_102875712_T_C/PhaseC_TRAF3_chr14_102875712_T_C_Evidence_Review_Decision_2026-06-02.md",
        "manuscript_use": "Bounded exact-variant support statement.",
        "blocked_claims": "No functional validation, causality, selection or ancestry-effect claim.",
    },
    {
        "result_package_id": "PHASEC_RESULT_005",
        "result_topic": "BIRC3_evidence",
        "allowed_statement_id": "ALLOWED_EVID_002",
        "primary_result_file": "results/qc/phaseC_BIRC3_chr11_102324736_G_T_evidence_review_v1.json",
        "supporting_tables_or_dirs": "11_results/phaseC/genotype_only_qc/annotation_supported_evidence_review/BIRC3_chr11_102324736_G_T",
        "audit_decision": "15_logs_and_audit/phaseC/genotype_only_qc/annotation_supported_evidence_review/BIRC3_chr11_102324736_G_T/PhaseC_BIRC3_chr11_102324736_G_T_Evidence_Review_Decision_2026-06-02.md",
        "manuscript_use": "Bounded local VEP identity/consequence statement.",
        "blocked_claims": "No external validation or functional validation claim.",
    },
    {
        "result_package_id": "PHASEC_RESULT_006",
        "result_topic": "TNFAIP3_paired_evidence",
        "allowed_statement_id": "ALLOWED_EVID_003",
        "primary_result_file": "results/qc/phaseC_TNFAIP3_paired_evidence_review_v1.json",
        "supporting_tables_or_dirs": "11_results/phaseC/genotype_only_qc/annotation_supported_evidence_review/TNFAIP3_paired_candidates",
        "audit_decision": "15_logs_and_audit/phaseC/genotype_only_qc/annotation_supported_evidence_review/TNFAIP3_paired_candidates/PhaseC_TNFAIP3_Paired_Evidence_Review_Decision_2026-06-02.md",
        "manuscript_use": "Bounded exact-variant identity plus paired/haplotype-level functional-context statement.",
        "blocked_claims": "No independent single-variant functional validation, project-level selection, or ancestry-effect claim.",
    },
    {
        "result_package_id": "PHASEC_RESULT_007",
        "result_topic": "manuscript_readiness_gate",
        "allowed_statement_id": "GATE_CONTROL",
        "primary_result_file": "results/qc/phaseC_manuscript_readiness_gate_v1.json",
        "supporting_tables_or_dirs": "11_results/phaseC/genotype_only_qc/manuscript_readiness_gate",
        "audit_decision": "15_logs_and_audit/phaseC/genotype_only_qc/manuscript_readiness_gate/PhaseC_Manuscript_Readiness_Gate_Decision_2026-06-02.md",
        "manuscript_use": "Controls allowed language, prohibited claims and required limitations.",
        "blocked_claims": "All statements exceeding the readiness gate.",
    },
]

missing_paths = []
for row in package_rows:
    for field in ["primary_result_file", "supporting_tables_or_dirs", "audit_decision"]:
        path = Path(row[field])
        if not path.exists():
            missing_paths.append(str(path))

if missing_paths:
    failures.extend(f"Missing package path: {path}" for path in missing_paths)

index_path = OUTDIR / "phaseC_results_package_index_v1.tsv"
write_rows(
    index_path,
    [
        "result_package_id",
        "result_topic",
        "allowed_statement_id",
        "primary_result_file",
        "supporting_tables_or_dirs",
        "audit_decision",
        "manuscript_use",
        "blocked_claims",
    ],
    package_rows,
)

summary_rows = [
    {"metric": "result_package_rows", "value": len(package_rows)},
    {"metric": "allowed_statement_rows_available", "value": len(allowed_rows)},
    {"metric": "missing_package_paths", "value": len(missing_paths)},
    {"metric": "new_analysis_performed", "value": "NO"},
    {"metric": "manuscript_prose_written", "value": "NO"},
    {"metric": "selection_or_population_claim_authorised", "value": "NO"},
    {"metric": "broad_biological_interpretation_authorised", "value": "NO"},
]
summary_path = OUTDIR / "phaseC_results_package_index_summary_v1.tsv"
write_rows(summary_path, ["metric", "value"], summary_rows)

payload = {
    "status": "PASS" if not failures else "FAIL",
    "operation": "Prepare Phase C results package index.",
    "result_package_rows": len(package_rows),
    "allowed_statement_rows_available": len(allowed_rows),
    "missing_package_paths": len(missing_paths),
    "new_analysis_performed": "NO",
    "manuscript_prose_written": "NO",
    "selection_or_population_claim_authorised": "NO",
    "broad_biological_interpretation_authorised": "NO",
    "outputs": {
        "results_package_index": str(index_path),
        "summary": str(summary_path),
    },
    "failures": failures,
    "boundary": (
        "This checkpoint maps Phase C manuscript-ready result statements to supporting files only. "
        "It does not perform new analysis, write manuscript prose, validate function, infer selection, "
        "make population-differentiation claims, assign ancestry effects or authorise broad biological interpretation."
    ),
}

OUTJSON.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
print(json.dumps(payload, indent=2))

if payload["status"] != "PASS":
    raise SystemExit(1)
