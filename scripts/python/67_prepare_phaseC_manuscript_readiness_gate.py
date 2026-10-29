#!/usr/bin/env python3
"""Prepare Phase C interpretation-boundary and manuscript-readiness gate.

This checkpoint converts completed technical/evidence-review outputs into a
controlled manuscript-readiness gate. It defines what can and cannot be said in
results/manuscript language. It does not write manuscript Results prose, does
not perform new evidence review, does not validate function, and does not infer
selection or population effects.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path

COMBINED_JSON = Path("results/qc/phaseC_annotation_supported_combined_evidence_summary_v1.json")
COMBINED_TABLE = Path(
    "11_results/phaseC/genotype_only_qc/annotation_supported_evidence_review/"
    "combined_summary/phaseC_annotation_supported_combined_evidence_review_status_v1.tsv"
)
PRIORITISATION_JSON = Path("results/qc/phaseC_candidate_prioritisation_applied_v1.json")
TOP_AUDIT_JSON = Path("results/qc/phaseC_top_candidate_table_audit_v1.json")
GROUP_FREQ_JSON = Path("results/qc/phaseC_strict_genotype_only_group_frequency_qc_v1.json")

OUTDIR = Path("11_results/phaseC/genotype_only_qc/manuscript_readiness_gate")
OUTJSON = Path("results/qc/phaseC_manuscript_readiness_gate_v1.json")

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


for path in [COMBINED_JSON, COMBINED_TABLE, PRIORITISATION_JSON, TOP_AUDIT_JSON, GROUP_FREQ_JSON]:
    if not path.exists():
        failures.append(f"Required input missing: {path}")

combined_payload = json.loads(COMBINED_JSON.read_text(encoding="utf-8")) if COMBINED_JSON.exists() else {}
combined_rows = read_tsv(COMBINED_TABLE) if COMBINED_TABLE.exists() else []

if combined_payload.get("status") != "PASS":
    failures.append("Combined annotation-supported evidence summary does not have PASS status.")
if len(combined_rows) != 4:
    failures.append(f"Expected 4 combined annotation-supported candidates; observed {len(combined_rows)}.")
if combined_payload.get("selection_or_population_claim_authorised") != "NO":
    failures.append("Combined summary does not explicitly block selection/population claims.")
if combined_payload.get("broad_biological_interpretation_authorised") != "NO":
    failures.append("Combined summary does not explicitly block broad biological interpretation.")

allowed_rows = [
    {
        "statement_id": "ALLOWED_TECH_001",
        "statement_class": "technical_qc",
        "allowed_statement": "Strict genotype-only Phase C analyses generated QC-controlled allele-frequency and group-frequency tables from validated strict genotype matrices.",
        "required_caveat": "This is a technical preparation/QC statement only.",
        "source_checkpoint": "strict_group_frequency_qc",
    },
    {
        "statement_id": "ALLOWED_TECH_002",
        "statement_class": "candidate_prioritisation",
        "allowed_statement": "The candidate-prioritisation workflow identified 500 applied candidates and 173 follow-up candidates under predeclared technical rules.",
        "required_caveat": "Prioritisation is a follow-up triage framework, not biological validation.",
        "source_checkpoint": "candidate_prioritisation_applied",
    },
    {
        "statement_id": "ALLOWED_TECH_003",
        "statement_class": "annotation_supported_candidates",
        "allowed_statement": "Four annotation-supported candidates were reviewed in first-pass evidence review across TRAF3, BIRC3 and TNFAIP3.",
        "required_caveat": "Candidate-specific evidence grades and boundaries must be preserved.",
        "source_checkpoint": "combined_annotation_supported_evidence_summary",
    },
    {
        "statement_id": "ALLOWED_EVID_001",
        "statement_class": "TRAF3_evidence",
        "allowed_statement": "TRAF3 chr14:102875712:T:C maps to rs1131877/p.Met129Thr with exact-variant database/protein-coordinate support.",
        "required_caveat": "No direct functional validation or evolutionary inference is established.",
        "source_checkpoint": "TRAF3_evidence_review",
    },
    {
        "statement_id": "ALLOWED_EVID_002",
        "statement_class": "BIRC3_evidence",
        "allowed_statement": "BIRC3 chr11:102324736:G:T has local VEP identity support as rs924735278 with W/L missense context.",
        "required_caveat": "Exact external variant evidence and direct functional validation were not established.",
        "source_checkpoint": "BIRC3_evidence_review",
    },
    {
        "statement_id": "ALLOWED_EVID_003",
        "statement_class": "TNFAIP3_evidence",
        "allowed_statement": "TNFAIP3 chr6:137874871:A:G and chr6:137875820:A:C map to rs376205580/T108A and rs141807543/I207L, respectively, with paired T108A;I207L haplotype-level functional context.",
        "required_caveat": "Do not treat haplotype-level evidence as independent single-variant functional validation.",
        "source_checkpoint": "TNFAIP3_paired_evidence_review",
    },
]

prohibited_rows = [
    {
        "claim_id": "PROHIBIT_001",
        "prohibited_claim": "The Phase C candidates are causal variants.",
        "reason": "Causality was not tested or validated.",
        "replacement_language": "technical candidates or evidence-review candidates",
    },
    {
        "claim_id": "PROHIBIT_002",
        "prohibited_claim": "Frequency differences show selection or adaptation.",
        "reason": "Formal population-genetic tests were not performed in this checkpoint.",
        "replacement_language": "descriptive frequency variation under strict QC",
    },
    {
        "claim_id": "PROHIBIT_003",
        "prohibited_claim": "Candidates are population-specific or ancestry-determining.",
        "reason": "The workflow does not authorise ancestry-effect claims.",
        "replacement_language": "group-frequency patterns requiring formal follow-up",
    },
    {
        "claim_id": "PROHIBIT_004",
        "prohibited_claim": "TNFAIP3 T108A and I207L are independently validated as single-variant functional effects by this project.",
        "reason": "Functional evidence is paired/haplotype-level from prior literature, not independent project validation.",
        "replacement_language": "paired T108A;I207L haplotype-level functional context",
    },
    {
        "claim_id": "PROHIBIT_005",
        "prohibited_claim": "BIRC3 rs924735278 is externally validated or functionally validated.",
        "reason": "Only local VEP identity/consequence context and gene-level context were established.",
        "replacement_language": "local VEP identity support with indirect variant context",
    },
]

candidate_language_rows = []
for row in combined_rows:
    candidate_language_rows.append(
        {
            "allele_key": row["allele_key"],
            "gene": row["gene"],
            "mapped_identity": row["mapped_identity"],
            "overall_evidence_grade": row["overall_evidence_grade"],
            "allowed_language": row["allowed_language"],
            "prohibited_language": row["prohibited_language"],
            "interpretation_boundary": row["interpretation_boundary"],
            "manuscript_readiness_label": "bounded_results_ready",
        }
    )

manuscript_gate_rows = [
    {
        "section": "Results",
        "status": "READY_WITH_BOUNDARIES",
        "permitted_content": "Technical QC outputs, candidate prioritisation summaries and bounded annotation-supported evidence summary.",
        "blocked_content": "Selection claims, adaptive interpretation, causal claims, ancestry effects, independent functional validation claims.",
    },
    {
        "section": "Discussion",
        "status": "READY_WITH_BOUNDARIES",
        "permitted_content": "Interpretive framing as hypothesis-generating and evidence-tiered.",
        "blocked_content": "Strong mechanistic conclusions or population-specific biological claims.",
    },
    {
        "section": "Limitations",
        "status": "REQUIRED",
        "permitted_content": "Explicit limitations on genotype-only data, descriptive frequencies, evidence grading and lack of independent functional validation.",
        "blocked_content": "Omission of claim boundaries.",
    },
    {
        "section": "Figures_Tables",
        "status": "READY_FOR_DESIGN",
        "permitted_content": "QC flow, prioritisation flow, evidence-grade table and bounded candidate summary.",
        "blocked_content": "Visual ranking implying biological importance or population superiority.",
    },
]

allowed_path = OUTDIR / "phaseC_allowed_manuscript_statements_v1.tsv"
prohibited_path = OUTDIR / "phaseC_prohibited_manuscript_claims_v1.tsv"
candidate_path = OUTDIR / "phaseC_candidate_manuscript_language_boundaries_v1.tsv"
gate_path = OUTDIR / "phaseC_manuscript_readiness_gate_v1.tsv"

write_rows(
    allowed_path,
    ["statement_id", "statement_class", "allowed_statement", "required_caveat", "source_checkpoint"],
    allowed_rows,
)
write_rows(
    prohibited_path,
    ["claim_id", "prohibited_claim", "reason", "replacement_language"],
    prohibited_rows,
)
write_rows(
    candidate_path,
    [
        "allele_key",
        "gene",
        "mapped_identity",
        "overall_evidence_grade",
        "allowed_language",
        "prohibited_language",
        "interpretation_boundary",
        "manuscript_readiness_label",
    ],
    candidate_language_rows,
)
write_rows(
    gate_path,
    ["section", "status", "permitted_content", "blocked_content"],
    manuscript_gate_rows,
)

summary_rows = [
    {"metric": "allowed_statement_rows", "value": len(allowed_rows)},
    {"metric": "prohibited_claim_rows", "value": len(prohibited_rows)},
    {"metric": "candidate_language_boundary_rows", "value": len(candidate_language_rows)},
    {"metric": "manuscript_gate_rows", "value": len(manuscript_gate_rows)},
    {"metric": "results_ready_with_boundaries", "value": "YES"},
    {"metric": "discussion_ready_with_boundaries", "value": "YES"},
    {"metric": "limitations_required", "value": "YES"},
    {"metric": "selection_or_population_claim_authorised", "value": "NO"},
    {"metric": "broad_biological_interpretation_authorised", "value": "NO"},
]
summary_path = OUTDIR / "phaseC_manuscript_readiness_gate_summary_v1.tsv"
write_rows(summary_path, ["metric", "value"], summary_rows)

payload = {
    "status": "PASS" if not failures else "FAIL",
    "operation": "Prepare Phase C interpretation-boundary and manuscript-readiness gate.",
    "allowed_statement_rows": len(allowed_rows),
    "prohibited_claim_rows": len(prohibited_rows),
    "candidate_language_boundary_rows": len(candidate_language_rows),
    "manuscript_gate_rows": len(manuscript_gate_rows),
    "results_ready_with_boundaries": "YES",
    "discussion_ready_with_boundaries": "YES",
    "limitations_required": "YES",
    "selection_or_population_claim_authorised": "NO",
    "broad_biological_interpretation_authorised": "NO",
    "outputs": {
        "allowed_manuscript_statements": str(allowed_path),
        "prohibited_manuscript_claims": str(prohibited_path),
        "candidate_language_boundaries": str(candidate_path),
        "manuscript_readiness_gate": str(gate_path),
        "summary": str(summary_path),
    },
    "failures": failures,
    "boundary": (
        "This checkpoint prepares controlled manuscript-readiness boundaries only. It does not write Results prose, "
        "perform new evidence review, validate function, infer selection, make population-differentiation claims, "
        "assign ancestry effects or authorise broad biological interpretation."
    ),
}

OUTJSON.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
print(json.dumps(payload, indent=2))

if payload["status"] != "PASS":
    raise SystemExit(1)
