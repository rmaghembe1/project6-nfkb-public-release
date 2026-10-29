#!/usr/bin/env python3
"""Draft bounded Phase C Results text.

This checkpoint drafts manuscript-facing Results text from the completed Phase C
display package and manuscript-readiness boundaries. It does not perform new
analysis, validate function, infer selection, assign ancestry effects or
authorise broad biological interpretation.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path

OUTDIR = Path("manuscript/phaseC/results_text")
AUDITDIR = Path("manuscript/phaseC/audit")
OUTJSON = Path("results/qc/phaseC_bounded_results_text_draft_v1.json")

for path in [OUTDIR, AUDITDIR, OUTJSON.parent]:
    path.mkdir(parents=True, exist_ok=True)

required_inputs = [
    Path("manuscript/phaseC/audit/phaseC_display_package_audit_summary_v1.tsv"),
    Path("manuscript/phaseC/tables/main/formatted/main_table1_annotation_supported_candidate_evidence_summary_formatted_v1.tsv"),
    Path("manuscript/phaseC/tables/main/formatted/main_table2_phaseC_candidate_language_boundaries_formatted_v1.tsv"),
    Path("manuscript/phaseC/captions/figure1_caption_v1.md"),
    Path("manuscript/phaseC/captions/figure2_caption_v1.md"),
    Path("manuscript/phaseC/captions/figure3_caption_v1.md"),
]

failures = []
for path in required_inputs:
    if not path.exists():
        failures.append(f"Missing required input: {path}")

results_text = """# Phase C Results text draft v1

## Phase C genotype-only QC and manuscript-readiness workflow

Phase C established a controlled genotype-only QC and manuscript-output workflow for the NF-κB regulatory-circuitry project. The workflow moved from strict genotype-only QC through candidate prioritisation, returned VEP annotation review, annotation-supported evidence review, manuscript-readiness gating, results-package indexing and display-package generation. Figure 1 summarises this audit-controlled workflow and its interpretation boundary. The workflow is descriptive and does not, by itself, support causal, adaptive, ancestry-effect, population-specific or functional-validation claims.

## Candidate prioritisation and evidence-review flow

The Phase C prioritisation framework was applied to 500 candidates and retained 173 follow-up candidates for downstream annotation and evidence triage. Four annotation-supported candidates were then included in the first evidence-review pass. Figure 2 summarises this transition from applied candidates to follow-up candidates and annotation-supported candidates. These counts are audit and prioritisation outputs only; they should not be interpreted as a ranking of biological importance, evidence of selection, population specificity, ancestry effects or functional causality.

## Annotation-supported candidate evidence landscape

The first annotation-supported evidence-review pass covered four candidates across TRAF3, BIRC3 and paired TNFAIP3. Figure 3 and Table 1 summarise candidate identity, evidence grade, functional-validation status, selection/population-claim status and interpretation boundaries. Table 2 provides candidate-specific language boundaries for manuscript reporting. Across these candidates, Phase C supports bounded reporting of annotation-supported and evidence-reviewed candidate status, but it does not establish causal effects, adaptive evolution, ancestry-determining effects, population-specific effects or independent single-variant functional validation.

## Display package completeness

The Phase C manuscript display package contains three figures, two main tables, supplementary tables, captions and audit files. The display-package audit confirmed that all 23 required display files were present, non-empty and free from SVG trailing-whitespace failures. This audit supports use of the display package for bounded manuscript drafting under the Phase C manuscript-readiness gate.

## Results-text boundary

All Phase C Results statements should remain tied to the committed audit outputs, figure/table package and manuscript-readiness gate. The current Results text is limited to genotype-only QC, candidate-prioritisation flow, annotation-supported evidence-review status and display-package completeness. It does not interpret the candidates as causal, adaptive, population-specific, ancestry-determining or independently functionally validated.
"""

results_path = OUTDIR / "phaseC_bounded_results_text_draft_v1.md"
results_path.write_text(results_text, encoding="utf-8")

forbidden_positive_phrases = [
    "causal variant",
    "selected variant",
    "ancestry-determining variant",
    "population-specific variant",
    "validated functional variant",
    "adaptive variant",
]

lower_text = results_text.lower()
for phrase in forbidden_positive_phrases:
    if phrase in lower_text:
        failures.append(f"Forbidden positive phrase detected in Results text: {phrase}")

# Allow "does not establish causal effects" etc., but block strong affirmative wording.
affirmative_blocks = [
    "establishes causal",
    "demonstrates selection",
    "proves adaptation",
    "confirms ancestry",
    "functionally validated by phase c",
]
for phrase in affirmative_blocks:
    if phrase in lower_text:
        failures.append(f"Unsupported affirmative claim detected: {phrase}")

summary_path = AUDITDIR / "phaseC_bounded_results_text_draft_summary_v1.tsv"
summary_rows = [
    {"metric": "results_text_written", "value": "YES"},
    {"metric": "results_sections", "value": 5},
    {"metric": "required_inputs_checked", "value": len(required_inputs)},
    {"metric": "missing_inputs", "value": len([f for f in failures if f.startswith("Missing required input")])},
    {"metric": "new_analysis_performed", "value": "NO"},
    {"metric": "figures_generated", "value": "NO"},
    {"metric": "tables_generated", "value": "NO"},
    {"metric": "selection_or_population_claim_authorised", "value": "NO"},
    {"metric": "broad_biological_interpretation_authorised", "value": "NO"},
]
with summary_path.open("w", newline="", encoding="utf-8") as handle:
    writer = csv.DictWriter(handle, fieldnames=["metric", "value"], delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(summary_rows)

payload = {
    "status": "PASS" if not failures else "FAIL",
    "operation": "Draft bounded Phase C Results text.",
    "results_text_written": "YES",
    "results_sections": 5,
    "required_inputs_checked": len(required_inputs),
    "missing_inputs": len([f for f in failures if f.startswith("Missing required input")]),
    "new_analysis_performed": "NO",
    "figures_generated": "NO",
    "tables_generated": "NO",
    "selection_or_population_claim_authorised": "NO",
    "broad_biological_interpretation_authorised": "NO",
    "outputs": {
        "results_text": str(results_path),
        "summary": str(summary_path),
    },
    "failures": failures,
    "boundary": (
        "This checkpoint drafts bounded Results text only. It does not perform new analysis, generate figures, "
        "generate tables, validate function, infer selection, make population-differentiation claims, assign ancestry "
        "effects or authorise broad biological interpretation."
    ),
}

OUTJSON.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
print(json.dumps(payload, indent=2))

if payload["status"] != "PASS":
    raise SystemExit(1)
