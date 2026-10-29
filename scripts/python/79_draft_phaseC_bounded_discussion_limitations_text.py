#!/usr/bin/env python3
"""Draft bounded Phase C Discussion and limitations text.

This checkpoint drafts manuscript-facing Discussion and limitations text from
the completed Phase C display package and bounded Results text. It does not
perform new analysis, validate function, infer selection, assign ancestry effects
or authorise broad biological interpretation.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path

OUTDIR = Path("manuscript/phaseC/discussion_text")
AUDITDIR = Path("manuscript/phaseC/audit")
OUTJSON = Path("results/qc/phaseC_bounded_discussion_limitations_text_draft_v1.json")

for path in [OUTDIR, AUDITDIR, OUTJSON.parent]:
    path.mkdir(parents=True, exist_ok=True)

required_inputs = [
    Path("manuscript/phaseC/results_text/phaseC_bounded_results_text_draft_v1.md"),
    Path("manuscript/phaseC/audit/phaseC_display_package_audit_summary_v1.tsv"),
    Path("manuscript/phaseC/tables/main/formatted/main_table2_phaseC_candidate_language_boundaries_formatted_v1.tsv"),
    Path("manuscript/phaseC/tables/supplementary/formatted/supp_table2A_allowed_phaseC_manuscript_statements_formatted_v1.tsv"),
    Path("manuscript/phaseC/tables/supplementary/formatted/supp_table2B_prohibited_phaseC_manuscript_claims_formatted_v1.tsv"),
]

failures = []
for path in required_inputs:
    if not path.exists():
        failures.append(f"Missing required input: {path}")

discussion_text = """# Phase C Discussion and limitations text draft v1

## Summary of Phase C contribution

Phase C provides a controlled genotype-only analysis and manuscript-display layer for the NF-κB regulatory-circuitry project. Its main contribution is not a claim of functional validation or adaptive differentiation, but a reproducible framework for moving from strict genotype-only QC to candidate prioritisation, annotation-supported evidence review, manuscript-readiness gating and auditable figure/table generation. The display package and bounded Results text provide a traceable basis for reporting candidate status while preserving explicit interpretation limits.

## Interpretation of candidate prioritisation

The Phase C prioritisation framework should be interpreted as a technical triage system. It helped organise candidates for annotation review and evidence checking, but the resulting scores and follow-up labels do not measure biological importance, causal effect size, adaptive value or population specificity. The prioritisation outputs are therefore most useful for transparent manuscript organisation and for identifying which candidates require future annotation, functional or population-genetic follow-up.

## Annotation-supported evidence review

The first annotation-supported review pass identified four candidates across TRAF3, BIRC3 and paired TNFAIP3 for bounded reporting. These candidates can be described as annotation-supported and evidence-reviewed within the limits recorded by the manuscript-readiness gate. However, Phase C does not demonstrate that these variants alter NF-κB signalling, inflammatory-response phenotypes, disease risk or treatment response. Where external evidence exists for gene-level or broader pathway relevance, it should be separated from direct evidence for the specific candidate variant or paired candidate state.

## Manuscript-readiness boundaries

A central strength of Phase C is that it explicitly records allowed and prohibited manuscript language. This reduces the risk of overinterpreting genotype-only descriptive outputs as causal, adaptive, ancestry-determining or functionally validated findings. The manuscript should preserve these boundaries by reporting candidate prioritisation, annotation identity, evidence grade and review status as audit-controlled outputs rather than as definitive biological conclusions.

## Limitations

Phase C is limited by its genotype-only design. It does not include direct experimental perturbation, expression-QTL testing, chromatin-response validation, protein-function assays, haplotype-resolved molecular phenotyping or clinical association testing. The current analysis also does not establish selection, adaptation, ancestry effects, population specificity or independent single-variant functional validation. For paired TNFAIP3 candidates, wording should remain paired or haplotype-level unless a future checkpoint explicitly resolves independent single-variant effects.

## Future work

Future work should use the Phase C outputs as a prioritisation and reporting framework rather than as final biological proof. Useful next steps include independent annotation verification, haplotype-aware analysis, stimulus-specific regulatory annotation, population-genetic neutrality or selection testing where appropriate, and functional assays for candidates with sufficient evidence support. Any future claim about functional effect, population differentiation, adaptation or disease relevance should be introduced only after a dedicated checkpoint supports that claim.

## Discussion-text boundary

The present Discussion draft is intentionally conservative. It explains what Phase C contributes, why the manuscript-readiness gate matters and what limitations must remain visible. It does not convert candidate prioritisation or annotation support into causal, adaptive, population-specific, ancestry-effect or independently functionally validated claims.
"""

discussion_path = OUTDIR / "phaseC_bounded_discussion_limitations_text_draft_v1.md"
discussion_path.write_text(discussion_text, encoding="utf-8")

lower_text = discussion_text.lower()

unsupported_affirmative_phrases = [
    "demonstrates selection",
    "proves adaptation",
    "establishes adaptation",
    "confirms ancestry effects",
    "establishes ancestry effects",
    "is a causal variant",
    "are causal variants",
    "functionally validated by phase c",
    "phase c validates function",
    "population-specific variant",
    "ancestry-determining variant",
]

for phrase in unsupported_affirmative_phrases:
    if phrase in lower_text:
        failures.append(f"Unsupported affirmative phrase detected: {phrase}")

summary_path = AUDITDIR / "phaseC_bounded_discussion_limitations_text_draft_summary_v1.tsv"
summary_rows = [
    {"metric": "discussion_text_written", "value": "YES"},
    {"metric": "discussion_sections", "value": 7},
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
    "operation": "Draft bounded Phase C Discussion and limitations text.",
    "discussion_text_written": "YES",
    "discussion_sections": 7,
    "required_inputs_checked": len(required_inputs),
    "missing_inputs": len([f for f in failures if f.startswith("Missing required input")]),
    "new_analysis_performed": "NO",
    "figures_generated": "NO",
    "tables_generated": "NO",
    "selection_or_population_claim_authorised": "NO",
    "broad_biological_interpretation_authorised": "NO",
    "outputs": {
        "discussion_text": str(discussion_path),
        "summary": str(summary_path),
    },
    "failures": failures,
    "boundary": (
        "This checkpoint drafts bounded Discussion and limitations text only. It does not perform new analysis, "
        "generate figures, generate tables, validate function, infer selection, make population-differentiation claims, "
        "assign ancestry effects or authorise broad biological interpretation."
    ),
}

OUTJSON.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
print(json.dumps(payload, indent=2))

if payload["status"] != "PASS":
    raise SystemExit(1)
