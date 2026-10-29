#!/usr/bin/env python3
"""Assemble bounded Phase C Methods section.

This checkpoint drafts a Methods section from existing Phase C audit records and
manuscript-output files. It does not perform new analysis, validate function,
infer selection, assign ancestry effects or authorise broad biological interpretation.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path

BASE = Path("manuscript/phaseC")
ASSEMBLY = BASE / "assembly"
SECTIONS = ASSEMBLY / "sections"
AUDIT = ASSEMBLY / "audit"
OUTJSON = Path("results/qc/phaseC_methods_assembly_v1.json")

for path in [SECTIONS, AUDIT, OUTJSON.parent]:
    path.mkdir(parents=True, exist_ok=True)

required_inputs = [
    Path("docs/phaseC/PhaseC_branch_closure_and_next_step_note.md"),
    Path("docs/manuscript/phaseC/PhaseC_manuscript_output_branch_closure_note.md"),
    Path("manuscript/phaseC/audit/phaseC_manuscript_output_readiness_summary_v1.tsv"),
    Path("manuscript/phaseC/audit/phaseC_display_package_audit_summary_v1.tsv"),
    Path("results/qc/phaseC_manuscript_source_data_export_v1.json"),
    Path("results/qc/phaseC_main_tables_formatting_v1.json"),
    Path("results/qc/phaseC_supplementary_tables_formatting_v1.json"),
    Path("results/qc/phaseC_figure1_generation_v1.json"),
    Path("results/qc/phaseC_figure2_generation_v1.json"),
    Path("results/qc/phaseC_figure3_generation_v1.json"),
    Path("results/qc/phaseC_display_package_audit_v1.json"),
]

failures = []
for path in required_inputs:
    if not path.exists():
        failures.append(f"Missing required input: {path}")

methods_text = """# Methods

## Phase C analysis scope

Phase C was designed as a controlled genotype-only QC, candidate-prioritisation and manuscript-output activation stage for the NF-κB regulatory-circuitry project. The phase was separated from downstream biological interpretation by explicit audit gates. Its purpose was to organise genotype-derived candidate information, apply predefined prioritisation and evidence-review boundaries, prepare manuscript-facing displays and generate bounded manuscript text. Phase C did not include new wet-lab validation, expression analysis, chromatin-response modelling, association testing, selection testing or ancestry-effect inference.

## Genotype-only QC and candidate-prioritisation workflow

The Phase C technical workflow began from committed genotype-only outputs and proceeded through strict allele-count, frequency and group-frequency QC gates. Candidate prioritisation was then applied as a technical triage framework. Prioritisation outputs were used to define applied candidates and follow-up candidates for annotation and evidence review. These outputs were treated as audit-controlled organisational categories rather than measures of biological importance, causality, adaptive value or population specificity.

## Annotation import and evidence-review gating

Returned VEP annotations were imported into the Phase C workflow and profiled to support consequence-context review. Annotation-supported candidates were then routed through evidence-review gates. Evidence review was performed for TRAF3, BIRC3 and paired TNFAIP3 candidates, followed by a combined annotation-supported evidence summary. The evidence-review process recorded candidate identity, evidence grade, functional-validation status, selection/population-claim status and interpretation boundaries. Candidate language boundaries were explicitly separated into allowed and prohibited manuscript statements to prevent unsupported causal, adaptive, ancestry-effect, population-specific or independent functional-validation claims.

## Manuscript-readiness gate

A manuscript-readiness gate was prepared before generating manuscript-facing outputs. This gate defined allowed manuscript statements, prohibited claims and candidate-specific wording boundaries. It restricted Phase C reporting to genotype-only QC, candidate-prioritisation, annotation-supported evidence-review status and display-package completeness. The gate did not authorise biological claims beyond the committed evidence level.

## Manuscript source-data export and display generation

Committed Phase C source tables were exported into the manuscript-output display package. Main Table 1 summarised annotation-supported candidate evidence-review status, and Main Table 2 summarised candidate-specific manuscript language boundaries. Supplementary Table 1 indexed the Phase C results package, while Supplementary Table 2A/2B documented allowed and prohibited manuscript statements. Figures 1–3 were generated as manuscript-facing workflow, prioritisation/evidence-flow and evidence-landscape displays. Each figure was exported as an editable SVG and high-resolution PNG, and each figure and table received a bounded caption.

## Display-package and manuscript-output audit

The display package was audited for completeness, file presence, non-empty outputs and SVG whitespace status. A manuscript-output readiness audit then checked the full package, including figures, tables, captions, bounded Results text, bounded Discussion/limitations text, QC JSONs and branch status. These audits verified readiness for manuscript assembly while preserving the Phase C manuscript-readiness boundaries.

## Manuscript assembly

The manuscript assembly branch was created from the closed manuscript-output branch. Modular Markdown section files were generated for title/scope, introduction placeholder, Methods, Results, Discussion/limitations, figure/table callouts and manuscript interpretation boundaries. The assembly process reused bounded Results and Discussion/limitations drafts from the completed manuscript-output package. No new analysis or new biological claim was introduced during manuscript assembly.
"""

methods_path = SECTIONS / "02_methods_phaseC_assembly_v1.md"
methods_path.write_text(methods_text, encoding="utf-8")

# Preserve a deprecated placeholder copy if it exists.
placeholder = SECTIONS / "02_methods_phaseC_assembly_placeholder_v1.md"
if placeholder.exists():
    placeholder_archive = SECTIONS / "02_methods_phaseC_assembly_placeholder_archived_v1.md"
    placeholder_archive.write_text(
        placeholder.read_text(encoding="utf-8") +
        "\n\nArchived after replacement by 02_methods_phaseC_assembly_v1.md.\n",
        encoding="utf-8",
    )

section_files = [
    SECTIONS / "00_title_and_scope_placeholder_v1.md",
    SECTIONS / "01_introduction_scope_placeholder_v1.md",
    methods_path,
    SECTIONS / "03_results_phaseC_bounded_v1.md",
    SECTIONS / "04_discussion_limitations_phaseC_bounded_v1.md",
    SECTIONS / "05_figures_tables_callouts_v1.md",
    SECTIONS / "06_manuscript_boundaries_v1.md",
]

for path in section_files:
    if not path.exists():
        failures.append(f"Missing assembly section: {path}")

master_path = ASSEMBLY / "phaseC_manuscript_assembly_master_v1.md"
with master_path.open("w", encoding="utf-8") as out:
    out.write("# Phase C manuscript assembly master v1\n\n")
    out.write("This master file assembles controlled Phase C manuscript sections under the manuscript-readiness gate.\n\n")
    for path in section_files:
        out.write(f"\n\n<!-- BEGIN {path.name} -->\n\n")
        if path.exists():
            out.write(path.read_text(encoding="utf-8"))
        else:
            out.write(f"[Missing section file: {path}]\n")
        out.write(f"\n\n<!-- END {path.name} -->\n")

lower_text = methods_text.lower()
for phrase in [
    "demonstrates selection",
    "proves adaptation",
    "causal variant",
    "ancestry-determining variant",
    "population-specific variant",
    "functionally validated by phase c",
]:
    if phrase in lower_text:
        failures.append(f"Unsupported Methods phrase detected: {phrase}")

summary_path = AUDIT / "phaseC_methods_assembly_summary_v1.tsv"
summary_rows = [
    {"metric": "methods_section_written", "value": "YES"},
    {"metric": "required_inputs_checked", "value": len(required_inputs)},
    {"metric": "missing_inputs", "value": len([f for f in failures if f.startswith("Missing required input")])},
    {"metric": "master_assembly_regenerated", "value": "YES"},
    {"metric": "new_analysis_performed", "value": "NO"},
    {"metric": "new_biological_claims_added", "value": "NO"},
    {"metric": "selection_or_population_claim_authorised", "value": "NO"},
    {"metric": "broad_biological_interpretation_authorised", "value": "NO"},
]
with summary_path.open("w", newline="", encoding="utf-8") as handle:
    writer = csv.DictWriter(handle, fieldnames=["metric", "value"], delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(summary_rows)

payload = {
    "status": "PASS" if not failures else "FAIL",
    "operation": "Assemble bounded Phase C Methods section.",
    "methods_section_written": "YES",
    "required_inputs_checked": len(required_inputs),
    "missing_inputs": len([f for f in failures if f.startswith("Missing required input")]),
    "master_assembly_regenerated": "YES",
    "new_analysis_performed": "NO",
    "new_biological_claims_added": "NO",
    "selection_or_population_claim_authorised": "NO",
    "broad_biological_interpretation_authorised": "NO",
    "outputs": {
        "methods_section": str(methods_path),
        "master": str(master_path),
        "summary": str(summary_path),
    },
    "failures": failures,
    "boundary": (
        "This checkpoint assembles a bounded Methods section only. It does not perform new analysis, validate function, "
        "infer selection, make population-differentiation claims, assign ancestry effects or authorise broad biological interpretation."
    ),
}

OUTJSON.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
print(json.dumps(payload, indent=2))

if payload["status"] != "PASS":
    raise SystemExit(1)
