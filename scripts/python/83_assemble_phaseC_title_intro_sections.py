#!/usr/bin/env python3
"""Assemble bounded Phase C title/scope and introduction sections.

This checkpoint replaces title/scope and introduction placeholders with
controlled manuscript-ready text. It does not perform new analysis, validate
function, infer selection, assign ancestry effects or authorise broad biological
interpretation.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path

BASE = Path("manuscript/phaseC")
ASSEMBLY = BASE / "assembly"
SECTIONS = ASSEMBLY / "sections"
AUDIT = ASSEMBLY / "audit"
OUTJSON = Path("results/qc/phaseC_title_intro_assembly_v1.json")

for path in [SECTIONS, AUDIT, OUTJSON.parent]:
    path.mkdir(parents=True, exist_ok=True)

required_inputs = [
    SECTIONS / "02_methods_phaseC_assembly_v1.md",
    SECTIONS / "03_results_phaseC_bounded_v1.md",
    SECTIONS / "04_discussion_limitations_phaseC_bounded_v1.md",
    SECTIONS / "06_manuscript_boundaries_v1.md",
    Path("docs/manuscript/phaseC/PhaseC_assembly_branch_activation_note.md"),
    Path("docs/manuscript/phaseC/PhaseC_manuscript_output_branch_closure_note.md"),
]

failures = []
for path in required_inputs:
    if not path.exists():
        failures.append(f"Missing required input: {path}")

title_text = """# Population-resolved NF-κB regulatory-circuitry variation: a genotype-only, evidence-gated Phase C manuscript assembly

## Scope statement

This manuscript assembly focuses on the Phase C genotype-only QC, candidate-prioritisation, annotation-supported evidence-review and manuscript-readiness outputs generated for the population-resolved NF-κB regulatory-circuitry project. The assembly reports audit-controlled workflow outputs, candidate triage, evidence-review status, figures, tables and interpretation boundaries. It does not present experimental validation, disease association testing, selection inference, ancestry-effect inference or causal functional claims.
"""

intro_text = """# Introduction

NF-κB signalling is a central regulatory axis in immune activation, inflammatory restraint, stress responses and host defence. The pathway integrates canonical and noncanonical signalling components, ubiquitin-editing machinery, inhibitory regulators and transcriptional outputs. Because NF-κB biology is embedded in complex regulatory networks, candidate genetic variation in this circuitry requires careful interpretation: a variant may be computationally prioritised or annotation-supported without being experimentally validated, causal, adaptive or population-specific.

The broader project was designed to examine population-resolved variation in human NF-κB regulatory circuitry using a structured genomic and annotation framework. Phase C represents a controlled manuscript-output stage of that project. Its role is to organise genotype-only outputs, apply candidate-prioritisation and evidence-review gates, prepare manuscript-facing figures and tables, and define what language can and cannot be used in a manuscript. This evidence-gated approach is intended to reduce overinterpretation of descriptive genotype-derived findings.

The present assembly therefore reports Phase C as a bounded analysis and reporting layer. It summarises the genotype-only QC and manuscript-readiness workflow, candidate-prioritisation flow, annotation-supported candidate evidence landscape, formatted tables and interpretation boundaries. The emphasis is on transparency, reproducibility and claim control. The analysis does not test functional consequences directly, does not establish selection or adaptation, does not assign ancestry effects and does not infer population-specific biological differences.

This framing is important because NF-κB genes are biologically consequential and easily overinterpreted. Genes such as TRAF3, BIRC3 and TNFAIP3 have known relevance to immune signalling at the gene or pathway level, but gene-level biological importance is not equivalent to direct evidence for a specific candidate variant. Phase C therefore separates candidate identity, annotation support, evidence-review status and manuscript language boundaries. This separation allows candidates to be reported as prioritised or evidence-reviewed while preserving uncertainty about function, causality and population interpretation.

The objective of this manuscript assembly is to present the Phase C workflow and outputs in a format suitable for downstream manuscript development. The assembly integrates bounded Results and Discussion/limitations text with Methods, figure/table callouts and interpretation boundaries. All sections remain under the Phase C manuscript-readiness gate.
"""

title_path = SECTIONS / "00_title_and_scope_v1.md"
intro_path = SECTIONS / "01_introduction_phaseC_scope_v1.md"
title_path.write_text(title_text, encoding="utf-8")
intro_path.write_text(intro_text, encoding="utf-8")

# Archive placeholders if present.
placeholder_pairs = [
    (SECTIONS / "00_title_and_scope_placeholder_v1.md", SECTIONS / "00_title_and_scope_placeholder_archived_v1.md"),
    (SECTIONS / "01_introduction_scope_placeholder_v1.md", SECTIONS / "01_introduction_scope_placeholder_archived_v1.md"),
]
for src, dst in placeholder_pairs:
    if src.exists():
        dst.write_text(src.read_text(encoding="utf-8") + f"\n\nArchived after replacement by {dst.name.replace('_placeholder_archived', '')}.\n", encoding="utf-8")

section_files = [
    title_path,
    intro_path,
    SECTIONS / "02_methods_phaseC_assembly_v1.md",
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
        out.write(path.read_text(encoding="utf-8") if path.exists() else f"[Missing section file: {path}]\n")
        out.write(f"\n\n<!-- END {path.name} -->\n")

combined = (title_text + "\n" + intro_text).lower()
for phrase in [
    "demonstrates selection",
    "proves adaptation",
    "establishes ancestry effects",
    "causal variant",
    "population-specific variant",
    "functionally validated by phase c",
]:
    if phrase in combined:
        failures.append(f"Unsupported title/introduction phrase detected: {phrase}")

summary_path = AUDIT / "phaseC_title_intro_assembly_summary_v1.tsv"
summary_rows = [
    {"metric": "title_scope_written", "value": "YES"},
    {"metric": "introduction_written", "value": "YES"},
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
    "operation": "Assemble bounded Phase C title/scope and introduction sections.",
    "title_scope_written": "YES",
    "introduction_written": "YES",
    "required_inputs_checked": len(required_inputs),
    "missing_inputs": len([f for f in failures if f.startswith("Missing required input")]),
    "master_assembly_regenerated": "YES",
    "new_analysis_performed": "NO",
    "new_biological_claims_added": "NO",
    "selection_or_population_claim_authorised": "NO",
    "broad_biological_interpretation_authorised": "NO",
    "outputs": {
        "title_scope": str(title_path),
        "introduction": str(intro_path),
        "master": str(master_path),
        "summary": str(summary_path),
    },
    "failures": failures,
    "boundary": (
        "This checkpoint assembles bounded title/scope and introduction sections only. It does not perform new analysis, "
        "validate function, infer selection, assign ancestry effects or authorise broad biological interpretation."
    ),
}

OUTJSON.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
print(json.dumps(payload, indent=2))

if payload["status"] != "PASS":
    raise SystemExit(1)
