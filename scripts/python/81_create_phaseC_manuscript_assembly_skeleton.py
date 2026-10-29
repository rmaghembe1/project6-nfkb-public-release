#!/usr/bin/env python3
"""Create Phase C manuscript assembly skeleton.

This checkpoint creates modular manuscript assembly files from the completed
Phase C manuscript-output package. It does not perform new analysis, validate
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
OUTJSON = Path("results/qc/phaseC_manuscript_assembly_skeleton_v1.json")

for path in [ASSEMBLY, SECTIONS, AUDIT, OUTJSON.parent]:
    path.mkdir(parents=True, exist_ok=True)

required_inputs = [
    BASE / "results_text" / "phaseC_bounded_results_text_draft_v1.md",
    BASE / "discussion_text" / "phaseC_bounded_discussion_limitations_text_draft_v1.md",
    BASE / "captions" / "figure1_caption_v1.md",
    BASE / "captions" / "figure2_caption_v1.md",
    BASE / "captions" / "figure3_caption_v1.md",
    BASE / "captions" / "main_table1_caption_v1.md",
    BASE / "captions" / "main_table2_caption_v1.md",
    BASE / "captions" / "supp_table1_caption_v1.md",
    BASE / "captions" / "supp_table2_caption_v1.md",
    BASE / "audit" / "phaseC_manuscript_output_readiness_summary_v1.tsv",
]

failures = []
for path in required_inputs:
    if not path.exists():
        failures.append(f"Missing required input: {path}")

title_path = SECTIONS / "00_title_and_scope_placeholder_v1.md"
title_path.write_text(
    "# Manuscript assembly placeholder: Phase C genotype-only NF-κB regulatory-circuitry analysis\n\n"
    "This file is a placeholder for manuscript title, author information and scope framing. "
    "It does not add new scientific interpretation.\n",
    encoding="utf-8",
)

intro_path = SECTIONS / "01_introduction_scope_placeholder_v1.md"
intro_path.write_text(
    "# Introduction scope placeholder\n\n"
    "This section remains a placeholder. The current assembly checkpoint does not draft a full introduction. "
    "Future introduction text should frame the project as a population-resolved, genotype-only and evidence-gated analysis of NF-κB regulatory circuitry without implying causal, adaptive, ancestry-effect or population-specific claims from Phase C outputs.\n",
    encoding="utf-8",
)

methods_path = SECTIONS / "02_methods_phaseC_assembly_placeholder_v1.md"
methods_path.write_text(
    "# Methods assembly placeholder\n\n"
    "This section will later assemble the Phase C genotype-only QC, candidate prioritisation, VEP annotation import, evidence-review gate, manuscript-readiness gate and display-package generation methods. "
    "No new methods prose is finalised in this checkpoint.\n",
    encoding="utf-8",
)

results_source = BASE / "results_text" / "phaseC_bounded_results_text_draft_v1.md"
results_path = SECTIONS / "03_results_phaseC_bounded_v1.md"
if results_source.exists():
    results_path.write_text(results_source.read_text(encoding="utf-8"), encoding="utf-8")

discussion_source = BASE / "discussion_text" / "phaseC_bounded_discussion_limitations_text_draft_v1.md"
discussion_path = SECTIONS / "04_discussion_limitations_phaseC_bounded_v1.md"
if discussion_source.exists():
    discussion_path.write_text(discussion_source.read_text(encoding="utf-8"), encoding="utf-8")

display_path = SECTIONS / "05_figures_tables_callouts_v1.md"
display_path.write_text(
    "# Figures and tables\n\n"
    "## Main figures\n\n"
    "- Figure 1. Phase C genotype-only QC and manuscript-readiness workflow.\n"
    "- Figure 2. Phase C candidate prioritisation and evidence-review flow.\n"
    "- Figure 3. Phase C annotation-supported evidence landscape.\n\n"
    "## Main tables\n\n"
    "- Table 1. Annotation-supported Phase C candidates and evidence-review status.\n"
    "- Table 2. Candidate-specific manuscript language boundaries for Phase C.\n\n"
    "## Supplementary tables\n\n"
    "- Supplementary Table 1. Phase C results package index.\n"
    "- Supplementary Table 2. Allowed and prohibited Phase C manuscript statements.\n",
    encoding="utf-8",
)

boundary_path = SECTIONS / "06_manuscript_boundaries_v1.md"
boundary_path.write_text(
    "# Manuscript interpretation boundaries\n\n"
    "All manuscript sections assembled from Phase C outputs must preserve the manuscript-readiness gate. "
    "Phase C outputs may be reported as genotype-only QC, prioritisation, annotation-supported evidence-review and manuscript-display outputs. "
    "They must not be converted into causal, adaptive, ancestry-effect, population-specific, independent functional-validation or broad biological claims unless a future checkpoint explicitly supports such claims.\n",
    encoding="utf-8",
)

master_path = ASSEMBLY / "phaseC_manuscript_assembly_master_v1.md"
section_files = [
    title_path,
    intro_path,
    methods_path,
    results_path,
    discussion_path,
    display_path,
    boundary_path,
]
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

manifest_path = AUDIT / "phaseC_manuscript_assembly_skeleton_manifest_v1.tsv"
manifest_rows = []
for path in section_files + [master_path]:
    manifest_rows.append(
        {
            "file": str(path),
            "exists": "YES" if path.exists() else "NO",
            "size_bytes": path.stat().st_size if path.exists() else 0,
            "role": "master" if path == master_path else "section",
        }
    )
with manifest_path.open("w", newline="", encoding="utf-8") as handle:
    writer = csv.DictWriter(handle, fieldnames=["file", "exists", "size_bytes", "role"], delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(manifest_rows)

summary_path = AUDIT / "phaseC_manuscript_assembly_skeleton_summary_v1.tsv"
summary_rows = [
    {"metric": "section_files_declared", "value": len(section_files)},
    {"metric": "master_file_written", "value": "YES" if master_path.exists() else "NO"},
    {"metric": "required_inputs_checked", "value": len(required_inputs)},
    {"metric": "missing_inputs", "value": len([f for f in failures if f.startswith("Missing required input")])},
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
    "operation": "Create Phase C manuscript assembly skeleton.",
    "section_files_declared": len(section_files),
    "master_file_written": "YES" if master_path.exists() else "NO",
    "required_inputs_checked": len(required_inputs),
    "missing_inputs": len([f for f in failures if f.startswith("Missing required input")]),
    "new_analysis_performed": "NO",
    "new_biological_claims_added": "NO",
    "selection_or_population_claim_authorised": "NO",
    "broad_biological_interpretation_authorised": "NO",
    "outputs": {
        "master": str(master_path),
        "manifest": str(manifest_path),
        "summary": str(summary_path),
    },
    "failures": failures,
    "boundary": (
        "This checkpoint creates a manuscript assembly skeleton only. It does not perform new analysis, validate function, "
        "infer selection, make population-differentiation claims, assign ancestry effects or authorise broad biological interpretation."
    ),
}

OUTJSON.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
print(json.dumps(payload, indent=2))

if payload["status"] != "PASS":
    raise SystemExit(1)
