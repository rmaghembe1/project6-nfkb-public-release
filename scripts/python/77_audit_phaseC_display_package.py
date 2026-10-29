#!/usr/bin/env python3
"""Audit Phase C manuscript display package completeness.

This checkpoint audits the manuscript-facing Phase C figure/table package after
generation of planned figures and formatted tables. It does not perform new
analysis, write Results prose, validate function, infer selection, assign
ancestry effects or authorise broad biological interpretation.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path

OUTJSON = Path("results/qc/phaseC_display_package_audit_v1.json")
AUDIT_DIR = Path("manuscript/phaseC/audit")
OUTJSON.parent.mkdir(parents=True, exist_ok=True)
AUDIT_DIR.mkdir(parents=True, exist_ok=True)

required_files = [
    {
        "item": "Figure 1 SVG",
        "path": "manuscript/phaseC/figures/svg/figure1_phaseC_genotype_only_qc_workflow_v1.svg",
        "type": "figure_svg",
    },
    {
        "item": "Figure 1 PNG",
        "path": "manuscript/phaseC/figures/png/figure1_phaseC_genotype_only_qc_workflow_v1.png",
        "type": "figure_png",
    },
    {
        "item": "Figure 1 caption",
        "path": "manuscript/phaseC/captions/figure1_caption_v1.md",
        "type": "caption",
    },
    {
        "item": "Figure 2 SVG",
        "path": "manuscript/phaseC/figures/svg/figure2_phaseC_candidate_prioritisation_evidence_flow_v1.svg",
        "type": "figure_svg",
    },
    {
        "item": "Figure 2 PNG",
        "path": "manuscript/phaseC/figures/png/figure2_phaseC_candidate_prioritisation_evidence_flow_v1.png",
        "type": "figure_png",
    },
    {
        "item": "Figure 2 caption",
        "path": "manuscript/phaseC/captions/figure2_caption_v1.md",
        "type": "caption",
    },
    {
        "item": "Figure 3 SVG",
        "path": "manuscript/phaseC/figures/svg/figure3_phaseC_annotation_supported_evidence_landscape_v1.svg",
        "type": "figure_svg",
    },
    {
        "item": "Figure 3 PNG",
        "path": "manuscript/phaseC/figures/png/figure3_phaseC_annotation_supported_evidence_landscape_v1.png",
        "type": "figure_png",
    },
    {
        "item": "Figure 3 caption",
        "path": "manuscript/phaseC/captions/figure3_caption_v1.md",
        "type": "caption",
    },
    {
        "item": "Main Table 1 TSV",
        "path": "manuscript/phaseC/tables/main/formatted/main_table1_annotation_supported_candidate_evidence_summary_formatted_v1.tsv",
        "type": "main_table",
    },
    {
        "item": "Main Table 1 Markdown",
        "path": "manuscript/phaseC/tables/main/formatted/main_table1_annotation_supported_candidate_evidence_summary_formatted_v1.md",
        "type": "main_table",
    },
    {
        "item": "Main Table 1 caption",
        "path": "manuscript/phaseC/captions/main_table1_caption_v1.md",
        "type": "caption",
    },
    {
        "item": "Main Table 2 TSV",
        "path": "manuscript/phaseC/tables/main/formatted/main_table2_phaseC_candidate_language_boundaries_formatted_v1.tsv",
        "type": "main_table",
    },
    {
        "item": "Main Table 2 Markdown",
        "path": "manuscript/phaseC/tables/main/formatted/main_table2_phaseC_candidate_language_boundaries_formatted_v1.md",
        "type": "main_table",
    },
    {
        "item": "Main Table 2 caption",
        "path": "manuscript/phaseC/captions/main_table2_caption_v1.md",
        "type": "caption",
    },
    {
        "item": "Supplementary Table 1 TSV",
        "path": "manuscript/phaseC/tables/supplementary/formatted/supp_table1_phaseC_results_package_index_formatted_v1.tsv",
        "type": "supplementary_table",
    },
    {
        "item": "Supplementary Table 1 Markdown",
        "path": "manuscript/phaseC/tables/supplementary/formatted/supp_table1_phaseC_results_package_index_formatted_v1.md",
        "type": "supplementary_table",
    },
    {
        "item": "Supplementary Table 1 caption",
        "path": "manuscript/phaseC/captions/supp_table1_caption_v1.md",
        "type": "caption",
    },
    {
        "item": "Supplementary Table 2A TSV",
        "path": "manuscript/phaseC/tables/supplementary/formatted/supp_table2A_allowed_phaseC_manuscript_statements_formatted_v1.tsv",
        "type": "supplementary_table",
    },
    {
        "item": "Supplementary Table 2A Markdown",
        "path": "manuscript/phaseC/tables/supplementary/formatted/supp_table2A_allowed_phaseC_manuscript_statements_formatted_v1.md",
        "type": "supplementary_table",
    },
    {
        "item": "Supplementary Table 2B TSV",
        "path": "manuscript/phaseC/tables/supplementary/formatted/supp_table2B_prohibited_phaseC_manuscript_claims_formatted_v1.tsv",
        "type": "supplementary_table",
    },
    {
        "item": "Supplementary Table 2B Markdown",
        "path": "manuscript/phaseC/tables/supplementary/formatted/supp_table2B_prohibited_phaseC_manuscript_claims_formatted_v1.md",
        "type": "supplementary_table",
    },
    {
        "item": "Supplementary Table 2 caption",
        "path": "manuscript/phaseC/captions/supp_table2_caption_v1.md",
        "type": "caption",
    },
]

status_rows = []
missing = []
empty = []
trailing_svg_whitespace = []

for entry in required_files:
    path = Path(entry["path"])
    exists = path.exists()
    size = path.stat().st_size if exists else 0

    if not exists:
        missing.append(str(path))
    elif size == 0:
        empty.append(str(path))

    trailing_ws = "NOT_APPLICABLE"
    if exists and path.suffix.lower() == ".svg":
        bad_lines = []
        with path.open(encoding="utf-8") as handle:
            for i, line in enumerate(handle, start=1):
                if line.rstrip("\n").endswith((" ", "\t")):
                    bad_lines.append(i)
        trailing_ws = "PASS" if not bad_lines else "FAIL"
        if bad_lines:
            trailing_svg_whitespace.append(f"{path}: {len(bad_lines)} lines")

    status_rows.append(
        {
            "item": entry["item"],
            "type": entry["type"],
            "path": str(path),
            "exists": "YES" if exists else "NO",
            "size_bytes": size,
            "nonempty": "YES" if size > 0 else "NO",
            "svg_trailing_whitespace": trailing_ws,
        }
    )

status_path = AUDIT_DIR / "phaseC_display_package_file_status_v1.tsv"
with status_path.open("w", newline="", encoding="utf-8") as handle:
    fieldnames = ["item", "type", "path", "exists", "size_bytes", "nonempty", "svg_trailing_whitespace"]
    writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(status_rows)

summary_rows = [
    {"metric": "required_files", "value": len(required_files)},
    {"metric": "missing_files", "value": len(missing)},
    {"metric": "empty_files", "value": len(empty)},
    {"metric": "svg_trailing_whitespace_failures", "value": len(trailing_svg_whitespace)},
    {"metric": "figures_svg_expected", "value": 3},
    {"metric": "figures_png_expected", "value": 3},
    {"metric": "main_tables_expected", "value": 4},
    {"metric": "supplementary_tables_expected", "value": 6},
    {"metric": "captions_expected", "value": 8},
    {"metric": "new_analysis_performed", "value": "NO"},
    {"metric": "results_prose_written", "value": "NO"},
    {"metric": "selection_or_population_claim_authorised", "value": "NO"},
    {"metric": "broad_biological_interpretation_authorised", "value": "NO"},
]

summary_path = AUDIT_DIR / "phaseC_display_package_audit_summary_v1.tsv"
with summary_path.open("w", newline="", encoding="utf-8") as handle:
    writer = csv.DictWriter(handle, fieldnames=["metric", "value"], delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(summary_rows)

failures = []
if missing:
    failures.append(f"Missing files: {missing}")
if empty:
    failures.append(f"Empty files: {empty}")
if trailing_svg_whitespace:
    failures.append(f"SVG trailing whitespace failures: {trailing_svg_whitespace}")

payload = {
    "status": "PASS" if not failures else "FAIL",
    "operation": "Audit Phase C manuscript display package completeness.",
    "required_files": len(required_files),
    "missing_files": len(missing),
    "empty_files": len(empty),
    "svg_trailing_whitespace_failures": len(trailing_svg_whitespace),
    "new_analysis_performed": "NO",
    "results_prose_written": "NO",
    "selection_or_population_claim_authorised": "NO",
    "broad_biological_interpretation_authorised": "NO",
    "outputs": {
        "file_status": str(status_path),
        "summary": str(summary_path),
    },
    "failures": failures,
    "boundary": (
        "This checkpoint audits completeness of the Phase C manuscript display package only. "
        "It does not perform new analysis, write Results prose, validate function, infer selection, "
        "make population-differentiation claims, assign ancestry effects or authorise broad biological interpretation."
    ),
}

OUTJSON.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
print(json.dumps(payload, indent=2))

if payload["status"] != "PASS":
    raise SystemExit(1)
