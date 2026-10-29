#!/usr/bin/env python3
"""Audit Phase C manuscript-output package readiness.

This checkpoint audits readiness of the Phase C manuscript-output branch after
figures, tables, captions, bounded Results text and bounded Discussion/
limitations text have been generated. It does not perform new analysis,
validate function, infer selection, assign ancestry effects or authorise broad
biological interpretation.
"""

from __future__ import annotations

import csv
import json
import subprocess
from pathlib import Path

OUTJSON = Path("results/qc/phaseC_manuscript_output_readiness_audit_v1.json")
AUDIT_DIR = Path("manuscript/phaseC/audit")
OUTJSON.parent.mkdir(parents=True, exist_ok=True)
AUDIT_DIR.mkdir(parents=True, exist_ok=True)

EXPECTED_BRANCH = "manuscript/phaseC-figures-tables-and-results"

required_files = [
    # Display package audit
    "manuscript/phaseC/audit/phaseC_display_package_audit_summary_v1.tsv",
    "manuscript/phaseC/audit/phaseC_display_package_file_status_v1.tsv",
    "results/qc/phaseC_display_package_audit_v1.json",

    # Figures and captions
    "manuscript/phaseC/figures/svg/figure1_phaseC_genotype_only_qc_workflow_v1.svg",
    "manuscript/phaseC/figures/png/figure1_phaseC_genotype_only_qc_workflow_v1.png",
    "manuscript/phaseC/captions/figure1_caption_v1.md",
    "manuscript/phaseC/figures/svg/figure2_phaseC_candidate_prioritisation_evidence_flow_v1.svg",
    "manuscript/phaseC/figures/png/figure2_phaseC_candidate_prioritisation_evidence_flow_v1.png",
    "manuscript/phaseC/captions/figure2_caption_v1.md",
    "manuscript/phaseC/figures/svg/figure3_phaseC_annotation_supported_evidence_landscape_v1.svg",
    "manuscript/phaseC/figures/png/figure3_phaseC_annotation_supported_evidence_landscape_v1.png",
    "manuscript/phaseC/captions/figure3_caption_v1.md",

    # Main tables and captions
    "manuscript/phaseC/tables/main/formatted/main_table1_annotation_supported_candidate_evidence_summary_formatted_v1.tsv",
    "manuscript/phaseC/tables/main/formatted/main_table1_annotation_supported_candidate_evidence_summary_formatted_v1.md",
    "manuscript/phaseC/captions/main_table1_caption_v1.md",
    "manuscript/phaseC/tables/main/formatted/main_table2_phaseC_candidate_language_boundaries_formatted_v1.tsv",
    "manuscript/phaseC/tables/main/formatted/main_table2_phaseC_candidate_language_boundaries_formatted_v1.md",
    "manuscript/phaseC/captions/main_table2_caption_v1.md",

    # Supplementary tables and captions
    "manuscript/phaseC/tables/supplementary/formatted/supp_table1_phaseC_results_package_index_formatted_v1.tsv",
    "manuscript/phaseC/tables/supplementary/formatted/supp_table1_phaseC_results_package_index_formatted_v1.md",
    "manuscript/phaseC/captions/supp_table1_caption_v1.md",
    "manuscript/phaseC/tables/supplementary/formatted/supp_table2A_allowed_phaseC_manuscript_statements_formatted_v1.tsv",
    "manuscript/phaseC/tables/supplementary/formatted/supp_table2A_allowed_phaseC_manuscript_statements_formatted_v1.md",
    "manuscript/phaseC/tables/supplementary/formatted/supp_table2B_prohibited_phaseC_manuscript_claims_formatted_v1.tsv",
    "manuscript/phaseC/tables/supplementary/formatted/supp_table2B_prohibited_phaseC_manuscript_claims_formatted_v1.md",
    "manuscript/phaseC/captions/supp_table2_caption_v1.md",

    # Text drafts
    "manuscript/phaseC/results_text/phaseC_bounded_results_text_draft_v1.md",
    "manuscript/phaseC/discussion_text/phaseC_bounded_discussion_limitations_text_draft_v1.md",
    "results/qc/phaseC_bounded_results_text_draft_v1.json",
    "results/qc/phaseC_bounded_discussion_limitations_text_draft_v1.json",
]

required_qc_jsons = [
    "results/qc/phaseC_manuscript_source_data_export_v1.json",
    "results/qc/phaseC_main_tables_formatting_v1.json",
    "results/qc/phaseC_supplementary_tables_formatting_v1.json",
    "results/qc/phaseC_figure1_generation_v1.json",
    "results/qc/phaseC_figure2_generation_v1.json",
    "results/qc/phaseC_figure3_generation_v1.json",
    "results/qc/phaseC_display_package_audit_v1.json",
    "results/qc/phaseC_bounded_results_text_draft_v1.json",
    "results/qc/phaseC_bounded_discussion_limitations_text_draft_v1.json",
]

failures = []

def run_git(args: list[str]) -> str:
    return subprocess.check_output(["git", *args], text=True).strip()

branch = run_git(["branch", "--show-current"])
head = run_git(["rev-parse", "--short", "HEAD"])
status_short = subprocess.check_output(
    ["git", "status", "--short", "--untracked-files=all"],
    text=True,
)

file_rows = []
missing = []
empty = []

for file_path in required_files:
    path = Path(file_path)
    exists = path.exists()
    size = path.stat().st_size if exists else 0
    if not exists:
        missing.append(file_path)
    elif size == 0:
        empty.append(file_path)
    file_rows.append(
        {
            "file": file_path,
            "exists": "YES" if exists else "NO",
            "size_bytes": size,
            "nonempty": "YES" if size > 0 else "NO",
        }
    )

qc_rows = []
qc_failures = []
for file_path in required_qc_jsons:
    path = Path(file_path)
    exists = path.exists()
    status = "MISSING"
    if exists:
        try:
            payload = json.loads(path.read_text(encoding="utf-8"))
            status = str(payload.get("status", "NO_STATUS"))
            if status != "PASS":
                qc_failures.append(f"{file_path}: {status}")
        except Exception as exc:
            status = f"READ_ERROR:{exc}"
            qc_failures.append(f"{file_path}: {status}")
    else:
        qc_failures.append(f"{file_path}: MISSING")
    qc_rows.append({"qc_json": file_path, "exists": "YES" if exists else "NO", "status": status})

if branch != EXPECTED_BRANCH:
    failures.append(f"Unexpected branch: {branch}; expected {EXPECTED_BRANCH}")
if missing:
    failures.append(f"Missing required files: {missing}")
if empty:
    failures.append(f"Empty required files: {empty}")
if qc_failures:
    failures.append(f"QC JSON failures: {qc_failures}")
if status_short.strip():
    failures.append("Working tree is not clean at audit time; commit audit outputs after review.")

# The working tree will become non-clean after this audit writes files, so we record pre-audit state.
pre_audit_working_tree_clean = "YES" if not status_short.strip() else "NO"

file_status_path = AUDIT_DIR / "phaseC_manuscript_output_readiness_file_status_v1.tsv"
with file_status_path.open("w", newline="", encoding="utf-8") as handle:
    writer = csv.DictWriter(handle, fieldnames=["file", "exists", "size_bytes", "nonempty"], delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(file_rows)

qc_status_path = AUDIT_DIR / "phaseC_manuscript_output_readiness_qc_status_v1.tsv"
with qc_status_path.open("w", newline="", encoding="utf-8") as handle:
    writer = csv.DictWriter(handle, fieldnames=["qc_json", "exists", "status"], delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(qc_rows)

summary_rows = [
    {"metric": "branch", "value": branch},
    {"metric": "head_commit", "value": head},
    {"metric": "expected_branch", "value": EXPECTED_BRANCH},
    {"metric": "required_files", "value": len(required_files)},
    {"metric": "missing_files", "value": len(missing)},
    {"metric": "empty_files", "value": len(empty)},
    {"metric": "required_qc_jsons", "value": len(required_qc_jsons)},
    {"metric": "qc_json_failures", "value": len(qc_failures)},
    {"metric": "pre_audit_working_tree_clean", "value": pre_audit_working_tree_clean},
    {"metric": "figures_ready", "value": "YES" if not missing and not empty else "CHECK"},
    {"metric": "tables_ready", "value": "YES" if not missing and not empty else "CHECK"},
    {"metric": "bounded_results_ready", "value": "YES" if Path("manuscript/phaseC/results_text/phaseC_bounded_results_text_draft_v1.md").exists() else "NO"},
    {"metric": "bounded_discussion_ready", "value": "YES" if Path("manuscript/phaseC/discussion_text/phaseC_bounded_discussion_limitations_text_draft_v1.md").exists() else "NO"},
    {"metric": "new_analysis_performed", "value": "NO"},
    {"metric": "selection_or_population_claim_authorised", "value": "NO"},
    {"metric": "broad_biological_interpretation_authorised", "value": "NO"},
]

summary_path = AUDIT_DIR / "phaseC_manuscript_output_readiness_summary_v1.tsv"
with summary_path.open("w", newline="", encoding="utf-8") as handle:
    writer = csv.DictWriter(handle, fieldnames=["metric", "value"], delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(summary_rows)

# Do not fail solely because pre-audit tree had expected untracked audit outputs.
hard_failures = [f for f in failures if not f.startswith("Working tree is not clean")]

payload = {
    "status": "PASS" if not hard_failures else "FAIL",
    "operation": "Audit Phase C manuscript-output package readiness.",
    "branch": branch,
    "head_commit": head,
    "required_files": len(required_files),
    "missing_files": len(missing),
    "empty_files": len(empty),
    "required_qc_jsons": len(required_qc_jsons),
    "qc_json_failures": len(qc_failures),
    "pre_audit_working_tree_clean": pre_audit_working_tree_clean,
    "new_analysis_performed": "NO",
    "selection_or_population_claim_authorised": "NO",
    "broad_biological_interpretation_authorised": "NO",
    "outputs": {
        "file_status": str(file_status_path),
        "qc_status": str(qc_status_path),
        "summary": str(summary_path),
    },
    "failures": failures,
    "hard_failures": hard_failures,
    "boundary": (
        "This checkpoint audits readiness of the Phase C manuscript-output package only. It does not perform new analysis, "
        "validate function, infer selection, make population-differentiation claims, assign ancestry effects or authorise broad biological interpretation."
    ),
}

OUTJSON.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
print(json.dumps(payload, indent=2))

if payload["status"] != "PASS":
    raise SystemExit(1)
