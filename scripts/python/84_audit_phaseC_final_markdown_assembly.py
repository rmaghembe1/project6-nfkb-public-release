#!/usr/bin/env python3
"""Audit final Phase C Markdown manuscript assembly.

This checkpoint audits the Phase C Markdown manuscript assembly after title,
introduction, Methods, Results, Discussion/limitations, figure/table callouts
and boundaries have been assembled. It does not perform new analysis, validate
function, infer selection, assign ancestry effects or authorise broad biological
interpretation.
"""

from __future__ import annotations

import csv
import json
import subprocess
from pathlib import Path

EXPECTED_BRANCH = "manuscript/phaseC-assembly-v1"

BASE = Path("manuscript/phaseC")
ASSEMBLY = BASE / "assembly"
SECTIONS = ASSEMBLY / "sections"
AUDIT = ASSEMBLY / "audit"
OUTJSON = Path("results/qc/phaseC_final_markdown_assembly_audit_v1.json")

for path in [AUDIT, OUTJSON.parent]:
    path.mkdir(parents=True, exist_ok=True)

required_sections = [
    SECTIONS / "00_title_and_scope_v1.md",
    SECTIONS / "01_introduction_phaseC_scope_v1.md",
    SECTIONS / "02_methods_phaseC_assembly_v1.md",
    SECTIONS / "03_results_phaseC_bounded_v1.md",
    SECTIONS / "04_discussion_limitations_phaseC_bounded_v1.md",
    SECTIONS / "05_figures_tables_callouts_v1.md",
    SECTIONS / "06_manuscript_boundaries_v1.md",
]

required_files = [
    ASSEMBLY / "phaseC_manuscript_assembly_master_v1.md",
    *required_sections,
    AUDIT / "phaseC_manuscript_assembly_skeleton_summary_v1.tsv",
    AUDIT / "phaseC_methods_assembly_summary_v1.tsv",
    AUDIT / "phaseC_title_intro_assembly_summary_v1.tsv",
    Path("results/qc/phaseC_manuscript_assembly_skeleton_v1.json"),
    Path("results/qc/phaseC_methods_assembly_v1.json"),
    Path("results/qc/phaseC_title_intro_assembly_v1.json"),
]

required_qc_jsons = [
    Path("results/qc/phaseC_manuscript_assembly_skeleton_v1.json"),
    Path("results/qc/phaseC_methods_assembly_v1.json"),
    Path("results/qc/phaseC_title_intro_assembly_v1.json"),
]

failures = []

def git_output(args: list[str]) -> str:
    return subprocess.check_output(["git", *args], text=True).strip()

branch = git_output(["branch", "--show-current"])
head_commit = git_output(["rev-parse", "--short", "HEAD"])
status_short = subprocess.check_output(
    ["git", "status", "--short", "--untracked-files=all"],
    text=True,
)

if branch != EXPECTED_BRANCH:
    failures.append(f"Unexpected branch: {branch}; expected {EXPECTED_BRANCH}")

file_rows = []
missing = []
empty = []
for path in required_files:
    exists = path.exists()
    size = path.stat().st_size if exists else 0
    if not exists:
        missing.append(str(path))
    elif size == 0:
        empty.append(str(path))
    file_rows.append(
        {
            "file": str(path),
            "exists": "YES" if exists else "NO",
            "size_bytes": size,
            "nonempty": "YES" if size > 0 else "NO",
        }
    )

if missing:
    failures.append(f"Missing required files: {missing}")
if empty:
    failures.append(f"Empty required files: {empty}")

qc_rows = []
qc_failures = []
for path in required_qc_jsons:
    exists = path.exists()
    status = "MISSING"
    if exists:
        try:
            payload = json.loads(path.read_text(encoding="utf-8"))
            status = str(payload.get("status", "NO_STATUS"))
            if status != "PASS":
                qc_failures.append(f"{path}: {status}")
        except Exception as exc:
            status = f"READ_ERROR:{exc}"
            qc_failures.append(f"{path}: {status}")
    else:
        qc_failures.append(f"{path}: MISSING")
    qc_rows.append({"qc_json": str(path), "exists": "YES" if exists else "NO", "status": status})

if qc_failures:
    failures.append(f"QC JSON failures: {qc_failures}")

master_path = ASSEMBLY / "phaseC_manuscript_assembly_master_v1.md"
master_text = master_path.read_text(encoding="utf-8") if master_path.exists() else ""

required_markers = [f"<!-- BEGIN {path.name} -->" for path in required_sections]
missing_markers = [marker for marker in required_markers if marker not in master_text]
if missing_markers:
    failures.append(f"Missing master assembly section markers: {missing_markers}")

unsupported_phrases = [
    "demonstrates selection",
    "proves adaptation",
    "establishes ancestry effects",
    "causal variant",
    "population-specific variant",
    "functionally validated by phase c",
]
unsupported_hits = []
lower_master = master_text.lower()
for phrase in unsupported_phrases:
    if phrase in lower_master:
        unsupported_hits.append(phrase)

if unsupported_hits:
    failures.append(f"Unsupported manuscript phrases detected: {unsupported_hits}")

pre_audit_working_tree_clean = "YES" if not status_short.strip() else "NO"

file_status_path = AUDIT / "phaseC_final_markdown_assembly_file_status_v1.tsv"
with file_status_path.open("w", newline="", encoding="utf-8") as handle:
    writer = csv.DictWriter(handle, fieldnames=["file", "exists", "size_bytes", "nonempty"], delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(file_rows)

qc_status_path = AUDIT / "phaseC_final_markdown_assembly_qc_status_v1.tsv"
with qc_status_path.open("w", newline="", encoding="utf-8") as handle:
    writer = csv.DictWriter(handle, fieldnames=["qc_json", "exists", "status"], delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(qc_rows)

summary_path = AUDIT / "phaseC_final_markdown_assembly_audit_summary_v1.tsv"
summary_rows = [
    {"metric": "branch", "value": branch},
    {"metric": "head_commit", "value": head_commit},
    {"metric": "expected_branch", "value": EXPECTED_BRANCH},
    {"metric": "required_files", "value": len(required_files)},
    {"metric": "missing_files", "value": len(missing)},
    {"metric": "empty_files", "value": len(empty)},
    {"metric": "required_qc_jsons", "value": len(required_qc_jsons)},
    {"metric": "qc_json_failures", "value": len(qc_failures)},
    {"metric": "required_section_markers", "value": len(required_markers)},
    {"metric": "missing_section_markers", "value": len(missing_markers)},
    {"metric": "unsupported_phrase_hits", "value": len(unsupported_hits)},
    {"metric": "pre_audit_working_tree_clean", "value": pre_audit_working_tree_clean},
    {"metric": "new_analysis_performed", "value": "NO"},
    {"metric": "new_biological_claims_added", "value": "NO"},
    {"metric": "selection_or_population_claim_authorised", "value": "NO"},
    {"metric": "broad_biological_interpretation_authorised", "value": "NO"},
]
with summary_path.open("w", newline="", encoding="utf-8") as handle:
    writer = csv.DictWriter(handle, fieldnames=["metric", "value"], delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(summary_rows)

# Newly generated audit files make the tree non-clean after the script; do not treat that as a hard failure.
hard_failures = [f for f in failures if not f.startswith("Working tree")]

payload = {
    "status": "PASS" if not hard_failures else "FAIL",
    "operation": "Audit final Phase C Markdown manuscript assembly.",
    "branch": branch,
    "head_commit": head_commit,
    "required_files": len(required_files),
    "missing_files": len(missing),
    "empty_files": len(empty),
    "required_qc_jsons": len(required_qc_jsons),
    "qc_json_failures": len(qc_failures),
    "required_section_markers": len(required_markers),
    "missing_section_markers": len(missing_markers),
    "unsupported_phrase_hits": len(unsupported_hits),
    "pre_audit_working_tree_clean": pre_audit_working_tree_clean,
    "new_analysis_performed": "NO",
    "new_biological_claims_added": "NO",
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
        "This checkpoint audits the final Phase C Markdown manuscript assembly only. "
        "It does not perform new analysis, validate function, infer selection, assign ancestry effects "
        "or authorise broad biological interpretation."
    ),
}

OUTJSON.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
print(json.dumps(payload, indent=2))

if payload["status"] != "PASS":
    raise SystemExit(1)
