#!/usr/bin/env python3
"""Audit final Phase C branch readiness.

This checkpoint verifies that the Phase C genotype-only QC and analysis branch
has completed its technical, evidence-review, manuscript-readiness, results-index
and figure/table-blueprint gates. It does not perform new analysis, write
manuscript prose, generate figures, infer selection, validate function or
authorise broad biological interpretation.
"""

from __future__ import annotations

import csv
import json
import subprocess
from pathlib import Path

REQUIRED_QC_JSONS = {
    "strict_allele_frequency_qc": "results/qc/phaseC_strict_genotype_only_allele_count_frequency_qc_v1.json",
    "strict_group_frequency_qc": "results/qc/phaseC_strict_genotype_only_group_frequency_qc_v1.json",
    "candidate_prioritisation_applied": "results/qc/phaseC_candidate_prioritisation_applied_v1.json",
    "top_candidate_audit": "results/qc/phaseC_top_candidate_table_audit_v1.json",
    "combined_evidence_summary": "results/qc/phaseC_annotation_supported_combined_evidence_summary_v1.json",
    "manuscript_readiness_gate": "results/qc/phaseC_manuscript_readiness_gate_v1.json",
    "results_package_index": "results/qc/phaseC_results_package_index_v1.json",
    "figure_table_blueprint": "results/qc/phaseC_figure_table_blueprint_v1.json",
}

REQUIRED_OUTPUT_DIRS = {
    "candidate_prioritisation_applied": "11_results/phaseC/genotype_only_qc/candidate_prioritisation_applied",
    "annotation_supported_evidence_review": "11_results/phaseC/genotype_only_qc/annotation_supported_evidence_review",
    "manuscript_readiness_gate": "11_results/phaseC/genotype_only_qc/manuscript_readiness_gate",
    "results_package_index": "11_results/phaseC/genotype_only_qc/results_package_index",
    "figure_table_blueprint": "11_results/phaseC/genotype_only_qc/figure_table_blueprint",
}

OUTDIR = Path("11_results/phaseC/genotype_only_qc/final_branch_readiness_audit")
OUTJSON = Path("results/qc/phaseC_final_branch_readiness_audit_v1.json")

OUTDIR.mkdir(parents=True, exist_ok=True)
OUTJSON.parent.mkdir(parents=True, exist_ok=True)

failures: list[str] = []


def write_rows(path: Path, fieldnames: list[str], rows: list[dict[str, object]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row.get(key, ".") for key in fieldnames})


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


json_rows = []
for label, path_text in REQUIRED_QC_JSONS.items():
    path = Path(path_text)
    exists = path.exists()
    status = "."
    if exists:
        payload = load_json(path)
        status = payload.get("status", ".")
        if status != "PASS":
            failures.append(f"Required QC JSON is not PASS: {path_text}")
    else:
        failures.append(f"Required QC JSON missing: {path_text}")

    json_rows.append(
        {
            "checkpoint": label,
            "path": path_text,
            "exists": "YES" if exists else "NO",
            "status": status,
        }
    )

dir_rows = []
for label, path_text in REQUIRED_OUTPUT_DIRS.items():
    path = Path(path_text)
    exists = path.exists()
    file_count = sum(1 for item in path.rglob("*") if item.is_file()) if exists else 0
    if not exists:
        failures.append(f"Required output directory missing: {path_text}")
    dir_rows.append(
        {
            "checkpoint": label,
            "path": path_text,
            "exists": "YES" if exists else "NO",
            "file_count": file_count,
        }
    )

# Git cleanliness and protected-artifact staging are external workflow checks,
# but we record current branch/commit for traceability.
branch = subprocess.check_output(["git", "rev-parse", "--abbrev-ref", "HEAD"], text=True).strip()
commit = subprocess.check_output(["git", "rev-parse", "--short", "HEAD"], text=True).strip()
status_short = subprocess.check_output(["git", "status", "--short", "--untracked-files=all"], text=True).strip()

branch_rows = [
    {"field": "branch", "value": branch},
    {"field": "commit", "value": commit},
    {"field": "working_tree_status_before_this_audit", "value": status_short if status_short else "clean"},
]

readiness_rows = [
    {
        "readiness_domain": "technical_qc",
        "status": "READY",
        "allowed_next_step": "Report QC-controlled genotype-only Phase C results with technical caveats.",
        "blocked_next_step": "Do not infer selection, adaptation or biological causality from QC outputs.",
    },
    {
        "readiness_domain": "candidate_prioritisation",
        "status": "READY_WITH_BOUNDARIES",
        "allowed_next_step": "Use prioritisation as a transparent triage framework.",
        "blocked_next_step": "Do not treat prioritisation score as biological importance or causal probability.",
    },
    {
        "readiness_domain": "evidence_review",
        "status": "READY_WITH_BOUNDARIES",
        "allowed_next_step": "Use candidate-specific evidence grades and boundaries.",
        "blocked_next_step": "Do not overstate functional validation, especially for BIRC3 or single TNFAIP3 variants.",
    },
    {
        "readiness_domain": "manuscript_readiness",
        "status": "READY_WITH_BOUNDARIES",
        "allowed_next_step": "Draft Results and Discussion only under the readiness gate.",
        "blocked_next_step": "Do not omit limitations or prohibited-claim controls.",
    },
    {
        "readiness_domain": "figures_tables",
        "status": "BLUEPRINT_READY_NOT_FINAL",
        "allowed_next_step": "Generate final figures/tables from the blueprint.",
        "blocked_next_step": "Do not create visual rankings implying selection, ancestry effect or biological superiority.",
    },
]

json_table = OUTDIR / "phaseC_final_branch_required_qc_json_status_v1.tsv"
dir_table = OUTDIR / "phaseC_final_branch_required_output_dir_status_v1.tsv"
branch_table = OUTDIR / "phaseC_final_branch_git_trace_v1.tsv"
readiness_table = OUTDIR / "phaseC_final_branch_readiness_domains_v1.tsv"

write_rows(json_table, ["checkpoint", "path", "exists", "status"], json_rows)
write_rows(dir_table, ["checkpoint", "path", "exists", "file_count"], dir_rows)
write_rows(branch_table, ["field", "value"], branch_rows)
write_rows(readiness_table, ["readiness_domain", "status", "allowed_next_step", "blocked_next_step"], readiness_rows)

summary_rows = [
    {"metric": "required_qc_jsons", "value": len(REQUIRED_QC_JSONS)},
    {"metric": "required_qc_jsons_pass", "value": sum(row["status"] == "PASS" for row in json_rows)},
    {"metric": "required_output_dirs", "value": len(REQUIRED_OUTPUT_DIRS)},
    {"metric": "required_output_dirs_present", "value": sum(row["exists"] == "YES" for row in dir_rows)},
    {"metric": "readiness_domains", "value": len(readiness_rows)},
    {"metric": "new_analysis_performed", "value": "NO"},
    {"metric": "manuscript_prose_written", "value": "NO"},
    {"metric": "final_figures_generated", "value": "NO"},
    {"metric": "selection_or_population_claim_authorised", "value": "NO"},
    {"metric": "broad_biological_interpretation_authorised", "value": "NO"},
    {"metric": "final_phaseC_branch_readiness", "value": "READY_WITH_BOUNDARIES" if not failures else "NOT_READY"},
]
summary_table = OUTDIR / "phaseC_final_branch_readiness_summary_v1.tsv"
write_rows(summary_table, ["metric", "value"], summary_rows)

payload = {
    "status": "PASS" if not failures else "FAIL",
    "operation": "Audit final Phase C branch readiness.",
    "branch": branch,
    "commit": commit,
    "required_qc_jsons": len(REQUIRED_QC_JSONS),
    "required_qc_jsons_pass": sum(row["status"] == "PASS" for row in json_rows),
    "required_output_dirs": len(REQUIRED_OUTPUT_DIRS),
    "required_output_dirs_present": sum(row["exists"] == "YES" for row in dir_rows),
    "readiness_domains": len(readiness_rows),
    "new_analysis_performed": "NO",
    "manuscript_prose_written": "NO",
    "final_figures_generated": "NO",
    "selection_or_population_claim_authorised": "NO",
    "broad_biological_interpretation_authorised": "NO",
    "final_phaseC_branch_readiness": "READY_WITH_BOUNDARIES" if not failures else "NOT_READY",
    "outputs": {
        "required_qc_json_status": str(json_table),
        "required_output_dir_status": str(dir_table),
        "git_trace": str(branch_table),
        "readiness_domains": str(readiness_table),
        "summary": str(summary_table),
    },
    "failures": failures,
    "boundary": (
        "This checkpoint audits final Phase C branch readiness only. It does not perform new analysis, "
        "write manuscript prose, generate figures, validate function, infer selection, make population-differentiation "
        "claims, assign ancestry effects or authorise broad biological interpretation."
    ),
}

OUTJSON.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
print(json.dumps(payload, indent=2))

if payload["status"] != "PASS":
    raise SystemExit(1)
