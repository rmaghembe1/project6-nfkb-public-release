#!/usr/bin/env python3
"""Prepare Phase C figure/table blueprint package.

This checkpoint designs manuscript figure and table blueprints for Phase C.
It maps each planned figure/table to result-package files and interpretation
boundaries. It does not generate final figures, write manuscript prose, perform
new analysis, infer selection, validate function or authorise broad biological
interpretation.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path

RESULTS_INDEX_JSON = Path("results/qc/phaseC_results_package_index_v1.json")
RESULTS_INDEX = Path(
    "11_results/phaseC/genotype_only_qc/results_package_index/"
    "phaseC_results_package_index_v1.tsv"
)
MANUSCRIPT_GATE_JSON = Path("results/qc/phaseC_manuscript_readiness_gate_v1.json")

OUTDIR = Path("11_results/phaseC/genotype_only_qc/figure_table_blueprint")
OUTJSON = Path("results/qc/phaseC_figure_table_blueprint_v1.json")

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


for path in [RESULTS_INDEX_JSON, RESULTS_INDEX, MANUSCRIPT_GATE_JSON]:
    if not path.exists():
        failures.append(f"Required input missing: {path}")

results_payload = json.loads(RESULTS_INDEX_JSON.read_text(encoding="utf-8")) if RESULTS_INDEX_JSON.exists() else {}
gate_payload = json.loads(MANUSCRIPT_GATE_JSON.read_text(encoding="utf-8")) if MANUSCRIPT_GATE_JSON.exists() else {}
index_rows = read_tsv(RESULTS_INDEX) if RESULTS_INDEX.exists() else []

if results_payload.get("status") != "PASS":
    failures.append("Results package index does not have PASS status.")
if gate_payload.get("status") != "PASS":
    failures.append("Manuscript-readiness gate does not have PASS status.")
if gate_payload.get("selection_or_population_claim_authorised") != "NO":
    failures.append("Manuscript-readiness gate does not block selection/population claims.")
if len(index_rows) != 7:
    failures.append(f"Expected 7 result-package rows; observed {len(index_rows)}.")

index_by_topic = {row["result_topic"]: row for row in index_rows}

blueprint_rows = [
    {
        "item_id": "FIGURE_1",
        "item_type": "main_figure",
        "proposed_title": "Phase C genotype-only QC and analysis workflow",
        "result_package_source": "manuscript_readiness_gate",
        "primary_source_file": index_by_topic.get("manuscript_readiness_gate", {}).get("primary_result_file", "."),
        "supporting_source_dir": index_by_topic.get("manuscript_readiness_gate", {}).get("supporting_tables_or_dirs", "."),
        "main_message": "Phase C proceeded through strict genotype-only QC, candidate prioritisation, evidence review and manuscript-readiness gates.",
        "allowed_interpretation": "Workflow and audit-control figure only.",
        "blocked_interpretation": "Do not imply biological discovery, selection, adaptation or population specificity.",
        "status": "blueprint_ready",
    },
    {
        "item_id": "FIGURE_2",
        "item_type": "main_figure",
        "proposed_title": "Candidate prioritisation and evidence-review flow",
        "result_package_source": "candidate_prioritisation;annotation_supported_combined_evidence",
        "primary_source_file": "results/qc/phaseC_candidate_prioritisation_applied_v1.json;results/qc/phaseC_annotation_supported_combined_evidence_summary_v1.json",
        "supporting_source_dir": "11_results/phaseC/genotype_only_qc/candidate_prioritisation_applied;11_results/phaseC/genotype_only_qc/annotation_supported_evidence_review/combined_summary",
        "main_message": "Predeclared prioritisation narrowed technical candidates to annotation-supported candidates for bounded evidence review.",
        "allowed_interpretation": "Technical triage and evidence-tiering summary.",
        "blocked_interpretation": "Do not present ranking as biological importance or causal probability.",
        "status": "blueprint_ready",
    },
    {
        "item_id": "FIGURE_3",
        "item_type": "main_figure",
        "proposed_title": "Evidence-grade landscape of annotation-supported candidates",
        "result_package_source": "annotation_supported_combined_evidence",
        "primary_source_file": index_by_topic.get("annotation_supported_combined_evidence", {}).get("primary_result_file", "."),
        "supporting_source_dir": index_by_topic.get("annotation_supported_combined_evidence", {}).get("supporting_tables_or_dirs", "."),
        "main_message": "TRAF3, BIRC3 and TNFAIP3 candidates differ in evidence grade and functional-context boundaries.",
        "allowed_interpretation": "Evidence-grade comparison with candidate-specific caveats.",
        "blocked_interpretation": "Do not imply all candidates have functional validation or equal biological support.",
        "status": "blueprint_ready",
    },
    {
        "item_id": "TABLE_1",
        "item_type": "main_table",
        "proposed_title": "Annotation-supported candidate evidence summary",
        "result_package_source": "annotation_supported_combined_evidence",
        "primary_source_file": "11_results/phaseC/genotype_only_qc/annotation_supported_evidence_review/combined_summary/phaseC_annotation_supported_combined_evidence_review_status_v1.tsv",
        "supporting_source_dir": index_by_topic.get("annotation_supported_combined_evidence", {}).get("supporting_tables_or_dirs", "."),
        "main_message": "Four annotation-supported candidates are summarised with mapped identity, evidence grade, functional-validation status and prohibited claims.",
        "allowed_interpretation": "Bounded evidence-summary table.",
        "blocked_interpretation": "Do not describe as causal, selected, ancestry-determining or independently functionally validated unless explicitly supported.",
        "status": "blueprint_ready",
    },
    {
        "item_id": "TABLE_2",
        "item_type": "main_table",
        "proposed_title": "Manuscript language boundaries for Phase C",
        "result_package_source": "manuscript_readiness_gate",
        "primary_source_file": "11_results/phaseC/genotype_only_qc/manuscript_readiness_gate/phaseC_candidate_manuscript_language_boundaries_v1.tsv",
        "supporting_source_dir": index_by_topic.get("manuscript_readiness_gate", {}).get("supporting_tables_or_dirs", "."),
        "main_message": "Each candidate has allowed language, prohibited language and interpretation boundaries.",
        "allowed_interpretation": "Claim-control and transparency table.",
        "blocked_interpretation": "No overclaiming beyond readiness gate.",
        "status": "blueprint_ready",
    },
    {
        "item_id": "SUPP_TABLE_1",
        "item_type": "supplementary_table",
        "proposed_title": "Phase C results package index",
        "result_package_source": "results_package_index",
        "primary_source_file": "11_results/phaseC/genotype_only_qc/results_package_index/phaseC_results_package_index_v1.tsv",
        "supporting_source_dir": "11_results/phaseC/genotype_only_qc/results_package_index",
        "main_message": "Maps each manuscript-ready statement to supporting files and audit decisions.",
        "allowed_interpretation": "Reproducibility and traceability table.",
        "blocked_interpretation": "No new analysis or biological interpretation.",
        "status": "blueprint_ready",
    },
    {
        "item_id": "SUPP_TABLE_2",
        "item_type": "supplementary_table",
        "proposed_title": "Allowed and prohibited Phase C manuscript statements",
        "result_package_source": "manuscript_readiness_gate",
        "primary_source_file": "11_results/phaseC/genotype_only_qc/manuscript_readiness_gate/phaseC_allowed_manuscript_statements_v1.tsv;11_results/phaseC/genotype_only_qc/manuscript_readiness_gate/phaseC_prohibited_manuscript_claims_v1.tsv",
        "supporting_source_dir": "11_results/phaseC/genotype_only_qc/manuscript_readiness_gate",
        "main_message": "Defines what can and cannot be stated from Phase C.",
        "allowed_interpretation": "Transparency and limitation-control supplement.",
        "blocked_interpretation": "No causal, adaptive, ancestry-effect or population-specific language.",
        "status": "blueprint_ready",
    },
]

for row in blueprint_rows:
    for path_text in row["primary_source_file"].split(";"):
        if path_text and path_text != "." and not Path(path_text).exists():
            failures.append(f"Missing primary source file for {row['item_id']}: {path_text}")
    for path_text in row["supporting_source_dir"].split(";"):
        if path_text and path_text != "." and not Path(path_text).exists():
            failures.append(f"Missing supporting source dir for {row['item_id']}: {path_text}")

figure_rows = [row for row in blueprint_rows if "figure" in row["item_type"]]
table_rows = [row for row in blueprint_rows if "table" in row["item_type"]]
main_rows = [row for row in blueprint_rows if row["item_type"].startswith("main")]
supp_rows = [row for row in blueprint_rows if row["item_type"].startswith("supplementary")]

blueprint_path = OUTDIR / "phaseC_figure_table_blueprint_v1.tsv"
write_rows(
    blueprint_path,
    [
        "item_id",
        "item_type",
        "proposed_title",
        "result_package_source",
        "primary_source_file",
        "supporting_source_dir",
        "main_message",
        "allowed_interpretation",
        "blocked_interpretation",
        "status",
    ],
    blueprint_rows,
)

summary_rows = [
    {"metric": "blueprint_rows", "value": len(blueprint_rows)},
    {"metric": "main_figure_rows", "value": len(figure_rows)},
    {"metric": "table_rows", "value": len(table_rows)},
    {"metric": "main_text_rows", "value": len(main_rows)},
    {"metric": "supplementary_rows", "value": len(supp_rows)},
    {"metric": "missing_source_paths", "value": len(failures)},
    {"metric": "final_figures_generated", "value": "NO"},
    {"metric": "manuscript_prose_written", "value": "NO"},
    {"metric": "selection_or_population_claim_authorised", "value": "NO"},
    {"metric": "broad_biological_interpretation_authorised", "value": "NO"},
]
summary_path = OUTDIR / "phaseC_figure_table_blueprint_summary_v1.tsv"
write_rows(summary_path, ["metric", "value"], summary_rows)

payload = {
    "status": "PASS" if not failures else "FAIL",
    "operation": "Prepare Phase C figure/table blueprint package.",
    "blueprint_rows": len(blueprint_rows),
    "main_figure_rows": len(figure_rows),
    "table_rows": len(table_rows),
    "main_text_rows": len(main_rows),
    "supplementary_rows": len(supp_rows),
    "missing_source_paths": len(failures),
    "final_figures_generated": "NO",
    "manuscript_prose_written": "NO",
    "selection_or_population_claim_authorised": "NO",
    "broad_biological_interpretation_authorised": "NO",
    "outputs": {
        "figure_table_blueprint": str(blueprint_path),
        "summary": str(summary_path),
    },
    "failures": failures,
    "boundary": (
        "This checkpoint designs Phase C figure/table blueprints only. It does not generate final figures, "
        "write manuscript prose, perform new analysis, validate function, infer selection, make "
        "population-differentiation claims, assign ancestry effects or authorise broad biological interpretation."
    ),
}

OUTJSON.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
print(json.dumps(payload, indent=2))

if payload["status"] != "PASS":
    raise SystemExit(1)
