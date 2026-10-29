#!/usr/bin/env python3
"""Export Phase C manuscript display source data.

This checkpoint copies/curates committed Phase C source tables into the
manuscript-output branch display package. It does not generate final figures,
write manuscript prose, perform new analysis, infer selection, validate function
or authorise broad biological interpretation.
"""

from __future__ import annotations

import csv
import json
import shutil
from pathlib import Path

OUTBASE = Path("manuscript/phaseC")
FIG_SOURCE = OUTBASE / "figures" / "source_data"
TABLE_MAIN = OUTBASE / "tables" / "main"
TABLE_SUPP = OUTBASE / "tables" / "supplementary"
AUDIT = OUTBASE / "audit"

OUTJSON = Path("results/qc/phaseC_manuscript_source_data_export_v1.json")

for path in [FIG_SOURCE, TABLE_MAIN, TABLE_SUPP, AUDIT, OUTJSON.parent]:
    path.mkdir(parents=True, exist_ok=True)

SOURCE_MAP = [
    {
        "display_item": "FIGURE_1",
        "source": "11_results/phaseC/genotype_only_qc/manuscript_readiness_gate/phaseC_manuscript_readiness_gate_v1.tsv",
        "destination": FIG_SOURCE / "figure1_phaseC_workflow_gate_source_v1.tsv",
        "purpose": "Workflow and audit-control figure source.",
    },
    {
        "display_item": "FIGURE_2",
        "source": "11_results/phaseC/genotype_only_qc/candidate_prioritisation_applied/phaseC_candidate_prioritisation_claim_label_counts_v1.tsv",
        "destination": FIG_SOURCE / "figure2_candidate_prioritisation_claim_label_counts_v1.tsv",
        "purpose": "Candidate prioritisation flow/count source.",
    },
    {
        "display_item": "FIGURE_2",
        "source": "11_results/phaseC/genotype_only_qc/annotation_supported_evidence_review/combined_summary/phaseC_annotation_supported_overall_evidence_grade_counts_v1.tsv",
        "destination": FIG_SOURCE / "figure2_annotation_supported_evidence_grade_counts_v1.tsv",
        "purpose": "Evidence-review flow/count source.",
    },
    {
        "display_item": "FIGURE_3",
        "source": "11_results/phaseC/genotype_only_qc/annotation_supported_evidence_review/combined_summary/phaseC_annotation_supported_combined_evidence_review_status_v1.tsv",
        "destination": FIG_SOURCE / "figure3_annotation_supported_evidence_landscape_source_v1.tsv",
        "purpose": "Evidence-grade landscape source.",
    },
    {
        "display_item": "TABLE_1",
        "source": "11_results/phaseC/genotype_only_qc/annotation_supported_evidence_review/combined_summary/phaseC_annotation_supported_combined_evidence_review_status_v1.tsv",
        "destination": TABLE_MAIN / "table1_annotation_supported_candidate_evidence_summary_v1.tsv",
        "purpose": "Main Table 1 source table.",
    },
    {
        "display_item": "TABLE_2",
        "source": "11_results/phaseC/genotype_only_qc/manuscript_readiness_gate/phaseC_candidate_manuscript_language_boundaries_v1.tsv",
        "destination": TABLE_MAIN / "table2_phaseC_candidate_language_boundaries_v1.tsv",
        "purpose": "Main Table 2 source table.",
    },
    {
        "display_item": "SUPP_TABLE_1",
        "source": "11_results/phaseC/genotype_only_qc/results_package_index/phaseC_results_package_index_v1.tsv",
        "destination": TABLE_SUPP / "supp_table1_phaseC_results_package_index_v1.tsv",
        "purpose": "Supplementary Table 1 source table.",
    },
    {
        "display_item": "SUPP_TABLE_2A",
        "source": "11_results/phaseC/genotype_only_qc/manuscript_readiness_gate/phaseC_allowed_manuscript_statements_v1.tsv",
        "destination": TABLE_SUPP / "supp_table2A_allowed_phaseC_manuscript_statements_v1.tsv",
        "purpose": "Supplementary Table 2A source table.",
    },
    {
        "display_item": "SUPP_TABLE_2B",
        "source": "11_results/phaseC/genotype_only_qc/manuscript_readiness_gate/phaseC_prohibited_manuscript_claims_v1.tsv",
        "destination": TABLE_SUPP / "supp_table2B_prohibited_phaseC_manuscript_claims_v1.tsv",
        "purpose": "Supplementary Table 2B source table.",
    },
]

failures = []
manifest_rows = []

for item in SOURCE_MAP:
    source = Path(item["source"])
    destination = Path(item["destination"])

    if not source.exists():
        failures.append(f"Missing source: {source}")
        continue

    shutil.copy2(source, destination)

    with destination.open(newline="", encoding="utf-8") as handle:
        row_count = sum(1 for _ in csv.DictReader(handle, delimiter="\t"))

    manifest_rows.append(
        {
            "display_item": item["display_item"],
            "source": str(source),
            "destination": str(destination),
            "purpose": item["purpose"],
            "row_count": row_count,
            "status": "copied",
        }
    )

manifest_path = AUDIT / "phaseC_manuscript_source_data_export_manifest_v1.tsv"
with manifest_path.open("w", newline="", encoding="utf-8") as handle:
    fieldnames = ["display_item", "source", "destination", "purpose", "row_count", "status"]
    writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(manifest_rows)

summary_path = AUDIT / "phaseC_manuscript_source_data_export_summary_v1.tsv"
summary_rows = [
    {"metric": "source_items_declared", "value": len(SOURCE_MAP)},
    {"metric": "source_items_exported", "value": len(manifest_rows)},
    {"metric": "missing_sources", "value": len(failures)},
    {"metric": "final_figures_generated", "value": "NO"},
    {"metric": "manuscript_prose_written", "value": "NO"},
    {"metric": "new_analysis_performed", "value": "NO"},
    {"metric": "selection_or_population_claim_authorised", "value": "NO"},
    {"metric": "broad_biological_interpretation_authorised", "value": "NO"},
]
with summary_path.open("w", newline="", encoding="utf-8") as handle:
    writer = csv.DictWriter(handle, fieldnames=["metric", "value"], delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(summary_rows)

payload = {
    "status": "PASS" if not failures else "FAIL",
    "operation": "Export Phase C manuscript display source data.",
    "source_items_declared": len(SOURCE_MAP),
    "source_items_exported": len(manifest_rows),
    "missing_sources": len(failures),
    "final_figures_generated": "NO",
    "manuscript_prose_written": "NO",
    "new_analysis_performed": "NO",
    "selection_or_population_claim_authorised": "NO",
    "broad_biological_interpretation_authorised": "NO",
    "outputs": {
        "manifest": str(manifest_path),
        "summary": str(summary_path),
    },
    "failures": failures,
    "boundary": (
        "This checkpoint exports manuscript display source data only. It does not generate final figures, "
        "write manuscript prose, perform new analysis, validate function, infer selection, make population-differentiation "
        "claims, assign ancestry effects or authorise broad biological interpretation."
    ),
}
OUTJSON.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
print(json.dumps(payload, indent=2))

if payload["status"] != "PASS":
    raise SystemExit(1)
