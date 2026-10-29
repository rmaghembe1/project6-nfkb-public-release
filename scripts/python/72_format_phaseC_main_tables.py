#!/usr/bin/env python3
"""Format Phase C main manuscript tables.

This checkpoint formats Main Table 1 and Main Table 2 from already exported
manuscript source data. It does not perform new analysis, write Results prose,
generate figures, infer selection, validate function or authorise broad
biological interpretation.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path

SRC_TABLE1 = Path("manuscript/phaseC/tables/main/table1_annotation_supported_candidate_evidence_summary_v1.tsv")
SRC_TABLE2 = Path("manuscript/phaseC/tables/main/table2_phaseC_candidate_language_boundaries_v1.tsv")

OUTDIR = Path("manuscript/phaseC/tables/main/formatted")
CAPDIR = Path("manuscript/phaseC/captions")
OUTJSON = Path("results/qc/phaseC_main_tables_formatting_v1.json")

OUTDIR.mkdir(parents=True, exist_ok=True)
CAPDIR.mkdir(parents=True, exist_ok=True)
OUTJSON.parent.mkdir(parents=True, exist_ok=True)

failures = []


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, fieldnames: list[str], rows: list[dict[str, str]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row.get(key, ".") for key in fieldnames})


def write_markdown(path: Path, fieldnames: list[str], rows: list[dict[str, str]]) -> None:
    with path.open("w", encoding="utf-8") as handle:
        handle.write("| " + " | ".join(fieldnames) + " |\n")
        handle.write("| " + " | ".join(["---"] * len(fieldnames)) + " |\n")
        for row in rows:
            values = [str(row.get(key, ".")).replace("|", "\\|") for key in fieldnames]
            handle.write("| " + " | ".join(values) + " |\n")


for path in [SRC_TABLE1, SRC_TABLE2]:
    if not path.exists():
        failures.append(f"Missing source table: {path}")

table1_rows = read_tsv(SRC_TABLE1) if SRC_TABLE1.exists() else []
table2_rows = read_tsv(SRC_TABLE2) if SRC_TABLE2.exists() else []

table1_fields = [
    "Candidate",
    "Gene",
    "Mapped identity",
    "Evidence grade",
    "Functional validation status",
    "Selection/population claim status",
    "Manuscript boundary",
]

table1_formatted = []
for row in table1_rows:
    table1_formatted.append(
        {
            "Candidate": row.get("allele_key", "."),
            "Gene": row.get("gene", "."),
            "Mapped identity": row.get("mapped_identity", "."),
            "Evidence grade": row.get("overall_evidence_grade", "."),
            "Functional validation status": row.get("functional_validation_available", "."),
            "Selection/population claim status": row.get("selection_or_population_claim_supported", "."),
            "Manuscript boundary": row.get("interpretation_boundary", "."),
        }
    )

table2_fields = [
    "Candidate",
    "Gene",
    "Allowed wording",
    "Prohibited wording",
    "Interpretation boundary",
]

table2_formatted = []
for row in table2_rows:
    table2_formatted.append(
        {
            "Candidate": row.get("allele_key", "."),
            "Gene": row.get("gene", "."),
            "Allowed wording": row.get("allowed_language", "."),
            "Prohibited wording": row.get("prohibited_language", "."),
            "Interpretation boundary": row.get("interpretation_boundary", "."),
        }
    )

if len(table1_formatted) != 4:
    failures.append(f"Expected 4 formatted Table 1 rows; observed {len(table1_formatted)}.")
if len(table2_formatted) != 4:
    failures.append(f"Expected 4 formatted Table 2 rows; observed {len(table2_formatted)}.")

# Scan only affirmative/interpretive columns for forbidden positive claims.
# Prohibited-wording columns are expected to contain blocked terms.
claim_scan_specs = [
    (
        table1_formatted,
        "Table 1",
        [
            "Candidate",
            "Gene",
            "Mapped identity",
            "Evidence grade",
            "Functional validation status",
            "Selection/population claim status",
            "Manuscript boundary",
        ],
    ),
    (
        table2_formatted,
        "Table 2",
        [
            "Candidate",
            "Gene",
            "Allowed wording",
            "Interpretation boundary",
        ],
    ),
]

for rows, label, scan_fields in claim_scan_specs:
    joined = " ".join(
        " ".join(str(row.get(field, "")) for field in scan_fields)
        for row in rows
    ).lower()
    for forbidden in [
        "causal variant",
        "selected variant",
        "ancestry-determining",
        "population-specific variant",
    ]:
        if forbidden in joined:
            failures.append(f"Forbidden positive wording detected in {label}: {forbidden}")

table1_tsv = OUTDIR / "main_table1_annotation_supported_candidate_evidence_summary_formatted_v1.tsv"
table1_md = OUTDIR / "main_table1_annotation_supported_candidate_evidence_summary_formatted_v1.md"
table2_tsv = OUTDIR / "main_table2_phaseC_candidate_language_boundaries_formatted_v1.tsv"
table2_md = OUTDIR / "main_table2_phaseC_candidate_language_boundaries_formatted_v1.md"

write_tsv(table1_tsv, table1_fields, table1_formatted)
write_markdown(table1_md, table1_fields, table1_formatted)
write_tsv(table2_tsv, table2_fields, table2_formatted)
write_markdown(table2_md, table2_fields, table2_formatted)

caption1 = CAPDIR / "main_table1_caption_v1.md"
caption1.write_text(
    "Table 1. Annotation-supported Phase C candidates and evidence-review status. "
    "The table summarises mapped candidate identity, evidence grade, functional-validation status, "
    "selection/population-claim status and manuscript interpretation boundaries. Evidence grades are "
    "used for bounded reporting only and do not imply causality, adaptation, population specificity or "
    "independent functional validation.\n",
    encoding="utf-8",
)

caption2 = CAPDIR / "main_table2_caption_v1.md"
caption2.write_text(
    "Table 2. Candidate-specific manuscript language boundaries for Phase C. "
    "Allowed and prohibited wording is shown for each annotation-supported candidate to preserve the "
    "Phase C manuscript-readiness gate. These boundaries prevent causal, adaptive, ancestry-effect, "
    "population-specific and unsupported functional-validation claims.\n",
    encoding="utf-8",
)

summary_path = Path("manuscript/phaseC/audit/phaseC_main_tables_formatting_summary_v1.tsv")
summary_path.parent.mkdir(parents=True, exist_ok=True)
summary_rows = [
    {"metric": "table1_rows", "value": len(table1_formatted)},
    {"metric": "table2_rows", "value": len(table2_formatted)},
    {"metric": "captions_written", "value": 2},
    {"metric": "new_analysis_performed", "value": "NO"},
    {"metric": "results_prose_written", "value": "NO"},
    {"metric": "figures_generated", "value": "NO"},
    {"metric": "selection_or_population_claim_authorised", "value": "NO"},
    {"metric": "broad_biological_interpretation_authorised", "value": "NO"},
]
with summary_path.open("w", newline="", encoding="utf-8") as handle:
    writer = csv.DictWriter(handle, fieldnames=["metric", "value"], delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(summary_rows)

payload = {
    "status": "PASS" if not failures else "FAIL",
    "operation": "Format Phase C main manuscript tables.",
    "table1_rows": len(table1_formatted),
    "table2_rows": len(table2_formatted),
    "captions_written": 2,
    "new_analysis_performed": "NO",
    "results_prose_written": "NO",
    "figures_generated": "NO",
    "selection_or_population_claim_authorised": "NO",
    "broad_biological_interpretation_authorised": "NO",
    "outputs": {
        "table1_tsv": str(table1_tsv),
        "table1_md": str(table1_md),
        "table2_tsv": str(table2_tsv),
        "table2_md": str(table2_md),
        "table1_caption": str(caption1),
        "table2_caption": str(caption2),
        "summary": str(summary_path),
    },
    "failures": failures,
    "boundary": (
        "This checkpoint formats main manuscript tables only. It does not perform new analysis, "
        "write Results prose, generate figures, validate function, infer selection, make population-differentiation "
        "claims, assign ancestry effects or authorise broad biological interpretation."
    ),
}

OUTJSON.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
print(json.dumps(payload, indent=2))

if payload["status"] != "PASS":
    raise SystemExit(1)
