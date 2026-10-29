#!/usr/bin/env python3
"""Format Phase C supplementary manuscript tables.

This checkpoint formats Supplementary Table 1 and Supplementary Table 2A/2B
from already exported manuscript source data. It does not perform new analysis,
write Results prose, generate figures, infer selection, validate function or
authorise broad biological interpretation.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path

SRC_SUPP1 = Path("manuscript/phaseC/tables/supplementary/supp_table1_phaseC_results_package_index_v1.tsv")
SRC_SUPP2A = Path("manuscript/phaseC/tables/supplementary/supp_table2A_allowed_phaseC_manuscript_statements_v1.tsv")
SRC_SUPP2B = Path("manuscript/phaseC/tables/supplementary/supp_table2B_prohibited_phaseC_manuscript_claims_v1.tsv")

OUTDIR = Path("manuscript/phaseC/tables/supplementary/formatted")
CAPDIR = Path("manuscript/phaseC/captions")
AUDITDIR = Path("manuscript/phaseC/audit")
OUTJSON = Path("results/qc/phaseC_supplementary_tables_formatting_v1.json")

for path in [OUTDIR, CAPDIR, AUDITDIR, OUTJSON.parent]:
    path.mkdir(parents=True, exist_ok=True)

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


for path in [SRC_SUPP1, SRC_SUPP2A, SRC_SUPP2B]:
    if not path.exists():
        failures.append(f"Missing source table: {path}")

supp1_rows = read_tsv(SRC_SUPP1) if SRC_SUPP1.exists() else []
supp2a_rows = read_tsv(SRC_SUPP2A) if SRC_SUPP2A.exists() else []
supp2b_rows = read_tsv(SRC_SUPP2B) if SRC_SUPP2B.exists() else []

supp1_fields = [
    "Result package",
    "Topic",
    "Allowed statement ID",
    "Primary result file",
    "Supporting source",
    "Audit decision",
    "Manuscript use",
    "Blocked claims",
]

supp1_formatted = []
for row in supp1_rows:
    supp1_formatted.append(
        {
            "Result package": row.get("result_package_id", "."),
            "Topic": row.get("result_topic", "."),
            "Allowed statement ID": row.get("allowed_statement_id", "."),
            "Primary result file": row.get("primary_result_file", "."),
            "Supporting source": row.get("supporting_tables_or_dirs", "."),
            "Audit decision": row.get("audit_decision", "."),
            "Manuscript use": row.get("manuscript_use", "."),
            "Blocked claims": row.get("blocked_claims", "."),
        }
    )

supp2a_fields = [
    "Statement ID",
    "Class",
    "Allowed statement",
    "Required caveat",
    "Source checkpoint",
]

supp2a_formatted = []
for row in supp2a_rows:
    supp2a_formatted.append(
        {
            "Statement ID": row.get("statement_id", "."),
            "Class": row.get("statement_class", "."),
            "Allowed statement": row.get("allowed_statement", "."),
            "Required caveat": row.get("required_caveat", "."),
            "Source checkpoint": row.get("source_checkpoint", "."),
        }
    )

supp2b_fields = [
    "Claim ID",
    "Prohibited claim",
    "Reason",
    "Replacement language",
]

supp2b_formatted = []
for row in supp2b_rows:
    supp2b_formatted.append(
        {
            "Claim ID": row.get("claim_id", "."),
            "Prohibited claim": row.get("prohibited_claim", "."),
            "Reason": row.get("reason", "."),
            "Replacement language": row.get("replacement_language", "."),
        }
    )

if len(supp1_formatted) != 7:
    failures.append(f"Expected 7 Supplementary Table 1 rows; observed {len(supp1_formatted)}.")
if len(supp2a_formatted) != 6:
    failures.append(f"Expected 6 Supplementary Table 2A rows; observed {len(supp2a_formatted)}.")
if len(supp2b_formatted) != 5:
    failures.append(f"Expected 5 Supplementary Table 2B rows; observed {len(supp2b_formatted)}.")

# In allowed statements, prevent affirmative overclaims.
allowed_joined = " ".join(
    " ".join(row.values()) for row in supp2a_formatted
).lower()
for forbidden in [
    "causal variant",
    "selected variant",
    "ancestry-determining",
    "population-specific variant",
]:
    if forbidden in allowed_joined:
        failures.append(f"Forbidden positive wording detected in Supplementary Table 2A: {forbidden}")

outputs = {
    "supp1_tsv": OUTDIR / "supp_table1_phaseC_results_package_index_formatted_v1.tsv",
    "supp1_md": OUTDIR / "supp_table1_phaseC_results_package_index_formatted_v1.md",
    "supp2a_tsv": OUTDIR / "supp_table2A_allowed_phaseC_manuscript_statements_formatted_v1.tsv",
    "supp2a_md": OUTDIR / "supp_table2A_allowed_phaseC_manuscript_statements_formatted_v1.md",
    "supp2b_tsv": OUTDIR / "supp_table2B_prohibited_phaseC_manuscript_claims_formatted_v1.tsv",
    "supp2b_md": OUTDIR / "supp_table2B_prohibited_phaseC_manuscript_claims_formatted_v1.md",
}

write_tsv(outputs["supp1_tsv"], supp1_fields, supp1_formatted)
write_markdown(outputs["supp1_md"], supp1_fields, supp1_formatted)
write_tsv(outputs["supp2a_tsv"], supp2a_fields, supp2a_formatted)
write_markdown(outputs["supp2a_md"], supp2a_fields, supp2a_formatted)
write_tsv(outputs["supp2b_tsv"], supp2b_fields, supp2b_formatted)
write_markdown(outputs["supp2b_md"], supp2b_fields, supp2b_formatted)

caption_supp1 = CAPDIR / "supp_table1_caption_v1.md"
caption_supp1.write_text(
    "Supplementary Table 1. Phase C results package index. The table maps manuscript-ready Phase C result statements to primary result files, supporting source directories, audit decisions, manuscript use and blocked claims.\n",
    encoding="utf-8",
)

caption_supp2 = CAPDIR / "supp_table2_caption_v1.md"
caption_supp2.write_text(
    "Supplementary Table 2. Allowed and prohibited Phase C manuscript statements. Panel A lists allowed bounded statements and required caveats; Panel B lists prohibited claims and replacement language to prevent unsupported causal, adaptive, ancestry-effect, population-specific or functional-validation interpretations.\n",
    encoding="utf-8",
)

summary_path = AUDITDIR / "phaseC_supplementary_tables_formatting_summary_v1.tsv"
summary_rows = [
    {"metric": "supp_table1_rows", "value": len(supp1_formatted)},
    {"metric": "supp_table2A_rows", "value": len(supp2a_formatted)},
    {"metric": "supp_table2B_rows", "value": len(supp2b_formatted)},
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
    "operation": "Format Phase C supplementary manuscript tables.",
    "supp_table1_rows": len(supp1_formatted),
    "supp_table2A_rows": len(supp2a_formatted),
    "supp_table2B_rows": len(supp2b_formatted),
    "captions_written": 2,
    "new_analysis_performed": "NO",
    "results_prose_written": "NO",
    "figures_generated": "NO",
    "selection_or_population_claim_authorised": "NO",
    "broad_biological_interpretation_authorised": "NO",
    "outputs": {
        **{key: str(value) for key, value in outputs.items()},
        "supp_table1_caption": str(caption_supp1),
        "supp_table2_caption": str(caption_supp2),
        "summary": str(summary_path),
    },
    "failures": failures,
    "boundary": (
        "This checkpoint formats supplementary manuscript tables only. It does not perform new analysis, "
        "write Results prose, generate figures, validate function, infer selection, make population-differentiation "
        "claims, assign ancestry effects or authorise broad biological interpretation."
    ),
}

OUTJSON.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
print(json.dumps(payload, indent=2))

if payload["status"] != "PASS":
    raise SystemExit(1)
