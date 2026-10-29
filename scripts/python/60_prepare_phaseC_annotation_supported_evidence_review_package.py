#!/usr/bin/env python3
"""Prepare first-pass evidence-review package for annotation-supported Phase C candidates.

This checkpoint prepares structured evidence-review templates for the four
annotation-supported candidates. It does not perform external evidence review,
does not assign evidence grades beyond not_reviewed, and does not make biological
claims.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path

TEMPLATE = Path(
    "11_results/phaseC/genotype_only_qc/evidence_review_gate/"
    "phaseC_annotation_supported_evidence_review_template_v1.tsv"
)
GATE_JSON = Path("results/qc/phaseC_evidence_review_gate_v1.json")
TOP_AUDIT_JSON = Path("results/qc/phaseC_top_candidate_table_audit_v1.json")

OUTDIR = Path("11_results/phaseC/genotype_only_qc/annotation_supported_evidence_review")
OUTJSON = Path("results/qc/phaseC_annotation_supported_evidence_review_package_v1.json")

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


for path in [TEMPLATE, GATE_JSON, TOP_AUDIT_JSON]:
    if not path.exists():
        failures.append(f"Required input missing: {path}")

rows = read_tsv(TEMPLATE) if TEMPLATE.exists() else []

if len(rows) != 4:
    failures.append(f"Expected 4 annotation-supported candidates; observed {len(rows)}.")

def parse_allele_key(allele_key: str) -> dict[str, str]:
    # Expected form: chr14:102875712:T:C
    parts = allele_key.split(":")
    if len(parts) != 4:
        return {"chrom": ".", "pos": ".", "ref": ".", "alt": "."}
    return {"chrom": parts[0], "pos": parts[1], "ref": parts[2], "alt": parts[3]}


review_rows = []
query_rows = []
source_check_rows = []

source_categories = [
    "curated_variant_database",
    "protein_domain_annotation",
    "functional_literature",
    "regulatory_QTL_evidence",
    "population_genetic_follow_up",
]

for row in rows:
    parsed = parse_allele_key(row["allele_key"])
    chrom_no_prefix = parsed["chrom"].replace("chr", "")
    hg38_short = f'{chrom_no_prefix}-{parsed["pos"]}-{parsed["ref"]}-{parsed["alt"]}'
    hg38_colon = f'{parsed["chrom"]}:{parsed["pos"]}:{parsed["ref"]}:{parsed["alt"]}'

    review_rows.append(
        {
            "allele_key": row["allele_key"],
            "gene": row["gene"],
            "chrom": parsed["chrom"],
            "pos_hg38": parsed["pos"],
            "ref": parsed["ref"],
            "alt": parsed["alt"],
            "hg38_short_query_id": hg38_short,
            "hg38_colon_query_id": hg38_colon,
            "claim_boundary_label": row["claim_boundary_label"],
            "technical_priority_rank": row["technical_priority_rank"],
            "technical_priority_score": row["technical_priority_score"],
            "canonical_impact_label": row["canonical_impact_label"],
            "variant_level_evidence_grade": "not_reviewed",
            "gene_level_evidence_grade": "not_reviewed",
            "regulatory_or_qtl_evidence_grade": "not_reviewed",
            "population_genetic_evidence_grade": "not_reviewed",
            "overall_evidence_grade": "not_reviewed",
            "review_status": "pending_external_evidence_review",
            "reviewer_notes": ".",
            "interpretation_boundary": "Package prepared only; no external evidence reviewed or interpreted.",
        }
    )

    query_templates = [
        {
            "query_type": "exact_variant_general",
            "query": f'"{row["allele_key"]}" OR "{hg38_colon}" OR "{hg38_short}"',
        },
        {
            "query_type": "gene_variant_general",
            "query": f'{row["gene"]} "{parsed["pos"]}" "{parsed["ref"]}>{parsed["alt"]}"',
        },
        {
            "query_type": "gene_missense_function",
            "query": f'{row["gene"]} missense variant NF-kB functional study',
        },
        {
            "query_type": "gene_domain_context",
            "query": f'{row["gene"]} protein domain {parsed["pos"]} variant',
        },
        {
            "query_type": "regulatory_qtl_context",
            "query": f'{row["gene"]} eQTL sQTL immune cell stimulation variant',
        },
    ]

    for item in query_templates:
        query_rows.append(
            {
                "allele_key": row["allele_key"],
                "gene": row["gene"],
                "query_type": item["query_type"],
                "query": item["query"],
                "review_status": "not_run",
                "interpretation_boundary": "Search query only; no evidence reviewed.",
            }
        )

    for category in source_categories:
        source_check_rows.append(
            {
                "allele_key": row["allele_key"],
                "gene": row["gene"],
                "source_category": category,
                "source_checked": "NO",
                "evidence_found": "not_reviewed",
                "evidence_level": "not_reviewed",
                "source_reference_or_identifier": ".",
                "review_notes": ".",
                "interpretation_boundary": "Checklist row only; no external evidence reviewed.",
            }
        )

review_path = OUTDIR / "phaseC_annotation_supported_evidence_review_workbook_v1.tsv"
write_rows(
    review_path,
    [
        "allele_key",
        "gene",
        "chrom",
        "pos_hg38",
        "ref",
        "alt",
        "hg38_short_query_id",
        "hg38_colon_query_id",
        "claim_boundary_label",
        "technical_priority_rank",
        "technical_priority_score",
        "canonical_impact_label",
        "variant_level_evidence_grade",
        "gene_level_evidence_grade",
        "regulatory_or_qtl_evidence_grade",
        "population_genetic_evidence_grade",
        "overall_evidence_grade",
        "review_status",
        "reviewer_notes",
        "interpretation_boundary",
    ],
    review_rows,
)

query_path = OUTDIR / "phaseC_annotation_supported_evidence_query_templates_v1.tsv"
write_rows(
    query_path,
    ["allele_key", "gene", "query_type", "query", "review_status", "interpretation_boundary"],
    query_rows,
)

source_check_path = OUTDIR / "phaseC_annotation_supported_evidence_source_checklist_v1.tsv"
write_rows(
    source_check_path,
    [
        "allele_key",
        "gene",
        "source_category",
        "source_checked",
        "evidence_found",
        "evidence_level",
        "source_reference_or_identifier",
        "review_notes",
        "interpretation_boundary",
    ],
    source_check_rows,
)

summary_rows = [
    {"metric": "annotation_supported_candidates", "value": len(review_rows)},
    {"metric": "query_template_rows", "value": len(query_rows)},
    {"metric": "source_checklist_rows", "value": len(source_check_rows)},
    {"metric": "source_categories_per_candidate", "value": len(source_categories)},
    {"metric": "external_evidence_review_performed", "value": "NO"},
]
summary_path = OUTDIR / "phaseC_annotation_supported_evidence_review_package_summary_v1.tsv"
write_rows(summary_path, ["metric", "value"], summary_rows)

payload = {
    "status": "PASS" if not failures else "FAIL",
    "operation": "Prepare first-pass evidence-review package for annotation-supported Phase C candidates.",
    "annotation_supported_candidates": len(review_rows),
    "query_template_rows": len(query_rows),
    "source_checklist_rows": len(source_check_rows),
    "source_categories_per_candidate": len(source_categories),
    "external_evidence_review_performed": "NO",
    "outputs": {
        "review_workbook": str(review_path),
        "query_templates": str(query_path),
        "source_checklist": str(source_check_path),
        "summary": str(summary_path),
    },
    "failures": failures,
    "boundary": (
        "This checkpoint prepares an evidence-review package only. It does not perform external evidence review, "
        "assign evidence grades, validate function, infer selection, make population-differentiation claims, "
        "or authorise biological interpretation."
    ),
}

OUTJSON.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
print(json.dumps(payload, indent=2))

if payload["status"] != "PASS":
    raise SystemExit(1)
