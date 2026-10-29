#!/usr/bin/env python3
"""Prepare Phase D candidate input inventory.

This checkpoint builds a Phase D regulatory/functional evidence input inventory
from committed Phase C manuscript/evidence outputs. It does not query external
resources, perform new biological analysis, validate function, infer selection,
assign ancestry effects or authorise broad biological interpretation.
"""

from __future__ import annotations

import csv
import json
import re
from pathlib import Path

SOURCE = Path("manuscript/phaseC/tables/main/table1_annotation_supported_candidate_evidence_summary_v1.tsv")
OUTDIR = Path("11_results/phaseD/regulatory_functional_evidence/candidate_input_inventory")
AUDITDIR = Path("15_logs_and_audit/phaseD/candidate_input_inventory")
OUTJSON = Path("results/qc/phaseD_candidate_input_inventory_v1.json")

OUTDIR.mkdir(parents=True, exist_ok=True)
AUDITDIR.mkdir(parents=True, exist_ok=True)
OUTJSON.parent.mkdir(parents=True, exist_ok=True)

failures = []

if not SOURCE.exists():
    failures.append(f"Missing Phase C source table: {SOURCE}")

def parse_allele_key(allele_key: str) -> dict[str, str]:
    """Parse allele key formatted as chr:pos:ref:alt."""
    parts = allele_key.split(":")
    if len(parts) != 4:
        return {
            "chromosome": "UNPARSED",
            "position": "UNPARSED",
            "reference_allele": "UNPARSED",
            "alternate_allele": "UNPARSED",
            "genome_build": "UNCONFIRMED_FROM_PHASEC",
            "coordinate_parse_status": "FAIL",
        }
    chrom, pos, ref, alt = parts
    return {
        "chromosome": chrom,
        "position": pos,
        "reference_allele": ref,
        "alternate_allele": alt,
        "genome_build": "UNCONFIRMED_FROM_PHASEC",
        "coordinate_parse_status": "PASS",
    }

rows = []
if SOURCE.exists():
    with SOURCE.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            rows.append(row)

inventory_rows = []
for idx, row in enumerate(rows, start=1):
    allele_key = row.get("allele_key", "").strip()
    parsed = parse_allele_key(allele_key)
    gene = row.get("gene", "").strip()
    mapped_identity = row.get("mapped_identity", "").strip()
    mapped_protein = row.get("mapped_protein_or_hgvs_p", "").strip()
    external_variant_name = row.get("external_variant_name", "").strip()

    has_rsid = "YES" if re.match(r"^rs\d+$", mapped_identity) else "NO"
    has_coordinate = "YES" if parsed["coordinate_parse_status"] == "PASS" else "NO"
    has_protein_label = "YES" if mapped_protein and mapped_protein != "." else "NO"

    if has_coordinate == "YES" and has_rsid == "YES":
        query_readiness = "READY_FOR_DIRECT_AND_RSID_QUERY_WITH_BUILD_CHECK"
    elif has_coordinate == "YES":
        query_readiness = "READY_FOR_COORDINATE_QUERY_WITH_BUILD_CHECK"
    elif has_rsid == "YES":
        query_readiness = "READY_FOR_RSID_QUERY_ONLY"
    else:
        query_readiness = "NOT_READY_NEEDS_IDENTIFIER_REVIEW"

    caution_flags = [
        "genome_build_unconfirmed",
        "do_not_infer_selection_or_ancestry_effect",
        "do_not_treat_regulatory_overlap_as_functional_validation",
    ]

    if gene == "TNFAIP3":
        caution_flags.append("preserve_paired_haplotype_context")
        caution_flags.append("do_not_claim_single_variant_functional_validation")
    if gene == "BIRC3":
        caution_flags.append("indirect_variant_context_requires_extra_identity_check")
    if row.get("functional_validation_available", "").strip() in {"NO", "NO_FOR_THIS_PROJECT_DATA"}:
        caution_flags.append("no_project_level_functional_validation")

    inventory_rows.append(
        {
            "phaseD_candidate_index": idx,
            "allele_key": allele_key,
            "gene": gene,
            "mapped_identity": mapped_identity,
            "external_variant_name": external_variant_name,
            "mapped_protein_or_hgvs_p": mapped_protein,
            "chromosome": parsed["chromosome"],
            "position": parsed["position"],
            "reference_allele": parsed["reference_allele"],
            "alternate_allele": parsed["alternate_allele"],
            "genome_build": parsed["genome_build"],
            "coordinate_parse_status": parsed["coordinate_parse_status"],
            "has_coordinate": has_coordinate,
            "has_rsid": has_rsid,
            "has_protein_label": has_protein_label,
            "phaseC_variant_level_evidence_grade": row.get("variant_level_evidence_grade", "").strip(),
            "phaseC_paired_haplotype_functional_evidence_grade": row.get("paired_haplotype_functional_evidence_grade", "").strip(),
            "phaseC_gene_level_evidence_grade": row.get("gene_level_evidence_grade", "").strip(),
            "phaseC_overall_evidence_grade": row.get("overall_evidence_grade", "").strip(),
            "phaseC_functional_validation_available": row.get("functional_validation_available", "").strip(),
            "phaseC_selection_or_population_claim_supported": row.get("selection_or_population_claim_supported", "").strip(),
            "phaseC_review_status": row.get("review_status", "").strip(),
            "required_phaseD_query_forms": "rsid;coordinate;gene_context;ld_proxy_if_direct_absent",
            "phaseD_query_readiness": query_readiness,
            "phaseD_direct_query_allowed": "YES_WITH_BUILD_CHECK" if has_coordinate == "YES" or has_rsid == "YES" else "NO",
            "phaseD_ld_proxy_query_allowed": "YES_ONLY_IF_DIRECT_ABSENT_AND_RULES_PREDECLARED",
            "phaseD_gene_level_context_allowed": "YES_BUT_NOT_VARIANT_FUNCTION_PROOF",
            "candidate_specific_cautions": ";".join(caution_flags),
            "allowed_phaseC_language": row.get("allowed_language", "").strip(),
            "prohibited_phaseC_language": row.get("prohibited_language", "").strip(),
            "phaseC_interpretation_boundary": row.get("interpretation_boundary", "").strip(),
        }
    )

if not inventory_rows:
    failures.append("No candidate rows were loaded from Phase C source table.")

required_columns = [
    "allele_key",
    "gene",
    "mapped_identity",
    "overall_evidence_grade",
]
if rows:
    missing_cols = [col for col in required_columns if col not in rows[0]]
    if missing_cols:
        failures.append(f"Missing required source columns: {missing_cols}")

inventory_path = OUTDIR / "phaseD_candidate_input_inventory_v1.tsv"
with inventory_path.open("w", encoding="utf-8", newline="") as handle:
    fieldnames = list(inventory_rows[0].keys()) if inventory_rows else [
        "phaseD_candidate_index", "allele_key", "gene", "mapped_identity"
    ]
    writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(inventory_rows)

gene_counts = {}
readiness_counts = {}
build_counts = {}
for row in inventory_rows:
    gene_counts[row["gene"]] = gene_counts.get(row["gene"], 0) + 1
    readiness_counts[row["phaseD_query_readiness"]] = readiness_counts.get(row["phaseD_query_readiness"], 0) + 1
    build_counts[row["genome_build"]] = build_counts.get(row["genome_build"], 0) + 1

gene_counts_path = OUTDIR / "phaseD_candidate_input_inventory_gene_counts_v1.tsv"
with gene_counts_path.open("w", encoding="utf-8", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=["gene", "candidate_count"], delimiter="\t", lineterminator="\n")
    writer.writeheader()
    for gene, count in sorted(gene_counts.items()):
        writer.writerow({"gene": gene, "candidate_count": count})

readiness_counts_path = OUTDIR / "phaseD_candidate_input_inventory_query_readiness_counts_v1.tsv"
with readiness_counts_path.open("w", encoding="utf-8", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=["phaseD_query_readiness", "candidate_count"], delimiter="\t", lineterminator="\n")
    writer.writeheader()
    for label, count in sorted(readiness_counts.items()):
        writer.writerow({"phaseD_query_readiness": label, "candidate_count": count})

summary_path = OUTDIR / "phaseD_candidate_input_inventory_summary_v1.tsv"
summary_rows = [
    {"metric": "source_table", "value": str(SOURCE)},
    {"metric": "candidate_rows", "value": len(inventory_rows)},
    {"metric": "genes_represented", "value": len(gene_counts)},
    {"metric": "coordinates_parsed", "value": sum(1 for r in inventory_rows if r["coordinate_parse_status"] == "PASS")},
    {"metric": "rsid_available", "value": sum(1 for r in inventory_rows if r["has_rsid"] == "YES")},
    {"metric": "genome_build_confirmed", "value": "NO"},
    {"metric": "external_queries_performed", "value": "NO"},
    {"metric": "new_analysis_performed", "value": "NO"},
    {"metric": "new_biological_claims_added", "value": "NO"},
    {"metric": "selection_or_population_claim_authorised", "value": "NO"},
    {"metric": "broad_biological_interpretation_authorised", "value": "NO"},
]
with summary_path.open("w", encoding="utf-8", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=["metric", "value"], delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(summary_rows)

payload = {
    "status": "PASS" if not failures else "FAIL",
    "operation": "Prepare Phase D candidate input inventory from committed Phase C outputs.",
    "source_table": str(SOURCE),
    "candidate_rows": len(inventory_rows),
    "genes_represented": len(gene_counts),
    "coordinates_parsed": sum(1 for r in inventory_rows if r["coordinate_parse_status"] == "PASS"),
    "rsid_available": sum(1 for r in inventory_rows if r["has_rsid"] == "YES"),
    "genome_build_confirmed": "NO",
    "external_queries_performed": "NO",
    "new_analysis_performed": "NO",
    "new_biological_claims_added": "NO",
    "selection_or_population_claim_authorised": "NO",
    "broad_biological_interpretation_authorised": "NO",
    "outputs": {
        "inventory": str(inventory_path),
        "gene_counts": str(gene_counts_path),
        "query_readiness_counts": str(readiness_counts_path),
        "summary": str(summary_path),
    },
    "failures": failures,
    "boundary": (
        "This checkpoint prepares a Phase D candidate input inventory only. "
        "It does not query external regulatory resources, perform new biological analysis, validate function, "
        "infer selection, assign ancestry effects or authorise broad biological interpretation."
    ),
}

OUTJSON.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
print(json.dumps(payload, indent=2))

if payload["status"] != "PASS":
    raise SystemExit(1)
