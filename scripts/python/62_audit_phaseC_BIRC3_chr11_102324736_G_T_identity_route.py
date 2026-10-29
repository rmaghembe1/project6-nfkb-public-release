#!/usr/bin/env python3
"""Audit identity route for BIRC3 chr11:102324736:G:T before evidence grading.

This checkpoint extracts local VEP and prioritisation context for the BIRC3
annotation-supported candidate. It does not perform external evidence review,
does not assign final evidence grades, and does not make biological claims.
"""

from __future__ import annotations

import csv
import json
from collections import Counter
from pathlib import Path

TARGET = "chr11:102324736:G:T"
GENE = "BIRC3"

WORKBOOK = Path(
    "11_results/phaseC/genotype_only_qc/annotation_supported_evidence_review/"
    "phaseC_annotation_supported_evidence_review_workbook_v1.tsv"
)
APPLIED = Path(
    "11_results/phaseC/genotype_only_qc/candidate_prioritisation_applied/"
    "phaseC_candidate_prioritisation_applied_v1.tsv"
)
ALL_VEP = Path(
    "11_results/phaseC/genotype_only_qc/returned_vep_import_audit/"
    "phaseC_returned_vep_all_transcript_annotations_v1.tsv"
)
CANONICAL_VEP = Path(
    "11_results/phaseC/genotype_only_qc/returned_vep_import_audit/"
    "phaseC_returned_vep_MANE_or_APPRIS_canonical_like_annotations_v1.tsv"
)
LANDSCAPE = Path(
    "11_results/phaseC/genotype_only_qc/vep_consequence_landscape/"
    "phaseC_vep_landscape_per_candidate_technical_summary_v1.tsv"
)
PACKAGE_JSON = Path("results/qc/phaseC_annotation_supported_evidence_review_package_v1.json")

OUTDIR = Path(
    "11_results/phaseC/genotype_only_qc/annotation_supported_evidence_review/"
    "BIRC3_chr11_102324736_G_T"
)
OUTJSON = Path("results/qc/phaseC_BIRC3_chr11_102324736_G_T_identity_route_audit_v1.json")

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


def pick(row: dict[str, str], names: list[str]) -> str:
    for name in names:
        value = row.get(name, "")
        if value not in {"", "."}:
            return value
    return "."


for path in [WORKBOOK, APPLIED, ALL_VEP, CANONICAL_VEP, LANDSCAPE, PACKAGE_JSON]:
    if not path.exists():
        failures.append(f"Required input missing: {path}")

workbook = read_tsv(WORKBOOK) if WORKBOOK.exists() else []
applied = read_tsv(APPLIED) if APPLIED.exists() else []
all_vep = read_tsv(ALL_VEP) if ALL_VEP.exists() else []
canonical_vep = read_tsv(CANONICAL_VEP) if CANONICAL_VEP.exists() else []
landscape = read_tsv(LANDSCAPE) if LANDSCAPE.exists() else []

workbook_target = [row for row in workbook if row.get("allele_key") == TARGET]
applied_target = [row for row in applied if row.get("allele_key") == TARGET]
landscape_target = [row for row in landscape if row.get("allele_key") == TARGET]
all_target = [row for row in all_vep if row.get("Uploaded_variation") == TARGET]
canonical_target = [row for row in canonical_vep if row.get("Uploaded_variation") == TARGET]

if len(workbook_target) != 1:
    failures.append(f"Expected 1 workbook target row; observed {len(workbook_target)}.")
if len(applied_target) != 1:
    failures.append(f"Expected 1 applied-prioritisation target row; observed {len(applied_target)}.")
if len(landscape_target) != 1:
    failures.append(f"Expected 1 VEP landscape target row; observed {len(landscape_target)}.")
if not all_target:
    failures.append("No all-transcript VEP rows found for target.")
if not canonical_target:
    failures.append("No canonical-like VEP rows found for target.")

vep_fields = [
    "Uploaded_variation",
    "Location",
    "Allele",
    "Gene",
    "Feature",
    "Feature_type",
    "Consequence",
    "cDNA_position",
    "CDS_position",
    "Protein_position",
    "Amino_acids",
    "Codons",
    "Existing_variation",
    "IMPACT_from_Extra",
    "SYMBOL_from_Extra",
    "MANE_SELECT_from_Extra",
    "MANE_PLUS_CLINICAL_from_Extra",
    "APPRIS_from_Extra",
    "CANONICAL_from_Extra",
    "HGVSc_from_Extra",
    "HGVSp_from_Extra",
]


def reduce_vep_rows(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    return [{field: row.get(field, ".") for field in vep_fields} for row in rows]


all_extract = OUTDIR / "phaseC_BIRC3_chr11_102324736_G_T_all_transcript_VEP_extract_v1.tsv"
canonical_extract = OUTDIR / "phaseC_BIRC3_chr11_102324736_G_T_canonical_like_VEP_extract_v1.tsv"

write_rows(all_extract, vep_fields, reduce_vep_rows(all_target))
write_rows(canonical_extract, vep_fields, reduce_vep_rows(canonical_target))

existing_variants = sorted(
    {
        row.get("Existing_variation", ".")
        for row in all_target + canonical_target
        if row.get("Existing_variation", ".") not in {"", "."}
    }
)
hgvsc_values = sorted(
    {
        row.get("HGVSc_from_Extra", ".")
        for row in canonical_target
        if row.get("HGVSc_from_Extra", ".") not in {"", "."}
    }
)
hgvsp_values = sorted(
    {
        row.get("HGVSp_from_Extra", ".")
        for row in canonical_target
        if row.get("HGVSp_from_Extra", ".") not in {"", "."}
    }
)

consequence_counts = Counter(row.get("Consequence", ".") for row in canonical_target)
impact_counts = Counter(row.get("IMPACT_from_Extra", ".") for row in canonical_target)

identity_rows = [
    {"field": "allele_key", "value": TARGET},
    {"field": "gene", "value": GENE},
    {"field": "workbook_rows", "value": len(workbook_target)},
    {"field": "applied_prioritisation_rows", "value": len(applied_target)},
    {"field": "landscape_rows", "value": len(landscape_target)},
    {"field": "all_transcript_vep_rows", "value": len(all_target)},
    {"field": "canonical_like_vep_rows", "value": len(canonical_target)},
    {"field": "existing_variation_values", "value": ",".join(existing_variants) if existing_variants else "."},
    {"field": "canonical_like_HGVSc_values", "value": ",".join(hgvsc_values) if hgvsc_values else "."},
    {"field": "canonical_like_HGVSp_values", "value": ",".join(hgvsp_values) if hgvsp_values else "."},
    {"field": "canonical_like_consequence_values", "value": ",".join(sorted(consequence_counts)) if consequence_counts else "."},
    {"field": "canonical_like_impact_values", "value": ",".join(sorted(impact_counts)) if impact_counts else "."},
    {"field": "external_evidence_review_performed", "value": "NO"},
    {"field": "identity_route_status", "value": "ready_for_manual_external_lookup" if not failures else "failed"},
]

identity_path = OUTDIR / "phaseC_BIRC3_chr11_102324736_G_T_identity_route_summary_v1.tsv"
write_rows(identity_path, ["field", "value"], identity_rows)

query_rows = [
    {
        "query_type": "exact_allele_key",
        "query": f'"{TARGET}"',
        "purpose": "Check exact coordinate string route.",
        "status": "not_run",
    },
    {
        "query_type": "existing_variation",
        "query": " OR ".join(f'"{item}"' for item in existing_variants) if existing_variants else ".",
        "purpose": "Check rsID or other exact variation identifiers from local VEP.",
        "status": "not_run",
    },
    {
        "query_type": "hgvsc",
        "query": " OR ".join(f'"{item}"' for item in hgvsc_values) if hgvsc_values else ".",
        "purpose": "Check transcript consequence identifiers from local VEP.",
        "status": "not_run",
    },
    {
        "query_type": "hgvsp",
        "query": " OR ".join(f'"{item}"' for item in hgvsp_values) if hgvsp_values else ".",
        "purpose": "Check protein consequence identifiers from local VEP.",
        "status": "not_run",
    },
    {
        "query_type": "gene_context",
        "query": "BIRC3 cIAP2 NF-kB apoptosis ubiquitin inflammation",
        "purpose": "Gene-level context only; not variant-specific evidence.",
        "status": "not_run",
    },
]

query_path = OUTDIR / "phaseC_BIRC3_chr11_102324736_G_T_external_lookup_queries_v1.tsv"
write_rows(query_path, ["query_type", "query", "purpose", "status"], query_rows)

payload = {
    "status": "PASS" if not failures else "FAIL",
    "operation": "Audit local identity route for BIRC3 chr11:102324736:G:T before evidence grading.",
    "allele_key": TARGET,
    "gene": GENE,
    "workbook_rows": len(workbook_target),
    "applied_prioritisation_rows": len(applied_target),
    "all_transcript_vep_rows": len(all_target),
    "canonical_like_vep_rows": len(canonical_target),
    "existing_variation_values": existing_variants,
    "canonical_like_HGVSc_values": hgvsc_values,
    "canonical_like_HGVSp_values": hgvsp_values,
    "external_evidence_review_performed": "NO",
    "outputs": {
        "identity_route_summary": str(identity_path),
        "all_transcript_VEP_extract": str(all_extract),
        "canonical_like_VEP_extract": str(canonical_extract),
        "external_lookup_queries": str(query_path),
    },
    "failures": failures,
    "boundary": (
        "This checkpoint audits local identity information only. It does not perform external evidence review, "
        "assign evidence grades, validate function, infer selection, make population-differentiation claims, "
        "or authorise biological interpretation."
    ),
}

OUTJSON.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
print(json.dumps(payload, indent=2))

if payload["status"] != "PASS":
    raise SystemExit(1)
