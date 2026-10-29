#!/usr/bin/env python3
"""Audit local identity routes for paired TNFAIP3 annotation-supported candidates.

This checkpoint extracts local VEP and prioritisation context for:
- TNFAIP3 chr6:137874871:A:G
- TNFAIP3 chr6:137875820:A:C

It does not perform external evidence review, does not assign final evidence
grades, and does not make biological, functional, population-genetic or
evolutionary claims.
"""

from __future__ import annotations

import csv
import json
from collections import Counter
from pathlib import Path

TARGETS = ["chr6:137874871:A:G", "chr6:137875820:A:C"]
GENE = "TNFAIP3"

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
    "TNFAIP3_paired_candidates"
)
OUTJSON = Path("results/qc/phaseC_TNFAIP3_paired_identity_route_audit_v1.json")

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


for path in [WORKBOOK, APPLIED, ALL_VEP, CANONICAL_VEP, LANDSCAPE, PACKAGE_JSON]:
    if not path.exists():
        failures.append(f"Required input missing: {path}")

workbook = read_tsv(WORKBOOK) if WORKBOOK.exists() else []
applied = read_tsv(APPLIED) if APPLIED.exists() else []
all_vep = read_tsv(ALL_VEP) if ALL_VEP.exists() else []
canonical_vep = read_tsv(CANONICAL_VEP) if CANONICAL_VEP.exists() else []
landscape = read_tsv(LANDSCAPE) if LANDSCAPE.exists() else []

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


summary_rows = []
query_rows = []
all_extract_rows = []
canonical_extract_rows = []

for target in TARGETS:
    workbook_target = [row for row in workbook if row.get("allele_key") == target]
    applied_target = [row for row in applied if row.get("allele_key") == target]
    landscape_target = [row for row in landscape if row.get("allele_key") == target]
    all_target = [row for row in all_vep if row.get("Uploaded_variation") == target]
    canonical_target = [row for row in canonical_vep if row.get("Uploaded_variation") == target]

    if len(workbook_target) != 1:
        failures.append(f"Expected 1 workbook target row for {target}; observed {len(workbook_target)}.")
    if len(applied_target) != 1:
        failures.append(f"Expected 1 applied-prioritisation target row for {target}; observed {len(applied_target)}.")
    if len(landscape_target) != 1:
        failures.append(f"Expected 1 VEP landscape target row for {target}; observed {len(landscape_target)}.")
    if not all_target:
        failures.append(f"No all-transcript VEP rows found for target {target}.")
    if not canonical_target:
        failures.append(f"No canonical-like VEP rows found for target {target}.")

    for row in reduce_vep_rows(all_target):
        row["target_allele_key"] = target
        all_extract_rows.append(row)
    for row in reduce_vep_rows(canonical_target):
        row["target_allele_key"] = target
        canonical_extract_rows.append(row)

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
    protein_positions = sorted(
        {
            row.get("Protein_position", ".")
            for row in canonical_target
            if row.get("Protein_position", ".") not in {"", "."}
        }
    )
    aa_changes = sorted(
        {
            row.get("Amino_acids", ".")
            for row in canonical_target
            if row.get("Amino_acids", ".") not in {"", "."}
        }
    )
    codon_changes = sorted(
        {
            row.get("Codons", ".")
            for row in canonical_target
            if row.get("Codons", ".") not in {"", "."}
        }
    )

    summary_rows.append(
        {
            "allele_key": target,
            "gene": GENE,
            "workbook_rows": len(workbook_target),
            "applied_prioritisation_rows": len(applied_target),
            "landscape_rows": len(landscape_target),
            "all_transcript_vep_rows": len(all_target),
            "canonical_like_vep_rows": len(canonical_target),
            "existing_variation_values": ",".join(existing_variants) if existing_variants else ".",
            "canonical_like_HGVSc_values": ",".join(hgvsc_values) if hgvsc_values else ".",
            "canonical_like_HGVSp_values": ",".join(hgvsp_values) if hgvsp_values else ".",
            "canonical_like_consequence_values": ",".join(sorted(consequence_counts)) if consequence_counts else ".",
            "canonical_like_impact_values": ",".join(sorted(impact_counts)) if impact_counts else ".",
            "canonical_like_protein_positions": ",".join(protein_positions) if protein_positions else ".",
            "canonical_like_amino_acid_changes": ",".join(aa_changes) if aa_changes else ".",
            "canonical_like_codon_changes": ",".join(codon_changes) if codon_changes else ".",
            "external_evidence_review_performed": "NO",
            "identity_route_status": "ready_for_manual_external_lookup" if not failures else "failed",
        }
    )

    query_rows.extend(
        [
            {
                "allele_key": target,
                "query_type": "exact_allele_key",
                "query": f'"{target}"',
                "purpose": "Check exact coordinate string route.",
                "status": "not_run",
            },
            {
                "allele_key": target,
                "query_type": "existing_variation",
                "query": " OR ".join(f'"{item}"' for item in existing_variants) if existing_variants else ".",
                "purpose": "Check rsID or other exact variation identifiers from local VEP.",
                "status": "not_run",
            },
            {
                "allele_key": target,
                "query_type": "hgvsc",
                "query": " OR ".join(f'"{item}"' for item in hgvsc_values) if hgvsc_values else ".",
                "purpose": "Check transcript consequence identifiers from local VEP.",
                "status": "not_run",
            },
            {
                "allele_key": target,
                "query_type": "hgvsp",
                "query": " OR ".join(f'"{item}"' for item in hgvsp_values) if hgvsp_values else ".",
                "purpose": "Check protein consequence identifiers from local VEP.",
                "status": "not_run",
            },
            {
                "allele_key": target,
                "query_type": "gene_context",
                "query": "TNFAIP3 A20 ubiquitin editing NF-kB inflammation immune regulation",
                "purpose": "Gene-level context only; not variant-specific evidence.",
                "status": "not_run",
            },
        ]
    )

all_extract_path = OUTDIR / "phaseC_TNFAIP3_paired_all_transcript_VEP_extract_v1.tsv"
canonical_extract_path = OUTDIR / "phaseC_TNFAIP3_paired_canonical_like_VEP_extract_v1.tsv"
summary_path = OUTDIR / "phaseC_TNFAIP3_paired_identity_route_summary_v1.tsv"
query_path = OUTDIR / "phaseC_TNFAIP3_paired_external_lookup_queries_v1.tsv"

write_rows(all_extract_path, ["target_allele_key"] + vep_fields, all_extract_rows)
write_rows(canonical_extract_path, ["target_allele_key"] + vep_fields, canonical_extract_rows)
write_rows(
    summary_path,
    [
        "allele_key",
        "gene",
        "workbook_rows",
        "applied_prioritisation_rows",
        "landscape_rows",
        "all_transcript_vep_rows",
        "canonical_like_vep_rows",
        "existing_variation_values",
        "canonical_like_HGVSc_values",
        "canonical_like_HGVSp_values",
        "canonical_like_consequence_values",
        "canonical_like_impact_values",
        "canonical_like_protein_positions",
        "canonical_like_amino_acid_changes",
        "canonical_like_codon_changes",
        "external_evidence_review_performed",
        "identity_route_status",
    ],
    summary_rows,
)
write_rows(query_path, ["allele_key", "query_type", "query", "purpose", "status"], query_rows)

payload = {
    "status": "PASS" if not failures else "FAIL",
    "operation": "Audit local identity routes for paired TNFAIP3 annotation-supported candidates before evidence grading.",
    "targets": TARGETS,
    "gene": GENE,
    "target_count": len(TARGETS),
    "summary_rows": len(summary_rows),
    "all_transcript_extract_rows": len(all_extract_rows),
    "canonical_like_extract_rows": len(canonical_extract_rows),
    "external_lookup_query_rows": len(query_rows),
    "external_evidence_review_performed": "NO",
    "outputs": {
        "identity_route_summary": str(summary_path),
        "all_transcript_VEP_extract": str(all_extract_path),
        "canonical_like_VEP_extract": str(canonical_extract_path),
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
