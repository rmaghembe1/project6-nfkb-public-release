#!/usr/bin/env python3
"""Summarise completed annotation-supported evidence reviews.

This checkpoint combines first-pass evidence-review outputs for the four
annotation-supported candidates: TRAF3, BIRC3 and paired TNFAIP3 candidates.

It summarises evidence status only. It does not perform new evidence review,
validate function, infer selection, make population-differentiation claims,
assign ancestry effects or authorise biological interpretation.
"""

from __future__ import annotations

import csv
import json
from collections import Counter
from pathlib import Path

INPUTS = {
    "TRAF3": Path(
        "11_results/phaseC/genotype_only_qc/annotation_supported_evidence_review/"
        "TRAF3_chr14_102875712_T_C/phaseC_TRAF3_chr14_102875712_T_C_review_status_v1.tsv"
    ),
    "BIRC3": Path(
        "11_results/phaseC/genotype_only_qc/annotation_supported_evidence_review/"
        "BIRC3_chr11_102324736_G_T/phaseC_BIRC3_chr11_102324736_G_T_review_status_v1.tsv"
    ),
    "TNFAIP3": Path(
        "11_results/phaseC/genotype_only_qc/annotation_supported_evidence_review/"
        "TNFAIP3_paired_candidates/phaseC_TNFAIP3_paired_review_status_v1.tsv"
    ),
}

QC_JSONS = {
    "TRAF3": Path("results/qc/phaseC_TRAF3_chr14_102875712_T_C_evidence_review_v1.json"),
    "BIRC3": Path("results/qc/phaseC_BIRC3_chr11_102324736_G_T_evidence_review_v1.json"),
    "TNFAIP3": Path("results/qc/phaseC_TNFAIP3_paired_evidence_review_v1.json"),
}

OUTDIR = Path(
    "11_results/phaseC/genotype_only_qc/annotation_supported_evidence_review/"
    "combined_summary"
)
OUTJSON = Path("results/qc/phaseC_annotation_supported_combined_evidence_summary_v1.json")

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


for label, path in INPUTS.items():
    if not path.exists():
        failures.append(f"Missing review-status input for {label}: {path}")
for label, path in QC_JSONS.items():
    if not path.exists():
        failures.append(f"Missing QC JSON for {label}: {path}")

combined_rows = []

for label, path in INPUTS.items():
    if not path.exists():
        continue

    rows = read_tsv(path)
    for row in rows:
        allele_key = row.get("allele_key", ".")
        gene = row.get("gene", label)

        mapped_identity = (
            row.get("mapped_rsid")
            or row.get("mapped_rsid_from_local_VEP")
            or "."
        )
        external_variant_name = row.get("external_variant_name", ".")
        mapped_protein = (
            row.get("mapped_protein_consequence")
            or row.get("external_hgvs_p")
            or "."
        )

        combined_rows.append(
            {
                "allele_key": allele_key,
                "gene": gene,
                "mapped_identity": mapped_identity,
                "external_variant_name": external_variant_name,
                "mapped_protein_or_hgvs_p": mapped_protein,
                "variant_level_evidence_grade": row.get("variant_level_evidence_grade", "."),
                "paired_haplotype_functional_evidence_grade": row.get(
                    "paired_haplotype_functional_evidence_grade", "."
                ),
                "gene_level_evidence_grade": row.get("gene_level_evidence_grade", "."),
                "overall_evidence_grade": row.get("overall_evidence_grade", "."),
                "functional_validation_available": row.get("functional_validation_available", "."),
                "selection_or_population_claim_supported": row.get(
                    "selection_or_population_claim_supported", "."
                ),
                "review_status": row.get("review_status", "."),
                "allowed_language": row.get("allowed_language", "."),
                "prohibited_language": row.get("prohibited_language", "."),
                "interpretation_boundary": row.get("interpretation_boundary", "."),
            }
        )

if len(combined_rows) != 4:
    failures.append(f"Expected 4 combined annotation-supported candidates; observed {len(combined_rows)}.")

expected_keys = {
    "chr14:102875712:T:C",
    "chr11:102324736:G:T",
    "chr6:137874871:A:G",
    "chr6:137875820:A:C",
}
observed_keys = {row["allele_key"] for row in combined_rows}
if observed_keys != expected_keys:
    failures.append(f"Combined allele-key set mismatch: {sorted(observed_keys)}")

# Ensure all review rows remain bounded.
for row in combined_rows:
    if "Do not describe" not in row["prohibited_language"]:
        failures.append(f"Missing prohibited-language boundary for {row['allele_key']}.")
    boundary_lower = row["interpretation_boundary"].lower()
    if "selection" not in boundary_lower and "functional validation" not in boundary_lower and "exact external variant evidence" not in boundary_lower:
        failures.append(f"Weak interpretation boundary for {row['allele_key']}.")

combined_path = OUTDIR / "phaseC_annotation_supported_combined_evidence_review_status_v1.tsv"
write_rows(
    combined_path,
    [
        "allele_key",
        "gene",
        "mapped_identity",
        "external_variant_name",
        "mapped_protein_or_hgvs_p",
        "variant_level_evidence_grade",
        "paired_haplotype_functional_evidence_grade",
        "gene_level_evidence_grade",
        "overall_evidence_grade",
        "functional_validation_available",
        "selection_or_population_claim_supported",
        "review_status",
        "allowed_language",
        "prohibited_language",
        "interpretation_boundary",
    ],
    combined_rows,
)

grade_counts = Counter(row["overall_evidence_grade"] for row in combined_rows)
gene_counts = Counter(row["gene"] for row in combined_rows)
functional_counts = Counter(row["functional_validation_available"] for row in combined_rows)
selection_counts = Counter(row["selection_or_population_claim_supported"] for row in combined_rows)

def counter_to_rows(counter: Counter, key_name: str, value_name: str) -> list[dict[str, object]]:
    return [{key_name: key, value_name: value} for key, value in sorted(counter.items())]

grade_counts_path = OUTDIR / "phaseC_annotation_supported_overall_evidence_grade_counts_v1.tsv"
gene_counts_path = OUTDIR / "phaseC_annotation_supported_gene_counts_v1.tsv"
functional_counts_path = OUTDIR / "phaseC_annotation_supported_functional_validation_counts_v1.tsv"
selection_counts_path = OUTDIR / "phaseC_annotation_supported_selection_claim_counts_v1.tsv"

write_rows(grade_counts_path, ["overall_evidence_grade", "candidate_count"], counter_to_rows(grade_counts, "overall_evidence_grade", "candidate_count"))
write_rows(gene_counts_path, ["gene", "candidate_count"], counter_to_rows(gene_counts, "gene", "candidate_count"))
write_rows(functional_counts_path, ["functional_validation_available", "candidate_count"], counter_to_rows(functional_counts, "functional_validation_available", "candidate_count"))
write_rows(selection_counts_path, ["selection_or_population_claim_supported", "candidate_count"], counter_to_rows(selection_counts, "selection_or_population_claim_supported", "candidate_count"))

summary_rows = [
    {"metric": "combined_annotation_supported_candidates", "value": len(combined_rows)},
    {"metric": "genes_represented", "value": len(gene_counts)},
    {"metric": "traf3_review_rows", "value": gene_counts.get("TRAF3", 0)},
    {"metric": "birc3_review_rows", "value": gene_counts.get("BIRC3", 0)},
    {"metric": "tnfaip3_review_rows", "value": gene_counts.get("TNFAIP3", 0)},
    {"metric": "new_evidence_review_performed", "value": "NO"},
    {"metric": "selection_or_population_claim_authorised", "value": "NO"},
    {"metric": "broad_biological_interpretation_authorised", "value": "NO"},
]
summary_path = OUTDIR / "phaseC_annotation_supported_combined_evidence_summary_v1.tsv"
write_rows(summary_path, ["metric", "value"], summary_rows)

payload = {
    "status": "PASS" if not failures else "FAIL",
    "operation": "Summarise completed annotation-supported evidence reviews.",
    "combined_annotation_supported_candidates": len(combined_rows),
    "genes_represented": len(gene_counts),
    "overall_evidence_grade_counts": dict(grade_counts),
    "functional_validation_counts": dict(functional_counts),
    "selection_claim_counts": dict(selection_counts),
    "new_evidence_review_performed": "NO",
    "selection_or_population_claim_authorised": "NO",
    "broad_biological_interpretation_authorised": "NO",
    "outputs": {
        "combined_review_status": str(combined_path),
        "summary": str(summary_path),
        "overall_evidence_grade_counts": str(grade_counts_path),
        "gene_counts": str(gene_counts_path),
        "functional_validation_counts": str(functional_counts_path),
        "selection_claim_counts": str(selection_counts_path),
    },
    "failures": failures,
    "boundary": (
        "This checkpoint summarises already completed annotation-supported evidence reviews only. "
        "It does not perform new evidence review, validate function, infer selection, make "
        "population-differentiation claims, assign ancestry effects, or authorise broad biological interpretation."
    ),
}

OUTJSON.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
print(json.dumps(payload, indent=2))

if payload["status"] != "PASS":
    raise SystemExit(1)
