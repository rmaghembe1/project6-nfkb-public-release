#!/usr/bin/env python3
"""First-pass evidence review for TRAF3 chr14:102875712:T:C.

This checkpoint reviews external evidence for one annotation-supported candidate:
TRAF3 chr14:102875712:T:C. It separates exact-variant evidence from gene-level
context and does not infer functional validation, selection, adaptation,
population specificity, ancestry effects or causality.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path

PACKAGE_JSON = Path("results/qc/phaseC_annotation_supported_evidence_review_package_v1.json")
WORKBOOK = Path(
    "11_results/phaseC/genotype_only_qc/annotation_supported_evidence_review/"
    "phaseC_annotation_supported_evidence_review_workbook_v1.tsv"
)

OUTDIR = Path(
    "11_results/phaseC/genotype_only_qc/annotation_supported_evidence_review/"
    "TRAF3_chr14_102875712_T_C"
)
OUTJSON = Path("results/qc/phaseC_TRAF3_chr14_102875712_T_C_evidence_review_v1.json")

OUTDIR.mkdir(parents=True, exist_ok=True)
OUTJSON.parent.mkdir(parents=True, exist_ok=True)

TARGET = "chr14:102875712:T:C"
TARGET_GENE = "TRAF3"

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


for path in [PACKAGE_JSON, WORKBOOK]:
    if not path.exists():
        failures.append(f"Required input missing: {path}")

workbook_rows = read_tsv(WORKBOOK) if WORKBOOK.exists() else []
target_rows = [row for row in workbook_rows if row.get("allele_key") == TARGET]

if len(target_rows) != 1:
    failures.append(f"Expected exactly one target workbook row for {TARGET}; observed {len(target_rows)}.")

# Evidence items are deliberately source-level and conservative.
# URLs are stored for reproducibility; claims remain bounded in review fields.
evidence_rows = [
    {
        "evidence_id": "TRAF3_M129T_EV001",
        "allele_key": TARGET,
        "gene": TARGET_GENE,
        "evidence_scope": "exact_variant",
        "source_category": "curated_variant_database",
        "source_name": "ClinVar RCV003487382.4",
        "source_url": "https://www.ncbi.nlm.nih.gov/clinvar/RCV003487382.4/",
        "evidence_observation": (
            "ClinVar record for NM_145725.3(TRAF3):c.386T>C (p.Met129Thr), "
            "dbSNP rs1131877; germline classification benign, multiple submitters/no conflicts."
        ),
        "evidence_grade": "variant_level_support",
        "supports_functional_validation": "NO",
        "supports_selection_or_population_claim": "NO",
        "interpretation_boundary": (
            "Exact-variant database support for identity and benign clinical assertion; "
            "not functional validation and not evolutionary evidence."
        ),
    },
    {
        "evidence_id": "TRAF3_M129T_EV002",
        "allele_key": TARGET,
        "gene": TARGET_GENE,
        "evidence_scope": "exact_variant",
        "source_category": "curated_variant_database",
        "source_name": "ClinVar RCV001518952.9",
        "source_url": "https://www.ncbi.nlm.nih.gov/clinvar/RCV001518952/",
        "evidence_observation": (
            "ClinVar condition-specific record for NM_145725.3(TRAF3):c.386T>C "
            "(p.Met129Thr) and herpes simplex encephalitis susceptibility; "
            "germline classification benign, single submitter."
        ),
        "evidence_grade": "variant_level_support",
        "supports_functional_validation": "NO",
        "supports_selection_or_population_claim": "NO",
        "interpretation_boundary": (
            "Exact-variant database support for identity and benign assertion in a condition record; "
            "not functional validation and not evidence of disease causality."
        ),
    },
    {
        "evidence_id": "TRAF3_M129T_EV003",
        "allele_key": TARGET,
        "gene": TARGET_GENE,
        "evidence_scope": "exact_variant",
        "source_category": "protein_domain_annotation",
        "source_name": "UniProtKB/Swiss-Prot variant VAR_052149",
        "source_url": "https://web.expasy.org/variant_pages/VAR_052149.html",
        "evidence_observation": (
            "UniProtKB/Swiss-Prot variant page maps TRAF3 p.Met129Thr to rs1131877 "
            "and places the residue at position 129 in the 568-aa canonical protein sequence."
        ),
        "evidence_grade": "variant_level_support",
        "supports_functional_validation": "NO",
        "supports_selection_or_population_claim": "NO",
        "interpretation_boundary": (
            "Exact variant/protein-coordinate support; contextual only and not proof of molecular effect."
        ),
    },
    {
        "evidence_id": "TRAF3_M129T_EV004",
        "allele_key": TARGET,
        "gene": TARGET_GENE,
        "evidence_scope": "gene_level",
        "source_category": "gene_function_database",
        "source_name": "NCBI Gene TRAF3",
        "source_url": "https://www.ncbi.nlm.nih.gov/gene/7187",
        "evidence_observation": (
            "NCBI Gene describes TRAF3 as a TNF receptor-associated factor involved in "
            "TNFR-family signal transduction, CD40 signalling, LTbetaR signalling, NF-kB activation context, "
            "and antiviral response."
        ),
        "evidence_grade": "gene_level_support_only",
        "supports_functional_validation": "NO",
        "supports_selection_or_population_claim": "NO",
        "interpretation_boundary": (
            "Gene-level pathway context only; must not be treated as variant-specific evidence."
        ),
    },
    {
        "evidence_id": "TRAF3_M129T_EV005",
        "allele_key": TARGET,
        "gene": TARGET_GENE,
        "evidence_scope": "gene_level",
        "source_category": "protein_function_database",
        "source_name": "UniProtKB/Swiss-Prot Q13114",
        "source_url": "https://www.uniprot.org/uniprotkb/Q13114/entry",
        "evidence_observation": (
            "UniProt describes TRAF3 as acting in TLR/RLR signalling and as an E3 ubiquitin ligase "
            "involved in Lys-63-linked ubiquitin-chain synthesis on pathway substrates."
        ),
        "evidence_grade": "gene_level_support_only",
        "supports_functional_validation": "NO",
        "supports_selection_or_population_claim": "NO",
        "interpretation_boundary": (
            "Gene/protein-level pathway context only; does not validate p.Met129Thr function."
        ),
    },
]

source_table = OUTDIR / "phaseC_TRAF3_chr14_102875712_T_C_evidence_sources_v1.tsv"
write_rows(
    source_table,
    [
        "evidence_id",
        "allele_key",
        "gene",
        "evidence_scope",
        "source_category",
        "source_name",
        "source_url",
        "evidence_observation",
        "evidence_grade",
        "supports_functional_validation",
        "supports_selection_or_population_claim",
        "interpretation_boundary",
    ],
    evidence_rows,
)

summary_rows = [
    {"metric": "allele_key", "value": TARGET},
    {"metric": "gene", "value": TARGET_GENE},
    {"metric": "mapped_rsid", "value": "rs1131877"},
    {"metric": "mapped_transcript_consequence", "value": "NM_145725.3:c.386T>C"},
    {"metric": "mapped_protein_consequence", "value": "p.Met129Thr"},
    {"metric": "external_evidence_items", "value": len(evidence_rows)},
    {"metric": "exact_variant_evidence_items", "value": sum(row["evidence_scope"] == "exact_variant" for row in evidence_rows)},
    {"metric": "gene_level_evidence_items", "value": sum(row["evidence_scope"] == "gene_level" for row in evidence_rows)},
    {"metric": "functional_validation_available", "value": "NO"},
    {"metric": "selection_or_population_claim_supported", "value": "NO"},
    {"metric": "overall_variant_level_evidence_grade", "value": "variant_level_support"},
    {"metric": "overall_gene_level_evidence_grade", "value": "gene_level_support_only"},
    {"metric": "overall_review_status", "value": "first_pass_review_complete"},
]
summary_table = OUTDIR / "phaseC_TRAF3_chr14_102875712_T_C_evidence_summary_v1.tsv"
write_rows(summary_table, ["metric", "value"], summary_rows)

review_row = {
    "allele_key": TARGET,
    "gene": TARGET_GENE,
    "mapped_rsid": "rs1131877",
    "mapped_transcript_consequence": "NM_145725.3:c.386T>C",
    "mapped_protein_consequence": "p.Met129Thr",
    "variant_level_evidence_grade": "variant_level_support",
    "gene_level_evidence_grade": "gene_level_support_only",
    "regulatory_or_qtl_evidence_grade": "not_reviewed",
    "population_genetic_evidence_grade": "not_reviewed",
    "overall_evidence_grade": "variant_level_support",
    "functional_validation_available": "NO",
    "selection_or_population_claim_supported": "NO",
    "review_status": "first_pass_review_complete",
    "allowed_language": (
        "TRAF3 chr14:102875712:T:C is an annotation-supported technical candidate "
        "mapped to rs1131877/p.Met129Thr with exact-variant database support."
    ),
    "prohibited_language": (
        "Do not describe as functionally validated, causal, adaptive, selected, "
        "population-specific or ancestry-determining."
    ),
    "interpretation_boundary": (
        "Exact-variant database/protein-coordinate support and gene-level pathway context are present, "
        "but no direct functional validation or evolutionary inference is established in this checkpoint."
    ),
}

review_table = OUTDIR / "phaseC_TRAF3_chr14_102875712_T_C_review_status_v1.tsv"
write_rows(
    review_table,
    [
        "allele_key",
        "gene",
        "mapped_rsid",
        "mapped_transcript_consequence",
        "mapped_protein_consequence",
        "variant_level_evidence_grade",
        "gene_level_evidence_grade",
        "regulatory_or_qtl_evidence_grade",
        "population_genetic_evidence_grade",
        "overall_evidence_grade",
        "functional_validation_available",
        "selection_or_population_claim_supported",
        "review_status",
        "allowed_language",
        "prohibited_language",
        "interpretation_boundary",
    ],
    [review_row],
)

payload = {
    "status": "PASS" if not failures else "FAIL",
    "operation": "First-pass evidence review for TRAF3 chr14:102875712:T:C.",
    "allele_key": TARGET,
    "gene": TARGET_GENE,
    "mapped_rsid": "rs1131877",
    "mapped_transcript_consequence": "NM_145725.3:c.386T>C",
    "mapped_protein_consequence": "p.Met129Thr",
    "external_evidence_items": len(evidence_rows),
    "exact_variant_evidence_items": sum(row["evidence_scope"] == "exact_variant" for row in evidence_rows),
    "gene_level_evidence_items": sum(row["evidence_scope"] == "gene_level" for row in evidence_rows),
    "variant_level_evidence_grade": "variant_level_support",
    "gene_level_evidence_grade": "gene_level_support_only",
    "overall_evidence_grade": "variant_level_support",
    "functional_validation_available": "NO",
    "selection_or_population_claim_supported": "NO",
    "outputs": {
        "evidence_sources": str(source_table),
        "evidence_summary": str(summary_table),
        "review_status": str(review_table),
    },
    "failures": failures,
    "boundary": (
        "This checkpoint performs first-pass evidence review for one candidate. It identifies exact-variant "
        "database/protein-coordinate support and gene-level pathway context, but does not validate function, "
        "infer selection, make population-differentiation claims, or authorise biological interpretation."
    ),
}

OUTJSON.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
print(json.dumps(payload, indent=2))

if payload["status"] != "PASS":
    raise SystemExit(1)
