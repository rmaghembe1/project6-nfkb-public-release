#!/usr/bin/env python3
"""First-pass evidence review for BIRC3 chr11:102324736:G:T.

This checkpoint reviews external and local identity/context evidence for one
annotation-supported candidate: BIRC3 chr11:102324736:G:T. It separates local
exact-variant identity support from gene-level context. It does not validate
function, infer selection, make population-differentiation claims, assign
ancestry effects, or authorise biological interpretation.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path

IDENTITY_JSON = Path("results/qc/phaseC_BIRC3_chr11_102324736_G_T_identity_route_audit_v1.json")
IDENTITY_SUMMARY = Path(
    "11_results/phaseC/genotype_only_qc/annotation_supported_evidence_review/"
    "BIRC3_chr11_102324736_G_T/"
    "phaseC_BIRC3_chr11_102324736_G_T_identity_route_summary_v1.tsv"
)
CANONICAL_VEP = Path(
    "11_results/phaseC/genotype_only_qc/annotation_supported_evidence_review/"
    "BIRC3_chr11_102324736_G_T/"
    "phaseC_BIRC3_chr11_102324736_G_T_canonical_like_VEP_extract_v1.tsv"
)

OUTDIR = Path(
    "11_results/phaseC/genotype_only_qc/annotation_supported_evidence_review/"
    "BIRC3_chr11_102324736_G_T"
)
OUTJSON = Path("results/qc/phaseC_BIRC3_chr11_102324736_G_T_evidence_review_v1.json")

OUTDIR.mkdir(parents=True, exist_ok=True)
OUTJSON.parent.mkdir(parents=True, exist_ok=True)

TARGET = "chr11:102324736:G:T"
GENE = "BIRC3"
RSID = "rs924735278"

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


for path in [IDENTITY_JSON, IDENTITY_SUMMARY, CANONICAL_VEP]:
    if not path.exists():
        failures.append(f"Required input missing: {path}")

identity_payload = json.loads(IDENTITY_JSON.read_text(encoding="utf-8")) if IDENTITY_JSON.exists() else {}
identity_summary = (
    {row["field"]: row["value"] for row in read_tsv(IDENTITY_SUMMARY)}
    if IDENTITY_SUMMARY.exists()
    else {}
)
canonical_rows = read_tsv(CANONICAL_VEP) if CANONICAL_VEP.exists() else []

if identity_payload.get("status") != "PASS":
    failures.append("BIRC3 identity-route audit JSON does not have PASS status.")
if identity_payload.get("allele_key") != TARGET:
    failures.append("BIRC3 identity-route audit does not match target allele key.")
if identity_summary.get("existing_variation_values") != RSID:
    failures.append(f"Expected local VEP existing variation {RSID}; observed {identity_summary.get('existing_variation_values')}.")
if len(canonical_rows) != 3:
    failures.append(f"Expected 3 canonical-like VEP rows; observed {len(canonical_rows)}.")

# This first-pass review is conservative:
# - exact-variant identity is supported locally by VEP rsID and canonical-like missense context;
# - no direct external exact-variant functional/clinical record is asserted here;
# - gene-level context is supported by official NCBI Gene and UniProt records.
evidence_rows = [
    {
        "evidence_id": "BIRC3_W76L_EV001",
        "allele_key": TARGET,
        "gene": GENE,
        "evidence_scope": "exact_variant_local_identity",
        "source_category": "local_VEP_identity_route",
        "source_name": "Project local VEP canonical-like extract",
        "source_url_or_path": str(CANONICAL_VEP),
        "evidence_observation": (
            "Local VEP canonical-like extract maps BIRC3 chr11:102324736:G:T to "
            "Existing_variation rs924735278 across three canonical-like transcript rows, "
            "with missense_variant consequence, MODERATE impact, protein position 76, "
            "amino-acid change W/L and codon change tGg/tTg."
        ),
        "evidence_grade": "indirect_variant_context",
        "supports_functional_validation": "NO",
        "supports_selection_or_population_claim": "NO",
        "interpretation_boundary": (
            "Local exact-variant identity and consequence context only; not external functional validation."
        ),
    },
    {
        "evidence_id": "BIRC3_W76L_EV002",
        "allele_key": TARGET,
        "gene": GENE,
        "evidence_scope": "gene_level",
        "source_category": "gene_function_database",
        "source_name": "NCBI Gene BIRC3",
        "source_url_or_path": "https://www.ncbi.nlm.nih.gov/gene/330",
        "evidence_observation": (
            "NCBI Gene describes BIRC3 as an IAP-family protein that binds TRAF1/TRAF2, "
            "inhibits apoptosis and contains three baculovirus IAP repeats and a RING finger domain."
        ),
        "evidence_grade": "gene_level_support_only",
        "supports_functional_validation": "NO",
        "supports_selection_or_population_claim": "NO",
        "interpretation_boundary": (
            "Gene-level pathway/protein-domain context only; must not be treated as variant-specific evidence."
        ),
    },
    {
        "evidence_id": "BIRC3_W76L_EV003",
        "allele_key": TARGET,
        "gene": GENE,
        "evidence_scope": "gene_level",
        "source_category": "protein_function_database",
        "source_name": "UniProtKB/Swiss-Prot Q13489",
        "source_url_or_path": "https://www.uniprot.org/uniprotkb/Q13489/entry",
        "evidence_observation": (
            "UniProt describes BIRC3 as a multifunctional protein regulating caspases/apoptosis, "
            "inflammatory signalling, immunity and NF-kappa-B-related signalling."
        ),
        "evidence_grade": "gene_level_support_only",
        "supports_functional_validation": "NO",
        "supports_selection_or_population_claim": "NO",
        "interpretation_boundary": (
            "Gene/protein-level pathway context only; does not validate the W76L-like candidate function."
        ),
    },
]

source_table = OUTDIR / "phaseC_BIRC3_chr11_102324736_G_T_evidence_sources_v1.tsv"
write_rows(
    source_table,
    [
        "evidence_id",
        "allele_key",
        "gene",
        "evidence_scope",
        "source_category",
        "source_name",
        "source_url_or_path",
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
    {"metric": "gene", "value": GENE},
    {"metric": "mapped_rsid_from_local_VEP", "value": RSID},
    {"metric": "local_canonical_like_consequence", "value": identity_summary.get("canonical_like_consequence_values", ".")},
    {"metric": "local_canonical_like_impact", "value": identity_summary.get("canonical_like_impact_values", ".")},
    {"metric": "local_protein_position", "value": "76"},
    {"metric": "local_amino_acid_change", "value": "W/L"},
    {"metric": "external_evidence_items", "value": len(evidence_rows)},
    {"metric": "exact_variant_external_database_items", "value": 0},
    {"metric": "local_exact_variant_identity_items", "value": 1},
    {"metric": "gene_level_evidence_items", "value": 2},
    {"metric": "functional_validation_available", "value": "NO"},
    {"metric": "selection_or_population_claim_supported", "value": "NO"},
    {"metric": "overall_variant_level_evidence_grade", "value": "indirect_variant_context"},
    {"metric": "overall_gene_level_evidence_grade", "value": "gene_level_support_only"},
    {"metric": "overall_review_status", "value": "first_pass_review_complete"},
]
summary_table = OUTDIR / "phaseC_BIRC3_chr11_102324736_G_T_evidence_summary_v1.tsv"
write_rows(summary_table, ["metric", "value"], summary_rows)

review_row = {
    "allele_key": TARGET,
    "gene": GENE,
    "mapped_rsid_from_local_VEP": RSID,
    "mapped_transcript_consequence": "not_confirmed_externally",
    "mapped_protein_consequence": "W76L_like_from_local_VEP",
    "variant_level_evidence_grade": "indirect_variant_context",
    "gene_level_evidence_grade": "gene_level_support_only",
    "regulatory_or_qtl_evidence_grade": "not_reviewed",
    "population_genetic_evidence_grade": "not_reviewed",
    "overall_evidence_grade": "indirect_variant_context",
    "functional_validation_available": "NO",
    "selection_or_population_claim_supported": "NO",
    "review_status": "first_pass_review_complete",
    "allowed_language": (
        "BIRC3 chr11:102324736:G:T is an annotation-supported technical candidate with "
        "local VEP identity support as rs924735278 and W/L missense context."
    ),
    "prohibited_language": (
        "Do not describe as externally validated, functionally validated, causal, adaptive, "
        "selected, population-specific or ancestry-determining."
    ),
    "interpretation_boundary": (
        "Local VEP identity/consequence context and gene-level pathway support are present, "
        "but exact external variant evidence and direct functional validation were not established in this checkpoint."
    ),
}

review_table = OUTDIR / "phaseC_BIRC3_chr11_102324736_G_T_review_status_v1.tsv"
write_rows(
    review_table,
    [
        "allele_key",
        "gene",
        "mapped_rsid_from_local_VEP",
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
    "operation": "First-pass evidence review for BIRC3 chr11:102324736:G:T.",
    "allele_key": TARGET,
    "gene": GENE,
    "mapped_rsid_from_local_VEP": RSID,
    "local_canonical_like_consequence": identity_summary.get("canonical_like_consequence_values", "."),
    "local_canonical_like_impact": identity_summary.get("canonical_like_impact_values", "."),
    "local_protein_position": "76",
    "local_amino_acid_change": "W/L",
    "external_evidence_items": len(evidence_rows),
    "exact_variant_external_database_items": 0,
    "local_exact_variant_identity_items": 1,
    "gene_level_evidence_items": 2,
    "variant_level_evidence_grade": "indirect_variant_context",
    "gene_level_evidence_grade": "gene_level_support_only",
    "overall_evidence_grade": "indirect_variant_context",
    "functional_validation_available": "NO",
    "selection_or_population_claim_supported": "NO",
    "outputs": {
        "evidence_sources": str(source_table),
        "evidence_summary": str(summary_table),
        "review_status": str(review_table),
    },
    "failures": failures,
    "boundary": (
        "This checkpoint performs first-pass evidence review for one candidate. It identifies local VEP "
        "identity/consequence context and gene-level pathway support, but does not establish exact external "
        "variant evidence, validate function, infer selection, make population-differentiation claims, "
        "or authorise biological interpretation."
    ),
}

OUTJSON.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
print(json.dumps(payload, indent=2))

if payload["status"] != "PASS":
    raise SystemExit(1)
