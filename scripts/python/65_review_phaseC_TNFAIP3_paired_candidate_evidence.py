#!/usr/bin/env python3
"""First-pass evidence review for paired TNFAIP3 annotation-supported candidates.

Targets:
- TNFAIP3 chr6:137874871:A:G -> rs376205580, T108A-like local VEP context
- TNFAIP3 chr6:137875820:A:C -> rs141807543, I207L-like local VEP context

This checkpoint separates exact-variant identity evidence, paired/haplotype-level
functional evidence and gene-level TNFAIP3/A20 context. It does not infer
selection from this project data, make population-differentiation claims, assign
ancestry effects, or authorise biological interpretation beyond bounded evidence
review.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path

IDENTITY_JSON = Path("results/qc/phaseC_TNFAIP3_paired_identity_route_audit_v1.json")
IDENTITY_SUMMARY = Path(
    "11_results/phaseC/genotype_only_qc/annotation_supported_evidence_review/"
    "TNFAIP3_paired_candidates/phaseC_TNFAIP3_paired_identity_route_summary_v1.tsv"
)
CANONICAL_VEP = Path(
    "11_results/phaseC/genotype_only_qc/annotation_supported_evidence_review/"
    "TNFAIP3_paired_candidates/phaseC_TNFAIP3_paired_canonical_like_VEP_extract_v1.tsv"
)

OUTDIR = Path(
    "11_results/phaseC/genotype_only_qc/annotation_supported_evidence_review/"
    "TNFAIP3_paired_candidates"
)
OUTJSON = Path("results/qc/phaseC_TNFAIP3_paired_evidence_review_v1.json")

OUTDIR.mkdir(parents=True, exist_ok=True)
OUTJSON.parent.mkdir(parents=True, exist_ok=True)

TARGETS = {
    "chr6:137874871:A:G": {
        "gene": "TNFAIP3",
        "rsid": "rs376205580",
        "local_amino_acid_change": "T/A",
        "local_protein_position": "108",
        "external_variant_name": "T108A",
        "external_hgvs_c": "NM_001270508.2:c.322A>G",
        "external_hgvs_p": "p.(Thr108Ala)",
    },
    "chr6:137875820:A:C": {
        "gene": "TNFAIP3",
        "rsid": "rs141807543",
        "local_amino_acid_change": "I/L",
        "local_protein_position": "207",
        "external_variant_name": "I207L",
        "external_hgvs_c": "NM_001270508.2:c.619A>C",
        "external_hgvs_p": "p.(Ile207Leu)",
    },
}

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
summary_rows_input = read_tsv(IDENTITY_SUMMARY) if IDENTITY_SUMMARY.exists() else []
canonical_rows = read_tsv(CANONICAL_VEP) if CANONICAL_VEP.exists() else []

if identity_payload.get("status") != "PASS":
    failures.append("TNFAIP3 paired identity-route audit JSON does not have PASS status.")
if set(identity_payload.get("targets", [])) != set(TARGETS):
    failures.append("TNFAIP3 identity-route targets do not match expected paired targets.")

summary_by_key = {row["allele_key"]: row for row in summary_rows_input}
for allele_key, expected in TARGETS.items():
    row = summary_by_key.get(allele_key)
    if row is None:
        failures.append(f"Missing identity summary row for {allele_key}.")
        continue
    if expected["rsid"] not in row.get("existing_variation_values", ""):
        failures.append(f"Expected {expected['rsid']} in local VEP existing variation for {allele_key}.")
    if row.get("canonical_like_consequence_values") != "missense_variant":
        failures.append(f"Expected missense_variant for {allele_key}.")
    if row.get("canonical_like_impact_values") != "MODERATE":
        failures.append(f"Expected MODERATE impact for {allele_key}.")
    if row.get("canonical_like_protein_positions") != expected["local_protein_position"]:
        failures.append(f"Unexpected protein position for {allele_key}.")
    if row.get("canonical_like_amino_acid_changes") != expected["local_amino_acid_change"]:
        failures.append(f"Unexpected amino-acid change for {allele_key}.")

# Evidence rows: exact variant identity, paired haplotype-level functional evidence,
# and gene-level context are kept distinct.
evidence_rows = []

for allele_key, info in TARGETS.items():
    evidence_rows.append(
        {
            "evidence_id": f"TNFAIP3_{info['external_variant_name']}_EV001",
            "allele_key": allele_key,
            "gene": info["gene"],
            "evidence_scope": "exact_variant_identity",
            "source_category": "curated_variant_database",
            "source_name": "Infevers variant detail",
            "source_url_or_path": (
                "https://infevers.umai-montpellier.fr/web/detail_mutation.php"
            ),
            "evidence_observation": (
                f"Infevers records {info['external_variant_name']} for TNFAIP3 with "
                f"{info['external_hgvs_c']}, {info['external_hgvs_p']} and rs number {info['rsid']}."
            ),
            "evidence_grade": "variant_level_support",
            "supports_functional_validation": "NO",
            "supports_selection_or_population_claim": "NO",
            "interpretation_boundary": (
                "Exact-variant identity support only; not functional validation by itself."
            ),
        }
    )
    evidence_rows.append(
        {
            "evidence_id": f"TNFAIP3_{info['external_variant_name']}_EV002",
            "allele_key": allele_key,
            "gene": info["gene"],
            "evidence_scope": "exact_variant_local_identity",
            "source_category": "local_VEP_identity_route",
            "source_name": "Project local VEP canonical-like extract",
            "source_url_or_path": str(CANONICAL_VEP),
            "evidence_observation": (
                f"Local VEP maps {allele_key} to {info['rsid']} with missense_variant, "
                f"MODERATE impact, protein position {info['local_protein_position']} and "
                f"amino-acid change {info['local_amino_acid_change']}."
            ),
            "evidence_grade": "variant_level_support",
            "supports_functional_validation": "NO",
            "supports_selection_or_population_claim": "NO",
            "interpretation_boundary": (
                "Local exact-variant identity and consequence context only."
            ),
        }
    )

# Paired/haplotype evidence is listed for both targets but explicitly scoped as paired evidence.
for allele_key, info in TARGETS.items():
    evidence_rows.append(
        {
            "evidence_id": f"TNFAIP3_{info['external_variant_name']}_EV003",
            "allele_key": allele_key,
            "gene": info["gene"],
            "evidence_scope": "paired_haplotype_functional_context",
            "source_category": "primary_functional_literature",
            "source_name": "Zammit et al. 2019 Nature Immunology / PubMed PMID 31534238",
            "source_url_or_path": "https://pubmed.ncbi.nlm.nih.gov/31534238/",
            "evidence_observation": (
                "Primary literature reports the paired A20 T108A;I207L allele as a Denisovan-introgressed "
                "TNFAIP3 allele and describes partial phosphorylation/activation effects for A20 allelic series."
            ),
            "evidence_grade": "functional_validation_available_haplotype_level",
            "supports_functional_validation": "YES_HAPLOTYPE_LEVEL_ONLY",
            "supports_selection_or_population_claim": "NO_FOR_THIS_PROJECT_DATA",
            "interpretation_boundary": (
                "Functional evidence applies to paired/haplotype context; do not treat as independent "
                "single-variant functional validation or as a new project-level selection claim."
            ),
        }
    )

gene_level_rows = [
    {
        "evidence_id": "TNFAIP3_GENE_EV001",
        "allele_key": "paired_TNFAIP3_candidates",
        "gene": "TNFAIP3",
        "evidence_scope": "gene_level",
        "source_category": "protein_function_database",
        "source_name": "UniProtKB/Swiss-Prot P21580",
        "source_url_or_path": "https://www.uniprot.org/uniprotkb/P21580/entry",
        "evidence_observation": (
            "UniProt describes TNFAIP3/A20 as a ubiquitin-editing enzyme with ubiquitin ligase "
            "and deubiquitinase activities involved in immune and inflammatory responses and "
            "termination of NF-kappa-B activity."
        ),
        "evidence_grade": "gene_level_support_only",
        "supports_functional_validation": "NO_VARIANT_SPECIFIC",
        "supports_selection_or_population_claim": "NO",
        "interpretation_boundary": (
            "Gene-level pathway context only; must not be promoted to variant-specific evidence."
        ),
    },
    {
        "evidence_id": "TNFAIP3_GENE_EV002",
        "allele_key": "paired_TNFAIP3_candidates",
        "gene": "TNFAIP3",
        "evidence_scope": "gene_level",
        "source_category": "review_literature",
        "source_name": "Vereecke et al. 2009 PubMed PMID 19643665",
        "source_url_or_path": "https://pubmed.ncbi.nlm.nih.gov/19643665/",
        "evidence_observation": (
            "Review literature describes A20/TNFAIP3 as a key negative feedback regulator of NF-kappa-B signalling."
        ),
        "evidence_grade": "gene_level_support_only",
        "supports_functional_validation": "NO_VARIANT_SPECIFIC",
        "supports_selection_or_population_claim": "NO",
        "interpretation_boundary": (
            "Gene-level pathway context only; not exact-variant functional evidence."
        ),
    },
]

evidence_rows.extend(gene_level_rows)

source_table = OUTDIR / "phaseC_TNFAIP3_paired_evidence_sources_v1.tsv"
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

review_rows = []
for allele_key, info in TARGETS.items():
    review_rows.append(
        {
            "allele_key": allele_key,
            "gene": info["gene"],
            "mapped_rsid": info["rsid"],
            "external_variant_name": info["external_variant_name"],
            "external_hgvs_c": info["external_hgvs_c"],
            "external_hgvs_p": info["external_hgvs_p"],
            "local_protein_position": info["local_protein_position"],
            "local_amino_acid_change": info["local_amino_acid_change"],
            "variant_level_evidence_grade": "variant_level_support",
            "paired_haplotype_functional_evidence_grade": "functional_validation_available_haplotype_level",
            "gene_level_evidence_grade": "gene_level_support_only",
            "regulatory_or_qtl_evidence_grade": "not_reviewed",
            "population_genetic_evidence_grade": "not_reviewed",
            "overall_evidence_grade": "variant_level_support_with_haplotype_level_functional_context",
            "functional_validation_available": "YES_HAPLOTYPE_LEVEL_ONLY",
            "selection_or_population_claim_supported": "NO_FOR_THIS_PROJECT_DATA",
            "review_status": "first_pass_review_complete",
            "allowed_language": (
                f"TNFAIP3 {allele_key} is an annotation-supported technical candidate mapping to "
                f"{info['rsid']} / {info['external_variant_name']}, with exact-variant identity support "
                "and paired T108A;I207L haplotype-level functional context."
            ),
            "prohibited_language": (
                "Do not describe as independently functionally validated as a single variant, causal, "
                "selected in this dataset, population-specific, ancestry-determining or biologically proven by this project."
            ),
            "interpretation_boundary": (
                "Exact-variant identity evidence and paired/haplotype-level functional literature are present, "
                "but this project does not independently infer selection, adaptation, population specificity or single-variant causality."
            ),
        }
    )

review_table = OUTDIR / "phaseC_TNFAIP3_paired_review_status_v1.tsv"
write_rows(
    review_table,
    [
        "allele_key",
        "gene",
        "mapped_rsid",
        "external_variant_name",
        "external_hgvs_c",
        "external_hgvs_p",
        "local_protein_position",
        "local_amino_acid_change",
        "variant_level_evidence_grade",
        "paired_haplotype_functional_evidence_grade",
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
    review_rows,
)

summary_rows = [
    {"metric": "target_candidates_reviewed", "value": len(TARGETS)},
    {"metric": "evidence_items_total", "value": len(evidence_rows)},
    {"metric": "exact_variant_identity_items", "value": 4},
    {"metric": "paired_haplotype_functional_items", "value": 2},
    {"metric": "gene_level_evidence_items", "value": 2},
    {"metric": "functional_validation_available", "value": "YES_HAPLOTYPE_LEVEL_ONLY"},
    {"metric": "selection_or_population_claim_supported", "value": "NO_FOR_THIS_PROJECT_DATA"},
    {"metric": "overall_review_status", "value": "first_pass_review_complete"},
]
summary_table = OUTDIR / "phaseC_TNFAIP3_paired_evidence_summary_v1.tsv"
write_rows(summary_table, ["metric", "value"], summary_rows)

payload = {
    "status": "PASS" if not failures else "FAIL",
    "operation": "First-pass evidence review for paired TNFAIP3 annotation-supported candidates.",
    "targets": list(TARGETS),
    "gene": "TNFAIP3",
    "target_candidates_reviewed": len(TARGETS),
    "evidence_items_total": len(evidence_rows),
    "exact_variant_identity_items": 4,
    "paired_haplotype_functional_items": 2,
    "gene_level_evidence_items": 2,
    "variant_level_evidence_grade": "variant_level_support",
    "paired_haplotype_functional_evidence_grade": "functional_validation_available_haplotype_level",
    "gene_level_evidence_grade": "gene_level_support_only",
    "overall_evidence_grade": "variant_level_support_with_haplotype_level_functional_context",
    "functional_validation_available": "YES_HAPLOTYPE_LEVEL_ONLY",
    "selection_or_population_claim_supported": "NO_FOR_THIS_PROJECT_DATA",
    "outputs": {
        "evidence_sources": str(source_table),
        "evidence_summary": str(summary_table),
        "review_status": str(review_table),
    },
    "failures": failures,
    "boundary": (
        "This checkpoint performs first-pass evidence review for paired TNFAIP3 candidates. It identifies "
        "exact-variant identity support and paired/haplotype-level functional context from prior literature, "
        "but does not independently validate function, infer selection, make population-differentiation claims, "
        "or authorise biological interpretation from this project data."
    ),
}

OUTJSON.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
print(json.dumps(payload, indent=2))

if payload["status"] != "PASS":
    raise SystemExit(1)
