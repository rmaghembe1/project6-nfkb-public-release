#!/usr/bin/env python3
"""Prepare Phase C evidence-review gate.

This checkpoint defines rules for later evidence review of technically prioritised
Phase C candidates. It does not perform literature review, functional annotation
interpretation, selection inference, or biological claim-making.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path

TOP_AUDIT_JSON = Path("results/qc/phaseC_top_candidate_table_audit_v1.json")
FOLLOW_UP = Path(
    "11_results/phaseC/genotype_only_qc/candidate_prioritisation_applied/"
    "phaseC_candidate_prioritisation_follow_up_candidates_v1.tsv"
)
ANNOTATION_SUPPORTED = Path(
    "11_results/phaseC/genotype_only_qc/top_candidate_audit/"
    "phaseC_annotation_supported_candidate_audit_table_v1.tsv"
)

OUTDIR = Path("11_results/phaseC/genotype_only_qc/evidence_review_gate")
OUTJSON = Path("results/qc/phaseC_evidence_review_gate_v1.json")

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


for path in [TOP_AUDIT_JSON, FOLLOW_UP, ANNOTATION_SUPPORTED]:
    if not path.exists():
        failures.append(f"Required input missing: {path}")

follow_up_rows = read_tsv(FOLLOW_UP) if FOLLOW_UP.exists() else []
annotation_rows = read_tsv(ANNOTATION_SUPPORTED) if ANNOTATION_SUPPORTED.exists() else []

if len(follow_up_rows) != 173:
    failures.append(f"Expected 173 follow-up candidate rows; observed {len(follow_up_rows)}.")
if len(annotation_rows) != 4:
    failures.append(f"Expected 4 annotation-supported candidate rows; observed {len(annotation_rows)}.")

gate_rows = [
    {
        "gate_domain": "entry_scope",
        "rule_id": "ERG_001",
        "rule": "Only non-technical-only candidates enter the evidence-review gate.",
        "applies_to": "173 follow-up candidate rows",
        "required_action": "Review evidence only after preserving technical claim-boundary labels.",
        "boundary": "Entry into review is not evidence of function.",
    },
    {
        "gate_domain": "entry_scope",
        "rule_id": "ERG_002",
        "rule": "Annotation-supported candidates receive first-pass evidence review priority.",
        "applies_to": "4 annotation-supported candidates",
        "required_action": "Review variant-level and gene-level evidence separately.",
        "boundary": "VEP MODERATE consequence is not functional validation.",
    },
    {
        "gate_domain": "evidence_type",
        "rule_id": "ERG_003",
        "rule": "Variant-level evidence must be separated from gene-level pathway evidence.",
        "applies_to": "all reviewed candidates",
        "required_action": "Record whether evidence is exact variant, same residue/domain, same gene, or pathway-only.",
        "boundary": "Gene/pathway evidence cannot be treated as variant-specific evidence.",
    },
    {
        "gate_domain": "evidence_type",
        "rule_id": "ERG_004",
        "rule": "Allowed evidence types include curated variant databases, primary functional studies, protein/domain annotations, eQTL/sQTL/regulatory evidence and disease association records.",
        "applies_to": "all reviewed candidates",
        "required_action": "Classify every evidence item by source type and evidence level.",
        "boundary": "Database presence alone does not establish functional effect.",
    },
    {
        "gate_domain": "population_interpretation",
        "rule_id": "ERG_005",
        "rule": "Frequency variation must remain descriptive unless formal population-genetic tests are performed later.",
        "applies_to": "frequency-variation candidates",
        "required_action": "Retain wording as descriptive broad-group AF variation under strict QC.",
        "boundary": "No selection, adaptation, ancestry-effect or population-specific claim is allowed.",
    },
    {
        "gate_domain": "claim_control",
        "rule_id": "ERG_006",
        "rule": "Every reviewed candidate must receive an evidence-grade label.",
        "applies_to": "all reviewed candidates",
        "required_action": "Use: no_external_support, gene_level_support_only, indirect_variant_context, variant_level_support, functional_validation_available.",
        "boundary": "Only functional_validation_available can support strong functional language.",
    },
    {
        "gate_domain": "claim_control",
        "rule_id": "ERG_007",
        "rule": "Manuscript-facing language must not exceed the weakest evidence boundary.",
        "applies_to": "all reviewed candidates",
        "required_action": "Use cautious terms such as technical candidate, evidence-review candidate, or follow-up candidate.",
        "boundary": "Avoid causal, adaptive, selected, disease-causing or functionally validated language unless independently supported.",
    },
]

gate_path = OUTDIR / "phaseC_evidence_review_gate_rules_v1.tsv"
write_rows(
    gate_path,
    ["gate_domain", "rule_id", "rule", "applies_to", "required_action", "boundary"],
    gate_rows,
)

evidence_grade_rows = [
    {
        "evidence_grade": "no_external_support",
        "definition": "No relevant external evidence found beyond this project technical workflow.",
        "allowed_language": "technical follow-up candidate",
        "prohibited_language": "functional or disease-relevant variant",
    },
    {
        "evidence_grade": "gene_level_support_only",
        "definition": "External evidence supports the gene/pathway role but not the exact variant.",
        "allowed_language": "candidate in a biologically relevant gene context",
        "prohibited_language": "variant has functional effect",
    },
    {
        "evidence_grade": "indirect_variant_context",
        "definition": "Evidence involves nearby residue/domain/regulatory region or computational prediction but not direct functional validation.",
        "allowed_language": "variant with indirect contextual support",
        "prohibited_language": "validated functional variant",
    },
    {
        "evidence_grade": "variant_level_support",
        "definition": "Exact variant has curated, association, regulatory or prior reported evidence, but not necessarily experimental validation.",
        "allowed_language": "variant with external supporting evidence",
        "prohibited_language": "causal or adaptive variant unless evidence supports that exact claim",
    },
    {
        "evidence_grade": "functional_validation_available",
        "definition": "Exact variant has direct experimental evidence for molecular, cellular or organismal effect.",
        "allowed_language": "functionally supported variant, with precise context",
        "prohibited_language": "generalized causal or population-level claims without appropriate evidence",
    },
]

grade_path = OUTDIR / "phaseC_evidence_review_grade_labels_v1.tsv"
write_rows(
    grade_path,
    ["evidence_grade", "definition", "allowed_language", "prohibited_language"],
    evidence_grade_rows,
)

review_template_rows = []
for row in annotation_rows:
    review_template_rows.append(
        {
            "allele_key": row["allele_key"],
            "gene": row["gene"],
            "claim_boundary_label": row["claim_boundary_label"],
            "technical_priority_rank": row["technical_priority_rank"],
            "technical_priority_score": row["technical_priority_score"],
            "canonical_impact_label": row["canonical_impact_label"],
            "review_priority": "first_pass_annotation_supported",
            "variant_level_evidence_grade": "not_reviewed",
            "gene_level_evidence_grade": "not_reviewed",
            "regulatory_or_qtl_evidence_grade": "not_reviewed",
            "population_genetic_evidence_grade": "not_reviewed",
            "review_status": "pending_evidence_review",
            "interpretation_boundary": "Evidence gate prepared only; no external evidence reviewed yet.",
        }
    )

template_path = OUTDIR / "phaseC_annotation_supported_evidence_review_template_v1.tsv"
write_rows(
    template_path,
    [
        "allele_key",
        "gene",
        "claim_boundary_label",
        "technical_priority_rank",
        "technical_priority_score",
        "canonical_impact_label",
        "review_priority",
        "variant_level_evidence_grade",
        "gene_level_evidence_grade",
        "regulatory_or_qtl_evidence_grade",
        "population_genetic_evidence_grade",
        "review_status",
        "interpretation_boundary",
    ],
    review_template_rows,
)

source_plan_rows = [
    {
        "source_category": "curated_variant_database",
        "examples": "ClinVar, dbSNP annotations, gnomAD variant pages where appropriate",
        "use": "variant identity, prior reports, clinical assertions where present",
        "limitation": "Presence or assertion is not functional proof.",
    },
    {
        "source_category": "protein_domain_annotation",
        "examples": "UniProt, InterPro, Pfam, AlphaFold/PDB context where available",
        "use": "domain/residue context",
        "limitation": "Domain location is contextual, not proof of effect.",
    },
    {
        "source_category": "functional_literature",
        "examples": "primary studies on exact variant, same domain/residue, or gene mechanism",
        "use": "evidence grading",
        "limitation": "Gene-level studies must not be promoted to variant-level proof.",
    },
    {
        "source_category": "regulatory_QTL_evidence",
        "examples": "GTEx, immune stimulation QTL/chromatin response resources where available",
        "use": "regulatory plausibility and expression/splicing context",
        "limitation": "QTL evidence requires allele matching and context matching.",
    },
    {
        "source_category": "population_genetic_follow_up",
        "examples": "FST, PBS, haplotype/LD, archaic-introgression analyses if later performed",
        "use": "future formal evolutionary analysis",
        "limitation": "Not part of current evidence gate; no selection claim now.",
    },
]

source_path = OUTDIR / "phaseC_evidence_review_source_plan_v1.tsv"
write_rows(
    source_path,
    ["source_category", "examples", "use", "limitation"],
    source_plan_rows,
)

payload = {
    "status": "PASS" if not failures else "FAIL",
    "operation": "Prepare Phase C evidence-review gate for technically prioritised candidates.",
    "follow_up_candidate_rows": len(follow_up_rows),
    "annotation_supported_candidate_rows": len(annotation_rows),
    "gate_rules": len(gate_rows),
    "evidence_grade_labels": len(evidence_grade_rows),
    "annotation_supported_review_template_rows": len(review_template_rows),
    "source_plan_rows": len(source_plan_rows),
    "outputs": {
        "gate_rules": str(gate_path),
        "evidence_grade_labels": str(grade_path),
        "annotation_supported_review_template": str(template_path),
        "source_plan": str(source_path),
    },
    "failures": failures,
    "boundary": (
        "This checkpoint prepares evidence-review rules only. It does not perform evidence review, "
        "validate function, infer selection, make population-differentiation claims, or authorise biological interpretation."
    ),
}

OUTJSON.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
print(json.dumps(payload, indent=2))

if payload["status"] != "PASS":
    raise SystemExit(1)
