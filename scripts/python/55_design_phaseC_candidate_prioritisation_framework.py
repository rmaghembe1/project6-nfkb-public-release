#!/usr/bin/env python3
"""Design Phase C candidate-prioritisation framework.

This checkpoint defines transparent prioritisation rules before applying them.
It does not rank variants, select candidates, interpret biology, infer selection,
or make population-differentiation claims.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path

LANDSCAPE_JSON = Path("results/qc/phaseC_vep_consequence_landscape_profile_v1.json")
FREQ_PROFILE_JSON = Path("results/qc/phaseC_strict_frequency_variation_candidate_profile_v1.json")
POP_STRATA_JSON = Path("results/qc/phaseC_genotype_only_population_strata_validation_v1.json")
GROUP_FREQ_JSON = Path("results/qc/phaseC_strict_genotype_only_group_frequency_qc_v1.json")
TECH_RULE_JSON = Path("results/qc/phaseC_genotype_only_technical_rule_freeze_v1.json")

OUTDIR = Path("11_results/phaseC/genotype_only_qc/candidate_prioritisation_framework")
OUTJSON = Path("results/qc/phaseC_candidate_prioritisation_framework_design_v1.json")

OUTDIR.mkdir(parents=True, exist_ok=True)
OUTJSON.parent.mkdir(parents=True, exist_ok=True)

failures: list[str] = []

required_inputs = [
    LANDSCAPE_JSON,
    FREQ_PROFILE_JSON,
    POP_STRATA_JSON,
    GROUP_FREQ_JSON,
    TECH_RULE_JSON,
]

for path in required_inputs:
    if not path.exists():
        failures.append(f"Required upstream QC JSON missing: {path}")

def write_rows(path: Path, fieldnames: list[str], rows: list[dict[str, object]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row.get(key, ".") for key in fieldnames})

framework_rows = [
    {
        "domain": "technical_inclusion",
        "criterion": "candidate_must_derive_from_frozen_strict_genotype_only_inventory",
        "purpose": "Prevent prioritisation of records outside the approved Phase C technical rule.",
        "allowed_values_or_rule": "YES required",
        "score_direction": "gate_only",
        "weight": 0,
        "boundary": "Technical gate only; not evidence of biology.",
    },
    {
        "domain": "technical_inclusion",
        "criterion": "candidate_must_have_valid_imported_VEP_annotation",
        "purpose": "Ensure candidate has audited consequence annotation.",
        "allowed_values_or_rule": "YES required",
        "score_direction": "gate_only",
        "weight": 0,
        "boundary": "Annotation availability only; not functional proof.",
    },
    {
        "domain": "frequency_variation",
        "criterion": "broad_group_AF_delta_rank",
        "purpose": "Capture descriptive broad-group allele-frequency variation among strict candidates.",
        "allowed_values_or_rule": "Use precomputed strict broad AF-delta candidate table; higher delta receives higher technical priority.",
        "score_direction": "higher_is_more_prioritised",
        "weight": 3,
        "boundary": "Descriptive frequency contrast only; not selection or ancestry effect.",
    },
    {
        "domain": "frequency_variation",
        "criterion": "observed_only_or_group_restricted_flag",
        "purpose": "Flag variants observed in limited broad groups for technical follow-up.",
        "allowed_values_or_rule": "Flag only; must be interpreted with sample size and missingness context.",
        "score_direction": "flag_not_standalone_score",
        "weight": 1,
        "boundary": "Not evidence of population specificity without additional modelling.",
    },
    {
        "domain": "consequence",
        "criterion": "canonical_like_VEP_impact",
        "purpose": "Use MANE/APPRIS canonical-like consequence impact as a technical annotation severity layer.",
        "allowed_values_or_rule": "MODERATE > LOW > MODIFIER; HIGH would be highest if present.",
        "score_direction": "higher_is_more_prioritised",
        "weight": 3,
        "boundary": "VEP consequence impact is annotation severity, not functional validation.",
    },
    {
        "domain": "consequence",
        "criterion": "any_transcript_moderate_or_high_flag",
        "purpose": "Retain variants with non-canonical moderate/high signal for review without over-weighting.",
        "allowed_values_or_rule": "YES/NO flag",
        "score_direction": "yes_adds_small_review_weight",
        "weight": 1,
        "boundary": "All-transcript signal may reflect noncanonical isoforms; requires transcript-policy review.",
    },
    {
        "domain": "gene_context",
        "criterion": "gene_module_role",
        "purpose": "Record pathway-module context within the predefined NF-kB regulatory-circuit panel.",
        "allowed_values_or_rule": "Core TF, IKK/IκB, noncanonical branch, ubiquitin editing/restraint, or other predefined panel module.",
        "score_direction": "descriptive_context",
        "weight": 1,
        "boundary": "Module membership supports contextual review only, not causality.",
    },
    {
        "domain": "robustness",
        "criterion": "small_n_or_unstable_group_flag",
        "purpose": "Prevent over-prioritisation of candidates driven by insufficient group counts.",
        "allowed_values_or_rule": "Candidates depending on group n < 20 or technical-review records must be downgraded or held.",
        "score_direction": "downgrade_or_hold",
        "weight": -3,
        "boundary": "Technical robustness gate.",
    },
    {
        "domain": "interpretation_boundary",
        "criterion": "claim_level",
        "purpose": "Force every candidate output to carry an interpretation boundary.",
        "allowed_values_or_rule": "technical_candidate_only, annotation_supported_candidate, follow_up_candidate; no causal/selection claim.",
        "score_direction": "mandatory_label",
        "weight": 0,
        "boundary": "No biological claim authorised by prioritisation alone.",
    },
]

framework_path = OUTDIR / "phaseC_candidate_prioritisation_framework_rules_v1.tsv"
write_rows(
    framework_path,
    ["domain", "criterion", "purpose", "allowed_values_or_rule", "score_direction", "weight", "boundary"],
    framework_rows,
)

score_rows = [
    {
        "score_component": "frequency_delta_score",
        "source_table": "phaseC_strict_top_broad_AF_delta_candidates_v1.tsv",
        "planned_rule": "Assign ordinal bins from broad AF delta distribution after confirming denominators and missingness.",
        "provisional_range": "0-3",
        "interpretation_limit": "Descriptive frequency variation only.",
    },
    {
        "score_component": "canonical_consequence_score",
        "source_table": "phaseC_vep_landscape_per_candidate_technical_summary_v1.tsv",
        "planned_rule": "MODIFIER=0, LOW=1, MODERATE=3, HIGH=4 if ever present; canonical-like only.",
        "provisional_range": "0-4",
        "interpretation_limit": "Annotation severity only.",
    },
    {
        "score_component": "any_transcript_review_score",
        "source_table": "phaseC_vep_landscape_per_candidate_technical_summary_v1.tsv",
        "planned_rule": "Add 1 if any transcript has MODERATE/HIGH but canonical-like does not.",
        "provisional_range": "0-1",
        "interpretation_limit": "Review flag only.",
    },
    {
        "score_component": "gene_context_score",
        "source_table": "predefined Project 6 NF-kB module map to be generated in the next checkpoint",
        "planned_rule": "Add contextual module label and modest score only after module-map file is frozen.",
        "provisional_range": "0-1",
        "interpretation_limit": "Pathway context only.",
    },
    {
        "score_component": "robustness_penalty",
        "source_table": "population strata, group frequency and technical-rule QC outputs",
        "planned_rule": "Downgrade or hold candidates dependent on small-n, missingness or technical-review conditions.",
        "provisional_range": "0 to -3",
        "interpretation_limit": "Technical caution only.",
    },
]

score_path = OUTDIR / "phaseC_candidate_prioritisation_score_components_v1.tsv"
write_rows(
    score_path,
    ["score_component", "source_table", "planned_rule", "provisional_range", "interpretation_limit"],
    score_rows,
)

claim_rows = [
    {
        "label": "technical_candidate_only",
        "definition": "Candidate passed technical discovery and annotation workflow but lacks functional or inferential support.",
        "allowed_language": "technical candidate; descriptive candidate; annotation-ready candidate",
        "prohibited_language": "functional variant; adaptive variant; ancestry-associated causal variant",
    },
    {
        "label": "annotation_supported_candidate",
        "definition": "Candidate has moderate/high canonical-like VEP consequence or clear transcript consequence flag.",
        "allowed_language": "annotation-supported follow-up candidate",
        "prohibited_language": "damaging variant unless independently validated by functional predictors and evidence review",
    },
    {
        "label": "frequency_variation_candidate",
        "definition": "Candidate shows descriptive broad-group AF variation under strict QC.",
        "allowed_language": "frequency-variation candidate under strict QC",
        "prohibited_language": "population-specific, selected, adaptive, or ancestry-determining variant",
    },
    {
        "label": "follow_up_candidate",
        "definition": "Candidate may be advanced to later functional/regulatory/literature review.",
        "allowed_language": "candidate for follow-up annotation and evidence review",
        "prohibited_language": "validated mechanism or biological effect",
    },
]

claim_path = OUTDIR / "phaseC_candidate_prioritisation_claim_boundaries_v1.tsv"
write_rows(
    claim_path,
    ["label", "definition", "allowed_language", "prohibited_language"],
    claim_rows,
)

next_step_rows = [
    {
        "order": 1,
        "next_checkpoint": "freeze_gene_module_context_map",
        "reason": "Gene/module context must be frozen before scoring to avoid post-hoc interpretation.",
    },
    {
        "order": 2,
        "next_checkpoint": "apply_prioritisation_framework",
        "reason": "Apply the predeclared framework to candidate tables after module context is frozen.",
    },
    {
        "order": 3,
        "next_checkpoint": "audit_top_candidate_table",
        "reason": "Check that top candidates are not driven by small-n, missingness or metadata artefacts.",
    },
    {
        "order": 4,
        "next_checkpoint": "evidence_review_gate",
        "reason": "Only after technical ranking should literature/functional evidence be reviewed.",
    },
]

next_path = OUTDIR / "phaseC_candidate_prioritisation_next_checkpoint_plan_v1.tsv"
write_rows(
    next_path,
    ["order", "next_checkpoint", "reason"],
    next_step_rows,
)

payload = {
    "status": "PASS" if not failures else "FAIL",
    "operation": "Design Phase C candidate-prioritisation framework before applying ranking.",
    "framework_rules": len(framework_rows),
    "score_components": len(score_rows),
    "claim_boundary_labels": len(claim_rows),
    "next_checkpoints": len(next_step_rows),
    "required_upstream_inputs_checked": len(required_inputs),
    "outputs": {
        "framework_rules": str(framework_path),
        "score_components": str(score_path),
        "claim_boundaries": str(claim_path),
        "next_checkpoint_plan": str(next_path),
    },
    "failures": failures,
    "boundary": (
        "This checkpoint designs prioritisation rules only. It does not apply scores, rank variants, "
        "select candidates, infer function, infer selection, make population-differentiation claims, "
        "or authorise biological interpretation."
    ),
}

OUTJSON.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
print(json.dumps(payload, indent=2))

if payload["status"] != "PASS":
    raise SystemExit(1)
