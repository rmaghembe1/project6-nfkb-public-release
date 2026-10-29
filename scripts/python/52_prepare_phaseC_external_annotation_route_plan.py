#!/usr/bin/env python3
"""Prepare external annotation route plan for Phase C candidates.

This checkpoint documents how the external annotation package should be used
and how returned annotation results must be imported. It does not submit files,
run annotation, import annotation, assign consequence, or interpret biology.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path

PACKAGE_JSON = Path("results/qc/phaseC_external_annotation_package_v1.json")
PACKAGE_ZIP = Path(
    "11_results/phaseC/genotype_only_qc/external_annotation_package/"
    "phaseC_external_annotation_package_v1.zip"
)
PACKAGE_DIR = Path(
    "11_results/phaseC/genotype_only_qc/external_annotation_package/"
    "phaseC_external_annotation_package_v1"
)
VCF = PACKAGE_DIR / "candidate_variants_for_annotation.vcf"
METADATA = PACKAGE_DIR / "candidate_variant_metadata.tsv"
MANIFEST = PACKAGE_DIR / "package_manifest.tsv"

OUTDIR = Path("11_results/phaseC/genotype_only_qc/external_annotation_route")
OUTJSON = Path("results/qc/phaseC_external_annotation_route_plan_v1.json")

EXPECTED_RECORDS = 500

OUTDIR.mkdir(parents=True, exist_ok=True)
OUTJSON.parent.mkdir(parents=True, exist_ok=True)

failures: list[str] = []


def count_vcf_records(path: Path) -> int:
    return sum(1 for line in path.read_text(encoding="utf-8").splitlines() if line and not line.startswith("#"))


def count_tsv_records(path: Path) -> int:
    with path.open(newline="", encoding="utf-8") as handle:
        return sum(1 for _ in csv.DictReader(handle, delimiter="\t"))


vcf_records = count_vcf_records(VCF)
metadata_records = count_tsv_records(METADATA)

if vcf_records != EXPECTED_RECORDS:
    failures.append(f"Expected {EXPECTED_RECORDS} VCF records; observed {vcf_records}.")
if metadata_records != EXPECTED_RECORDS:
    failures.append(f"Expected {EXPECTED_RECORDS} metadata records; observed {metadata_records}.")
if not PACKAGE_ZIP.exists():
    failures.append(f"Package ZIP missing: {PACKAGE_ZIP}.")
if not MANIFEST.exists():
    failures.append(f"Package manifest missing: {MANIFEST}.")

route_rows = [
    {
        "route": "Ensembl_VEP_web_or_API_GRCh38",
        "input_file": "candidate_variants_for_annotation.vcf",
        "genome_build": "GRCh38",
        "recommended_outputs": (
            "variant_id, uploaded_variation, location, allele, gene_id, gene_symbol, feature_id, "
            "feature_type, consequence_terms, impact, canonical_flag, MANE_Select_or_MANE_Plus_Clinical_if_available, "
            "HGVSc, HGVSp, cDNA_position, CDS_position, protein_position, existing_variation_if_available"
        ),
        "required_boundary": "Use for consequence annotation only; do not infer function or biological importance at this stage.",
    },
    {
        "route": "Local_VEP_future",
        "input_file": "candidate_variants_for_annotation.vcf",
        "genome_build": "GRCh38",
        "recommended_outputs": (
            "tabular VEP output plus run parameters, cache/database version, command line, and warning log"
        ),
        "required_boundary": "Import only through a separate audited checkpoint.",
    },
    {
        "route": "Collaborator_annotation_environment",
        "input_file": "phaseC_external_annotation_package_v1.zip",
        "genome_build": "GRCh38",
        "recommended_outputs": (
            "annotated TSV/VCF, tool name, tool version, genome build, database/cache version, command or settings, run date"
        ),
        "required_boundary": "Returned annotations must be audited before use.",
    },
]

route_plan_path = OUTDIR / "phaseC_external_annotation_route_plan_v1.tsv"
with route_plan_path.open("w", newline="", encoding="utf-8") as handle:
    writer = csv.DictWriter(handle, fieldnames=list(route_rows[0]), delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(route_rows)

expected_return_rows = [
    {"required_return_item": "annotated_candidate_file", "description": "Annotated TSV, VCF or VEP output covering the 500 submitted candidate records."},
    {"required_return_item": "tool_name", "description": "Name of annotation tool used, for example Ensembl VEP."},
    {"required_return_item": "tool_version", "description": "Tool/software version."},
    {"required_return_item": "genome_build", "description": "Must be GRCh38 or clearly documented equivalent."},
    {"required_return_item": "database_or_cache_version", "description": "Annotation database/cache/release version."},
    {"required_return_item": "run_settings_or_command", "description": "Exact settings, command line, or web/API options used."},
    {"required_return_item": "run_date", "description": "Date annotation was run."},
    {"required_return_item": "warning_or_error_log", "description": "Warnings, failed records, unrecognised alleles, or empty output statement."},
]

return_spec_path = OUTDIR / "phaseC_external_annotation_expected_return_spec_v1.tsv"
with return_spec_path.open("w", newline="", encoding="utf-8") as handle:
    writer = csv.DictWriter(handle, fieldnames=list(expected_return_rows[0]), delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(expected_return_rows)

import_checklist_rows = [
    {"audit_step": "record_count_check", "required_action": "Confirm returned annotation has exactly 500 submitted candidate IDs or explicitly account for missing/failed records."},
    {"audit_step": "allele_identity_check", "required_action": "Verify CHROM, POS, REF and ALT match the submitted VCF/candidate metadata."},
    {"audit_step": "build_check", "required_action": "Confirm annotation was run on GRCh38."},
    {"audit_step": "tool_version_check", "required_action": "Record tool name, version and annotation database/cache version."},
    {"audit_step": "consequence_schema_check", "required_action": "Map returned consequence fields to a stable internal schema before biological review."},
    {"audit_step": "multi_transcript_policy", "required_action": "Separate canonical/MANE consequence from all-transcript consequences; do not collapse without documented rule."},
    {"audit_step": "boundary_check", "required_action": "Do not interpret functional importance, selection or population differentiation during import."},
]

checklist_path = OUTDIR / "phaseC_returned_annotation_import_audit_checklist_v1.tsv"
with checklist_path.open("w", newline="", encoding="utf-8") as handle:
    writer = csv.DictWriter(handle, fieldnames=list(import_checklist_rows[0]), delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(import_checklist_rows)

readme_path = OUTDIR / "README_external_annotation_route_v1.md"
readme_path.write_text(
    """# Phase C external annotation route plan

## Package to use

Use the committed package:

`11_results/phaseC/genotype_only_qc/external_annotation_package/phaseC_external_annotation_package_v1.zip`

For annotation, the primary file inside the package is:

`candidate_variants_for_annotation.vcf`

## Recommended route

Run consequence annotation against GRCh38, preferably using Ensembl VEP or an equivalent controlled annotation environment.

## Required returned metadata

Returned annotation must include or be accompanied by:

- tool name;
- tool version;
- genome build;
- annotation database/cache/release version;
- command line, web settings or API settings;
- run date;
- warning/error log or explicit no-warning statement.

## Import rule

Returned annotation results must not overwrite this package. They must be imported in a separate audited checkpoint with record-count, allele-identity, genome-build and schema checks.

## Boundary

This route plan does not run annotation, assign consequence, prioritise variants, infer selection, make population-differentiation claims, or authorise biological interpretation.
""",
    encoding="utf-8",
)

payload = {
    "status": "PASS" if not failures else "FAIL",
    "operation": "Prepare external annotation route plan and returned-result import rules.",
    "package_zip": str(PACKAGE_ZIP),
    "submitted_vcf_records": vcf_records,
    "metadata_records": metadata_records,
    "outputs": {
        "route_plan": str(route_plan_path),
        "expected_return_spec": str(return_spec_path),
        "returned_annotation_import_checklist": str(checklist_path),
        "readme": str(readme_path),
    },
    "failures": failures,
    "boundary": (
        "This prepares the external annotation route only. It does not submit files, run annotation, "
        "import annotation, assign consequence, prioritise function, infer selection or authorise biological interpretation."
    ),
}

OUTJSON.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
print(json.dumps(payload, indent=2))

if payload["status"] != "PASS":
    raise SystemExit(1)
