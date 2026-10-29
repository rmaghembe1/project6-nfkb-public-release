#!/usr/bin/env python3
"""Prepare an external consequence-annotation package for Phase C candidates.

This checkpoint packages the 500 candidate VCF and supporting metadata for
external or future local annotation. It does not run annotation and does not
assign consequences or biological interpretation.
"""

from __future__ import annotations

import csv
import json
import shutil
import zipfile
from pathlib import Path

CANDIDATE_UNION = Path(
    "11_results/phaseC/genotype_only_qc/functional_annotation_readiness/"
    "phaseC_frequency_variation_candidate_union_for_annotation_v1.tsv"
)
VCF_INPUT = Path(
    "11_results/phaseC/genotype_only_qc/consequence_annotation_context/"
    "phaseC_frequency_variation_candidate_union_for_annotation_v1.vcf"
)
BED_INPUT = Path(
    "11_results/phaseC/genotype_only_qc/functional_annotation_readiness/"
    "phaseC_frequency_variation_candidate_union_for_annotation_v1.bed"
)
ANNOTATION_CONTEXT_JSON = Path(
    "results/qc/phaseC_candidate_consequence_annotation_context_audit_v1.json"
)
ANNOTATION_CONTEXT_SUMMARY = Path(
    "11_results/phaseC/genotype_only_qc/consequence_annotation_context/"
    "phaseC_candidate_consequence_annotation_context_summary_v1.tsv"
)
TOOL_AVAILABILITY = Path(
    "11_results/phaseC/genotype_only_qc/consequence_annotation_context/"
    "phaseC_candidate_annotation_tool_availability_v1.tsv"
)
VARIANT_CLASS_COUNTS = Path(
    "11_results/phaseC/genotype_only_qc/consequence_annotation_context/"
    "phaseC_candidate_variant_class_counts_v1.tsv"
)
GENE_COUNTS = Path(
    "11_results/phaseC/genotype_only_qc/consequence_annotation_context/"
    "phaseC_candidate_gene_counts_v1.tsv"
)

OUTDIR = Path("11_results/phaseC/genotype_only_qc/external_annotation_package")
OUTJSON = Path("results/qc/phaseC_external_annotation_package_v1.json")

PACKAGE_DIR = OUTDIR / "phaseC_external_annotation_package_v1"
ZIP_PATH = OUTDIR / "phaseC_external_annotation_package_v1.zip"

EXPECTED_CANDIDATES = 500

OUTDIR.mkdir(parents=True, exist_ok=True)
OUTJSON.parent.mkdir(parents=True, exist_ok=True)

if PACKAGE_DIR.exists():
    shutil.rmtree(PACKAGE_DIR)
PACKAGE_DIR.mkdir(parents=True, exist_ok=True)

failures: list[str] = []


def count_vcf_records(path: Path) -> int:
    return sum(1 for line in path.read_text(encoding="utf-8").splitlines() if line and not line.startswith("#"))


def count_tsv_records(path: Path) -> int:
    with path.open(newline="", encoding="utf-8") as handle:
        return sum(1 for _ in csv.DictReader(handle, delimiter="\t"))


def copy_into_package(src: Path, dst_name: str) -> Path:
    dst = PACKAGE_DIR / dst_name
    shutil.copyfile(src, dst)
    return dst


vcf_records = count_vcf_records(VCF_INPUT)
union_records = count_tsv_records(CANDIDATE_UNION)

if vcf_records != EXPECTED_CANDIDATES:
    failures.append(f"Expected {EXPECTED_CANDIDATES} VCF records; observed {vcf_records}.")
if union_records != EXPECTED_CANDIDATES:
    failures.append(f"Expected {EXPECTED_CANDIDATES} union records; observed {union_records}.")

copied_files = [
    copy_into_package(VCF_INPUT, "candidate_variants_for_annotation.vcf"),
    copy_into_package(CANDIDATE_UNION, "candidate_variant_metadata.tsv"),
    copy_into_package(BED_INPUT, "candidate_variants_for_interval_tools.bed"),
    copy_into_package(ANNOTATION_CONTEXT_JSON, "annotation_context_audit.json"),
    copy_into_package(ANNOTATION_CONTEXT_SUMMARY, "annotation_context_summary.tsv"),
    copy_into_package(TOOL_AVAILABILITY, "local_tool_availability.tsv"),
    copy_into_package(VARIANT_CLASS_COUNTS, "candidate_variant_class_counts.tsv"),
    copy_into_package(GENE_COUNTS, "candidate_gene_counts.tsv"),
]

readme = PACKAGE_DIR / "README_external_annotation_package.md"
readme.write_text(
    """# Phase C external consequence-annotation package

## Purpose

This package contains 500 strict frequency-variation candidate variants prepared for downstream consequence annotation.

## Contents

- `candidate_variants_for_annotation.vcf`: true VCF input with 500 candidate records.
- `candidate_variant_metadata.tsv`: union candidate table with gene, allele key, parsed coordinates and descriptive source labels.
- `candidate_variants_for_interval_tools.bed`: BED-like coordinate input.
- `annotation_context_audit.json`: QC audit for candidate representation and local tool availability.
- `annotation_context_summary.tsv`: compact audit summary.
- `local_tool_availability.tsv`: local annotation tool availability at package creation.
- `candidate_variant_class_counts.tsv`: SNV/indel class counts.
- `candidate_gene_counts.tsv`: candidate counts by gene.
- `package_manifest.tsv`: package file manifest.

## Recommended annotation route

Use Ensembl VEP against GRCh38 with equivalent settings documented by the analyst. The recommended minimal outputs are:

- consequence terms;
- gene/transcript ID;
- canonical or MANE transcript flag where available;
- variant class;
- coding/protein consequence if applicable;
- impact/severity labels if produced by the annotation tool;
- nearest gene or regulatory consequence only if explicitly configured.

## Required import rule

Returned annotation outputs must be imported as a separate audited checkpoint. Do not overwrite this package.

## Boundary

The candidate labels in this package are descriptive QC/prioritisation-support labels only. This package does not assign functional consequence, infer selection, make population-differentiation claims, or authorise biological interpretation.
""",
    encoding="utf-8",
)

copied_files.append(readme)

purpose_map = {
    "candidate_variants_for_annotation.vcf": "Primary VCF input for downstream consequence annotation.",
    "candidate_variant_metadata.tsv": "Candidate metadata and descriptive source labels.",
    "candidate_variants_for_interval_tools.bed": "BED-like coordinate input for interval tools.",
    "annotation_context_audit.json": "Machine-readable audit context.",
    "annotation_context_summary.tsv": "Compact audit summary.",
    "local_tool_availability.tsv": "Local tool availability at package creation.",
    "candidate_variant_class_counts.tsv": "Variant class composition.",
    "candidate_gene_counts.tsv": "Candidate gene composition.",
    "README_external_annotation_package.md": "Human-readable package instructions and boundaries.",
    "package_manifest.tsv": "Package file manifest.",
}

manifest_path = PACKAGE_DIR / "package_manifest.tsv"

# First write a placeholder-free manifest including all files plus the manifest itself.
all_manifest_files = sorted(copied_files + [manifest_path], key=lambda p: p.name)

manifest_rows = []
for path in all_manifest_files:
    size_bytes = path.stat().st_size if path.exists() else 0
    manifest_rows.append({
        "package_file": path.name,
        "relative_path": path.name,
        "size_bytes": size_bytes,
        "purpose": purpose_map.get(path.name, "Package file."),
    })

with manifest_path.open("w", newline="", encoding="utf-8") as handle:
    writer = csv.DictWriter(handle, fieldnames=list(manifest_rows[0]), delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(manifest_rows)

# Recompute the manifest row size now that the manifest exists.
manifest_rows = []
all_package_files = sorted([p for p in PACKAGE_DIR.iterdir() if p.is_file()], key=lambda p: p.name)
for path in all_package_files:
    manifest_rows.append({
        "package_file": path.name,
        "relative_path": path.name,
        "size_bytes": path.stat().st_size,
        "purpose": purpose_map.get(path.name, "Package file."),
    })

with manifest_path.open("w", newline="", encoding="utf-8") as handle:
    writer = csv.DictWriter(handle, fieldnames=list(manifest_rows[0]), delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(manifest_rows)

all_package_files = sorted([p for p in PACKAGE_DIR.iterdir() if p.is_file()], key=lambda p: p.name)

with zipfile.ZipFile(ZIP_PATH, "w", compression=zipfile.ZIP_DEFLATED) as zf:
    for path in all_package_files:
        zf.write(path, arcname=path.name)

payload = {
    "status": "PASS" if not failures else "FAIL",
    "operation": "Prepare external consequence-annotation package for Phase C candidates.",
    "candidate_records": EXPECTED_CANDIDATES,
    "observed_vcf_records": vcf_records,
    "observed_union_records": union_records,
    "package_directory": str(PACKAGE_DIR),
    "package_zip": str(ZIP_PATH),
    "package_files": [p.name for p in all_package_files],
    "failures": failures,
    "boundary": (
        "This packages candidates for external/future annotation only. It does not run annotation, "
        "assign consequence, prioritise function, infer selection or authorise biological interpretation."
    ),
}

OUTJSON.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
print(json.dumps(payload, indent=2))

if payload["status"] != "PASS":
    raise SystemExit(1)
