#!/usr/bin/env python3
"""Build strict Phase C genotype-only analysis matrices from frozen technical-rule inventories."""

from __future__ import annotations

import csv
import json
import re
import subprocess
from collections import Counter
from pathlib import Path

STRICT = Path(
    "11_results/phaseC/genotype_only_qc/final_technical_rule/"
    "phaseC_strict_primary_genotype_only_records_v1.tsv"
)
MANIFEST = Path(
    "config/phaseB/full_panel_extraction/"
    "nfkb_full_panel_gene_level_extraction_manifest_GRCh38.tsv"
)
AUTOSOMAL_COHORT = Path(
    "11_results/phaseB/harmonisation_qc/"
    "cohort_primary_genotype_metadata_linked_n4150_v1.txt"
)
XX_SAMPLES = Path(
    "11_results/phaseC/genotype_only_qc/"
    "phaseC_IKBKG_candidate_primary_standard_XX_unflagged_samples_v1.txt"
)
XY_SAMPLES = Path(
    "11_results/phaseC/genotype_only_qc/"
    "phaseC_IKBKG_candidate_primary_standard_XY_unflagged_samples_v1.txt"
)

OUTDIR = Path("11_results/phaseC/genotype_only_qc/strict_analysis_matrices")
OUTJSON = Path("results/qc/phaseC_strict_genotype_only_analysis_matrices_v1.json")

EXPECTED_AUTOSOMAL_RECORDS = 39438
EXPECTED_CHRX_RECORDS = 383
EXPECTED_AUTOSOMAL_SAMPLES = 4150
EXPECTED_XX_SAMPLES = 1904
EXPECTED_XY_SAMPLES = 2216

OUTDIR.mkdir(parents=True, exist_ok=True)
OUTJSON.parent.mkdir(parents=True, exist_ok=True)


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def sample_count(path: Path) -> int:
    return sum(1 for line in path.read_text(encoding="utf-8").splitlines() if line)


def query_records(vcf: str, sample_file: Path, key_set: set[str]) -> list[list[str]]:
    cmd = [
        "bcftools",
        "query",
        "-S",
        str(sample_file),
        "-f",
        "%CHROM:%POS:%REF:%ALT[\t%GT]\n",
        vcf,
    ]
    run = subprocess.run(cmd, check=True, text=True, capture_output=True)
    rows: list[list[str]] = []
    for line in run.stdout.splitlines():
        fields = line.split("\t")
        key = fields[0]
        if key in key_set:
            rows.append(fields)
    return rows


strict = read_tsv(STRICT)
manifest = read_tsv(MANIFEST)
vcf_by_gene = {row["gene_symbol"]: row["genotype_output"] for row in manifest}

autosomal_strict = [row for row in strict if row["analysis_layer"] == "autosomal"]
chrx_strict = [row for row in strict if row["analysis_layer"] == "chrX_sex_aware"]

failures: list[str] = []

if len(autosomal_strict) != EXPECTED_AUTOSOMAL_RECORDS:
    failures.append(f"Expected {EXPECTED_AUTOSOMAL_RECORDS} autosomal strict records; observed {len(autosomal_strict)}.")
if len(chrx_strict) != EXPECTED_CHRX_RECORDS:
    failures.append(f"Expected {EXPECTED_CHRX_RECORDS} chrX strict records; observed {len(chrx_strict)}.")

if sample_count(AUTOSOMAL_COHORT) != EXPECTED_AUTOSOMAL_SAMPLES:
    failures.append("Autosomal cohort does not contain 4,150 samples.")
if sample_count(XX_SAMPLES) != EXPECTED_XX_SAMPLES:
    failures.append("XX chrX primary sample list does not contain 1,904 samples.")
if sample_count(XY_SAMPLES) != EXPECTED_XY_SAMPLES:
    failures.append("XY chrX primary sample list does not contain 2,216 samples.")

autosomal_keys_by_gene: dict[str, set[str]] = {}
for row in autosomal_strict:
    autosomal_keys_by_gene.setdefault(row["gene_symbol"], set()).add(row["allele_key"])

chrx_keys = {row["allele_key"] for row in chrx_strict}

autosomal_matrix = OUTDIR / "phaseC_strict_autosomal_genotype_GT_matrix_v1.tsv"
chrx_xx_matrix = OUTDIR / "phaseC_strict_IKBKG_chrX_XX_genotype_GT_matrix_v1.tsv"
chrx_xy_matrix = OUTDIR / "phaseC_strict_IKBKG_chrX_XY_genotype_GT_matrix_v1.tsv"

autosomal_record_count = 0
autosomal_unexpected = 0
autosomal_missing = 0
autosomal_called = 0
autosomal_by_gene = Counter()

with autosomal_matrix.open("w", newline="", encoding="utf-8") as handle:
    writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
    writer.writerow(["gene_symbol", "allele_key", "GT_values_in_frozen_4150_sample_order"])

    for gene in sorted(autosomal_keys_by_gene):
        rows = query_records(vcf_by_gene[gene], AUTOSOMAL_COHORT, autosomal_keys_by_gene[gene])
        observed_keys = {row[0] for row in rows}
        missing_keys = autosomal_keys_by_gene[gene] - observed_keys
        if missing_keys:
            failures.append(f"{gene}: {len(missing_keys)} strict autosomal keys missing from VCF query.")

        for fields in rows:
            key = fields[0]
            gts = fields[1:]
            if len(gts) != EXPECTED_AUTOSOMAL_SAMPLES:
                failures.append(f"{gene} {key}: expected 4,150 GT cells; observed {len(gts)}.")
                continue

            for gt in gts:
                if "." in gt:
                    autosomal_missing += 1
                elif gt in {"0/0", "0/1", "1/0", "1/1"}:
                    autosomal_called += 1
                else:
                    autosomal_unexpected += 1

            writer.writerow([gene, key, ",".join(gts)])
            autosomal_record_count += 1
            autosomal_by_gene[gene] += 1

def write_chrx_matrix(sample_file: Path, output: Path, stratum: str, expected_samples: int) -> tuple[int, int, int, int]:
    record_count = 0
    called = 0
    missing = 0
    unexpected = 0
    rows = query_records(vcf_by_gene["IKBKG"], sample_file, chrx_keys)
    observed_keys = {row[0] for row in rows}
    missing_keys = chrx_keys - observed_keys
    if missing_keys:
        failures.append(f"IKBKG {stratum}: {len(missing_keys)} strict chrX keys missing from VCF query.")

    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["gene_symbol", "allele_key", f"GT_values_in_frozen_{stratum}_sample_order"])

        for fields in rows:
            key = fields[0]
            gts = fields[1:]
            if len(gts) != expected_samples:
                failures.append(f"IKBKG {stratum} {key}: expected {expected_samples} GT cells; observed {len(gts)}.")
                continue

            for gt in gts:
                if "." in gt:
                    missing += 1
                elif stratum == "XX" and re.fullmatch(r"\d+/\d+", gt):
                    called += 1
                elif stratum == "XY" and re.fullmatch(r"\d+", gt):
                    called += 1
                else:
                    unexpected += 1

            writer.writerow(["IKBKG", key, ",".join(gts)])
            record_count += 1

    return record_count, called, missing, unexpected


xx_records, xx_called, xx_missing, xx_unexpected = write_chrx_matrix(
    XX_SAMPLES, chrx_xx_matrix, "XX", EXPECTED_XX_SAMPLES
)
xy_records, xy_called, xy_missing, xy_unexpected = write_chrx_matrix(
    XY_SAMPLES, chrx_xy_matrix, "XY", EXPECTED_XY_SAMPLES
)

summary_rows = [
    {
        "matrix": "strict_autosomal",
        "records": autosomal_record_count,
        "samples": EXPECTED_AUTOSOMAL_SAMPLES,
        "GT_cells": autosomal_record_count * EXPECTED_AUTOSOMAL_SAMPLES,
        "called_GT_cells": autosomal_called,
        "missing_GT_cells": autosomal_missing,
        "unexpected_nonmissing_GT_cells": autosomal_unexpected,
    },
    {
        "matrix": "strict_IKBKG_chrX_XX",
        "records": xx_records,
        "samples": EXPECTED_XX_SAMPLES,
        "GT_cells": xx_records * EXPECTED_XX_SAMPLES,
        "called_GT_cells": xx_called,
        "missing_GT_cells": xx_missing,
        "unexpected_nonmissing_GT_cells": xx_unexpected,
    },
    {
        "matrix": "strict_IKBKG_chrX_XY",
        "records": xy_records,
        "samples": EXPECTED_XY_SAMPLES,
        "GT_cells": xy_records * EXPECTED_XY_SAMPLES,
        "called_GT_cells": xy_called,
        "missing_GT_cells": xy_missing,
        "unexpected_nonmissing_GT_cells": xy_unexpected,
    },
]

summary_path = OUTDIR / "phaseC_strict_genotype_only_GT_matrix_summary_v1.tsv"
with summary_path.open("w", newline="", encoding="utf-8") as handle:
    writer = csv.DictWriter(handle, fieldnames=list(summary_rows[0]), delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(summary_rows)

gene_counts_path = OUTDIR / "phaseC_strict_autosomal_records_by_gene_v1.tsv"
with gene_counts_path.open("w", newline="", encoding="utf-8") as handle:
    writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
    writer.writerow(["gene_symbol", "strict_autosomal_record_count"])
    for gene, count in sorted(autosomal_by_gene.items()):
        writer.writerow([gene, count])

if autosomal_record_count != EXPECTED_AUTOSOMAL_RECORDS:
    failures.append(f"Autosomal matrix contains {autosomal_record_count}; expected {EXPECTED_AUTOSOMAL_RECORDS}.")
if xx_records != EXPECTED_CHRX_RECORDS or xy_records != EXPECTED_CHRX_RECORDS:
    failures.append(f"chrX matrices contain XX={xx_records}, XY={xy_records}; expected {EXPECTED_CHRX_RECORDS} each.")
if autosomal_unexpected != 0 or xx_unexpected != 0 or xy_unexpected != 0:
    failures.append("Unexpected nonmissing GT encodings were observed in strict matrices.")

payload = {
    "status": "PASS" if not failures else "FAIL",
    "operation": "Build strict genotype-only GT matrices from frozen Phase C technical rule.",
    "strict_record_counts": {
        "autosomal": autosomal_record_count,
        "IKBKG_chrX": EXPECTED_CHRX_RECORDS,
    },
    "sample_counts": {
        "autosomal": EXPECTED_AUTOSOMAL_SAMPLES,
        "IKBKG_chrX_XX": EXPECTED_XX_SAMPLES,
        "IKBKG_chrX_XY": EXPECTED_XY_SAMPLES,
    },
    "matrix_outputs": {
        "autosomal": str(autosomal_matrix),
        "IKBKG_chrX_XX": str(chrx_xx_matrix),
        "IKBKG_chrX_XY": str(chrx_xy_matrix),
    },
    "GT_cell_summary": summary_rows,
    "failures": failures,
    "boundary": (
        "Matrices contain GT strings only. No allele frequencies, population comparisons, "
        "selection inference, functional prioritisation or biological interpretation are authorised."
    ),
}

OUTJSON.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
print(json.dumps(payload, indent=2))

if payload["status"] != "PASS":
    raise SystemExit(1)
