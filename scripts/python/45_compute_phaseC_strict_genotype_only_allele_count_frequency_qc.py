#!/usr/bin/env python3
"""Compute strict genotype-only allele-count and frequency QC tables.

This is a technical QC step only. It computes AC/AN/AF from frozen strict GT
matrices but does not authorise population comparison or biological interpretation.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path

MATRIX_DIR = Path("11_results/phaseC/genotype_only_qc/strict_analysis_matrices")
AUTO_MATRIX = MATRIX_DIR / "phaseC_strict_autosomal_genotype_GT_matrix_v1.tsv"
CHRX_XX_MATRIX = MATRIX_DIR / "phaseC_strict_IKBKG_chrX_XX_genotype_GT_matrix_v1.tsv"
CHRX_XY_MATRIX = MATRIX_DIR / "phaseC_strict_IKBKG_chrX_XY_genotype_GT_matrix_v1.tsv"

OUTDIR = Path("11_results/phaseC/genotype_only_qc/strict_frequency_qc")
OUTJSON = Path("results/qc/phaseC_strict_genotype_only_allele_count_frequency_qc_v1.json")

EXPECTED_AUTO_RECORDS = 39438
EXPECTED_CHRX_RECORDS = 383
EXPECTED_AUTO_SAMPLES = 4150
EXPECTED_XX_SAMPLES = 1904
EXPECTED_XY_SAMPLES = 2216

OUTDIR.mkdir(parents=True, exist_ok=True)
OUTJSON.parent.mkdir(parents=True, exist_ok=True)


def parse_autosomal_gt(gt: str) -> tuple[int, int, int]:
    """Return called_alleles, alt_alleles, missing_cells for diploid autosomal GT."""
    if "." in gt:
        return 0, 0, 1
    if gt not in {"0/0", "0/1", "1/0", "1/1"}:
        raise ValueError(f"Unexpected autosomal GT: {gt}")
    a, b = gt.split("/")
    return 2, int(a) + int(b), 0


def parse_xx_gt(gt: str) -> tuple[int, int, int]:
    """Return called_alleles, alt_alleles, missing_cells for diploid XX chrX GT."""
    if "." in gt:
        return 0, 0, 1
    if "/" not in gt:
        raise ValueError(f"Unexpected XX chrX GT: {gt}")
    a, b = gt.split("/")
    return 2, int(a) + int(b), 0


def parse_xy_gt(gt: str) -> tuple[int, int, int]:
    """Return called_alleles, alt_alleles, missing_cells for haploid XY chrX GT."""
    if "." in gt:
        return 0, 0, 1
    if "/" in gt or "|" in gt:
        raise ValueError(f"Unexpected XY chrX GT: {gt}")
    return 1, int(gt), 0


def read_matrix(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


failures: list[str] = []

auto_rows = read_matrix(AUTO_MATRIX)
xx_rows = read_matrix(CHRX_XX_MATRIX)
xy_rows = read_matrix(CHRX_XY_MATRIX)

if len(auto_rows) != EXPECTED_AUTO_RECORDS:
    failures.append(f"Expected {EXPECTED_AUTO_RECORDS} autosomal matrix records; observed {len(auto_rows)}.")
if len(xx_rows) != EXPECTED_CHRX_RECORDS:
    failures.append(f"Expected {EXPECTED_CHRX_RECORDS} XX chrX matrix records; observed {len(xx_rows)}.")
if len(xy_rows) != EXPECTED_CHRX_RECORDS:
    failures.append(f"Expected {EXPECTED_CHRX_RECORDS} XY chrX matrix records; observed {len(xy_rows)}.")

auto_out = []
auto_total_an = auto_total_ac = auto_missing_cells = 0
auto_unexpected = 0

for row in auto_rows:
    gts = row["GT_values_in_frozen_4150_sample_order"].split(",")
    if len(gts) != EXPECTED_AUTO_SAMPLES:
        failures.append(f"{row['gene_symbol']} {row['allele_key']}: expected 4,150 autosomal GT cells; observed {len(gts)}.")
        continue

    an = ac = missing = 0
    for gt in gts:
        try:
            called, alt, miss = parse_autosomal_gt(gt)
        except ValueError:
            auto_unexpected += 1
            continue
        an += called
        ac += alt
        missing += miss

    af = ac / an if an else ""
    auto_total_an += an
    auto_total_ac += ac
    auto_missing_cells += missing

    auto_out.append({
        "analysis_layer": "autosomal",
        "gene_symbol": row["gene_symbol"],
        "allele_key": row["allele_key"],
        "sample_count": EXPECTED_AUTO_SAMPLES,
        "missing_GT_cells": missing,
        "AN": an,
        "AC": ac,
        "AF": af,
    })

with (OUTDIR / "phaseC_strict_autosomal_allele_count_frequency_qc_v1.tsv").open("w", newline="", encoding="utf-8") as handle:
    writer = csv.DictWriter(handle, fieldnames=list(auto_out[0]), delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(auto_out)

xx_by_key = {row["allele_key"]: row for row in xx_rows}
xy_by_key = {row["allele_key"]: row for row in xy_rows}

if set(xx_by_key) != set(xy_by_key):
    failures.append("XX and XY chrX matrices do not contain identical allele-key sets.")

chrx_out = []
chrx_total_an = chrx_total_ac = chrx_missing_cells = 0
chrx_unexpected = 0

for key in sorted(xx_by_key):
    xx = xx_by_key[key]
    xy = xy_by_key[key]

    xx_gts = xx["GT_values_in_frozen_XX_sample_order"].split(",")
    xy_gts = xy["GT_values_in_frozen_XY_sample_order"].split(",")

    if len(xx_gts) != EXPECTED_XX_SAMPLES:
        failures.append(f"IKBKG {key}: expected 1,904 XX GT cells; observed {len(xx_gts)}.")
        continue
    if len(xy_gts) != EXPECTED_XY_SAMPLES:
        failures.append(f"IKBKG {key}: expected 2,216 XY GT cells; observed {len(xy_gts)}.")
        continue

    xx_an = xx_ac = xx_missing = 0
    xy_an = xy_ac = xy_missing = 0

    for gt in xx_gts:
        try:
            called, alt, miss = parse_xx_gt(gt)
        except ValueError:
            chrx_unexpected += 1
            continue
        xx_an += called
        xx_ac += alt
        xx_missing += miss

    for gt in xy_gts:
        try:
            called, alt, miss = parse_xy_gt(gt)
        except ValueError:
            chrx_unexpected += 1
            continue
        xy_an += called
        xy_ac += alt
        xy_missing += miss

    an = xx_an + xy_an
    ac = xx_ac + xy_ac
    af = ac / an if an else ""

    chrx_total_an += an
    chrx_total_ac += ac
    chrx_missing_cells += xx_missing + xy_missing

    chrx_out.append({
        "analysis_layer": "chrX_sex_aware",
        "gene_symbol": "IKBKG",
        "allele_key": key,
        "XX_sample_count": EXPECTED_XX_SAMPLES,
        "XY_sample_count": EXPECTED_XY_SAMPLES,
        "XX_missing_GT_cells": xx_missing,
        "XY_missing_GT_cells": xy_missing,
        "XX_AN": xx_an,
        "XY_AN": xy_an,
        "AN": an,
        "XX_AC": xx_ac,
        "XY_AC": xy_ac,
        "AC": ac,
        "AF": af,
    })

with (OUTDIR / "phaseC_strict_IKBKG_chrX_allele_count_frequency_qc_v1.tsv").open("w", newline="", encoding="utf-8") as handle:
    writer = csv.DictWriter(handle, fieldnames=list(chrx_out[0]), delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(chrx_out)

summary_rows = [
    {
        "analysis_layer": "autosomal",
        "records": len(auto_out),
        "sample_model": "diploid_4150_samples",
        "total_AN": auto_total_an,
        "total_AC": auto_total_ac,
        "missing_GT_cells": auto_missing_cells,
        "unexpected_GT_cells": auto_unexpected,
    },
    {
        "analysis_layer": "chrX_sex_aware",
        "records": len(chrx_out),
        "sample_model": "diploid_XX_1904_plus_haploid_XY_2216",
        "total_AN": chrx_total_an,
        "total_AC": chrx_total_ac,
        "missing_GT_cells": chrx_missing_cells,
        "unexpected_GT_cells": chrx_unexpected,
    },
]

with (OUTDIR / "phaseC_strict_allele_count_frequency_qc_summary_v1.tsv").open("w", newline="", encoding="utf-8") as handle:
    writer = csv.DictWriter(handle, fieldnames=list(summary_rows[0]), delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(summary_rows)

# Hard denominator sanity checks from previous matrix summary.
expected_auto_an = (EXPECTED_AUTO_RECORDS * EXPECTED_AUTO_SAMPLES * 2) - (78697 * 2)
expected_chrx_an = ((EXPECTED_CHRX_RECORDS * EXPECTED_XX_SAMPLES * 2) - (45 * 2)) + (
    (EXPECTED_CHRX_RECORDS * EXPECTED_XY_SAMPLES) - 3940
)

if auto_total_an != expected_auto_an:
    failures.append(f"Autosomal total AN mismatch: expected {expected_auto_an}; observed {auto_total_an}.")
if chrx_total_an != expected_chrx_an:
    failures.append(f"chrX total AN mismatch: expected {expected_chrx_an}; observed {chrx_total_an}.")
if auto_unexpected != 0 or chrx_unexpected != 0:
    failures.append("Unexpected GT encodings observed during allele-count QC.")

payload = {
    "status": "PASS" if not failures else "FAIL",
    "operation": "Compute strict genotype-only allele-count/frequency QC tables.",
    "record_counts": {
        "autosomal": len(auto_out),
        "chrX_sex_aware": len(chrx_out),
    },
    "denominator_model": {
        "autosomal": "AN = 2 * called diploid autosomal genotypes",
        "chrX_sex_aware": "AN = 2 * called XX chrX genotypes + 1 * called XY chrX genotypes",
    },
    "summary": summary_rows,
    "failures": failures,
    "boundary": (
        "AC/AN/AF are computed for technical QC only. No population comparison, selection inference, "
        "functional prioritisation or biological interpretation is authorised."
    ),
}

OUTJSON.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
print(json.dumps(payload, indent=2))

if payload["status"] != "PASS":
    raise SystemExit(1)
