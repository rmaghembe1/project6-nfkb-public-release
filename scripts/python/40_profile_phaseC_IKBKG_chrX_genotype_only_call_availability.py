#!/usr/bin/env python3
"""Profile Phase C genotype-only IKBKG/chrX call availability using sex-aware primary-cohort handling."""

from __future__ import annotations

import csv
import json
import re
import subprocess
from collections import Counter
from pathlib import Path

META = Path("02_reference_data/gnomad_HGDP1kGP_metadata/gnomad_meta_updated.tsv")
COHORT = Path(
    "11_results/phaseB/harmonisation_qc/"
    "cohort_primary_genotype_metadata_linked_n4150_v1.txt"
)
IKBKG = Path(
    "03_raw_genomic_resources/interval_extraction/full_panel/genotype/"
    "IKBKG_GRCh38_gnomad_hgdp_tgp.vcf.gz"
)
ROOT = Path("11_results/phaseC/genotype_only_qc")
OUTDIR = ROOT / "chrX_IKBKG_call_availability"
OUTJSON = Path("results/qc/phaseC_IKBKG_chrX_genotype_only_call_availability_profile_v1.json")

SEX_META = ROOT / "phaseC_genotype_only_cohort_sex_ploidy_metadata_v1.tsv"
XX_LIST = ROOT / "phaseC_IKBKG_candidate_primary_standard_XX_unflagged_samples_v1.txt"
XY_LIST = ROOT / "phaseC_IKBKG_candidate_primary_standard_XY_unflagged_samples_v1.txt"
REVIEW_REGISTER = ROOT / "phaseC_IKBKG_nonstandard_or_ambiguous_chrX_samples_v1.tsv"
EXACT_GT = ROOT / "phaseC_IKBKG_primary_genotype_exact_GT_string_counts_by_karyotype_v1.tsv"

EXPECTED_KARYOTYPES = {
    "XX": 1904,
    "XY": 2216,
    "X": 20,
    "ambiguous": 6,
    "XXY": 2,
    "XYY": 2,
}
EXPECTED_PRIMARY = 4120
EXPECTED_REVIEW = 30
EXPECTED_RECORDS = 937

SEX_FIELDS = [
    "sample_id", "sex_karyotype", "chrX_ploidy", "chrY_ploidy",
    "sex_aneuploidy", "ambiguous_sex", "reported_sex"
]


def variant_class(key: str) -> str:
    _, _, ref, alt = key.split(":", 3)
    if len(ref) == 1 and len(alt) == 1:
        return "SNV"
    if len(ref) != len(alt):
        return "INDEL"
    return "MNV_OR_COMPLEX"


def query_records(sample_file: Path) -> list[tuple[str, list[str]]]:
    cmd = [
        "bcftools", "query", "-S", str(sample_file),
        "-f", "%CHROM:%POS:%REF:%ALT[\t%GT]\n", str(IKBKG)
    ]
    run = subprocess.run(cmd, check=True, text=True, capture_output=True)
    records = []
    for line in run.stdout.splitlines():
        fields = line.split("\t")
        records.append((fields[0], fields[1:]))
    return records


OUTDIR.mkdir(parents=True, exist_ok=True)
OUTJSON.parent.mkdir(parents=True, exist_ok=True)
failures: list[str] = []

cohort = [x for x in COHORT.read_text(encoding="utf-8").splitlines() if x]
cohort_set = set(cohort)
if len(cohort) != 4150 or len(cohort_set) != 4150:
    failures.append("Frozen genotype-only cohort is not exactly 4,150 unique samples.")

with META.open(newline="", encoding="utf-8") as handle:
    reader = csv.DictReader(handle, delimiter="\t")
    metadata = {row["s"]: row for row in reader if row["s"] in cohort_set}

if len(metadata) != 4150:
    failures.append(f"Expected metadata for 4,150 cohort samples; observed {len(metadata)}.")

sex_rows = []
for sample in cohort:
    row = metadata[sample]
    sex_rows.append({
        "sample_id": sample,
        "sex_karyotype": row["sex_imputation.sex_karyotype"],
        "chrX_ploidy": row["sex_imputation.chrX_ploidy"],
        "chrY_ploidy": row["sex_imputation.chrY_ploidy"],
        "sex_aneuploidy": row["sample_filters.sex_aneuploidy"],
        "ambiguous_sex": row["sample_filters.ambiguous_sex"],
        "reported_sex": row["sex"],
    })

with SEX_META.open("w", newline="", encoding="utf-8") as handle:
    writer = csv.DictWriter(handle, fieldnames=SEX_FIELDS, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(sex_rows)

karyotype_counts = Counter(row["sex_karyotype"] for row in sex_rows)
if dict(karyotype_counts) != EXPECTED_KARYOTYPES:
    failures.append(f"Unexpected sex-karyotype counts: {dict(karyotype_counts)}.")

xx = [
    row["sample_id"] for row in sex_rows
    if row["sex_karyotype"] == "XX"
    and row["sex_aneuploidy"] == "false"
    and row["ambiguous_sex"] == "false"
]
xy = [
    row["sample_id"] for row in sex_rows
    if row["sex_karyotype"] == "XY"
    and row["sex_aneuploidy"] == "false"
    and row["ambiguous_sex"] == "false"
]
review = [row for row in sex_rows if row["sex_karyotype"] not in {"XX", "XY"}]

if len(xx) != 1904 or len(xy) != 2216 or len(xx) + len(xy) != EXPECTED_PRIMARY:
    failures.append("Primary standard unflagged XX/XY chrX cohort does not contain 4,120 samples.")
if len(review) != EXPECTED_REVIEW:
    failures.append("Technical-review chrX register does not contain 30 samples.")

XX_LIST.write_text("\n".join(xx) + "\n", encoding="utf-8")
XY_LIST.write_text("\n".join(xy) + "\n", encoding="utf-8")

with REVIEW_REGISTER.open("w", newline="", encoding="utf-8") as handle:
    writer = csv.DictWriter(handle, fieldnames=SEX_FIELDS, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(review)

records = {"XX": query_records(XX_LIST), "XY": query_records(XY_LIST)}
if len(records["XX"]) != EXPECTED_RECORDS or len(records["XY"]) != EXPECTED_RECORDS:
    failures.append("IKBKG query does not contain the expected 937 records.")
if [key for key, _ in records["XX"]] != [key for key, _ in records["XY"]]:
    failures.append("IKBKG allele-key order differs between XX and XY subsets.")

exact_counts = {"XX": Counter(), "XY": Counter()}
site_rows = []
unexpected_total = 0

for i, (key, xx_gts) in enumerate(records["XX"]):
    xy_key, xy_gts = records["XY"][i]
    if xy_key != key:
        failures.append(f"Allele-key mismatch at record {i}: {key} versus {xy_key}.")
        continue

    row = {"allele_key": key, "variant_class": variant_class(key)}

    for stratum, gts, expected in [("XX", xx_gts, len(xx)), ("XY", xy_gts, len(xy))]:
        if len(gts) != expected:
            failures.append(f"{key} {stratum}: expected {expected} GT cells, observed {len(gts)}.")

        called = missing = unexpected = 0
        for gt in gts:
            exact_counts[stratum][gt] += 1
            if "." in gt:
                missing += 1
            elif stratum == "XX" and re.fullmatch(r"\d+/\d+", gt):
                called += 1
            elif stratum == "XY" and re.fullmatch(r"\d+", gt):
                called += 1
            else:
                unexpected += 1

        unexpected_total += unexpected
        row[f"{stratum}_expected_GT_cells"] = expected
        row[f"{stratum}_called_GT_cells"] = called
        row[f"{stratum}_missing_GT_cells"] = missing
        row[f"{stratum}_unexpected_nonmissing_GT_cells"] = unexpected
        row[f"{stratum}_missingness_fraction"] = missing / expected

    row["combined_expected_GT_cells"] = len(xx) + len(xy)
    row["combined_called_GT_cells"] = row["XX_called_GT_cells"] + row["XY_called_GT_cells"]
    row["combined_missing_GT_cells"] = row["XX_missing_GT_cells"] + row["XY_missing_GT_cells"]
    row["combined_unexpected_nonmissing_GT_cells"] = (
        row["XX_unexpected_nonmissing_GT_cells"] + row["XY_unexpected_nonmissing_GT_cells"]
    )
    row["combined_missingness_fraction"] = (
        row["combined_missing_GT_cells"] / row["combined_expected_GT_cells"]
    )
    site_rows.append(row)

with EXACT_GT.open("w", newline="", encoding="utf-8") as handle:
    writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
    writer.writerow(["sex_karyotype", "GT_string", "count"])
    for stratum in ["XX", "XY"]:
        for gt, count in sorted(exact_counts[stratum].items()):
            writer.writerow([stratum, gt, count])

profile = OUTDIR / "phaseC_IKBKG_genotype_only_call_availability_by_allele_key_v1.tsv"
with profile.open("w", newline="", encoding="utf-8") as handle:
    writer = csv.DictWriter(handle, fieldnames=list(site_rows[0]), delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(site_rows)

candidate_rules = [
    ("combined_missingness_gt_0.001", lambda r: r["combined_missingness_fraction"] > 0.001),
    ("combined_missingness_gt_0.01", lambda r: r["combined_missingness_fraction"] > 0.01),
    ("either_stratum_missingness_gt_0.001", lambda r: r["XX_missingness_fraction"] > 0.001 or r["XY_missingness_fraction"] > 0.001),
    ("either_stratum_missingness_gt_0.01", lambda r: r["XX_missingness_fraction"] > 0.01 or r["XY_missingness_fraction"] > 0.01),
]

threshold_rows = []
for label, rule in candidate_rules:
    flagged = sum(rule(row) for row in site_rows)
    threshold_rows.append({
        "candidate_rule_component": label,
        "flagged_allele_keys": flagged,
        "retained_allele_keys_if_excluded": len(site_rows) - flagged,
    })

threshold_file = OUTDIR / "phaseC_IKBKG_candidate_missingness_threshold_counts_v1.tsv"
with threshold_file.open("w", newline="", encoding="utf-8") as handle:
    writer = csv.DictWriter(handle, fieldnames=list(threshold_rows[0]), delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(threshold_rows)

total_expected = sum(int(row["combined_expected_GT_cells"]) for row in site_rows)
total_missing = sum(int(row["combined_missing_GT_cells"]) for row in site_rows)

payload = {
    "status": "PASS" if not failures and unexpected_total == 0 else "FAIL",
    "operation": "Descriptive sex-aware genotype-only IKBKG/chrX GT call-availability profiling only.",
    "full_genotype_only_cohort_count": len(cohort),
    "sex_karyotype_counts": dict(karyotype_counts),
    "primary_chrX_QC_cohort": {"XX": len(xx), "XY": len(xy), "total": len(xx) + len(xy)},
    "chrX_technical_review_sample_count": len(review),
    "chrX_technical_review_karyotype_counts": dict(Counter(row["sex_karyotype"] for row in review)),
    "IKBKG_variant_records_profiled": len(site_rows),
    "expected_GT_cells_in_primary_chrX_QC_cohort": total_expected,
    "missing_GT_cells": total_missing,
    "unexpected_nonmissing_GT_cells": unexpected_total,
    "overall_missingness_fraction": total_missing / total_expected,
    "failures": failures,
    "boundary": (
        "No genotype-only missingness rule, allele-frequency calculation, "
        "population comparison or biological interpretation is authorised."
    ),
}
OUTJSON.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
print(json.dumps(payload, indent=2))

if payload["status"] != "PASS":
    raise SystemExit(1)
