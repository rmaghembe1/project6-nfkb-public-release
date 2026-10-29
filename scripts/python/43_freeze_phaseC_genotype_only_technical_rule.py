#!/usr/bin/env python3
"""Freeze the Phase C genotype-only technical inclusion rule from reviewed candidate inventories."""

from __future__ import annotations

import csv
import json
import shutil
from collections import Counter
from pathlib import Path

CANDIDATE_DIR = Path("11_results/phaseC/genotype_only_qc/rule_design")
OUTDIR = Path("11_results/phaseC/genotype_only_qc/final_technical_rule")
OUTJSON = Path("results/qc/phaseC_genotype_only_technical_rule_freeze_v1.json")

STRICT_IN = CANDIDATE_DIR / "phaseC_candidate_strict_genotype_only_records_FILTER_PASS_missingness_rule_v1.tsv"
REVIEW_IN = CANDIDATE_DIR / "phaseC_candidate_technical_review_genotype_only_records_FILTER_PASS_missingness_rule_v1.tsv"
REASON_IN = CANDIDATE_DIR / "phaseC_candidate_technical_review_reason_counts_v1.tsv"

STRICT_OUT = OUTDIR / "phaseC_strict_primary_genotype_only_records_v1.tsv"
REVIEW_OUT = OUTDIR / "phaseC_technical_review_genotype_only_records_v1.tsv"
REASON_OUT = OUTDIR / "phaseC_technical_review_reason_counts_v1.tsv"
SUMMARY_OUT = OUTDIR / "phaseC_genotype_only_technical_rule_summary_v1.tsv"

EXPECTED_TOTAL = 56529
EXPECTED_STRICT = 39821
EXPECTED_REVIEW = 16708
EXPECTED_AUTOSOMAL_TOTAL = 55592
EXPECTED_CHRX_TOTAL = 937
EXPECTED_CHRX_STRICT = 383
EXPECTED_AUTOSOMAL_STRICT = 39438

OUTDIR.mkdir(parents=True, exist_ok=True)
OUTJSON.parent.mkdir(parents=True, exist_ok=True)


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


strict = read_tsv(STRICT_IN)
review = read_tsv(REVIEW_IN)
reason = read_tsv(REASON_IN)

failures: list[str] = []

if len(strict) != EXPECTED_STRICT:
    failures.append(f"Expected {EXPECTED_STRICT} strict records; observed {len(strict)}.")
if len(review) != EXPECTED_REVIEW:
    failures.append(f"Expected {EXPECTED_REVIEW} review records; observed {len(review)}.")
if len(strict) + len(review) != EXPECTED_TOTAL:
    failures.append(f"Expected total {EXPECTED_TOTAL}; observed {len(strict) + len(review)}.")

strict_keys = {(r["analysis_layer"], r["gene_symbol"], r["allele_key"]) for r in strict}
review_keys = {(r["analysis_layer"], r["gene_symbol"], r["allele_key"]) for r in review}

if strict_keys & review_keys:
    failures.append("Strict and technical-review inventories are not disjoint.")

all_keys = strict_keys | review_keys
if len(all_keys) != EXPECTED_TOTAL:
    failures.append("Combined strict/review key inventory is not unique at expected total size.")

strict_layer_counts = Counter(r["analysis_layer"] for r in strict)
review_layer_counts = Counter(r["analysis_layer"] for r in review)

if strict_layer_counts.get("autosomal", 0) != EXPECTED_AUTOSOMAL_STRICT:
    failures.append("Autosomal strict count does not match expected 39,438.")
if strict_layer_counts.get("chrX_sex_aware", 0) != EXPECTED_CHRX_STRICT:
    failures.append("chrX strict count does not match expected 383.")
if strict_layer_counts.get("autosomal", 0) + review_layer_counts.get("autosomal", 0) != EXPECTED_AUTOSOMAL_TOTAL:
    failures.append("Autosomal strict plus review count does not equal 55,592.")
if strict_layer_counts.get("chrX_sex_aware", 0) + review_layer_counts.get("chrX_sex_aware", 0) != EXPECTED_CHRX_TOTAL:
    failures.append("chrX strict plus review count does not equal 937.")

# Validate the frozen strict rule itself.
for row in strict:
    if row["site_FILTER"] != "PASS":
        failures.append(f"Strict record has non-PASS FILTER: {row['analysis_layer']} {row['allele_key']}")
        break
    if row["analysis_layer"] == "autosomal":
        if float(row["missingness_fraction"]) > 0.01:
            failures.append(f"Autosomal strict record exceeds missingness rule: {row['allele_key']}")
            break
    elif row["analysis_layer"] == "chrX_sex_aware":
        if float(row["XX_missingness_fraction"]) > 0.01 or float(row["XY_missingness_fraction"]) > 0.01:
            failures.append(f"chrX strict record exceeds sex-stratum missingness rule: {row['allele_key']}")
            break
    else:
        failures.append(f"Unexpected analysis layer in strict inventory: {row['analysis_layer']}")
        break

# Validate review records have reason codes.
for row in review:
    if not row.get("technical_review_reason", ""):
        failures.append(f"Review record lacks reason code: {row['analysis_layer']} {row['allele_key']}")
        break

shutil.copyfile(STRICT_IN, STRICT_OUT)
shutil.copyfile(REVIEW_IN, REVIEW_OUT)
shutil.copyfile(REASON_IN, REASON_OUT)

summary_rows = [
    {
        "inventory": "strict_primary",
        "analysis_layer": "autosomal",
        "record_count": strict_layer_counts.get("autosomal", 0),
        "rule": "FILTER=PASS and missingness_fraction <= 0.01",
    },
    {
        "inventory": "strict_primary",
        "analysis_layer": "chrX_sex_aware",
        "record_count": strict_layer_counts.get("chrX_sex_aware", 0),
        "rule": "FILTER=PASS and XX_missingness_fraction <= 0.01 and XY_missingness_fraction <= 0.01",
    },
    {
        "inventory": "technical_review",
        "analysis_layer": "autosomal",
        "record_count": review_layer_counts.get("autosomal", 0),
        "rule": "Preserved for technical review/sensitivity with reason codes",
    },
    {
        "inventory": "technical_review",
        "analysis_layer": "chrX_sex_aware",
        "record_count": review_layer_counts.get("chrX_sex_aware", 0),
        "rule": "Preserved for technical review/sensitivity with reason codes",
    },
]

with SUMMARY_OUT.open("w", newline="", encoding="utf-8") as handle:
    writer = csv.DictWriter(handle, fieldnames=list(summary_rows[0]), delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(summary_rows)

payload = {
    "status": "PASS" if not failures else "FAIL",
    "operation": "Freeze Phase C genotype-only technical inclusion rule.",
    "frozen_rule": {
        "autosomal_strict_primary": "FILTER=PASS and missingness_fraction <= 0.01",
        "chrX_IKBKG_strict_primary": (
            "FILTER=PASS and XX_missingness_fraction <= 0.01 and "
            "XY_missingness_fraction <= 0.01"
        ),
        "technical_review": "All non-PASS or high-missingness records preserved with reason codes.",
    },
    "record_counts": {
        "strict_primary_total": len(strict),
        "strict_primary_autosomal": strict_layer_counts.get("autosomal", 0),
        "strict_primary_chrX_sex_aware": strict_layer_counts.get("chrX_sex_aware", 0),
        "technical_review_total": len(review),
        "technical_review_autosomal": review_layer_counts.get("autosomal", 0),
        "technical_review_chrX_sex_aware": review_layer_counts.get("chrX_sex_aware", 0),
        "combined_total": len(strict) + len(review),
    },
    "technical_review_reason_counts": {
        row["technical_review_reason"]: int(row["record_count"]) for row in reason
    },
    "failures": failures,
    "boundary": (
        "This freezes technical record inclusion inventories only. It does not compute allele frequencies, "
        "perform population comparisons, infer selection, prioritise function or authorise biological interpretation."
    ),
}

OUTJSON.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
print(json.dumps(payload, indent=2))

if payload["status"] != "PASS":
    raise SystemExit(1)
