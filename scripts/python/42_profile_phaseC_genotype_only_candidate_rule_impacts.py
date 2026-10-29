#!/usr/bin/env python3
"""Profile candidate genotype-only technical inclusion rules before freezing a rule."""

from __future__ import annotations

import csv
import json
from collections import Counter, defaultdict
from pathlib import Path

JOINED = Path(
    "11_results/phaseC/genotype_only_qc/site_quality_context/"
    "phaseC_genotype_only_filter_missingness_joined_inventory_v1.tsv"
)
CHRX_PROFILE = Path(
    "11_results/phaseC/genotype_only_qc/chrX_IKBKG_call_availability/"
    "phaseC_IKBKG_genotype_only_call_availability_by_allele_key_v1.tsv"
)
OUTDIR = Path("11_results/phaseC/genotype_only_qc/rule_design")
OUTJSON = Path("results/qc/phaseC_genotype_only_candidate_rule_impact_profile_v1.json")

OUTDIR.mkdir(parents=True, exist_ok=True)
OUTJSON.parent.mkdir(parents=True, exist_ok=True)


def load_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


joined = load_tsv(JOINED)
chrx = {row["allele_key"]: row for row in load_tsv(CHRX_PROFILE)}

failures: list[str] = []

if len(joined) != 56529:
    failures.append(f"Expected 56,529 joined records; observed {len(joined)}.")

autosomal = [r for r in joined if r["analysis_layer"] == "autosomal"]
chrx_rows = [r for r in joined if r["analysis_layer"] == "chrX_sex_aware"]

if len(autosomal) != 55592:
    failures.append(f"Expected 55,592 autosomal records; observed {len(autosomal)}.")
if len(chrx_rows) != 937:
    failures.append(f"Expected 937 chrX records; observed {len(chrx_rows)}.")

candidate_rows = []

rules = [
    ("FILTER_PASS_only", "site_FILTER == PASS"),
    ("FILTER_PASS_and_missingness_le_0.05", "site_FILTER == PASS and combined/autosomal missingness <= 0.05"),
    ("FILTER_PASS_and_missingness_le_0.01", "site_FILTER == PASS and combined/autosomal missingness <= 0.01"),
    ("FILTER_PASS_and_missingness_le_0.001", "site_FILTER == PASS and combined/autosomal missingness <= 0.001"),
    ("FILTER_PASS_SNV_missingness_le_0.01_INDEL_missingness_le_0.05",
     "site_FILTER == PASS and SNV missingness <= 0.01 and INDEL/MNV/complex missingness <= 0.05"),
]

for rule_name, description in rules:
    for layer, records in [("autosomal", autosomal), ("chrX_sex_aware", chrx_rows)]:
        retained = 0
        review = 0

        for row in records:
            miss = float(row["missingness_fraction"])
            filt = row["site_FILTER"]
            klass = row["variant_class"]

            keep = False
            if rule_name == "FILTER_PASS_only":
                keep = filt == "PASS"
            elif rule_name == "FILTER_PASS_and_missingness_le_0.05":
                keep = filt == "PASS" and miss <= 0.05
            elif rule_name == "FILTER_PASS_and_missingness_le_0.01":
                keep = filt == "PASS" and miss <= 0.01
            elif rule_name == "FILTER_PASS_and_missingness_le_0.001":
                keep = filt == "PASS" and miss <= 0.001
            elif rule_name == "FILTER_PASS_SNV_missingness_le_0.01_INDEL_missingness_le_0.05":
                if filt == "PASS":
                    if klass == "SNV":
                        keep = miss <= 0.01
                    else:
                        keep = miss <= 0.05

            if keep:
                retained += 1
            else:
                review += 1

        candidate_rows.append({
            "candidate_rule": rule_name,
            "rule_description": description,
            "analysis_layer": layer,
            "input_records": len(records),
            "retained_records": retained,
            "technical_review_records": review,
            "retained_fraction": retained / len(records) if records else "",
        })

# Additional chrX sex-stratum-aware candidate rules.
chrx_sex_rules = [
    ("chrX_FILTER_PASS_and_combined_missingness_le_0.01",
     lambda r: r["site_FILTER"] == "PASS" and float(r["missingness_fraction"]) <= 0.01),
    ("chrX_FILTER_PASS_and_XX_XY_missingness_each_le_0.01",
     lambda r: r["site_FILTER"] == "PASS"
     and float(chrx[r["allele_key"]]["XX_missingness_fraction"]) <= 0.01
     and float(chrx[r["allele_key"]]["XY_missingness_fraction"]) <= 0.01),
    ("chrX_FILTER_PASS_and_XX_XY_missingness_each_le_0.05",
     lambda r: r["site_FILTER"] == "PASS"
     and float(chrx[r["allele_key"]]["XX_missingness_fraction"]) <= 0.05
     and float(chrx[r["allele_key"]]["XY_missingness_fraction"]) <= 0.05),
]

for rule_name, predicate in chrx_sex_rules:
    retained = sum(predicate(row) for row in chrx_rows)
    candidate_rows.append({
        "candidate_rule": rule_name,
        "rule_description": "IKBKG chrX sex-stratum-aware candidate rule",
        "analysis_layer": "chrX_sex_aware",
        "input_records": len(chrx_rows),
        "retained_records": retained,
        "technical_review_records": len(chrx_rows) - retained,
        "retained_fraction": retained / len(chrx_rows),
    })

with (OUTDIR / "phaseC_genotype_only_candidate_rule_impact_summary_v1.tsv").open("w", newline="", encoding="utf-8") as handle:
    writer = csv.DictWriter(handle, fieldnames=list(candidate_rows[0]), delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(candidate_rows)

# Build reason-coded inventories for the most conservative interpretable candidate:
# Autosomal: FILTER=PASS and missingness <= 0.01.
# chrX: FILTER=PASS and XX and XY missingness each <= 0.01.
strict_rows = []
review_rows = []

for row in autosomal:
    miss = float(row["missingness_fraction"])
    reasons = []
    if row["site_FILTER"] != "PASS":
        reasons.append("NON_PASS_FILTER")
    if miss > 0.01:
        reasons.append("MISSINGNESS_GT_0.01")

    out = {
        "analysis_layer": row["analysis_layer"],
        "gene_symbol": row["gene_symbol"],
        "allele_key": row["allele_key"],
        "variant_class": row["variant_class"],
        "site_FILTER": row["site_FILTER"],
        "missingness_fraction": row["missingness_fraction"],
    }
    if reasons:
        out["technical_review_reason"] = ";".join(reasons)
        review_rows.append(out)
    else:
        strict_rows.append(out)

for row in chrx_rows:
    key = row["allele_key"]
    sexrow = chrx[key]
    xx_miss = float(sexrow["XX_missingness_fraction"])
    xy_miss = float(sexrow["XY_missingness_fraction"])
    reasons = []
    if row["site_FILTER"] != "PASS":
        reasons.append("NON_PASS_FILTER")
    if xx_miss > 0.01:
        reasons.append("XX_MISSINGNESS_GT_0.01")
    if xy_miss > 0.01:
        reasons.append("XY_MISSINGNESS_GT_0.01")

    out = {
        "analysis_layer": row["analysis_layer"],
        "gene_symbol": row["gene_symbol"],
        "allele_key": row["allele_key"],
        "variant_class": row["variant_class"],
        "site_FILTER": row["site_FILTER"],
        "combined_missingness_fraction": row["missingness_fraction"],
        "XX_missingness_fraction": sexrow["XX_missingness_fraction"],
        "XY_missingness_fraction": sexrow["XY_missingness_fraction"],
    }
    if reasons:
        out["technical_review_reason"] = ";".join(reasons)
        review_rows.append(out)
    else:
        strict_rows.append(out)

strict_path = OUTDIR / "phaseC_candidate_strict_genotype_only_records_FILTER_PASS_missingness_rule_v1.tsv"
review_path = OUTDIR / "phaseC_candidate_technical_review_genotype_only_records_FILTER_PASS_missingness_rule_v1.tsv"

strict_fields = sorted({k for row in strict_rows for k in row})
review_fields = sorted({k for row in review_rows for k in row})

with strict_path.open("w", newline="", encoding="utf-8") as handle:
    writer = csv.DictWriter(handle, fieldnames=strict_fields, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(strict_rows)

with review_path.open("w", newline="", encoding="utf-8") as handle:
    writer = csv.DictWriter(handle, fieldnames=review_fields, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(review_rows)

reason_counts = Counter(row["technical_review_reason"] for row in review_rows)
with (OUTDIR / "phaseC_candidate_technical_review_reason_counts_v1.tsv").open("w", newline="", encoding="utf-8") as handle:
    writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
    writer.writerow(["technical_review_reason", "record_count"])
    for reason, count in sorted(reason_counts.items(), key=lambda x: (-x[1], x[0])):
        writer.writerow([reason, count])

payload = {
    "status": "PASS" if not failures else "FAIL",
    "operation": "Candidate genotype-only technical inclusion rule-impact profiling; no rule is frozen.",
    "records_profiled": {
        "autosomal": len(autosomal),
        "chrX_sex_aware": len(chrx_rows),
        "combined": len(joined),
    },
    "candidate_default_strict_rule": {
        "autosomal": "FILTER=PASS and missingness <= 0.01",
        "chrX_sex_aware": "FILTER=PASS and XX_missingness <= 0.01 and XY_missingness <= 0.01",
        "strict_records": len(strict_rows),
        "technical_review_records": len(review_rows),
    },
    "technical_review_reason_counts": dict(reason_counts),
    "failures": failures,
    "boundary": (
        "These are candidate rule-impact outputs only. No genotype-only technical rule, "
        "frequency calculation, population comparison or biological interpretation is authorised."
    ),
}

OUTJSON.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
print(json.dumps(payload, indent=2))

if payload["status"] != "PASS":
    raise SystemExit(1)
