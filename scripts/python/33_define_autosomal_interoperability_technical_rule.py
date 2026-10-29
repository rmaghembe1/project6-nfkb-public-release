#!/usr/bin/env python3
"""Define analysis-specific autosomal genotype-phased interoperability technical sets."""

from __future__ import annotations

import csv
import json
from collections import Counter
from pathlib import Path

INPUT = Path(
    "11_results/phaseB/genotype_call_concordance_qc/autosomal_filter_rule_development/"
    "autosomal_shared_allele_key_discordance_missingness_profile_v1.tsv"
)
OUTDIR = Path(
    "11_results/phaseB/genotype_call_concordance_qc/autosomal_filter_rule"
)
OUTJSON = Path(
    "results/qc/phaseB_autosomal_interoperability_technical_rule_v1.json"
)

MISSINGNESS_MAX = 0.01

with INPUT.open(newline="", encoding="utf-8") as handle:
    rows = list(csv.DictReader(handle, delimiter="\t"))

strict_rows = []
review_rows = []
reason_counts = Counter()

for row in rows:
    discordant = int(row["discordant_GT_calls"])
    missingness = float(row["missingness_fraction_of_expected_cells"])

    reasons = []
    if discordant > 0:
        reasons.append("ANY_GT_DISCORDANCE")
    if missingness > MISSINGNESS_MAX:
        reasons.append("MISSINGNESS_GT_0.01")

    out = dict(row)
    out["technical_rule_status"] = "STRICT_INTEROPERABILITY_RETAIN" if not reasons else "TECHNICAL_REVIEW_SET"
    out["technical_rule_reason"] = "PASS" if not reasons else ";".join(reasons)

    if reasons:
        review_rows.append(out)
        if len(reasons) == 2:
            reason_counts["BOTH_DISCORDANCE_AND_MISSINGNESS_GT_0.01"] += 1
        elif reasons[0] == "ANY_GT_DISCORDANCE":
            reason_counts["DISCORDANCE_ONLY"] += 1
        else:
            reason_counts["MISSINGNESS_GT_0.01_ONLY"] += 1
    else:
        strict_rows.append(out)

OUTDIR.mkdir(parents=True, exist_ok=True)

fields = list(strict_rows[0].keys())

with (OUTDIR / "autosomal_strict_interoperability_retained_allele_keys_v1.tsv").open(
    "w", newline="", encoding="utf-8"
) as handle:
    writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(strict_rows)

with (OUTDIR / "autosomal_interoperability_technical_review_allele_keys_v1.tsv").open(
    "w", newline="", encoding="utf-8"
) as handle:
    writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(review_rows)

with (OUTDIR / "autosomal_interoperability_rule_summary_v1.tsv").open(
    "w", newline="", encoding="utf-8"
) as handle:
    writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
    writer.writerow(["category", "allele_key_count"])
    writer.writerow(["shared_autosomal_input_keys", len(rows)])
    writer.writerow(["strict_interoperability_retained_keys", len(strict_rows)])
    writer.writerow(["technical_review_set_keys", len(review_rows)])
    writer.writerow(["discordance_only_review_keys", reason_counts["DISCORDANCE_ONLY"]])
    writer.writerow(["missingness_gt_0.01_only_review_keys", reason_counts["MISSINGNESS_GT_0.01_ONLY"]])
    writer.writerow(["both_discordance_and_missingness_gt_0.01_review_keys", reason_counts["BOTH_DISCORDANCE_AND_MISSINGNESS_GT_0.01"]])

payload = {
    "status": "PASS",
    "operation": "Definition of autosomal phased/interoperability technical key sets only.",
    "rule_scope": "Shared autosomal genotype-phased allele keys in the frozen 4091-sample interoperability cohort.",
    "strict_interoperability_rule": {
        "discordant_GT_calls": "must equal 0",
        "missingness_fraction_of_expected_cells": "must be <= 0.01",
    },
    "shared_autosomal_input_keys": len(rows),
    "strict_interoperability_retained_keys": len(strict_rows),
    "technical_review_set_keys": len(review_rows),
    "review_reason_counts": dict(reason_counts),
    "excluded_from_rule_scope": ["IKBKG_chrX", "genotype_only_analysis_layer"],
    "policy_boundary": (
        "This strict rule governs autosomal genotype-phased interoperability and later phased analyses only. "
        "It must not be applied as a quality exclusion rule for genotype-only analyses."
    ),
    "interpretation_boundary": "No biological, allele-frequency, haplotype or selection interpretation is authorised.",
}
OUTJSON.parent.mkdir(parents=True, exist_ok=True)
OUTJSON.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
print(json.dumps(payload, indent=2))

if len(rows) != 22846 or len(strict_rows) != 20553 or len(review_rows) != 2293:
    raise SystemExit("ERROR: technical-set counts do not match the profiled candidate-rule result.")
