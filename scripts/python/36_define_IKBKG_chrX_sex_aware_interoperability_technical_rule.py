#!/usr/bin/env python3
"""Define the sex-aware IKBKG chrX phased-interoperability technical key sets."""

from __future__ import annotations

import csv
import json
from collections import Counter
from pathlib import Path

INPUT = Path(
    "11_results/phaseB/genotype_call_concordance_qc/chrX_IKBKG_filter_rule_development/"
    "IKBKG_chrX_shared_allele_key_discordance_unavailable_profile_v1.tsv"
)
OUTDIR = Path(
    "11_results/phaseB/genotype_call_concordance_qc/chrX_IKBKG_filter_rule"
)
OUTJSON = Path("results/qc/phaseB_IKBKG_chrX_interoperability_technical_rule_v1.json")

MAX_UNAVAILABLE_FRACTION = 0.01

with INPUT.open(newline="", encoding="utf-8") as handle:
    rows = list(csv.DictReader(handle, delimiter="\t"))

strict_rows = []
review_rows = []
reason_pattern_counts = Counter()

for row in rows:
    discordant = int(row["combined_discordant_GT_calls"])
    xx_unavailable = float(row["XX_unavailable_fraction"])
    xy_unavailable = float(row["XY_unavailable_fraction"])

    reasons = []
    if discordant > 0:
        reasons.append("ANY_GT_DISCORDANCE")
    if xx_unavailable > MAX_UNAVAILABLE_FRACTION:
        reasons.append("XX_UNAVAILABLE_GT_0.01")
    if xy_unavailable > MAX_UNAVAILABLE_FRACTION:
        reasons.append("XY_UNAVAILABLE_GT_0.01")

    out = dict(row)
    if reasons:
        out["technical_rule_status"] = "TECHNICAL_REVIEW_SET"
        out["technical_rule_reason"] = ";".join(reasons)
        review_rows.append(out)
        reason_pattern_counts[";".join(reasons)] += 1
    else:
        out["technical_rule_status"] = "STRICT_INTEROPERABILITY_RETAIN"
        out["technical_rule_reason"] = "PASS"
        strict_rows.append(out)

OUTDIR.mkdir(parents=True, exist_ok=True)
fields = list(strict_rows[0].keys())

strict_path = OUTDIR / "IKBKG_chrX_strict_interoperability_retained_allele_keys_v1.tsv"
review_path = OUTDIR / "IKBKG_chrX_interoperability_technical_review_allele_keys_v1.tsv"
summary_path = OUTDIR / "IKBKG_chrX_interoperability_rule_summary_v1.tsv"
reason_path = OUTDIR / "IKBKG_chrX_interoperability_review_reason_patterns_v1.tsv"

with strict_path.open("w", newline="", encoding="utf-8") as handle:
    writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(strict_rows)

with review_path.open("w", newline="", encoding="utf-8") as handle:
    writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(review_rows)

with summary_path.open("w", newline="", encoding="utf-8") as handle:
    writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
    writer.writerow(["category", "allele_key_count"])
    writer.writerow(["shared_IKBKG_input_keys", len(rows)])
    writer.writerow(["strict_interoperability_retained_keys", len(strict_rows)])
    writer.writerow(["technical_review_set_keys", len(review_rows)])

with reason_path.open("w", newline="", encoding="utf-8") as handle:
    writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
    writer.writerow(["technical_rule_reason_pattern", "allele_key_count"])
    for reason, count in sorted(reason_pattern_counts.items()):
        writer.writerow([reason, count])

payload = {
    "status": "PASS",
    "operation": "Definition of IKBKG chrX sex-aware phased/interoperability technical key sets only.",
    "rule_scope": "Shared IKBKG genotype-phased allele keys evaluated separately in clean XX and XY strata.",
    "strict_interoperability_rule": {
        "combined_discordant_GT_calls": "must equal 0",
        "XX_unavailable_fraction": "must be <= 0.01",
        "XY_unavailable_fraction": "must be <= 0.01",
    },
    "shared_IKBKG_input_keys": len(rows),
    "strict_interoperability_retained_keys": len(strict_rows),
    "technical_review_set_keys": len(review_rows),
    "review_reason_pattern_counts": dict(reason_pattern_counts),
    "policy_boundary": (
        "This rule governs sex-aware IKBKG genotype-phased interoperability and later phased workflows only. "
        "It must not be treated as a genotype-only quality exclusion rule."
    ),
    "interpretation_boundary": "No biological, allele-frequency, haplotype or selection interpretation is authorised.",
}
OUTJSON.parent.mkdir(parents=True, exist_ok=True)
OUTJSON.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
print(json.dumps(payload, indent=2))

if len(rows) != 340 or len(strict_rows) != 240 or len(review_rows) != 100:
    raise SystemExit("ERROR: IKBKG technical-set counts do not match the profiled candidate-rule result.")
