#!/usr/bin/env python3
"""Profile IKBKG chrX sex-aware GT discordance and unavailable cells before rule definition."""

from __future__ import annotations

import csv
import json
from collections import Counter, defaultdict
from pathlib import Path

INPUT = Path(
    "11_results/phaseB/genotype_call_concordance_qc/chrX_IKBKG_sex_aware/"
    "IKBKG_chrX_sex_aware_direct_GT_concordance_by_shared_allele_key_v1.tsv"
)
OUTDIR = Path(
    "11_results/phaseB/genotype_call_concordance_qc/chrX_IKBKG_filter_rule_development"
)
OUTJSON = Path("results/qc/phaseB_IKBKG_chrX_filter_rule_profile_v1.json")

EXPECTED = {"XX": 1896, "XY": 2195}
TOTAL_EXPECTED = sum(EXPECTED.values())


def variant_class(allele_key: str) -> str:
    _, _, ref, alt = allele_key.split(":", 3)
    if len(ref) == 1 and len(alt) == 1:
        return "SNV"
    if len(ref) != len(alt):
        return "INDEL"
    return "MNV_OR_COMPLEX"


with INPUT.open(newline="", encoding="utf-8") as handle:
    rows = list(csv.DictReader(handle, delimiter="\t"))

by_key: dict[str, dict[str, dict[str, int]]] = defaultdict(dict)
for row in rows:
    stratum = row["sex_karyotype"]
    key = row["allele_key"]
    by_key[key][stratum] = {
        "compared": int(row["nonmissing_comparisons"]),
        "matches": int(row["matching_GT_calls"]),
        "discordant": int(row["discordant_GT_calls"]),
        "unavailable": int(row["unavailable_or_unexpected_cells"]),
    }

failures = []
if len(by_key) != 340:
    failures.append(f"Expected 340 shared IKBKG keys, observed {len(by_key)}.")

profiles = []
for key in sorted(by_key):
    strata = by_key[key]
    if set(strata) != {"XX", "XY"}:
        failures.append(f"{key}: does not contain both XX and XY strata.")
        continue

    xx = strata["XX"]
    xy = strata["XY"]
    compared = xx["compared"] + xy["compared"]
    matches = xx["matches"] + xy["matches"]
    discordant = xx["discordant"] + xy["discordant"]
    unavailable = xx["unavailable"] + xy["unavailable"]

    profiles.append({
        "allele_key": key,
        "variant_class": variant_class(key),
        "XX_nonmissing_comparisons": xx["compared"],
        "XX_matching_GT_calls": xx["matches"],
        "XX_discordant_GT_calls": xx["discordant"],
        "XX_unavailable_cells": xx["unavailable"],
        "XX_unavailable_fraction": xx["unavailable"] / EXPECTED["XX"],
        "XY_nonmissing_comparisons": xy["compared"],
        "XY_matching_GT_calls": xy["matches"],
        "XY_discordant_GT_calls": xy["discordant"],
        "XY_unavailable_cells": xy["unavailable"],
        "XY_unavailable_fraction": xy["unavailable"] / EXPECTED["XY"],
        "combined_nonmissing_comparisons": compared,
        "combined_matching_GT_calls": matches,
        "combined_discordant_GT_calls": discordant,
        "combined_unavailable_cells": unavailable,
        "combined_unavailable_fraction": unavailable / TOTAL_EXPECTED,
        "combined_concordance_fraction": matches / compared if compared else "",
    })

OUTDIR.mkdir(parents=True, exist_ok=True)

profile_path = OUTDIR / "IKBKG_chrX_shared_allele_key_discordance_unavailable_profile_v1.tsv"
with profile_path.open("w", newline="", encoding="utf-8") as handle:
    writer = csv.DictWriter(handle, fieldnames=list(profiles[0]), delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(profiles)

class_counts = defaultdict(Counter)
for row in profiles:
    klass = row["variant_class"]
    class_counts[klass]["allele_keys"] += 1
    class_counts[klass]["keys_with_any_discordance"] += int(row["combined_discordant_GT_calls"] > 0)
    class_counts[klass]["discordant_GT_calls"] += row["combined_discordant_GT_calls"]
    class_counts[klass]["nonmissing_comparisons"] += row["combined_nonmissing_comparisons"]
    class_counts[klass]["unavailable_cells"] += row["combined_unavailable_cells"]

class_path = OUTDIR / "IKBKG_chrX_concordance_summary_by_variant_class_v1.tsv"
with class_path.open("w", newline="", encoding="utf-8") as handle:
    fields = [
        "variant_class", "allele_keys", "keys_with_any_discordance",
        "discordant_GT_calls", "nonmissing_comparisons",
        "discordance_fraction_among_nonmissing", "unavailable_cells"
    ]
    writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    for klass in sorted(class_counts):
        x = class_counts[klass]
        writer.writerow({
            "variant_class": klass,
            "allele_keys": x["allele_keys"],
            "keys_with_any_discordance": x["keys_with_any_discordance"],
            "discordant_GT_calls": x["discordant_GT_calls"],
            "nonmissing_comparisons": x["nonmissing_comparisons"],
            "discordance_fraction_among_nonmissing": (
                x["discordant_GT_calls"] / x["nonmissing_comparisons"]
                if x["nonmissing_comparisons"] else ""
            ),
            "unavailable_cells": x["unavailable_cells"],
        })

rules = [
    ("any_discordance", lambda r: r["combined_discordant_GT_calls"] > 0),
    ("combined_unavailable_gt_0.01", lambda r: r["combined_unavailable_fraction"] > 0.01),
    ("XX_unavailable_gt_0.01", lambda r: r["XX_unavailable_fraction"] > 0.01),
    ("XY_unavailable_gt_0.01", lambda r: r["XY_unavailable_fraction"] > 0.01),
    ("either_stratum_unavailable_gt_0.01", lambda r: r["XX_unavailable_fraction"] > 0.01 or r["XY_unavailable_fraction"] > 0.01),
    ("discordance_any_or_combined_unavailable_gt_0.01", lambda r: r["combined_discordant_GT_calls"] > 0 or r["combined_unavailable_fraction"] > 0.01),
    ("discordance_any_or_either_stratum_unavailable_gt_0.01", lambda r: r["combined_discordant_GT_calls"] > 0 or r["XX_unavailable_fraction"] > 0.01 or r["XY_unavailable_fraction"] > 0.01),
]

threshold_rows = []
for label, rule in rules:
    flagged = [row for row in profiles if rule(row)]
    threshold_rows.append({
        "candidate_rule_component": label,
        "flagged_allele_keys": len(flagged),
        "retained_allele_keys_if_excluded": len(profiles) - len(flagged),
    })

threshold_path = OUTDIR / "IKBKG_chrX_candidate_discordance_unavailable_threshold_counts_v1.tsv"
with threshold_path.open("w", newline="", encoding="utf-8") as handle:
    writer = csv.DictWriter(handle, fieldnames=list(threshold_rows[0]), delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(threshold_rows)

top = sorted(
    profiles,
    key=lambda r: (r["combined_discordant_GT_calls"], r["combined_unavailable_cells"]),
    reverse=True,
)[:40]

top_path = OUTDIR / "IKBKG_chrX_top40_discordant_shared_allele_keys_v1.tsv"
with top_path.open("w", newline="", encoding="utf-8") as handle:
    writer = csv.DictWriter(handle, fieldnames=list(top[0]), delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(top)

payload = {
    "status": "PASS" if not failures else "FAIL",
    "operation": "Technical profiling of sex-aware IKBKG chrX concordance and unavailable cells before rule definition.",
    "shared_IKBKG_allele_keys_profiled": len(profiles),
    "allele_keys_with_any_discordance": sum(r["combined_discordant_GT_calls"] > 0 for r in profiles),
    "total_discordant_GT_calls": sum(r["combined_discordant_GT_calls"] for r in profiles),
    "total_unavailable_cells": sum(r["combined_unavailable_cells"] for r in profiles),
    "strata": EXPECTED,
    "failures": failures,
    "interpretation_boundary": (
        "This profile informs IKBKG-specific technical-rule development only; "
        "no chrX exclusion rule or biological interpretation is authorised."
    ),
}
OUTJSON.parent.mkdir(parents=True, exist_ok=True)
OUTJSON.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
print(json.dumps(payload, indent=2))

if payload["status"] != "PASS":
    raise SystemExit(1)
