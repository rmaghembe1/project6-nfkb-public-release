#!/usr/bin/env python3
"""Profile autosomal genotype-phased GT discordance and missingness before defining filtering rules."""

from __future__ import annotations

import csv
import json
from collections import Counter, defaultdict
from pathlib import Path

INPUT = Path(
    "11_results/phaseB/genotype_call_concordance_qc/autosomal_full_panel/"
    "autosomal_direct_GT_concordance_by_shared_allele_key_v1.tsv"
)
OUTDIR = Path(
    "11_results/phaseB/genotype_call_concordance_qc/autosomal_filter_rule_development"
)
OUTJSON = Path("results/qc/phaseB_autosomal_concordance_filter_rule_profile_v1.json")

SAMPLES = 4091


def variant_class(allele_key: str) -> str:
    _, _, ref, alt = allele_key.split(":", 3)
    if len(ref) == 1 and len(alt) == 1:
        return "SNV"
    if len(ref) != len(alt):
        return "INDEL"
    return "MNV_OR_COMPLEX"


with INPUT.open(newline="", encoding="utf-8") as handle:
    rows = list(csv.DictReader(handle, delimiter="\t"))

profiles = []
for row in rows:
    compared = int(row["nonmissing_GT_comparisons"])
    discordant = int(row["discordant_GT_calls"])
    missing = int(row["missing_or_nondiploid_in_either_resource"])
    profiles.append({
        "gene_symbol": row["gene_symbol"],
        "allele_key": row["allele_key"],
        "variant_class": variant_class(row["allele_key"]),
        "nonmissing_GT_comparisons": compared,
        "discordant_GT_calls": discordant,
        "missing_or_nondiploid_in_either_resource": missing,
        "discordance_fraction_among_nonmissing": discordant / compared if compared else "",
        "missingness_fraction_of_expected_cells": missing / SAMPLES,
    })

OUTDIR.mkdir(parents=True, exist_ok=True)

profile_tsv = OUTDIR / "autosomal_shared_allele_key_discordance_missingness_profile_v1.tsv"
with profile_tsv.open("w", newline="", encoding="utf-8") as handle:
    fields = list(profiles[0].keys())
    writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(profiles)

class_counts = defaultdict(lambda: Counter())
for row in profiles:
    klass = row["variant_class"]
    class_counts[klass]["allele_keys"] += 1
    class_counts[klass]["keys_with_any_discordance"] += int(row["discordant_GT_calls"] > 0)
    class_counts[klass]["discordant_GT_calls"] += row["discordant_GT_calls"]
    class_counts[klass]["nonmissing_GT_comparisons"] += row["nonmissing_GT_comparisons"]
    class_counts[klass]["missing_cells"] += row["missing_or_nondiploid_in_either_resource"]

class_tsv = OUTDIR / "autosomal_concordance_summary_by_variant_class_v1.tsv"
with class_tsv.open("w", newline="", encoding="utf-8") as handle:
    fields = [
        "variant_class", "allele_keys", "keys_with_any_discordance",
        "discordant_GT_calls", "nonmissing_GT_comparisons",
        "discordance_fraction_among_nonmissing", "missing_cells"
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
            "nonmissing_GT_comparisons": x["nonmissing_GT_comparisons"],
            "discordance_fraction_among_nonmissing": (
                x["discordant_GT_calls"] / x["nonmissing_GT_comparisons"]
                if x["nonmissing_GT_comparisons"] else ""
            ),
            "missing_cells": x["missing_cells"],
        })

threshold_rows = []
thresholds = [
    ("any_discordance", lambda r: r["discordant_GT_calls"] > 0),
    ("discordance_gt_0.001", lambda r: r["discordance_fraction_among_nonmissing"] != "" and r["discordance_fraction_among_nonmissing"] > 0.001),
    ("discordance_gt_0.01", lambda r: r["discordance_fraction_among_nonmissing"] != "" and r["discordance_fraction_among_nonmissing"] > 0.01),
    ("missingness_gt_0.01", lambda r: r["missingness_fraction_of_expected_cells"] > 0.01),
    ("missingness_gt_0.05", lambda r: r["missingness_fraction_of_expected_cells"] > 0.05),
    ("discordance_any_or_missingness_gt_0.01", lambda r: r["discordant_GT_calls"] > 0 or r["missingness_fraction_of_expected_cells"] > 0.01),
]
for label, rule in thresholds:
    selected = [row for row in profiles if rule(row)]
    threshold_rows.append({
        "candidate_rule_component": label,
        "flagged_allele_keys": len(selected),
        "retained_allele_keys_if_excluded": len(profiles) - len(selected),
    })

threshold_tsv = OUTDIR / "autosomal_candidate_discordance_missingness_threshold_counts_v1.tsv"
with threshold_tsv.open("w", newline="", encoding="utf-8") as handle:
    fields = list(threshold_rows[0].keys())
    writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(threshold_rows)

top = sorted(
    profiles,
    key=lambda r: (r["discordant_GT_calls"], r["missing_or_nondiploid_in_either_resource"]),
    reverse=True,
)[:50]
top_tsv = OUTDIR / "autosomal_top50_discordant_shared_allele_keys_v1.tsv"
with top_tsv.open("w", newline="", encoding="utf-8") as handle:
    fields = list(top[0].keys())
    writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(top)

payload = {
    "status": "PASS",
    "operation": "Technical profiling of autosomal shared-key GT discordance and missingness before filtering-rule definition.",
    "shared_autosomal_allele_keys_profiled": len(profiles),
    "allele_keys_with_any_discordance": sum(row["discordant_GT_calls"] > 0 for row in profiles),
    "total_discordant_GT_calls": sum(row["discordant_GT_calls"] for row in profiles),
    "total_missing_or_nondiploid_cells": sum(row["missing_or_nondiploid_in_either_resource"] for row in profiles),
    "variant_classes_present": sorted(class_counts.keys()),
    "interpretation_boundary": (
        "This profile informs technical filtering-rule development only; "
        "no downstream analysis exclusions are yet authorised."
    ),
}
OUTJSON.parent.mkdir(parents=True, exist_ok=True)
OUTJSON.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
print(json.dumps(payload, indent=2))
