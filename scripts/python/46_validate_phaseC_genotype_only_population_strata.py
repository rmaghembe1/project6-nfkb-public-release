#!/usr/bin/env python3
"""Validate population/sample strata before strict genotype-only group-frequency analysis.

This checkpoint validates metadata availability and sample group structure only.
It does not compute group allele frequencies or biological comparisons.
"""

from __future__ import annotations

import csv
import json
from collections import Counter
from pathlib import Path

META = Path("02_reference_data/gnomad_HGDP1kGP_metadata/gnomad_meta_updated.tsv")
COHORT = Path("11_results/phaseB/harmonisation_qc/cohort_primary_genotype_metadata_linked_n4150_v1.txt")
XX = Path("11_results/phaseC/genotype_only_qc/phaseC_IKBKG_candidate_primary_standard_XX_unflagged_samples_v1.txt")
XY = Path("11_results/phaseC/genotype_only_qc/phaseC_IKBKG_candidate_primary_standard_XY_unflagged_samples_v1.txt")

OUTDIR = Path("11_results/phaseC/genotype_only_qc/population_strata_validation")
OUTJSON = Path("results/qc/phaseC_genotype_only_population_strata_validation_v1.json")

EXPECTED_COHORT = 4150
EXPECTED_CHRX_PRIMARY = 4120
MIN_GROUP_SIZE_CANDIDATE = 20

BROAD_GROUP_COLUMN = "hgdp_tgp_meta.Genetic.region"
POPULATION_GROUP_COLUMN = "population"
SENSITIVITY_GROUP_COLUMN = "population_inference.pop"

OUTDIR.mkdir(parents=True, exist_ok=True)
OUTJSON.parent.mkdir(parents=True, exist_ok=True)


def read_lines(path: Path) -> list[str]:
    return [x for x in path.read_text(encoding="utf-8").splitlines() if x]


def normalise(value: str) -> str:
    value = (value or "").strip()
    return value if value else "MISSING"


def write_group_counts(rows_by_sample, ordered_samples, col: str, filename: str):
    counts = Counter(normalise(rows_by_sample[s].get(col, "")) for s in ordered_samples)
    rows = []
    for group, count in sorted(counts.items(), key=lambda x: (-x[1], x[0])):
        rows.append({
            "metadata_column": col,
            "group": group,
            "sample_count": count,
            "meets_candidate_min_n_20": "YES" if count >= MIN_GROUP_SIZE_CANDIDATE else "NO",
        })
    with (OUTDIR / filename).open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    return rows


cohort = read_lines(COHORT)
xx = read_lines(XX)
xy = read_lines(XY)
chrx_primary = set(xx) | set(xy)
chrx_ordered = [s for s in cohort if s in chrx_primary]

failures = []

if len(cohort) != EXPECTED_COHORT or len(set(cohort)) != EXPECTED_COHORT:
    failures.append(f"Expected {EXPECTED_COHORT} unique genotype-only cohort samples; observed {len(cohort)} lines and {len(set(cohort))} unique IDs.")
if len(xx) + len(xy) != EXPECTED_CHRX_PRIMARY:
    failures.append(f"Expected {EXPECTED_CHRX_PRIMARY} primary chrX samples; observed {len(xx) + len(xy)}.")
if len(chrx_ordered) != EXPECTED_CHRX_PRIMARY:
    failures.append(f"Expected {EXPECTED_CHRX_PRIMARY} chrX ordered samples; observed {len(chrx_ordered)}.")

with META.open(newline="", encoding="utf-8") as handle:
    reader = csv.DictReader(handle, delimiter="\t")
    fieldnames = reader.fieldnames or []
    rows_by_sample = {row["s"]: row for row in reader if row.get("s") in set(cohort)}

if len(rows_by_sample) != EXPECTED_COHORT:
    failures.append(f"Expected metadata for {EXPECTED_COHORT} cohort samples; observed {len(rows_by_sample)}.")

required_columns = [
    BROAD_GROUP_COLUMN,
    POPULATION_GROUP_COLUMN,
    SENSITIVITY_GROUP_COLUMN,
    "sex_imputation.sex_karyotype",
    "sample_filters.sex_aneuploidy",
    "sample_filters.ambiguous_sex",
    "sex",
]
for col in required_columns:
    if col not in fieldnames:
        failures.append(f"Required metadata column missing: {col}")

available_population_related_columns = [
    c for c in fieldnames
    if any(token in c.lower() for token in ["population", "pop", "ancestry", "region", "project", "country", "sex"])
]

column_summary_rows = []
for col in available_population_related_columns:
    values = Counter(normalise(rows_by_sample[s].get(col, "")) for s in cohort)
    non_missing = sum(v for k, v in values.items() if k != "MISSING")
    nonmissing_sizes = [v for k, v in values.items() if k != "MISSING"]
    column_summary_rows.append({
        "metadata_column": col,
        "distinct_values": len(values),
        "non_missing_samples": non_missing,
        "missing_samples": values.get("MISSING", 0),
        "largest_group_size": max(values.values()) if values else 0,
        "smallest_nonmissing_group_size": min(nonmissing_sizes) if nonmissing_sizes else 0,
    })

with (OUTDIR / "phaseC_population_metadata_column_availability_v1.tsv").open("w", newline="", encoding="utf-8") as handle:
    writer = csv.DictWriter(handle, fieldnames=list(column_summary_rows[0]), delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(column_summary_rows)

broad_rows = write_group_counts(rows_by_sample, cohort, BROAD_GROUP_COLUMN, "phaseC_genotype_only_broad_group_counts_v1.tsv")
pop_rows = write_group_counts(rows_by_sample, cohort, POPULATION_GROUP_COLUMN, "phaseC_genotype_only_population_group_counts_v1.tsv")
sens_rows = write_group_counts(rows_by_sample, cohort, SENSITIVITY_GROUP_COLUMN, "phaseC_genotype_only_population_inference_group_counts_v1.tsv")

chrx_broad_rows = write_group_counts(rows_by_sample, chrx_ordered, BROAD_GROUP_COLUMN, "phaseC_chrX_primary_broad_group_counts_v1.tsv")
chrx_pop_rows = write_group_counts(rows_by_sample, chrx_ordered, POPULATION_GROUP_COLUMN, "phaseC_chrX_primary_population_group_counts_v1.tsv")
chrx_sens_rows = write_group_counts(rows_by_sample, chrx_ordered, SENSITIVITY_GROUP_COLUMN, "phaseC_chrX_primary_population_inference_group_counts_v1.tsv")

with (OUTDIR / "phaseC_genotype_only_population_strata_sample_order_v1.txt").open("w", encoding="utf-8") as handle:
    handle.write("\n".join(cohort) + "\n")

with (OUTDIR / "phaseC_chrX_primary_population_strata_sample_order_v1.txt").open("w", encoding="utf-8") as handle:
    handle.write("\n".join(chrx_ordered) + "\n")

payload = {
    "status": "PASS" if not failures else "FAIL",
    "operation": "Validate metadata population strata before strict genotype-only group-frequency analysis.",
    "sample_counts": {
        "genotype_only_cohort": len(cohort),
        "chrX_primary_cohort": len(chrx_ordered),
        "XX_chrX_primary": len(xx),
        "XY_chrX_primary": len(xy),
    },
    "selected_grouping_columns": {
        "broad_group": BROAD_GROUP_COLUMN,
        "population_group": POPULATION_GROUP_COLUMN,
        "sensitivity_inferred_population_group": SENSITIVITY_GROUP_COLUMN,
    },
    "candidate_min_group_size": MIN_GROUP_SIZE_CANDIDATE,
    "broad_group_count": len(broad_rows),
    "population_group_count": len(pop_rows),
    "sensitivity_inferred_population_group_count": len(sens_rows),
    "chrX_broad_group_count": len(chrx_broad_rows),
    "chrX_population_group_count": len(chrx_pop_rows),
    "chrX_sensitivity_inferred_population_group_count": len(chrx_sens_rows),
    "broad_groups_below_min_n_20": sum(row["meets_candidate_min_n_20"] == "NO" for row in broad_rows),
    "population_groups_below_min_n_20": sum(row["meets_candidate_min_n_20"] == "NO" for row in pop_rows),
    "chrX_population_groups_below_min_n_20": sum(row["meets_candidate_min_n_20"] == "NO" for row in chrx_pop_rows),
    "failures": failures,
    "boundary": "This validates metadata strata only. No group-level allele frequencies, population comparisons, ancestry claims or biological interpretation are authorised.",
}

OUTJSON.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
print(json.dumps(payload, indent=2))

if payload["status"] != "PASS":
    raise SystemExit(1)
