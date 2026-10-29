#!/usr/bin/env python3
"""Compute strict genotype-only group-frequency QC tables.

This is a technical preparation/QC step only. It computes group-level AC/AN/AF
from strict genotype-only GT matrices using validated population strata.
It does not authorise population comparison or biological interpretation.
"""

from __future__ import annotations

import csv
import json
from collections import Counter, defaultdict
from pathlib import Path

META = Path("02_reference_data/gnomad_HGDP1kGP_metadata/gnomad_meta_updated.tsv")

MATRIX_DIR = Path("11_results/phaseC/genotype_only_qc/strict_analysis_matrices")
AUTO_MATRIX = MATRIX_DIR / "phaseC_strict_autosomal_genotype_GT_matrix_v1.tsv"
CHRX_XX_MATRIX = MATRIX_DIR / "phaseC_strict_IKBKG_chrX_XX_genotype_GT_matrix_v1.tsv"
CHRX_XY_MATRIX = MATRIX_DIR / "phaseC_strict_IKBKG_chrX_XY_genotype_GT_matrix_v1.tsv"

STRATA_DIR = Path("11_results/phaseC/genotype_only_qc/population_strata_validation")
AUTO_ORDER = STRATA_DIR / "phaseC_genotype_only_population_strata_sample_order_v1.txt"
CHRX_ORDER = STRATA_DIR / "phaseC_chrX_primary_population_strata_sample_order_v1.txt"

OUTDIR = Path("11_results/phaseC/genotype_only_qc/strict_group_frequency_qc")
OUTJSON = Path("results/qc/phaseC_strict_genotype_only_group_frequency_qc_v1.json")

BROAD_COL = "hgdp_tgp_meta.Genetic.region"
POP_COL = "population"
INFERRED_COL = "population_inference.pop"
MIN_N = 20

EXPECTED_AUTO_RECORDS = 39438
EXPECTED_CHRX_RECORDS = 383
EXPECTED_AUTO_SAMPLES = 4150
EXPECTED_CHRX_SAMPLES = 4120
EXPECTED_XX_SAMPLES = 1904
EXPECTED_XY_SAMPLES = 2216

OUTDIR.mkdir(parents=True, exist_ok=True)
OUTJSON.parent.mkdir(parents=True, exist_ok=True)


def read_lines(path: Path) -> list[str]:
    return [x for x in path.read_text(encoding="utf-8").splitlines() if x]


def normalise(value: str) -> str:
    value = (value or "").strip()
    return value if value else "MISSING"


def parse_auto_gt(gt: str) -> tuple[int, int, int]:
    if "." in gt:
        return 0, 0, 1
    if gt not in {"0/0", "0/1", "1/0", "1/1"}:
        raise ValueError(f"Unexpected autosomal GT: {gt}")
    a, b = gt.split("/")
    return 2, int(a) + int(b), 0


def parse_xx_gt(gt: str) -> tuple[int, int, int]:
    if "." in gt:
        return 0, 0, 1
    if "/" not in gt:
        raise ValueError(f"Unexpected XX chrX GT: {gt}")
    a, b = gt.split("/")
    return 2, int(a) + int(b), 0


def parse_xy_gt(gt: str) -> tuple[int, int, int]:
    if "." in gt:
        return 0, 0, 1
    if "/" in gt or "|" in gt:
        raise ValueError(f"Unexpected XY chrX GT: {gt}")
    return 1, int(gt), 0


def load_metadata(cohort: set[str]) -> dict[str, dict[str, str]]:
    with META.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return {row["s"]: row for row in reader if row.get("s") in cohort}


def build_group_index(samples: list[str], meta: dict[str, dict[str, str]], col: str) -> dict[str, list[int]]:
    groups: dict[str, list[int]] = defaultdict(list)
    for i, sample in enumerate(samples):
        groups[normalise(meta[sample].get(col, ""))].append(i)
    return dict(groups)


def group_sizes_from_index(groups: dict[str, list[int]]) -> dict[str, int]:
    return {group: len(indices) for group, indices in groups.items()}


auto_samples = read_lines(AUTO_ORDER)
chrx_samples = read_lines(CHRX_ORDER)

failures: list[str] = []

if len(auto_samples) != EXPECTED_AUTO_SAMPLES:
    failures.append(f"Expected {EXPECTED_AUTO_SAMPLES} autosomal samples; observed {len(auto_samples)}.")
if len(chrx_samples) != EXPECTED_CHRX_SAMPLES:
    failures.append(f"Expected {EXPECTED_CHRX_SAMPLES} chrX primary samples; observed {len(chrx_samples)}.")

meta = load_metadata(set(auto_samples))
if len(meta) != EXPECTED_AUTO_SAMPLES:
    failures.append(f"Expected metadata for {EXPECTED_AUTO_SAMPLES} samples; observed {len(meta)}.")

# chrX ordered samples are XX followed by XY? No: validation order follows cohort order.
# For dosage, the matrices remain separate and each has its own sample order from the original sample lists.
xx_samples = read_lines(Path("11_results/phaseC/genotype_only_qc/phaseC_IKBKG_candidate_primary_standard_XX_unflagged_samples_v1.txt"))
xy_samples = read_lines(Path("11_results/phaseC/genotype_only_qc/phaseC_IKBKG_candidate_primary_standard_XY_unflagged_samples_v1.txt"))

if len(xx_samples) != EXPECTED_XX_SAMPLES:
    failures.append(f"Expected {EXPECTED_XX_SAMPLES} XX samples; observed {len(xx_samples)}.")
if len(xy_samples) != EXPECTED_XY_SAMPLES:
    failures.append(f"Expected {EXPECTED_XY_SAMPLES} XY samples; observed {len(xy_samples)}.")

grouping_specs = [
    ("broad", BROAD_COL, "primary_broad"),
    ("population", POP_COL, "primary_population"),
    ("inferred_population", INFERRED_COL, "sensitivity_inferred_population"),
]

auto_group_indices = {
    label: build_group_index(auto_samples, meta, col)
    for label, col, _purpose in grouping_specs
}
xx_group_indices = {
    label: build_group_index(xx_samples, meta, col)
    for label, col, _purpose in grouping_specs
}
xy_group_indices = {
    label: build_group_index(xy_samples, meta, col)
    for label, col, _purpose in grouping_specs
}

def write_group_size_table() -> None:
    rows = []
    for label, col, purpose in grouping_specs:
        for cohort_label, group_index in [
            ("autosomal", auto_group_indices[label]),
            ("IKBKG_chrX_XX", xx_group_indices[label]),
            ("IKBKG_chrX_XY", xy_group_indices[label]),
        ]:
            for group, indices in sorted(group_index.items(), key=lambda x: (-len(x[1]), x[0])):
                rows.append({
                    "grouping_layer": label,
                    "grouping_purpose": purpose,
                    "metadata_column": col,
                    "cohort_layer": cohort_label,
                    "group": group,
                    "sample_count": len(indices),
                    "meets_min_n_20": "YES" if len(indices) >= MIN_N else "NO",
                })
    with (OUTDIR / "phaseC_strict_group_frequency_sample_counts_v1.tsv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


write_group_size_table()


def process_autosomal_grouping(label: str, col: str, purpose: str, group_index: dict[str, list[int]]) -> tuple[Path, dict[str, int]]:
    output = OUTDIR / f"phaseC_strict_autosomal_{label}_group_frequency_qc_v1.tsv"
    summary = Counter()
    with AUTO_MATRIX.open(newline="", encoding="utf-8") as in_handle, output.open("w", newline="", encoding="utf-8") as out_handle:
        reader = csv.DictReader(in_handle, delimiter="\t")
        fields = [
            "analysis_layer", "grouping_layer", "grouping_purpose", "metadata_column",
            "group", "group_sample_count", "meets_min_n_20",
            "gene_symbol", "allele_key", "missing_GT_cells", "AN", "AC", "AF",
        ]
        writer = csv.DictWriter(out_handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()

        record_count = 0
        for row in reader:
            record_count += 1
            gts = row["GT_values_in_frozen_4150_sample_order"].split(",")
            if len(gts) != EXPECTED_AUTO_SAMPLES:
                failures.append(f"{row['gene_symbol']} {row['allele_key']}: autosomal GT length mismatch.")
                continue

            for group, indices in group_index.items():
                an = ac = missing = 0
                for i in indices:
                    called, alt, miss = parse_auto_gt(gts[i])
                    an += called
                    ac += alt
                    missing += miss
                writer.writerow({
                    "analysis_layer": "autosomal",
                    "grouping_layer": label,
                    "grouping_purpose": purpose,
                    "metadata_column": col,
                    "group": group,
                    "group_sample_count": len(indices),
                    "meets_min_n_20": "YES" if len(indices) >= MIN_N else "NO",
                    "gene_symbol": row["gene_symbol"],
                    "allele_key": row["allele_key"],
                    "missing_GT_cells": missing,
                    "AN": an,
                    "AC": ac,
                    "AF": ac / an if an else "",
                })
            summary["input_records"] = record_count

    summary["output_rows"] = summary["input_records"] * len(group_index)
    summary["groups"] = len(group_index)
    summary["groups_below_min_n_20"] = sum(len(v) < MIN_N for v in group_index.values())
    return output, dict(summary)


def process_chrx_grouping(label: str, col: str, purpose: str, xx_index: dict[str, list[int]], xy_index: dict[str, list[int]]) -> tuple[Path, dict[str, int]]:
    output = OUTDIR / f"phaseC_strict_IKBKG_chrX_{label}_group_frequency_qc_v1.tsv"

    with CHRX_XX_MATRIX.open(newline="", encoding="utf-8") as handle:
        xx_rows = {row["allele_key"]: row for row in csv.DictReader(handle, delimiter="\t")}
    with CHRX_XY_MATRIX.open(newline="", encoding="utf-8") as handle:
        xy_rows = {row["allele_key"]: row for row in csv.DictReader(handle, delimiter="\t")}

    if set(xx_rows) != set(xy_rows):
        failures.append(f"chrX XX/XY key mismatch for grouping {label}.")

    all_groups = sorted(set(xx_index) | set(xy_index))
    with output.open("w", newline="", encoding="utf-8") as out_handle:
        fields = [
            "analysis_layer", "grouping_layer", "grouping_purpose", "metadata_column",
            "group", "XX_sample_count", "XY_sample_count", "group_sample_count",
            "meets_min_n_20", "gene_symbol", "allele_key",
            "XX_missing_GT_cells", "XY_missing_GT_cells", "missing_GT_cells",
            "XX_AN", "XY_AN", "AN", "XX_AC", "XY_AC", "AC", "AF",
        ]
        writer = csv.DictWriter(out_handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()

        record_count = 0
        for key in sorted(xx_rows):
            record_count += 1
            xx_gts = xx_rows[key]["GT_values_in_frozen_XX_sample_order"].split(",")
            xy_gts = xy_rows[key]["GT_values_in_frozen_XY_sample_order"].split(",")

            if len(xx_gts) != EXPECTED_XX_SAMPLES:
                failures.append(f"IKBKG {key}: XX GT length mismatch.")
                continue
            if len(xy_gts) != EXPECTED_XY_SAMPLES:
                failures.append(f"IKBKG {key}: XY GT length mismatch.")
                continue

            for group in all_groups:
                xx_indices = xx_index.get(group, [])
                xy_indices = xy_index.get(group, [])

                xx_an = xx_ac = xx_missing = 0
                xy_an = xy_ac = xy_missing = 0

                for i in xx_indices:
                    called, alt, miss = parse_xx_gt(xx_gts[i])
                    xx_an += called
                    xx_ac += alt
                    xx_missing += miss

                for i in xy_indices:
                    called, alt, miss = parse_xy_gt(xy_gts[i])
                    xy_an += called
                    xy_ac += alt
                    xy_missing += miss

                an = xx_an + xy_an
                ac = xx_ac + xy_ac
                n = len(xx_indices) + len(xy_indices)

                writer.writerow({
                    "analysis_layer": "chrX_sex_aware",
                    "grouping_layer": label,
                    "grouping_purpose": purpose,
                    "metadata_column": col,
                    "group": group,
                    "XX_sample_count": len(xx_indices),
                    "XY_sample_count": len(xy_indices),
                    "group_sample_count": n,
                    "meets_min_n_20": "YES" if n >= MIN_N else "NO",
                    "gene_symbol": "IKBKG",
                    "allele_key": key,
                    "XX_missing_GT_cells": xx_missing,
                    "XY_missing_GT_cells": xy_missing,
                    "missing_GT_cells": xx_missing + xy_missing,
                    "XX_AN": xx_an,
                    "XY_AN": xy_an,
                    "AN": an,
                    "XX_AC": xx_ac,
                    "XY_AC": xy_ac,
                    "AC": ac,
                    "AF": ac / an if an else "",
                })

    summary = {
        "input_records": record_count,
        "groups": len(all_groups),
        "output_rows": record_count * len(all_groups),
        "groups_below_min_n_20": sum(
            (len(xx_index.get(g, [])) + len(xy_index.get(g, []))) < MIN_N for g in all_groups
        ),
    }
    return output, summary


output_summaries = []

for label, col, purpose in grouping_specs:
    auto_path, auto_summary = process_autosomal_grouping(label, col, purpose, auto_group_indices[label])
    chrx_path, chrx_summary = process_chrx_grouping(label, col, purpose, xx_group_indices[label], xy_group_indices[label])

    output_summaries.append({
        "output_layer": f"autosomal_{label}",
        "path": str(auto_path),
        **auto_summary,
    })
    output_summaries.append({
        "output_layer": f"IKBKG_chrX_{label}",
        "path": str(chrx_path),
        **chrx_summary,
    })

summary_path = OUTDIR / "phaseC_strict_group_frequency_qc_output_summary_v1.tsv"
with summary_path.open("w", newline="", encoding="utf-8") as handle:
    fields = ["output_layer", "path", "input_records", "groups", "output_rows", "groups_below_min_n_20"]
    writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(output_summaries)

expected_rows = {
    "autosomal_broad": EXPECTED_AUTO_RECORDS * 7,
    "IKBKG_chrX_broad": EXPECTED_CHRX_RECORDS * 7,
    "autosomal_population": EXPECTED_AUTO_RECORDS * 80,
    "IKBKG_chrX_population": EXPECTED_CHRX_RECORDS * 80,
    "autosomal_inferred_population": EXPECTED_AUTO_RECORDS * 9,
    "IKBKG_chrX_inferred_population": EXPECTED_CHRX_RECORDS * 9,
}

for item in output_summaries:
    layer = item["output_layer"]
    if layer in expected_rows and int(item["output_rows"]) != expected_rows[layer]:
        failures.append(f"{layer}: expected {expected_rows[layer]} rows; observed {item['output_rows']}.")

payload = {
    "status": "PASS" if not failures else "FAIL",
    "operation": "Compute strict genotype-only group-frequency QC tables.",
    "grouping_columns": {
        "broad": BROAD_COL,
        "population": POP_COL,
        "inferred_population_sensitivity": INFERRED_COL,
    },
    "minimum_group_size_flag": MIN_N,
    "output_summaries": output_summaries,
    "failures": failures,
    "boundary": (
        "Group-level AC/AN/AF are computed for technical preparation only. "
        "No population comparison, ancestry claim, selection inference, functional prioritisation "
        "or biological interpretation is authorised."
    ),
}

OUTJSON.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
print(json.dumps(payload, indent=2))

if payload["status"] != "PASS":
    raise SystemExit(1)
