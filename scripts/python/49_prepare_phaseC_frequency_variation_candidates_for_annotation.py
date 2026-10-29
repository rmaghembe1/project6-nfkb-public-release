#!/usr/bin/env python3
"""Prepare strict frequency-variation candidates for functional annotation.

This is an annotation-readiness checkpoint only. It prepares candidate input
tables but does not run annotation, assign consequence, prioritise function,
or make biological claims.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path

CAND_DIR = Path("11_results/phaseC/genotype_only_qc/frequency_variation_candidates")
TOP_DELTA = CAND_DIR / "phaseC_strict_top_broad_AF_delta_candidates_v1.tsv"
OBSERVED_ONE = CAND_DIR / "phaseC_strict_observed_only_one_broad_group_candidates_v1.tsv"

OUTDIR = Path("11_results/phaseC/genotype_only_qc/functional_annotation_readiness")
OUTJSON = Path("results/qc/phaseC_functional_annotation_readiness_v1.json")

EXPECTED_TOP_N = 250
EXPECTED_OBSERVED_ONE_N = 250

OUTDIR.mkdir(parents=True, exist_ok=True)
OUTJSON.parent.mkdir(parents=True, exist_ok=True)


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def parse_allele_key(key: str) -> dict[str, str | int]:
    parts = key.split(":")
    if len(parts) != 4:
        raise ValueError(f"Allele key is not CHROM:POS:REF:ALT: {key}")
    chrom, pos, ref, alt = parts
    return {
        "CHROM": chrom,
        "POS": int(pos),
        "REF": ref,
        "ALT": alt,
    }


def bed_coordinates(pos: int, ref: str) -> tuple[int, int]:
    # BED is 0-based half-open. For simple variants, span the REF allele.
    start = pos - 1
    end = start + max(len(ref), 1)
    return start, end


top_rows = read_tsv(TOP_DELTA)
observed_rows = read_tsv(OBSERVED_ONE)

failures: list[str] = []

if len(top_rows) != EXPECTED_TOP_N:
    failures.append(f"Expected {EXPECTED_TOP_N} top AF-delta candidates; observed {len(top_rows)}.")
if len(observed_rows) != EXPECTED_OBSERVED_ONE_N:
    failures.append(f"Expected {EXPECTED_OBSERVED_ONE_N} observed-only candidates; observed {len(observed_rows)}.")

candidate_map: dict[tuple[str, str, str], dict[str, object]] = {}

def add_candidates(rows: list[dict[str, str]], source_label: str) -> None:
    for rank, row in enumerate(rows, start=1):
        key = (row["analysis_layer"], row["gene_symbol"], row["allele_key"])
        parsed = parse_allele_key(row["allele_key"])

        if key not in candidate_map:
            candidate_map[key] = {
                "analysis_layer": row["analysis_layer"],
                "gene_symbol": row["gene_symbol"],
                "allele_key": row["allele_key"],
                "CHROM": parsed["CHROM"],
                "POS": parsed["POS"],
                "REF": parsed["REF"],
                "ALT": parsed["ALT"],
                "candidate_sources": [],
                "source_ranks": [],
                "max_AF_group": row["max_AF_group"],
                "max_AF": row["max_AF"],
                "min_AF_group": row["min_AF_group"],
                "min_AF": row["min_AF"],
                "max_minus_min_AF": row["max_minus_min_AF"],
                "broad_groups_with_AC_gt_0": row["broad_groups_with_AC_gt_0"],
                "observed_only_one_broad_group": row["observed_only_one_broad_group"],
                "groups_with_missing_GT": row["groups_with_missing_GT"],
            }

        candidate_map[key]["candidate_sources"].append(source_label)
        candidate_map[key]["source_ranks"].append(f"{source_label}:{rank}")


add_candidates(top_rows, "top_broad_AF_delta")
add_candidates(observed_rows, "observed_only_one_broad_group")

union_rows = []
for item in candidate_map.values():
    bed_start, bed_end = bed_coordinates(int(item["POS"]), str(item["REF"]))
    union_rows.append({
        "analysis_layer": item["analysis_layer"],
        "gene_symbol": item["gene_symbol"],
        "allele_key": item["allele_key"],
        "CHROM": item["CHROM"],
        "POS": item["POS"],
        "BED_start_0based": bed_start,
        "BED_end_0based_exclusive": bed_end,
        "REF": item["REF"],
        "ALT": item["ALT"],
        "candidate_sources": ";".join(item["candidate_sources"]),
        "source_ranks": ";".join(item["source_ranks"]),
        "max_AF_group": item["max_AF_group"],
        "max_AF": item["max_AF"],
        "min_AF_group": item["min_AF_group"],
        "min_AF": item["min_AF"],
        "max_minus_min_AF": item["max_minus_min_AF"],
        "broad_groups_with_AC_gt_0": item["broad_groups_with_AC_gt_0"],
        "observed_only_one_broad_group": item["observed_only_one_broad_group"],
        "groups_with_missing_GT": item["groups_with_missing_GT"],
    })

union_rows = sorted(
    union_rows,
    key=lambda r: (str(r["CHROM"]), int(r["POS"]), str(r["REF"]), str(r["ALT"]), str(r["gene_symbol"])),
)

if not union_rows:
    failures.append("Union candidate table is empty.")

union_path = OUTDIR / "phaseC_frequency_variation_candidate_union_for_annotation_v1.tsv"
with union_path.open("w", newline="", encoding="utf-8") as handle:
    writer = csv.DictWriter(handle, fieldnames=list(union_rows[0]), delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(union_rows)

bed_path = OUTDIR / "phaseC_frequency_variation_candidate_union_for_annotation_v1.bed"
with bed_path.open("w", encoding="utf-8") as handle:
    for row in union_rows:
        name = f"{row['gene_symbol']}|{row['allele_key']}|{row['candidate_sources']}"
        handle.write(
            f"{row['CHROM']}\t{row['BED_start_0based']}\t{row['BED_end_0based_exclusive']}\t{name}\n"
        )

vcf_like_path = OUTDIR / "phaseC_frequency_variation_candidate_union_for_annotation_v1.vcf.tsv"
with vcf_like_path.open("w", newline="", encoding="utf-8") as handle:
    writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
    writer.writerow(["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"])
    for row in union_rows:
        info = (
            f"GENE={row['gene_symbol']};"
            f"ANALYSIS_LAYER={row['analysis_layer']};"
            f"CANDIDATE_SOURCES={row['candidate_sources']};"
            f"MAX_AF_GROUP={row['max_AF_group']};"
            f"MAX_MINUS_MIN_AF={row['max_minus_min_AF']}"
        )
        writer.writerow([
            row["CHROM"],
            row["POS"],
            row["allele_key"],
            row["REF"],
            row["ALT"],
            ".",
            "PASS",
            info,
        ])

source_summary = {
    "top_broad_AF_delta_input": len(top_rows),
    "observed_only_one_broad_group_input": len(observed_rows),
    "union_candidate_records": len(union_rows),
    "overlap_between_candidate_sources": len(top_rows) + len(observed_rows) - len(union_rows),
    "autosomal_union_records": sum(r["analysis_layer"] == "autosomal" for r in union_rows),
    "chrX_union_records": sum(r["analysis_layer"] == "chrX_sex_aware" for r in union_rows),
}

summary_path = OUTDIR / "phaseC_functional_annotation_readiness_summary_v1.tsv"
with summary_path.open("w", newline="", encoding="utf-8") as handle:
    writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
    writer.writerow(["metric", "value"])
    for key, value in source_summary.items():
        writer.writerow([key, value])

payload = {
    "status": "PASS" if not failures else "FAIL",
    "operation": "Prepare strict frequency-variation candidates for functional annotation input.",
    "candidate_sources": {
        "top_broad_AF_delta": str(TOP_DELTA),
        "observed_only_one_broad_group": str(OBSERVED_ONE),
    },
    "summary": source_summary,
    "outputs": {
        "union_tsv": str(union_path),
        "bed_like": str(bed_path),
        "vcf_like_tsv": str(vcf_like_path),
        "summary": str(summary_path),
    },
    "failures": failures,
    "boundary": (
        "This prepares annotation inputs only. It does not run functional annotation, assign consequence, "
        "prioritise function, infer selection or authorise biological interpretation."
    ),
}

OUTJSON.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
print(json.dumps(payload, indent=2))

if payload["status"] != "PASS":
    raise SystemExit(1)
