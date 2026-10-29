#!/usr/bin/env python3
"""Generate chromosome inventory for later interval-restricted genotype/haplotype querying."""

from __future__ import annotations

import argparse
import csv
import json
from collections import defaultdict
from pathlib import Path


def chromosome_key(chrom: str) -> int:
    label = chrom.removeprefix("chr")
    return 23 if label == "X" else 24 if label == "Y" else int(label)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--coordinates", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--qc", required=True, type=Path)
    args = parser.parse_args()

    by_chromosome: dict[str, list[str]] = defaultdict(list)

    with args.coordinates.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))

    for row in rows:
        by_chromosome[row["chromosome"]].append(row["gene_symbol"])

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.qc.parent.mkdir(parents=True, exist_ok=True)

    with args.out.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow([
            "chromosome", "primary_gene_count", "primary_genes",
            "genotype_vcf_required", "phased_bcf_required_later",
            "extraction_status"
        ])
        for chrom in sorted(by_chromosome, key=chromosome_key):
            genes = sorted(by_chromosome[chrom])
            writer.writerow([
                chrom,
                len(genes),
                ",".join(genes),
                "yes",
                "yes_after_sample_compatibility_check",
                "not_started_coordinate_freeze_only",
            ])

    payload = {
        "status": "PASS" if len(rows) == 24 and len(by_chromosome) == 14 else "FAIL",
        "primary_gene_count": len(rows),
        "required_chromosome_count": len(by_chromosome),
        "required_chromosomes": sorted(by_chromosome, key=chromosome_key),
        "policy": "Inventory only; no genotype or haplotype extraction initiated.",
    }
    args.qc.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(payload, indent=2))

    if payload["status"] != "PASS":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
