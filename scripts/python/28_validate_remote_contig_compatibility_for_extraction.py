#!/usr/bin/env python3
"""Validate frozen BED chromosome names against explicit Phase B remote-resource routes."""

from __future__ import annotations

import argparse
import csv
import json
import subprocess
from pathlib import Path


def read_bed_chromosomes(path: Path) -> list[str]:
    chromosomes = []
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            chrom = line.split("\t", 1)[0]
            if chrom not in chromosomes:
                chromosomes.append(chrom)
    return chromosomes


def read_remote_header(url: str) -> tuple[set[str], str]:
    result = subprocess.run(
        ["bcftools", "view", "-h", url],
        capture_output=True,
        text=True,
        timeout=300,
        check=False,
    )
    if result.returncode != 0:
        raise RuntimeError(
            f"bcftools header read failed for {url}; stderr={result.stderr.strip()}"
        )

    contigs = set()
    for line in result.stdout.splitlines():
        if line.startswith("##contig=<ID="):
            contig = line.split("##contig=<ID=", 1)[1].split(",", 1)[0].rstrip(">")
            contigs.add(contig)

    return contigs, result.stderr.strip()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--bed", required=True, type=Path)
    parser.add_argument("--routes", required=True, type=Path)
    parser.add_argument("--out-table", required=True, type=Path)
    parser.add_argument("--out-json", required=True, type=Path)
    args = parser.parse_args()

    required_chromosomes = read_bed_chromosomes(args.bed)

    with args.routes.open(newline="", encoding="utf-8") as handle:
        route_rows = list(csv.DictReader(handle, delimiter="\t"))

    expected_fields = {
        "chromosome", "primary_genes", "genotype_vcf_url", "genotype_index_url",
        "phased_bcf_url", "phased_index_url", "route_status", "intended_use"
    }
    observed_fields = set(route_rows[0].keys()) if route_rows else set()
    missing_fields = sorted(expected_fields - observed_fields)
    if missing_fields:
        raise SystemExit(f"ERROR: route manifest missing fields: {missing_fields}")

    route_by_chromosome = {row["chromosome"]: row for row in route_rows}
    failures = []
    warnings = []
    output = []

    missing_routes = sorted(set(required_chromosomes) - set(route_by_chromosome))
    extra_routes = sorted(set(route_by_chromosome) - set(required_chromosomes))
    if missing_routes:
        failures.append(f"Missing routes for frozen BED chromosomes: {missing_routes}")
    if extra_routes:
        failures.append(f"Routes present for chromosomes outside frozen BED: {extra_routes}")

    for chrom in required_chromosomes:
        route = route_by_chromosome.get(chrom)
        if route is None:
            continue

        try:
            genotype_contigs, genotype_warning = read_remote_header(route["genotype_vcf_url"])
            phased_contigs, phased_warning = read_remote_header(route["phased_bcf_url"])
        except RuntimeError as exc:
            failures.append(f"{chrom}: {exc}")
            continue

        genotype_match = chrom in genotype_contigs
        phased_match = chrom in phased_contigs

        if genotype_warning:
            warnings.append({
                "chromosome": chrom,
                "resource": "genotype",
                "stderr": genotype_warning,
            })
        if phased_warning:
            warnings.append({
                "chromosome": chrom,
                "resource": "phased",
                "stderr": phased_warning,
            })

        output.append({
            "chromosome": chrom,
            "primary_genes": route["primary_genes"],
            "genotype_vcf_url": route["genotype_vcf_url"],
            "phased_bcf_url": route["phased_bcf_url"],
            "route_status": route["route_status"],
            "genotype_contig_match": "yes" if genotype_match else "no",
            "phased_contig_match": "yes" if phased_match else "no",
        })

        if not genotype_match:
            failures.append(f"{chrom}: BED contig absent from genotype header")
        if not phased_match:
            failures.append(f"{chrom}: BED contig absent from phased header")

    fields = [
        "chromosome", "primary_genes", "genotype_vcf_url", "phased_bcf_url",
        "route_status", "genotype_contig_match", "phased_contig_match"
    ]

    args.out_table.parent.mkdir(parents=True, exist_ok=True)
    with args.out_table.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(output)

    payload = {
        "status": (
            "PASS"
            if not failures and len(output) == len(required_chromosomes)
            else "FAIL"
        ),
        "operation": (
            "Remote-header contig compatibility check using explicit Phase B "
            "resource routes before interval extraction."
        ),
        "frozen_bed_chromosome_count": len(required_chromosomes),
        "route_manifest_chromosome_count": len(route_rows),
        "validated_resource_pairs": len(output),
        "failures": failures,
        "nonfatal_header_warnings": warnings,
        "interpretation_boundary": (
            "This check reads remote headers only. It extracts no variant records "
            "and authorises no biological interpretation."
        ),
    }

    args.out_json.parent.mkdir(parents=True, exist_ok=True)
    args.out_json.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(payload, indent=2))

    if payload["status"] != "PASS":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
