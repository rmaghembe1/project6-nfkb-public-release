#!/usr/bin/env python3
"""Validate sample identity consistency across all required genotype and phased resources."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path


def read_samples(path: Path) -> list[str]:
    samples = [
        line.strip()
        for line in path.read_text(encoding="utf-8").splitlines()
        if line.strip()
    ]
    if len(samples) != len(set(samples)):
        raise ValueError(f"Duplicate sample identifiers found in {path}")
    return samples


def read_metadata_samples(path: Path, id_field: str) -> set[str]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if id_field not in (reader.fieldnames or []):
            raise ValueError(
                f"Metadata identifier field '{id_field}' not found in {path}"
            )
        samples = [
            row[id_field].strip()
            for row in reader
            if row[id_field].strip()
        ]
    if len(samples) != len(set(samples)):
        raise ValueError("Duplicate sample identifiers found in metadata.")
    return set(samples)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--inventory", required=True, type=Path)
    parser.add_argument("--sample-dir", required=True, type=Path)
    parser.add_argument("--metadata", required=True, type=Path)
    parser.add_argument("--metadata-id-field", default="s")
    parser.add_argument("--out-summary", required=True, type=Path)
    parser.add_argument("--out-json", required=True, type=Path)
    args = parser.parse_args()

    with args.inventory.open(newline="", encoding="utf-8") as handle:
        chromosomes = [
            row["chromosome"]
            for row in csv.DictReader(handle, delimiter="\t")
        ]

    metadata = read_metadata_samples(args.metadata, args.metadata_id_field)

    genotype_sets: dict[str, set[str]] = {}
    phased_sets: dict[str, set[str]] = {}
    summary_rows: list[dict[str, object]] = []

    for chrom in chromosomes:
        genotype_file = (
            args.sample_dir /
            f"gnomad_v3.1.2_hgdp_tgp_{chrom}_genotype_samples.txt"
        )
        phased_label = "chrX_nonPAR" if chrom == "chrX" else chrom
        phased_file = (
            args.sample_dir /
            f"hgdp1kgp_{phased_label}_phased_haplotype_samples.txt"
        )

        if not genotype_file.exists():
            raise FileNotFoundError(f"Missing genotype sample list: {genotype_file}")
        if not phased_file.exists():
            raise FileNotFoundError(f"Missing phased sample list: {phased_file}")

        genotype = set(read_samples(genotype_file))
        phased = set(read_samples(phased_file))
        genotype_sets[chrom] = genotype
        phased_sets[chrom] = phased

        summary_rows.append({
            "chromosome": chrom,
            "genotype_sample_count": len(genotype),
            "phased_sample_count": len(phased),
            "phased_not_in_genotype_count": len(phased - genotype),
            "genotype_without_metadata_count": len(genotype - metadata),
            "phased_without_metadata_count": len(phased - metadata),
        })

    baseline_chromosome = "chr6" if "chr6" in chromosomes else chromosomes[0]
    genotype_baseline = genotype_sets[baseline_chromosome]
    phased_baseline = phased_sets[baseline_chromosome]

    genotype_identity_failures: list[dict[str, int | str]] = []
    phased_identity_failures: list[dict[str, int | str]] = []

    for chrom in chromosomes:
        genotype = genotype_sets[chrom]
        phased = phased_sets[chrom]

        if genotype != genotype_baseline:
            genotype_identity_failures.append({
                "chromosome": chrom,
                "baseline_only_count": len(genotype_baseline - genotype),
                "chromosome_only_count": len(genotype - genotype_baseline),
            })

        if phased != phased_baseline:
            phased_identity_failures.append({
                "chromosome": chrom,
                "baseline_only_count": len(phased_baseline - phased),
                "chromosome_only_count": len(phased - phased_baseline),
            })

    phased_outside_genotype = sorted(phased_baseline - genotype_baseline)
    genotype_without_metadata = sorted(genotype_baseline - metadata)
    phased_without_metadata = sorted(phased_baseline - metadata)
    genotype_not_in_phased = sorted(genotype_baseline - phased_baseline)

    status = (
        "PASS"
        if len(chromosomes) == 14
        and not genotype_identity_failures
        and not phased_identity_failures
        and not phased_outside_genotype
        and not phased_without_metadata
        and len(genotype_without_metadata) == 1
        and genotype_without_metadata == ["CHMI_CHMI3_WGS2"]
        and len(genotype_not_in_phased) == 60
        else "FAIL"
    )

    args.out_summary.parent.mkdir(parents=True, exist_ok=True)
    args.out_json.parent.mkdir(parents=True, exist_ok=True)

    fields = [
        "chromosome",
        "genotype_sample_count",
        "phased_sample_count",
        "phased_not_in_genotype_count",
        "genotype_without_metadata_count",
        "phased_without_metadata_count",
    ]
    with args.out_summary.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fields, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(summary_rows)

    payload = {
        "status": status,
        "validation_scope": (
            "Identity-level sample-set consistency across all 14 chromosomes "
            "required by the frozen NF-kB primary panel."
        ),
        "baseline_chromosome_for_identity_comparison": baseline_chromosome,
        "required_chromosome_count": len(chromosomes),
        "genotype_sample_count": len(genotype_baseline),
        "phased_sample_count": len(phased_baseline),
        "metadata_sample_count": len(metadata),
        "fully_interoperable_phased_sample_count": len(
            genotype_baseline & phased_baseline & metadata
        ),
        "genotype_identity_failures_across_chromosomes": genotype_identity_failures,
        "phased_identity_failures_across_chromosomes": phased_identity_failures,
        "phased_samples_absent_from_genotype": phased_outside_genotype,
        "phased_samples_absent_from_metadata": phased_without_metadata,
        "genotype_samples_absent_from_metadata": genotype_without_metadata,
        "genotype_not_in_phased_count": len(genotype_not_in_phased),
        "interpretation": (
            "The phased resource is a resource-defined filtered subset of the "
            "primary genotype resource and is fully linked to metadata; "
            "CHMI_CHMI3_WGS2 must be excluded from metadata-dependent genotype analyses."
        ),
    }

    args.out_json.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(payload, indent=2))

    if status != "PASS":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
