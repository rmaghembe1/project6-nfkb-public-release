#!/usr/bin/env python3
"""Compare sample sets across gnomAD genotype VCF, phased BCF and metadata resources."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path


def read_sample_list(path: Path) -> list[str]:
    values = [
        line.strip()
        for line in path.read_text(encoding="utf-8").splitlines()
        if line.strip()
    ]
    if len(values) != len(set(values)):
        raise ValueError(f"Duplicate sample IDs detected in {path}")
    return values


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--genotype-samples", required=True, type=Path)
    parser.add_argument("--phased-samples", required=True, type=Path)
    parser.add_argument("--metadata", required=True, type=Path)
    parser.add_argument("--metadata-id-field", default="s")
    parser.add_argument("--out-json", required=True, type=Path)
    parser.add_argument("--out-membership", required=True, type=Path)
    parser.add_argument("--out-genotype-only", required=True, type=Path)
    parser.add_argument("--out-phased-only", required=True, type=Path)
    parser.add_argument("--out-metadata-only", required=True, type=Path)
    args = parser.parse_args()

    genotype_order = read_sample_list(args.genotype_samples)
    phased_order = read_sample_list(args.phased_samples)

    with args.metadata.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = reader.fieldnames or []
        if args.metadata_id_field not in fields:
            raise SystemExit(
                f"Metadata ID field '{args.metadata_id_field}' not found. "
                f"Available fields: {fields}"
            )
        metadata_order = [
            row[args.metadata_id_field].strip()
            for row in reader
            if row[args.metadata_id_field].strip()
        ]

    if len(metadata_order) != len(set(metadata_order)):
        raise ValueError("Duplicate metadata sample IDs detected.")

    genotype = set(genotype_order)
    phased = set(phased_order)
    metadata = set(metadata_order)

    shared_all = genotype & phased & metadata
    genotype_only_vs_phased = genotype - phased
    phased_only_vs_genotype = phased - genotype
    metadata_only_vs_genotype = metadata - genotype
    genotype_missing_metadata = genotype - metadata
    phased_missing_metadata = phased - metadata

    args.out_json.parent.mkdir(parents=True, exist_ok=True)
    args.out_membership.parent.mkdir(parents=True, exist_ok=True)

    all_ids = sorted(genotype | phased | metadata)
    with args.out_membership.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["sample_id", "in_genotype", "in_phased", "in_metadata"])
        for sample in all_ids:
            writer.writerow([
                sample,
                "yes" if sample in genotype else "no",
                "yes" if sample in phased else "no",
                "yes" if sample in metadata else "no",
            ])

    for path, values in [
        (args.out_genotype_only, genotype_only_vs_phased),
        (args.out_phased_only, phased_only_vs_genotype),
        (args.out_metadata_only, metadata_only_vs_genotype),
    ]:
        path.write_text(
            "".join(f"{sample}\n" for sample in sorted(values)),
            encoding="utf-8"
        )

    payload = {
        "status": "PASS" if not phased_only_vs_genotype else "REVIEW",
        "interpretation_rule": (
            "Phased samples may be a filtered subset of genotype samples; "
            "any phased sample absent from the primary genotype set requires review."
        ),
        "genotype_sample_count": len(genotype),
        "phased_sample_count": len(phased),
        "metadata_sample_count": len(metadata),
        "shared_across_all_three_count": len(shared_all),
        "genotype_not_in_phased_count": len(genotype_only_vs_phased),
        "phased_not_in_genotype_count": len(phased_only_vs_genotype),
        "metadata_not_in_genotype_count": len(metadata_only_vs_genotype),
        "genotype_without_metadata_count": len(genotype_missing_metadata),
        "phased_without_metadata_count": len(phased_missing_metadata),
        "metadata_id_field": args.metadata_id_field,
        "outputs": {
            "membership_table": str(args.out_membership),
            "genotype_not_in_phased": str(args.out_genotype_only),
            "phased_not_in_genotype": str(args.out_phased_only),
            "metadata_not_in_genotype": str(args.out_metadata_only),
        },
    }

    args.out_json.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(payload, indent=2))

    if phased_only_vs_genotype:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
