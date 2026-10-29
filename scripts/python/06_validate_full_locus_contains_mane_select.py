#!/usr/bin/env python3
"""Validate that MANE Select transcript intervals are contained within matching full gene loci."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path


def read_by_gene(path: Path) -> dict[str, dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return {
            row["gene_symbol"]: row
            for row in csv.DictReader(handle, delimiter="\t")
        }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--full-loci", required=True, type=Path)
    parser.add_argument("--mane-select", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args()

    full_loci = read_by_gene(args.full_loci)
    mane_select = read_by_gene(args.mane_select)

    missing_or_unmatched = sorted(set(full_loci) ^ set(mane_select))
    failures: list[dict[str, str]] = []

    for gene in sorted(set(full_loci) & set(mane_select)):
        locus = full_loci[gene]
        transcript = mane_select[gene]
        passed = (
            locus["chromosome"] == transcript["chromosome"]
            and locus["strand"] == transcript["strand"]
            and int(locus["start_1based"]) <= int(transcript["start_1based"])
            and int(locus["end_1based"]) >= int(transcript["end_1based"])
        )
        if not passed:
            failures.append({
                "gene_symbol": gene,
                "full_locus_chromosome": locus["chromosome"],
                "mane_chromosome": transcript["chromosome"],
                "full_locus_start": locus["start_1based"],
                "mane_start": transcript["start_1based"],
                "full_locus_end": locus["end_1based"],
                "mane_end": transcript["end_1based"],
                "full_locus_strand": locus["strand"],
                "mane_strand": transcript["strand"],
            })

    status = (
        "PASS"
        if len(full_loci) == 24
        and len(mane_select) == 24
        and not missing_or_unmatched
        and not failures
        else "FAIL"
    )

    payload = {
        "status": status,
        "validation_rule": (
            "Each MANE Select transcript must be contained within the matching "
            "full annotated gene locus on the same chromosome and strand."
        ),
        "full_gene_locus_count": len(full_loci),
        "mane_select_transcript_count": len(mane_select),
        "missing_or_unmatched_genes": missing_or_unmatched,
        "containment_failures": failures,
    }

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(payload, indent=2))

    if status != "PASS":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
