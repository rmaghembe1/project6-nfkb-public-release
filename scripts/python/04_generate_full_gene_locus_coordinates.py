#!/usr/bin/env python3
"""Generate full annotated GRCh38 gene-locus spans for the primary NF-kB panel from the MANE RefSeq GTF."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import re
from pathlib import Path


def sha256sum(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def parse_attributes(text: str) -> dict[str, str]:
    attributes: dict[str, str] = {}
    for key, value in re.findall(r'(\S+)\s+"([^"]*)"', text):
        attributes[key] = value
    return attributes


def chromosome_sort_key(chrom: str) -> tuple[int, str]:
    cleaned = chrom.removeprefix("chr")
    if cleaned.isdigit():
        return (int(cleaned), "")
    order = {"X": 23, "Y": 24, "M": 25, "MT": 25}
    return (order.get(cleaned, 99), cleaned)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--panel", required=True, type=Path)
    parser.add_argument("--gtf", required=True, type=Path)
    parser.add_argument("--out-tsv", required=True, type=Path)
    parser.add_argument("--out-bed", required=True, type=Path)
    parser.add_argument("--qc", required=True, type=Path)
    args = parser.parse_args()

    with args.panel.open(encoding="utf-8", newline="") as handle:
        panel_rows = list(csv.DictReader(handle, delimiter="\t"))

    primary_rows = [
        row for row in panel_rows
        if row["analysis_tier"].strip().lower() == "primary"
    ]
    primary_symbols = {row["gene_symbol"] for row in primary_rows}
    module_by_gene = {row["gene_symbol"]: row["module_id"] for row in primary_rows}
    module_name_by_gene = {row["gene_symbol"]: row["module_name"] for row in primary_rows}

    found: dict[str, list[dict[str, object]]] = {gene: [] for gene in primary_symbols}

    with gzip.open(args.gtf, "rt", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2] != "gene":
                continue

            chrom, source, _, start, end, _, strand, _, attribute_text = fields
            attributes = parse_attributes(attribute_text)
            gene = attributes.get("gene")

            if gene in primary_symbols:
                found[gene].append({
                    "gene_symbol": gene,
                    "module_id": module_by_gene[gene],
                    "module_name": module_name_by_gene[gene],
                    "chromosome": chrom,
                    "start_1based": int(start),
                    "end_1based": int(end),
                    "strand": strand,
                    "source": source,
                    "gene_biotype": attributes.get("gene_biotype", ""),
                    "gene_id": attributes.get("gene_id", ""),
                })

    missing = sorted(gene for gene, records in found.items() if len(records) == 0)
    duplicated = sorted(gene for gene, records in found.items() if len(records) > 1)
    coordinates = [records[0] for records in found.values() if len(records) == 1]
    coordinates.sort(key=lambda row: (chromosome_sort_key(str(row["chromosome"])),
                                      int(row["start_1based"])))

    status = "PASS" if len(primary_symbols) == 24 and not missing and not duplicated else "FAIL"

    args.out_tsv.parent.mkdir(parents=True, exist_ok=True)
    args.out_bed.parent.mkdir(parents=True, exist_ok=True)
    args.qc.parent.mkdir(parents=True, exist_ok=True)

    tsv_fields = [
        "module_id", "module_name", "gene_symbol", "chromosome",
        "start_1based", "end_1based", "strand", "source",
        "gene_biotype", "gene_id", "coordinate_definition"
    ]

    with args.out_tsv.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=tsv_fields, delimiter="\t",
                                lineterminator="\n")
        writer.writeheader()
        for row in coordinates:
            output = dict(row)
            output["coordinate_definition"] = "MANE_v1.5_RefSeq_GTF_full_gene_locus_GRCh38_1based_inclusive"
            writer.writerow(output)

    with args.out_bed.open("w", encoding="utf-8", newline="") as handle:
        for row in coordinates:
            handle.write(
                f'{row["chromosome"]}\t{int(row["start_1based"]) - 1}\t'
                f'{row["end_1based"]}\t{row["gene_symbol"]}\t'
                f'{row["module_id"]}\t{row["strand"]}\n'
            )

    qc_payload = {
        "status": status,
        "coordinate_policy": "Full annotated gene-locus spans from the MANE v1.5 RefSeq GTF gene feature; no external regulatory flanks added.",
        "coordinate_system_tsv": "GRCh38, 1-based inclusive",
        "coordinate_system_bed": "GRCh38, 0-based half-open",
        "panel_primary_gene_count": len(primary_symbols),
        "genes_with_single_mane_gene_record": len(coordinates),
        "missing_primary_genes": missing,
        "duplicated_primary_gene_records": duplicated,
        "source_gtf": str(args.gtf),
        "source_gtf_sha256": sha256sum(args.gtf),
    }

    args.qc.write_text(json.dumps(qc_payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(qc_payload, indent=2))

    if status != "PASS":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
