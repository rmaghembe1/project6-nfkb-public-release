#!/usr/bin/env python3
"""Generate GRCh38 MANE Select transcript intervals for the primary NF-kB panel."""

from __future__ import annotations

import argparse
import csv
import gzip
import json
from pathlib import Path


ACCESSION_TO_CHROM = {
    "NC_000001.11": "chr1", "NC_000002.12": "chr2",
    "NC_000003.12": "chr3", "NC_000004.12": "chr4",
    "NC_000005.10": "chr5", "NC_000006.12": "chr6",
    "NC_000007.14": "chr7", "NC_000008.11": "chr8",
    "NC_000009.12": "chr9", "NC_000010.11": "chr10",
    "NC_000011.10": "chr11", "NC_000012.12": "chr12",
    "NC_000013.11": "chr13", "NC_000014.9": "chr14",
    "NC_000015.10": "chr15", "NC_000016.10": "chr16",
    "NC_000017.11": "chr17", "NC_000018.10": "chr18",
    "NC_000019.10": "chr19", "NC_000020.11": "chr20",
    "NC_000021.9": "chr21", "NC_000022.11": "chr22",
    "NC_000023.11": "chrX", "NC_000024.10": "chrY",
}


def chrom_key(chrom: str) -> int:
    label = chrom.removeprefix("chr")
    return 23 if label == "X" else 24 if label == "Y" else int(label)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--panel", required=True, type=Path)
    parser.add_argument("--mane-summary", required=True, type=Path)
    parser.add_argument("--out-tsv", required=True, type=Path)
    parser.add_argument("--out-bed", required=True, type=Path)
    parser.add_argument("--qc", required=True, type=Path)
    args = parser.parse_args()

    with args.panel.open(newline="", encoding="utf-8") as handle:
        primary = {
            row["gene_symbol"]: row
            for row in csv.DictReader(handle, delimiter="\t")
            if row["analysis_tier"].strip().lower() == "primary"
        }

    select_by_gene: dict[str, list[dict[str, str]]] = {
        gene: [] for gene in primary
    }

    with gzip.open(args.mane_summary, "rt", newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            gene = row["symbol"]
            if gene in primary and row["MANE_status"] == "MANE Select":
                select_by_gene[gene].append(row)

    missing = sorted(gene for gene, rows in select_by_gene.items() if not rows)
    nonunique = sorted(gene for gene, rows in select_by_gene.items() if len(rows) != 1)

    rows_out: list[dict[str, str]] = []
    for gene, rows in select_by_gene.items():
        if len(rows) != 1:
            continue
        mane = rows[0]
        chrom = ACCESSION_TO_CHROM.get(mane["GRCh38_chr"], mane["GRCh38_chr"])
        rows_out.append({
            "module_id": primary[gene]["module_id"],
            "module_name": primary[gene]["module_name"],
            "gene_symbol": gene,
            "chromosome": chrom,
            "start_1based": mane["chr_start"],
            "end_1based": mane["chr_end"],
            "strand": mane["chr_strand"],
            "RefSeq_nuc": mane["RefSeq_nuc"],
            "RefSeq_prot": mane["RefSeq_prot"],
            "Ensembl_Gene": mane["Ensembl_Gene"],
            "Ensembl_nuc": mane["Ensembl_nuc"],
            "Ensembl_prot": mane["Ensembl_prot"],
            "HGNC_ID": mane["HGNC_ID"],
            "MANE_status": mane["MANE_status"],
            "coordinate_definition": "MANE_v1.5_Select_transcript_GRCh38_1based_inclusive",
        })

    rows_out.sort(key=lambda row: (chrom_key(row["chromosome"]),
                                    int(row["start_1based"])))

    status = "PASS" if len(rows_out) == 24 and not missing and not nonunique else "FAIL"

    args.out_tsv.parent.mkdir(parents=True, exist_ok=True)
    args.out_bed.parent.mkdir(parents=True, exist_ok=True)
    args.qc.parent.mkdir(parents=True, exist_ok=True)

    fields = list(rows_out[0].keys()) if rows_out else []
    with args.out_tsv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t",
                                lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows_out)

    with args.out_bed.open("w", encoding="utf-8") as handle:
        for row in rows_out:
            handle.write(
                f'{row["chromosome"]}\t{int(row["start_1based"]) - 1}\t'
                f'{row["end_1based"]}\t{row["gene_symbol"]}\t'
                f'{row["module_id"]}\t{row["strand"]}\n'
            )

    payload = {
        "status": status,
        "mane_select_transcript_count": len(rows_out),
        "missing_primary_genes": missing,
        "genes_with_nonunique_mane_select": nonunique,
        "coordinate_policy": "MANE Select transcript coordinates derived from the MANE v1.5 summary file.",
        "coordinate_system_tsv": "GRCh38, 1-based inclusive",
        "coordinate_system_bed": "GRCh38, 0-based half-open",
    }

    args.qc.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(payload, indent=2))

    if status != "PASS":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
