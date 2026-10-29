#!/usr/bin/env python3
"""Validate MANE-linked transcript/protein consequences for TNFAIP3/A20 target anchor alleles."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path


EXPECTED = {
    "rs376205580": {
        "variant": "chr6:137874871:A>G",
        "protein_change": "A20_T108A",
        "ensembl_hgvsc": "ENST00000612899.5:c.322A>G",
        "ensembl_hgvsp": "ENSP00000481570.1:p.Thr108Ala",
        "refseq_hgvsc": "NM_001270508.2:c.322A>G",
        "refseq_hgvsp": "NP_001257437.1:p.Thr108Ala",
        "amino_acids": "T/A",
        "codons": "Act/Gct",
        "protein_position": 108,
    },
    "rs141807543": {
        "variant": "chr6:137875820:A>C",
        "protein_change": "A20_I207L",
        "ensembl_hgvsc": "ENST00000612899.5:c.619A>C",
        "ensembl_hgvsp": "ENSP00000481570.1:p.Ile207Leu",
        "refseq_hgvsc": "NM_001270508.2:c.619A>C",
        "refseq_hgvsp": "NP_001257437.1:p.Ile207Leu",
        "amino_acids": "I/L",
        "codons": "Atc/Ctc",
        "protein_position": 207,
    },
}


def read_target_anchor(path: Path) -> dict[str, dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return {
            row["rsid"]: row
            for row in csv.DictReader(handle, delimiter="\t")
        }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--target-anchor", required=True, type=Path)
    parser.add_argument("--vep-response", required=True, type=Path)
    parser.add_argument("--out-table", required=True, type=Path)
    parser.add_argument("--out-qc", required=True, type=Path)
    args = parser.parse_args()

    anchor = read_target_anchor(args.target_anchor)
    response = json.loads(args.vep_response.read_text(encoding="utf-8"))

    failures: list[str] = []
    output_rows: list[dict[str, str | int]] = []

    if set(anchor) != set(EXPECTED):
        failures.append(
            f"Target-anchor rsID set differs from expected set: {sorted(anchor)}"
        )

    records_by_rsid: dict[str, dict] = {}
    for record in response:
        input_text = record.get("input", "")
        for rsid in EXPECTED:
            if rsid in input_text:
                records_by_rsid[rsid] = record

    for rsid, expected in EXPECTED.items():
        if rsid not in records_by_rsid:
            failures.append(f"{rsid}: absent from VEP response")
            continue

        target = anchor.get(rsid)
        if target is None:
            failures.append(f"{rsid}: absent from parsed genomic anchor table")
            continue

        observed_variant = (
            f'{target["chromosome"]}:{target["position_1based"]}:'
            f'{target["ref"]}>{target["target_alt"]}'
        )
        if observed_variant != expected["variant"]:
            failures.append(
                f"{rsid}: genomic target allele mismatch; observed {observed_variant}"
            )

        record = records_by_rsid[rsid]
        if record.get("most_severe_consequence") != "missense_variant":
            failures.append(f"{rsid}: most severe consequence is not missense_variant")

        consequences = record.get("transcript_consequences", [])
        ensembl = [
            item for item in consequences
            if item.get("hgvsc") == expected["ensembl_hgvsc"]
            and item.get("hgvsp") == expected["ensembl_hgvsp"]
        ]
        refseq = [
            item for item in consequences
            if item.get("hgvsc") == expected["refseq_hgvsc"]
            and item.get("hgvsp") == expected["refseq_hgvsp"]
        ]

        if len(ensembl) != 1:
            failures.append(f"{rsid}: expected exactly one Ensembl MANE-linked consequence")
            continue
        if len(refseq) != 1:
            failures.append(f"{rsid}: expected exactly one RefSeq MANE-linked consequence")
            continue

        ens = ensembl[0]
        ref = refseq[0]

        for label, item in [("Ensembl", ens), ("RefSeq", ref)]:
            if item.get("consequence_terms") != ["missense_variant"]:
                failures.append(f"{rsid}: {label} consequence term mismatch")
            if item.get("amino_acids") != expected["amino_acids"]:
                failures.append(f"{rsid}: {label} amino-acid change mismatch")
            if item.get("codons") != expected["codons"]:
                failures.append(f"{rsid}: {label} codon change mismatch")
            if item.get("protein_start") != expected["protein_position"]:
                failures.append(f"{rsid}: {label} protein position mismatch")
            if item.get("canonical") != 1:
                failures.append(f"{rsid}: {label} consequence is not canonical/MANE-linked")

        output_rows.append({
            "anchor_id": "A1",
            "gene_symbol": "TNFAIP3",
            "anchor_haplotype": "A20_T108A-I207L_Denisovan-derived",
            "protein_change": expected["protein_change"],
            "rsid": rsid,
            "grch38_target_variant": expected["variant"],
            "vep_most_severe_consequence": record.get("most_severe_consequence", ""),
            "mane_refseq_transcript": "NM_001270508.2",
            "mane_refseq_hgvsc": expected["refseq_hgvsc"],
            "mane_refseq_protein": "NP_001257437.1",
            "mane_refseq_hgvsp": expected["refseq_hgvsp"],
            "mane_ensembl_transcript": "ENST00000612899.5",
            "mane_ensembl_hgvsc": expected["ensembl_hgvsc"],
            "mane_ensembl_protein": "ENSP00000481570.1",
            "mane_ensembl_hgvsp": expected["ensembl_hgvsp"],
            "amino_acid_change": expected["amino_acids"],
            "codon_change": expected["codons"],
            "protein_position": expected["protein_position"],
            "confirmation_status": "PASS_MANE_linked_missense_consequence_confirmed",
        })

    output_rows.sort(key=lambda row: int(str(row["grch38_target_variant"]).split(":")[1]))

    status = "PASS" if len(output_rows) == 2 and not failures else "FAIL"

    args.out_table.parent.mkdir(parents=True, exist_ok=True)
    args.out_qc.parent.mkdir(parents=True, exist_ok=True)

    fields = list(output_rows[0].keys()) if output_rows else [
        "anchor_id", "gene_symbol", "anchor_haplotype", "protein_change",
        "rsid", "grch38_target_variant", "confirmation_status"
    ]
    with args.out_table.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fields, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(output_rows)

    payload = {
        "status": status,
        "anchor_id": "A1",
        "gene_symbol": "TNFAIP3",
        "validated_target_variant_count": len(output_rows),
        "mane_refseq_transcript": "NM_001270508.2",
        "mane_ensembl_transcript": "ENST00000612899.5",
        "validated_consequences": [
            {
                "rsid": row["rsid"],
                "grch38_target_variant": row["grch38_target_variant"],
                "refseq_hgvsc": row["mane_refseq_hgvsc"],
                "refseq_hgvsp": row["mane_refseq_hgvsp"],
            }
            for row in output_rows
        ],
        "interpretation": (
            "The target TNFAIP3 alleles defining the published A20 T108A-I207L "
            "anchor are confirmed as MANE-linked missense variants; non-target "
            "dbSNP alternate alleles remain outside the target anchor definition."
        ),
        "failures": failures,
    }

    args.out_qc.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(payload, indent=2))

    if status != "PASS":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
