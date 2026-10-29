#!/usr/bin/env python3
"""Resolve GRCh38 genomic placement and target alleles for the TNFAIP3/A20 anchor."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path


TARGETS = {
    "rs376205580": {
        "protein_change": "A20_T108A",
        "target_ref": "A",
        "target_alt": "G",
        "anchor_component": "Denisovan-derived A20 T108A component",
    },
    "rs141807543": {
        "protein_change": "A20_I207L",
        "target_ref": "A",
        "target_alt": "C",
        "anchor_component": "Denisovan-derived A20 I207L component",
    },
}


def read_tnfaip3_interval(path: Path) -> dict[str, str]:
    with path.open(newline="", encoding="utf-8") as handle:
        rows = [
            row for row in csv.DictReader(handle, delimiter="\t")
            if row["gene_symbol"] == "TNFAIP3"
        ]
    if len(rows) != 1:
        raise ValueError(f"Expected one TNFAIP3 interval in {path}; found {len(rows)}")
    return rows[0]


def grch38_alleles_from_refsnp(path: Path) -> dict[str, object]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    rsid = f"rs{payload['refsnp_id']}"

    placements = payload["primary_snapshot_data"].get("placements_with_allele", [])
    chr6_placements = [
        placement for placement in placements
        if placement.get("seq_id") == "NC_000006.12"
    ]

    if len(chr6_placements) != 1:
        raise ValueError(
            f"Expected exactly one GRCh38 chr6 placement for {rsid}; "
            f"found {len(chr6_placements)}"
        )

    placement = chr6_placements[0]
    assemblies = [
        item.get("assembly_name", "")
        for item in placement.get("placement_annot", {})
        .get("seq_id_traits_by_assembly", [])
    ]

    non_reference_alleles = []
    for allele_record in placement.get("alleles", []):
        spdi = allele_record.get("allele", {}).get("spdi", {})
        if not spdi:
            continue

        ref = spdi.get("deleted_sequence")
        alt = spdi.get("inserted_sequence")
        pos0 = spdi.get("position")

        if ref == alt:
            continue

        non_reference_alleles.append({
            "position_1based": int(pos0) + 1,
            "ref": ref,
            "alt": alt,
        })

    unique_positions = {row["position_1based"] for row in non_reference_alleles}
    unique_refs = {row["ref"] for row in non_reference_alleles}

    if len(unique_positions) != 1 or len(unique_refs) != 1:
        raise ValueError(f"Unexpected multi-position or multi-reference record for {rsid}")

    return {
        "rsid": rsid,
        "seq_id": placement.get("seq_id"),
        "assemblies": assemblies,
        "is_ptlp": placement.get("is_ptlp"),
        "position_1based": next(iter(unique_positions)),
        "ref": next(iter(unique_refs)),
        "non_reference_alleles": sorted(
            row["alt"] for row in non_reference_alleles
        ),
        "source_json": str(path),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--dbsnp-dir", required=True, type=Path)
    parser.add_argument("--full-locus", required=True, type=Path)
    parser.add_argument("--mane-select", required=True, type=Path)
    parser.add_argument("--out-all-alleles", required=True, type=Path)
    parser.add_argument("--out-target-anchor", required=True, type=Path)
    parser.add_argument("--out-qc", required=True, type=Path)
    args = parser.parse_args()

    full_locus = read_tnfaip3_interval(args.full_locus)
    mane_select = read_tnfaip3_interval(args.mane_select)

    all_rows = []
    target_rows = []
    failures = []

    for rsid, expected in TARGETS.items():
        record = grch38_alleles_from_refsnp(args.dbsnp_dir / f"{rsid}.json")

        pos = int(record["position_1based"])
        ref = str(record["ref"])
        alts = list(record["non_reference_alleles"])
        target_alt = expected["target_alt"]

        target_present = ref == expected["target_ref"] and target_alt in alts

        in_full_locus = (
            record["seq_id"] == "NC_000006.12"
            and full_locus["chromosome"] == "chr6"
            and int(full_locus["start_1based"]) <= pos <= int(full_locus["end_1based"])
        )
        in_mane_select = (
            record["seq_id"] == "NC_000006.12"
            and mane_select["chromosome"] == "chr6"
            and int(mane_select["start_1based"]) <= pos <= int(mane_select["end_1based"])
        )

        for alt in alts:
            all_rows.append({
                "gene_symbol": "TNFAIP3",
                "rsid": rsid,
                "protein_anchor_component": expected["protein_change"],
                "chromosome": "chr6",
                "position_1based": pos,
                "ref": ref,
                "alt": alt,
                "is_target_anchor_allele": "yes" if alt == target_alt else "no",
                "assembly": ";".join(record["assemblies"]),
                "seq_id": record["seq_id"],
                "is_ptlp": str(record["is_ptlp"]),
                "interpretation": (
                    "Published A20 T108A-I207L target allele"
                    if alt == target_alt
                    else "Additional dbSNP alternate allele; not assigned to target anchor"
                ),
            })

        target_rows.append({
            "anchor_id": "A1",
            "gene_symbol": "TNFAIP3",
            "anchor_haplotype": "A20_T108A-I207L_Denisovan-derived",
            "protein_change": expected["protein_change"],
            "rsid": rsid,
            "chromosome": "chr6",
            "position_1based": pos,
            "ref": ref,
            "target_alt": target_alt,
            "other_dbsnp_alts": ",".join(alt for alt in alts if alt != target_alt),
            "assembly": ";".join(record["assemblies"]),
            "seq_id": record["seq_id"],
            "contained_in_full_tnfaip3_locus": "yes" if in_full_locus else "no",
            "contained_in_mane_select_transcript": "yes" if in_mane_select else "no",
            "resolution_status": (
                "genomic_position_and_target_allele_resolved;"
                "transcript_consequence_confirmation_pending"
            ),
        })

        if not target_present:
            failures.append(f"{rsid}: expected target allele {ref}>{target_alt} not found")
        if not in_full_locus:
            failures.append(f"{rsid}: variant is outside frozen full TNFAIP3 locus")
        if not in_mane_select:
            failures.append(f"{rsid}: variant is outside frozen MANE Select TNFAIP3 transcript")

    all_rows.sort(key=lambda row: (int(row["position_1based"]), row["alt"]))
    target_rows.sort(key=lambda row: int(row["position_1based"]))

    args.out_all_alleles.parent.mkdir(parents=True, exist_ok=True)
    args.out_target_anchor.parent.mkdir(parents=True, exist_ok=True)
    args.out_qc.parent.mkdir(parents=True, exist_ok=True)

    all_fields = list(all_rows[0].keys())
    with args.out_all_alleles.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=all_fields, delimiter="\t",
                                lineterminator="\n")
        writer.writeheader()
        writer.writerows(all_rows)

    target_fields = list(target_rows[0].keys())
    with args.out_target_anchor.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=target_fields, delimiter="\t",
                                lineterminator="\n")
        writer.writeheader()
        writer.writerows(target_rows)

    payload = {
        "status": "PASS" if not failures else "FAIL",
        "anchor_id": "A1",
        "gene_symbol": "TNFAIP3",
        "anchor_haplotype": "A20_T108A-I207L_Denisovan-derived",
        "resolved_variant_count": len(target_rows),
        "target_anchor_alleles": [
            {
                "rsid": row["rsid"],
                "protein_change": row["protein_change"],
                "grch38_variant": (
                    f'{row["chromosome"]}:{row["position_1based"]}:'
                    f'{row["ref"]}>{row["target_alt"]}'
                ),
            }
            for row in target_rows
        ],
        "all_dbsnp_nonreference_allele_count": len(all_rows),
        "containment_policy": (
            "Both target variants must fall within the frozen full TNFAIP3 locus "
            "and MANE Select transcript interval."
        ),
        "interpretation_boundary": (
            "The target genomic alleles are selected for the published A20 "
            "T108A-I207L anchor; additional dbSNP alternate alleles are retained "
            "for audit but are not assigned to the target haplotype."
        ),
        "pending_requirement": (
            "Confirm target transcript/protein consequences before anchor-based "
            "frequency or haplotype interpretation."
        ),
        "failures": failures,
    }

    args.out_qc.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(payload, indent=2))

    if failures:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
