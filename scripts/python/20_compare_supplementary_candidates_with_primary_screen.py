#!/usr/bin/env python3
"""Compare supplementary direct-design candidates against the primary classified competitor set."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
from pathlib import Path


def normalise_doi(value: str) -> str:
    return value.strip().lower().replace("https://doi.org/", "").replace("http://doi.org/", "")


def normalise_title(value: str) -> str:
    value = re.sub(r"<[^>]+>", "", value)
    value = value.lower().replace("&lt;i&gt;", "").replace("&lt;/i&gt;", "")
    value = re.sub(r"[^a-z0-9]+", " ", value)
    return re.sub(r"\s+", " ", value).strip()


def record_keys(row: dict[str, str]) -> set[str]:
    keys = set()
    doi = normalise_doi(row.get("doi", ""))
    pmid = row.get("pmid", "").strip()
    title = normalise_title(row.get("title", ""))
    if doi:
        keys.add(f"doi:{doi}")
    if pmid:
        keys.add(f"pmid:{pmid}")
    if title:
        keys.add(f"title:{hashlib.sha256(title.encode('utf-8')).hexdigest()}")
    return keys


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--primary-classified", required=True, type=Path)
    parser.add_argument("--supplementary-unique", required=True, type=Path)
    parser.add_argument("--out-table", required=True, type=Path)
    parser.add_argument("--out-new", required=True, type=Path)
    parser.add_argument("--out-json", required=True, type=Path)
    args = parser.parse_args()

    with args.primary_classified.open(newline="", encoding="utf-8") as handle:
        primary = list(csv.DictReader(handle, delimiter="\t"))
    with args.supplementary_unique.open(newline="", encoding="utf-8") as handle:
        supplementary = list(csv.DictReader(handle, delimiter="\t"))

    primary_by_key: dict[str, dict[str, str]] = {}
    for row in primary:
        for key in record_keys(row):
            primary_by_key[key] = row

    comparison_rows = []
    new_rows = []

    for index, row in enumerate(supplementary, start=1):
        matched_primary = None
        matched_key = ""
        for key in record_keys(row):
            if key in primary_by_key:
                matched_primary = primary_by_key[key]
                matched_key = key
                break

        comparison = {
            "supplementary_id": f"SP{index:03d}",
            "supplementary_query_ids": row["query_ids"],
            "title": row["title"],
            "year": row["year"],
            "doi": row["doi"],
            "pmid": row["pmid"],
            "already_in_primary_screen": "yes" if matched_primary else "no",
            "matched_primary_screening_id": (
                matched_primary["screening_id"] if matched_primary else ""
            ),
            "matched_primary_classification": (
                matched_primary["screening_classification"] if matched_primary else ""
            ),
            "matched_by_key": matched_key,
            "supplementary_screening_classification": (
                matched_primary["screening_classification"] if matched_primary else ""
            ),
            "supplementary_screening_reason": (
                "Previously classified in primary formal screening set."
                if matched_primary else ""
            ),
            "direct_competitor": (
                "yes"
                if matched_primary
                and matched_primary["screening_classification"] == "DIRECT_COMPETITOR"
                else "no" if matched_primary else ""
            ),
            "abstract": row["abstract"],
        }
        comparison_rows.append(comparison)
        if not matched_primary:
            new_rows.append(comparison)

    fields = list(comparison_rows[0].keys())
    args.out_table.parent.mkdir(parents=True, exist_ok=True)
    args.out_new.parent.mkdir(parents=True, exist_ok=True)
    args.out_json.parent.mkdir(parents=True, exist_ok=True)

    with args.out_table.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(comparison_rows)

    with args.out_new.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(new_rows)

    payload = {
        "status": "SUPPLEMENTARY_PRIMARY_OVERLAP_ASSESSED_NEW_RECORD_SCREENING_PENDING",
        "supplementary_unique_candidate_count": len(supplementary),
        "already_screened_in_primary_set_count": len(comparison_rows) - len(new_rows),
        "new_supplementary_candidate_count": len(new_rows),
        "already_screened_records": [
            {
                "supplementary_id": row["supplementary_id"],
                "matched_primary_screening_id": row["matched_primary_screening_id"],
                "classification": row["matched_primary_classification"],
                "title": row["title"],
            }
            for row in comparison_rows
            if row["already_in_primary_screen"] == "yes"
        ],
        "new_supplementary_records": [
            {
                "supplementary_id": row["supplementary_id"],
                "query_ids": row["supplementary_query_ids"],
                "title": row["title"],
                "doi": row["doi"],
            }
            for row in new_rows
        ],
        "interpretation_boundary": (
            "Novelty assessment cannot be closed until genuinely new supplementary "
            "records have been screened and classified."
        ),
    }

    args.out_json.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
