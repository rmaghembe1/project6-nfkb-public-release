#!/usr/bin/env python3
"""Compare priority citation-chain candidates against previously classified competitor records."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
from pathlib import Path


def normalise_doi(value: str) -> str:
    return (
        (value or "")
        .strip()
        .lower()
        .replace("https://doi.org/", "")
        .replace("http://doi.org/", "")
    )


def normalise_title(value: str) -> str:
    value = re.sub(r"<[^>]+>", "", value or "")
    value = value.replace("&lt;i&gt;", "").replace("&lt;/i&gt;", "")
    value = value.replace("&gt;", ">").replace("&lt;", "<")
    value = value.lower()
    value = re.sub(r"[^a-z0-9]+", " ", value)
    return re.sub(r"\s+", " ", value).strip()


def keys(row: dict[str, str]) -> set[str]:
    output = set()
    doi = normalise_doi(row.get("doi", ""))
    title = normalise_title(row.get("title", ""))
    if doi:
        output.add(f"doi:{doi}")
    if title:
        output.add(f"title:{hashlib.sha256(title.encode('utf-8')).hexdigest()}")
    return output


def load_table(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--citation-priority", required=True, type=Path)
    parser.add_argument("--primary-classified", required=True, type=Path)
    parser.add_argument("--supplementary-classified", required=True, type=Path)
    parser.add_argument("--out-comparison", required=True, type=Path)
    parser.add_argument("--out-new", required=True, type=Path)
    parser.add_argument("--out-json", required=True, type=Path)
    args = parser.parse_args()

    citation_rows = load_table(args.citation_priority)
    primary_rows = load_table(args.primary_classified)
    supplementary_rows = load_table(args.supplementary_classified)

    prior_index: dict[str, dict[str, str]] = {}

    for row in primary_rows:
        record = {
            "prior_source": "PRIMARY_FORMAL_SCREEN",
            "prior_id": row["screening_id"],
            "prior_classification": row["screening_classification"],
            "prior_reason": row["screening_reason"],
            "title": row["title"],
            "doi": row["doi"],
        }
        for key in keys(row):
            prior_index[key] = record

    for row in supplementary_rows:
        record = {
            "prior_source": "SUPPLEMENTARY_DIRECT_DESIGN_SCREEN",
            "prior_id": row["supplementary_id"],
            "prior_classification": row["supplementary_screening_classification"],
            "prior_reason": row["supplementary_screening_reason"],
            "title": row["title"],
            "doi": row["doi"],
        }
        for key in keys(row):
            prior_index[key] = record

    comparison = []
    new_records = []

    for i, row in enumerate(citation_rows, start=1):
        match = None
        matched_key = ""
        for key in keys(row):
            if key in prior_index:
                match = prior_index[key]
                matched_key = key
                break

        out = dict(row)
        out["citation_screening_id"] = f"CC{i:03d}"
        out["already_screened"] = "yes" if match else "no"
        out["matched_prior_source"] = match["prior_source"] if match else ""
        out["matched_prior_id"] = match["prior_id"] if match else ""
        out["matched_prior_classification"] = (
            match["prior_classification"] if match else ""
        )
        out["matched_prior_reason"] = match["prior_reason"] if match else ""
        out["matched_by_key"] = matched_key
        out["citation_manual_classification"] = (
            match["prior_classification"] if match else ""
        )
        out["citation_manual_reason"] = (
            "Previously adjudicated in formal or supplementary screening."
            if match else ""
        )
        out["direct_competitor"] = (
            "yes"
            if match and match["prior_classification"] == "DIRECT_COMPETITOR"
            else "no" if match else ""
        )
        comparison.append(out)
        if not match:
            new_records.append(out)

    fields = list(comparison[0].keys())
    args.out_comparison.parent.mkdir(parents=True, exist_ok=True)
    args.out_new.parent.mkdir(parents=True, exist_ok=True)
    args.out_json.parent.mkdir(parents=True, exist_ok=True)

    with args.out_comparison.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fields, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(comparison)

    with args.out_new.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fields, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(new_records)

    payload = {
        "status": "CITATION_PRIORITY_COMPARISON_COMPLETE_NEW_RECORD_SCREENING_PENDING",
        "priority_citation_candidate_count": len(citation_rows),
        "already_classified_in_prior_screening_count": len(comparison) - len(new_records),
        "genuinely_new_priority_citation_candidate_count": len(new_records),
        "already_classified_records": [
            {
                "citation_screening_id": row["citation_screening_id"],
                "matched_prior_id": row["matched_prior_id"],
                "matched_prior_classification": row["matched_prior_classification"],
                "title": row["title"],
            }
            for row in comparison
            if row["already_screened"] == "yes"
        ],
        "new_priority_citation_records": [
            {
                "citation_screening_id": row["citation_screening_id"],
                "tier": row["high_specificity_tier"],
                "title": row["title"],
                "doi": row["doi"],
            }
            for row in new_records
        ],
        "interpretation_boundary": (
            "Only genuinely new Tier 1/Tier 2 citation-neighbour records require "
            "additional manual classification before citation-chaining closure."
        ),
    }

    args.out_json.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
