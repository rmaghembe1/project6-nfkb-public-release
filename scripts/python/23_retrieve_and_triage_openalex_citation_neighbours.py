#!/usr/bin/env python3
"""Retrieve and triage first-degree OpenAlex citation neighbours for Project 6 seeds."""

from __future__ import annotations

import argparse
import csv
import json
import os
import re
import time
import urllib.parse
import urllib.request
from collections import defaultdict
from pathlib import Path


OPENALEX_WORKS = "https://api.openalex.org/works"
MAX_OBSERVED_LINKS = 1000
PER_PAGE = 200

NFKB_TERMS = [
    "nf-kb", "nf-κb", "nf kappa b", "nuclear factor kappa b",
    "nfkb1", "nfkb2", "rela", "relb", "rel", "ikbkg", "ikk",
    "chuk", "ikbkb", "nfkbia", "tnfaip3", "a20", "tnip1", "cyld", "otulin",
]

DESIGN_TERMS = [
    "population", "ancestry", "ancestral", "evolution", "evolutionary",
    "selection", "adaptation", "adaptive", "haplotype", "introgress",
    "denisovan", "archaic", "polymorphism", "variant", "genetic interaction",
    "genetic association", "genomics", "genome", "eqtl", "regulatory",
    "copy number", "cnv", "structural variant",
]


def fetch_json(url: str, api_key: str) -> dict:
    if api_key:
        delimiter = "&" if "?" in url else "?"
        url = f"{url}{delimiter}api_key={urllib.parse.quote(api_key)}"
    request = urllib.request.Request(
        url,
        headers={"User-Agent": "Project6_NFkB_citation_chain_audit/1.0"},
    )
    with urllib.request.urlopen(request, timeout=180) as response:
        return json.loads(response.read().decode("utf-8"))


def abstract_from_inverted_index(index: dict | None) -> str:
    if not index:
        return ""
    positions = []
    for word, locations in index.items():
        for location in locations:
            positions.append((location, word))
    return " ".join(word for _, word in sorted(positions))


def normalise_doi(doi: str | None) -> str:
    if not doi:
        return ""
    return (
        doi.lower()
        .replace("https://doi.org/", "")
        .replace("http://doi.org/", "")
        .strip()
    )


def source_name(work: dict) -> str:
    location = work.get("primary_location") or {}
    source = location.get("source") or {}
    return str(source.get("display_name", ""))


def retrieve_direction(
    seed_id: str,
    work_id: str,
    direction: str,
    raw_dir: Path,
    api_key: str,
) -> list[dict]:
    if direction == "backward_reference":
        # Works cited by the seed paper.
        filter_value = f"cited_by:{work_id}"
    elif direction == "forward_citation":
        # Works that cite the seed paper.
        filter_value = f"cites:{work_id}"
    else:
        raise ValueError(direction)

    records = []
    cursor = "*"
    page = 0

    while True:
        page += 1
        params = urllib.parse.urlencode({
            "filter": filter_value,
            "per-page": str(PER_PAGE),
            "cursor": cursor,
        })
        payload = fetch_json(f"{OPENALEX_WORKS}?{params}", api_key)

        raw_path = raw_dir / f"{seed_id}_{direction}_page_{page:03d}.json"
        raw_path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")

        page_records = payload.get("results", []) or []
        records.extend(page_records)

        print(
            f"{seed_id}\t{direction}\tpage={page}\t"
            f"page_records={len(page_records)}\tcumulative={len(records)}"
        )

        next_cursor = (payload.get("meta") or {}).get("next_cursor")
        if not page_records or not next_cursor:
            break
        cursor = next_cursor
        time.sleep(0.2)

    return records


def is_plausible_screening_candidate(title: str, abstract: str) -> tuple[str, str]:
    text = f"{title} {abstract}".lower()
    matched_nfkb = sorted({term for term in NFKB_TERMS if term in text})
    matched_design = sorted({term for term in DESIGN_TERMS if term in text})

    if matched_nfkb and matched_design:
        return (
            "yes",
            f"NF-kB terms={';'.join(matched_nfkb)} | design terms={';'.join(matched_design)}"
        )
    return ("no", "")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--count-table", required=True, type=Path)
    parser.add_argument("--raw-dir", required=True, type=Path)
    parser.add_argument("--out-links", required=True, type=Path)
    parser.add_argument("--out-unique", required=True, type=Path)
    parser.add_argument("--out-plausible", required=True, type=Path)
    parser.add_argument("--out-json", required=True, type=Path)
    args = parser.parse_args()

    api_key = os.environ.get("OPENALEX_API_KEY", "").strip()

    with args.count_table.open(newline="", encoding="utf-8") as handle:
        seeds = list(csv.DictReader(handle, delimiter="\t"))

    expected_links = sum(
        int(row["total_first_degree_neighbour_count_before_deduplication"])
        for row in seeds
    )
    if expected_links > 800:
        raise SystemExit(
            f"ERROR: count-gated citation network has {expected_links} links, "
            "which exceeds the approved threshold of 800."
        )

    args.raw_dir.mkdir(parents=True, exist_ok=True)
    args.out_links.parent.mkdir(parents=True, exist_ok=True)
    args.out_unique.parent.mkdir(parents=True, exist_ok=True)
    args.out_plausible.parent.mkdir(parents=True, exist_ok=True)
    args.out_json.parent.mkdir(parents=True, exist_ok=True)

    link_rows = []
    work_records: dict[str, dict] = {}
    observed_by_seed_direction = []

    for seed in seeds:
        for direction in ["backward_reference", "forward_citation"]:
            records = retrieve_direction(
                seed["seed_id"],
                seed["openalex_work_id"],
                direction,
                args.raw_dir,
                api_key,
            )
            observed_by_seed_direction.append({
                "seed_id": seed["seed_id"],
                "direction": direction,
                "retrieved_link_count": len(records),
            })

            for work in records:
                work_id_url = str(work.get("id", ""))
                work_id = work_id_url.rsplit("/", 1)[-1] if work_id_url else ""
                if not work_id:
                    continue

                work_records[work_id] = work
                link_rows.append({
                    "seed_id": seed["seed_id"],
                    "seed_source_screening_id": seed["source_screening_id"],
                    "seed_classification": seed["classification"],
                    "relation_to_seed": direction,
                    "neighbour_openalex_work_id": work_id,
                    "neighbour_doi": normalise_doi(work.get("doi")),
                    "neighbour_title": str(work.get("display_name", "")),
                    "neighbour_year": str(work.get("publication_year", "")),
                    "neighbour_source": source_name(work),
                })

                if len(link_rows) > MAX_OBSERVED_LINKS:
                    raise SystemExit(
                        "ERROR: observed citation links exceed the retrieval safeguard "
                        f"of {MAX_OBSERVED_LINKS}; stop and refine the citation-chain scope."
                    )

            time.sleep(0.2)

    unique_context: dict[str, dict[str, set[str]]] = defaultdict(
        lambda: {"seed_ids": set(), "seed_screening_ids": set(), "directions": set()}
    )
    for row in link_rows:
        work_id = row["neighbour_openalex_work_id"]
        unique_context[work_id]["seed_ids"].add(row["seed_id"])
        unique_context[work_id]["seed_screening_ids"].add(
            row["seed_source_screening_id"]
        )
        unique_context[work_id]["directions"].add(row["relation_to_seed"])

    link_fields = [
        "seed_id", "seed_source_screening_id", "seed_classification",
        "relation_to_seed", "neighbour_openalex_work_id", "neighbour_doi",
        "neighbour_title", "neighbour_year", "neighbour_source",
    ]
    with args.out_links.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=link_fields, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(link_rows)

    unique_rows = []
    plausible_rows = []

    for work_id, context in unique_context.items():
        work = work_records[work_id]
        title = str(work.get("display_name", ""))
        abstract = abstract_from_inverted_index(work.get("abstract_inverted_index"))
        plausible, reason = is_plausible_screening_candidate(title, abstract)

        row = {
            "neighbour_openalex_work_id": work_id,
            "doi": normalise_doi(work.get("doi")),
            "title": title,
            "publication_year": str(work.get("publication_year", "")),
            "source": source_name(work),
            "work_type": str(work.get("type", "")),
            "cited_by_count": str(work.get("cited_by_count", "")),
            "linked_seed_ids": ";".join(sorted(context["seed_ids"])),
            "linked_seed_screening_ids": ";".join(
                sorted(context["seed_screening_ids"])
            ),
            "relations_to_seeds": ";".join(sorted(context["directions"])),
            "plausible_screening_candidate": plausible,
            "triage_reason": reason,
            "abstract": abstract,
            "screening_classification": "",
            "screening_reason": "",
            "direct_competitor": "",
        }
        unique_rows.append(row)
        if plausible == "yes":
            plausible_rows.append(row)

    unique_rows.sort(
        key=lambda row: (row["publication_year"], row["title"]),
        reverse=True,
    )
    plausible_rows.sort(
        key=lambda row: (row["publication_year"], row["title"]),
        reverse=True,
    )

    unique_fields = list(unique_rows[0].keys()) if unique_rows else [
        "neighbour_openalex_work_id", "doi", "title", "publication_year",
        "source", "work_type", "cited_by_count", "linked_seed_ids",
        "linked_seed_screening_ids", "relations_to_seeds",
        "plausible_screening_candidate", "triage_reason", "abstract",
        "screening_classification", "screening_reason", "direct_competitor",
    ]

    for output_path, rows in [
        (args.out_unique, unique_rows),
        (args.out_plausible, plausible_rows),
    ]:
        with output_path.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(
                handle, fieldnames=unique_fields, delimiter="\t", lineterminator="\n"
            )
            writer.writeheader()
            writer.writerows(rows)

    payload = {
        "status": "CITATION_NEIGHBOUR_RETRIEVAL_COMPLETE_SCREENING_PENDING",
        "citation_source": "OpenAlex",
        "approved_count_gate_links_before_deduplication": expected_links,
        "observed_links_retrieved_before_deduplication": len(link_rows),
        "maximum_observed_link_safeguard": MAX_OBSERVED_LINKS,
        "unique_neighbour_records_after_openalex_id_deduplication": len(unique_rows),
        "plausible_nfkb_design_overlap_candidates_for_screening": len(plausible_rows),
        "observed_counts_by_seed_and_direction": observed_by_seed_direction,
        "interpretation_boundary": (
            "Automated term triage identifies records for manual screening only. "
            "The competitor-search requirement remains open until plausible citation "
            "neighbours are classified and checked against already screened records."
        ),
    }

    args.out_json.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
