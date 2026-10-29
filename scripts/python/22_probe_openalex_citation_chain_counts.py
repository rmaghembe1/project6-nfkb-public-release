#!/usr/bin/env python3
"""Count backward-reference and forward-citation neighbours for Project 6 seed papers using OpenAlex."""

from __future__ import annotations

import argparse
import csv
import json
import os
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path


OPENALEX = "https://api.openalex.org/works"


def request_json(url: str, api_key: str) -> dict:
    separator = "&" if "?" in url else "?"
    if api_key:
        url = f"{url}{separator}api_key={urllib.parse.quote(api_key)}"
    req = urllib.request.Request(
        url,
        headers={"User-Agent": "Project6_NFkB_competitor_citation_chain/1.0"},
    )
    try:
        with urllib.request.urlopen(req, timeout=120) as response:
            return json.loads(response.read().decode("utf-8"))
    except urllib.error.HTTPError as exc:
        if exc.code in {401, 403} and not api_key:
            raise SystemExit(
                "ERROR: OpenAlex requires an API key for this request. "
                "Create a free OpenAlex API key, export OPENALEX_API_KEY, "
                "and rerun this count-only probe."
            )
        raise


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--seeds", required=True, type=Path)
    parser.add_argument("--raw-dir", required=True, type=Path)
    parser.add_argument("--out-table", required=True, type=Path)
    parser.add_argument("--out-json", required=True, type=Path)
    args = parser.parse_args()

    api_key = os.environ.get("OPENALEX_API_KEY", "").strip()

    with args.seeds.open(newline="", encoding="utf-8") as handle:
        seeds = list(csv.DictReader(handle, delimiter="\t"))

    args.raw_dir.mkdir(parents=True, exist_ok=True)
    args.out_table.parent.mkdir(parents=True, exist_ok=True)
    args.out_json.parent.mkdir(parents=True, exist_ok=True)

    output = []
    unresolved = []

    for seed in seeds:
        doi = seed["doi"].strip().lower()
        lookup = f"{OPENALEX}/doi:{urllib.parse.quote(doi, safe='/')}"
        payload = request_json(lookup, api_key)

        raw_path = args.raw_dir / f'{seed["seed_id"]}_openalex_seed_record.json'
        raw_path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")

        work_id_url = str(payload.get("id", ""))
        work_id = work_id_url.rsplit("/", 1)[-1] if work_id_url else ""
        if not work_id:
            unresolved.append(seed["seed_id"])
            continue

        referenced_works = payload.get("referenced_works", []) or []
        cited_by_count = int(payload.get("cited_by_count", 0) or 0)

        row = {
            "seed_id": seed["seed_id"],
            "source_screening_id": seed["source_screening_id"],
            "classification": seed["classification"],
            "title": seed["title"],
            "doi": doi,
            "openalex_work_id": work_id,
            "backward_reference_count": len(referenced_works),
            "forward_citation_count": cited_by_count,
            "total_first_degree_neighbour_count_before_deduplication": (
                len(referenced_works) + cited_by_count
            ),
        }
        output.append(row)
        print(
            f'{seed["seed_id"]}\t{work_id}\t'
            f'backward={len(referenced_works)}\tforward={cited_by_count}\t'
            f'total={len(referenced_works) + cited_by_count}'
        )

    fields = [
        "seed_id",
        "source_screening_id",
        "classification",
        "title",
        "doi",
        "openalex_work_id",
        "backward_reference_count",
        "forward_citation_count",
        "total_first_degree_neighbour_count_before_deduplication",
    ]
    with args.out_table.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fields, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(output)

    total_backward = sum(int(row["backward_reference_count"]) for row in output)
    total_forward = sum(int(row["forward_citation_count"]) for row in output)
    total_neighbours = total_backward + total_forward

    if unresolved:
        status = "SEED_RESOLUTION_INCOMPLETE_REVIEW_REQUIRED"
        retrieval_decision = "Do not retrieve neighbours until all seed works resolve."
    elif total_neighbours <= 800:
        status = "CITATION_CHAIN_COUNT_GATE_PASS"
        retrieval_decision = "Proceed to bounded neighbour-metadata retrieval and relevance screening."
    else:
        status = "CITATION_CHAIN_COUNT_GATE_REFINEMENT_REQUIRED"
        retrieval_decision = (
            "Citation network exceeds the prespecified 800-neighbour bound; "
            "restrict chaining to forward citations and selected backward seeds."
        )

    payload = {
        "status": status,
        "citation_source": "OpenAlex",
        "seed_record_count": len(seeds),
        "resolved_seed_record_count": len(output),
        "unresolved_seed_ids": unresolved,
        "total_backward_reference_links_before_deduplication": total_backward,
        "total_forward_citation_links_before_deduplication": total_forward,
        "total_first_degree_neighbour_links_before_deduplication": total_neighbours,
        "prespecified_full_neighbour_retrieval_threshold": 800,
        "retrieval_decision": retrieval_decision,
        "interpretation_boundary": (
            "Citation-link counts are not a novelty conclusion. "
            "Neighbour metadata must be retrieved, deduplicated and screened "
            "for direct-design overlap before the competitor audit can close."
        ),
    }
    args.out_json.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
