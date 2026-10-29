#!/usr/bin/env python3
"""Count PubMed and Europe PMC hits for formal competitor queries before retrieval."""

from __future__ import annotations

import argparse
import csv
import json
import os
import time
import urllib.parse
import urllib.request
from pathlib import Path


PUBMED_EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
EPMC_API = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"


def fetch_json(url: str) -> dict:
    request = urllib.request.Request(
        url,
        headers={"User-Agent": "Project6_NFkB_competitor_audit/1.0"},
    )
    with urllib.request.urlopen(request, timeout=120) as response:
        return json.loads(response.read().decode("utf-8"))


def classify_count(count: int) -> str:
    if count == 0:
        return "ZERO_RETRIEVAL_REVIEW_QUERY"
    if count <= 100:
        return "MANAGEABLE_FOR_FULL_RETRIEVAL"
    if count <= 300:
        return "RETRIEVABLE_WITH_CAUTION"
    return "TOO_BROAD_REFINE_BEFORE_RETRIEVAL"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--queries", required=True, type=Path)
    parser.add_argument("--out-table", required=True, type=Path)
    parser.add_argument("--out-json", required=True, type=Path)
    args = parser.parse_args()

    email = os.environ.get("NCBI_EMAIL", "").strip()
    if not email:
        raise SystemExit("ERROR: export NCBI_EMAIL before running.")
    tool = "Project6_NFkB_competitor_audit"

    with args.queries.open(newline="", encoding="utf-8") as handle:
        queries = list(csv.DictReader(handle, delimiter="\t"))

    output_rows = []

    for row in queries:
        qid = row["query_id"]
        purpose = row["screening_purpose"]

        pubmed_params = urllib.parse.urlencode({
            "db": "pubmed",
            "term": row["pubmed_query"],
            "retmode": "json",
            "retmax": "0",
            "tool": tool,
            "email": email,
        })
        pubmed_payload = fetch_json(f"{PUBMED_EUTILS}?{pubmed_params}")
        pubmed_count = int(pubmed_payload["esearchresult"]["count"])
        print(f"{qid}\tPubMed\t{pubmed_count}\t{classify_count(pubmed_count)}")
        output_rows.append({
            "query_id": qid,
            "source_database": "PubMed",
            "screening_purpose": purpose,
            "record_count": pubmed_count,
            "retrieval_decision": classify_count(pubmed_count),
            "query": row["pubmed_query"],
        })
        time.sleep(0.5)

        epmc_params = urllib.parse.urlencode({
            "query": row["europepmc_query"],
            "format": "json",
            "resultType": "lite",
            "pageSize": "1",
        })
        epmc_payload = fetch_json(f"{EPMC_API}?{epmc_params}")
        epmc_count = int(epmc_payload.get("hitCount", 0))
        print(f"{qid}\tEurope PMC\t{epmc_count}\t{classify_count(epmc_count)}")
        output_rows.append({
            "query_id": qid,
            "source_database": "Europe PMC",
            "screening_purpose": purpose,
            "record_count": epmc_count,
            "retrieval_decision": classify_count(epmc_count),
            "query": row["europepmc_query"],
        })
        time.sleep(0.5)

    args.out_table.parent.mkdir(parents=True, exist_ok=True)
    args.out_json.parent.mkdir(parents=True, exist_ok=True)

    fields = [
        "query_id",
        "source_database",
        "screening_purpose",
        "record_count",
        "retrieval_decision",
        "query",
    ]
    with args.out_table.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fields, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(output_rows)

    too_broad = [
        {
            "query_id": row["query_id"],
            "source_database": row["source_database"],
            "record_count": row["record_count"],
        }
        for row in output_rows
        if row["retrieval_decision"] == "TOO_BROAD_REFINE_BEFORE_RETRIEVAL"
    ]

    payload = {
        "status": "COUNT_GATE_COMPLETE_RETRIEVAL_NOT_YET_PERMITTED",
        "search_sources": ["PubMed", "Europe PMC"],
        "query_count": len(queries),
        "source_query_combinations_counted": len(output_rows),
        "threshold_policy": {
            "0": "review query sensitivity",
            "1-100": "manageable for full retrieval",
            "101-300": "retrieve only with caution and explicit cap",
            ">300": "refine query before full retrieval",
        },
        "too_broad_source_query_combinations": too_broad,
        "full_retrieval_permitted_without_revision": len(too_broad) == 0,
        "interpretation_boundary": (
            "Count-only reconnaissance is not evidence of novelty or redundancy. "
            "Full candidate retrieval and manual classification remain required."
        ),
    }

    args.out_json.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
