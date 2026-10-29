#!/usr/bin/env python3
"""Probe Europe PMC access status for the Project 6 full-text verification subset."""

from __future__ import annotations

import argparse
import csv
import json
import urllib.parse
import urllib.request
from pathlib import Path


EPMC_SEARCH = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"


def fetch_json(url: str) -> dict:
    req = urllib.request.Request(
        url,
        headers={"User-Agent": "Project6_NFkB_competitor_audit/1.0"},
    )
    with urllib.request.urlopen(req, timeout=120) as response:
        return json.loads(response.read().decode("utf-8"))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--register", required=True, type=Path)
    parser.add_argument("--out-table", required=True, type=Path)
    parser.add_argument("--out-json", required=True, type=Path)
    args = parser.parse_args()

    with args.register.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))

    output = []
    for row in rows:
        doi = row["doi"]
        query = f'DOI:"{doi}"'
        params = urllib.parse.urlencode({
            "query": query,
            "format": "json",
            "resultType": "core",
            "pageSize": "5",
        })
        payload = fetch_json(f"{EPMC_SEARCH}?{params}")
        matches = payload.get("resultList", {}).get("result", [])

        exact = [
            item for item in matches
            if str(item.get("doi", "")).lower() == doi.lower()
        ]
        item = exact[0] if exact else (matches[0] if matches else {})

        pmcid = str(item.get("pmcid", ""))
        is_open_access = str(item.get("isOpenAccess", "")).lower()
        has_pdf = str(item.get("hasPDF", "")).lower()
        has_full_text_xml = str(item.get("hasTextMinedTerms", "")).lower()

        if pmcid and is_open_access == "y":
            access_status = "OPEN_ACCESS_PMC_RECORD_AVAILABLE"
        elif pmcid:
            access_status = "PMC_RECORD_PRESENT_ACCESS_REVIEW_REQUIRED"
        elif item:
            access_status = "INDEXED_NO_PMC_FULLTEXT_IDENTIFIED"
        else:
            access_status = "NO_EXACT_EPMC_DOI_MATCH"

        output.append({
            "screening_id": row["screening_id"],
            "classification": row["classification"],
            "doi": doi,
            "title": row["title"],
            "epmc_hit_count": payload.get("hitCount", 0),
            "pmid": str(item.get("pmid", "")),
            "pmcid": pmcid,
            "is_open_access": is_open_access,
            "has_pdf": has_pdf,
            "access_status": access_status,
            "fulltext_xml_endpoint_if_pmc": (
                f"https://www.ebi.ac.uk/europepmc/webservices/rest/{pmcid}/fullTextXML"
                if pmcid else ""
            ),
        })

    args.out_table.parent.mkdir(parents=True, exist_ok=True)
    args.out_json.parent.mkdir(parents=True, exist_ok=True)

    fields = list(output[0].keys())
    with args.out_table.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fields, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(output)

    counts = {}
    for row in output:
        counts[row["access_status"]] = counts.get(row["access_status"], 0) + 1

    summary = {
        "status": "FULLTEXT_ACCESS_PROBE_COMPLETE",
        "queried_records": len(output),
        "access_status_counts": counts,
        "records": output,
        "interpretation_boundary": (
            "Open-access probing establishes retrievability only. "
            "It does not itself perform full-text novelty verification."
        ),
    }
    args.out_json.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
