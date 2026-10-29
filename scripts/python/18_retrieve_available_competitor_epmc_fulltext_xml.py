#!/usr/bin/env python3
"""Retrieve official Europe PMC full-text XML for competitor-verification records with PMCIDs."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import urllib.error
import urllib.request
from pathlib import Path


def retrieve_xml(pmcid: str) -> tuple[str, bytes | None, str]:
    url = f"https://www.ebi.ac.uk/europepmc/webservices/rest/{pmcid}/fullTextXML"
    request = urllib.request.Request(
        url,
        headers={"User-Agent": "Project6_NFkB_competitor_audit/1.0"},
    )
    try:
        with urllib.request.urlopen(request, timeout=180) as response:
            content = response.read()
            content_type = response.headers.get("Content-Type", "")
        if content.lstrip().startswith(b"<"):
            return "FULLTEXT_XML_RETRIEVED", content, content_type
        return "ENDPOINT_RETURNED_NON_XML_CONTENT", content, content_type
    except urllib.error.HTTPError as exc:
        return f"HTTP_ERROR_{exc.code}", None, ""
    except urllib.error.URLError as exc:
        return f"URL_ERROR_{exc.reason}", None, ""


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--access-table", required=True, type=Path)
    parser.add_argument("--raw-dir", required=True, type=Path)
    parser.add_argument("--out-table", required=True, type=Path)
    parser.add_argument("--out-json", required=True, type=Path)
    parser.add_argument("--out-checksums", required=True, type=Path)
    args = parser.parse_args()

    with args.access_table.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))

    args.raw_dir.mkdir(parents=True, exist_ok=True)
    args.out_table.parent.mkdir(parents=True, exist_ok=True)
    args.out_json.parent.mkdir(parents=True, exist_ok=True)
    args.out_checksums.parent.mkdir(parents=True, exist_ok=True)

    output = []
    checksum_lines = []

    for row in rows:
        pmcid = row["pmcid"].strip()
        if not pmcid:
            output.append({
                "screening_id": row["screening_id"],
                "classification": row["classification"],
                "doi": row["doi"],
                "pmcid": "",
                "xml_retrieval_status": "NO_PMCID_ABSTRACT_METADATA_ONLY",
                "xml_bytes": "0",
                "sha256": "",
                "local_raw_xml": "",
            })
            continue

        status, content, content_type = retrieve_xml(pmcid)
        local_path = args.raw_dir / f'{row["screening_id"]}_{pmcid}_fullTextXML.xml'
        sha256 = ""

        if content is not None and status == "FULLTEXT_XML_RETRIEVED":
            local_path.write_bytes(content)
            sha256 = hashlib.sha256(content).hexdigest()
            checksum_lines.append(f"{sha256}  {local_path}")
            byte_count = len(content)
        else:
            byte_count = len(content) if content is not None else 0

        output.append({
            "screening_id": row["screening_id"],
            "classification": row["classification"],
            "doi": row["doi"],
            "pmcid": pmcid,
            "xml_retrieval_status": status,
            "xml_bytes": str(byte_count),
            "sha256": sha256,
            "local_raw_xml": str(local_path) if sha256 else "",
        })
        print(f'{row["screening_id"]}\t{pmcid}\t{status}\t{byte_count}')

    fields = list(output[0].keys())
    with args.out_table.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fields, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(output)

    args.out_checksums.write_text(
        "\n".join(checksum_lines) + ("\n" if checksum_lines else ""),
        encoding="utf-8",
    )

    retrieved = [
        row["screening_id"] for row in output
        if row["xml_retrieval_status"] == "FULLTEXT_XML_RETRIEVED"
    ]
    metadata_only = [
        row["screening_id"] for row in output
        if row["xml_retrieval_status"] == "NO_PMCID_ABSTRACT_METADATA_ONLY"
    ]
    failures = [
        row for row in output
        if row["xml_retrieval_status"] not in {
            "FULLTEXT_XML_RETRIEVED",
            "NO_PMCID_ABSTRACT_METADATA_ONLY",
        }
    ]

    payload = {
        "status": (
            "FULLTEXT_XML_RETRIEVAL_COMPLETE"
            if not failures else "FULLTEXT_XML_RETRIEVAL_REVIEW_REQUIRED"
        ),
        "verification_subset_records": len(output),
        "fulltext_xml_retrieved_count": len(retrieved),
        "fulltext_xml_retrieved_screening_ids": retrieved,
        "abstract_metadata_only_screening_ids": metadata_only,
        "endpoint_failures": failures,
        "interpretation_boundary": (
            "Retrieved Europe PMC XML will be used for design-scope verification only. "
            "Records without retrievable XML remain assessed from indexed abstract and "
            "metadata unless additional legal full text is supplied."
        ),
    }
    args.out_json.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
