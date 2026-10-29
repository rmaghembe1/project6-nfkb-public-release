#!/usr/bin/env python3
"""Retrieve bounded PubMed and Europe PMC competitor candidates after count-gate approval."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import time
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from collections import defaultdict
from pathlib import Path


PUBMED_ESEARCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
PUBMED_EFETCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
EPMC_SEARCH = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"
MAX_RECORDS_PER_SOURCE_QUERY = 100


def fetch_bytes(url: str) -> bytes:
    request = urllib.request.Request(
        url,
        headers={"User-Agent": "Project6_NFkB_competitor_audit/1.0"},
    )
    with urllib.request.urlopen(request, timeout=180) as response:
        return response.read()


def write_raw(path: Path, content: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(content)


def xml_text(node: ET.Element | None) -> str:
    if node is None:
        return ""
    return "".join(node.itertext()).strip()


def stable_key(record: dict[str, str]) -> str:
    doi = record["doi"].strip().lower()
    pmid = record["pmid"].strip()
    title = " ".join(record["title"].lower().split())
    if doi:
        return f"doi:{doi}"
    if pmid:
        return f"pmid:{pmid}"
    return f"title_sha256:{hashlib.sha256(title.encode('utf-8')).hexdigest()}"


def read_count_gate(path: Path) -> dict[tuple[str, str], int]:
    counts: dict[tuple[str, str], int] = {}
    with path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            count = int(row["record_count"])
            if count > MAX_RECORDS_PER_SOURCE_QUERY:
                raise SystemExit(
                    "ERROR: bounded retrieval blocked because "
                    f'{row["query_id"]} / {row["source_database"]} has '
                    f"{count} records, exceeding the cap of "
                    f"{MAX_RECORDS_PER_SOURCE_QUERY}."
                )
            counts[(row["query_id"], row["source_database"])] = count
    return counts


def retrieve_pubmed(
    query_id: str,
    query: str,
    expected_count: int,
    raw_dir: Path,
    email: str,
) -> list[dict[str, str]]:
    if expected_count == 0:
        print(f"{query_id}\tPubMed\tSKIP_ZERO_HITS")
        return []

    params = urllib.parse.urlencode({
        "db": "pubmed",
        "term": query,
        "retmode": "json",
        "retmax": str(MAX_RECORDS_PER_SOURCE_QUERY),
        "tool": "Project6_NFkB_competitor_audit",
        "email": email,
    })
    raw_search = fetch_bytes(f"{PUBMED_ESEARCH}?{params}")
    write_raw(raw_dir / f"{query_id}_pubmed_esearch.json", raw_search)
    search = json.loads(raw_search.decode("utf-8"))
    observed_count = int(search["esearchresult"]["count"])
    ids = search["esearchresult"].get("idlist", [])

    if observed_count != expected_count:
        raise SystemExit(
            f"ERROR: PubMed count changed for {query_id}: "
            f"count gate={expected_count}, retrieval={observed_count}. "
            "Re-run the count gate before retrieval."
        )
    if len(ids) != expected_count:
        raise SystemExit(
            f"ERROR: PubMed ID retrieval mismatch for {query_id}: "
            f"expected {expected_count}, retrieved {len(ids)} IDs."
        )

    print(f"{query_id}\tPubMed\tRETRIEVING\t{observed_count}")

    records: list[dict[str, str]] = []
    for batch_number, start in enumerate(range(0, len(ids), 100), start=1):
        batch = ids[start:start + 100]
        fetch_params = urllib.parse.urlencode({
            "db": "pubmed",
            "id": ",".join(batch),
            "retmode": "xml",
            "tool": "Project6_NFkB_competitor_audit",
            "email": email,
        })
        raw_fetch = fetch_bytes(f"{PUBMED_EFETCH}?{fetch_params}")
        write_raw(
            raw_dir / f"{query_id}_pubmed_efetch_{batch_number:03d}.xml",
            raw_fetch,
        )
        root = ET.fromstring(raw_fetch)

        for article in root.findall(".//PubmedArticle"):
            pmid = xml_text(article.find(".//PMID"))
            title = xml_text(article.find(".//ArticleTitle"))
            abstract = " ".join(
                xml_text(item)
                for item in article.findall(".//Abstract/AbstractText")
            ).strip()
            journal = xml_text(article.find(".//Journal/Title"))
            year = (
                xml_text(article.find(".//PubDate/Year"))
                or xml_text(article.find(".//ArticleDate/Year"))
                or xml_text(article.find(".//PubDate/MedlineDate"))[:4]
            )
            doi = ""
            for aid in article.findall(".//ArticleId"):
                if aid.attrib.get("IdType") == "doi":
                    doi = xml_text(aid)
                    break
            records.append({
                "query_id": query_id,
                "source_database": "PubMed",
                "pmid": pmid,
                "pmcid": "",
                "doi": doi,
                "title": title,
                "abstract": abstract,
                "journal": journal,
                "year": year,
                "publication_type": "",
            })
        time.sleep(0.4)

    if len(records) != expected_count:
        raise SystemExit(
            f"ERROR: PubMed parsed-record mismatch for {query_id}: "
            f"expected {expected_count}, parsed {len(records)}."
        )
    return records


def retrieve_epmc(
    query_id: str,
    query: str,
    expected_count: int,
    raw_dir: Path,
) -> list[dict[str, str]]:
    if expected_count == 0:
        print(f"{query_id}\tEurope PMC\tSKIP_ZERO_HITS")
        return []

    params = urllib.parse.urlencode({
        "query": query,
        "format": "json",
        "resultType": "core",
        "pageSize": str(MAX_RECORDS_PER_SOURCE_QUERY),
    })
    raw = fetch_bytes(f"{EPMC_SEARCH}?{params}")
    write_raw(raw_dir / f"{query_id}_epmc_page_001.json", raw)
    payload = json.loads(raw.decode("utf-8"))

    observed_count = int(payload.get("hitCount", 0))
    if observed_count != expected_count:
        raise SystemExit(
            f"ERROR: Europe PMC count changed for {query_id}: "
            f"count gate={expected_count}, retrieval={observed_count}. "
            "Re-run the count gate before retrieval."
        )
    if observed_count > MAX_RECORDS_PER_SOURCE_QUERY:
        raise SystemExit(
            f"ERROR: Europe PMC retrieval cap exceeded for {query_id}: "
            f"{observed_count} records."
        )

    print(f"{query_id}\tEurope PMC\tRETRIEVING\t{observed_count}")

    records: list[dict[str, str]] = []
    for item in payload.get("resultList", {}).get("result", []):
        records.append({
            "query_id": query_id,
            "source_database": "Europe PMC",
            "pmid": str(item.get("pmid", "")),
            "pmcid": str(item.get("pmcid", "")),
            "doi": str(item.get("doi", "")),
            "title": str(item.get("title", "")),
            "abstract": str(item.get("abstractText", "")),
            "journal": str(item.get("journalTitle", "")),
            "year": str(item.get("pubYear", "")),
            "publication_type": str(item.get("pubType", "")),
        })

    if len(records) != expected_count:
        raise SystemExit(
            f"ERROR: Europe PMC parsed-record mismatch for {query_id}: "
            f"expected {expected_count}, parsed {len(records)}."
        )
    return records


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--queries", required=True, type=Path)
    parser.add_argument("--count-table", required=True, type=Path)
    parser.add_argument("--raw-dir", required=True, type=Path)
    parser.add_argument("--out-all", required=True, type=Path)
    parser.add_argument("--out-unique", required=True, type=Path)
    parser.add_argument("--out-summary", required=True, type=Path)
    args = parser.parse_args()

    email = os.environ.get("NCBI_EMAIL", "").strip()
    if not email:
        raise SystemExit("ERROR: export NCBI_EMAIL before running.")

    counts = read_count_gate(args.count_table)

    with args.queries.open(newline="", encoding="utf-8") as handle:
        queries = list(csv.DictReader(handle, delimiter="\t"))

    args.raw_dir.mkdir(parents=True, exist_ok=True)
    records: list[dict[str, str]] = []
    per_query_counts: list[dict[str, str | int]] = []

    for row in queries:
        query_id = row["query_id"]

        pubmed_expected = counts[(query_id, "PubMed")]
        pubmed = retrieve_pubmed(
            query_id,
            row["pubmed_query"],
            pubmed_expected,
            args.raw_dir,
            email,
        )
        records.extend(pubmed)
        per_query_counts.append({
            "query_id": query_id,
            "source_database": "PubMed",
            "retrieved_records": len(pubmed),
        })
        time.sleep(0.5)

        epmc_expected = counts[(query_id, "Europe PMC")]
        epmc = retrieve_epmc(
            query_id,
            row["europepmc_query"],
            epmc_expected,
            args.raw_dir,
        )
        records.extend(epmc)
        per_query_counts.append({
            "query_id": query_id,
            "source_database": "Europe PMC",
            "retrieved_records": len(epmc),
        })
        time.sleep(0.5)

    fields = [
        "query_id",
        "source_database",
        "pmid",
        "pmcid",
        "doi",
        "title",
        "abstract",
        "journal",
        "year",
        "publication_type",
    ]
    args.out_all.parent.mkdir(parents=True, exist_ok=True)
    with args.out_all.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fields, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(records)

    grouped: dict[str, list[dict[str, str]]] = defaultdict(list)
    for record in records:
        grouped[stable_key(record)].append(record)

    unique_rows: list[dict[str, str]] = []
    for key, grouped_records in grouped.items():
        representative = max(
            grouped_records,
            key=lambda row: len(row.get("abstract", "")),
        )
        unique_rows.append({
            "record_key": key,
            "query_ids": ";".join(sorted({row["query_id"] for row in grouped_records})),
            "source_databases": ";".join(
                sorted({row["source_database"] for row in grouped_records})
            ),
            "pmid": representative["pmid"],
            "pmcid": representative["pmcid"],
            "doi": representative["doi"],
            "title": representative["title"],
            "abstract": representative["abstract"],
            "journal": representative["journal"],
            "year": representative["year"],
            "publication_type": representative["publication_type"],
            "screening_classification": "",
            "screening_reason": "",
        })

    unique_rows.sort(
        key=lambda row: (row["year"], row["title"]),
        reverse=True,
    )
    unique_fields = [
        "record_key",
        "query_ids",
        "source_databases",
        "pmid",
        "pmcid",
        "doi",
        "title",
        "abstract",
        "journal",
        "year",
        "publication_type",
        "screening_classification",
        "screening_reason",
    ]
    with args.out_unique.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=unique_fields, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(unique_rows)

    positive_controls = {}
    for qid in sorted({row["query_id"] for row in queries}):
        if qid.endswith("07A") or qid.endswith("07B"):
            positive_controls[qid] = {
                "pubmed_retrieved": counts[(qid, "PubMed")] == 1,
                "europepmc_retrieved": counts[(qid, "Europe PMC")] == 1,
            }

    payload = {
        "status": "BOUNDED_RETRIEVAL_COMPLETE_SCREENING_PENDING",
        "query_register": str(args.queries),
        "count_gate_table": str(args.count_table),
        "maximum_allowed_records_per_source_query": MAX_RECORDS_PER_SOURCE_QUERY,
        "total_source_query_records_retrieved": len(records),
        "unique_candidate_records_after_deduplication": len(unique_rows),
        "per_query_source_retrieval_counts": per_query_counts,
        "positive_control_retrieval": positive_controls,
        "interpretation_boundary": (
            "Bounded retrieval is not a novelty conclusion. Unique candidates "
            "must be screened and classified before direct-competitor or "
            "distinctness conclusions are recorded."
        ),
    }
    args.out_summary.parent.mkdir(parents=True, exist_ok=True)
    args.out_summary.write_text(
        json.dumps(payload, indent=2) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
