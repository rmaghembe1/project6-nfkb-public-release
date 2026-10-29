#!/usr/bin/env python3
"""Resolve bibliographic metadata for selected Project 6 novelty-boundary papers using Crossref."""

from __future__ import annotations

import csv
import json
import os
import urllib.parse
import urllib.request
from pathlib import Path


SELECTED = [
    {
        "id": "L08",
        "doi": "10.1371/journal.pgen.1012067",
        "role_in_project": "Archaic regulatory-introgression boundary including predicted TNFAIP3 effects",
        "phaseA_assessment": "Partial overlap verified; not NF-kB circuitry-wide",
    },
    {
        "id": "L09",
        "doi": "10.1186/1471-2148-13-19",
        "role_in_project": "Population-differentiated CNV pathway-evolution boundary implicating NF-kB",
        "phaseA_assessment": "Partial overlap verified by official full-text XML; distinct variant class and design",
    },
    {
        "id": "L10",
        "doi": "10.2217/pme.09.54",
        "role_in_project": "Broad human population-pathway evolutionary boundary implicating NF-kB",
        "phaseA_assessment": "Near competitor boundary verified from abstract/metadata; official XML endpoint returned HTTP 404",
    },
    {
        "id": "L11",
        "doi": "10.4049/jimmunol.1000324",
        "role_in_project": "African-derived TNFAIP3 population-functional boundary",
        "phaseA_assessment": "Partial overlap verified; single-locus rather than circuitry-wide",
    },
    {
        "id": "L12",
        "doi": "10.1016/j.molimm.2013.07.006",
        "role_in_project": "Multi-gene human NF-kB clinical-genetic claim boundary",
        "phaseA_assessment": "Partial overlap verified; prohibits generic first-multi-gene claim",
    },
]


def fetch_crossref(doi: str) -> dict:
    url = "https://api.crossref.org/works/" + urllib.parse.quote(doi, safe="")
    contact_email = os.environ.get("PROJECT6_CROSSREF_CONTACT_EMAIL", "").strip()
    if (not contact_email or "@" not in contact_email
            or any(ch in contact_email for ch in ("\r", "\n", " "))):
        raise ValueError("Set PROJECT6_CROSSREF_CONTACT_EMAIL to a valid contact email for Crossref requests")
    request = urllib.request.Request(
        url,
        headers={"User-Agent": f"Project6_NFkB_manifest_metadata/1.0 (mailto:{contact_email})"},
    )
    with urllib.request.urlopen(request, timeout=120) as response:
        return json.loads(response.read().decode("utf-8"))["message"]


def author_string(message: dict) -> str:
    authors = message.get("author", []) or []
    if not authors:
        return ""
    first = authors[0]
    lead = " ".join(
        part for part in [first.get("given", ""), first.get("family", "")]
        if part
    ).strip()
    return f"{lead} et al." if len(authors) > 1 else lead


def get_year(message: dict) -> str:
    for field in ["published-print", "published-online", "published", "issued"]:
        value = message.get(field, {})
        parts = value.get("date-parts", [])
        if parts and parts[0]:
            return str(parts[0][0])
    return ""


def main() -> None:
    raw_dir = Path("02_reference_data/competitor_search_raw/manifest_metadata")
    out = Path("11_results/competitor_search/manifest_updates/Project6_literature_manifest_boundary_additions_v1.tsv")
    audit = Path("15_logs_and_audit/competitor_search/manifest_updates/Project6_boundary_metadata_resolution_v1.json")
    raw_dir.mkdir(parents=True, exist_ok=True)
    out.parent.mkdir(parents=True, exist_ok=True)
    audit.parent.mkdir(parents=True, exist_ok=True)

    rows = []
    audit_rows = []

    for selected in SELECTED:
        message = fetch_crossref(selected["doi"])
        raw_path = raw_dir / f'{selected["id"]}_crossref.json'
        raw_path.write_text(json.dumps(message, indent=2) + "\n", encoding="utf-8")

        title_list = message.get("title", []) or []
        journal_list = message.get("container-title", []) or []

        row = {
            "id": selected["id"],
            "authors": author_string(message),
            "title": title_list[0] if title_list else "",
            "journal": journal_list[0] if journal_list else "",
            "year": get_year(message),
            "doi": selected["doi"],
            "role_in_project": selected["role_in_project"],
            "phaseA_assessment": selected["phaseA_assessment"],
        }
        rows.append(row)
        audit_rows.append({
            "id": selected["id"],
            "doi": selected["doi"],
            "resolved_title": row["title"],
            "resolved_journal": row["journal"],
            "resolved_year": row["year"],
            "resolved_authors": row["authors"],
        })

    fields = [
        "id", "authors", "title", "journal", "year", "doi",
        "role_in_project", "phaseA_assessment"
    ]
    with out.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fields, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)

    audit.write_text(
        json.dumps(
            {
                "status": "SELECTED_BOUNDARY_METADATA_RESOLVED_REVIEW_BEFORE_MANIFEST_APPEND",
                "metadata_source": "Crossref DOI metadata API",
                "selected_boundary_records": audit_rows,
            },
            indent=2,
        ) + "\n",
        encoding="utf-8",
    )

    print(out.read_text(encoding="utf-8"))
    print(audit.read_text(encoding="utf-8"))


if __name__ == "__main__":
    main()
