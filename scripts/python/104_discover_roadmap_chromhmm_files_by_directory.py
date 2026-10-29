#!/usr/bin/env python3
"""Discover Roadmap ChromHMM files from actual directory links.

This checkpoint discovers real Roadmap ChromHMM file links from directory listings
and maps priority immune/blood EIDs to observed filenames where possible.

It does not download full chromatin-state files, select final files, perform
candidate overlap, retrieve candidate chromatin-state evidence, validate function,
infer causality, infer selection, assign ancestry effects or authorise broad
biological interpretation.
"""

from __future__ import annotations

import csv
import json
import re
import urllib.error
import urllib.request
from html.parser import HTMLParser
from pathlib import Path
from urllib.parse import urljoin

PRIORITY_EIDS = Path("11_results/phaseD/regulatory_functional_evidence/query_pass1/roadmap_chromhmm/eid_refinement/phaseD_roadmap_immune_blood_eid_priority_for_file_selection_v1.tsv")
OUTDIR = Path("11_results/phaseD/regulatory_functional_evidence/query_pass1/roadmap_chromhmm/file_discovery")
OUTJSON = Path("results/qc/phaseD_roadmap_chromhmm_file_discovery_v1.json")

BASE_URL = "https://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/"
SEED_URLS = [
    BASE_URL,
    urljoin(BASE_URL, "imputed12marks/"),
]

OUTDIR.mkdir(parents=True, exist_ok=True)
OUTJSON.parent.mkdir(parents=True, exist_ok=True)

class LinkParser(HTMLParser):
    def __init__(self):
        super().__init__()
        self.links: list[str] = []

    def handle_starttag(self, tag, attrs):
        if tag.lower() == "a":
            href = dict(attrs).get("href")
            if href:
                self.links.append(href)

def fetch_links(url: str, timeout: int = 45) -> tuple[dict[str, str], list[str]]:
    req = urllib.request.Request(
        url,
        headers={"User-Agent": "Project6-PhaseD-Roadmap-ChromHMM-directory-discovery/1.0"},
    )
    try:
        with urllib.request.urlopen(req, timeout=timeout) as resp:
            text = resp.read(200000).decode("utf-8", errors="replace")
            parser = LinkParser()
            parser.feed(text)
            return {
                "url": url,
                "http_status": str(resp.status),
                "reachable_2xx": "YES" if 200 <= resp.status < 300 else "NO",
                "content_type": resp.headers.get("content-type", ""),
                "bytes_read": str(len(text.encode("utf-8"))),
                "error_type": "",
                "error_message": "",
            }, parser.links
    except urllib.error.HTTPError as exc:
        return {
            "url": url,
            "http_status": str(exc.code),
            "reachable_2xx": "NO",
            "content_type": exc.headers.get("content-type", "") if exc.headers else "",
            "bytes_read": "0",
            "error_type": "HTTPError",
            "error_message": str(exc),
        }, []
    except Exception as exc:
        return {
            "url": url,
            "http_status": "",
            "reachable_2xx": "NO",
            "content_type": "",
            "bytes_read": "0",
            "error_type": type(exc).__name__,
            "error_message": str(exc)[:300],
        }, []

with PRIORITY_EIDS.open("r", encoding="utf-8", newline="") as handle:
    priority_rows = list(csv.DictReader(handle, delimiter="\t"))

priority_eids = {r["eid"] for r in priority_rows}

visited = set()
to_visit = list(SEED_URLS)
route_rows = []
link_rows = []

# Shallow recursive directory discovery: enough to resolve model/subdirectory structure,
# but does not download full data files.
while to_visit and len(visited) < 50:
    url = to_visit.pop(0)
    if url in visited:
        continue
    visited.add(url)

    route_status, links = fetch_links(url)
    route_rows.append({
        **route_status,
        "links_found": str(len(links)),
        "candidate_overlap_performed": "NO",
        "candidate_chromatin_state_evidence_retrieved": "NO",
        "biological_interpretation_performed": "NO",
    })

    for href in links:
        if href.startswith("?") or href in {"../", "/"}:
            continue
        resolved = urljoin(url, href)
        name = href.rstrip("/").split("/")[-1]
        is_dir = "YES" if href.endswith("/") else "NO"
        is_bed_like = "YES" if re.search(r"\.(bed|bed\.gz|txt|gz)$", name.lower()) else "NO"
        eid_hits = sorted(eid for eid in priority_eids if eid in name or f"/{eid}" in resolved)

        link_rows.append({
            "source_url": url,
            "href": href,
            "resolved_url": resolved,
            "link_name": name,
            "is_directory": is_dir,
            "is_bed_like": is_bed_like,
            "priority_eid_hits": ";".join(eid_hits),
            "priority_eid_hit_count": str(len(eid_hits)),
            "candidate_file_for_priority_eid": "YES" if eid_hits and is_bed_like == "YES" else "NO",
            "specific_chromhmm_file_selected": "NO",
            "candidate_overlap_performed": "NO",
            "candidate_chromatin_state_evidence_retrieved": "NO",
            "biological_interpretation_performed": "NO",
        })

        # Follow directories under the Roadmap ChromHMM model tree only.
        if is_dir == "YES" and resolved.startswith(BASE_URL) and resolved not in visited and len(visited) + len(to_visit) < 50:
            to_visit.append(resolved)

candidate_links = [r for r in link_rows if r["candidate_file_for_priority_eid"] == "YES"]

mapped_rows = []
for prow in priority_rows:
    eid = prow["eid"]
    hits = [r for r in candidate_links if eid in r["priority_eid_hits"].split(";")]
    if hits:
        for hit in hits:
            mapped_rows.append({
                "eid": eid,
                "broad_group": prow.get("broad_group", ""),
                "mnemonic": prow.get("mnemonic", ""),
                "label": prow.get("label", ""),
                "harmonised_name": prow.get("harmonised_name", ""),
                "structured_category": prow.get("structured_category", ""),
                "sample_type": prow.get("sample_type", ""),
                "standard_status": prow.get("standard_status", ""),
                "observed_file_name": hit["link_name"],
                "observed_file_url": hit["resolved_url"],
                "apparent_build": "UNCONFIRMED",
                "file_mapping_status": "OBSERVED_CANDIDATE_FILE_NOT_FINAL",
                "specific_chromhmm_file_selected": "NO",
                "candidate_overlap_performed": "NO",
                "candidate_chromatin_state_evidence_retrieved": "NO",
                "biological_interpretation_performed": "NO",
            })
    else:
        mapped_rows.append({
            "eid": eid,
            "broad_group": prow.get("broad_group", ""),
            "mnemonic": prow.get("mnemonic", ""),
            "label": prow.get("label", ""),
            "harmonised_name": prow.get("harmonised_name", ""),
            "structured_category": prow.get("structured_category", ""),
            "sample_type": prow.get("sample_type", ""),
            "standard_status": prow.get("standard_status", ""),
            "observed_file_name": "",
            "observed_file_url": "",
            "apparent_build": "UNCONFIRMED",
            "file_mapping_status": "NO_OBSERVED_FILE_MATCH_IN_DISCOVERED_LINKS",
            "specific_chromhmm_file_selected": "NO",
            "candidate_overlap_performed": "NO",
            "candidate_chromatin_state_evidence_retrieved": "NO",
            "biological_interpretation_performed": "NO",
        })

route_path = OUTDIR / "phaseD_roadmap_chromhmm_directory_routes_discovered_v1.tsv"
links_path = OUTDIR / "phaseD_roadmap_chromhmm_directory_links_inventory_v1.tsv"
candidate_path = OUTDIR / "phaseD_roadmap_chromhmm_priority_eid_candidate_file_links_v1.tsv"
mapping_path = OUTDIR / "phaseD_roadmap_chromhmm_priority_eid_file_discovery_mapping_v1.tsv"
summary_path = OUTDIR / "phaseD_roadmap_chromhmm_file_discovery_summary_v1.tsv"

for path, rows in [
    (route_path, route_rows),
    (links_path, link_rows),
    (candidate_path, candidate_links),
    (mapping_path, mapped_rows),
]:
    with path.open("w", encoding="utf-8", newline="") as handle:
        fieldnames = list(rows[0].keys()) if rows else []
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)

mapped_eids = {r["eid"] for r in mapped_rows if r["file_mapping_status"] == "OBSERVED_CANDIDATE_FILE_NOT_FINAL"}

summary_rows = [
    {"metric": "priority_eids_tested", "value": len(priority_rows)},
    {"metric": "directory_routes_checked", "value": len(route_rows)},
    {"metric": "directory_routes_reachable_2xx", "value": sum(1 for r in route_rows if r["reachable_2xx"] == "YES")},
    {"metric": "links_inventory_rows", "value": len(link_rows)},
    {"metric": "priority_eid_candidate_file_links", "value": len(candidate_links)},
    {"metric": "priority_eids_with_observed_candidate_file", "value": len(mapped_eids)},
    {"metric": "priority_eids_without_observed_candidate_file", "value": len(priority_eids - mapped_eids)},
    {"metric": "specific_chromhmm_file_selected", "value": "NO"},
    {"metric": "build_confirmed_for_file", "value": "NO"},
    {"metric": "candidate_overlap_performed", "value": "NO"},
    {"metric": "candidate_chromatin_state_evidence_retrieved", "value": "NO"},
    {"metric": "biological_interpretation_performed", "value": "NO"},
    {"metric": "functional_validation_claim_authorised", "value": "NO"},
    {"metric": "broad_biological_interpretation_authorised", "value": "NO"},
]

with summary_path.open("w", encoding="utf-8", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=["metric", "value"], delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(summary_rows)

payload = {
    "status": "PASS",
    "operation": "Discover Roadmap ChromHMM files from actual directory links.",
    "priority_eids_tested": len(priority_rows),
    "directory_routes_checked": len(route_rows),
    "directory_routes_reachable_2xx": sum(1 for r in route_rows if r["reachable_2xx"] == "YES"),
    "links_inventory_rows": len(link_rows),
    "priority_eid_candidate_file_links": len(candidate_links),
    "priority_eids_with_observed_candidate_file": len(mapped_eids),
    "priority_eids_without_observed_candidate_file": len(priority_eids - mapped_eids),
    "specific_chromhmm_file_selected": "NO",
    "build_confirmed_for_file": "NO",
    "candidate_overlap_performed": "NO",
    "candidate_chromatin_state_evidence_retrieved": "NO",
    "biological_interpretation_performed": "NO",
    "functional_validation_claim_authorised": "NO",
    "broad_biological_interpretation_authorised": "NO",
    "outputs": {
        "routes": str(route_path),
        "links": str(links_path),
        "candidate_links": str(candidate_path),
        "mapping": str(mapping_path),
        "summary": str(summary_path),
    },
    "boundary": (
        "This checkpoint discovers actual Roadmap ChromHMM directory links only. It does not select final files, "
        "download full chromatin-state data, perform candidate overlap, retrieve candidate evidence, validate "
        "function or authorise biological interpretation."
    ),
}

OUTJSON.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
print(json.dumps(payload, indent=2))
