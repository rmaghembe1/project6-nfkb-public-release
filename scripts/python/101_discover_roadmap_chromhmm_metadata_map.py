#!/usr/bin/env python3
"""Discover Roadmap ChromHMM metadata/file-map routes.

This checkpoint searches Roadmap portal pages and likely metadata/download
routes to identify sample/cell/tissue labels and file naming patterns needed
before selecting immune/blood-relevant ChromHMM files.

It does not select a final file, download full chromatin-state datasets, perform
candidate overlap, retrieve candidate chromatin-state evidence, validate
function, infer causality, infer selection, assign ancestry effects or authorise
broad biological interpretation.
"""

from __future__ import annotations

import csv
import json
import re
import urllib.error
import urllib.request
from html.parser import HTMLParser
from pathlib import Path

OUTDIR = Path("11_results/phaseD/regulatory_functional_evidence/query_pass1/roadmap_chromhmm/metadata_map")
OUTJSON = Path("results/qc/phaseD_roadmap_chromhmm_metadata_map_discovery_v1.json")

OUTDIR.mkdir(parents=True, exist_ok=True)
OUTJSON.parent.mkdir(parents=True, exist_ok=True)

ROUTES = [
    {
        "route_id": "roadmap_metadata_portal",
        "url": "https://egg2.wustl.edu/roadmap/web_portal/meta.html",
        "route_type": "metadata_page",
    },
    {
        "route_id": "roadmap_metadata_csv_candidate",
        "url": "https://egg2.wustl.edu/roadmap/data/byFileType/metadata/metadata.csv",
        "route_type": "metadata_file_candidate",
    },
    {
        "route_id": "roadmap_jul2013_sample_metadata_archive",
        "url": "https://egg2.wustl.edu/roadmap/data/byFileType/metadata/July2013/sample_metadata.csv",
        "route_type": "metadata_file_candidate",
    },
    {
        "route_id": "roadmap_roadmap_epigenomics_metadata",
        "url": "https://egg2.wustl.edu/roadmap/data/byFileType/metadata/",
        "route_type": "metadata_directory",
    },
    {
        "route_id": "roadmap_chromhmm_imputed12marks_dir",
        "url": "https://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/imputed12marks/",
        "route_type": "chromhmm_directory",
    },
]

IMMUNE_BLOOD_TERMS = [
    "blood", "monocyte", "macrophage", "neutrophil", "t cell", "t-cell", "tcell",
    "b cell", "b-cell", "bcell", "nk", "cd4", "cd8", "pbmc", "lymph",
    "thymus", "spleen", "hematopoietic", "haematopoietic", "immune",
    "leukocyte", "granulocyte", "cd14", "cd19", "cd34", "cd56"
]

class LinkParser(HTMLParser):
    def __init__(self):
        super().__init__()
        self.links: list[str] = []

    def handle_starttag(self, tag, attrs):
        if tag.lower() == "a":
            href = dict(attrs).get("href")
            if href:
                self.links.append(href)

def fetch(url: str, timeout: int = 45) -> dict[str, str]:
    req = urllib.request.Request(
        url,
        headers={
            "User-Agent": "Project6-PhaseD-Roadmap-metadata-map-discovery/1.0",
            "Accept": "text/html,text/plain,text/csv,application/json,*/*",
        },
    )
    try:
        with urllib.request.urlopen(req, timeout=timeout) as resp:
            raw = resp.read(10000)
            text = raw.decode("utf-8", errors="replace")
            return {
                "http_status": str(resp.status),
                "reachable_2xx": "YES" if 200 <= resp.status < 300 else "NO",
                "content_type": resp.headers.get("content-type", ""),
                "bytes_read": str(len(raw)),
                "text_excerpt": text.replace("\n", " ")[:1000],
                "raw_text": text,
                "error_type": "",
                "error_message": "",
            }
    except urllib.error.HTTPError as exc:
        raw = b""
        try:
            raw = exc.read(3000)
        except Exception:
            pass
        text = raw.decode("utf-8", errors="replace")
        return {
            "http_status": str(exc.code),
            "reachable_2xx": "NO",
            "content_type": exc.headers.get("content-type", "") if exc.headers else "",
            "bytes_read": str(len(raw)),
            "text_excerpt": text.replace("\n", " ")[:1000],
            "raw_text": text,
            "error_type": "HTTPError",
            "error_message": str(exc),
        }
    except Exception as exc:
        return {
            "http_status": "",
            "reachable_2xx": "NO",
            "content_type": "",
            "bytes_read": "0",
            "text_excerpt": "",
            "raw_text": "",
            "error_type": type(exc).__name__,
            "error_message": str(exc)[:300],
        }

route_rows = []
link_rows = []
candidate_metadata_hits = []

for route in ROUTES:
    result = fetch(route["url"])
    raw_text = result.pop("raw_text")
    route_rows.append({
        **route,
        **result,
        "metadata_or_file_map_selected": "NO",
        "candidate_overlap_performed": "NO",
        "candidate_chromatin_state_evidence_retrieved": "NO",
        "biological_interpretation_performed": "NO",
    })

    parser = LinkParser()
    parser.feed(raw_text)
    base = route["url"].rsplit("/", 1)[0] + "/"
    for href in parser.links:
        href_l = href.lower()
        full = href if href.startswith("http") else base + href
        name = href.rstrip("/").split("/")[-1]
        immune_flag = "YES" if any(term in href_l for term in IMMUNE_BLOOD_TERMS) else "NO"
        metadata_like = "YES" if any(x in href_l for x in ["meta", "sample", "mnemonic", "chromhmm", "imputed", ".bed", ".txt", ".csv", ".xls"]) else "NO"
        link_rows.append({
            "source_route_id": route["route_id"],
            "href": href,
            "resolved_url": full,
            "link_name": name,
            "metadata_or_chromhmm_like": metadata_like,
            "immune_blood_relevant_name_flag": immune_flag,
            "candidate_overlap_performed": "NO",
            "candidate_chromatin_state_evidence_retrieved": "NO",
        })

    lower_text = raw_text.lower()
    immune_term_hits = sorted({term for term in IMMUNE_BLOOD_TERMS if term in lower_text})
    if immune_term_hits:
        candidate_metadata_hits.append({
            "route_id": route["route_id"],
            "url": route["url"],
            "immune_or_blood_terms_detected": ";".join(immune_term_hits),
            "metadata_map_selected": "NO",
            "candidate_overlap_performed": "NO",
        })

routes_path = OUTDIR / "phaseD_roadmap_chromhmm_metadata_route_diagnostics_v1.tsv"
with routes_path.open("w", encoding="utf-8", newline="") as handle:
    fieldnames = list(route_rows[0].keys())
    writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(route_rows)

links_path = OUTDIR / "phaseD_roadmap_chromhmm_metadata_links_inventory_v1.tsv"
with links_path.open("w", encoding="utf-8", newline="") as handle:
    fieldnames = list(link_rows[0].keys()) if link_rows else [
        "source_route_id", "href", "resolved_url", "link_name",
        "metadata_or_chromhmm_like", "immune_blood_relevant_name_flag",
        "candidate_overlap_performed", "candidate_chromatin_state_evidence_retrieved"
    ]
    writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(link_rows)

immune_links = [r for r in link_rows if r["immune_blood_relevant_name_flag"] == "YES" or r["metadata_or_chromhmm_like"] == "YES"]
immune_links_path = OUTDIR / "phaseD_roadmap_chromhmm_metadata_candidate_links_v1.tsv"
with immune_links_path.open("w", encoding="utf-8", newline="") as handle:
    fieldnames = list(link_rows[0].keys()) if link_rows else [
        "source_route_id", "href", "resolved_url", "link_name",
        "metadata_or_chromhmm_like", "immune_blood_relevant_name_flag",
        "candidate_overlap_performed", "candidate_chromatin_state_evidence_retrieved"
    ]
    writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(immune_links)

hits_path = OUTDIR / "phaseD_roadmap_chromhmm_metadata_immune_term_hits_v1.tsv"
with hits_path.open("w", encoding="utf-8", newline="") as handle:
    fieldnames = ["route_id", "url", "immune_or_blood_terms_detected", "metadata_map_selected", "candidate_overlap_performed"]
    writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(candidate_metadata_hits)

summary_path = OUTDIR / "phaseD_roadmap_chromhmm_metadata_map_discovery_summary_v1.tsv"
summary_rows = [
    {"metric": "routes_checked", "value": len(route_rows)},
    {"metric": "routes_reachable_2xx", "value": sum(1 for r in route_rows if r["reachable_2xx"] == "YES")},
    {"metric": "links_inventory_rows", "value": len(link_rows)},
    {"metric": "candidate_metadata_or_chromhmm_links", "value": len(immune_links)},
    {"metric": "routes_with_immune_blood_text_hits", "value": len(candidate_metadata_hits)},
    {"metric": "metadata_map_selected", "value": "NO"},
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
    "operation": "Discover Roadmap ChromHMM metadata/file-map routes before file selection.",
    "routes_checked": len(route_rows),
    "routes_reachable_2xx": sum(1 for r in route_rows if r["reachable_2xx"] == "YES"),
    "links_inventory_rows": len(link_rows),
    "candidate_metadata_or_chromhmm_links": len(immune_links),
    "routes_with_immune_blood_text_hits": len(candidate_metadata_hits),
    "metadata_map_selected": "NO",
    "specific_chromhmm_file_selected": "NO",
    "build_confirmed_for_file": "NO",
    "candidate_overlap_performed": "NO",
    "candidate_chromatin_state_evidence_retrieved": "NO",
    "biological_interpretation_performed": "NO",
    "functional_validation_claim_authorised": "NO",
    "broad_biological_interpretation_authorised": "NO",
    "outputs": {
        "routes": str(routes_path),
        "links_inventory": str(links_path),
        "candidate_links": str(immune_links_path),
        "immune_term_hits": str(hits_path),
        "summary": str(summary_path),
    },
    "boundary": (
        "This checkpoint discovers Roadmap metadata/file-map routes only. It does not select a final file, "
        "download full chromatin-state data, perform candidate overlap, retrieve candidate evidence, validate "
        "function or authorise biological interpretation."
    ),
}

OUTJSON.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
print(json.dumps(payload, indent=2))
