#!/usr/bin/env python3
"""Inventory Roadmap ChromHMM files before candidate overlap.

This checkpoint inventories candidate Roadmap ChromHMM files from the selected
imputed12marks directory and flags immune/blood-relevant file-name patterns.

It does not download full chromatin-state files, perform candidate overlap,
retrieve candidate chromatin-state evidence, validate function, infer causality,
infer selection, assign ancestry effects or authorise broad biological interpretation.
"""

from __future__ import annotations

import csv
import json
import re
import urllib.request
from html.parser import HTMLParser
from pathlib import Path

OUTDIR = Path("11_results/phaseD/regulatory_functional_evidence/query_pass1/roadmap_chromhmm/file_inventory")
OUTJSON = Path("results/qc/phaseD_roadmap_chromhmm_file_inventory_v1.json")

OUTDIR.mkdir(parents=True, exist_ok=True)
OUTJSON.parent.mkdir(parents=True, exist_ok=True)

BASE_URL = "https://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/imputed12marks/"

IMMUNE_BLOOD_PATTERNS = [
    "blood", "monocyte", "macrophage", "neutrophil", "tcell", "t-cell",
    "bcell", "b-cell", "nk", "cd4", "cd8", "pbmc", "lymph", "thymus",
    "spleen", "hematopoietic", "haematopoietic", "immune"
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

req = urllib.request.Request(
    BASE_URL,
    headers={"User-Agent": "Project6-PhaseD-Roadmap-ChromHMM-file-inventory/1.0"},
)
with urllib.request.urlopen(req, timeout=60) as resp:
    html = resp.read().decode("utf-8", errors="replace")

parser = LinkParser()
parser.feed(html)

rows = []
for href in parser.links:
    if href.startswith("?") or href in {"../", "/"}:
        continue
    file_name = href.rstrip("/")
    full_url = BASE_URL + href
    lower = file_name.lower()

    file_type = "directory" if href.endswith("/") else "file"
    apparent_bed_like = "YES" if re.search(r"\.(bed|bed\.gz|mnemonics\.bed|segments\.bed|txt|gz)$", lower) else "UNKNOWN"
    immune_flag = "YES" if any(pat in lower for pat in IMMUNE_BLOOD_PATTERNS) else "NO"

    rows.append({
        "file_name": file_name,
        "url": full_url,
        "file_type": file_type,
        "apparent_bed_like": apparent_bed_like,
        "immune_blood_relevant_name_flag": immune_flag,
        "apparent_build": "UNCONFIRMED",
        "cell_or_tissue_context": "UNRESOLVED_FROM_FILENAME",
        "candidate_overlap_performed": "NO",
        "candidate_chromatin_state_evidence_retrieved": "NO",
        "biological_interpretation_performed": "NO",
    })

inventory_path = OUTDIR / "phaseD_roadmap_chromhmm_file_inventory_v1.tsv"
with inventory_path.open("w", encoding="utf-8", newline="") as handle:
    fieldnames = list(rows[0].keys()) if rows else [
        "file_name", "url", "file_type", "apparent_bed_like",
        "immune_blood_relevant_name_flag", "apparent_build",
        "cell_or_tissue_context", "candidate_overlap_performed",
        "candidate_chromatin_state_evidence_retrieved",
        "biological_interpretation_performed"
    ]
    writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)

immune_rows = [r for r in rows if r["immune_blood_relevant_name_flag"] == "YES"]
immune_path = OUTDIR / "phaseD_roadmap_chromhmm_immune_blood_candidate_files_v1.tsv"
with immune_path.open("w", encoding="utf-8", newline="") as handle:
    fieldnames = list(rows[0].keys()) if rows else list(csv.DictReader(open(inventory_path, encoding="utf-8"), delimiter="\t").fieldnames or [])
    writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(immune_rows)

summary_path = OUTDIR / "phaseD_roadmap_chromhmm_file_inventory_summary_v1.tsv"
summary_rows = [
    {"metric": "selected_route_url", "value": BASE_URL},
    {"metric": "directory_links_parsed", "value": len(parser.links)},
    {"metric": "inventory_rows", "value": len(rows)},
    {"metric": "immune_blood_candidate_rows_by_filename", "value": len(immune_rows)},
    {"metric": "specific_file_selected_for_overlap", "value": "NO"},
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
    "operation": "Inventory Roadmap ChromHMM files before candidate overlap.",
    "selected_route_url": BASE_URL,
    "directory_links_parsed": len(parser.links),
    "inventory_rows": len(rows),
    "immune_blood_candidate_rows_by_filename": len(immune_rows),
    "specific_file_selected_for_overlap": "NO",
    "build_confirmed_for_file": "NO",
    "candidate_overlap_performed": "NO",
    "candidate_chromatin_state_evidence_retrieved": "NO",
    "biological_interpretation_performed": "NO",
    "functional_validation_claim_authorised": "NO",
    "broad_biological_interpretation_authorised": "NO",
    "outputs": {
        "inventory": str(inventory_path),
        "immune_blood_candidates": str(immune_path),
        "summary": str(summary_path),
    },
    "boundary": (
        "This checkpoint inventories Roadmap ChromHMM files only. It does not select a final file, "
        "download full chromatin-state data, perform candidate overlap, retrieve candidate evidence, "
        "validate function or authorise biological interpretation."
    ),
}

OUTJSON.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
print(json.dumps(payload, indent=2))
