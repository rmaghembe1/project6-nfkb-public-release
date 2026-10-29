#!/usr/bin/env python3
"""Select an ENCODE cCRE resource route before candidate overlap.

This checkpoint evaluates candidate ENCODE/SCREEN/UCSC cCRE access routes for
future GRCh38 exact-position overlap. It does not download full cCRE datasets,
perform candidate overlap, retrieve candidate regulatory evidence, validate
function, infer causality, infer selection, assign ancestry effects or authorise
broad biological interpretation.
"""

from __future__ import annotations

import csv
import json
import shutil
import urllib.error
import urllib.request
from pathlib import Path

OUTDIR = Path("11_results/phaseD/regulatory_functional_evidence/query_pass1/encode_ccre/resource_route_selection")
OUTJSON = Path("results/qc/phaseD_encode_ccre_resource_route_selection_v1.json")

OUTDIR.mkdir(parents=True, exist_ok=True)
OUTJSON.parent.mkdir(parents=True, exist_ok=True)

ROUTES = [
    {
        "route_id": "screen_v4_downloads_page",
        "resource_name": "SCREEN_v4_downloads",
        "url": "https://screen-v4.wenglab.org/downloads",
        "route_type": "download_page",
        "candidate_overlap_ready": "NO",
        "reason": "page route only; exact file URL still needs confirmation",
    },
    {
        "route_id": "encode_portal_v4_ccre_analysis_json",
        "resource_name": "ENCODE_portal_cCRE_analysis",
        "url": "https://www.encodeproject.org/analyses/ENCAN463NBN/?format=json",
        "route_type": "portal_json",
        "candidate_overlap_ready": "NO",
        "reason": "metadata route; selected file accession/download URL must be extracted before overlap",
    },
    {
        "route_id": "ucsc_encode3_ccre_bigbed",
        "resource_name": "UCSC_ENCODE3_cCRE_bigBed",
        "url": "https://hgdownload.soe.ucsc.edu/gbdb/hg38/encode3/ccre/encodeCcreCombined.bb",
        "route_type": "bigBed",
        "candidate_overlap_ready": "YES_IF_BIGBEDTOBED_AVAILABLE",
        "reason": "direct hg38/GRCh38 BigBed route; requires bigBedToBed or equivalent interval reader",
    },
]

def check_url(url: str, timeout: int = 45) -> dict[str, str]:
    req = urllib.request.Request(
        url,
        headers={
            "User-Agent": "Project6-PhaseD-ENCODE-cCRE-route-selection/1.0",
            "Accept": "application/json,text/html,application/octet-stream,*/*",
            "Range": "bytes=0-1023",
        },
    )
    try:
        with urllib.request.urlopen(req, timeout=timeout) as resp:
            raw = resp.read(1024)
            return {
                "http_status": str(resp.status),
                "reachable": "YES" if 200 <= resp.status < 400 else "NO",
                "content_type": resp.headers.get("content-type", ""),
                "bytes_read": str(len(raw)),
                "error_type": "",
                "error_message": "",
            }
    except urllib.error.HTTPError as exc:
        return {
            "http_status": str(exc.code),
            "reachable": "NO",
            "content_type": exc.headers.get("content-type", "") if exc.headers else "",
            "bytes_read": "0",
            "error_type": "HTTPError",
            "error_message": str(exc),
        }
    except Exception as exc:
        return {
            "http_status": "",
            "reachable": "NO",
            "content_type": "",
            "bytes_read": "0",
            "error_type": type(exc).__name__,
            "error_message": str(exc)[:300],
        }

bigbedtobed_path = shutil.which("bigBedToBed")
bedtools_path = shutil.which("bedtools")

rows = []
for route in ROUTES:
    status = check_url(route["url"])
    rows.append({
        **route,
        **status,
        "bigBedToBed_available": "YES" if bigbedtobed_path else "NO",
        "bedtools_available": "YES" if bedtools_path else "NO",
        "candidate_overlap_performed": "NO",
        "candidate_regulatory_evidence_retrieved": "NO",
        "biological_interpretation_performed": "NO",
    })

# Selection rule: choose the first reachable route that is closest to reproducible overlap.
selected_route = "NONE"
selection_status = "NO_ROUTE_SELECTED"
selection_reason = "No route selected for overlap yet."

for row in rows:
    if row["route_id"] == "ucsc_encode3_ccre_bigbed" and row["reachable"] == "YES":
        selected_route = row["route_id"]
        if bigbedtobed_path:
            selection_status = "SELECTED_READY_FOR_FUTURE_OVERLAP_PROTOCOL"
            selection_reason = "UCSC ENCODE3 cCRE BigBed is reachable and bigBedToBed is available."
        else:
            selection_status = "SELECTED_BUT_TOOL_MISSING"
            selection_reason = "UCSC ENCODE3 cCRE BigBed is reachable, but bigBedToBed is not available."
        break

if selected_route == "NONE":
    for row in rows:
        if row["route_id"] == "encode_portal_v4_ccre_analysis_json" and row["reachable"] == "YES":
            selected_route = row["route_id"]
            selection_status = "METADATA_ROUTE_SELECTED_NEEDS_FILE_EXTRACTION"
            selection_reason = "ENCODE portal cCRE analysis JSON is reachable; file extraction route needed before overlap."
            break

route_table = OUTDIR / "phaseD_encode_ccre_resource_route_selection_v1.tsv"
with route_table.open("w", encoding="utf-8", newline="") as handle:
    fieldnames = list(rows[0].keys())
    writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)

summary_path = OUTDIR / "phaseD_encode_ccre_resource_route_selection_summary_v1.tsv"
summary_rows = [
    {"metric": "routes_evaluated", "value": len(rows)},
    {"metric": "routes_reachable", "value": sum(1 for r in rows if r["reachable"] == "YES")},
    {"metric": "selected_route", "value": selected_route},
    {"metric": "selection_status", "value": selection_status},
    {"metric": "selection_reason", "value": selection_reason},
    {"metric": "bigBedToBed_available", "value": "YES" if bigbedtobed_path else "NO"},
    {"metric": "bedtools_available", "value": "YES" if bedtools_path else "NO"},
    {"metric": "candidate_overlap_performed", "value": "NO"},
    {"metric": "candidate_regulatory_evidence_retrieved", "value": "NO"},
    {"metric": "biological_interpretation_performed", "value": "NO"},
    {"metric": "functional_validation_claim_authorised", "value": "NO"},
    {"metric": "selection_or_population_claim_authorised", "value": "NO"},
    {"metric": "broad_biological_interpretation_authorised", "value": "NO"},
]
with summary_path.open("w", encoding="utf-8", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=["metric", "value"], delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(summary_rows)

payload = {
    "status": "PASS",
    "operation": "Select ENCODE cCRE resource route before candidate overlap.",
    "routes_evaluated": len(rows),
    "routes_reachable": sum(1 for r in rows if r["reachable"] == "YES"),
    "selected_route": selected_route,
    "selection_status": selection_status,
    "selection_reason": selection_reason,
    "bigBedToBed_available": "YES" if bigbedtobed_path else "NO",
    "bedtools_available": "YES" if bedtools_path else "NO",
    "candidate_overlap_performed": "NO",
    "candidate_regulatory_evidence_retrieved": "NO",
    "biological_interpretation_performed": "NO",
    "functional_validation_claim_authorised": "NO",
    "selection_or_population_claim_authorised": "NO",
    "broad_biological_interpretation_authorised": "NO",
    "outputs": {
        "route_selection": str(route_table),
        "summary": str(summary_path),
    },
    "boundary": (
        "This checkpoint selects or defers an ENCODE cCRE resource route only. "
        "It does not perform candidate overlap, retrieve regulatory evidence, validate function, "
        "infer selection, assign ancestry effects or authorise broad biological interpretation."
    ),
}

OUTJSON.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
print(json.dumps(payload, indent=2))
