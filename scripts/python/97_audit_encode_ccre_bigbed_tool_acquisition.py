#!/usr/bin/env python3
"""Audit bigBedToBed acquisition routes before ENCODE cCRE overlap.

This checkpoint tests local availability and feasible acquisition routes for
bigBedToBed. It does not install tools, download full cCRE data, convert BigBed,
perform candidate overlap, retrieve regulatory evidence, validate function or
authorise biological interpretation.
"""

from __future__ import annotations

import csv
import json
import shutil
import urllib.error
import urllib.request
from pathlib import Path

OUTDIR = Path("11_results/phaseD/regulatory_functional_evidence/query_pass1/encode_ccre/tooling")
OUTJSON = Path("results/qc/phaseD_encode_ccre_bigbed_tool_acquisition_v1.json")
OUTDIR.mkdir(parents=True, exist_ok=True)
OUTJSON.parent.mkdir(parents=True, exist_ok=True)

TOOLS = ["bigBedToBed", "conda", "mamba", "micromamba", "wget", "curl", "bedtools"]

ROUTES = [
    {
        "route": "bioconda_ucsc_bigbedtobed",
        "description": "Install ucsc-bigbedtobed using conda/mamba/micromamba with Bioconda.",
        "url": "https://bioconda.github.io/recipes/ucsc-bigbedtobed/README.html",
    },
    {
        "route": "ucsc_linux_x86_64_binary_directory",
        "description": "Download UCSC standalone Linux x86_64 binary directly if allowed.",
        "url": "https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/",
    },
    {
        "route": "ucsc_bigbed_help",
        "description": "UCSC bigBed help page documenting bigBedToBed.",
        "url": "https://genome.ucsc.edu/goldenpath/help/bigBed.html",
    },
]

def check_url(url: str, timeout: int = 30) -> dict[str, str]:
    req = urllib.request.Request(
        url,
        headers={
            "User-Agent": "Project6-PhaseD-bigBedToBed-acquisition-audit/1.0",
            "Accept": "text/html,text/plain,*/*",
        },
    )
    try:
        with urllib.request.urlopen(req, timeout=timeout) as resp:
            raw = resp.read(1024)
            return {
                "http_status": str(resp.status),
                "reachable_2xx": "YES" if 200 <= resp.status < 300 else "NO",
                "content_type": resp.headers.get("content-type", ""),
                "bytes_read": str(len(raw)),
                "error_type": "",
                "error_message": "",
            }
    except urllib.error.HTTPError as exc:
        return {
            "http_status": str(exc.code),
            "reachable_2xx": "NO",
            "content_type": exc.headers.get("content-type", "") if exc.headers else "",
            "bytes_read": "0",
            "error_type": "HTTPError",
            "error_message": str(exc),
        }
    except Exception as exc:
        return {
            "http_status": "",
            "reachable_2xx": "NO",
            "content_type": "",
            "bytes_read": "0",
            "error_type": type(exc).__name__,
            "error_message": str(exc)[:300],
        }

tool_rows = []
for tool in TOOLS:
    path = shutil.which(tool)
    tool_rows.append({
        "tool": tool,
        "available": "YES" if path else "NO",
        "path": path or "",
        "tool_install_performed": "NO",
        "candidate_overlap_performed": "NO",
    })

route_rows = []
for route in ROUTES:
    status = check_url(route["url"])
    route_rows.append({
        **route,
        **status,
        "tool_install_performed": "NO",
        "candidate_overlap_performed": "NO",
        "candidate_regulatory_evidence_retrieved": "NO",
    })

tool_table = OUTDIR / "phaseD_encode_ccre_bigbed_tool_availability_v1.tsv"
with tool_table.open("w", encoding="utf-8", newline="") as handle:
    fieldnames = list(tool_rows[0].keys())
    writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(tool_rows)

route_table = OUTDIR / "phaseD_encode_ccre_bigbed_tool_acquisition_routes_v1.tsv"
with route_table.open("w", encoding="utf-8", newline="") as handle:
    fieldnames = list(route_rows[0].keys())
    writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(route_rows)

available_tools = {r["tool"]: r["available"] for r in tool_rows}
conda_like_available = any(available_tools.get(t) == "YES" for t in ["conda", "mamba", "micromamba"])
download_tool_available = any(available_tools.get(t) == "YES" for t in ["wget", "curl"])
route_reachable = sum(1 for r in route_rows if r["reachable_2xx"] == "YES")

recommended_next_action = "NONE"
if available_tools.get("bigBedToBed") == "YES":
    recommended_next_action = "VERIFY_BIGBEDTOBED_AND_DEFINE_OVERLAP_EXECUTION"
elif conda_like_available:
    recommended_next_action = "INSTALL_UCSC_BIGBEDTOBED_FROM_BIOCONDA_IN_CONTROLLED_CHECKPOINT"
elif download_tool_available:
    recommended_next_action = "DOWNLOAD_UCSC_STANDALONE_BIGBEDTOBED_IN_CONTROLLED_CHECKPOINT"
else:
    recommended_next_action = "DEFER_AND_DEFINE_ENCODE_PORTAL_FILE_EXTRACTION_ROUTE"

summary_path = OUTDIR / "phaseD_encode_ccre_bigbed_tool_acquisition_summary_v1.tsv"
summary_rows = [
    {"metric": "bigBedToBed_available", "value": available_tools.get("bigBedToBed", "NO")},
    {"metric": "conda_like_available", "value": "YES" if conda_like_available else "NO"},
    {"metric": "download_tool_available", "value": "YES" if download_tool_available else "NO"},
    {"metric": "acquisition_routes_checked", "value": len(route_rows)},
    {"metric": "acquisition_routes_reachable_2xx", "value": route_reachable},
    {"metric": "recommended_next_action", "value": recommended_next_action},
    {"metric": "tool_install_performed", "value": "NO"},
    {"metric": "candidate_overlap_performed", "value": "NO"},
    {"metric": "candidate_regulatory_evidence_retrieved", "value": "NO"},
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
    "operation": "Audit bigBedToBed acquisition routes before ENCODE cCRE overlap.",
    "bigBedToBed_available": available_tools.get("bigBedToBed", "NO"),
    "conda_like_available": "YES" if conda_like_available else "NO",
    "download_tool_available": "YES" if download_tool_available else "NO",
    "acquisition_routes_checked": len(route_rows),
    "acquisition_routes_reachable_2xx": route_reachable,
    "recommended_next_action": recommended_next_action,
    "tool_install_performed": "NO",
    "candidate_overlap_performed": "NO",
    "candidate_regulatory_evidence_retrieved": "NO",
    "biological_interpretation_performed": "NO",
    "functional_validation_claim_authorised": "NO",
    "broad_biological_interpretation_authorised": "NO",
    "outputs": {
        "tool_availability": str(tool_table),
        "acquisition_routes": str(route_table),
        "summary": str(summary_path),
    },
    "boundary": (
        "This checkpoint audits bigBedToBed acquisition feasibility only. It does not install tools, "
        "download full cCRE data, perform candidate overlap, retrieve regulatory evidence, validate "
        "function or authorise biological interpretation."
    ),
}

OUTJSON.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
print(json.dumps(payload, indent=2))
