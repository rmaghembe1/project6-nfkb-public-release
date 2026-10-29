#!/usr/bin/env python3
"""Diagnose Roadmap ChromHMM access before candidate overlap.

This checkpoint tests Roadmap ChromHMM route/access feasibility only. It does
not download full chromatin-state datasets, perform candidate overlap, retrieve
candidate chromatin-state evidence, validate function, infer causality, infer
selection, assign ancestry effects or authorise broad biological interpretation.
"""

from __future__ import annotations

import csv
import json
import shutil
import urllib.error
import urllib.request
from pathlib import Path

OUTDIR = Path("11_results/phaseD/regulatory_functional_evidence/query_pass1/roadmap_chromhmm/access_diagnostics")
OUTJSON = Path("results/qc/phaseD_roadmap_chromhmm_access_diagnostics_v1.json")

OUTDIR.mkdir(parents=True, exist_ok=True)
OUTJSON.parent.mkdir(parents=True, exist_ok=True)

ROUTES = [
    {
        "route_id": "roadmap_main_portal",
        "resource_name": "Roadmap_Epigenomics_web_portal",
        "url": "https://egg2.wustl.edu/roadmap/web_portal/",
        "route_type": "portal",
        "candidate_overlap_ready": "NO",
        "reason": "portal route only; specific ChromHMM file route still required",
    },
    {
        "route_id": "roadmap_chromhmm_15state_browser",
        "resource_name": "Roadmap_ChromHMM_15state_browser",
        "url": "https://egg2.wustl.edu/roadmap/web_portal/chr_state_learning.html",
        "route_type": "documentation_page",
        "candidate_overlap_ready": "NO",
        "reason": "documentation route; specific file route and build status still required",
    },
    {
        "route_id": "roadmap_chromhmm_imputed12marks_dir",
        "resource_name": "Roadmap_ChromHMM_imputed12marks_directory",
        "url": "https://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/imputed12marks/",
        "route_type": "directory",
        "candidate_overlap_ready": "NO",
        "reason": "directory route; specific cell/tissue BED files must be selected before overlap",
    },
]

def check_url(url: str, timeout: int = 45) -> dict[str, str]:
    req = urllib.request.Request(
        url,
        headers={
            "User-Agent": "Project6-PhaseD-Roadmap-ChromHMM-access-diagnostics/1.0",
            "Accept": "text/html,text/plain,application/json,*/*",
        },
    )
    try:
        with urllib.request.urlopen(req, timeout=timeout) as resp:
            raw = resp.read(1500)
            return {
                "http_status": str(resp.status),
                "reachable_2xx": "YES" if 200 <= resp.status < 300 else "NO",
                "content_type": resp.headers.get("content-type", ""),
                "bytes_read": str(len(raw)),
                "response_excerpt": raw.decode("utf-8", errors="replace").replace("\n", " ")[:500],
                "error_type": "",
                "error_message": "",
            }
    except urllib.error.HTTPError as exc:
        raw = b""
        try:
            raw = exc.read(1500)
        except Exception:
            pass
        return {
            "http_status": str(exc.code),
            "reachable_2xx": "NO",
            "content_type": exc.headers.get("content-type", "") if exc.headers else "",
            "bytes_read": str(len(raw)),
            "response_excerpt": raw.decode("utf-8", errors="replace").replace("\n", " ")[:500],
            "error_type": "HTTPError",
            "error_message": str(exc),
        }
    except Exception as exc:
        return {
            "http_status": "",
            "reachable_2xx": "NO",
            "content_type": "",
            "bytes_read": "0",
            "response_excerpt": "",
            "error_type": type(exc).__name__,
            "error_message": str(exc)[:300],
        }

tool_rows = []
for tool in ["bedtools", "wget", "curl", "python"]:
    path = shutil.which(tool)
    tool_rows.append({
        "tool": tool,
        "available": "YES" if path else "NO",
        "path": path or "",
        "candidate_overlap_performed": "NO",
    })

route_rows = []
for route in ROUTES:
    status = check_url(route["url"])
    route_rows.append({
        **route,
        **status,
        "candidate_overlap_performed": "NO",
        "candidate_chromatin_state_evidence_retrieved": "NO",
        "biological_interpretation_performed": "NO",
    })

routes_path = OUTDIR / "phaseD_roadmap_chromhmm_access_routes_v1.tsv"
with routes_path.open("w", encoding="utf-8", newline="") as handle:
    fieldnames = list(route_rows[0].keys())
    writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(route_rows)

tools_path = OUTDIR / "phaseD_roadmap_chromhmm_tool_availability_v1.tsv"
with tools_path.open("w", encoding="utf-8", newline="") as handle:
    fieldnames = list(tool_rows[0].keys())
    writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(tool_rows)

reachable = sum(1 for r in route_rows if r["reachable_2xx"] == "YES")
bedtools_available = "YES" if any(r["tool"] == "bedtools" and r["available"] == "YES" for r in tool_rows) else "NO"

summary_path = OUTDIR / "phaseD_roadmap_chromhmm_access_diagnostics_summary_v1.tsv"
summary_rows = [
    {"metric": "routes_checked", "value": len(route_rows)},
    {"metric": "routes_reachable_2xx", "value": reachable},
    {"metric": "bedtools_available", "value": bedtools_available},
    {"metric": "specific_chromhmm_file_selected", "value": "NO"},
    {"metric": "cell_or_tissue_context_selected", "value": "NO"},
    {"metric": "candidate_overlap_performed", "value": "NO"},
    {"metric": "candidate_chromatin_state_evidence_retrieved", "value": "NO"},
    {"metric": "biological_interpretation_performed", "value": "NO"},
    {"metric": "functional_validation_claim_authorised", "value": "NO"},
    {"metric": "causality_claim_authorised", "value": "NO"},
    {"metric": "selection_or_population_claim_authorised", "value": "NO"},
    {"metric": "broad_biological_interpretation_authorised", "value": "NO"},
]
with summary_path.open("w", encoding="utf-8", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=["metric", "value"], delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(summary_rows)

payload = {
    "status": "PASS",
    "operation": "Diagnose Roadmap ChromHMM access before candidate overlap.",
    "routes_checked": len(route_rows),
    "routes_reachable_2xx": reachable,
    "bedtools_available": bedtools_available,
    "specific_chromhmm_file_selected": "NO",
    "cell_or_tissue_context_selected": "NO",
    "candidate_overlap_performed": "NO",
    "candidate_chromatin_state_evidence_retrieved": "NO",
    "biological_interpretation_performed": "NO",
    "functional_validation_claim_authorised": "NO",
    "causality_claim_authorised": "NO",
    "selection_or_population_claim_authorised": "NO",
    "broad_biological_interpretation_authorised": "NO",
    "outputs": {
        "routes": str(routes_path),
        "tools": str(tools_path),
        "summary": str(summary_path),
    },
    "boundary": (
        "This checkpoint diagnoses Roadmap ChromHMM access only. It does not select a specific "
        "cell/tissue file, perform candidate overlap, retrieve chromatin-state evidence, validate "
        "function, infer selection, assign ancestry effects or authorise broad biological interpretation."
    ),
}

OUTJSON.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
print(json.dumps(payload, indent=2))
