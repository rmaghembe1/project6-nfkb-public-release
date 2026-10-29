#!/usr/bin/env python3
"""Diagnose ENCODE SCREEN/cCRE access before overlap execution.

This checkpoint tests whether ENCODE SCREEN/cCRE resources are reachable and
whether a future GRCh38 exact-position overlap can be implemented reproducibly.
It does not perform candidate cCRE overlap, retrieve candidate regulatory
evidence, validate function, infer causality, infer selection, assign ancestry
effects or authorise broad biological interpretation.
"""

from __future__ import annotations

import csv
import json
import urllib.error
import urllib.request
from pathlib import Path

OUTDIR = Path("11_results/phaseD/regulatory_functional_evidence/query_pass1/encode_ccre/access_diagnostics")
OUTJSON = Path("results/qc/phaseD_encode_ccre_access_diagnostics_v1.json")

OUTDIR.mkdir(parents=True, exist_ok=True)
OUTJSON.parent.mkdir(parents=True, exist_ok=True)

TEST_URLS = [
    {
        "label": "encode_screen_home",
        "url": "https://screen.encodeproject.org/",
        "route_type": "portal",
    },
    {
        "label": "encode_screen_search",
        "url": "https://screen.encodeproject.org/search/",
        "route_type": "portal_search",
    },
    {
        "label": "encode_project_search_ccre",
        "url": "https://www.encodeproject.org/search/?type=Annotation&annotation_type=candidate+Cis-Regulatory+Elements",
        "route_type": "encode_portal_search",
    },
]

def check_url(url: str, timeout: int = 45) -> dict[str, str]:
    req = urllib.request.Request(
        url,
        headers={
            "User-Agent": "Project6-PhaseD-ENCODE-cCRE-access-diagnostics/1.0",
            "Accept": "text/html,application/json,text/plain,*/*",
        },
    )
    try:
        with urllib.request.urlopen(req, timeout=timeout) as resp:
            raw = resp.read(1200).decode("utf-8", errors="replace")
            return {
                "http_status": str(resp.status),
                "reachable_2xx": "YES" if 200 <= resp.status < 300 else "NO",
                "content_type": resp.headers.get("content-type", ""),
                "bytes_read": str(len(raw.encode("utf-8"))),
                "response_excerpt": raw.replace("\n", " ")[:400],
                "error_type": "",
                "error_message": "",
            }
    except urllib.error.HTTPError as exc:
        raw = ""
        try:
            raw = exc.read(1200).decode("utf-8", errors="replace")
        except Exception:
            pass
        return {
            "http_status": str(exc.code),
            "reachable_2xx": "NO",
            "content_type": exc.headers.get("content-type", "") if exc.headers else "",
            "bytes_read": str(len(raw.encode("utf-8"))),
            "response_excerpt": raw.replace("\n", " ")[:400],
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

rows = []
for item in TEST_URLS:
    status = check_url(item["url"])
    rows.append({
        "label": item["label"],
        "route_type": item["route_type"],
        "url": item["url"],
        **status,
        "candidate_overlap_performed": "NO",
        "candidate_regulatory_evidence_retrieved": "NO",
        "biological_interpretation_performed": "NO",
    })

diagnostics_path = OUTDIR / "phaseD_encode_ccre_access_diagnostics_v1.tsv"
with diagnostics_path.open("w", encoding="utf-8", newline="") as handle:
    fieldnames = list(rows[0].keys())
    writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)

reachable = sum(1 for r in rows if r["reachable_2xx"] == "YES")

summary_path = OUTDIR / "phaseD_encode_ccre_access_diagnostics_summary_v1.tsv"
summary_rows = [
    {"metric": "access_routes_checked", "value": len(rows)},
    {"metric": "routes_reachable_2xx", "value": reachable},
    {"metric": "candidate_overlap_performed", "value": "NO"},
    {"metric": "candidate_regulatory_evidence_retrieved", "value": "NO"},
    {"metric": "resource_version_confirmed", "value": "NO"},
    {"metric": "download_or_api_route_selected", "value": "NO"},
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
    "operation": "Diagnose ENCODE SCREEN/cCRE access before candidate overlap execution.",
    "access_routes_checked": len(rows),
    "routes_reachable_2xx": reachable,
    "candidate_overlap_performed": "NO",
    "candidate_regulatory_evidence_retrieved": "NO",
    "resource_version_confirmed": "NO",
    "download_or_api_route_selected": "NO",
    "biological_interpretation_performed": "NO",
    "functional_validation_claim_authorised": "NO",
    "selection_or_population_claim_authorised": "NO",
    "broad_biological_interpretation_authorised": "NO",
    "outputs": {
        "diagnostics": str(diagnostics_path),
        "summary": str(summary_path),
    },
    "boundary": (
        "This checkpoint diagnoses ENCODE cCRE access only. It does not perform candidate overlap, "
        "retrieve regulatory evidence, validate function, infer selection, assign ancestry effects or "
        "authorise broad biological interpretation."
    ),
}

OUTJSON.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
print(json.dumps(payload, indent=2))
