#!/usr/bin/env python3
"""Audit Phase D resource-level access.

This checkpoint tests resource-level URL accessibility only. It does not query
Phase D candidate variants, retrieve candidate-level regulatory evidence,
perform biological interpretation, validate function, infer selection, assign
ancestry effects or authorise broad biological interpretation.
"""

from __future__ import annotations

import csv
import json
import socket
import ssl
import time
import urllib.error
import urllib.request
from pathlib import Path

OUTDIR = Path("11_results/phaseD/regulatory_functional_evidence/resource_access_feasibility")
OUTJSON = Path("results/qc/phaseD_resource_level_access_audit_v1.json")

OUTDIR.mkdir(parents=True, exist_ok=True)
OUTJSON.parent.mkdir(parents=True, exist_ok=True)

resources = [
    {
        "resource": "GTEx",
        "url": "https://gtexportal.org/home/apiPage",
        "access_route": "public_api_documentation_page",
        "candidate_query_performed": "NO",
    },
    {
        "resource": "eQTL_Catalogue",
        "url": "https://www.ebi.ac.uk/eqtl/api-docs/",
        "access_route": "public_api_documentation_page",
        "candidate_query_performed": "NO",
    },
    {
        "resource": "DICE",
        "url": "https://dice-database.org/",
        "access_route": "public_portal_page",
        "candidate_query_performed": "NO",
    },
    {
        "resource": "BLUEPRINT",
        "url": "https://blueprint-data.bsc.es/",
        "access_route": "public_portal_page",
        "candidate_query_performed": "NO",
    },
    {
        "resource": "ENCODE_SCREEN_cCRE",
        "url": "https://screen.encodeproject.org/",
        "access_route": "public_portal_page",
        "candidate_query_performed": "NO",
    },
    {
        "resource": "Roadmap_ChromHMM",
        "url": "https://egg2.wustl.edu/roadmap/web_portal/",
        "access_route": "public_portal_page",
        "candidate_query_performed": "NO",
    },
    {
        "resource": "ABC_maps",
        "url": "https://www.engreitzlab.org/resources",
        "access_route": "public_resource_page",
        "candidate_query_performed": "NO",
    },
    {
        "resource": "GWAS_Catalog",
        "url": "https://www.ebi.ac.uk/gwas/rest/api",
        "access_route": "public_rest_api_base",
        "candidate_query_performed": "NO",
    },
    {
        "resource": "LDlink",
        "url": "https://ldlink.nih.gov/apiaccess",
        "access_route": "public_api_access_page_token_required_for_api_calls",
        "candidate_query_performed": "NO",
    },
]

def fetch_status(url: str, timeout: int = 20) -> dict[str, str]:
    req = urllib.request.Request(
        url,
        method="GET",
        headers={
            "User-Agent": "Project6-PhaseD-resource-access-audit/1.0",
        },
    )
    start = time.time()
    try:
        with urllib.request.urlopen(req, timeout=timeout, context=ssl.create_default_context()) as resp:
            elapsed_ms = int((time.time() - start) * 1000)
            content_type = resp.headers.get("content-type", "")
            # Read only a small header chunk to avoid downloading full resources.
            chunk = resp.read(512)
            return {
                "http_status": str(resp.status),
                "reachable": "YES" if 200 <= int(resp.status) < 400 else "NO",
                "content_type": content_type,
                "elapsed_ms": str(elapsed_ms),
                "bytes_read": str(len(chunk)),
                "error_type": "",
                "error_message": "",
            }
    except urllib.error.HTTPError as exc:
        elapsed_ms = int((time.time() - start) * 1000)
        return {
            "http_status": str(exc.code),
            "reachable": "YES" if 200 <= int(exc.code) < 400 else "NO",
            "content_type": exc.headers.get("content-type", "") if exc.headers else "",
            "elapsed_ms": str(elapsed_ms),
            "bytes_read": "0",
            "error_type": "HTTPError",
            "error_message": str(exc),
        }
    except (urllib.error.URLError, TimeoutError, socket.timeout) as exc:
        elapsed_ms = int((time.time() - start) * 1000)
        return {
            "http_status": "",
            "reachable": "NO",
            "content_type": "",
            "elapsed_ms": str(elapsed_ms),
            "bytes_read": "0",
            "error_type": type(exc).__name__,
            "error_message": str(exc)[:300],
        }

rows = []
for item in resources:
    status = fetch_status(item["url"])
    row = {
        "resource": item["resource"],
        "url": item["url"],
        "access_route": item["access_route"],
        **status,
        "candidate_query_performed": item["candidate_query_performed"],
        "candidate_level_evidence_retrieved": "NO",
        "biological_interpretation_performed": "NO",
    }
    rows.append(row)

access_path = OUTDIR / "phaseD_resource_level_access_audit_v1.tsv"
fieldnames = [
    "resource",
    "url",
    "access_route",
    "http_status",
    "reachable",
    "content_type",
    "elapsed_ms",
    "bytes_read",
    "error_type",
    "error_message",
    "candidate_query_performed",
    "candidate_level_evidence_retrieved",
    "biological_interpretation_performed",
]
with access_path.open("w", encoding="utf-8", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)

reachable_count = sum(1 for r in rows if r["reachable"] == "YES")
not_reachable_count = sum(1 for r in rows if r["reachable"] != "YES")

summary_path = OUTDIR / "phaseD_resource_level_access_audit_summary_v1.tsv"
summary_rows = [
    {"metric": "resources_tested", "value": len(rows)},
    {"metric": "resources_reachable", "value": reachable_count},
    {"metric": "resources_not_reachable", "value": not_reachable_count},
    {"metric": "candidate_query_performed", "value": "NO"},
    {"metric": "candidate_level_evidence_retrieved", "value": "NO"},
    {"metric": "new_analysis_performed", "value": "NO"},
    {"metric": "new_biological_claims_added", "value": "NO"},
    {"metric": "selection_or_population_claim_authorised", "value": "NO"},
    {"metric": "broad_biological_interpretation_authorised", "value": "NO"},
]
with summary_path.open("w", encoding="utf-8", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=["metric", "value"], delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(summary_rows)

payload = {
    "status": "PASS",
    "operation": "Audit Phase D resource-level access without candidate-level queries.",
    "resources_tested": len(rows),
    "resources_reachable": reachable_count,
    "resources_not_reachable": not_reachable_count,
    "candidate_query_performed": "NO",
    "candidate_level_evidence_retrieved": "NO",
    "new_analysis_performed": "NO",
    "new_biological_claims_added": "NO",
    "selection_or_population_claim_authorised": "NO",
    "broad_biological_interpretation_authorised": "NO",
    "outputs": {
        "access_audit": str(access_path),
        "summary": str(summary_path),
    },
    "boundary": (
        "This checkpoint tests resource-level access only. It does not query candidate variants, retrieve "
        "candidate-level regulatory evidence, validate function, infer selection, assign ancestry effects or "
        "authorise broad biological interpretation."
    ),
}

OUTJSON.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
print(json.dumps(payload, indent=2))
