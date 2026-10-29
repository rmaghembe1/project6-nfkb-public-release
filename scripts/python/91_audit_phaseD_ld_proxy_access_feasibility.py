#!/usr/bin/env python3
"""Audit LD-proxy metadata access feasibility.

This checkpoint checks whether LD proxy metadata retrieval can proceed under
the Phase D LD-proxy protocol. It does not retrieve proxies, does not query
candidate LD, does not interpret population biology, and does not authorise
selection, ancestry-effect or functional claims.
"""

from __future__ import annotations

import csv
import json
import os
import urllib.error
import urllib.request
from pathlib import Path

OUTDIR = Path("11_results/phaseD/regulatory_functional_evidence/query_pass1/ld_proxy/access_feasibility")
OUTJSON = Path("results/qc/phaseD_ld_proxy_access_feasibility_v1.json")

OUTDIR.mkdir(parents=True, exist_ok=True)
OUTJSON.parent.mkdir(parents=True, exist_ok=True)

token_names = ["LDLINK_TOKEN", "LDLINK_API_TOKEN"]
token_present = any(os.environ.get(name, "").strip() for name in token_names)

def check_url(url: str, timeout: int = 30) -> dict[str, str]:
    req = urllib.request.Request(
        url,
        headers={"User-Agent": "Project6-PhaseD-LDlink-access-feasibility/1.0"},
    )
    try:
        with urllib.request.urlopen(req, timeout=timeout) as resp:
            raw = resp.read(512).decode("utf-8", errors="replace")
            return {
                "url": url,
                "http_status": str(resp.status),
                "reachable_2xx": "YES" if 200 <= resp.status < 300 else "NO",
                "content_type": resp.headers.get("content-type", ""),
                "bytes_read": str(len(raw.encode("utf-8"))),
                "error_type": "",
                "error_message": "",
            }
    except urllib.error.HTTPError as exc:
        return {
            "url": url,
            "http_status": str(exc.code),
            "reachable_2xx": "NO",
            "content_type": exc.headers.get("content-type", "") if exc.headers else "",
            "bytes_read": "0",
            "error_type": "HTTPError",
            "error_message": str(exc),
        }
    except Exception as exc:
        return {
            "url": url,
            "http_status": "",
            "reachable_2xx": "NO",
            "content_type": "",
            "bytes_read": "0",
            "error_type": type(exc).__name__,
            "error_message": str(exc)[:300],
        }

rows = []
for label, url in [
    ("ldlink_api_access_page", "https://ldlink.nih.gov/apiaccess"),
    ("ldlink_home", "https://ldlink.nih.gov/"),
]:
    status = check_url(url)
    rows.append({
        "resource": "LDlink",
        "check_label": label,
        **status,
        "ld_proxy_retrieval_performed": "NO",
        "candidate_ld_query_performed": "NO",
        "population_interpretation_performed": "NO",
    })

access_path = OUTDIR / "phaseD_ld_proxy_access_feasibility_v1.tsv"
with access_path.open("w", encoding="utf-8", newline="") as handle:
    fieldnames = list(rows[0].keys())
    writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)

summary_path = OUTDIR / "phaseD_ld_proxy_access_feasibility_summary_v1.tsv"
summary_rows = [
    {"metric": "ldlink_pages_checked", "value": len(rows)},
    {"metric": "ldlink_pages_reachable_2xx", "value": sum(1 for r in rows if r["reachable_2xx"] == "YES")},
    {"metric": "ldlink_token_present_in_environment", "value": "YES" if token_present else "NO"},
    {"metric": "token_environment_names_checked", "value": ";".join(token_names)},
    {"metric": "ld_proxy_retrieval_performed", "value": "NO"},
    {"metric": "candidate_ld_query_performed", "value": "NO"},
    {"metric": "population_panel_selected", "value": "NO"},
    {"metric": "population_interpretation_performed", "value": "NO"},
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
    "operation": "Audit LD-proxy metadata access feasibility without candidate LD retrieval.",
    "ldlink_pages_checked": len(rows),
    "ldlink_pages_reachable_2xx": sum(1 for r in rows if r["reachable_2xx"] == "YES"),
    "ldlink_token_present_in_environment": "YES" if token_present else "NO",
    "ld_proxy_retrieval_performed": "NO",
    "candidate_ld_query_performed": "NO",
    "population_panel_selected": "NO",
    "population_interpretation_performed": "NO",
    "new_biological_claims_added": "NO",
    "selection_or_population_claim_authorised": "NO",
    "broad_biological_interpretation_authorised": "NO",
    "outputs": {
        "access_feasibility": str(access_path),
        "summary": str(summary_path),
    },
    "boundary": (
        "This checkpoint audits LD-proxy access feasibility only. It does not retrieve LD proxies, "
        "query candidate LD, select a population panel, interpret population biology, infer selection "
        "or authorise broad biological interpretation."
    ),
}

OUTJSON.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
print(json.dumps(payload, indent=2))
