#!/usr/bin/env python3
"""Diagnose GTEx direct rsID endpoint behavior.

This checkpoint tests GTEx endpoint formatting and response behavior before
candidate-level GTEx results are accepted. Candidate rsIDs are used only to
classify HTTP response behavior and variant-record availability. This checkpoint
does not interpret eQTLs, tissue specificity, causality, function, selection,
ancestry effects or population-specific biology.
"""

from __future__ import annotations

import csv
import json
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

OUTDIR = Path("11_results/phaseD/regulatory_functional_evidence/query_pass1/gtex/endpoint_diagnostics")
OUTJSON = Path("results/qc/phaseD_gtex_direct_rsid_endpoint_diagnostics_v1.json")

OUTDIR.mkdir(parents=True, exist_ok=True)
OUTJSON.parent.mkdir(parents=True, exist_ok=True)

candidate_rsids = ["rs1131877", "rs924735278", "rs376205580", "rs141807543"]

test_urls = [
    {
        "label": "redoc",
        "url": "https://gtexportal.org/api/v2/redoc",
        "candidate_specific": "NO",
        "query_class": "documentation",
    },
    {
        "label": "openapi_json",
        "url": "https://gtexportal.org/api/v2/openapi.json",
        "candidate_specific": "NO",
        "query_class": "documentation",
    },
]

for dataset in ["gtex_v8", "gtex_v10"]:
    for rsid in candidate_rsids:
        encoded = urllib.parse.quote(rsid, safe="")
        test_urls.extend([
            {
                "label": f"{rsid}_{dataset}_variant_snpId",
                "url": f"https://gtexportal.org/api/v2/dataset/variant?snpId={encoded}&datasetId={dataset}",
                "candidate_specific": "YES",
                "query_class": "variant_snpId_lookup",
            },
            {
                "label": f"{rsid}_{dataset}_singleTissueEqtl_snpId",
                "url": f"https://gtexportal.org/api/v2/association/singleTissueEqtl?snpId={encoded}&datasetId={dataset}&itemsPerPage=250",
                "candidate_specific": "YES",
                "query_class": "eqtl_snpId_lookup",
            },
        ])

def get_url(url: str, timeout: int = 45) -> dict[str, str]:
    req = urllib.request.Request(
        url,
        headers={
            "User-Agent": "Project6-PhaseD-GTEx-endpoint-diagnostics/1.0",
            "Accept": "application/json,text/html,text/plain,*/*",
        },
    )
    try:
        with urllib.request.urlopen(req, timeout=timeout) as resp:
            text = resp.read(2000).decode("utf-8", errors="replace")
            return {
                "http_status": str(resp.status),
                "reachable_2xx": "YES" if 200 <= resp.status < 300 else "NO",
                "content_type": resp.headers.get("content-type", ""),
                "bytes_read": str(len(text.encode("utf-8"))),
                "response_excerpt": text.replace("\n", " ")[:500],
                "error_type": "",
                "error_message": "",
            }
    except urllib.error.HTTPError as exc:
        text = ""
        try:
            text = exc.read(2000).decode("utf-8", errors="replace")
        except Exception:
            text = ""
        return {
            "http_status": str(exc.code),
            "reachable_2xx": "NO",
            "content_type": exc.headers.get("content-type", "") if exc.headers else "",
            "bytes_read": str(len(text.encode("utf-8"))),
            "response_excerpt": text.replace("\n", " ")[:500],
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
for item in test_urls:
    status = get_url(item["url"])
    rows.append({
        "label": item["label"],
        "query_class": item["query_class"],
        "url": item["url"],
        "candidate_specific": item["candidate_specific"],
        **status,
        "candidate_result_accepted": "NO",
        "eqtl_interpretation_performed": "NO",
        "tissue_specific_interpretation_performed": "NO",
    })

diagnostics_path = OUTDIR / "phaseD_gtex_direct_rsid_endpoint_diagnostics_v1.tsv"
with diagnostics_path.open("w", encoding="utf-8", newline="") as handle:
    fieldnames = list(rows[0].keys())
    writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)

candidate_rows = [r for r in rows if r["candidate_specific"] == "YES"]
candidate_2xx = sum(1 for r in candidate_rows if r["reachable_2xx"] == "YES")
candidate_non2xx = len(candidate_rows) - candidate_2xx

summary_path = OUTDIR / "phaseD_gtex_direct_rsid_endpoint_diagnostics_summary_v1.tsv"
summary_rows = [
    {"metric": "endpoint_tests", "value": len(rows)},
    {"metric": "candidate_endpoint_tests", "value": len(candidate_rows)},
    {"metric": "candidate_endpoint_2xx", "value": candidate_2xx},
    {"metric": "candidate_endpoint_non2xx", "value": candidate_non2xx},
    {"metric": "candidate_results_accepted", "value": "NO"},
    {"metric": "eqtl_interpretation_performed", "value": "NO"},
    {"metric": "tissue_specific_interpretation_performed", "value": "NO"},
    {"metric": "new_biological_claims_added", "value": "NO"},
    {"metric": "causality_claim_authorised", "value": "NO"},
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
    "operation": "Diagnose GTEx direct rsID endpoint behavior before accepting candidate-level GTEx outputs.",
    "endpoint_tests": len(rows),
    "candidate_endpoint_tests": len(candidate_rows),
    "candidate_endpoint_2xx": candidate_2xx,
    "candidate_endpoint_non2xx": candidate_non2xx,
    "candidate_results_accepted": "NO",
    "eqtl_interpretation_performed": "NO",
    "tissue_specific_interpretation_performed": "NO",
    "new_biological_claims_added": "NO",
    "causality_claim_authorised": "NO",
    "functional_validation_claim_authorised": "NO",
    "selection_or_population_claim_authorised": "NO",
    "broad_biological_interpretation_authorised": "NO",
    "outputs": {
        "diagnostics": str(diagnostics_path),
        "summary": str(summary_path),
    },
    "boundary": (
        "This checkpoint diagnoses GTEx endpoint behavior only. Candidate-specific tests classify HTTP response "
        "behavior and are not accepted as biological evidence, eQTL interpretation or tissue-specific interpretation."
    ),
}

OUTJSON.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
print(json.dumps(payload, indent=2))
