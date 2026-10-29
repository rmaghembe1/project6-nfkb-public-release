#!/usr/bin/env python3
"""Diagnose eQTL Catalogue association endpoint behavior.

This checkpoint tests endpoint formatting and response behavior for eQTL
Catalogue before candidate-level evidence retrieval is accepted. It includes
candidate rsIDs only to test HTTP response classes, not to interpret QTL
evidence. It does not retrieve or interpret candidate-level QTL biology,
validate function, infer causality, selection, ancestry effects or authorise
broad biological interpretation.
"""

from __future__ import annotations

import csv
import json
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

OUTDIR = Path("11_results/phaseD/regulatory_functional_evidence/query_pass1/eqtl_catalogue/endpoint_diagnostics")
OUTJSON = Path("results/qc/phaseD_eqtl_catalogue_endpoint_diagnostics_v1.json")

OUTDIR.mkdir(parents=True, exist_ok=True)
OUTJSON.parent.mkdir(parents=True, exist_ok=True)

test_queries = [
    {"label": "docs_page", "url": "https://www.ebi.ac.uk/eqtl/api-docs/", "candidate_specific": "NO"},
    {"label": "studies_endpoint", "url": "https://www.ebi.ac.uk/eqtl/api/studies", "candidate_specific": "NO"},
    {"label": "associations_base", "url": "https://www.ebi.ac.uk/eqtl/api/associations", "candidate_specific": "NO"},
    {"label": "example_variant_path", "url": "https://www.ebi.ac.uk/eqtl/api/associations/chr19_226776_C_T", "candidate_specific": "NO"},
    {"label": "example_variant_param", "url": "https://www.ebi.ac.uk/eqtl/api/associations/?variant_id=chr19_226776_C_T", "candidate_specific": "NO"},
]

candidate_rsids = ["rs1131877", "rs924735278", "rs376205580", "rs141807543"]
for rsid in candidate_rsids:
    encoded = urllib.parse.quote(rsid, safe="")
    test_queries.extend([
        {"label": f"{rsid}_path", "url": f"https://www.ebi.ac.uk/eqtl/api/associations/{encoded}", "candidate_specific": "YES"},
        {"label": f"{rsid}_param_variant_id", "url": f"https://www.ebi.ac.uk/eqtl/api/associations/?variant_id={encoded}", "candidate_specific": "YES"},
        {"label": f"{rsid}_param_rsid", "url": f"https://www.ebi.ac.uk/eqtl/api/associations/?rsid={encoded}", "candidate_specific": "YES"},
    ])

def request_url(url: str, timeout: int = 45) -> dict[str, str]:
    req = urllib.request.Request(
        url,
        headers={
            "User-Agent": "Project6-PhaseD-eQTLCatalogue-endpoint-diagnostics/1.0",
            "Accept": "application/json,text/plain,*/*",
        },
    )
    try:
        with urllib.request.urlopen(req, timeout=timeout) as resp:
            raw = resp.read(1000).decode("utf-8", errors="replace")
            return {
                "http_status": str(resp.status),
                "reachable_2xx": "YES" if 200 <= resp.status < 300 else "NO",
                "content_type": resp.headers.get("content-type", ""),
                "bytes_read": str(len(raw.encode("utf-8"))),
                "response_excerpt": raw.replace("\n", " ")[:300],
                "error_type": "",
                "error_message": "",
            }
    except urllib.error.HTTPError as exc:
        raw = ""
        try:
            raw = exc.read(1000).decode("utf-8", errors="replace")
        except Exception:
            raw = ""
        return {
            "http_status": str(exc.code),
            "reachable_2xx": "NO",
            "content_type": exc.headers.get("content-type", "") if exc.headers else "",
            "bytes_read": str(len(raw.encode("utf-8"))),
            "response_excerpt": raw.replace("\n", " ")[:300],
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
for item in test_queries:
    status = request_url(item["url"])
    rows.append({
        "label": item["label"],
        "url": item["url"],
        "candidate_specific": item["candidate_specific"],
        **status,
        "candidate_level_interpretation_performed": "NO",
        "qtl_interpretation_performed": "NO",
    })

diagnostics_path = OUTDIR / "phaseD_eqtl_catalogue_endpoint_diagnostics_v1.tsv"
with diagnostics_path.open("w", encoding="utf-8", newline="") as handle:
    fieldnames = list(rows[0].keys())
    writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)

candidate_rows = [r for r in rows if r["candidate_specific"] == "YES"]
candidate_2xx = sum(1 for r in candidate_rows if r["reachable_2xx"] == "YES")
candidate_400 = sum(1 for r in candidate_rows if r["http_status"] == "400")

summary_path = OUTDIR / "phaseD_eqtl_catalogue_endpoint_diagnostics_summary_v1.tsv"
summary_rows = [
    {"metric": "endpoint_tests", "value": len(rows)},
    {"metric": "candidate_endpoint_tests", "value": len(candidate_rows)},
    {"metric": "candidate_endpoint_2xx", "value": candidate_2xx},
    {"metric": "candidate_endpoint_400", "value": candidate_400},
    {"metric": "candidate_query_results_accepted", "value": "NO"},
    {"metric": "endpoint_diagnostics_required_before_query_commit", "value": "YES"},
    {"metric": "qtl_interpretation_performed", "value": "NO"},
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
    "operation": "Diagnose eQTL Catalogue endpoint behavior before accepting candidate-level query outputs.",
    "endpoint_tests": len(rows),
    "candidate_endpoint_tests": len(candidate_rows),
    "candidate_endpoint_2xx": candidate_2xx,
    "candidate_endpoint_400": candidate_400,
    "candidate_query_results_accepted": "NO",
    "qtl_interpretation_performed": "NO",
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
        "This checkpoint diagnoses eQTL Catalogue endpoint behavior only. Candidate-specific tests are used "
        "to classify HTTP response behavior and are not accepted as biological evidence or QTL interpretation."
    ),
}

OUTJSON.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
print(json.dumps(payload, indent=2))
