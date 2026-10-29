#!/usr/bin/env python3
"""Query GWAS Catalog for Phase D candidate rsIDs.

This checkpoint performs controlled direct rsID lookups against the GWAS
Catalog REST API. It records direct rsID availability and lightweight record
counts only. It does not interpret associations, validate function, infer
causality, infer selection, assign ancestry effects or authorise broad
biological interpretation.
"""

from __future__ import annotations

import csv
import json
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

CANDIDATE_INVENTORY = Path("11_results/phaseD/regulatory_functional_evidence/candidate_input_inventory/phaseD_candidate_input_inventory_v1.tsv")
OUTDIR = Path("11_results/phaseD/regulatory_functional_evidence/query_pass1/gwas_catalog")
OUTJSON = Path("results/qc/phaseD_gwas_catalog_direct_rsid_query_v1.json")

OUTDIR.mkdir(parents=True, exist_ok=True)
OUTJSON.parent.mkdir(parents=True, exist_ok=True)

BASE_ENDPOINTS = [
    # v1 endpoint retained as a practical direct-SNP endpoint while v2 migration is ongoing.
    "https://www.ebi.ac.uk/gwas/rest/api/singleNucleotidePolymorphisms/{rsid}",
]

failures: list[str] = []

if not CANDIDATE_INVENTORY.exists():
    failures.append(f"Missing candidate inventory: {CANDIDATE_INVENTORY}")

candidates = []
if CANDIDATE_INVENTORY.exists():
    with CANDIDATE_INVENTORY.open("r", encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            candidates.append({
                "phaseD_candidate_index": row.get("phaseD_candidate_index", ""),
                "allele_key": row.get("allele_key", ""),
                "gene": row.get("gene", ""),
                "rsid": row.get("mapped_identity", ""),
                "candidate_specific_cautions": row.get("candidate_specific_cautions", ""),
            })

if not candidates:
    failures.append("No candidates loaded from inventory.")

def get_json(url: str, timeout: int = 30) -> tuple[int | None, dict | list | None, str]:
    req = urllib.request.Request(
        url,
        headers={"User-Agent": "Project6-PhaseD-GWASCatalog-rsID-query/1.0"},
    )
    try:
        with urllib.request.urlopen(req, timeout=timeout) as resp:
            text = resp.read().decode("utf-8", errors="replace")
            try:
                payload = json.loads(text)
            except json.JSONDecodeError:
                return resp.status, None, "JSONDecodeError"
            return resp.status, payload, ""
    except urllib.error.HTTPError as exc:
        return exc.code, None, f"HTTPError: {exc}"
    except Exception as exc:
        return None, None, f"{type(exc).__name__}: {exc}"

def count_hal_embedded(payload: dict | list | None) -> tuple[int, str]:
    if not isinstance(payload, dict):
        return 0, ""
    embedded = payload.get("_embedded", {})
    if not isinstance(embedded, dict):
        return 0, ""
    counts = []
    keys = []
    for key, value in embedded.items():
        keys.append(key)
        if isinstance(value, list):
            counts.append(len(value))
        elif isinstance(value, dict):
            counts.append(len(value))
    return sum(counts), ";".join(keys)

rows = []
raw_summary = {}

for cand in candidates:
    rsid = cand["rsid"]
    encoded = urllib.parse.quote(rsid, safe="")
    direct_record_found = "NO"
    endpoint_used = ""
    http_status = ""
    error_message = ""
    embedded_record_count = 0
    embedded_keys = ""
    payload_keys = ""

    for template in BASE_ENDPOINTS:
        url = template.format(rsid=encoded)
        status, payload, err = get_json(url)
        endpoint_used = url
        http_status = "" if status is None else str(status)
        error_message = err

        if payload is not None:
            direct_record_found = "YES"
            if isinstance(payload, dict):
                payload_keys = ";".join(sorted(payload.keys()))
            embedded_record_count, embedded_keys = count_hal_embedded(payload)
            raw_summary[rsid] = {
                "endpoint": url,
                "http_status": http_status,
                "payload_keys": payload_keys,
                "embedded_keys": embedded_keys,
                "embedded_record_count": embedded_record_count,
            }
            break

        # Conservative delay; API allows more, but this keeps the checkpoint gentle.
        time.sleep(0.25)

    rows.append({
        "phaseD_candidate_index": cand["phaseD_candidate_index"],
        "allele_key": cand["allele_key"],
        "gene": cand["gene"],
        "rsid": rsid,
        "resource": "GWAS_Catalog",
        "query_type": "direct_rsid_lookup",
        "endpoint_used": endpoint_used,
        "http_status": http_status,
        "direct_record_found": direct_record_found,
        "payload_keys": payload_keys,
        "embedded_keys": embedded_keys,
        "embedded_record_count": embedded_record_count,
        "candidate_query_executed": "YES",
        "candidate_level_evidence_retrieved": "DIRECT_RECORD_METADATA_ONLY" if direct_record_found == "YES" else "NO_RECORD_RETURNED",
        "association_interpretation_performed": "NO",
        "causality_claim_made": "NO",
        "functional_validation_claim_made": "NO",
        "selection_or_population_claim_made": "NO",
        "ancestry_effect_claim_made": "NO",
        "candidate_specific_cautions": cand["candidate_specific_cautions"],
        "error_message": error_message,
    })

access_path = OUTDIR / "phaseD_gwas_catalog_direct_rsid_query_results_v1.tsv"
with access_path.open("w", encoding="utf-8", newline="") as handle:
    fieldnames = list(rows[0].keys()) if rows else ["rsid"]
    writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)

found = sum(1 for r in rows if r["direct_record_found"] == "YES")
not_found = sum(1 for r in rows if r["direct_record_found"] != "YES")

summary_path = OUTDIR / "phaseD_gwas_catalog_direct_rsid_query_summary_v1.tsv"
summary_rows = [
    {"metric": "candidate_rsids_queried", "value": len(rows)},
    {"metric": "direct_records_found", "value": found},
    {"metric": "direct_records_not_found", "value": not_found},
    {"metric": "coordinate_query_used", "value": "NO"},
    {"metric": "ld_proxy_query_used", "value": "NO"},
    {"metric": "association_interpretation_performed", "value": "NO"},
    {"metric": "new_analysis_performed", "value": "NO"},
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

raw_summary_path = OUTDIR / "phaseD_gwas_catalog_direct_rsid_query_raw_payload_summary_v1.json"
raw_summary_path.write_text(json.dumps(raw_summary, indent=2) + "\n", encoding="utf-8")

payload = {
    "status": "PASS" if not failures else "FAIL",
    "operation": "Controlled GWAS Catalog direct rsID lookup for Phase D candidates.",
    "candidate_rsids_queried": len(rows),
    "direct_records_found": found,
    "direct_records_not_found": not_found,
    "coordinate_query_used": "NO",
    "ld_proxy_query_used": "NO",
    "association_interpretation_performed": "NO",
    "new_analysis_performed": "NO",
    "new_biological_claims_added": "NO",
    "causality_claim_authorised": "NO",
    "functional_validation_claim_authorised": "NO",
    "selection_or_population_claim_authorised": "NO",
    "broad_biological_interpretation_authorised": "NO",
    "outputs": {
        "results": str(access_path),
        "summary": str(summary_path),
        "raw_payload_summary": str(raw_summary_path),
    },
    "failures": failures,
    "boundary": (
        "This checkpoint performs direct rsID availability lookup in GWAS Catalog only. "
        "Returned records are contextual catalog metadata and are not interpreted as proof of function, "
        "causality, selection, ancestry effect or population-specific biology."
    ),
}

OUTJSON.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
print(json.dumps(payload, indent=2))

if payload["status"] != "PASS":
    raise SystemExit(1)
