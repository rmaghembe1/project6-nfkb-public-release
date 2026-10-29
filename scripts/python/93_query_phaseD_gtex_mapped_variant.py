#!/usr/bin/env python3
"""Query GTEx mapped variant availability for Phase D candidate rs1131877.

This checkpoint performs a controlled GTEx mapped-variant query for the one
candidate with a confirmed GTEx variantId mapping from endpoint diagnostics:
TRAF3 rs1131877 -> chr14_102875712_T_C_b38.

It records endpoint status and returned record counts only. It does not
interpret eQTLs, tissue specificity, causality, function, selection, ancestry
effects or population-specific biology.
"""

from __future__ import annotations

import csv
import json
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

OUTDIR = Path("11_results/phaseD/regulatory_functional_evidence/query_pass1/gtex/mapped_variant_query")
OUTJSON = Path("results/qc/phaseD_gtex_mapped_variant_query_v1.json")

OUTDIR.mkdir(parents=True, exist_ok=True)
OUTJSON.parent.mkdir(parents=True, exist_ok=True)

candidate = {
    "gene": "TRAF3",
    "rsid": "rs1131877",
    "gtex_variant_id": "chr14_102875712_T_C_b38",
}

datasets = ["gtex_v8", "gtex_v10"]

def get_json(url: str, timeout: int = 60) -> tuple[int | None, object | None, str, str]:
    req = urllib.request.Request(
        url,
        headers={
            "User-Agent": "Project6-PhaseD-GTEx-mapped-variant-query/1.0",
            "Accept": "application/json,text/plain,*/*",
        },
    )
    try:
        with urllib.request.urlopen(req, timeout=timeout) as resp:
            text = resp.read().decode("utf-8", errors="replace")
            content_type = resp.headers.get("content-type", "")
            try:
                payload = json.loads(text)
            except json.JSONDecodeError:
                payload = text
            return resp.status, payload, content_type, ""
    except urllib.error.HTTPError as exc:
        text = ""
        try:
            text = exc.read(1000).decode("utf-8", errors="replace")
        except Exception:
            pass
        return exc.code, None, exc.headers.get("content-type", "") if exc.headers else "", f"HTTPError: {exc}; {text[:200]}"
    except Exception as exc:
        return None, None, "", f"{type(exc).__name__}: {exc}"

def payload_counts(payload: object | None) -> tuple[int, int, str, str]:
    """Return page count, total count, structure label and payload keys."""
    if payload is None:
        return 0, 0, "none", ""
    if isinstance(payload, dict):
        keys = ";".join(sorted(payload.keys()))
        data = payload.get("data")
        page_count = len(data) if isinstance(data, list) else 0
        total = page_count
        paging = payload.get("paging_info")
        if isinstance(paging, dict):
            total = int(paging.get("totalNumberOfItems", total) or 0)
        return page_count, total, "dict_data" if isinstance(data, list) else "dict_metadata", keys
    if isinstance(payload, list):
        return len(payload), len(payload), "list", "list"
    if isinstance(payload, str):
        return (1 if payload.strip() else 0), (1 if payload.strip() else 0), "text", "text"
    return 0, 0, type(payload).__name__, ""

rows = []
raw_summary = {}

for dataset in datasets:
    variant_encoded = urllib.parse.quote(candidate["gtex_variant_id"], safe="")
    endpoint = (
        "https://gtexportal.org/api/v2/association/singleTissueEqtl"
        f"?variantId={variant_encoded}&datasetId={dataset}&itemsPerPage=250"
    )

    status, payload, content_type, error_message = get_json(endpoint)
    page_count, total_count, payload_structure, payload_keys = payload_counts(payload)

    records_returned = "YES" if status is not None and 200 <= status < 300 and total_count > 0 else "NO"

    rows.append({
        "gene": candidate["gene"],
        "rsid": candidate["rsid"],
        "gtex_variant_id": candidate["gtex_variant_id"],
        "dataset_id": dataset,
        "resource": "GTEx",
        "query_type": "mapped_variant_single_tissue_eqtl_availability",
        "endpoint_used": endpoint,
        "http_status": "" if status is None else str(status),
        "content_type": content_type,
        "records_returned": records_returned,
        "page_record_count": page_count,
        "total_record_count": total_count,
        "payload_structure": payload_structure,
        "payload_keys": payload_keys,
        "candidate_query_executed": "YES",
        "candidate_level_evidence_retrieved": "RECORD_COUNT_METADATA_ONLY" if records_returned == "YES" else "NO_RECORD_RETURNED",
        "eqtl_interpretation_performed": "NO",
        "tissue_specific_interpretation_performed": "NO",
        "causality_claim_made": "NO",
        "functional_validation_claim_made": "NO",
        "selection_or_population_claim_made": "NO",
        "ancestry_effect_claim_made": "NO",
        "error_message": error_message,
    })

    raw_summary[dataset] = {
        "endpoint": endpoint,
        "http_status": "" if status is None else str(status),
        "records_returned": records_returned,
        "page_record_count": page_count,
        "total_record_count": total_count,
        "payload_structure": payload_structure,
        "payload_keys": payload_keys,
    }

results_path = OUTDIR / "phaseD_gtex_mapped_variant_query_results_v1.tsv"
with results_path.open("w", encoding="utf-8", newline="") as handle:
    fieldnames = list(rows[0].keys())
    writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)

datasets_with_records = sum(1 for r in rows if r["records_returned"] == "YES")

summary_path = OUTDIR / "phaseD_gtex_mapped_variant_query_summary_v1.tsv"
summary_rows = [
    {"metric": "candidate_rsids_queried", "value": 1},
    {"metric": "candidate_rsid", "value": candidate["rsid"]},
    {"metric": "gtex_variant_id", "value": candidate["gtex_variant_id"]},
    {"metric": "datasets_queried", "value": len(datasets)},
    {"metric": "datasets_with_records_returned", "value": datasets_with_records},
    {"metric": "unmapped_rsids_excluded", "value": "rs924735278;rs376205580;rs141807543"},
    {"metric": "coordinate_query_used", "value": "NO"},
    {"metric": "ld_proxy_query_used", "value": "NO"},
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

raw_summary_path = OUTDIR / "phaseD_gtex_mapped_variant_query_raw_payload_summary_v1.json"
raw_summary_path.write_text(json.dumps(raw_summary, indent=2) + "\n", encoding="utf-8")

payload = {
    "status": "PASS",
    "operation": "Controlled GTEx mapped-variant availability query for Phase D candidate rs1131877.",
    "candidate_rsids_queried": 1,
    "candidate_rsid": candidate["rsid"],
    "gtex_variant_id": candidate["gtex_variant_id"],
    "datasets_queried": len(datasets),
    "datasets_with_records_returned": datasets_with_records,
    "unmapped_rsids_excluded": ["rs924735278", "rs376205580", "rs141807543"],
    "coordinate_query_used": "NO",
    "ld_proxy_query_used": "NO",
    "eqtl_interpretation_performed": "NO",
    "tissue_specific_interpretation_performed": "NO",
    "new_biological_claims_added": "NO",
    "causality_claim_authorised": "NO",
    "functional_validation_claim_authorised": "NO",
    "selection_or_population_claim_authorised": "NO",
    "broad_biological_interpretation_authorised": "NO",
    "outputs": {
        "results": str(results_path),
        "summary": str(summary_path),
        "raw_payload_summary": str(raw_summary_path),
    },
    "boundary": (
        "This checkpoint records GTEx mapped-variant endpoint availability and record counts only. "
        "It does not interpret eQTLs, tissue specificity, causality, function, selection, ancestry effects "
        "or population-specific biology."
    ),
}

OUTJSON.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
print(json.dumps(payload, indent=2))
