#!/usr/bin/env python3
"""Create Phase D query pass 1 candidate-resource query matrix.

This checkpoint creates a planned query matrix only. It does not execute
candidate-level external queries, retrieve regulatory evidence, perform
biological interpretation, validate function, infer selection, assign ancestry
effects or authorise broad biological interpretation.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path

CANDIDATES = Path("11_results/phaseD/regulatory_functional_evidence/candidate_input_inventory/phaseD_candidate_input_inventory_v1.tsv")
PROTOCOL = Path("11_results/phaseD/regulatory_functional_evidence/query_protocol/phaseD_resource_specific_candidate_query_protocol_v1.tsv")
OUTDIR = Path("11_results/phaseD/regulatory_functional_evidence/query_pass1")
OUTJSON = Path("results/qc/phaseD_query_pass1_matrix_v1.json")

OUTDIR.mkdir(parents=True, exist_ok=True)
OUTJSON.parent.mkdir(parents=True, exist_ok=True)

failures = []
if not CANDIDATES.exists():
    failures.append(f"Missing candidate inventory: {CANDIDATES}")
if not PROTOCOL.exists():
    failures.append(f"Missing query protocol: {PROTOCOL}")

candidate_rows = []
protocol_rows = []

if CANDIDATES.exists():
    with CANDIDATES.open("r", encoding="utf-8", newline="") as handle:
        candidate_rows = list(csv.DictReader(handle, delimiter="\t"))

if PROTOCOL.exists():
    with PROTOCOL.open("r", encoding="utf-8", newline="") as handle:
        protocol_rows = list(csv.DictReader(handle, delimiter="\t"))

if not candidate_rows:
    failures.append("No candidate rows loaded.")
if not protocol_rows:
    failures.append("No protocol rows loaded.")

query_rows = []
for cand in candidate_rows:
    rsid = cand.get("mapped_identity", "").strip()
    gene = cand.get("gene", "").strip()
    allele_key = cand.get("allele_key", "").strip()
    phaseD_index = cand.get("phaseD_candidate_index", "").strip()
    caution = cand.get("candidate_specific_cautions", "").strip()

    for prot in protocol_rows:
        resource = prot.get("resource", "").strip()
        primary = prot.get("primary_query_input", "").strip()
        secondary = prot.get("secondary_query_input", "").strip()

        if primary == "rsID":
            planned_primary_query_value = rsid
        elif primary == "gene_symbol":
            planned_primary_query_value = gene
        elif primary == "coordinate_with_build_check":
            planned_primary_query_value = allele_key
        else:
            planned_primary_query_value = "PROTOCOL_REVIEW_REQUIRED"

        if secondary == "gene_symbol":
            planned_secondary_query_value = gene
        elif secondary == "coordinate_with_build_check":
            planned_secondary_query_value = allele_key
        elif secondary == "rsID_if_supported":
            planned_secondary_query_value = rsid
        elif secondary in {"NA", ""}:
            planned_secondary_query_value = "NA"
        else:
            planned_secondary_query_value = secondary

        direct_query_allowed = "YES" if primary == "rsID" else "YES_WITH_BUILD_CHECK" if primary == "coordinate_with_build_check" else "REVIEW"
        if resource == "LDlink_or_1000G_proxy":
            direct_query_allowed = "NO_PROXY_METADATA_ONLY"

        query_rows.append(
            {
                "phaseD_candidate_index": phaseD_index,
                "allele_key": allele_key,
                "gene": gene,
                "rsid": rsid,
                "resource": resource,
                "query_unit": prot.get("query_unit", "").strip(),
                "primary_query_input": primary,
                "planned_primary_query_value": planned_primary_query_value,
                "secondary_query_input": secondary,
                "planned_secondary_query_value": planned_secondary_query_value,
                "direct_hit_rule": prot.get("direct_hit_rule", "").strip(),
                "proxy_rule": prot.get("proxy_rule", "").strip(),
                "evidence_output_type": prot.get("evidence_output_type", "").strip(),
                "direct_query_allowed_pass1": direct_query_allowed,
                "candidate_query_executed": "NO",
                "candidate_level_evidence_retrieved": "NO",
                "interpretation_performed": "NO",
                "candidate_specific_cautions": caution,
                "boundary": prot.get("boundary", "").strip(),
            }
        )

query_matrix_path = OUTDIR / "phaseD_query_pass1_candidate_resource_matrix_v1.tsv"
with query_matrix_path.open("w", encoding="utf-8", newline="") as handle:
    fieldnames = list(query_rows[0].keys()) if query_rows else [
        "phaseD_candidate_index", "allele_key", "gene", "rsid", "resource"
    ]
    writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(query_rows)

resource_counts = {}
candidate_counts = {}
for row in query_rows:
    resource_counts[row["resource"]] = resource_counts.get(row["resource"], 0) + 1
    candidate_counts[row["rsid"]] = candidate_counts.get(row["rsid"], 0) + 1

resource_counts_path = OUTDIR / "phaseD_query_pass1_resource_counts_v1.tsv"
with resource_counts_path.open("w", encoding="utf-8", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=["resource", "planned_query_rows"], delimiter="\t", lineterminator="\n")
    writer.writeheader()
    for resource, count in sorted(resource_counts.items()):
        writer.writerow({"resource": resource, "planned_query_rows": count})

candidate_counts_path = OUTDIR / "phaseD_query_pass1_candidate_counts_v1.tsv"
with candidate_counts_path.open("w", encoding="utf-8", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=["rsid", "planned_resource_rows"], delimiter="\t", lineterminator="\n")
    writer.writeheader()
    for rsid, count in sorted(candidate_counts.items()):
        writer.writerow({"rsid": rsid, "planned_resource_rows": count})

summary_path = OUTDIR / "phaseD_query_pass1_matrix_summary_v1.tsv"
summary_rows = [
    {"metric": "candidate_rows", "value": len(candidate_rows)},
    {"metric": "resources_in_protocol", "value": len(protocol_rows)},
    {"metric": "planned_candidate_resource_rows", "value": len(query_rows)},
    {"metric": "candidate_query_executed", "value": "NO"},
    {"metric": "candidate_level_evidence_retrieved", "value": "NO"},
    {"metric": "interpretation_performed", "value": "NO"},
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
    "status": "PASS" if not failures else "FAIL",
    "operation": "Create Phase D query pass 1 candidate-resource query matrix.",
    "candidate_rows": len(candidate_rows),
    "resources_in_protocol": len(protocol_rows),
    "planned_candidate_resource_rows": len(query_rows),
    "candidate_query_executed": "NO",
    "candidate_level_evidence_retrieved": "NO",
    "interpretation_performed": "NO",
    "new_analysis_performed": "NO",
    "new_biological_claims_added": "NO",
    "selection_or_population_claim_authorised": "NO",
    "broad_biological_interpretation_authorised": "NO",
    "outputs": {
        "query_matrix": str(query_matrix_path),
        "resource_counts": str(resource_counts_path),
        "candidate_counts": str(candidate_counts_path),
        "summary": str(summary_path),
    },
    "failures": failures,
    "boundary": (
        "This checkpoint creates a planned query matrix only. It does not execute candidate-level external queries, "
        "retrieve regulatory evidence, perform biological interpretation, validate function, infer selection, assign ancestry effects "
        "or authorise broad biological interpretation."
    ),
}

OUTJSON.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
print(json.dumps(payload, indent=2))

if payload["status"] != "PASS":
    raise SystemExit(1)
