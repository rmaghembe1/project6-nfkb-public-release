#!/usr/bin/env python3
"""Confirm Phase D candidate coordinate builds using Ensembl rsID lookup.

This checkpoint performs coordinate-build confirmation only. It queries Ensembl
GRCh38 and GRCh37 REST variation endpoints by rsID, records mappings and compares
them with project coordinate labels. It does not query regulatory resources,
perform interval overlap, validate function, infer selection, assign ancestry
effects or authorise broad biological interpretation.
"""

from __future__ import annotations

import csv
import json
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

CANDIDATES = [
    {"candidate": "1", "gene": "TRAF3", "rsid": "rs1131877", "project_coordinate": "chr14:102875712:T:C"},
    {"candidate": "2", "gene": "BIRC3", "rsid": "rs924735278", "project_coordinate": "chr11:102324736:G:T"},
    {"candidate": "3", "gene": "TNFAIP3", "rsid": "rs376205580", "project_coordinate": "chr6:137874871:A:G"},
    {"candidate": "4", "gene": "TNFAIP3", "rsid": "rs141807543", "project_coordinate": "chr6:137875820:A:C"},
]

OUTDIR = Path("11_results/phaseD/regulatory_functional_evidence/query_pass1/coordinate_build_confirmation/query")
OUTJSON = Path("results/qc/phaseD_coordinate_build_confirmation_ensembl_v1.json")
OUTDIR.mkdir(parents=True, exist_ok=True)
OUTJSON.parent.mkdir(parents=True, exist_ok=True)

SERVERS = [
    {"build": "GRCh38", "server": "https://rest.ensembl.org"},
    {"build": "GRCh37", "server": "https://grch37.rest.ensembl.org"},
]

def parse_project_coordinate(coord: str) -> dict[str, str]:
    chrom_pos, ref, alt = coord.rsplit(":", 2)
    chrom, pos = chrom_pos.split(":")
    return {
        "project_chromosome": chrom,
        "project_position": pos,
        "project_ref": ref,
        "project_alt": alt,
    }

def request_json(url: str, timeout: int = 45) -> tuple[int | None, object | None, str]:
    req = urllib.request.Request(
        url,
        headers={
            "User-Agent": "Project6-PhaseD-coordinate-build-confirmation/1.0",
            "Accept": "application/json",
            "Content-Type": "application/json",
        },
    )
    try:
        with urllib.request.urlopen(req, timeout=timeout) as resp:
            text = resp.read().decode("utf-8", errors="replace")
            return resp.status, json.loads(text), ""
    except urllib.error.HTTPError as exc:
        text = ""
        try:
            text = exc.read(1000).decode("utf-8", errors="replace")
        except Exception:
            pass
        return exc.code, None, f"HTTPError: {exc}; {text[:300]}"
    except Exception as exc:
        return None, None, f"{type(exc).__name__}: {exc}"

def normalize_chrom(value: str) -> str:
    return value.replace("chr", "")

def extract_mappings(payload: object | None) -> list[dict[str, str]]:
    mappings = []
    if not isinstance(payload, dict):
        return mappings

    # Ensembl variation endpoint often exposes mappings under "mappings".
    raw_mappings = payload.get("mappings", [])
    if isinstance(raw_mappings, list):
        for m in raw_mappings:
            if not isinstance(m, dict):
                continue
            mappings.append({
                "mapping_source": "mappings",
                "seq_region_name": str(m.get("seq_region_name", "")),
                "start": str(m.get("start", "")),
                "end": str(m.get("end", "")),
                "strand": str(m.get("strand", "")),
                "allele_string": str(m.get("allele_string", "")),
                "location": str(m.get("location", "")),
            })

    # Some responses may include synonyms/source metadata but no mappings.
    return mappings

rows = []
summary_by_candidate = []

for cand in CANDIDATES:
    project = parse_project_coordinate(cand["project_coordinate"])

    build_hits = []
    for server in SERVERS:
        rsid_encoded = urllib.parse.quote(cand["rsid"], safe="")
        url = f"{server['server']}/variation/homo_sapiens/{rsid_encoded}?content-type=application/json"
        status, payload, error = request_json(url)
        mappings = extract_mappings(payload)

        if not mappings:
            rows.append({
                **cand,
                **project,
                "queried_build": server["build"],
                "endpoint_used": url,
                "http_status": "" if status is None else str(status),
                "mapping_returned": "NO",
                "mapped_chromosome": "",
                "mapped_position": "",
                "mapped_end": "",
                "allele_string": "",
                "location": "",
                "project_coordinate_match": "NO",
                "build_decision_candidate": "not_resolved",
                "coordinate_interpretation_performed": "NO",
                "regulatory_overlap_performed": "NO",
                "biological_interpretation_performed": "NO",
                "error_message": error,
            })
        else:
            for m in mappings:
                mapped_chrom = m["seq_region_name"]
                mapped_pos = m["start"]
                chrom_match = normalize_chrom(project["project_chromosome"]) == normalize_chrom(mapped_chrom)
                pos_match = project["project_position"] == mapped_pos
                coord_match = chrom_match and pos_match
                if coord_match:
                    build_hits.append(server["build"])

                rows.append({
                    **cand,
                    **project,
                    "queried_build": server["build"],
                    "endpoint_used": url,
                    "http_status": "" if status is None else str(status),
                    "mapping_returned": "YES",
                    "mapped_chromosome": mapped_chrom,
                    "mapped_position": mapped_pos,
                    "mapped_end": m["end"],
                    "allele_string": m["allele_string"],
                    "location": m["location"],
                    "project_coordinate_match": "YES" if coord_match else "NO",
                    "build_decision_candidate": f"candidate_{server['build']}" if coord_match else "coordinate_mismatch",
                    "coordinate_interpretation_performed": "NO",
                    "regulatory_overlap_performed": "NO",
                    "biological_interpretation_performed": "NO",
                    "error_message": error,
                })
        time.sleep(0.2)

    if len(build_hits) == 1:
        decision = f"confirmed_{build_hits[0]}"
    elif len(build_hits) > 1:
        decision = "build_conflict"
    else:
        decision = "not_resolved"

    summary_by_candidate.append({
        "candidate": cand["candidate"],
        "gene": cand["gene"],
        "rsid": cand["rsid"],
        "project_coordinate": cand["project_coordinate"],
        "matching_builds": ";".join(build_hits) if build_hits else "NONE",
        "build_decision": decision,
        "coordinate_overlap_authorised": "YES" if decision in {"confirmed_GRCh38", "confirmed_GRCh37"} else "NO",
        "biological_interpretation_performed": "NO",
    })

results_path = OUTDIR / "phaseD_coordinate_build_confirmation_ensembl_results_v1.tsv"
with results_path.open("w", encoding="utf-8", newline="") as handle:
    fieldnames = list(rows[0].keys())
    writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)

candidate_summary_path = OUTDIR / "phaseD_coordinate_build_confirmation_candidate_summary_v1.tsv"
with candidate_summary_path.open("w", encoding="utf-8", newline="") as handle:
    fieldnames = list(summary_by_candidate[0].keys())
    writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(summary_by_candidate)

confirmed = sum(1 for r in summary_by_candidate if r["build_decision"] in {"confirmed_GRCh38", "confirmed_GRCh37"})
not_resolved = sum(1 for r in summary_by_candidate if r["build_decision"] == "not_resolved")
conflicts = sum(1 for r in summary_by_candidate if r["build_decision"] == "build_conflict")

summary_path = OUTDIR / "phaseD_coordinate_build_confirmation_summary_v1.tsv"
summary_rows = [
    {"metric": "candidate_rsids_queried", "value": len(CANDIDATES)},
    {"metric": "builds_tested", "value": "GRCh38;GRCh37"},
    {"metric": "candidates_confirmed_one_build", "value": confirmed},
    {"metric": "candidates_not_resolved", "value": not_resolved},
    {"metric": "candidates_with_build_conflict", "value": conflicts},
    {"metric": "coordinate_overlap_performed", "value": "NO"},
    {"metric": "regulatory_resource_query_performed", "value": "NO"},
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
    "operation": "Confirm Phase D candidate coordinate builds using Ensembl rsID lookup.",
    "candidate_rsids_queried": len(CANDIDATES),
    "builds_tested": ["GRCh38", "GRCh37"],
    "candidates_confirmed_one_build": confirmed,
    "candidates_not_resolved": not_resolved,
    "candidates_with_build_conflict": conflicts,
    "coordinate_overlap_performed": "NO",
    "regulatory_resource_query_performed": "NO",
    "biological_interpretation_performed": "NO",
    "functional_validation_claim_authorised": "NO",
    "selection_or_population_claim_authorised": "NO",
    "broad_biological_interpretation_authorised": "NO",
    "outputs": {
        "results": str(results_path),
        "candidate_summary": str(candidate_summary_path),
        "summary": str(summary_path),
    },
    "boundary": (
        "This checkpoint performs coordinate-build confirmation only. It does not query regulatory resources, "
        "perform interval overlap, validate function, infer selection, assign ancestry effects or authorise broad biological interpretation."
    ),
}

OUTJSON.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
print(json.dumps(payload, indent=2))
