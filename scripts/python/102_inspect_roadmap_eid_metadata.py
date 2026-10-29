#!/usr/bin/env python3
"""Inspect Roadmap EID metadata before ChromHMM file selection.

This checkpoint retrieves/inspects Roadmap EID_metadata.tab and identifies
immune/blood-relevant Roadmap sample IDs by metadata text matching.

It does not select final ChromHMM files, download chromatin-state files, perform
candidate overlap, retrieve candidate chromatin-state evidence, validate function,
infer causality, infer selection, assign ancestry effects or authorise broad
biological interpretation.
"""

from __future__ import annotations

import csv
import hashlib
import json
import urllib.request
from pathlib import Path

OUTDIR = Path("11_results/phaseD/regulatory_functional_evidence/query_pass1/roadmap_chromhmm/eid_metadata")
OUTJSON = Path("results/qc/phaseD_roadmap_eid_metadata_inspection_v1.json")

OUTDIR.mkdir(parents=True, exist_ok=True)
OUTJSON.parent.mkdir(parents=True, exist_ok=True)

URL = "https://egg2.wustl.edu/roadmap/data/byFileType/metadata/EID_metadata.tab"
RAW = OUTDIR / "EID_metadata.tab"
ALL_TSV = OUTDIR / "phaseD_roadmap_eid_metadata_all_rows_v1.tsv"
IMMUNE_TSV = OUTDIR / "phaseD_roadmap_eid_metadata_immune_blood_candidate_rows_v1.tsv"
SUMMARY_TSV = OUTDIR / "phaseD_roadmap_eid_metadata_inspection_summary_v1.tsv"

IMMUNE_TERMS = [
    "blood", "monocyte", "macrophage", "neutrophil", "t cell", "t-cell", "tcell",
    "b cell", "b-cell", "bcell", "nk", "natural killer", "cd4", "cd8", "pbmc",
    "lymph", "thymus", "spleen", "hematopoietic", "haematopoietic", "immune",
    "leukocyte", "granulocyte", "cd14", "cd19", "cd34", "cd56"
]

def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()

req = urllib.request.Request(
    URL,
    headers={"User-Agent": "Project6-PhaseD-Roadmap-EID-metadata-inspection/1.0"},
)
with urllib.request.urlopen(req, timeout=60) as resp:
    raw = resp.read()
    http_status = resp.status
    content_type = resp.headers.get("content-type", "")

RAW.write_bytes(raw)

text = raw.decode("utf-8", errors="replace")
lines = [line for line in text.splitlines() if line.strip()]

# Try to parse as tabular metadata. Preserve unknown column structure safely.
reader = csv.reader(lines, delimiter="\t")
parsed = list(reader)

header = parsed[0] if parsed else []
rows = parsed[1:] if len(parsed) > 1 else []

max_cols = max([len(header)] + [len(r) for r in rows]) if parsed else 0
fieldnames = [f"col_{i+1}" for i in range(max_cols)]

normalised_rows = []
for row in rows:
    padded = row + [""] * (max_cols - len(row))
    row_text = " ".join(padded).lower()
    hits = sorted({term for term in IMMUNE_TERMS if term in row_text})
    normalised = {fieldnames[i]: padded[i] for i in range(max_cols)}
    normalised.update({
        "immune_blood_relevant_text_flag": "YES" if hits else "NO",
        "immune_blood_terms_detected": ";".join(hits),
        "specific_chromhmm_file_selected": "NO",
        "candidate_overlap_performed": "NO",
        "candidate_chromatin_state_evidence_retrieved": "NO",
        "biological_interpretation_performed": "NO",
    })
    normalised_rows.append(normalised)

immune_rows = [r for r in normalised_rows if r["immune_blood_relevant_text_flag"] == "YES"]

with ALL_TSV.open("w", encoding="utf-8", newline="") as handle:
    out_fields = fieldnames + [
        "immune_blood_relevant_text_flag",
        "immune_blood_terms_detected",
        "specific_chromhmm_file_selected",
        "candidate_overlap_performed",
        "candidate_chromatin_state_evidence_retrieved",
        "biological_interpretation_performed",
    ]
    writer = csv.DictWriter(handle, fieldnames=out_fields, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(normalised_rows)

with IMMUNE_TSV.open("w", encoding="utf-8", newline="") as handle:
    out_fields = fieldnames + [
        "immune_blood_relevant_text_flag",
        "immune_blood_terms_detected",
        "specific_chromhmm_file_selected",
        "candidate_overlap_performed",
        "candidate_chromatin_state_evidence_retrieved",
        "biological_interpretation_performed",
    ]
    writer = csv.DictWriter(handle, fieldnames=out_fields, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(immune_rows)

summary_rows = [
    {"metric": "metadata_url", "value": URL},
    {"metric": "http_status", "value": http_status},
    {"metric": "content_type", "value": content_type},
    {"metric": "raw_file_size_bytes", "value": RAW.stat().st_size},
    {"metric": "raw_file_sha256", "value": sha256(RAW)},
    {"metric": "parsed_rows", "value": len(normalised_rows)},
    {"metric": "parsed_columns_max", "value": max_cols},
    {"metric": "immune_blood_candidate_rows_by_text", "value": len(immune_rows)},
    {"metric": "metadata_map_selected", "value": "YES_AS_EID_METADATA_TAB"},
    {"metric": "specific_chromhmm_file_selected", "value": "NO"},
    {"metric": "build_confirmed_for_file", "value": "NO"},
    {"metric": "candidate_overlap_performed", "value": "NO"},
    {"metric": "candidate_chromatin_state_evidence_retrieved", "value": "NO"},
    {"metric": "biological_interpretation_performed", "value": "NO"},
    {"metric": "functional_validation_claim_authorised", "value": "NO"},
    {"metric": "broad_biological_interpretation_authorised", "value": "NO"},
]
with SUMMARY_TSV.open("w", encoding="utf-8", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=["metric", "value"], delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(summary_rows)

payload = {
    "status": "PASS",
    "operation": "Inspect Roadmap EID_metadata.tab before ChromHMM file selection.",
    "metadata_url": URL,
    "http_status": http_status,
    "raw_file_size_bytes": RAW.stat().st_size,
    "raw_file_sha256": sha256(RAW),
    "parsed_rows": len(normalised_rows),
    "parsed_columns_max": max_cols,
    "immune_blood_candidate_rows_by_text": len(immune_rows),
    "metadata_map_selected": "YES_AS_EID_METADATA_TAB",
    "specific_chromhmm_file_selected": "NO",
    "build_confirmed_for_file": "NO",
    "candidate_overlap_performed": "NO",
    "candidate_chromatin_state_evidence_retrieved": "NO",
    "biological_interpretation_performed": "NO",
    "functional_validation_claim_authorised": "NO",
    "broad_biological_interpretation_authorised": "NO",
    "outputs": {
        "raw_metadata": str(RAW),
        "all_rows": str(ALL_TSV),
        "immune_blood_candidates": str(IMMUNE_TSV),
        "summary": str(SUMMARY_TSV),
    },
    "boundary": (
        "This checkpoint inspects Roadmap EID metadata only. It does not select final ChromHMM files, "
        "perform candidate overlap, retrieve candidate evidence, validate function or authorise biological interpretation."
    ),
}

OUTJSON.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
print(json.dumps(payload, indent=2))
