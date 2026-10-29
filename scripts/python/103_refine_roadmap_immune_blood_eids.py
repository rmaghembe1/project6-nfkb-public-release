#!/usr/bin/env python3
"""Refine Roadmap immune/blood EID selection using structured metadata.

This checkpoint refines immune/blood-relevant Roadmap EIDs using structured
metadata columns from EID_metadata.tab rather than broad substring matching.

It does not select final ChromHMM files, download chromatin-state files, perform
candidate overlap, retrieve candidate chromatin-state evidence, validate function,
infer causality, infer selection, assign ancestry effects or authorise broad
biological interpretation.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path

INFILE = Path("11_results/phaseD/regulatory_functional_evidence/query_pass1/roadmap_chromhmm/eid_metadata/phaseD_roadmap_eid_metadata_all_rows_v1.tsv")
OUTDIR = Path("11_results/phaseD/regulatory_functional_evidence/query_pass1/roadmap_chromhmm/eid_refinement")
OUTJSON = Path("results/qc/phaseD_roadmap_immune_blood_eid_refinement_v1.json")

OUTDIR.mkdir(parents=True, exist_ok=True)
OUTJSON.parent.mkdir(parents=True, exist_ok=True)

# Based on parsed EID_metadata.tab structure observed in prior checkpoint:
# col_1 = EID
# col_2 = broad group
# col_4 = short mnemonic
# col_5 = tissue/cell label
# col_6 = harmonised cell/tissue name
# col_7 = tissue/category
# col_8 = sample type
# col_9 = age/source details
# col_13 = consolidated/standard mark status
PRIMARY_STRUCTURED_CATEGORIES = {"BLOOD", "SPLEEN", "THYMUS"}
PRIMARY_GROUP_TERMS = {"HSC & B-cell", "Blood & T-cell"}
IMMUNE_LABEL_TERMS = [
    "monocyte", "neutrophil", "b cell", "b-cell", "cd19", "t cell", "t-cell",
    "cd3", "cd4", "cd8", "treg", "th17", "natural killer", "cd56",
    "hematopoietic", "haematopoietic", "cd34", "mononuclear", "lymphoblastoid",
    "leukemia", "k562", "spleen", "thymus"
]

with INFILE.open("r", encoding="utf-8", newline="") as handle:
    rows = list(csv.DictReader(handle, delimiter="\t"))

refined_rows = []
for row in rows:
    eid = row.get("col_1", "")
    group = row.get("col_2", "")
    mnemonic = row.get("col_4", "")
    label = row.get("col_5", "")
    harmonised = row.get("col_6", "")
    category = row.get("col_7", "")
    sample_type = row.get("col_8", "")
    standard_status = row.get("col_13", "")

    text = " ".join([group, mnemonic, label, harmonised, category, sample_type]).lower()

    structured_category_match = category in PRIMARY_STRUCTURED_CATEGORIES
    group_match = group in PRIMARY_GROUP_TERMS
    immune_label_hits = sorted({term for term in IMMUNE_LABEL_TERMS if term in text})

    # Conservative selection: structured category/group support OR clear immune label.
    selected = structured_category_match or group_match or bool(immune_label_hits)

    # Tiering for future file selection.
    if category == "BLOOD" or group in PRIMARY_GROUP_TERMS:
        tier = "tier1_blood_immune"
    elif category in {"SPLEEN", "THYMUS"}:
        tier = "tier2_lymphoid_tissue"
    elif immune_label_hits:
        tier = "tier3_immune_label_only"
    else:
        tier = "not_selected"

    # Prefer consolidated entries where available, but do not exclude split datasets yet.
    future_priority = "PRIORITISE" if tier == "tier1_blood_immune" and standard_status in {"C", "SD"} else (
        "SECONDARY" if tier in {"tier2_lymphoid_tissue", "tier3_immune_label_only"} else "NOT_SELECTED"
    )

    refined_rows.append({
        "eid": eid,
        "broad_group": group,
        "mnemonic": mnemonic,
        "label": label,
        "harmonised_name": harmonised,
        "structured_category": category,
        "sample_type": sample_type,
        "standard_status": standard_status,
        "structured_category_match": "YES" if structured_category_match else "NO",
        "broad_group_match": "YES" if group_match else "NO",
        "immune_label_hits": ";".join(immune_label_hits),
        "refined_immune_blood_selection": "YES" if selected else "NO",
        "selection_tier": tier,
        "future_file_selection_priority": future_priority,
        "specific_chromhmm_file_selected": "NO",
        "candidate_overlap_performed": "NO",
        "candidate_chromatin_state_evidence_retrieved": "NO",
        "biological_interpretation_performed": "NO",
    })

selected_rows = [r for r in refined_rows if r["refined_immune_blood_selection"] == "YES"]
priority_rows = [r for r in refined_rows if r["future_file_selection_priority"] == "PRIORITISE"]

all_path = OUTDIR / "phaseD_roadmap_immune_blood_eid_refinement_all_rows_v1.tsv"
selected_path = OUTDIR / "phaseD_roadmap_immune_blood_eid_refined_selected_v1.tsv"
priority_path = OUTDIR / "phaseD_roadmap_immune_blood_eid_priority_for_file_selection_v1.tsv"
summary_path = OUTDIR / "phaseD_roadmap_immune_blood_eid_refinement_summary_v1.tsv"

fieldnames = list(refined_rows[0].keys()) if refined_rows else []
for path, data in [(all_path, refined_rows), (selected_path, selected_rows), (priority_path, priority_rows)]:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(data)

tier_counts = {}
for r in selected_rows:
    tier_counts[r["selection_tier"]] = tier_counts.get(r["selection_tier"], 0) + 1

summary_rows = [
    {"metric": "input_rows", "value": len(rows)},
    {"metric": "refined_selected_eids", "value": len(selected_rows)},
    {"metric": "priority_file_selection_eids", "value": len(priority_rows)},
    {"metric": "tier1_blood_immune_eids", "value": tier_counts.get("tier1_blood_immune", 0)},
    {"metric": "tier2_lymphoid_tissue_eids", "value": tier_counts.get("tier2_lymphoid_tissue", 0)},
    {"metric": "tier3_immune_label_only_eids", "value": tier_counts.get("tier3_immune_label_only", 0)},
    {"metric": "specific_chromhmm_file_selected", "value": "NO"},
    {"metric": "build_confirmed_for_file", "value": "NO"},
    {"metric": "candidate_overlap_performed", "value": "NO"},
    {"metric": "candidate_chromatin_state_evidence_retrieved", "value": "NO"},
    {"metric": "biological_interpretation_performed", "value": "NO"},
    {"metric": "functional_validation_claim_authorised", "value": "NO"},
    {"metric": "broad_biological_interpretation_authorised", "value": "NO"},
]
with summary_path.open("w", encoding="utf-8", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=["metric", "value"], delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(summary_rows)

payload = {
    "status": "PASS",
    "operation": "Refine Roadmap immune/blood EID selection using structured metadata.",
    "input_rows": len(rows),
    "refined_selected_eids": len(selected_rows),
    "priority_file_selection_eids": len(priority_rows),
    "tier1_blood_immune_eids": tier_counts.get("tier1_blood_immune", 0),
    "tier2_lymphoid_tissue_eids": tier_counts.get("tier2_lymphoid_tissue", 0),
    "tier3_immune_label_only_eids": tier_counts.get("tier3_immune_label_only", 0),
    "specific_chromhmm_file_selected": "NO",
    "build_confirmed_for_file": "NO",
    "candidate_overlap_performed": "NO",
    "candidate_chromatin_state_evidence_retrieved": "NO",
    "biological_interpretation_performed": "NO",
    "functional_validation_claim_authorised": "NO",
    "broad_biological_interpretation_authorised": "NO",
    "outputs": {
        "all_rows": str(all_path),
        "selected_eids": str(selected_path),
        "priority_eids": str(priority_path),
        "summary": str(summary_path),
    },
    "boundary": (
        "This checkpoint refines Roadmap immune/blood EID selection only. It does not select final "
        "ChromHMM files, perform candidate overlap, retrieve candidate evidence, validate function "
        "or authorise biological interpretation."
    ),
}

OUTJSON.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
print(json.dumps(payload, indent=2))
