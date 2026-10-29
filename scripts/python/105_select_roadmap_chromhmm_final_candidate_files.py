#!/usr/bin/env python3
"""Select final Roadmap ChromHMM candidate files under committed policy.

This checkpoint selects final candidate Roadmap ChromHMM files for the 27
priority immune/blood EIDs using the committed file-selection policy:
imputed12marks, hg38lift, mnemonics.bed.gz preferred, segments.bed.gz fallback.

It does not download full chromatin-state segmentation files, perform candidate
overlap, retrieve candidate chromatin-state evidence, validate function, infer
causality, infer selection, assign ancestry effects or authorise broad biological
interpretation.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path

DISCOVERY = Path("11_results/phaseD/regulatory_functional_evidence/query_pass1/roadmap_chromhmm/file_discovery/phaseD_roadmap_chromhmm_priority_eid_file_discovery_mapping_v1.tsv")
OUTDIR = Path("11_results/phaseD/regulatory_functional_evidence/query_pass1/roadmap_chromhmm/final_file_mapping")
OUTJSON = Path("results/qc/phaseD_roadmap_chromhmm_final_file_mapping_v1.json")

OUTDIR.mkdir(parents=True, exist_ok=True)
OUTJSON.parent.mkdir(parents=True, exist_ok=True)

PREFERRED_TERMS = ["imputed12marks", "hg38lift", "mnemonics.bed.gz"]
FALLBACK_TERMS = ["imputed12marks", "hg38lift", "segments.bed.gz"]

EXCLUDE_TERMS = [
    "dense",
    "expanded",
    "stateno",
    "conservation_overlap",
    "coremarks",
    "core_k27ac",
]

def has_all_terms(text: str, terms: list[str]) -> bool:
    lower = text.lower()
    return all(term.lower() in lower for term in terms)

def has_excluded(text: str) -> bool:
    lower = text.lower()
    return any(term.lower() in lower for term in EXCLUDE_TERMS)

with DISCOVERY.open("r", encoding="utf-8", newline="") as handle:
    rows = list(csv.DictReader(handle, delimiter="\t"))

by_eid: dict[str, list[dict[str, str]]] = {}
for row in rows:
    eid = row.get("eid", "")
    if not eid:
        continue
    by_eid.setdefault(eid, []).append(row)

selected_rows = []
candidate_review_rows = []

for eid, eid_rows in sorted(by_eid.items()):
    observed = [r for r in eid_rows if r.get("file_mapping_status") == "OBSERVED_CANDIDATE_FILE_NOT_FINAL" and r.get("observed_file_url")]
    for r in observed:
        url = r.get("observed_file_url", "")
        fname = r.get("observed_file_name", "")
        preferred = has_all_terms(url, PREFERRED_TERMS) and not has_excluded(url)
        fallback = has_all_terms(url, FALLBACK_TERMS) and not has_excluded(url)

        candidate_review_rows.append({
            "eid": eid,
            "observed_file_name": fname,
            "observed_file_url": url,
            "matches_preferred_policy": "YES" if preferred else "NO",
            "matches_fallback_policy": "YES" if fallback else "NO",
            "excluded_by_policy": "YES" if has_excluded(url) else "NO",
            "specific_chromhmm_file_selected": "NO",
            "candidate_overlap_performed": "NO",
            "candidate_chromatin_state_evidence_retrieved": "NO",
            "biological_interpretation_performed": "NO",
        })

    preferred_hits = [
        r for r in observed
        if has_all_terms(r.get("observed_file_url", ""), PREFERRED_TERMS)
        and not has_excluded(r.get("observed_file_url", ""))
    ]
    fallback_hits = [
        r for r in observed
        if has_all_terms(r.get("observed_file_url", ""), FALLBACK_TERMS)
        and not has_excluded(r.get("observed_file_url", ""))
    ]

    if preferred_hits:
        chosen = sorted(preferred_hits, key=lambda x: x.get("observed_file_url", ""))[0]
        selection_status = "FINAL_CANDIDATE_FILE_SELECTED_BY_POLICY"
        selection_class = "preferred_imputed12marks_hg38lift_mnemonics"
    elif fallback_hits:
        chosen = sorted(fallback_hits, key=lambda x: x.get("observed_file_url", ""))[0]
        selection_status = "FINAL_CANDIDATE_FILE_SELECTED_BY_FALLBACK_POLICY"
        selection_class = "fallback_imputed12marks_hg38lift_segments"
    else:
        chosen = observed[0] if observed else {}
        selection_status = "NO_POLICY_MATCH_REQUIRES_MANUAL_REVIEW"
        selection_class = "not_selected"

    selected_rows.append({
        "eid": eid,
        "broad_group": chosen.get("broad_group", ""),
        "mnemonic": chosen.get("mnemonic", ""),
        "label": chosen.get("label", ""),
        "harmonised_name": chosen.get("harmonised_name", ""),
        "structured_category": chosen.get("structured_category", ""),
        "sample_type": chosen.get("sample_type", ""),
        "standard_status": chosen.get("standard_status", ""),
        "selected_file_name": chosen.get("observed_file_name", "") if "SELECTED" in selection_status else "",
        "selected_file_url": chosen.get("observed_file_url", "") if "SELECTED" in selection_status else "",
        "selection_status": selection_status,
        "selection_class": selection_class,
        "selected_model": "imputed12marks" if "SELECTED" in selection_status else "",
        "selected_build": "hg38lift" if "SELECTED" in selection_status else "",
        "selected_output_type": "mnemonics.bed.gz" if "mnemonics" in selection_class else ("segments.bed.gz" if "segments" in selection_class else ""),
        "apparent_build": "GRCh38_lifted_hg38lift" if "SELECTED" in selection_status else "UNCONFIRMED",
        "build_compatible_with_candidate_coordinates": "YES" if "SELECTED" in selection_status else "NO",
        "specific_chromhmm_file_selected": "YES" if "SELECTED" in selection_status else "NO",
        "download_performed": "NO",
        "candidate_overlap_performed": "NO",
        "candidate_chromatin_state_evidence_retrieved": "NO",
        "biological_interpretation_performed": "NO",
    })

selected_path = OUTDIR / "phaseD_roadmap_chromhmm_final_candidate_files_v1.tsv"
review_path = OUTDIR / "phaseD_roadmap_chromhmm_policy_candidate_file_review_v1.tsv"
summary_path = OUTDIR / "phaseD_roadmap_chromhmm_final_file_mapping_summary_v1.tsv"

with selected_path.open("w", encoding="utf-8", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=list(selected_rows[0].keys()), delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(selected_rows)

with review_path.open("w", encoding="utf-8", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=list(candidate_review_rows[0].keys()) if candidate_review_rows else [
        "eid", "observed_file_name", "observed_file_url", "matches_preferred_policy",
        "matches_fallback_policy", "excluded_by_policy", "specific_chromhmm_file_selected",
        "candidate_overlap_performed", "candidate_chromatin_state_evidence_retrieved",
        "biological_interpretation_performed"
    ], delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(candidate_review_rows)

selected_count = sum(1 for r in selected_rows if r["specific_chromhmm_file_selected"] == "YES")
preferred_count = sum(1 for r in selected_rows if r["selection_class"] == "preferred_imputed12marks_hg38lift_mnemonics")
fallback_count = sum(1 for r in selected_rows if r["selection_class"] == "fallback_imputed12marks_hg38lift_segments")
manual_count = sum(1 for r in selected_rows if r["selection_status"] == "NO_POLICY_MATCH_REQUIRES_MANUAL_REVIEW")

summary_rows = [
    {"metric": "priority_eids_evaluated", "value": len(selected_rows)},
    {"metric": "final_candidate_files_selected", "value": selected_count},
    {"metric": "preferred_mnemonics_files_selected", "value": preferred_count},
    {"metric": "fallback_segments_files_selected", "value": fallback_count},
    {"metric": "manual_review_required_eids", "value": manual_count},
    {"metric": "all_selected_files_build_compatible_with_grch38_candidates", "value": "YES" if selected_count == len(selected_rows) else "NO"},
    {"metric": "download_performed", "value": "NO"},
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
    "operation": "Select final Roadmap ChromHMM candidate files under committed policy.",
    "priority_eids_evaluated": len(selected_rows),
    "final_candidate_files_selected": selected_count,
    "preferred_mnemonics_files_selected": preferred_count,
    "fallback_segments_files_selected": fallback_count,
    "manual_review_required_eids": manual_count,
    "download_performed": "NO",
    "candidate_overlap_performed": "NO",
    "candidate_chromatin_state_evidence_retrieved": "NO",
    "biological_interpretation_performed": "NO",
    "functional_validation_claim_authorised": "NO",
    "broad_biological_interpretation_authorised": "NO",
    "outputs": {
        "final_candidate_files": str(selected_path),
        "policy_review": str(review_path),
        "summary": str(summary_path),
    },
    "boundary": (
        "This checkpoint selects final Roadmap ChromHMM candidate files under policy only. It does not "
        "download full files, perform candidate overlap, retrieve candidate evidence, validate function or "
        "authorise biological interpretation."
    ),
}

OUTJSON.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
print(json.dumps(payload, indent=2))
