#!/usr/bin/env python3
"""Execute controlled ENCODE cCRE exact-position overlap.

This checkpoint downloads the selected UCSC hg38 ENCODE3 cCRE BigBed file,
converts it to BED using bigBedToBed, creates exact 1-bp GRCh38 candidate BED
intervals, and performs bedtools intersect.

It records overlap metadata only. It does not interpret regulatory biology,
validate function, infer causality, infer selection, assign ancestry effects
or authorise broad biological interpretation.
"""

from __future__ import annotations

import csv
import gzip
import hashlib
import json
import shutil
import subprocess
import urllib.request
from pathlib import Path

OUTDIR = Path("11_results/phaseD/regulatory_functional_evidence/query_pass1/encode_ccre/overlap_execution")
OUTJSON = Path("results/qc/phaseD_encode_ccre_exact_position_overlap_v1.json")

OUTDIR.mkdir(parents=True, exist_ok=True)
OUTJSON.parent.mkdir(parents=True, exist_ok=True)

BIGBED_URL = "https://hgdownload.soe.ucsc.edu/gbdb/hg38/encode3/ccre/encodeCcreCombined.bb"
BIGBED_PATH = OUTDIR / "encodeCcreCombined.hg38.bigBed"
CCRE_BED_RAW = OUTDIR / "encodeCcreCombined.hg38.raw.bed"
CANDIDATE_BED = OUTDIR / "phaseD_candidates_grch38_exact_1bp_v1.bed"
INTERSECT_RAW = OUTDIR / "phaseD_encode_ccre_exact_position_intersections_raw_v1.bed"
RESULTS_TSV = OUTDIR / "phaseD_encode_ccre_exact_position_overlap_results_v1.tsv"
SUMMARY_TSV = OUTDIR / "phaseD_encode_ccre_exact_position_overlap_summary_v1.tsv"
FILE_MANIFEST = OUTDIR / "phaseD_encode_ccre_exact_position_overlap_file_manifest_v1.tsv"

CANDIDATES = [
    {"candidate": "1", "gene": "TRAF3", "rsid": "rs1131877", "chromosome": "chr14", "position": 102875712, "ref": "T", "alt": "C", "caution": "none"},
    {"candidate": "2", "gene": "BIRC3", "rsid": "rs924735278", "chromosome": "chr11", "position": 102324736, "ref": "G", "alt": "T", "caution": "identity_caution_active"},
    {"candidate": "3", "gene": "TNFAIP3", "rsid": "rs376205580", "chromosome": "chr6", "position": 137874871, "ref": "A", "alt": "G", "caution": "preserve_paired_haplotype_context"},
    {"candidate": "4", "gene": "TNFAIP3", "rsid": "rs141807543", "chromosome": "chr6", "position": 137875820, "ref": "A", "alt": "C", "caution": "preserve_paired_haplotype_context"},
]

def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()

def run(cmd: list[str], log_name: str) -> subprocess.CompletedProcess[str]:
    log_path = OUTDIR / log_name
    proc = subprocess.run(cmd, text=True, capture_output=True)
    log_path.write_text(
        "COMMAND\t" + " ".join(cmd) + "\n"
        + "RETURN_CODE\t" + str(proc.returncode) + "\n"
        + "STDOUT\n" + proc.stdout + "\n"
        + "STDERR\n" + proc.stderr + "\n",
        encoding="utf-8",
    )
    if proc.returncode != 0:
        raise RuntimeError(f"Command failed: {' '.join(cmd)}; see {log_path}")
    return proc

bigbedtobed = shutil.which("bigBedToBed")
bedtools = shutil.which("bedtools")

if not bigbedtobed:
    raise RuntimeError("bigBedToBed not available.")
if not bedtools:
    raise RuntimeError("bedtools not available.")

# Download selected BigBed file if not already present.
if not BIGBED_PATH.exists() or BIGBED_PATH.stat().st_size == 0:
    urllib.request.urlretrieve(BIGBED_URL, BIGBED_PATH)

# Convert BigBed to BED.
run([bigbedtobed, str(BIGBED_PATH), str(CCRE_BED_RAW)], "phaseD_encode_ccre_bigbedtobed_log_v1.txt")

# Candidate BED: 0-based half-open interval [pos-1, pos)
with CANDIDATE_BED.open("w", encoding="utf-8", newline="") as handle:
    for c in CANDIDATES:
        start = c["position"] - 1
        end = c["position"]
        name = "|".join([c["candidate"], c["gene"], c["rsid"], c["ref"], c["alt"], c["caution"]])
        handle.write(f"{c['chromosome']}\t{start}\t{end}\t{name}\n")

# Exact overlap.
run(
    [bedtools, "intersect", "-a", str(CANDIDATE_BED), "-b", str(CCRE_BED_RAW), "-wa", "-wb"],
    "phaseD_encode_ccre_bedtools_intersect_log_v1.txt",
).stdout

# Capture intersect output directly.
proc = subprocess.run(
    [bedtools, "intersect", "-a", str(CANDIDATE_BED), "-b", str(CCRE_BED_RAW), "-wa", "-wb"],
    text=True,
    capture_output=True,
    check=True,
)
INTERSECT_RAW.write_text(proc.stdout, encoding="utf-8")

# Parse candidate metadata.
candidate_map = {}
for c in CANDIDATES:
    name = "|".join([c["candidate"], c["gene"], c["rsid"], c["ref"], c["alt"], c["caution"]])
    candidate_map[name] = c

overlap_counts = {c["rsid"]: 0 for c in CANDIDATES}
result_rows = []

if proc.stdout.strip():
    for line in proc.stdout.splitlines():
        fields = line.split("\t")
        a_chrom, a_start, a_end, a_name = fields[:4]
        b_fields = fields[4:]
        cand = candidate_map.get(a_name, {})
        rsid = cand.get("rsid", "")
        overlap_counts[rsid] = overlap_counts.get(rsid, 0) + 1
        result_rows.append({
            "candidate": cand.get("candidate", ""),
            "gene": cand.get("gene", ""),
            "rsid": rsid,
            "confirmed_build": "GRCh38",
            "candidate_chromosome": a_chrom,
            "candidate_bed_start_0based": a_start,
            "candidate_bed_end_0based": a_end,
            "variant_position_1based": str(cand.get("position", "")),
            "ref": cand.get("ref", ""),
            "alt": cand.get("alt", ""),
            "candidate_caution": cand.get("caution", ""),
            "ccre_overlap_found": "YES",
            "ccre_overlap_index_for_candidate": str(overlap_counts[rsid]),
            "ccre_raw_bed_fields": "|".join(b_fields),
            "candidate_overlap_performed": "YES",
            "candidate_regulatory_evidence_retrieved": "ANNOTATION_OVERLAP_METADATA_ONLY",
            "biological_interpretation_performed": "NO",
            "functional_validation_claim_made": "NO",
            "causality_claim_made": "NO",
            "selection_or_population_claim_made": "NO",
            "ancestry_effect_claim_made": "NO",
        })

# Add explicit no-overlap rows.
for c in CANDIDATES:
    if overlap_counts.get(c["rsid"], 0) == 0:
        result_rows.append({
            "candidate": c["candidate"],
            "gene": c["gene"],
            "rsid": c["rsid"],
            "confirmed_build": "GRCh38",
            "candidate_chromosome": c["chromosome"],
            "candidate_bed_start_0based": str(c["position"] - 1),
            "candidate_bed_end_0based": str(c["position"]),
            "variant_position_1based": str(c["position"]),
            "ref": c["ref"],
            "alt": c["alt"],
            "candidate_caution": c["caution"],
            "ccre_overlap_found": "NO",
            "ccre_overlap_index_for_candidate": "0",
            "ccre_raw_bed_fields": "",
            "candidate_overlap_performed": "YES",
            "candidate_regulatory_evidence_retrieved": "NO_OVERLAP_METADATA_ONLY",
            "biological_interpretation_performed": "NO",
            "functional_validation_claim_made": "NO",
            "causality_claim_made": "NO",
            "selection_or_population_claim_made": "NO",
            "ancestry_effect_claim_made": "NO",
        })

with RESULTS_TSV.open("w", encoding="utf-8", newline="") as handle:
    fieldnames = list(result_rows[0].keys())
    writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(result_rows)

candidates_with_overlap = sum(1 for c in CANDIDATES if overlap_counts.get(c["rsid"], 0) > 0)
total_overlaps = sum(overlap_counts.values())

summary_rows = [
    {"metric": "selected_file_url", "value": BIGBED_URL},
    {"metric": "bigBedToBed_available", "value": "YES"},
    {"metric": "bedtools_available", "value": "YES"},
    {"metric": "confirmed_build", "value": "GRCh38"},
    {"metric": "candidate_count", "value": len(CANDIDATES)},
    {"metric": "candidate_overlap_performed", "value": "YES"},
    {"metric": "candidate_regulatory_evidence_retrieved", "value": "ANNOTATION_OVERLAP_METADATA_ONLY"},
    {"metric": "candidates_with_ccre_overlap", "value": candidates_with_overlap},
    {"metric": "total_ccre_overlap_rows", "value": total_overlaps},
    {"metric": "interval_rule", "value": "exact_variant_position_1bp"},
    {"metric": "bed_interval_convention", "value": "0_based_half_open"},
    {"metric": "expanded_window_used", "value": "NO"},
    {"metric": "biological_interpretation_performed", "value": "NO"},
    {"metric": "functional_validation_claim_authorised", "value": "NO"},
    {"metric": "causality_claim_authorised", "value": "NO"},
    {"metric": "selection_or_population_claim_authorised", "value": "NO"},
    {"metric": "broad_biological_interpretation_authorised", "value": "NO"},
]
with SUMMARY_TSV.open("w", encoding="utf-8", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=["metric", "value"], delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(summary_rows)

files = [
    (BIGBED_PATH, "source_bigbed"),
    (CCRE_BED_RAW, "converted_bed"),
    (CANDIDATE_BED, "candidate_bed"),
    (INTERSECT_RAW, "raw_intersection"),
    (RESULTS_TSV, "results"),
    (SUMMARY_TSV, "summary"),
    (OUTDIR / "phaseD_encode_ccre_bigbedtobed_log_v1.txt", "bigbedtobed_log"),
    (OUTDIR / "phaseD_encode_ccre_bedtools_intersect_log_v1.txt", "bedtools_log"),
]
with FILE_MANIFEST.open("w", encoding="utf-8", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=["file", "role", "exists", "size_bytes", "sha256"], delimiter="\t", lineterminator="\n")
    writer.writeheader()
    for path, role in files:
        writer.writerow({
            "file": str(path),
            "role": role,
            "exists": "YES" if path.exists() else "NO",
            "size_bytes": path.stat().st_size if path.exists() else 0,
            "sha256": sha256(path) if path.exists() else "",
        })

payload = {
    "status": "PASS",
    "operation": "Execute ENCODE cCRE exact-position overlap using confirmed GRCh38 candidates.",
    "selected_file_url": BIGBED_URL,
    "bigBedToBed_available": "YES",
    "bedtools_available": "YES",
    "confirmed_build": "GRCh38",
    "candidate_count": len(CANDIDATES),
    "candidate_overlap_performed": "YES",
    "candidate_regulatory_evidence_retrieved": "ANNOTATION_OVERLAP_METADATA_ONLY",
    "candidates_with_ccre_overlap": candidates_with_overlap,
    "total_ccre_overlap_rows": total_overlaps,
    "interval_rule": "exact_variant_position_1bp",
    "bed_interval_convention": "0_based_half_open",
    "expanded_window_used": "NO",
    "biological_interpretation_performed": "NO",
    "functional_validation_claim_authorised": "NO",
    "causality_claim_authorised": "NO",
    "selection_or_population_claim_authorised": "NO",
    "broad_biological_interpretation_authorised": "NO",
    "outputs": {
        "results": str(RESULTS_TSV),
        "summary": str(SUMMARY_TSV),
        "file_manifest": str(FILE_MANIFEST),
    },
    "boundary": (
        "This checkpoint records ENCODE cCRE exact-position overlap metadata only. "
        "It does not interpret regulatory biology, validate function, infer causality, "
        "infer selection, assign ancestry effects or authorise broad biological interpretation."
    ),
}

OUTJSON.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
print(json.dumps(payload, indent=2))
