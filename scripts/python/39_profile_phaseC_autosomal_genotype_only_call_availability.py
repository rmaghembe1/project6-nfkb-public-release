#!/usr/bin/env python3
"""Profile autosomal genotype-only GT availability before defining Phase C filtering rules."""

from __future__ import annotations

import csv
import json
import subprocess
from pathlib import Path

MANIFEST = Path("config/phaseB/full_panel_extraction/nfkb_full_panel_gene_level_extraction_manifest_GRCh38.tsv")
COHORT = Path("11_results/phaseB/harmonisation_qc/cohort_primary_genotype_metadata_linked_n4150_v1.txt")
OUTDIR = Path("11_results/phaseC/genotype_only_qc/autosomal_call_availability")
OUTJSON = Path("results/qc/phaseC_autosomal_genotype_only_call_availability_profile_v1.json")

EXPECTED_SAMPLES = 4150


def variant_class(key: str) -> str:
    _, _, ref, alt = key.split(":", 3)
    if len(ref) == 1 and len(alt) == 1:
        return "SNV"
    if len(ref) != len(alt):
        return "INDEL"
    return "MNV_OR_COMPLEX"


def record_stream(path: str):
    cmd = [
        "bcftools", "query", "-S", str(COHORT),
        "-f", "%CHROM:%POS:%REF:%ALT[\t%GT]\n", path
    ]
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    assert proc.stdout is not None
    for line in proc.stdout:
        fields = line.rstrip("\n").split("\t")
        yield fields[0], fields[1:]
    stderr = proc.stderr.read() if proc.stderr else ""
    code = proc.wait()
    if code != 0:
        raise RuntimeError(f"bcftools query failed for {path}: {stderr}")


with MANIFEST.open(newline="", encoding="utf-8") as handle:
    manifest = list(csv.DictReader(handle, delimiter="\t"))

autosomal = [row for row in manifest if row["chromosome"] != "chrX"]
site_rows = []
locus_rows = []
failures = []

for row in autosomal:
    gene = row["gene_symbol"]
    path = row["genotype_output"]
    locus_records = locus_missing = locus_unexpected_called = locus_called = 0

    for key, gts in record_stream(path):
        if len(gts) != EXPECTED_SAMPLES:
            failures.append(f"{gene} {key}: expected {EXPECTED_SAMPLES} GT cells; observed {len(gts)}.")
            continue

        missing = 0
        called = 0
        unexpected_called = 0

        for gt in gts:
            if "." in gt:
                missing += 1
            elif gt in {"0/0", "0/1", "1/0", "1/1"}:
                called += 1
            else:
                unexpected_called += 1

        site_rows.append({
            "gene_symbol": gene,
            "allele_key": key,
            "variant_class": variant_class(key),
            "expected_GT_cells": EXPECTED_SAMPLES,
            "called_diploid_GT_cells": called,
            "missing_GT_cells": missing,
            "unexpected_nonmissing_GT_cells": unexpected_called,
            "missingness_fraction": missing / EXPECTED_SAMPLES,
        })

        locus_records += 1
        locus_missing += missing
        locus_called += called
        locus_unexpected_called += unexpected_called

    locus_rows.append({
        "gene_symbol": gene,
        "autosomal_variant_records": locus_records,
        "expected_GT_cells": locus_records * EXPECTED_SAMPLES,
        "called_diploid_GT_cells": locus_called,
        "missing_GT_cells": locus_missing,
        "unexpected_nonmissing_GT_cells": locus_unexpected_called,
        "missingness_fraction": locus_missing / (locus_records * EXPECTED_SAMPLES) if locus_records else "",
    })

OUTDIR.mkdir(parents=True, exist_ok=True)

with (OUTDIR / "phaseC_autosomal_genotype_only_call_availability_by_allele_key_v1.tsv").open("w", newline="", encoding="utf-8") as handle:
    writer = csv.DictWriter(handle, fieldnames=list(site_rows[0]), delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(site_rows)

with (OUTDIR / "phaseC_autosomal_genotype_only_call_availability_by_locus_v1.tsv").open("w", newline="", encoding="utf-8") as handle:
    writer = csv.DictWriter(handle, fieldnames=list(locus_rows[0]), delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(locus_rows)

thresholds = [0, 0.001, 0.01, 0.05]
with (OUTDIR / "phaseC_autosomal_genotype_only_candidate_missingness_threshold_counts_v1.tsv").open("w", newline="", encoding="utf-8") as handle:
    writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
    writer.writerow(["missingness_rule", "flagged_allele_keys", "retained_allele_keys_if_excluded"])
    for threshold in thresholds:
        flagged = sum(float(row["missingness_fraction"]) > threshold for row in site_rows)
        writer.writerow([f"missingness_gt_{threshold}", flagged, len(site_rows) - flagged])

total_expected = sum(int(row["expected_GT_cells"]) for row in locus_rows)
total_missing = sum(int(row["missing_GT_cells"]) for row in locus_rows)
total_unexpected = sum(int(row["unexpected_nonmissing_GT_cells"]) for row in locus_rows)

payload = {
    "status": "PASS" if not failures and len(locus_rows) == 23 and total_unexpected == 0 else "FAIL",
    "operation": "Descriptive autosomal genotype-only GT call-availability profiling only.",
    "genotype_only_sample_count": EXPECTED_SAMPLES,
    "autosomal_loci_profiled": len(locus_rows),
    "autosomal_variant_records_profiled": len(site_rows),
    "expected_GT_cells": total_expected,
    "missing_GT_cells": total_missing,
    "unexpected_nonmissing_GT_cells": total_unexpected,
    "overall_missingness_fraction": total_missing / total_expected if total_expected else None,
    "failures": failures,
    "boundary": "No genotype-only filtering threshold, allele-frequency analysis or biological interpretation is authorised.",
}
OUTJSON.parent.mkdir(parents=True, exist_ok=True)
OUTJSON.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
print(json.dumps(payload, indent=2))

if payload["status"] != "PASS":
    raise SystemExit(1)
