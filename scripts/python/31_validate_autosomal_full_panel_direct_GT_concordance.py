#!/usr/bin/env python3
"""Direct GT-only concordance QC across shared autosomal genotype/phased allele keys."""

from __future__ import annotations

import csv
import json
import subprocess
from pathlib import Path

MANIFEST = Path("config/phaseB/full_panel_extraction/nfkb_full_panel_gene_level_extraction_manifest_GRCh38.tsv")
COHORT = Path("11_results/phaseB/harmonisation_qc/cohort_direct_genotype_phased_interoperability_n4091_v1.txt")
OUTDIR = Path("11_results/phaseB/genotype_call_concordance_qc/autosomal_full_panel")
OUTJSON = Path("results/qc/phaseB_autosomal_full_panel_direct_GT_concordance_v1.json")

PER_LOCUS = OUTDIR / "autosomal_direct_GT_concordance_by_locus_v1.tsv"
PER_SITE = OUTDIR / "autosomal_direct_GT_concordance_by_shared_allele_key_v1.tsv"
DISCORDANT = OUTDIR / "autosomal_shared_allele_keys_with_discordant_GT_calls_v1.tsv"


def run_text(args: list[str]) -> str:
    result = subprocess.run(args, check=True, text=True, capture_output=True)
    return result.stdout


def sample_order(path: str) -> list[str]:
    text = run_text(["bcftools", "query", "-l", "-S", str(COHORT), path])
    return [line for line in text.splitlines() if line]


def gt_record_stream(path: str):
    cmd = [
        "bcftools", "query", "-S", str(COHORT),
        "-f", "%CHROM:%POS:%REF:%ALT[\t%GT]\n", path
    ]
    proc = subprocess.Popen(cmd, text=True, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
    assert proc.stdout is not None
    for line in proc.stdout:
        fields = line.rstrip("\n").split("\t")
        yield fields[0], fields[1:]
    return_code = proc.wait()
    if return_code != 0:
        raise RuntimeError(f"bcftools query failed for {path}, exit code {return_code}")


def canonical_diploid_gt(gt: str) -> str | None:
    if "." in gt:
        return None
    alleles = gt.replace("|", "/").split("/")
    if len(alleles) != 2:
        return None
    return "/".join(sorted(alleles))


def compare_vectors(genotype_gts: list[str], phased_gts: list[str]) -> tuple[int, int, int, int]:
    if len(genotype_gts) != len(phased_gts):
        raise RuntimeError("GT-vector length differs between genotype and phased resources.")

    compared = matches = discordant = unavailable = 0
    for genotype_gt, phased_gt in zip(genotype_gts, phased_gts):
        g = canonical_diploid_gt(genotype_gt)
        p = canonical_diploid_gt(phased_gt)
        if g is None or p is None:
            unavailable += 1
        else:
            compared += 1
            if g == p:
                matches += 1
            else:
                discordant += 1
    return compared, matches, discordant, unavailable


OUTDIR.mkdir(parents=True, exist_ok=True)
OUTJSON.parent.mkdir(parents=True, exist_ok=True)

with MANIFEST.open(newline="", encoding="utf-8") as handle:
    manifest_rows = list(csv.DictReader(handle, delimiter="\t"))

autosomal_rows = [row for row in manifest_rows if row["chromosome"] != "chrX"]

expected_samples = sum(1 for _ in COHORT.open(encoding="utf-8"))
failures: list[str] = []
locus_results: list[dict[str, object]] = []
site_results: list[dict[str, object]] = []

total_shared_keys = 0
total_genotype_keys = 0
total_phased_keys = 0
total_phased_only_keys = 0
total_potential_cells = 0
total_compared = 0
total_matches = 0
total_discordant = 0
total_unavailable = 0

for row in autosomal_rows:
    gene = row["gene_symbol"]
    genotype_path = row["genotype_output"]
    phased_path = row["phased_output"]

    print(f"[{gene}] validating sample order and direct GT concordance", flush=True)

    genotype_samples = sample_order(genotype_path)
    phased_samples = sample_order(phased_path)

    order_match = genotype_samples == phased_samples
    if not order_match:
        failures.append(f"{gene}: genotype/phased sample order mismatch after cohort restriction")
    if len(genotype_samples) != expected_samples or len(phased_samples) != expected_samples:
        failures.append(f"{gene}: restricted sample count differs from expected {expected_samples}")

    phased_records = {key: gts for key, gts in gt_record_stream(phased_path)}
    genotype_key_count = 0
    shared_count = 0
    locus_compared = locus_matches = locus_discordant = locus_unavailable = 0
    seen_shared: set[str] = set()

    for key, genotype_gts in gt_record_stream(genotype_path):
        genotype_key_count += 1
        phased_gts = phased_records.get(key)
        if phased_gts is None:
            continue

        shared_count += 1
        seen_shared.add(key)
        compared, matches, discordant, unavailable = compare_vectors(genotype_gts, phased_gts)

        locus_compared += compared
        locus_matches += matches
        locus_discordant += discordant
        locus_unavailable += unavailable

        site_results.append({
            "gene_symbol": gene,
            "allele_key": key,
            "nonmissing_GT_comparisons": compared,
            "matching_GT_calls": matches,
            "discordant_GT_calls": discordant,
            "missing_or_nondiploid_in_either_resource": unavailable,
        })

    phased_only = sorted(set(phased_records) - seen_shared)
    if phased_only:
        failures.append(f"{gene}: {len(phased_only)} phased-only allele keys detected")

    locus_fraction = locus_matches / locus_compared if locus_compared else None
    locus_results.append({
        "gene_symbol": gene,
        "chromosome": row["chromosome"],
        "genotype_allele_keys": genotype_key_count,
        "phased_allele_keys": len(phased_records),
        "shared_allele_keys_compared": shared_count,
        "phased_only_allele_keys": len(phased_only),
        "potential_sample_by_site_cells": shared_count * expected_samples,
        "nonmissing_GT_comparisons": locus_compared,
        "matching_GT_calls": locus_matches,
        "discordant_GT_calls": locus_discordant,
        "missing_or_nondiploid_in_either_resource_cells": locus_unavailable,
        "technical_concordance_fraction": locus_fraction,
        "sample_order_match": "PASS" if order_match else "FAIL",
    })

    total_genotype_keys += genotype_key_count
    total_phased_keys += len(phased_records)
    total_shared_keys += shared_count
    total_phased_only_keys += len(phased_only)
    total_potential_cells += shared_count * expected_samples
    total_compared += locus_compared
    total_matches += locus_matches
    total_discordant += locus_discordant
    total_unavailable += locus_unavailable

with PER_LOCUS.open("w", newline="", encoding="utf-8") as handle:
    fields = list(locus_results[0].keys())
    writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(locus_results)

with PER_SITE.open("w", newline="", encoding="utf-8") as handle:
    fields = list(site_results[0].keys())
    writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(site_results)

with DISCORDANT.open("w", newline="", encoding="utf-8") as handle:
    fields = list(site_results[0].keys())
    writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows([row for row in site_results if int(row["discordant_GT_calls"]) > 0])

payload = {
    "status": "PASS" if not failures and len(locus_results) == 23 and total_phased_only_keys == 0 else "FAIL",
    "method": "Direct GT-only comparison after exact CHROM:POS:REF:ALT matching in frozen 4091-sample interoperability cohort; autosomal loci only.",
    "autosomal_loci_validated": len(locus_results),
    "excluded_pending_sex_aware_qc": ["IKBKG_chrX"],
    "interoperability_sample_count": expected_samples,
    "genotype_allele_keys_in_autosomal_loci": total_genotype_keys,
    "phased_allele_keys_in_autosomal_loci": total_phased_keys,
    "shared_allele_keys_compared": total_shared_keys,
    "phased_only_allele_keys": total_phased_only_keys,
    "potential_sample_by_site_cells": total_potential_cells,
    "nonmissing_GT_comparisons": total_compared,
    "matching_GT_calls": total_matches,
    "discordant_GT_calls": total_discordant,
    "missing_or_nondiploid_in_either_resource_cells": total_unavailable,
    "technical_concordance_fraction": total_matches / total_compared if total_compared else None,
    "allele_keys_with_any_discordant_GT_calls": sum(
        int(row["discordant_GT_calls"]) > 0 for row in site_results
    ),
    "failures": failures,
    "boundary": "Phase separators are ignored for GT-call concordance; haplotype phase-direction concordance and chrX sex-aware concordance are not evaluated here."
}

OUTJSON.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
print(json.dumps(payload, indent=2))

if payload["status"] != "PASS":
    raise SystemExit(1)
