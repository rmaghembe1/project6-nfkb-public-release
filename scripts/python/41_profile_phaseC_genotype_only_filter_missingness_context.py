#!/usr/bin/env python3
from __future__ import annotations

import csv
import json
import subprocess
from collections import defaultdict, Counter
from pathlib import Path

MANIFEST = Path("config/phaseB/full_panel_extraction/nfkb_full_panel_gene_level_extraction_manifest_GRCh38.tsv")
AUTO = Path("11_results/phaseC/genotype_only_qc/autosomal_call_availability/phaseC_autosomal_genotype_only_call_availability_by_allele_key_v1.tsv")
CHRX = Path("11_results/phaseC/genotype_only_qc/chrX_IKBKG_call_availability/phaseC_IKBKG_genotype_only_call_availability_by_allele_key_v1.tsv")

OUT = Path("11_results/phaseC/genotype_only_qc/site_quality_context")
OUTJSON = Path("results/qc/phaseC_genotype_only_filter_missingness_context_profile_v1.json")

OUT.mkdir(parents=True, exist_ok=True)
OUTJSON.parent.mkdir(parents=True, exist_ok=True)


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def filter_map(vcf: str) -> dict[str, str]:
    cmd = ["bcftools", "query", "-f", "%CHROM:%POS:%REF:%ALT\t%FILTER\n", vcf]
    run = subprocess.run(cmd, check=True, text=True, capture_output=True)
    out = {}
    for line in run.stdout.splitlines():
        key, filt = line.split("\t", 1)
        out[key] = filt
    return out


manifest = read_tsv(MANIFEST)
auto = read_tsv(AUTO)
chrx = read_tsv(CHRX)

vcf_by_gene = {row["gene_symbol"]: row["genotype_output"] for row in manifest}
filters = {gene: filter_map(vcf) for gene, vcf in vcf_by_gene.items()}

failures = []
rows = []

for row in auto:
    gene = row["gene_symbol"]
    key = row["allele_key"]
    filt = filters[gene].get(key)
    if filt is None:
        failures.append(f"Missing FILTER record for {gene} {key}")
        continue
    rows.append({
        "analysis_layer": "autosomal",
        "gene_symbol": gene,
        "allele_key": key,
        "variant_class": row["variant_class"],
        "site_FILTER": filt,
        "missingness_fraction": float(row["missingness_fraction"]),
    })

for row in chrx:
    key = row["allele_key"]
    filt = filters["IKBKG"].get(key)
    if filt is None:
        failures.append(f"Missing FILTER record for IKBKG {key}")
        continue
    rows.append({
        "analysis_layer": "chrX_sex_aware",
        "gene_symbol": "IKBKG",
        "allele_key": key,
        "variant_class": row["variant_class"],
        "site_FILTER": filt,
        "missingness_fraction": float(row["combined_missingness_fraction"]),
    })

auto_n = sum(r["analysis_layer"] == "autosomal" for r in rows)
chrx_n = sum(r["analysis_layer"] == "chrX_sex_aware" for r in rows)

if auto_n != 55592:
    failures.append(f"Expected 55,592 autosomal records; observed {auto_n}")
if chrx_n != 937:
    failures.append(f"Expected 937 chrX records; observed {chrx_n}")

joined = OUT / "phaseC_genotype_only_filter_missingness_joined_inventory_v1.tsv"
with joined.open("w", newline="", encoding="utf-8") as handle:
    fields = ["analysis_layer", "gene_symbol", "allele_key", "variant_class", "site_FILTER", "missingness_fraction"]
    writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)

by_filter = defaultdict(list)
by_filter_class = defaultdict(list)
by_locus_filter = defaultdict(list)

for row in rows:
    miss = float(row["missingness_fraction"])
    by_filter[(row["analysis_layer"], row["site_FILTER"])].append(miss)
    by_filter_class[(row["analysis_layer"], row["site_FILTER"], row["variant_class"])].append(miss)
    by_locus_filter[(row["gene_symbol"], row["site_FILTER"])].append(miss)


def median(vals: list[float]) -> float:
    vals = sorted(vals)
    n = len(vals)
    return vals[n // 2] if n % 2 else (vals[n // 2 - 1] + vals[n // 2]) / 2


with (OUT / "phaseC_genotype_only_missingness_by_FILTER_status_v1.tsv").open("w", newline="", encoding="utf-8") as handle:
    writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
    writer.writerow(["analysis_layer", "site_FILTER", "record_count", "records_missingness_gt_0.01", "records_missingness_gt_0.05", "mean_missingness_fraction", "median_missingness_fraction"])
    for key, vals in sorted(by_filter.items()):
        writer.writerow([key[0], key[1], len(vals), sum(v > 0.01 for v in vals), sum(v > 0.05 for v in vals), sum(vals) / len(vals), median(vals)])

with (OUT / "phaseC_genotype_only_missingness_by_FILTER_and_variant_class_v1.tsv").open("w", newline="", encoding="utf-8") as handle:
    writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
    writer.writerow(["analysis_layer", "site_FILTER", "variant_class", "record_count", "records_missingness_gt_0.01", "records_missingness_gt_0.05", "mean_missingness_fraction"])
    for key, vals in sorted(by_filter_class.items()):
        writer.writerow([key[0], key[1], key[2], len(vals), sum(v > 0.01 for v in vals), sum(v > 0.05 for v in vals), sum(vals) / len(vals)])

with (OUT / "phaseC_genotype_only_missingness_by_locus_and_FILTER_status_v1.tsv").open("w", newline="", encoding="utf-8") as handle:
    writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
    writer.writerow(["gene_symbol", "site_FILTER", "record_count", "records_missingness_gt_0.01", "records_missingness_gt_0.05", "mean_missingness_fraction"])
    for key, vals in sorted(by_locus_filter.items()):
        writer.writerow([key[0], key[1], len(vals), sum(v > 0.01 for v in vals), sum(v > 0.05 for v in vals), sum(vals) / len(vals)])

filter_counts = Counter((r["analysis_layer"], r["site_FILTER"]) for r in rows)

payload = {
    "status": "PASS" if not failures else "FAIL",
    "operation": "Descriptive genotype-only FILTER and missingness-context audit.",
    "records_joined": {"autosomal": auto_n, "chrX_sex_aware": chrx_n},
    "site_FILTER_counts": {f"{k[0]}|{k[1]}": v for k, v in sorted(filter_counts.items())},
    "failures": failures,
    "boundary": "No FILTER-based exclusion, missingness threshold, frequency calculation or biological interpretation is authorised."
}

OUTJSON.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
print(json.dumps(payload, indent=2))

if payload["status"] != "PASS":
    raise SystemExit(1)
