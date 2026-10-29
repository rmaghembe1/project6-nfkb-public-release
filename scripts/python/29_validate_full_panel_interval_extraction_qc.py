#!/usr/bin/env python3
"""Validate completed Phase B full-panel interval extraction without biological interpretation."""

import csv
import json
import subprocess
from pathlib import Path

MANIFEST = Path("config/phaseB/full_panel_extraction/nfkb_full_panel_gene_level_extraction_manifest_GRCh38.tsv")
OUT_TSV = Path("11_results/phaseB/interval_extraction/full_panel/nfkb_full_panel_extraction_qc_summary_v1.tsv")
OUT_JSON = Path("results/qc/phaseB_full_panel_interval_extraction_validation_v1.json")

def run(args):
    return subprocess.check_output(args, text=True, stderr=subprocess.DEVNULL).strip()

def sample_count(path):
    text = run(["bcftools", "query", "-l", path])
    return 0 if not text else len(text.splitlines())

def positions(path):
    text = run(["bcftools", "query", "-f", "%POS\n", path])
    vals = [int(x) for x in text.splitlines() if x]
    return vals

with MANIFEST.open(newline="", encoding="utf-8") as h:
    rows = list(csv.DictReader(h, delimiter="\t"))

qc_rows = []
failures = []

for row in rows:
    gene = row["gene_symbol"]
    chrom = row["chromosome"]
    start, end = map(int, row["region_1based_inclusive"].split(":")[1].split("-"))

    for resource, field, suffix, expected_samples in [
        ("genotype", "genotype_output", ".tbi", 4151),
        ("phased", "phased_output", ".csi", 4091),
    ]:
        output = Path(row[field])
        index = Path(str(output) + suffix)

        if not output.exists() or not index.exists():
            failures.append(f"{gene}:{resource}: missing output or index")
            continue

        samples = sample_count(str(output))
        pos = positions(str(output))
        variants = len(pos)
        first = min(pos) if pos else None
        last = max(pos) if pos else None
        contained = bool(pos) and first >= start and last <= end

        if samples != expected_samples:
            failures.append(f"{gene}:{resource}: sample count {samples} != {expected_samples}")
        if variants == 0:
            failures.append(f"{gene}:{resource}: zero extracted variants")
        if not contained:
            failures.append(f"{gene}:{resource}: extracted positions outside requested interval or empty")

        qc_rows.append({
            "gene_symbol": gene,
            "chromosome": chrom,
            "region": row["region_1based_inclusive"],
            "resource": resource,
            "sample_count": samples,
            "variant_count": variants,
            "first_position": first if first is not None else "",
            "last_position": last if last is not None else "",
            "interval_containment": "PASS" if contained else "FAIL",
        })

def anchor_present(path, position, ref, alt):
    fmt = "%CHROM\t%POS\t%REF\t%ALT\n"
    text = run(["bcftools", "query", "-r", f"chr6:{position}-{position}", "-f", fmt, path])
    return any(
        fields[1] == str(position) and fields[2] == ref and alt in fields[3].split(",")
        for fields in (line.split("\t") for line in text.splitlines() if line)
    )

tnfaip3 = next(row for row in rows if row["gene_symbol"] == "TNFAIP3")
anchor_checks = {}
for resource, field in [("genotype", "genotype_output"), ("phased", "phased_output")]:
    path = tnfaip3[field]
    anchor_checks[f"{resource}_rs376205580_A_G"] = anchor_present(path, 137874871, "A", "G")
    anchor_checks[f"{resource}_rs141807543_A_C"] = anchor_present(path, 137875820, "A", "C")

if not all(anchor_checks.values()):
    failures.append("TNFAIP3/A20 anchor recovery failed in one or more resource extracts")

OUT_TSV.parent.mkdir(parents=True, exist_ok=True)
with OUT_TSV.open("w", newline="", encoding="utf-8") as h:
    fields = list(qc_rows[0].keys())
    w = csv.DictWriter(h, fieldnames=fields, delimiter="\t", lineterminator="\n")
    w.writeheader()
    w.writerows(qc_rows)

payload = {
    "status": "PASS" if not failures and len(qc_rows) == 48 else "FAIL",
    "operation": "Full-panel interval extraction file/index/sample/count/containment/anchor QC only.",
    "expected_gene_loci": 24,
    "validated_resource_extracts": len(qc_rows),
    "genotype_extracts_validated": sum(r["resource"] == "genotype" for r in qc_rows),
    "phased_extracts_validated": sum(r["resource"] == "phased" for r in qc_rows),
    "expected_sample_counts": {"genotype": 4151, "phased": 4091},
    "anchor_checks": anchor_checks,
    "failures": failures,
    "interpretation_boundary": "This QC authorises no biological, population-frequency, haplotype or selection interpretation.",
}
OUT_JSON.parent.mkdir(parents=True, exist_ok=True)
OUT_JSON.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
print(json.dumps(payload, indent=2))

if payload["status"] != "PASS":
    raise SystemExit(1)
