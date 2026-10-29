#!/usr/bin/env python3
import csv
import json
import subprocess
from pathlib import Path

GENO = "03_raw_genomic_resources/interval_extraction/full_panel/genotype/TNFAIP3_GRCh38_gnomad_hgdp_tgp.vcf.gz"
PHASE = "03_raw_genomic_resources/interval_extraction/full_panel/phased/TNFAIP3_GRCh38_hgdp1kgp_phased.bcf"
COHORT = "11_results/phaseB/harmonisation_qc/cohort_direct_genotype_phased_interoperability_n4091_v1.txt"
OUTDIR = Path("11_results/phaseB/genotype_call_concordance_qc/pilot_TNFAIP3")
JSON_OUT = Path("results/qc/phaseB_TNFAIP3_direct_GT_concordance_pilot_v1.json")

def run_query(path):
    cmd = ["bcftools", "query", "-S", COHORT, "-f", "%CHROM:%POS:%REF:%ALT[\t%GT]\n", path]
    result = subprocess.run(cmd, check=True, text=True, capture_output=True)
    records = {}
    for line in result.stdout.splitlines():
        fields = line.split("\t")
        records[fields[0]] = fields[1:]
    return records

def canonical_gt(gt):
    if gt in {".", "./.", ".|."} or "." in gt:
        return None
    alleles = gt.replace("|", "/").split("/")
    return "/".join(sorted(alleles))

genotype = run_query(GENO)
phased = run_query(PHASE)
shared = sorted(set(genotype) & set(phased))

per_site = []
total_compared = 0
total_matches = 0
total_mismatches = 0
total_missing_either = 0

for key in shared:
    gts_g = genotype[key]
    gts_p = phased[key]
    compared = matches = mismatches = missing = 0
    for g, p in zip(gts_g, gts_p):
        cg = canonical_gt(g)
        cp = canonical_gt(p)
        if cg is None or cp is None:
            missing += 1
        else:
            compared += 1
            if cg == cp:
                matches += 1
            else:
                mismatches += 1

    per_site.append([key, compared, matches, mismatches, missing])
    total_compared += compared
    total_matches += matches
    total_mismatches += mismatches
    total_missing_either += missing

OUTDIR.mkdir(parents=True, exist_ok=True)
with (OUTDIR / "TNFAIP3_direct_GT_concordance_by_shared_allele_key_v1.tsv").open("w", newline="") as handle:
    writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
    writer.writerow(["allele_key", "nonmissing_sample_comparisons", "matching_GT_calls", "discordant_GT_calls", "missing_in_either_resource"])
    writer.writerows(per_site)

payload = {
    "status": "PASS" if len(shared) == 363 else "FAIL",
    "method": "Direct GT-only comparison after exact CHROM:POS:REF:ALT matching in frozen 4091-sample interoperability cohort.",
    "genotype_allele_keys": len(genotype),
    "phased_allele_keys": len(phased),
    "shared_allele_keys_compared": len(shared),
    "potential_sample_by_site_cells": len(shared) * 4091,
    "nonmissing_GT_comparisons": total_compared,
    "matching_GT_calls": total_matches,
    "discordant_GT_calls": total_mismatches,
    "missing_in_either_resource_cells": total_missing_either,
    "technical_concordance_fraction": total_matches / total_compared if total_compared else None,
    "boundary": "Phase separators are ignored for GT-call concordance; haplotype phase-direction concordance is not evaluated here."
}
JSON_OUT.parent.mkdir(parents=True, exist_ok=True)
JSON_OUT.write_text(json.dumps(payload, indent=2) + "\n")
print(json.dumps(payload, indent=2))
