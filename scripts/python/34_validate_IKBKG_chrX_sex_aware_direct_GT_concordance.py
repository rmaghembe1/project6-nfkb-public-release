#!/usr/bin/env python3
"""Sex-aware direct GT concordance QC for chrX non-PAR IKBKG."""

from __future__ import annotations

import csv
import json
import subprocess
from pathlib import Path

GENO = "03_raw_genomic_resources/interval_extraction/full_panel/genotype/IKBKG_GRCh38_gnomad_hgdp_tgp.vcf.gz"
PHASE = "03_raw_genomic_resources/interval_extraction/full_panel/phased/IKBKG_GRCh38_hgdp1kgp_phased.bcf"
RECON = Path("11_results/phaseB/genotype_call_concordance_qc/chrX_IKBKG_reconnaissance")
OUTDIR = Path("11_results/phaseB/genotype_call_concordance_qc/chrX_IKBKG_sex_aware")
OUTJSON = Path("results/qc/phaseB_IKBKG_chrX_sex_aware_direct_GT_concordance_v1.json")

STRATA = {
    "XX": RECON / "IKBKG_clean_XX_samples_v1.txt",
    "XY": RECON / "IKBKG_clean_XY_samples_v1.txt",
}


def run_text(args: list[str]) -> str:
    return subprocess.run(args, check=True, text=True, capture_output=True).stdout


def sample_order(path: str, sample_file: Path) -> list[str]:
    text = run_text(["bcftools", "query", "-l", "-S", str(sample_file), path])
    return [line for line in text.splitlines() if line]


def records(path: str, sample_file: Path) -> dict[str, list[str]]:
    text = run_text([
        "bcftools", "query", "-S", str(sample_file),
        "-f", "%CHROM:%POS:%REF:%ALT[\t%GT]\n", path
    ])
    out = {}
    for line in text.splitlines():
        fields = line.split("\t")
        out[fields[0]] = fields[1:]
    return out


def canon_xx(gt: str) -> str | None:
    if "." in gt:
        return None
    alleles = gt.replace("|", "/").split("/")
    if len(alleles) != 2:
        return "UNEXPECTED_PLOIDY"
    return "/".join(sorted(alleles))


def canon_xy_genotype(gt: str) -> str | None:
    if "." in gt:
        return None
    if gt in {"0", "1"}:
        return gt
    return "UNEXPECTED_NONMISSING_DIPLOID_OR_OTHER"


def canon_xy_phased(gt: str) -> str | None:
    if "." in gt:
        return None
    normalized = gt.replace("/", "|")
    if normalized == "0|0":
        return "0"
    if normalized == "1|1":
        return "1"
    return "UNEXPECTED_HETEROZYGOUS_OR_OTHER"


OUTDIR.mkdir(parents=True, exist_ok=True)

summary_rows = []
site_rows = []
failures = []
shared_keys_reference = None
total_compared = total_matches = total_discordant = total_unavailable = 0
unexpected_xy_genotype = unexpected_xy_phased = 0

for stratum, sample_file in STRATA.items():
    n_samples = sum(1 for _ in sample_file.open(encoding="utf-8"))
    geno_order = sample_order(GENO, sample_file)
    phase_order = sample_order(PHASE, sample_file)
    order_match = geno_order == phase_order
    if not order_match or len(geno_order) != n_samples or len(phase_order) != n_samples:
        failures.append(f"{stratum}: sample-order or sample-count mismatch")

    genotype = records(GENO, sample_file)
    phased = records(PHASE, sample_file)
    shared = sorted(set(genotype) & set(phased))
    phased_only = sorted(set(phased) - set(genotype))

    if phased_only:
        failures.append(f"{stratum}: {len(phased_only)} phased-only keys")

    if shared_keys_reference is None:
        shared_keys_reference = shared
    elif shared != shared_keys_reference:
        failures.append("Shared allele-key inventory differs between XX and XY strata")

    stratum_compared = stratum_matches = stratum_discordant = stratum_unavailable = 0

    for key in shared:
        compared = matches = discordant = unavailable = 0
        for g, p in zip(genotype[key], phased[key]):
            if stratum == "XX":
                cg = canon_xx(g)
                cp = canon_xx(p)
            else:
                cg = canon_xy_genotype(g)
                cp = canon_xy_phased(p)
                unexpected_xy_genotype += int(cg == "UNEXPECTED_NONMISSING_DIPLOID_OR_OTHER")
                unexpected_xy_phased += int(cp == "UNEXPECTED_HETEROZYGOUS_OR_OTHER")

            if cg is None or cp is None:
                unavailable += 1
            elif cg.startswith("UNEXPECTED") or cp.startswith("UNEXPECTED"):
                unavailable += 1
            else:
                compared += 1
                if cg == cp:
                    matches += 1
                else:
                    discordant += 1

        site_rows.append({
            "sex_karyotype": stratum,
            "allele_key": key,
            "nonmissing_comparisons": compared,
            "matching_GT_calls": matches,
            "discordant_GT_calls": discordant,
            "unavailable_or_unexpected_cells": unavailable,
        })
        stratum_compared += compared
        stratum_matches += matches
        stratum_discordant += discordant
        stratum_unavailable += unavailable

    summary_rows.append({
        "sex_karyotype": stratum,
        "sample_count": n_samples,
        "genotype_keys": len(genotype),
        "phased_keys": len(phased),
        "shared_keys": len(shared),
        "phased_only_keys": len(phased_only),
        "potential_sample_by_site_cells": len(shared) * n_samples,
        "nonmissing_comparisons": stratum_compared,
        "matching_GT_calls": stratum_matches,
        "discordant_GT_calls": stratum_discordant,
        "unavailable_or_unexpected_cells": stratum_unavailable,
        "technical_concordance_fraction": stratum_matches / stratum_compared if stratum_compared else "",
        "sample_order_match": "PASS" if order_match else "FAIL",
    })

    total_compared += stratum_compared
    total_matches += stratum_matches
    total_discordant += stratum_discordant
    total_unavailable += stratum_unavailable

with (OUTDIR / "IKBKG_chrX_sex_aware_direct_GT_concordance_by_karyotype_v1.tsv").open("w", newline="", encoding="utf-8") as handle:
    writer = csv.DictWriter(handle, fieldnames=list(summary_rows[0]), delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(summary_rows)

with (OUTDIR / "IKBKG_chrX_sex_aware_direct_GT_concordance_by_shared_allele_key_v1.tsv").open("w", newline="", encoding="utf-8") as handle:
    writer = csv.DictWriter(handle, fieldnames=list(site_rows[0]), delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(site_rows)

with (OUTDIR / "IKBKG_chrX_shared_allele_keys_with_discordant_GT_calls_v1.tsv").open("w", newline="", encoding="utf-8") as handle:
    writer = csv.DictWriter(handle, fieldnames=list(site_rows[0]), delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows([row for row in site_rows if int(row["discordant_GT_calls"]) > 0])

payload = {
    "status": "PASS" if not failures and len(shared_keys_reference or []) == 340 and unexpected_xy_phased == 0 else "FAIL",
    "method": "Sex-aware direct GT comparison for chrX non-PAR IKBKG after exact CHROM:POS:REF:ALT matching.",
    "XX_sample_count": summary_rows[0]["sample_count"],
    "XY_sample_count": summary_rows[1]["sample_count"],
    "shared_IKBKG_allele_keys": len(shared_keys_reference or []),
    "phased_only_allele_keys": sum(int(row["phased_only_keys"]) for row in summary_rows),
    "potential_sample_by_site_cells": len(shared_keys_reference or []) * sum(int(row["sample_count"]) for row in summary_rows),
    "nonmissing_comparisons": total_compared,
    "matching_GT_calls": total_matches,
    "discordant_GT_calls": total_discordant,
    "unavailable_or_unexpected_cells": total_unavailable,
    "technical_concordance_fraction": total_matches / total_compared if total_compared else None,
    "unexpected_nonmissing_XY_genotype_encodings": unexpected_xy_genotype,
    "unexpected_heterozygous_or_other_XY_phased_encodings": unexpected_xy_phased,
    "allele_strata_rows_with_any_discordance": sum(int(row["discordant_GT_calls"]) > 0 for row in site_rows),
    "failures": failures,
    "boundary": "XX calls are compared as diploid genotypes; XY genotype haploid calls are compared with homozygous pseudo-diploid phased calls. No biological interpretation is authorised.",
}
OUTJSON.parent.mkdir(parents=True, exist_ok=True)
OUTJSON.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
print(json.dumps(payload, indent=2))

if payload["status"] != "PASS":
    raise SystemExit(1)
