#!/usr/bin/env python3
"""Validate Phase B extraction/harmonisation QC completion without authorising biological analysis."""

from __future__ import annotations

import csv
import json
from pathlib import Path

QC = Path("results/qc")
OUTDIR = Path("11_results/phaseB/gate_closure")
OUTJSON = QC / "phaseB_extraction_harmonisation_gate_closure_validation_v1.json"

AUTOSOMAL_STRICT = Path(
    "11_results/phaseB/genotype_call_concordance_qc/autosomal_filter_rule/"
    "autosomal_strict_interoperability_retained_allele_keys_v1.tsv"
)
AUTOSOMAL_REVIEW = Path(
    "11_results/phaseB/genotype_call_concordance_qc/autosomal_filter_rule/"
    "autosomal_interoperability_technical_review_allele_keys_v1.tsv"
)
CHRX_STRICT = Path(
    "11_results/phaseB/genotype_call_concordance_qc/chrX_IKBKG_filter_rule/"
    "IKBKG_chrX_strict_interoperability_retained_allele_keys_v1.tsv"
)
CHRX_REVIEW = Path(
    "11_results/phaseB/genotype_call_concordance_qc/chrX_IKBKG_filter_rule/"
    "IKBKG_chrX_interoperability_technical_review_allele_keys_v1.tsv"
)

JSON_INPUTS = {
    "full_panel_extraction_qc": QC / "phaseB_full_panel_interval_extraction_validation_v1.json",
    "cross_resource_allele_key_harmonisation": QC / "phaseB_cross_resource_allele_key_harmonisation_validation_v1.json",
    "autosomal_GT_concordance": QC / "phaseB_autosomal_full_panel_direct_GT_concordance_v1.json",
    "autosomal_interoperability_rule": QC / "phaseB_autosomal_interoperability_technical_rule_v1.json",
    "IKBKG_sex_aware_GT_concordance": QC / "phaseB_IKBKG_chrX_sex_aware_direct_GT_concordance_v1.json",
    "IKBKG_interoperability_rule": QC / "phaseB_IKBKG_chrX_interoperability_technical_rule_v1.json",
}


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def load_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


missing_files = [str(p) for p in [*JSON_INPUTS.values(), AUTOSOMAL_STRICT, AUTOSOMAL_REVIEW, CHRX_STRICT, CHRX_REVIEW] if not p.exists()]
failures: list[str] = []
checks: dict[str, bool] = {}

if missing_files:
    failures.extend(f"Missing required gate input: {p}" for p in missing_files)
    payload = {"status": "FAIL", "missing_files": missing_files, "failures": failures}
    OUTJSON.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(payload, indent=2))
    raise SystemExit(1)

j = {name: load_json(path) for name, path in JSON_INPUTS.items()}

checks["all_prerequisite_JSON_status_PASS"] = all(x.get("status") == "PASS" for x in j.values())
checks["full_panel_extracts_validated_48"] = j["full_panel_extraction_qc"]["validated_resource_extracts"] == 48
checks["all_phased_record_keys_genotype_linked"] = (
    j["cross_resource_allele_key_harmonisation"]["phased_only_record_keys"] == 0
    and j["cross_resource_allele_key_harmonisation"]["phased_unique_record_keys"] == 23186
)
checks["autosomal_loci_validated_23"] = j["autosomal_GT_concordance"]["autosomal_loci_validated"] == 23
checks["autosomal_strict_rule_counts"] = (
    j["autosomal_interoperability_rule"]["strict_interoperability_retained_keys"] == 20553
    and j["autosomal_interoperability_rule"]["technical_review_set_keys"] == 2293
)
checks["IKBKG_sex_aware_validation"] = (
    j["IKBKG_sex_aware_GT_concordance"]["XX_sample_count"] == 1896
    and j["IKBKG_sex_aware_GT_concordance"]["XY_sample_count"] == 2195
    and j["IKBKG_sex_aware_GT_concordance"]["phased_only_allele_keys"] == 0
    and j["IKBKG_sex_aware_GT_concordance"]["unexpected_heterozygous_or_other_XY_phased_encodings"] == 0
)
checks["IKBKG_strict_rule_counts"] = (
    j["IKBKG_interoperability_rule"]["strict_interoperability_retained_keys"] == 240
    and j["IKBKG_interoperability_rule"]["technical_review_set_keys"] == 100
)

autosomal_strict = load_tsv(AUTOSOMAL_STRICT)
autosomal_review = load_tsv(AUTOSOMAL_REVIEW)
chrx_strict = load_tsv(CHRX_STRICT)
chrx_review = load_tsv(CHRX_REVIEW)

autosomal_strict_keys = {row["allele_key"] for row in autosomal_strict}
autosomal_review_keys = {row["allele_key"] for row in autosomal_review}
chrx_strict_keys = {row["allele_key"] for row in chrx_strict}
chrx_review_keys = {row["allele_key"] for row in chrx_review}

checks["strict_inventory_lengths_match_rule_JSON"] = (
    len(autosomal_strict) == 20553 and len(chrx_strict) == 240
)
checks["review_inventory_lengths_match_rule_JSON"] = (
    len(autosomal_review) == 2293 and len(chrx_review) == 100
)
checks["strict_and_review_sets_nonoverlapping"] = (
    not (autosomal_strict_keys & autosomal_review_keys)
    and not (chrx_strict_keys & chrx_review_keys)
)
checks["combined_strict_key_count_20793"] = len(autosomal_strict) + len(chrx_strict) == 20793
checks["combined_review_key_count_2393"] = len(autosomal_review) + len(chrx_review) == 2393

for name, passed in checks.items():
    if not passed:
        failures.append(f"Failed gate check: {name}")

OUTDIR.mkdir(parents=True, exist_ok=True)

with (OUTDIR / "phaseB_strict_phased_interoperability_technical_inventory_v1.tsv").open("w", newline="", encoding="utf-8") as handle:
    writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
    writer.writerow(["technical_layer", "locus_scope", "allele_key", "rule_basis"])
    for row in autosomal_strict:
        writer.writerow(["autosomal", row["gene_symbol"], row["allele_key"], "zero_GT_discordance_and_missingness_le_0.01"])
    for row in chrx_strict:
        writer.writerow(["chrX_sex_aware", "IKBKG", row["allele_key"], "zero_GT_discordance_and_XX_XY_unavailable_le_0.01"])

with (OUTDIR / "phaseB_technical_review_sensitivity_inventory_v1.tsv").open("w", newline="", encoding="utf-8") as handle:
    writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
    writer.writerow(["technical_layer", "locus_scope", "allele_key", "technical_rule_reason"])
    for row in autosomal_review:
        writer.writerow(["autosomal", row["gene_symbol"], row["allele_key"], row["technical_rule_reason"]])
    for row in chrx_review:
        writer.writerow(["chrX_sex_aware", "IKBKG", row["allele_key"], row["technical_rule_reason"]])

with (OUTDIR / "phaseB_extraction_harmonisation_gate_closure_summary_v1.tsv").open("w", newline="", encoding="utf-8") as handle:
    writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
    writer.writerow(["gate_component", "result"])
    writer.writerow(["frozen_gene_loci_extracted", 24])
    writer.writerow(["validated_resource_extracts", 48])
    writer.writerow(["shared_phased_genotype_record_keys", 23186])
    writer.writerow(["phased_only_record_keys", 0])
    writer.writerow(["autosomal_strict_interoperability_keys", 20553])
    writer.writerow(["autosomal_technical_review_keys", 2293])
    writer.writerow(["IKBKG_chrX_strict_interoperability_keys", 240])
    writer.writerow(["IKBKG_chrX_technical_review_keys", 100])
    writer.writerow(["combined_strict_interoperability_keys", 20793])
    writer.writerow(["combined_technical_review_keys", 2393])

payload = {
    "status": "PASS" if not failures else "FAIL",
    "operation": "Phase B extraction and genotype-phased harmonisation technical gate-closure validation.",
    "checks": checks,
    "summary": {
        "frozen_gene_loci": 24,
        "validated_resource_extracts": 48,
        "shared_genotype_phased_allele_keys": 23186,
        "phased_only_allele_keys": 0,
        "autosomal_strict_interoperability_keys": 20553,
        "IKBKG_chrX_strict_interoperability_keys": 240,
        "combined_strict_interoperability_keys": 20793,
        "combined_technical_review_sensitivity_keys": 2393,
    },
    "remaining_boundaries": [
        "Genotype-only variant QC remains separate from phased/interoperability rules.",
        "Strict technical inventories do not themselves authorise biological interpretation.",
        "Allele-frequency, haplotype, selection and functional analyses require a subsequent authorised analysis phase."
    ],
    "failures": failures,
}
OUTJSON.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
print(json.dumps(payload, indent=2))

if payload["status"] != "PASS":
    raise SystemExit(1)
