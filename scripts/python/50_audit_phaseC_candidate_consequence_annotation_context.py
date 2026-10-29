#!/usr/bin/env python3
"""Audit consequence-annotation context for Phase C candidate variants.

This checkpoint does not run VEP/SnpEff/ANNOVAR and does not assign functional
consequences. It validates candidate representation, variant classes, gene/source
composition and creates tool-ready manifest files for downstream annotation.
"""

from __future__ import annotations

import csv
import json
import shutil
import subprocess
from collections import Counter
from pathlib import Path

CANDIDATE = Path(
    "11_results/phaseC/genotype_only_qc/functional_annotation_readiness/"
    "phaseC_frequency_variation_candidate_union_for_annotation_v1.tsv"
)
BED = Path(
    "11_results/phaseC/genotype_only_qc/functional_annotation_readiness/"
    "phaseC_frequency_variation_candidate_union_for_annotation_v1.bed"
)
VCF_LIKE = Path(
    "11_results/phaseC/genotype_only_qc/functional_annotation_readiness/"
    "phaseC_frequency_variation_candidate_union_for_annotation_v1.vcf.tsv"
)

OUTDIR = Path("11_results/phaseC/genotype_only_qc/consequence_annotation_context")
OUTJSON = Path("results/qc/phaseC_candidate_consequence_annotation_context_audit_v1.json")

EXPECTED_CANDIDATES = 500

OUTDIR.mkdir(parents=True, exist_ok=True)
OUTJSON.parent.mkdir(parents=True, exist_ok=True)

failures: list[str] = []


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def variant_class(ref: str, alt: str) -> str:
    if len(ref) == 1 and len(alt) == 1:
        return "SNV"
    if len(ref) == len(alt):
        return "MNV_or_complex_substitution"
    if len(ref) < len(alt):
        return "insertion_or_delins"
    if len(ref) > len(alt):
        return "deletion_or_delins"
    return "other"


def tool_version(tool: str) -> str:
    path = shutil.which(tool)
    if not path:
        return "NOT_FOUND"
    version_cmds = {
        "bcftools": [tool, "--version"],
        "bedtools": [tool, "--version"],
        "tabix": [tool, "--version"],
        "bgzip": [tool, "--version"],
        "vep": [tool, "--help"],
        "snpEff": [tool, "-version"],
    }
    try:
        run = subprocess.run(
            version_cmds.get(tool, [tool, "--version"]),
            text=True,
            capture_output=True,
            timeout=10,
        )
        text = (run.stdout or run.stderr or "").splitlines()
        first = text[0] if text else "FOUND_VERSION_UNREADABLE"
        return f"{path} | {first[:200]}"
    except Exception as exc:  # noqa: BLE001
        return f"{path} | VERSION_CHECK_FAILED: {exc}"


rows = read_tsv(CANDIDATE)

if len(rows) != EXPECTED_CANDIDATES:
    failures.append(f"Expected {EXPECTED_CANDIDATES} candidates; observed {len(rows)}.")

required = [
    "analysis_layer", "gene_symbol", "allele_key", "CHROM", "POS", "REF", "ALT",
    "BED_start_0based", "BED_end_0based_exclusive", "candidate_sources",
]
for col in required:
    if rows and col not in rows[0]:
        failures.append(f"Required candidate column missing: {col}")

seen_keys = set()
class_counts = Counter()
gene_counts = Counter()
source_counts = Counter()
layer_counts = Counter()
chrom_counts = Counter()
invalid_coordinate_rows = 0

for row in rows:
    key = (row["analysis_layer"], row["gene_symbol"], row["allele_key"])
    if key in seen_keys:
        failures.append(f"Duplicate candidate key observed: {key}")
        break
    seen_keys.add(key)

    pos = int(row["POS"])
    bed_start = int(row["BED_start_0based"])
    bed_end = int(row["BED_end_0based_exclusive"])

    if bed_start != pos - 1 or bed_end <= bed_start:
        invalid_coordinate_rows += 1

    cls = variant_class(row["REF"], row["ALT"])
    class_counts[cls] += 1
    gene_counts[row["gene_symbol"]] += 1
    layer_counts[row["analysis_layer"]] += 1
    chrom_counts[row["CHROM"]] += 1

    for source in row["candidate_sources"].split(";"):
        source_counts[source] += 1

if invalid_coordinate_rows:
    failures.append(f"Invalid BED coordinate rows observed: {invalid_coordinate_rows}.")

# Write compact context summaries.
def write_counter(counter: Counter, key_name: str, value_name: str, filename: str) -> None:
    with (OUTDIR / filename).open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow([key_name, value_name])
        for key, value in sorted(counter.items(), key=lambda x: (-x[1], x[0])):
            writer.writerow([key, value])


write_counter(class_counts, "variant_class", "candidate_count", "phaseC_candidate_variant_class_counts_v1.tsv")
write_counter(gene_counts, "gene_symbol", "candidate_count", "phaseC_candidate_gene_counts_v1.tsv")
write_counter(source_counts, "candidate_source", "candidate_count", "phaseC_candidate_source_counts_v1.tsv")
write_counter(layer_counts, "analysis_layer", "candidate_count", "phaseC_candidate_analysis_layer_counts_v1.tsv")
write_counter(chrom_counts, "CHROM", "candidate_count", "phaseC_candidate_chromosome_counts_v1.tsv")

# Create true VCF text with minimal headers, still as input only.
vcf_path = OUTDIR / "phaseC_frequency_variation_candidate_union_for_annotation_v1.vcf"
with vcf_path.open("w", encoding="utf-8") as handle:
    handle.write("##fileformat=VCFv4.2\n")
    handle.write("##source=Project6_PhaseC_annotation_readiness\n")
    handle.write("##INFO=<ID=GENE,Number=1,Type=String,Description=\"Candidate gene symbol\">\n")
    handle.write("##INFO=<ID=ANALYSIS_LAYER,Number=1,Type=String,Description=\"Project6 analysis layer\">\n")
    handle.write("##INFO=<ID=CANDIDATE_SOURCES,Number=.,Type=String,Description=\"Descriptive candidate source labels\">\n")
    handle.write("##INFO=<ID=MAX_AF_GROUP,Number=1,Type=String,Description=\"Group with maximum descriptive AF\">\n")
    handle.write("##INFO=<ID=MAX_MINUS_MIN_AF,Number=1,Type=Float,Description=\"Descriptive max-minus-min AF metric\">\n")
    handle.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
    for row in rows:
        info = (
            f"GENE={row['gene_symbol']};"
            f"ANALYSIS_LAYER={row['analysis_layer']};"
            f"CANDIDATE_SOURCES={row['candidate_sources'].replace(';', ',')};"
            f"MAX_AF_GROUP={row['max_AF_group']};"
            f"MAX_MINUS_MIN_AF={row['max_minus_min_AF']}"
        )
        handle.write(
            f"{row['CHROM']}\t{row['POS']}\t{row['allele_key']}\t"
            f"{row['REF']}\t{row['ALT']}\t.\tPASS\t{info}\n"
        )

# Tool audit.
tools = ["vep", "snpEff", "bcftools", "bedtools", "tabix", "bgzip"]
tool_rows = [{"tool": tool, "availability": tool_version(tool)} for tool in tools]
with (OUTDIR / "phaseC_candidate_annotation_tool_availability_v1.tsv").open("w", newline="", encoding="utf-8") as handle:
    writer = csv.DictWriter(handle, fieldnames=["tool", "availability"], delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(tool_rows)

summary = {
    "candidate_records": len(rows),
    "unique_candidate_keys": len(seen_keys),
    "variant_class_counts": dict(class_counts),
    "analysis_layer_counts": dict(layer_counts),
    "candidate_source_counts": dict(source_counts),
    "gene_count": len(gene_counts),
    "chromosome_count": len(chrom_counts),
    "invalid_coordinate_rows": invalid_coordinate_rows,
    "vep_available": shutil.which("vep") is not None,
    "snpeff_available": shutil.which("snpEff") is not None,
    "bcftools_available": shutil.which("bcftools") is not None,
    "bedtools_available": shutil.which("bedtools") is not None,
}

summary_path = OUTDIR / "phaseC_candidate_consequence_annotation_context_summary_v1.tsv"
with summary_path.open("w", newline="", encoding="utf-8") as handle:
    writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
    writer.writerow(["metric", "value"])
    for key, value in summary.items():
        writer.writerow([key, json.dumps(value) if isinstance(value, dict) else value])

payload = {
    "status": "PASS" if not failures else "FAIL",
    "operation": "Audit consequence-annotation context and create annotation-ready VCF input.",
    "input_candidate_table": str(CANDIDATE),
    "summary": summary,
    "outputs": {
        "candidate_variant_class_counts": str(OUTDIR / "phaseC_candidate_variant_class_counts_v1.tsv"),
        "candidate_gene_counts": str(OUTDIR / "phaseC_candidate_gene_counts_v1.tsv"),
        "candidate_source_counts": str(OUTDIR / "phaseC_candidate_source_counts_v1.tsv"),
        "candidate_analysis_layer_counts": str(OUTDIR / "phaseC_candidate_analysis_layer_counts_v1.tsv"),
        "candidate_chromosome_counts": str(OUTDIR / "phaseC_candidate_chromosome_counts_v1.tsv"),
        "annotation_tool_availability": str(OUTDIR / "phaseC_candidate_annotation_tool_availability_v1.tsv"),
        "vcf_input": str(vcf_path),
        "summary": str(summary_path),
    },
    "failures": failures,
    "boundary": (
        "This audits annotation context and prepares input only. It does not run annotation, "
        "assign consequence, prioritise function, infer selection or authorise biological interpretation."
    ),
}

OUTJSON.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
print(json.dumps(payload, indent=2))

if payload["status"] != "PASS":
    raise SystemExit(1)
