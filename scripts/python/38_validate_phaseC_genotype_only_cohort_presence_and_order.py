#!/usr/bin/env python3
"""Validate Phase C genotype-only cohort membership and cross-extract order invariance."""

from __future__ import annotations

import csv
import json
import subprocess
from pathlib import Path

COHORT = Path(
    "11_results/phaseB/harmonisation_qc/"
    "cohort_primary_genotype_metadata_linked_n4150_v1.txt"
)
GENODIR = Path("03_raw_genomic_resources/interval_extraction/full_panel/genotype")
REFERENCE = GENODIR / "TNFAIP3_GRCh38_gnomad_hgdp_tgp.vcf.gz"
OUTDIR = Path("11_results/phaseC/genotype_only_qc")
OUTJSON = Path("results/qc/phaseC_genotype_only_cohort_presence_order_validation_v1.json")
EXCLUDED = OUTDIR / "phaseC_genotype_resource_samples_excluded_from_primary_cohort_v1.txt"
TABLE = OUTDIR / "phaseC_genotype_extract_cohort_presence_and_order_v1.tsv"

EXPECTED_COHORT_COUNT = 4150
EXPECTED_FULL_RESOURCE_COUNT = 4151
EXPECTED_EXTRACT_COUNT = 24
EXPECTED_EXCLUDED_SAMPLE = "CHMI_CHMI3_WGS2"


def query_samples(path: Path, subset: Path | None = None) -> list[str]:
    cmd = ["bcftools", "query", "-l"]
    if subset is not None:
        cmd.extend(["-S", str(subset)])
    cmd.append(str(path))
    result = subprocess.run(cmd, check=True, text=True, capture_output=True)
    return [line for line in result.stdout.splitlines() if line]


OUTDIR.mkdir(parents=True, exist_ok=True)
OUTJSON.parent.mkdir(parents=True, exist_ok=True)

cohort_order = [line for line in COHORT.read_text(encoding="utf-8").splitlines() if line]
cohort_members = set(cohort_order)

reference_full_order = query_samples(REFERENCE)
reference_subset_order = query_samples(REFERENCE, COHORT)
reference_subset_members = set(reference_subset_order)
reference_excluded = sorted(set(reference_full_order) - cohort_members)
EXCLUDED.write_text("\n".join(reference_excluded) + ("\n" if reference_excluded else ""), encoding="utf-8")

files = sorted(GENODIR.glob("*.vcf.gz"))
failures: list[str] = []
rows: list[dict[str, object]] = []

if len(cohort_order) != EXPECTED_COHORT_COUNT or len(cohort_members) != EXPECTED_COHORT_COUNT:
    failures.append("Frozen genotype-only cohort does not contain exactly 4,150 unique samples.")
if len(files) != EXPECTED_EXTRACT_COUNT:
    failures.append(f"Expected 24 genotype extracts; observed {len(files)}.")
if reference_subset_members != cohort_members:
    failures.append("TNFAIP3 reference extract membership differs from the frozen genotype-only cohort.")
if reference_excluded != [EXPECTED_EXCLUDED_SAMPLE]:
    failures.append(
        f"Expected only excluded resource sample {EXPECTED_EXCLUDED_SAMPLE}; observed {reference_excluded}."
    )

for path in files:
    gene = path.name.split("_", 1)[0]
    full_order = query_samples(path)
    subset_order = query_samples(path, COHORT)
    subset_members = set(subset_order)

    membership_pass = subset_members == cohort_members and len(subset_order) == EXPECTED_COHORT_COUNT
    order_pass = subset_order == reference_subset_order
    status = (
        "PASS"
        if len(full_order) == EXPECTED_FULL_RESOURCE_COUNT and membership_pass and order_pass
        else "FAIL"
    )

    if status != "PASS":
        failures.append(f"{gene}: cohort presence/order validation failed.")

    rows.append({
        "gene_symbol": gene,
        "genotype_extract_samples_total": len(full_order),
        "cohort_samples_retrieved": len(subset_order),
        "membership_matches_frozen_cohort": "PASS" if membership_pass else "FAIL",
        "order_matches_genotype_reference": "PASS" if order_pass else "FAIL",
        "status": status,
    })

with TABLE.open("w", newline="", encoding="utf-8") as handle:
    writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)

payload = {
    "status": "PASS" if not failures else "FAIL",
    "operation": "Phase C genotype-only cohort membership and genotype-extract order-invariance validation.",
    "phaseC_scope_boundary": (
        "Genotype-only technical cohort validation only; phased/interoperability filters are not applied."
    ),
    "reference_genotype_extract": str(REFERENCE),
    "frozen_primary_genotype_cohort_count": len(cohort_members),
    "genotype_extracts_validated": len(files),
    "expected_resource_samples_per_extract": EXPECTED_FULL_RESOURCE_COUNT,
    "expected_primary_cohort_samples_per_extract": EXPECTED_COHORT_COUNT,
    "excluded_resource_samples": reference_excluded,
    "extracts_failing_validation": [row["gene_symbol"] for row in rows if row["status"] != "PASS"],
    "failures": failures,
    "interpretation_boundary": (
        "This validation does not authorise genotype-only variant filtering, "
        "allele-frequency analysis or biological interpretation."
    ),
}
OUTJSON.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
print(json.dumps(payload, indent=2))

if payload["status"] != "PASS":
    raise SystemExit(1)
