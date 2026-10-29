#!/usr/bin/env python3
"""Audit and import returned VEP annotations for Phase C candidates.

This checkpoint imports returned VEP annotation as structured technical outputs.
It validates record identity and schema, separates all VEP rows from MANE/canonical-
like subsets, and does not make biological interpretation.
"""

from __future__ import annotations

import csv
import json
from collections import Counter, defaultdict
from pathlib import Path

SUBMITTED = Path(
    "11_results/phaseC/genotype_only_qc/functional_annotation_readiness/"
    "phaseC_frequency_variation_candidate_union_for_annotation_v1.tsv"
)
VEP = Path(
    "10_external_returns/phaseC_candidate_vep_annotation_2026-06-02/"
    "x1BMV83ZybBYkQ1Z.vep.txt"
)
METADATA_NOTE = Path(
    "10_external_returns/phaseC_candidate_vep_annotation_2026-06-02/"
    "annotation_run_metadata_manual_entry.md"
)

OUTDIR = Path("11_results/phaseC/genotype_only_qc/returned_vep_import_audit")
OUTJSON = Path("results/qc/phaseC_returned_vep_annotation_import_audit_v1.json")

EXPECTED_SUBMITTED = 500

OUTDIR.mkdir(parents=True, exist_ok=True)
OUTJSON.parent.mkdir(parents=True, exist_ok=True)

failures: list[str] = []


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def parse_extra(extra: str) -> dict[str, str]:
    out = {}
    if not extra or extra == "-":
        return out
    for item in extra.split(";"):
        if not item:
            continue
        if "=" in item:
            key, value = item.split("=", 1)
            out[key] = value
        else:
            out[item] = "YES"
    return out


def clean_value(value) -> str:
    """Return a whitespace-safe TSV cell; use '.' for blank cells."""
    if value is None:
        return "."
    value = str(value).strip()
    return value if value else "."


def write_tsv_rows(path: Path, fieldnames: list[str], rows: list[dict[str, str]]) -> None:
    """Write rectangular TSV without trailing blank delimiters."""
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(fieldnames)
        for row in rows:
            writer.writerow([clean_value(row.get(field, ".")) for field in fieldnames])


submitted_rows = read_tsv(SUBMITTED)
submitted_ids = {row["allele_key"] for row in submitted_rows}

if len(submitted_rows) != EXPECTED_SUBMITTED:
    failures.append(f"Expected {EXPECTED_SUBMITTED} submitted rows; observed {len(submitted_rows)}.")
if len(submitted_ids) != EXPECTED_SUBMITTED:
    failures.append(f"Expected {EXPECTED_SUBMITTED} unique submitted allele keys; observed {len(submitted_ids)}.")

# Read VEP output. Header begins with #Uploaded_variation.
with VEP.open(encoding="utf-8") as handle:
    lines = [line.rstrip("\n") for line in handle if line.rstrip("\n")]

header_lines = [line for line in lines if line.startswith("#")]
data_lines = [line for line in lines if not line.startswith("#")]

if not header_lines:
    failures.append("No VEP header line found.")
    header = []
else:
    header = header_lines[-1].lstrip("#").split("\t")

required_cols = [
    "Uploaded_variation", "Location", "Allele", "Gene", "Feature", "Feature_type",
    "Consequence", "Existing_variation", "Extra",
]
for col in required_cols:
    if col not in header:
        failures.append(f"Required VEP column missing: {col}")

vep_rows = []
if header:
    for raw in data_lines:
        parts = raw.split("\t")
        if len(parts) != len(header):
            failures.append(f"VEP row field-count mismatch: expected {len(header)}, observed {len(parts)}.")
            break
        row = {key: value.strip() for key, value in zip(header, parts)}
        extra = parse_extra(row.get("Extra", ""))
        row["SYMBOL_from_Extra"] = extra.get("SYMBOL", "")
        row["IMPACT_from_Extra"] = extra.get("IMPACT", "")
        row["BIOTYPE_from_Extra"] = extra.get("BIOTYPE", "")
        row["MANE_from_Extra"] = extra.get("MANE", "")
        row["MANE_SELECT_from_Extra"] = extra.get("MANE_SELECT", "")
        row["APPRIS_from_Extra"] = extra.get("APPRIS", "")
        row["TSL_from_Extra"] = extra.get("TSL", "")
        row["HGVSc_from_Extra"] = extra.get("HGVSc", "")
        row["HGVSp_from_Extra"] = extra.get("HGVSp", "")
        row["CLIN_SIG_from_Extra"] = extra.get("CLIN_SIG", "")
        row["PHENO_from_Extra"] = extra.get("PHENO", "")
        row["AF_from_Extra"] = extra.get("AF", "")
        vep_rows.append(row)

vep_ids = {row["Uploaded_variation"] for row in vep_rows}
missing_ids = sorted(submitted_ids - vep_ids)
extra_ids = sorted(vep_ids - submitted_ids)

if missing_ids:
    failures.append(f"Submitted IDs missing from VEP output: {len(missing_ids)}.")
if extra_ids:
    failures.append(f"VEP output contains unexpected uploaded IDs: {len(extra_ids)}.")

# Write all imported rows with parsed Extra fields.
imported_path = OUTDIR / "phaseC_returned_vep_all_transcript_annotations_v1.tsv"
fields = header + [
    "SYMBOL_from_Extra",
    "IMPACT_from_Extra",
    "BIOTYPE_from_Extra",
    "MANE_from_Extra",
    "MANE_SELECT_from_Extra",
    "APPRIS_from_Extra",
    "TSL_from_Extra",
    "HGVSc_from_Extra",
    "HGVSp_from_Extra",
    "CLIN_SIG_from_Extra",
    "PHENO_from_Extra",
    "AF_from_Extra",
]
write_tsv_rows(imported_path, fields, vep_rows)

# MANE/canonical-like rows. We are not choosing “best” consequence; just extracting labelled rows.
mane_rows = [
    row for row in vep_rows
    if row.get("MANE_from_Extra") or row.get("MANE_SELECT_from_Extra")
]
canonical_like_rows = [
    row for row in vep_rows
    if row.get("MANE_from_Extra") or row.get("MANE_SELECT_from_Extra") or row.get("APPRIS_from_Extra") in {"P1", "A1"}
]

mane_path = OUTDIR / "phaseC_returned_vep_MANE_labelled_annotations_v1.tsv"
canonical_like_path = OUTDIR / "phaseC_returned_vep_MANE_or_APPRIS_canonical_like_annotations_v1.tsv"

write_tsv_rows(mane_path, fields, mane_rows)
write_tsv_rows(canonical_like_path, fields, canonical_like_rows)

# Per-submitted-variant summary.
by_id = defaultdict(list)
for row in vep_rows:
    by_id[row["Uploaded_variation"]].append(row)

summary_rows = []
for allele_key in sorted(submitted_ids):
    rows = by_id.get(allele_key, [])
    consequences = sorted({row["Consequence"] for row in rows})
    impacts = sorted({row["IMPACT_from_Extra"] for row in rows if row["IMPACT_from_Extra"]})
    symbols = sorted({row["SYMBOL_from_Extra"] for row in rows if row["SYMBOL_from_Extra"]})
    feature_types = sorted({row["Feature_type"] for row in rows})
    summary_rows.append({
        "allele_key": allele_key,
        "vep_row_count": len(rows),
        "feature_type_count": len(feature_types),
        "consequence_terms_union": ",".join(consequences),
        "impact_terms_union": ",".join(impacts),
        "symbol_terms_union": ",".join(symbols),
        "has_MANE_labelled_row": "YES" if any(r.get("MANE_from_Extra") or r.get("MANE_SELECT_from_Extra") for r in rows) else "NO",
        "has_APPRIS_P1_or_A1_row": "YES" if any(r.get("APPRIS_from_Extra") in {"P1", "A1"} for r in rows) else "NO",
        "has_existing_variation": "YES" if any(r.get("Existing_variation") and r.get("Existing_variation") != "-" for r in rows) else "NO",
    })

variant_summary_path = OUTDIR / "phaseC_returned_vep_per_candidate_annotation_summary_v1.tsv"
write_tsv_rows(variant_summary_path, list(summary_rows[0]), summary_rows)

def sanitise_tsv_file(path: Path) -> None:
    """Rewrite TSV with stripped cells and '.' for blank cells to avoid trailing tabs."""
    rows = []
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.reader(handle, delimiter="\t")
        for row in reader:
            rows.append(row)

    if not rows:
        return

    width = max(len(row) for row in rows)
    cleaned_rows = []
    for row in rows:
        padded = row + [""] * (width - len(row))
        cleaned_rows.append([(cell.strip() if cell.strip() else ".") for cell in padded])

    # Keep the header unchanged except stripping whitespace.
    cleaned_rows[0] = [cell.strip() for cell in rows[0]] + [""] * (width - len(rows[0]))

    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerows(cleaned_rows)


def sanitise_returned_vep_tsv_outputs() -> None:
    for path in OUTDIR.glob("phaseC_returned_vep_*_v1.tsv"):
        sanitise_tsv_file(path)


# Compact count tables.

def write_counter(counter: Counter, key_name: str, value_name: str, filename: str) -> None:
    with (OUTDIR / filename).open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow([key_name, value_name])
        for key, value in sorted(counter.items(), key=lambda x: (-x[1], x[0])):
            writer.writerow([key, value])


write_counter(Counter(row["Consequence"] for row in vep_rows), "consequence", "vep_row_count", "phaseC_returned_vep_consequence_counts_v1.tsv")
write_counter(Counter(row["IMPACT_from_Extra"] or "MISSING" for row in vep_rows), "impact", "vep_row_count", "phaseC_returned_vep_impact_counts_v1.tsv")
write_counter(Counter(row["Feature_type"] for row in vep_rows), "feature_type", "vep_row_count", "phaseC_returned_vep_feature_type_counts_v1.tsv")
write_counter(Counter(row["SYMBOL_from_Extra"] or "MISSING" for row in vep_rows), "symbol", "vep_row_count", "phaseC_returned_vep_symbol_counts_v1.tsv")

payload = {
    "status": "PASS" if not failures else "FAIL",
    "operation": "Audit and import returned VEP candidate annotations.",
    "input_vep_file": str(VEP),
    "submitted_candidate_records": len(submitted_rows),
    "returned_vep_rows": len(vep_rows),
    "unique_returned_uploaded_variation_ids": len(vep_ids),
    "missing_submitted_ids": len(missing_ids),
    "unexpected_returned_ids": len(extra_ids),
    "mane_labelled_rows": len(mane_rows),
    "mane_or_appris_canonical_like_rows": len(canonical_like_rows),
    "per_candidate_summary_rows": len(summary_rows),
    "metadata_note_exists": METADATA_NOTE.exists(),
    "outputs": {
        "all_transcript_annotations": str(imported_path),
        "MANE_labelled_annotations": str(mane_path),
        "MANE_or_APPRIS_canonical_like_annotations": str(canonical_like_path),
        "per_candidate_summary": str(variant_summary_path),
        "consequence_counts": str(OUTDIR / "phaseC_returned_vep_consequence_counts_v1.tsv"),
        "impact_counts": str(OUTDIR / "phaseC_returned_vep_impact_counts_v1.tsv"),
        "feature_type_counts": str(OUTDIR / "phaseC_returned_vep_feature_type_counts_v1.tsv"),
        "symbol_counts": str(OUTDIR / "phaseC_returned_vep_symbol_counts_v1.tsv"),
    },
    "failures": failures,
    "boundary": (
        "This imports returned VEP annotations as technical annotation outputs only. It does not select causal "
        "variants, prioritise functional importance, infer selection, make population-differentiation claims "
        "or authorise biological interpretation."
    ),
}

OUTJSON.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
print(json.dumps(payload, indent=2))

if payload["status"] != "PASS":
    raise SystemExit(1)
