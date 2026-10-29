#!/usr/bin/env python3
"""Profile imported VEP consequence landscape for Phase C candidates.

This is a technical annotation-summary checkpoint. It summarizes imported VEP
annotations across genes, consequence classes, VEP impact categories, feature
types, MANE/APPRIS labels and candidate source layers. It does not rank variants,
select candidates, infer selection, or make biological interpretation.
"""

from __future__ import annotations

import csv
import json
from collections import Counter, defaultdict
from pathlib import Path

ANNOT = Path(
    "11_results/phaseC/genotype_only_qc/returned_vep_import_audit/"
    "phaseC_returned_vep_all_transcript_annotations_v1.tsv"
)
CANONICAL = Path(
    "11_results/phaseC/genotype_only_qc/returned_vep_import_audit/"
    "phaseC_returned_vep_MANE_or_APPRIS_canonical_like_annotations_v1.tsv"
)
PER_CANDIDATE = Path(
    "11_results/phaseC/genotype_only_qc/returned_vep_import_audit/"
    "phaseC_returned_vep_per_candidate_annotation_summary_v1.tsv"
)
CANDIDATE_METADATA = Path(
    "11_results/phaseC/genotype_only_qc/functional_annotation_readiness/"
    "phaseC_frequency_variation_candidate_union_for_annotation_v1.tsv"
)

OUTDIR = Path("11_results/phaseC/genotype_only_qc/vep_consequence_landscape")
OUTJSON = Path("results/qc/phaseC_vep_consequence_landscape_profile_v1.json")

OUTDIR.mkdir(parents=True, exist_ok=True)
OUTJSON.parent.mkdir(parents=True, exist_ok=True)

failures: list[str] = []


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_rows(path: Path, fieldnames: list[str], rows: list[dict[str, object]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row.get(key, ".") for key in fieldnames})


def write_counter(path: Path, key_name: str, count_name: str, counter: Counter) -> None:
    rows = [
        {key_name: key, count_name: value}
        for key, value in sorted(counter.items(), key=lambda x: (-x[1], x[0]))
    ]
    write_rows(path, [key_name, count_name], rows)


annot_rows = read_tsv(ANNOT)
canonical_rows = read_tsv(CANONICAL)
candidate_rows = read_tsv(CANDIDATE_METADATA)
per_candidate_rows = read_tsv(PER_CANDIDATE)

if len(candidate_rows) != 500:
    failures.append(f"Expected 500 candidate metadata rows; observed {len(candidate_rows)}.")
if len(per_candidate_rows) != 500:
    failures.append(f"Expected 500 per-candidate annotation rows; observed {len(per_candidate_rows)}.")
if len(annot_rows) != 13070:
    failures.append(f"Expected 13,070 all-transcript VEP rows; observed {len(annot_rows)}.")
if len(canonical_rows) != 2426:
    failures.append(f"Expected 2,426 canonical-like VEP rows; observed {len(canonical_rows)}.")

candidate_by_id = {row["allele_key"]: row for row in candidate_rows}
all_ids = {row["Uploaded_variation"] for row in annot_rows}
candidate_ids = set(candidate_by_id)


def first_present(row: dict[str, str], names: list[str], default: str = ".") -> str:
    for name in names:
        value = row.get(name, "")
        if value and value != ".":
            return value
    return default

if all_ids != candidate_ids:
    failures.append(
        f"Candidate/VEP uploaded-variation ID mismatch: {len(candidate_ids - all_ids)} missing and {len(all_ids - candidate_ids)} unexpected."
    )

# Basic landscape counts.
write_counter(
    OUTDIR / "phaseC_vep_landscape_all_transcript_gene_symbol_counts_v1.tsv",
    "gene_symbol",
    "vep_row_count",
    Counter(row.get("SYMBOL_from_Extra", ".") for row in annot_rows),
)

write_counter(
    OUTDIR / "phaseC_vep_landscape_all_transcript_impact_counts_v1.tsv",
    "impact",
    "vep_row_count",
    Counter(row.get("IMPACT_from_Extra", ".") for row in annot_rows),
)

write_counter(
    OUTDIR / "phaseC_vep_landscape_all_transcript_feature_type_counts_v1.tsv",
    "feature_type",
    "vep_row_count",
    Counter(row.get("Feature_type", ".") for row in annot_rows),
)

write_counter(
    OUTDIR / "phaseC_vep_landscape_all_transcript_consequence_counts_v1.tsv",
    "consequence",
    "vep_row_count",
    Counter(row.get("Consequence", ".") for row in annot_rows),
)

write_counter(
    OUTDIR / "phaseC_vep_landscape_canonical_like_gene_symbol_counts_v1.tsv",
    "gene_symbol",
    "canonical_like_vep_row_count",
    Counter(row.get("SYMBOL_from_Extra", ".") for row in canonical_rows),
)

write_counter(
    OUTDIR / "phaseC_vep_landscape_canonical_like_impact_counts_v1.tsv",
    "impact",
    "canonical_like_vep_row_count",
    Counter(row.get("IMPACT_from_Extra", ".") for row in canonical_rows),
)

write_counter(
    OUTDIR / "phaseC_vep_landscape_canonical_like_consequence_counts_v1.tsv",
    "consequence",
    "canonical_like_vep_row_count",
    Counter(row.get("Consequence", ".") for row in canonical_rows),
)

# Gene x impact and gene x consequence tables.
gene_impact = Counter(
    (row.get("SYMBOL_from_Extra", "."), row.get("IMPACT_from_Extra", "."))
    for row in annot_rows
)
gene_impact_rows = [
    {"gene_symbol": gene, "impact": impact, "vep_row_count": count}
    for (gene, impact), count in sorted(gene_impact.items(), key=lambda x: (x[0][0], x[0][1]))
]
write_rows(
    OUTDIR / "phaseC_vep_landscape_all_transcript_gene_by_impact_v1.tsv",
    ["gene_symbol", "impact", "vep_row_count"],
    gene_impact_rows,
)

gene_consequence = Counter(
    (row.get("SYMBOL_from_Extra", "."), row.get("Consequence", "."))
    for row in annot_rows
)
gene_consequence_rows = [
    {"gene_symbol": gene, "consequence": consequence, "vep_row_count": count}
    for (gene, consequence), count in sorted(gene_consequence.items(), key=lambda x: (x[0][0], x[0][1]))
]
write_rows(
    OUTDIR / "phaseC_vep_landscape_all_transcript_gene_by_consequence_v1.tsv",
    ["gene_symbol", "consequence", "vep_row_count"],
    gene_consequence_rows,
)

# Per-candidate technical consequence landscape joined to candidate source metadata.
rows_by_id = defaultdict(list)
canonical_by_id = defaultdict(list)
for row in annot_rows:
    rows_by_id[row["Uploaded_variation"]].append(row)
for row in canonical_rows:
    canonical_by_id[row["Uploaded_variation"]].append(row)

candidate_summary_rows = []
for allele_key in sorted(candidate_ids):
    meta = candidate_by_id[allele_key]
    rows = rows_by_id[allele_key]
    canonical_like = canonical_by_id.get(allele_key, [])
    impacts = sorted({row.get("IMPACT_from_Extra", ".") for row in rows})
    consequences = sorted({row.get("Consequence", ".") for row in rows})
    canonical_impacts = sorted({row.get("IMPACT_from_Extra", ".") for row in canonical_like})
    canonical_consequences = sorted({row.get("Consequence", ".") for row in canonical_like})
    symbols = sorted({row.get("SYMBOL_from_Extra", ".") for row in rows if row.get("SYMBOL_from_Extra", ".") != "."})
    symbol_union = ",".join(symbols) if symbols else "."
    candidate_summary_rows.append(
        {
            "allele_key": allele_key,
            "gene": first_present(
                meta,
                ["gene", "target_gene", "gene_symbol", "symbol", "locus", "locus_gene"],
                symbol_union,
            ),
            "candidate_source": first_present(
                meta,
                ["candidate_source", "source", "candidate_set", "candidate_origin", "selection_source"],
                ".",
            ),
            "analysis_layer": first_present(meta, ["analysis_layer", "layer"], "."),
            "variant_class": first_present(meta, ["variant_class", "class", "variant_type"], "."),
            "all_transcript_row_count": len(rows),
            "canonical_like_row_count": len(canonical_like),
            "all_transcript_impact_union": ",".join(impacts),
            "all_transcript_consequence_union": ",".join(consequences),
            "canonical_like_impact_union": ",".join(canonical_impacts) if canonical_impacts else ".",
            "canonical_like_consequence_union": ",".join(canonical_consequences) if canonical_consequences else ".",
            "symbol_union": symbol_union,
            "has_moderate_or_high_impact_any_transcript": "YES" if any(i in {"MODERATE", "HIGH"} for i in impacts) else "NO",
            "has_moderate_or_high_impact_canonical_like": "YES" if any(i in {"MODERATE", "HIGH"} for i in canonical_impacts) else "NO",
            "has_mane_or_appris_canonical_like_annotation": "YES" if canonical_like else "NO",
        }
    )

candidate_summary_path = OUTDIR / "phaseC_vep_landscape_per_candidate_technical_summary_v1.tsv"
write_rows(
    candidate_summary_path,
    [
        "allele_key",
        "gene",
        "candidate_source",
        "analysis_layer",
        "variant_class",
        "all_transcript_row_count",
        "canonical_like_row_count",
        "all_transcript_impact_union",
        "all_transcript_consequence_union",
        "canonical_like_impact_union",
        "canonical_like_consequence_union",
        "symbol_union",
        "has_moderate_or_high_impact_any_transcript",
        "has_moderate_or_high_impact_canonical_like",
        "has_mane_or_appris_canonical_like_annotation",
    ],
    candidate_summary_rows,
)

# Candidate-level count summaries.
source_impact = Counter(
    (row["candidate_source"], row["has_moderate_or_high_impact_canonical_like"])
    for row in candidate_summary_rows
)
source_impact_rows = [
    {
        "candidate_source": source,
        "has_moderate_or_high_impact_canonical_like": flag,
        "candidate_count": count,
    }
    for (source, flag), count in sorted(source_impact.items(), key=lambda x: (x[0][0], x[0][1]))
]
write_rows(
    OUTDIR / "phaseC_vep_landscape_candidate_source_by_canonical_like_impact_flag_v1.tsv",
    ["candidate_source", "has_moderate_or_high_impact_canonical_like", "candidate_count"],
    source_impact_rows,
)

gene_candidate_flags = Counter(
    (row["gene"], row["has_moderate_or_high_impact_canonical_like"])
    for row in candidate_summary_rows
)
gene_candidate_flag_rows = [
    {
        "gene": gene,
        "has_moderate_or_high_impact_canonical_like": flag,
        "candidate_count": count,
    }
    for (gene, flag), count in sorted(gene_candidate_flags.items(), key=lambda x: (x[0][0], x[0][1]))
]
write_rows(
    OUTDIR / "phaseC_vep_landscape_gene_by_canonical_like_impact_flag_v1.tsv",
    ["gene", "has_moderate_or_high_impact_canonical_like", "candidate_count"],
    gene_candidate_flag_rows,
)

moderate_or_high_any = sum(
    row["has_moderate_or_high_impact_any_transcript"] == "YES" for row in candidate_summary_rows
)
moderate_or_high_canonical = sum(
    row["has_moderate_or_high_impact_canonical_like"] == "YES" for row in candidate_summary_rows
)
canonical_available = sum(
    row["has_mane_or_appris_canonical_like_annotation"] == "YES" for row in candidate_summary_rows
)

profile_summary_rows = [
    {"metric": "submitted_candidate_records", "value": len(candidate_rows)},
    {"metric": "all_transcript_vep_rows", "value": len(annot_rows)},
    {"metric": "canonical_like_vep_rows", "value": len(canonical_rows)},
    {"metric": "per_candidate_summary_rows", "value": len(candidate_summary_rows)},
    {"metric": "candidates_with_canonical_like_annotation", "value": canonical_available},
    {"metric": "candidates_with_moderate_or_high_impact_any_transcript", "value": moderate_or_high_any},
    {"metric": "candidates_with_moderate_or_high_impact_canonical_like", "value": moderate_or_high_canonical},
]
write_rows(
    OUTDIR / "phaseC_vep_consequence_landscape_profile_summary_v1.tsv",
    ["metric", "value"],
    profile_summary_rows,
)

payload = {
    "status": "PASS" if not failures else "FAIL",
    "operation": "Profile imported VEP consequence landscape for Phase C candidates.",
    "candidate_records": len(candidate_rows),
    "all_transcript_vep_rows": len(annot_rows),
    "canonical_like_vep_rows": len(canonical_rows),
    "per_candidate_summary_rows": len(candidate_summary_rows),
    "candidates_with_canonical_like_annotation": canonical_available,
    "candidates_with_moderate_or_high_impact_any_transcript": moderate_or_high_any,
    "candidates_with_moderate_or_high_impact_canonical_like": moderate_or_high_canonical,
    "outputs": {
        "profile_summary": str(OUTDIR / "phaseC_vep_consequence_landscape_profile_summary_v1.tsv"),
        "per_candidate_technical_summary": str(candidate_summary_path),
        "all_transcript_gene_counts": str(OUTDIR / "phaseC_vep_landscape_all_transcript_gene_symbol_counts_v1.tsv"),
        "all_transcript_consequence_counts": str(OUTDIR / "phaseC_vep_landscape_all_transcript_consequence_counts_v1.tsv"),
        "all_transcript_impact_counts": str(OUTDIR / "phaseC_vep_landscape_all_transcript_impact_counts_v1.tsv"),
        "canonical_like_consequence_counts": str(OUTDIR / "phaseC_vep_landscape_canonical_like_consequence_counts_v1.tsv"),
        "canonical_like_impact_counts": str(OUTDIR / "phaseC_vep_landscape_canonical_like_impact_counts_v1.tsv"),
        "gene_by_impact": str(OUTDIR / "phaseC_vep_landscape_all_transcript_gene_by_impact_v1.tsv"),
        "gene_by_consequence": str(OUTDIR / "phaseC_vep_landscape_all_transcript_gene_by_consequence_v1.tsv"),
    },
    "failures": failures,
    "boundary": (
        "This checkpoint profiles imported VEP consequences as technical annotation summaries only. "
        "It does not rank variants, select candidates, infer selection, make population-differentiation claims, "
        "or authorise biological interpretation."
    ),
}

OUTJSON.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
print(json.dumps(payload, indent=2))

if payload["status"] != "PASS":
    raise SystemExit(1)
