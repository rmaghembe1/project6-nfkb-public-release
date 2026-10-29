#!/usr/bin/env python3
"""Prepare a manual screening worklist for the formal Project 6 competitor audit."""

from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path


def normalize_title(title: str) -> str:
    title = title.lower()
    title = re.sub(r"<[^>]+>", "", title)
    title = re.sub(r"[^a-z0-9 ]+", " ", title)
    title = re.sub(r"\s+", " ", title).strip()
    return title


def priority_hint(title: str, abstract: str, query_ids: str) -> str:
    text = f"{title} {abstract}".lower()

    if "alternatively spliced nfkb1" in text and "andean" in text:
        return "HIGH_NEAR_COMPETITOR_NFKB1_POPULATION_FUNCTIONAL"

    if (
        "denisovan, modern human and mouse tnfaip3 alleles" in text
        or "phospho-tuning immunity through denisovan" in text
    ):
        return "HIGH_FOUNDATIONAL_TNFAIP3_ANCHOR"

    if "mapping the gene regulatory landscape of archaic hominin introgression" in text:
        return "HIGH_PARTIAL_OVERLAP_ARCHAIC_REGULATORY_CONTEXT"

    if "environmental and genetic disease modifiers of haploinsufficiency of a20" in text:
        return "HIGH_SUPPORTING_ANCHOR_EXTENSION"

    if "archaic hominin admixture facilitated adaptation" in text:
        return "HIGH_BROADER_EVOLUTIONARY_CONTEXT_REVIEW"

    if any(term in text for term in [
        "nfkb1", "nfkb2", "ikbkg", "tnfaip3", "a20", "tnip1", "cyld"
    ]) and any(term in text for term in [
        "population", "adaptation", "selection", "haplotype",
        "introgress", "denisovan", "evolution", "eqtl", "enhancer"
    ]):
        return "MEDIUM_REVIEW_FOR_PARTIAL_OVERLAP"

    if any(term in text for term in [
        "sheep", "partridge", "pig", "fox", "carassius",
        "eleutheronema", "mice", "mouse"
    ]):
        return "LOW_LIKELY_NONHUMAN_EXCLUSION"

    if "mendelian randomization" in text or "drug target" in text:
        return "LOW_LIKELY_FUNCTIONAL_BACKGROUND_OR_EXCLUSION"

    return "MEDIUM_SCREEN_ABSTRACT"

def family_hint(title: str) -> str:
    norm = normalize_title(title)

    if "mapping the gene regulatory landscape of archaic hominin introgression in modern papuans" in norm:
        return "PF001_PAPUAN_ARCHAIC_INTROGRESSION_PREPRINT_AND_PUBLICATION"
    if "alternatively spliced nfkb1 transcripts enriched in andean" in norm:
        return "PF002_NFKB1_ANDEAN_FUNCTIONAL_BOUNDARY"
    if "denisovan modern human and mouse tnfaip3 alleles tune a20 phosphorylation and immunity" in norm:
        return "PF003_TNFAIP3_A20_DENISOVAN_ANCHOR"
    if "phospho tuning immunity through denisovan modern human and mouse tnfaip3 gene variants" in norm:
        return "PF003_TNFAIP3_A20_DENISOVAN_ANCHOR"
    return ""

def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()

    with args.input.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))

    output_rows = []
    for index, row in enumerate(rows, start=1):
        output_rows.append({
            "screening_id": f"C{index:03d}",
            "publication_family_id": family_hint(row["title"]),
            "triage_priority_hint": priority_hint(
                row["title"], row["abstract"], row["query_ids"]
            ),
            "query_ids": row["query_ids"],
            "source_databases": row["source_databases"],
            "year": row["year"],
            "title": row["title"],
            "journal": row["journal"],
            "doi": row["doi"],
            "pmid": row["pmid"],
            "pmcid": row["pmcid"],
            "abstract": row["abstract"],
            "human_population_relevance": "",
            "nfkb_circuitry_scope": "",
            "population_genomic_or_evolutionary_scope": "",
            "functional_annotation_overlap": "",
            "full_text_check_required": "",
            "screening_classification": "",
            "screening_reason": "",
            "retain_in_literature_manifest": "",
        })

    fields = list(output_rows[0].keys())
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fields, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(output_rows)

    print(f"Prepared screening worklist for {len(output_rows)} unique candidate records.")
    print("Priority hints are triage aids only; manual screening classification remains required.")


if __name__ == "__main__":
    main()
