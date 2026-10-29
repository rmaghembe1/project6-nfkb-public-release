#!/usr/bin/env python3
"""Refine broadly triaged OpenAlex citation neighbours into a high-specificity screening set."""

from __future__ import annotations

import argparse
import csv
import json
import re
from collections import Counter
from pathlib import Path


GENE_PATTERNS = {
    "NFKB1": re.compile(r"\bNFKB1\b", re.IGNORECASE),
    "NFKB2": re.compile(r"\bNFKB2\b", re.IGNORECASE),
    "RELA": re.compile(r"\bRELA\b"),
    "RELB": re.compile(r"\bRELB\b"),
    "REL": re.compile(r"\bREL\b"),
    "IKBKG": re.compile(r"\bIKBKG\b", re.IGNORECASE),
    "CHUK": re.compile(r"\bCHUK\b", re.IGNORECASE),
    "IKBKB": re.compile(r"\bIKBKB\b", re.IGNORECASE),
    "NFKBIA": re.compile(r"\bNFKBIA\b", re.IGNORECASE),
    "NFKBIE": re.compile(r"\bNFKBIE\b", re.IGNORECASE),
    "IKBKE": re.compile(r"\bIKBKE\b", re.IGNORECASE),
    "TNFAIP3": re.compile(r"\bTNFAIP3\b", re.IGNORECASE),
    "TNIP1": re.compile(r"\bTNIP1\b", re.IGNORECASE),
    "CYLD": re.compile(r"\bCYLD\b", re.IGNORECASE),
    "OTULIN": re.compile(r"\bOTULIN\b", re.IGNORECASE),
}

PATHWAY_PATTERNS = {
    "NF-kB pathway": re.compile(
        r"\bNF\s*[-‐-–]?\s*[κkK](?:appa)?\s*B\b|\bNF\s*[-‐-–]?\s*kB\b",
        re.IGNORECASE,
    ),
    "nuclear factor kappa B": re.compile(
        r"\bnuclear factor\s*[-‐-–]?\s*(?:κ|kappa)\s*B\b",
        re.IGNORECASE,
    ),
}

EVOLUTIONARY_PATTERNS = {
    "population genetics": re.compile(r"\bpopulation genetic(?:s)?\b", re.IGNORECASE),
    "population genomics": re.compile(r"\bpopulation genomic(?:s)?\b", re.IGNORECASE),
    "human population": re.compile(r"\bhuman population(?:s)?\b", re.IGNORECASE),
    "global population": re.compile(r"\bglobal population(?:s)?\b", re.IGNORECASE),
    "ancestry": re.compile(r"\bancestr(?:y|al)\b", re.IGNORECASE),
    "human evolution": re.compile(r"\bhuman evolution(?:ary)?\b", re.IGNORECASE),
    "positive selection": re.compile(r"\bpositive selection\b", re.IGNORECASE),
    "natural selection": re.compile(r"\bnatural selection\b", re.IGNORECASE),
    "local adaptation": re.compile(r"\blocal adaptation\b", re.IGNORECASE),
    "adaptive introgression": re.compile(r"\badaptive introgression\b", re.IGNORECASE),
    "introgression": re.compile(r"\bintrogress(?:ion|ed|ed)?\b", re.IGNORECASE),
    "Denisovan": re.compile(r"\bDenisovan\b", re.IGNORECASE),
    "archaic": re.compile(r"\barchaic\b", re.IGNORECASE),
    "haplotype": re.compile(r"\bhaplotype(?:s)?\b", re.IGNORECASE),
    "eQTL": re.compile(r"\beQTL(?:s)?\b", re.IGNORECASE),
    "response QTL": re.compile(r"\bresponse\s+(?:eQTL|QTL)(?:s)?\b", re.IGNORECASE),
    "regulatory variation": re.compile(
        r"\bregulatory (?:variant|variants|variation|polymorphism|polymorphisms)\b",
        re.IGNORECASE,
    ),
    "copy number variation": re.compile(
        r"\bcopy number (?:variant|variants|variation|variations)\b|\bCNV(?:s)?\b",
        re.IGNORECASE,
    ),
    "structural variation": re.compile(
        r"\bstructural (?:variant|variants|variation|variations)\b",
        re.IGNORECASE,
    ),
}

CLINICAL_GENETIC_PATTERNS = {
    "genetic association": re.compile(r"\bgenetic association\b", re.IGNORECASE),
    "gene-gene interaction": re.compile(
        r"\bgene\s*[-–]?\s*gene interaction(?:s)?\b|\bgenetic interaction(?:s)?\b",
        re.IGNORECASE,
    ),
    "polymorphism": re.compile(r"\bpolymorphism(?:s)?\b", re.IGNORECASE),
    "case-control": re.compile(r"\bcase\s*[-–]?\s*control\b", re.IGNORECASE),
    "susceptibility": re.compile(r"\bsusceptibilit(?:y|ies)\b", re.IGNORECASE),
    "risk": re.compile(r"\brisk\b", re.IGNORECASE),
}


def compact_text(value: str) -> str:
    return re.sub(r"\s+", " ", value or "").strip()


def matched_terms(patterns: dict[str, re.Pattern], text: str) -> list[str]:
    return sorted(name for name, pattern in patterns.items() if pattern.search(text))


def classify_for_screening(text: str) -> tuple[str, list[str], list[str], list[str], str]:
    genes = matched_terms(GENE_PATTERNS, text)
    pathway = matched_terms(PATHWAY_PATTERNS, text)
    evolution = matched_terms(EVOLUTIONARY_PATTERNS, text)
    clinical = matched_terms(CLINICAL_GENETIC_PATTERNS, text)

    strong_nfkb = bool(genes or pathway)
    tnfaip3_anchor = "TNFAIP3" in genes
    multi_gene = len(genes) >= 2

    if strong_nfkb and evolution:
        tier = "TIER1_EVOLUTIONARY_POPULATION_FUNCTIONAL_OVERLAP"
        rationale = (
            "Specific NF-kB pathway/gene signal with explicit evolutionary, population, "
            "haplotype, regulatory-variation, eQTL or structural-variation term."
        )
    elif multi_gene and clinical:
        tier = "TIER2_MULTI_GENE_CLINICAL_GENETIC_BOUNDARY"
        rationale = (
            "At least two NF-kB circuit genes with clinical genetic-association or "
            "gene-interaction terminology."
        )
    elif tnfaip3_anchor and clinical:
        tier = "TIER3_TNFAIP3_FUNCTIONAL_DISEASE_BOUNDARY"
        rationale = (
            "TNFAIP3/A20 locus-specific genetic or functional disease evidence; "
            "retained only as an interpretative boundary."
        )
    else:
        tier = "NOT_RETAINED_BY_HIGH_SPECIFICITY_FILTER"
        rationale = "No retained high-specificity Project 6 competitor-boundary pattern."

    return tier, genes + pathway, evolution, clinical, rationale


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-unique", required=True, type=Path)
    parser.add_argument("--out-all-triaged", required=True, type=Path)
    parser.add_argument("--out-retained", required=True, type=Path)
    parser.add_argument("--out-tier1-tier2", required=True, type=Path)
    parser.add_argument("--out-json", required=True, type=Path)
    args = parser.parse_args()

    with args.input_unique.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))

    triaged = []
    retained = []
    tier1_tier2 = []

    for row in rows:
        title = compact_text(row["title"])
        abstract = compact_text(row["abstract"])
        text = f"{title} {abstract}"
        tier, nfkb_terms, evolutionary_terms, clinical_terms, rationale = (
            classify_for_screening(text)
        )

        out = dict(row)
        out["title"] = title
        out["abstract"] = abstract
        out["high_specificity_tier"] = tier
        out["matched_specific_nfkb_terms"] = ";".join(nfkb_terms)
        out["matched_evolutionary_or_population_terms"] = ";".join(evolutionary_terms)
        out["matched_clinical_genetic_terms"] = ";".join(clinical_terms)
        out["high_specificity_rationale"] = rationale
        out["manual_screening_classification"] = ""
        out["manual_screening_reason"] = ""
        out["direct_competitor"] = ""
        triaged.append(out)

        if tier != "NOT_RETAINED_BY_HIGH_SPECIFICITY_FILTER":
            retained.append(out)
        if tier in {
            "TIER1_EVOLUTIONARY_POPULATION_FUNCTIONAL_OVERLAP",
            "TIER2_MULTI_GENE_CLINICAL_GENETIC_BOUNDARY",
        }:
            tier1_tier2.append(out)

    fields = list(triaged[0].keys())
    for path, records in [
        (args.out_all_triaged, triaged),
        (args.out_retained, retained),
        (args.out_tier1_tier2, tier1_tier2),
    ]:
        path.parent.mkdir(parents=True, exist_ok=True)
        with path.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(
                handle, fieldnames=fields, delimiter="\t", lineterminator="\n"
            )
            writer.writeheader()
            writer.writerows(records)

    counts = Counter(row["high_specificity_tier"] for row in triaged)
    payload = {
        "status": "HIGH_SPECIFICITY_CITATION_TRIAGE_COMPLETE_MANUAL_SCREENING_PENDING",
        "input_unique_neighbour_records": len(rows),
        "high_specificity_tier_counts": dict(sorted(counts.items())),
        "retained_candidate_records_all_tiers": len(retained),
        "priority_manual_screening_records_tier1_tier2": len(tier1_tier2),
        "method_boundary": (
            "This is a deterministic specificity refinement of an already retrieved "
            "citation-neighbour set. It does not add or remove retrieved evidence and "
            "does not establish novelty without manual classification."
        ),
    }
    args.out_json.parent.mkdir(parents=True, exist_ok=True)
    args.out_json.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
