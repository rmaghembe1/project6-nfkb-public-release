#!/usr/bin/env python3
"""Classify genuinely new high-priority OpenAlex citation-chain candidates for Project 6."""

from __future__ import annotations

import argparse
import csv
import json
from collections import Counter
from pathlib import Path


DECISIONS = {
    "CC002": {
        "classification": "SUPPORTING_FUNCTIONAL_BACKGROUND",
        "direct_competitor": "no",
        "retain_in_boundary_manifest": "no",
        "full_text_check_required": "no",
        "reason": (
            "High-altitude adaptation commentary/review mentioning NFKB1 among multiple "
            "potential contributors to haemoglobin adaptation; not a primary NF-kB "
            "population-genomic circuitry analysis."
        ),
    },
    "CC004": {
        "classification": "PARTIAL_OVERLAP_CONTEXT",
        "direct_competitor": "no",
        "retain_in_boundary_manifest": "no",
        "full_text_check_required": "no",
        "reason": (
            "Human SLE case-control association study analysing IKBKE and NFKBIA with "
            "BANK1; represents clinical genetic overlap, not evolutionary circuitry."
        ),
    },
    "CC005": {
        "classification": "SUPPORTING_FUNCTIONAL_BACKGROUND",
        "direct_competitor": "no",
        "retain_in_boundary_manifest": "no",
        "full_text_check_required": "no",
        "reason": (
            "Review of HIF/NF-kB crosstalk and high-altitude inflammatory physiology; "
            "contextual relevance only, without germline NF-kB circuit analysis."
        ),
    },
    "CC008": {
        "classification": "EXCLUDE_NON_PRIMARY_COMMENTARY",
        "direct_competitor": "no",
        "retain_in_boundary_manifest": "no",
        "full_text_check_required": "no",
        "reason": (
            "Research Highlights commentary summarising the published TNFAIP3/A20 "
            "Denisovan anchor paper already captured as C033."
        ),
    },
    "CC010": {
        "classification": "EXCLUDE_DISEASE_ASSOCIATION_ONLY",
        "direct_competitor": "no",
        "retain_in_boundary_manifest": "no",
        "full_text_check_required": "no",
        "reason": (
            "Meta-analysis of TNFAIP3 and TNIP1 polymorphisms in SLE risk; "
            "clinical susceptibility evidence without evolutionary-genomics scope."
        ),
    },
    "CC011": {
        "classification": "EXCLUDE_DISEASE_ASSOCIATION_ONLY",
        "direct_competitor": "no",
        "retain_in_boundary_manifest": "no",
        "full_text_check_required": "no",
        "reason": (
            "TNIP1 gastric-carcinoma case-control SNP/haplotype association study; "
            "single-locus disease-risk design."
        ),
    },
    "CC012": {
        "classification": "EXCLUDE_DISEASE_ASSOCIATION_ONLY",
        "direct_competitor": "no",
        "retain_in_boundary_manifest": "no",
        "full_text_check_required": "no",
        "reason": (
            "TNFAIP3 case-control susceptibility/haplotype analysis in Korean RA/SLE; "
            "single-locus disease-association design."
        ),
    },
    "CC014": {
        "classification": "PARTIAL_OVERLAP_CONTEXT",
        "direct_competitor": "no",
        "retain_in_boundary_manifest": "yes",
        "full_text_check_required": "no",
        "reason": (
            "Human Behçet disease case-control analysis jointly evaluating NFKB1 and "
            "NFKBIA polymorphisms; multi-gene canonical NF-kB clinical-genetic boundary "
            "without evolutionary or global population-genomic design."
        ),
    },
    "CC015": {
        "classification": "PARTIAL_OVERLAP_CONTEXT",
        "direct_competitor": "no",
        "retain_in_boundary_manifest": "no",
        "full_text_check_required": "no",
        "reason": (
            "Multi-locus SLE risk-interaction study including TNFAIP3 and TNIP1 among "
            "broader immune loci; clinical association rather than NF-kB evolution."
        ),
    },
    "CC016": {
        "classification": "EXCLUDE_DISEASE_ASSOCIATION_ONLY",
        "direct_competitor": "no",
        "retain_in_boundary_manifest": "no",
        "full_text_check_required": "no",
        "reason": (
            "TNIP1 SNP/haplotype and expression study in Chinese Han SLE; "
            "single-locus clinical susceptibility evidence."
        ),
    },
    "CC017": {
        "classification": "EXCLUDE_NONHUMAN_OR_NONRELEVANT",
        "direct_competitor": "no",
        "retain_in_boundary_manifest": "no",
        "full_text_check_required": "no",
        "reason": (
            "Thymocyte-selection mechanism study; the term selection refers to immune "
            "cell development rather than population evolution."
        ),
    },
    "CC018": {
        "classification": "EXCLUDE_DISEASE_ASSOCIATION_ONLY",
        "direct_competitor": "no",
        "retain_in_boundary_manifest": "no",
        "full_text_check_required": "no",
        "reason": (
            "IKBKE and IFIH1 SLE risk-variant/haplotype study; functional disease "
            "genetics without NF-kB circuit evolutionary analysis."
        ),
    },
    "CC019": {
        "classification": "SUPPORTING_FUNCTIONAL_BACKGROUND",
        "direct_competitor": "no",
        "retain_in_boundary_manifest": "yes",
        "full_text_check_required": "no",
        "reason": (
            "Multi-ancestry TNIP1 functional risk-haplotype study in SLE; supports "
            "population-aware restraint-arm regulatory interpretation but is disease focused."
        ),
    },
    "CC020": {
        "classification": "PARTIAL_OVERLAP_CONTEXT",
        "direct_competitor": "no",
        "retain_in_boundary_manifest": "no",
        "full_text_check_required": "no",
        "reason": (
            "Chinese SLE replication study involving TNFAIP3, TNIP1 and other "
            "susceptibility genes; clinical association design rather than pathway evolution."
        ),
    },
    "CC021": {
        "classification": "PARTIAL_OVERLAP_CONTEXT",
        "direct_competitor": "no",
        "retain_in_boundary_manifest": "yes",
        "full_text_check_required": "no",
        "reason": (
            "Human SLE gene-gene interaction study including TNFAIP3 and REL within "
            "NF-kB signalling; important clinical multi-gene boundary, but not an "
            "evolutionary-population circuitry analysis."
        ),
    },
    "CC022": {
        "classification": "SUPPORTING_FUNCTIONAL_BACKGROUND",
        "direct_competitor": "no",
        "retain_in_boundary_manifest": "no",
        "full_text_check_required": "no",
        "reason": (
            "Japanese TNIP1 SLE association study assessing relation to TNFAIP3; "
            "supports NF-kB restraint biology but is clinical susceptibility research."
        ),
    },
    "CC025": {
        "classification": "NEAR_COMPETITOR_NOVELTY_BOUNDARY",
        "direct_competitor": "no",
        "retain_in_boundary_manifest": "yes",
        "full_text_check_required": "yes",
        "reason": (
            "Human HapMap population analysis of geographically differentiated SNPs "
            "and pathway-level evolutionary variation reporting NF-kB activation among "
            "implicated pathways; closest newly identified population-pathway boundary, "
            "but abstract does not describe a predefined NF-kB regulatory-circuit design."
        ),
    },
    "CC026": {
        "classification": "PARTIAL_OVERLAP_CONTEXT",
        "direct_competitor": "no",
        "retain_in_boundary_manifest": "yes",
        "full_text_check_required": "no",
        "reason": (
            "Coeliac-disease GWAS/follow-up identifies heritable variation at TNFAIP3 "
            "and REL within NF-kB signalling; multi-locus disease-risk boundary rather "
            "than evolutionary circuitry analysis."
        ),
    },
    "CC027": {
        "classification": "PARTIAL_OVERLAP_CONTEXT",
        "direct_competitor": "no",
        "retain_in_boundary_manifest": "yes",
        "full_text_check_required": "no",
        "reason": (
            "Human NHL case-control study genotyping TNF/LTA and nine canonical NF-kB "
            "genes, including CHUK, IKBKB, NFKB1 and REL; strong multi-gene clinical "
            "pathway-genetic boundary but not population-resolved evolutionary genomics."
        ),
    },
}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--summary", required=True, type=Path)
    args = parser.parse_args()

    with args.input.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))

    observed = {row["citation_screening_id"] for row in rows}
    expected = set(DECISIONS)
    if observed != expected:
        raise SystemExit(
            f"ERROR: citation screening-ID mismatch. "
            f"Missing={sorted(expected-observed)}; unexpected={sorted(observed-expected)}"
        )

    classified = []
    for row in rows:
        decision = DECISIONS[row["citation_screening_id"]]
        row["citation_manual_classification"] = decision["classification"]
        row["citation_manual_reason"] = decision["reason"]
        row["direct_competitor"] = decision["direct_competitor"]
        row["retain_in_boundary_manifest"] = decision["retain_in_boundary_manifest"]
        row["full_text_check_required"] = decision["full_text_check_required"]
        classified.append(row)

    fields = list(classified[0].keys())
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fields, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(classified)

    counts = Counter(row["citation_manual_classification"] for row in classified)
    direct = [row for row in classified if row["direct_competitor"] == "yes"]
    retained = [
        {
            "citation_screening_id": row["citation_screening_id"],
            "classification": row["citation_manual_classification"],
            "title": row["title"],
            "doi": row["doi"],
        }
        for row in classified
        if row["retain_in_boundary_manifest"] == "yes"
    ]
    full_text = [
        {
            "citation_screening_id": row["citation_screening_id"],
            "classification": row["citation_manual_classification"],
            "title": row["title"],
            "doi": row["doi"],
        }
        for row in classified
        if row["full_text_check_required"] == "yes"
    ]

    payload = {
        "status": "NEW_CITATION_PRIORITY_SCREENING_COMPLETE_TARGETED_VERIFICATION_PENDING",
        "new_priority_citation_records_screened": len(classified),
        "classification_counts": dict(sorted(counts.items())),
        "direct_competitor_count_in_new_citation_priority_records": len(direct),
        "direct_competitor_ids": [row["citation_screening_id"] for row in direct],
        "records_recommended_for_boundary_manifest": retained,
        "new_targeted_full_text_verification_records": full_text,
        "interpretation": (
            "No direct competitor was identified among the genuinely new citation-priority "
            "records. CC025 introduces an important human population/pathway-evolution "
            "boundary and requires targeted verification before formal closure. Several "
            "records establish prior multi-gene NF-kB clinical-genetic studies, refining "
            "the novelty claim away from any generic first-multi-gene assertion."
        ),
    }

    args.summary.parent.mkdir(parents=True, exist_ok=True)
    args.summary.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
