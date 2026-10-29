#!/usr/bin/env python3
"""Classify genuinely new supplementary competitor-search records for Project 6."""

from __future__ import annotations

import argparse
import csv
import json
from collections import Counter
from pathlib import Path


DECISIONS = {
    "SP002": {
        "classification": "EXCLUDE_NONHUMAN_OR_NONRELEVANT",
        "direct_competitor": "no",
        "retain_for_final_manifest": "no",
        "reason": (
            "KIRC tumour bioinformatics/prognostic study using NF-kB-related gene "
            "expression, CNV and somatic-variant data; not human population "
            "evolutionary or germline circuitry analysis."
        ),
    },
    "SP004": {
        "classification": "EXCLUDE_NONHUMAN_OR_NONRELEVANT",
        "direct_competitor": "no",
        "retain_for_final_manifest": "no",
        "reason": (
            "DLBCL tumour-genetic classification study in which NF-kB-pathway "
            "mutations occur as disease molecular features; not population evolution."
        ),
    },
    "SP005": {
        "classification": "PARTIAL_OVERLAP_CONTEXT",
        "direct_competitor": "no",
        "retain_for_final_manifest": "yes",
        "reason": (
            "Human Chinese case-control rheumatoid-arthritis study genotyping multiple "
            "NF-kB signalling-pathway SNPs and testing genetic interactions; establishes "
            "a multi-gene clinical-association boundary but lacks evolutionary, global "
            "population, regulatory-haplotype or selection analysis."
        ),
    },
    "SP006": {
        "classification": "EXCLUDE_DISEASE_ASSOCIATION_ONLY",
        "direct_competitor": "no",
        "retain_for_final_manifest": "no",
        "reason": (
            "Single-cell transcriptomic analysis of DADA2 monocytes identifying "
            "noncanonical NF-kB activation; no population-genetic or evolutionary design."
        ),
    },
    "SP008": {
        "classification": "EXCLUDE_NONHUMAN_OR_NONRELEVANT",
        "direct_competitor": "no",
        "retain_for_final_manifest": "no",
        "reason": (
            "Rectal-cancer exome study characterising somatic PCDHB3 effects on "
            "NF-kB transcriptional activity; unrelated to germline population evolution."
        ),
    },
    "SP009": {
        "classification": "EXCLUDE_DISEASE_ASSOCIATION_ONLY",
        "direct_competitor": "no",
        "retain_for_final_manifest": "no",
        "reason": (
            "Single HA20 case report with a pathogenic TNFAIP3 splice mutation and "
            "functional inflammation assays; not a population-resolved competitor."
        ),
    },
    "SP010": {
        "classification": "EXCLUDE_NONHUMAN_OR_NONRELEVANT",
        "direct_competitor": "no",
        "retain_for_final_manifest": "no",
        "reason": (
            "Conference poster-session compilation retrieved through term overlap; "
            "no identifiable competitor study relevant to Project 6."
        ),
    },
    "SP011": {
        "classification": "PARTIAL_OVERLAP_CONTEXT",
        "direct_competitor": "no",
        "retain_for_final_manifest": "yes",
        "reason": (
            "Human Taiwanese quantitative-trait association study assessing variants "
            "across several NF-kB pathway genes and sICAM1 levels; relevant multi-gene "
            "population-association boundary but not evolutionary genomics."
        ),
    },
    "SP012": {
        "classification": "SUPPORTING_FUNCTIONAL_BACKGROUND",
        "direct_competitor": "no",
        "retain_for_final_manifest": "no",
        "reason": (
            "Mouse/mechanistic study showing GIT2 recruits CYLD to terminate "
            "TLR-induced NF-kB signalling; supports restraint biology but is not a "
            "human population-genomic competitor."
        ),
    },
    "SP013": {
        "classification": "PARTIAL_OVERLAP_CONTEXT",
        "direct_competitor": "no",
        "retain_for_final_manifest": "yes",
        "reason": (
            "Human Chinese SLE case-control study testing multiple NF-kB signalling "
            "genes, including NFKB1, IKBKB, TNFAIP3 and TNIP1, with significant "
            "gene-gene interactions; important multi-gene clinical-genetic boundary "
            "but lacks population-evolutionary circuit analysis."
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

    observed = {row["supplementary_id"] for row in rows}
    expected = set(DECISIONS)
    if observed != expected:
        raise SystemExit(
            f"ERROR: supplementary screening-ID mismatch. "
            f"Missing={sorted(expected-observed)}; unexpected={sorted(observed-expected)}"
        )

    classified = []
    for row in rows:
        decision = DECISIONS[row["supplementary_id"]]
        row["supplementary_screening_classification"] = decision["classification"]
        row["supplementary_screening_reason"] = decision["reason"]
        row["direct_competitor"] = decision["direct_competitor"]
        row["retain_for_final_manifest"] = decision["retain_for_final_manifest"]
        classified.append(row)

    fields = list(classified[0].keys())
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fields, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(classified)

    counts = Counter(
        row["supplementary_screening_classification"] for row in classified
    )
    direct = [
        row for row in classified if row["direct_competitor"] == "yes"
    ]
    retained = [
        {
            "supplementary_id": row["supplementary_id"],
            "classification": row["supplementary_screening_classification"],
            "title": row["title"],
            "doi": row["doi"],
        }
        for row in classified
        if row["retain_for_final_manifest"] == "yes"
    ]

    payload = {
        "status": "SUPPLEMENTARY_NEW_RECORD_SCREENING_COMPLETE_CITATION_CHAINING_PENDING",
        "new_supplementary_records_screened": len(classified),
        "classification_counts": dict(sorted(counts.items())),
        "direct_competitor_count_in_new_supplementary_records": len(direct),
        "direct_competitor_screening_ids": [
            row["supplementary_id"] for row in direct
        ],
        "records_recommended_for_novelty_boundary_manifest": retained,
        "interpretation": (
            "No direct competitor was identified among the genuinely new "
            "supplementary records. Three studies provide additional boundaries "
            "showing multi-gene human NF-kB genetic association in clinical or "
            "quantitative-trait settings, without a global population-resolved "
            "evolutionary-functional circuitry design."
        ),
        "remaining_audit_requirement": (
            "Perform bounded citation-chaining/manual reference check from central "
            "novelty-boundary and anchor publications before formal Gate 1 closure."
        ),
    }

    args.summary.parent.mkdir(parents=True, exist_ok=True)
    args.summary.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
