#!/usr/bin/env python3
"""Assign manually curated abstract-level competitor classifications for Project 6."""

from __future__ import annotations

import argparse
import csv
import json
from collections import Counter
from pathlib import Path


DECISIONS = {
    "C001": ("PF001_PAPUAN_ARCHAIC_INTROGRESSION_PREPRINT_AND_PUBLICATION", "PARTIAL_OVERLAP_CONTEXT", "yes", "yes", "Published Papuan archaic-regulatory study; assays introgressed regulatory variants and identifies predicted TNFAIP3 effects, but does not describe pathway-wide NF-kB evolutionary circuitry analysis."),
    "C002": ("", "SUPPORTING_FUNCTIONAL_BACKGROUND", "no", "no", "Regulatory-variant prioritisation pipeline in autoimmune disease; includes TNFAIP3 context but is not population-evolutionary NF-kB circuitry analysis."),
    "C003": ("", "EXCLUDE_NONHUMAN_OR_NONRELEVANT", "no", "no", "Tibetan sheep adaptation study; nonhuman."),
    "C004": ("", "EXCLUDE_DISEASE_ASSOCIATION_ONLY", "no", "no", "Celiac-disease candidate polymorphism association study; no evolutionary NF-kB circuitry scope."),
    "C005": ("PF001_PAPUAN_ARCHAIC_INTROGRESSION_PREPRINT_AND_PUBLICATION", "SUPERSEDED_PREPRINT", "no", "no", "Preprint predecessor of the published Papuan archaic-regulatory study C001."),
    "C006": ("", "EXCLUDE_NONHUMAN_OR_NONRELEVANT", "no", "no", "Colorectal-cancer multi-omics study with pathway enrichment only; no NF-kB population-evolutionary design."),
    "C007": ("", "EXCLUDE_NONHUMAN_OR_NONRELEVANT", "no", "no", "Human CNS-disease MR/single-cell study; NF-kB appears as tumour-cell state context, not pathway evolution."),
    "C008": ("PF002_NFKB1_ANDEAN_FUNCTIONAL_BOUNDARY", "NEAR_COMPETITOR_NOVELTY_BOUNDARY", "yes", "yes", "Published human WGS/transcriptomic population-functional study of NFKB1 alternative splicing in Andean Aymara; closest locus-specific boundary but not circuitry-wide."),
    "C009": ("", "EXCLUDE_NONHUMAN_OR_NONRELEVANT", "no", "no", "Sepsis MR/trans-omics study; NF-kB is downstream inflammatory context only."),
    "C010": ("", "SUPPORTING_FUNCTIONAL_BACKGROUND", "no", "no", "TNFAIP3 emerges from disease-target MR/colocalisation; useful background but not competitor design."),
    "C011": ("", "EXCLUDE_DISEASE_ASSOCIATION_ONLY", "no", "no", "Clinical WES study of Behçet-like disease and TNFAIP3 loss-of-function variants; no population-evolutionary scope."),
    "C012": ("", "PARTIAL_OVERLAP_CONTEXT", "no", "yes", "Pathway-oriented regulatory-variant/eQTL workflow, but applied to BCR signalling in Europeans rather than NF-kB evolutionary circuitry."),
    "C013": ("", "EXCLUDE_NONHUMAN_OR_NONRELEVANT", "no", "no", "Atherosclerosis MR/transcriptomic study with NF-kB enrichment only."),
    "C014": ("", "EXCLUDE_NONHUMAN_OR_NONRELEVANT", "no", "no", "Tibetan partridge adaptation study; nonhuman."),
    "C015": ("", "EXCLUDE_DISEASE_ASSOCIATION_ONLY", "no", "no", "TNIP1 association study in infant RSV disease; no evolutionary circuitry scope."),
    "C016": ("", "EXCLUDE_NONHUMAN_OR_NONRELEVANT", "no", "no", "Fish population-genomics study; nonhuman."),
    "C017": ("", "EXCLUDE_DISEASE_ASSOCIATION_ONLY", "no", "no", "Psoriasis susceptibility and immunophenotype study of A20/CYLD/Cezanne; no evolutionary genomics."),
    "C018": ("PF004_TNFAIP3_SLE_ENHANCER_PREPRINT_AND_PUBLICATION", "SUPPORTING_FUNCTIONAL_BACKGROUND", "no", "yes", "Published functional dissection of a TNFAIP3 enhancer risk haplotype; directly supports regulatory-annotation rationale but is disease-locus focused."),
    "C019": ("", "EXCLUDE_NONHUMAN_OR_NONRELEVANT", "no", "no", "Bone mineral-density RNAm-SNP study with pathway enrichment only."),
    "C020": ("", "EXCLUDE_DISEASE_ASSOCIATION_ONLY", "no", "no", "TNFAIP3 rheumatoid-arthritis association study in Mexicans; no evolutionary circuit design."),
    "C021": ("", "SUPPORTING_FUNCTIONAL_BACKGROUND", "no", "yes", "A20 haploinsufficiency modifier study involving introgressed I207L; mechanistic anchor extension rather than competitor."),
    "C022": ("", "EXCLUDE_NONHUMAN_OR_NONRELEVANT", "no", "no", "Mouse behavioural knockdown study with NF-kB readouts; irrelevant to human pathway evolution."),
    "C023": ("", "SUPPORTING_FUNCTIONAL_BACKGROUND", "no", "no", "Cord-blood eQTL/asthma study referencing TNFAIP3/TNIP2 expression; stimulus-regulatory context only."),
    "C024": ("", "EXCLUDE_NONHUMAN_OR_NONRELEVANT", "no", "no", "Nutrigenomic review/bioinformatic synthesis; RELA/NFKB1 are inferred regulators, not population-evolution targets."),
    "C025": ("", "SUPPORTING_FUNCTIONAL_BACKGROUND", "no", "yes", "Single-cell eQTL analysis of SLE-risk loci highlights TNFAIP3 in monocytes; supports cell-context annotation but not evolutionary genomics."),
    "C026": ("", "EXCLUDE_NONHUMAN_OR_NONRELEVANT", "no", "no", "Arctic fox immunogenetic selection study; nonhuman."),
    "C027": ("", "EXCLUDE_DISEASE_ASSOCIATION_ONLY", "no", "no", "TNFAIP3 haplotype association with scleritis in Chinese Han; no evolutionary circuitry analysis."),
    "C028": ("PF004_TNFAIP3_SLE_ENHANCER_PREPRINT_AND_PUBLICATION", "SUPERSEDED_PREPRINT", "no", "no", "Preprint predecessor of published TNFAIP3 enhancer study C018."),
    "C029": ("", "EXCLUDE_DISEASE_ASSOCIATION_ONLY", "no", "no", "TNFAIP3/TNIP1 ankylosing-spondylitis association study; no evolutionary circuitry analysis."),
    "C030": ("PF003_TNFAIP3_A20_DENISOVAN_ANCHOR", "SUPERSEDED_PREPRINT", "no", "no", "Preprint predecessor of published TNFAIP3/A20 Denisovan anchor paper C033."),
    "C031": ("", "EXCLUDE_DISEASE_ASSOCIATION_ONLY", "no", "no", "SLE candidate-allele association study in a Dominican cohort; no pathway-evolution design."),
    "C032": ("", "SUPPORTING_FUNCTIONAL_BACKGROUND", "no", "no", "Genome-wide lncRNA-mRNA eQTL module study with NF-kB pathway enrichment; generic functional-annotation context."),
    "C033": ("PF003_TNFAIP3_A20_DENISOVAN_ANCHOR", "FOUNDATIONAL_ANCHOR", "yes", "yes", "Published TNFAIP3/A20 T108A-I207L Denisovan introgression and functional-immunity anchor paper."),
    "C034": ("", "EXCLUDE_DISEASE_ASSOCIATION_ONLY", "no", "no", "Clinical septic-shock outcome association across selected NF-kB genes; no population-evolutionary design."),
    "C035": ("", "SUPPORTING_FUNCTIONAL_BACKGROUND", "no", "yes", "Stimulus-responsive chromatin resource identifies condition-specific TNFAIP3 regulatory variation; supports functional annotation strategy."),
    "C036": ("", "EXCLUDE_NONHUMAN_OR_NONRELEVANT", "no", "no", "Meishan pig selective-sweep study; nonhuman."),
    "C037": ("", "EXCLUDE_NONHUMAN_OR_NONRELEVANT", "no", "no", "Carassius adaptation expression/selection study; nonhuman."),
    "C038": ("", "EXCLUDE_NONHUMAN_OR_NONRELEVANT", "no", "no", "NRAMP1 regulatory review; NF-kB/RELA are contextual regulators rather than project competitors."),
    "C039": ("", "EXCLUDE_DISEASE_ASSOCIATION_ONLY", "no", "no", "TNFAIP3 SLE association in Han Chinese; no evolutionary circuit design."),
    "C040": ("", "SUPPORTING_FUNCTIONAL_BACKGROUND", "no", "yes", "Stimulated-fibroblast gene-environment functional association at TNFAIP3; supports response-sensitive annotation context."),
    "C041": ("", "EXCLUDE_DISEASE_ASSOCIATION_ONLY", "no", "no", "Korean SLE GWAS/replication study with TNFAIP3 secondary signal; no evolutionary circuitry scope."),
    "C042": ("", "EXCLUDE_NONHUMAN_OR_NONRELEVANT", "no", "no", "Mouse stress-response Tnfaip3 cis-eQTL/network study; nonhuman."),
    "C043": ("", "PARTIAL_OVERLAP_CONTEXT", "yes", "yes", "Human adaptive-introgression survey includes TNFAIP3 among immune-related archaic haplotypes and eQTLs; broader evolutionary context but not NF-kB circuitry analysis."),
    "C044": ("", "EXCLUDE_DISEASE_ASSOCIATION_ONLY", "no", "no", "Behçet disease immunogenetics review; TNFAIP3 appears among disease susceptibility loci."),
    "C045": ("", "SUPPORTING_FUNCTIONAL_BACKGROUND", "no", "yes", "TNFAIP3 psoriasis risk-haplotype meta-analysis with regulatory annotation; locus-specific functional background."),
    "C046": ("", "EXCLUDE_DISEASE_ASSOCIATION_ONLY", "no", "no", "Psoriasis susceptibility fine-mapping; TNFAIP3 is one disease locus among several."),
    "C047": ("", "EXCLUDE_DISEASE_ASSOCIATION_ONLY", "no", "no", "Multiple-sclerosis case-control association at TNFAIP3/TNFRSF1A; no evolutionary circuitry scope."),
    "C048": ("", "EXCLUDE_DISEASE_ASSOCIATION_ONLY", "no", "no", "Rheumatoid-arthritis association at TNFAIP3 in Chinese Han; no evolutionary circuitry scope."),
    "C049": ("", "EXCLUDE_DISEASE_ASSOCIATION_ONLY", "no", "no", "Graves disease TNFAIP3 association study; no evolutionary circuitry scope."),
    "C050": ("", "PARTIAL_OVERLAP_CONTEXT", "yes", "yes", "Human population-differentiated CNV pathway-evolution study mentioning NF-kB; requires full-text verification as conceptual boundary."),
    "C051": ("", "PARTIAL_OVERLAP_CONTEXT", "yes", "yes", "Human multi-gene NF-kB pathway disease-association study including IKBKE, NFKB1, TNIP1 and TNFAIP3; not evolutionary by abstract but relevant design boundary."),
    "C052": ("", "EXCLUDE_DISEASE_ASSOCIATION_ONLY", "no", "no", "Human Graves disease FCRL locus study; NF-kB appears as transcription-factor binding mechanism only."),
    "C053": ("", "EXCLUDE_DISEASE_ASSOCIATION_ONLY", "no", "no", "TNFAIP3 pharmacogenetic response association in psoriasis; no evolutionary circuitry scope."),
    "C054": ("", "EXCLUDE_DISEASE_ASSOCIATION_ONLY", "no", "no", "TNFAIP3 rheumatoid-arthritis susceptibility association in African Americans; no circuit-level evolutionary design."),
    "C055": ("", "SUPPORTING_FUNCTIONAL_BACKGROUND", "no", "yes", "TNFAIP3 sequencing/coding-variant analysis across autoimmune disease; functional locus background, not competitor."),
    "C056": ("", "EXCLUDE_DISEASE_ASSOCIATION_ONLY", "no", "no", "SLE immunodysregulation study; TNFAIP3 expression is secondary context."),
    "C057": ("", "SUPPORTING_FUNCTIONAL_BACKGROUND", "yes", "yes", "Functional TNFAIP3 downstream regulatory haplotype defined across ethnically diverse SLE cohorts; important locus-regulatory background."),
    "C058": ("", "EXCLUDE_DISEASE_ASSOCIATION_ONLY", "no", "no", "TNFAIP3 systemic-sclerosis association study; no evolutionary circuitry scope."),
    "C059": ("", "PARTIAL_OVERLAP_CONTEXT", "yes", "yes", "African-derived TNFAIP3 coding haplotypes with functional and structural assessment; single-locus population-functional boundary."),
    "C060": ("", "EXCLUDE_DISEASE_ASSOCIATION_ONLY", "no", "no", "Wegener granulomatosis candidate association; TNFAIP3 signal nominal and not evolutionary."),
    "C061": ("", "EXCLUDE_DISEASE_ASSOCIATION_ONLY", "no", "no", "Rheumatoid-arthritis association study focused on TRAF1, with TNFAIP3 contextual only."),
    "C062": ("", "EXCLUDE_DISEASE_ASSOCIATION_ONLY", "no", "no", "TNFAIP3 SLE risk-haplotype meta-analysis; disease-locus association rather than evolutionary analysis."),
    "C063": ("", "EXCLUDE_DISEASE_ASSOCIATION_ONLY", "no", "no", "TNFAIP3/6q23 rheumatoid-arthritis association study; no evolutionary circuitry scope."),
}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--summary", required=True, type=Path)
    args = parser.parse_args()

    with args.input.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))

    observed_ids = {row["screening_id"] for row in rows}
    expected_ids = set(DECISIONS)
    if observed_ids != expected_ids:
        raise SystemExit(
            f"ERROR: screening ID mismatch. Missing={sorted(expected_ids-observed_ids)}; "
            f"unexpected={sorted(observed_ids-expected_ids)}"
        )

    classified = []
    for row in rows:
        family, classification, full_text, retain, reason = DECISIONS[row["screening_id"]]
        row["publication_family_id"] = family or row["publication_family_id"]
        row["human_population_relevance"] = (
            "yes" if classification in {
                "FOUNDATIONAL_ANCHOR",
                "NEAR_COMPETITOR_NOVELTY_BOUNDARY",
                "PARTIAL_OVERLAP_CONTEXT",
                "SUPPORTING_FUNCTIONAL_BACKGROUND",
                "SUPERSEDED_PREPRINT",
            } and "nonhuman" not in reason.lower() else "no"
        )
        row["nfkb_circuitry_scope"] = (
            "partial_or_single_locus"
            if classification in {
                "FOUNDATIONAL_ANCHOR",
                "NEAR_COMPETITOR_NOVELTY_BOUNDARY",
                "PARTIAL_OVERLAP_CONTEXT",
                "SUPPORTING_FUNCTIONAL_BACKGROUND",
            }
            else "absent_or_not_relevant"
        )
        row["population_genomic_or_evolutionary_scope"] = (
            "relevant_boundary"
            if classification in {
                "FOUNDATIONAL_ANCHOR",
                "NEAR_COMPETITOR_NOVELTY_BOUNDARY",
                "PARTIAL_OVERLAP_CONTEXT",
            }
            else "not_competitor_scope"
        )
        row["full_text_check_required"] = full_text
        row["screening_classification"] = classification
        row["screening_reason"] = reason
        row["retain_in_literature_manifest"] = retain
        classified.append(row)

    fields = list(classified[0].keys())
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fields, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(classified)

    counts = Counter(row["screening_classification"] for row in classified)
    full_text_subset = [
        {
            "screening_id": row["screening_id"],
            "title": row["title"],
            "doi": row["doi"],
            "classification": row["screening_classification"],
        }
        for row in classified
        if row["full_text_check_required"] == "yes"
    ]
    direct_competitors = [
        row["screening_id"] for row in classified
        if row["screening_classification"] == "DIRECT_COMPETITOR"
    ]

    payload = {
        "status": "ABSTRACT_SCREENING_COMPLETE_FULL_TEXT_VERIFICATION_PENDING",
        "screened_unique_records": len(classified),
        "classification_counts": dict(sorted(counts.items())),
        "direct_competitor_count_at_abstract_screening": len(direct_competitors),
        "direct_competitor_screening_ids": direct_competitors,
        "full_text_verification_subset_count": len(full_text_subset),
        "full_text_verification_subset": full_text_subset,
        "provisional_interpretation": (
            "No direct competitor was identified at abstract screening. "
            "A locus-specific NFKB1 human population-functional study is the closest "
            "near competitor; several TNFAIP3 and CNV/pathway studies define partial "
            "overlap or supporting boundaries. Final novelty assessment remains pending "
            "full-text verification and supplementary citation-chaining/manual checks."
        ),
    }
    args.summary.parent.mkdir(parents=True, exist_ok=True)
    args.summary.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
