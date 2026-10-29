#!/usr/bin/env python3
"""Extract scoped sentence-level evidence from retrieved competitor full-text XML records."""

from __future__ import annotations

import argparse
import csv
import html
import json
import re
import xml.etree.ElementTree as ET
from pathlib import Path


TERMS = [
    "tnfaip3", "nfkb1", "nf-kb", "nf-κb", "nuclear factor",
    "introgress", "population", "selection", "eqtl", "haplotype",
    "copy number", "cnv", "pathway", "circuit",
]

MAX_EXCERPTS_PER_RECORD = 14


def xml_to_text(path: Path) -> str:
    root = ET.fromstring(path.read_bytes())
    parts = []
    for node in root.iter():
        if node.text:
            parts.append(node.text)
        if node.tail:
            parts.append(node.tail)
    text = html.unescape(" ".join(parts))
    text = re.sub(r"\s+", " ", text).strip()
    return text


def sentences(text: str) -> list[str]:
    return re.split(r"(?<=[.!?])\s+(?=[A-Z0-9])", text)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--retrieval-table", required=True, type=Path)
    parser.add_argument("--out-report", required=True, type=Path)
    parser.add_argument("--out-summary", required=True, type=Path)
    args = parser.parse_args()

    with args.retrieval_table.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))

    summaries = []

    args.out_report.parent.mkdir(parents=True, exist_ok=True)
    with args.out_report.open("w", encoding="utf-8") as report:
        report.write("# Project 6 full-text scope-verification evidence report\n\n")
        report.write(
            "Purpose: retrieve sentence-level passages relevant to determining whether "
            "candidate publications duplicate or only partially overlap the Project 6 "
            "NF-kB regulatory-circuit evolutionary-genomics design.\n\n"
        )

        for row in rows:
            if row["xml_retrieval_status"] != "FULLTEXT_XML_RETRIEVED":
                continue

            text = xml_to_text(Path(row["local_raw_xml"]))
            candidate_sentences = []
            seen = set()

            for sentence in sentences(text):
                lower = sentence.lower()
                matched = [term for term in TERMS if term in lower]
                if matched and sentence not in seen:
                    seen.add(sentence)
                    candidate_sentences.append((matched, sentence))

            selected = candidate_sentences[:MAX_EXCERPTS_PER_RECORD]
            summaries.append({
                "screening_id": row["screening_id"],
                "pmcid": row["pmcid"],
                "candidate_sentence_count": len(candidate_sentences),
                "excerpted_sentence_count": len(selected),
            })

            report.write("=" * 100 + "\n")
            report.write(f'Screening ID: {row["screening_id"]}\n')
            report.write(f'Classification before full-text verification: {row["classification"]}\n')
            report.write(f'DOI: {row["doi"]}\n')
            report.write(f'PMCID: {row["pmcid"]}\n')
            report.write(f'Candidate matched sentences: {len(candidate_sentences)}\n\n')

            for i, (matched, sentence) in enumerate(selected, start=1):
                report.write(f'Excerpt {i} | matched terms: {", ".join(matched)}\n')
                report.write(sentence + "\n\n")

    payload = {
        "status": "SCOPED_FULLTEXT_EVIDENCE_EXTRACTION_COMPLETE",
        "records_processed": len(summaries),
        "records": summaries,
        "interpretation_boundary": (
            "Keyword-matched excerpts support manual scope verification; they do not "
            "replace full-text scientific judgement or independently prove novelty."
        ),
    }
    args.out_summary.parent.mkdir(parents=True, exist_ok=True)
    args.out_summary.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
