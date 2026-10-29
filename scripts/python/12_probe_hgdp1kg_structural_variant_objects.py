#!/usr/bin/env python3
"""Probe public gnomAD HGDP+1kGP object listings for structural-variant resources without downloading datasets."""

from __future__ import annotations

import argparse
import csv
import json
import re
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from pathlib import Path


BUCKET_ENDPOINT = "https://storage.googleapis.com/gcp-public-data--gnomad"

PREFIXES = {
    "secondary_analyses_hgdp_1kg_v2": "release/3.1/secondary_analyses/hgdp_1kg_v2/",
    "resources_hgdp_1kg": "resources/hgdp_1kg/",
    "v312_hgdp_tgp_vcf_objects": "release/3.1.2/vcf/genomes/gnomad.genomes.v3.1.2.hgdp_tgp",
}

SV_TERMS = re.compile(
    r"(^|[._/\-])(sv|svs|structural|cnv|copy[_-]?number|del|dup|inv|bnd|"
    r"breakend|mei|mobile[_-]?element)([._/\-]|$)",
    re.IGNORECASE,
)


def list_objects(prefix: str) -> list[dict[str, str | int]]:
    objects: list[dict[str, str | int]] = []
    continuation_token: str | None = None

    while True:
        params = {
            "list-type": "2",
            "prefix": prefix,
            "max-keys": "1000",
        }
        if continuation_token:
            params["continuation-token"] = continuation_token

        url = f"{BUCKET_ENDPOINT}?{urllib.parse.urlencode(params)}"
        with urllib.request.urlopen(url, timeout=120) as response:
            xml = response.read()

        root = ET.fromstring(xml)
        namespace = {"s3": "http://doc.s3.amazonaws.com/2006-03-01"}

        for content in root.findall("s3:Contents", namespace):
            key = content.findtext("s3:Key", default="", namespaces=namespace)
            size = int(content.findtext("s3:Size", default="0", namespaces=namespace))
            objects.append({"key": key, "size_bytes": size})

        truncated = root.findtext("s3:IsTruncated", default="false", namespaces=namespace)
        if truncated.lower() != "true":
            break

        continuation_token = root.findtext(
            "s3:NextContinuationToken", default="", namespaces=namespace
        )
        if not continuation_token:
            raise RuntimeError(f"Listing for prefix {prefix} was truncated but no token returned.")

    return objects


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out-all", required=True, type=Path)
    parser.add_argument("--out-sv-candidates", required=True, type=Path)
    parser.add_argument("--out-json", required=True, type=Path)
    args = parser.parse_args()

    all_rows: list[dict[str, str | int]] = []
    candidate_rows: list[dict[str, str | int]] = []

    for label, prefix in PREFIXES.items():
        objects = list_objects(prefix)
        for obj in objects:
            row = {
                "prefix_label": label,
                "queried_prefix": prefix,
                "object_key": obj["key"],
                "size_bytes": obj["size_bytes"],
                "sv_term_match": "yes" if SV_TERMS.search(str(obj["key"])) else "no",
            }
            all_rows.append(row)
            if row["sv_term_match"] == "yes":
                candidate_rows.append(row)

    args.out_all.parent.mkdir(parents=True, exist_ok=True)
    args.out_sv_candidates.parent.mkdir(parents=True, exist_ok=True)
    args.out_json.parent.mkdir(parents=True, exist_ok=True)

    fields = [
        "prefix_label",
        "queried_prefix",
        "object_key",
        "size_bytes",
        "sv_term_match",
    ]

    for output_path, rows in [
        (args.out_all, all_rows),
        (args.out_sv_candidates, candidate_rows),
    ]:
        with output_path.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(
                handle, fieldnames=fields, delimiter="\t", lineterminator="\n"
            )
            writer.writeheader()
            writer.writerows(rows)

    counts_by_prefix = {
        label: sum(1 for row in all_rows if row["prefix_label"] == label)
        for label in PREFIXES
    }
    candidates_by_prefix = {
        label: sum(1 for row in candidate_rows if row["prefix_label"] == label)
        for label in PREFIXES
    }

    payload = {
        "status": "RECONNAISSANCE_COMPLETE",
        "method": "Google Cloud public-bucket XML object listing only; no dataset object downloaded.",
        "queried_prefixes": PREFIXES,
        "object_counts_by_prefix": counts_by_prefix,
        "sv_candidate_counts_by_prefix": candidates_by_prefix,
        "sv_candidate_object_count": len(candidate_rows),
        "decision_rule": (
            "Any candidate SV object must still be verified for GRCh38 alignment, "
            "HGDP+1kGP sample linkage, genotype-level content and compatibility "
            "with the approved analytical cohorts before inclusion."
        ),
    }

    args.out_json.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
