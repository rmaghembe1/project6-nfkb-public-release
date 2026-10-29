#!/usr/bin/env python3
from __future__ import annotations

import csv
import hashlib
import json
import subprocess
from pathlib import Path

CANDIDATES = [
    ("TRAF3", "rs1131877", "chr14", 102875712, "T", "C",
     "03_raw_genomic_resources/interval_extraction/full_panel/phased/TRAF3_GRCh38_hgdp1kgp_phased.bcf"),
    ("BIRC3", "rs924735278", "chr11", 102324736, "G", "T",
     "03_raw_genomic_resources/interval_extraction/full_panel/phased/BIRC3_GRCh38_hgdp1kgp_phased.bcf"),
    ("TNFAIP3", "rs376205580", "chr6", 137874871, "A", "G",
     "03_raw_genomic_resources/interval_extraction/full_panel/phased/TNFAIP3_GRCh38_hgdp1kgp_phased.bcf"),
    ("TNFAIP3", "rs141807543", "chr6", 137875820, "A", "C",
     "03_raw_genomic_resources/interval_extraction/full_panel/phased/TNFAIP3_GRCh38_hgdp1kgp_phased.bcf"),
]

OUTDIR = Path("11_results/phaseD/regulatory_functional_evidence/query_pass1/local_ld_coordinate_presence_audit")
LOGDIR = Path("15_logs_and_audit/phaseD/query_pass1_local_ld_coordinate_presence_audit")
OUTDIR.mkdir(parents=True, exist_ok=True)
LOGDIR.mkdir(parents=True, exist_ok=True)

rows = []

for gene, rsid, chrom, pos, ref, alt, bcf in CANDIDATES:
    region = f"{chrom}:{pos}-{pos}"
    cmd = ["bcftools", "view", "-H", bcf, region]
    run = subprocess.run(cmd, text=True, capture_output=True)

    log_path = OUTDIR / f"{rsid}_coordinate_presence_bcftools.log"
    log_path.write_text(
        "COMMAND\t" + " ".join(cmd) + "\n\nSTDOUT\n" + run.stdout +
        "\nSTDERR\n" + run.stderr + f"\nRETURN_CODE\t{run.returncode}\n",
        encoding="utf-8"
    )

    matched_exact = "NO"
    record_count = 0
    ids_observed = []
    allele_strings = []

    for line in run.stdout.splitlines():
        fields = line.split("\t")
        if len(fields) < 5:
            continue
        record_count += 1
        obs_chrom, obs_pos, obs_id, obs_ref, obs_alt = fields[:5]
        ids_observed.append(obs_id)
        allele_strings.append(f"{obs_ref}>{obs_alt}")
        alt_list = obs_alt.split(",")
        if obs_chrom == chrom and obs_pos == str(pos) and obs_ref == ref and alt in alt_list:
            matched_exact = "YES"

    rows.append({
        "gene": gene,
        "rsid": rsid,
        "coordinate_GRCh38": f"{chrom}:{pos}:{ref}:{alt}",
        "source_file": bcf,
        "bcftools_return_code": str(run.returncode),
        "records_at_exact_position": str(record_count),
        "ids_observed": ";".join(sorted(set(ids_observed))) if ids_observed else "",
        "alleles_observed": ";".join(sorted(set(allele_strings))) if allele_strings else "",
        "exact_coordinate_ref_alt_match": matched_exact,
        "candidate_ld_extraction_performed": "NO",
        "proxy_evidence_generated": "NO",
        "absence_claim_authorised": "NO",
        "boundary": "coordinate_presence_audit_only_not_LD_proxy_evidence_not_functional_validation"
    })

out_tsv = OUTDIR / "phaseD_local_ld_candidate_coordinate_presence_audit_v1.tsv"
with out_tsv.open("w", newline="", encoding="utf-8") as handle:
    writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()), delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)

summary = {
    "status": "PASS",
    "checkpoint_type": "local_ld_candidate_coordinate_presence_audit",
    "candidate_count": len(rows),
    "candidate_ld_extraction_performed": "NO",
    "proxy_evidence_generated": "NO",
    "candidates_with_exact_coordinate_ref_alt_match": sum(1 for r in rows if r["exact_coordinate_ref_alt_match"] == "YES"),
    "absence_claim_authorised": "NO",
    "boundary": "This audit checks exact coordinate/ref/alt presence only. It does not extract LD, generate proxy evidence, validate function, infer causality, or authorise absence claims."
}

summary_json = Path("results/qc/phaseD_local_ld_candidate_coordinate_presence_audit_v1.json")
summary_json.parent.mkdir(parents=True, exist_ok=True)
summary_json.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")

decision = LOGDIR / "PhaseD_Local_LD_Candidate_Coordinate_Presence_Audit_Decision_2026-06-09.md"
decision.write_text(f"""# Phase D local LD candidate coordinate-presence audit decision: 2026-06-09

## Scope

This checkpoint audits whether the four Phase D candidate variants are present in local phased HGDP/1KGP BCF files by exact GRCh38 coordinate, reference allele and alternate allele.

No LD extraction is performed.

No candidate-centred LD window is extracted.

No proxy evidence is generated.

## Rationale

The previous local LD-proxy route-status checkpoint returned `candidate_not_found` for all candidates using rsID-based filtering. That result may reflect local BCF variant-ID representation and must not be interpreted as absence of LD proxies.

## Output state

Candidate count: {len(rows)}.

Candidates with exact coordinate/ref/alt match: {summary["candidates_with_exact_coordinate_ref_alt_match"]}.

Candidate LD extraction performed: NO.

Proxy evidence generated: NO.

Absence claim authorised: NO.

## Boundary

This checkpoint is a coordinate-presence audit only.

It does not demonstrate presence or absence of LD proxies.

It does not validate function, infer causality, infer selection or adaptation, infer ancestry effects, make population-specific biological claims, or make tissue-specific biological claims.

BIRC3 identity caution remains active.

TNFAIP3 paired/haplotype context remains active.
""", encoding="utf-8")

checksum = LOGDIR / "phaseD_local_ld_candidate_coordinate_presence_audit_outputs_sha256_2026-06-09.txt"
paths = [out_tsv, summary_json, decision] + sorted(OUTDIR.glob("*_coordinate_presence_bcftools.log"))
with checksum.open("w", encoding="utf-8") as handle:
    for path in paths:
        handle.write(f"{hashlib.sha256(path.read_bytes()).hexdigest()}  {path}\n")

print(json.dumps(summary, indent=2))
