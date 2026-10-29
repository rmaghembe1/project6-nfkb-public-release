#!/usr/bin/env python3
"""
Phase D local LD-proxy execution script.

STATUS AT COORDINATE-ANCHOR PATCH-REVIEW CHECKPOINT: PATCHED_NOT_RUN.

This script is designed for a future execution checkpoint only. It must not be
interpreted as having generated LD-proxy evidence unless it is explicitly run in
a later committed execution checkpoint.

Patch rationale:
- The rsID-based local route-status execution did not find candidate rsIDs in
  the phased BCFs.
- The coordinate/ref/alt audit found exact local phased BCF records for BIRC3
  rs924735278 and both TNFAIP3 candidates, but with ID='.'.
- TRAF3 rs1131877 was not found by exact coordinate/ref/alt in the local phased
  BCF and is therefore reported as a local-route status, not an absence claim.

Evidence boundaries:
- LD-proxy output is not direct candidate evidence.
- Empty output, candidate-not-found, coordinate-not-found, tool failure,
  parse error or no proxy returned must not be interpreted as absence of LD evidence.
- Population labels are technical reference LD contexts only, not ancestry,
  selection, adaptation or population-specific biological interpretation.
"""

from __future__ import annotations

import csv
import json
import shutil
import subprocess
from datetime import datetime, timezone
from pathlib import Path

CANDIDATES = [
    {
        "query_variant": "rs1131877",
        "anchor_id": "chr14:102875712:T:C",
        "query_gene": "TRAF3",
        "query_coordinate_GRCh38": "chr14:102875712:T:C",
        "chromosome": "14",
        "chromosome_label": "chr14",
        "position": 102875712,
        "ref": "T",
        "alt": "C",
        "source_file": "03_raw_genomic_resources/interval_extraction/full_panel/phased/TRAF3_GRCh38_hgdp1kgp_phased.bcf",
        "local_coordinate_presence_expected": "NO",
        "BIRC3_identity_caution": "NO",
        "TNFAIP3_paired_haplotype_context": "NO",
    },
    {
        "query_variant": "rs924735278",
        "anchor_id": "chr11:102324736:G:T",
        "query_gene": "BIRC3",
        "query_coordinate_GRCh38": "chr11:102324736:G:T",
        "chromosome": "11",
        "chromosome_label": "chr11",
        "position": 102324736,
        "ref": "G",
        "alt": "T",
        "source_file": "03_raw_genomic_resources/interval_extraction/full_panel/phased/BIRC3_GRCh38_hgdp1kgp_phased.bcf",
        "local_coordinate_presence_expected": "YES",
        "BIRC3_identity_caution": "YES",
        "TNFAIP3_paired_haplotype_context": "NO",
    },
    {
        "query_variant": "rs376205580",
        "anchor_id": "chr6:137874871:A:G",
        "query_gene": "TNFAIP3",
        "query_coordinate_GRCh38": "chr6:137874871:A:G",
        "chromosome": "6",
        "chromosome_label": "chr6",
        "position": 137874871,
        "ref": "A",
        "alt": "G",
        "source_file": "03_raw_genomic_resources/interval_extraction/full_panel/phased/TNFAIP3_GRCh38_hgdp1kgp_phased.bcf",
        "local_coordinate_presence_expected": "YES",
        "BIRC3_identity_caution": "NO",
        "TNFAIP3_paired_haplotype_context": "YES",
    },
    {
        "query_variant": "rs141807543",
        "anchor_id": "chr6:137875820:A:C",
        "query_gene": "TNFAIP3",
        "query_coordinate_GRCh38": "chr6:137875820:A:C",
        "chromosome": "6",
        "chromosome_label": "chr6",
        "position": 137875820,
        "ref": "A",
        "alt": "C",
        "source_file": "03_raw_genomic_resources/interval_extraction/full_panel/phased/TNFAIP3_GRCh38_hgdp1kgp_phased.bcf",
        "local_coordinate_presence_expected": "YES",
        "BIRC3_identity_caution": "NO",
        "TNFAIP3_paired_haplotype_context": "YES",
    },
]

OUTDIR = Path("11_results/phaseD/regulatory_functional_evidence/query_pass1/local_ld_proxy_coordinate_anchor_execution")
TMPDIR = OUTDIR / "tmp"
AUDIT_JSON = Path("results/qc/phaseD_local_ld_proxy_coordinate_anchor_execution_v1.json")

WINDOW_BP = 1_000_000
HIGH_LD_R2 = 0.8
GENOME_BUILD = "GRCh38"
POPULATION_PANEL = "HGDP_1KGP_all_available_samples_technical_reference"
SOURCE_ROUTE = "local_HGDP_1KGP_phased_reference_route_coordinate_anchor"

REQUIRED_TOOLS = ["bcftools", "plink2"]


def run_cmd(cmd: list[str], log_path: Path) -> subprocess.CompletedProcess[str]:
    with log_path.open("w", encoding="utf-8") as log:
        log.write("COMMAND\t" + " ".join(cmd) + "\n\n")
        proc = subprocess.run(cmd, text=True, capture_output=True)
        log.write("STDOUT\n")
        log.write(proc.stdout)
        log.write("\nSTDERR\n")
        log.write(proc.stderr)
        log.write(f"\nRETURN_CODE\t{proc.returncode}\n")
        return proc


def candidate_region(candidate: dict[str, object], window: int = WINDOW_BP) -> str:
    pos = int(candidate["position"])
    start = max(1, pos - window)
    end = pos + window
    return f"{candidate['chromosome_label']}:{start}-{end}"


def candidate_point(candidate: dict[str, object]) -> str:
    return f"{candidate['chromosome_label']}:{candidate['position']}-{candidate['position']}"


def verify_tools() -> list[str]:
    missing = []
    for tool in REQUIRED_TOOLS:
        if shutil.which(tool) is None:
            missing.append(tool)
    return missing


def write_status_row(candidate: dict[str, object], status: str, error_message: str = "") -> dict[str, str]:
    return {
        "query_variant": str(candidate["query_variant"]),
        "anchor_id": str(candidate["anchor_id"]),
        "query_gene": str(candidate["query_gene"]),
        "query_coordinate_GRCh38": str(candidate["query_coordinate_GRCh38"]),
        "proxy_variant": "",
        "proxy_coordinate_GRCh38": "",
        "chromosome": str(candidate["chromosome_label"]),
        "position": str(candidate["position"]),
        "alleles": f"{candidate['ref']}/{candidate['alt']}",
        "distance_bp": "",
        "D_prime": "",
        "r2": "",
        "population_panel": POPULATION_PANEL,
        "genome_build": GENOME_BUILD,
        "source_route": SOURCE_ROUTE,
        "source_file": str(candidate["source_file"]),
        "query_timestamp": datetime.now(timezone.utc).isoformat(),
        "query_status": status,
        "error_message": error_message[:500],
        "evidence_boundary": (
            "LD_proxy_status_only_not_direct_candidate_evidence_not_functional_validation_"
            "not_absence_of_LD_evidence"
        ),
    }


def record_matches_candidate(line: str, candidate: dict[str, object]) -> bool:
    fields = line.split("\t")
    if len(fields) < 5:
        return False
    chrom, pos, _vid, ref, alt = fields[:5]
    return (
        chrom == str(candidate["chromosome_label"])
        and pos == str(candidate["position"])
        and ref == str(candidate["ref"])
        and str(candidate["alt"]) in alt.split(",")
    )


def expected_coordinate_id(candidate: dict[str, object]) -> str:
    return f"{candidate['chromosome_label']}:{candidate['position']}:{candidate['ref']}:{candidate['alt']}"


def coordinate_id_setid_expression() -> str:
    return "%CHROM:%POS:%REF:%ALT"


def parse_ld_files(candidate: dict[str, object], ld_files: list[Path]) -> tuple[list[dict[str, str]], list[dict[str, str]], bool]:
    all_rows: list[dict[str, str]] = []
    high_rows: list[dict[str, str]] = []
    parsed_any = False
    query = str(candidate["query_variant"])
    anchor_id = str(candidate["anchor_id"])

    for ld_file in ld_files:
        with ld_file.open(encoding="utf-8", errors="replace") as handle:
            header = handle.readline().strip().split()
            if not header:
                continue

            for line in handle:
                parts = line.strip().split()
                if len(parts) != len(header):
                    continue
                rec = dict(zip(header, parts))

                r2_value = rec.get("R2") or rec.get("UNPHASED_R2") or rec.get("PHASED_R2") or ""
                try:
                    r2_float = float(r2_value)
                except Exception:
                    continue

                proxy_id = (
                    rec.get("ID_B")
                    or rec.get("ID")
                    or rec.get("SNP_B")
                    or rec.get("ID2")
                    or ""
                )

                if proxy_id == anchor_id:
                    continue

                proxy_chrom = rec.get("CHROM_B") or rec.get("#CHROM") or str(candidate["chromosome_label"])
                proxy_pos = rec.get("POS_B") or rec.get("POS") or ""
                proxy_ref = rec.get("REF_B") or rec.get("REF") or ""
                proxy_alt = rec.get("ALT_B") or rec.get("ALT") or ""

                distance = ""
                try:
                    if proxy_pos:
                        distance = str(abs(int(proxy_pos) - int(candidate["position"])))
                except Exception:
                    distance = ""

                row = {
                    "query_variant": query,
                    "anchor_id": anchor_id,
                    "query_gene": str(candidate["query_gene"]),
                    "query_coordinate_GRCh38": str(candidate["query_coordinate_GRCh38"]),
                    "proxy_variant": proxy_id,
                    "proxy_coordinate_GRCh38": f"{proxy_chrom}:{proxy_pos}" if proxy_pos else "",
                    "chromosome": proxy_chrom,
                    "position": proxy_pos,
                    "alleles": f"{proxy_ref}/{proxy_alt}".strip("/"),
                    "distance_bp": distance,
                    "D_prime": "",
                    "r2": str(r2_float),
                    "population_panel": POPULATION_PANEL,
                    "genome_build": GENOME_BUILD,
                    "source_route": SOURCE_ROUTE,
                    "source_file": str(candidate["source_file"]),
                    "query_timestamp": datetime.now(timezone.utc).isoformat(),
                    "query_status": "success",
                    "error_message": "",
                    "evidence_boundary": (
                        "LD_proxy_evidence_only_not_direct_candidate_evidence_not_functional_validation"
                    ),
                }
                all_rows.append(row)
                parsed_any = True
                if r2_float >= HIGH_LD_R2:
                    high_rows.append(row)

    return all_rows, high_rows, parsed_any


def main() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    TMPDIR.mkdir(parents=True, exist_ok=True)
    AUDIT_JSON.parent.mkdir(parents=True, exist_ok=True)

    missing_tools = verify_tools()
    all_rows: list[dict[str, str]] = []
    high_ld_rows: list[dict[str, str]] = []
    command_logs: list[str] = []

    if missing_tools:
        for candidate in CANDIDATES:
            all_rows.append(write_status_row(candidate, "tool_error", "Missing tools: " + ";".join(missing_tools)))
    else:
        for candidate in CANDIDATES:
            query = str(candidate["query_variant"])
            source = Path(str(candidate["source_file"]))
            region = candidate_region(candidate)

            if not source.exists():
                all_rows.append(write_status_row(candidate, "tool_error", f"Missing source BCF: {source}"))
                continue

            # Coordinate/ref/alt presence check at exact candidate position.
            presence_log = OUTDIR / f"{query}_coordinate_anchor_presence_check.log"
            presence_cmd = ["bcftools", "view", "-H", str(source), candidate_point(candidate)]
            presence_proc = run_cmd(presence_cmd, presence_log)
            command_logs.append(str(presence_log))

            if presence_proc.returncode != 0:
                all_rows.append(write_status_row(candidate, "tool_error", f"bcftools coordinate check failed; see {presence_log}"))
                continue

            if not any(record_matches_candidate(line, candidate) for line in presence_proc.stdout.splitlines()):
                all_rows.append(
                    write_status_row(
                        candidate,
                        "exact_candidate_coordinate_not_found_in_local_phased_BCF",
                        "Exact GRCh38 coordinate/ref/alt candidate record not found in local phased BCF. This is not absence of LD evidence.",
                    )
                )
                continue

            # Extract candidate-centred window as VCF.gz.
            window_vcf = TMPDIR / f"{query}_plusminus1Mb.raw.vcf.gz"
            window_log = OUTDIR / f"{query}_window_extraction.log"
            extract_cmd = ["bcftools", "view", "-Oz", "-o", str(window_vcf), str(source), region]
            extract_proc = run_cmd(extract_cmd, window_log)
            command_logs.append(str(window_log))

            if extract_proc.returncode != 0 or not window_vcf.exists():
                all_rows.append(write_status_row(candidate, "window_extraction_failed", f"Window extraction failed; see {window_log}"))
                continue

            index_raw_log = OUTDIR / f"{query}_raw_window_index.log"
            index_raw_cmd = ["bcftools", "index", "-t", "-f", str(window_vcf)]
            index_raw_proc = run_cmd(index_raw_cmd, index_raw_log)
            command_logs.append(str(index_raw_log))
            if index_raw_proc.returncode != 0:
                all_rows.append(write_status_row(candidate, "tool_error", f"Raw window indexing failed; see {index_raw_log}"))
                continue

            anchored_vcf = TMPDIR / f"{query}_plusminus1Mb.coordinate_id.vcf.gz"
            annotate_log = OUTDIR / f"{query}_coordinate_id_setid.log"
            annotate_cmd = [
                "bcftools", "annotate",
                "--set-id", coordinate_id_setid_expression(),
                "-Oz", "-o", str(anchored_vcf),
                str(window_vcf),
            ]
            annotate_proc = run_cmd(annotate_cmd, annotate_log)
            command_logs.append(str(annotate_log))

            if annotate_proc.returncode != 0 or not anchored_vcf.exists():
                all_rows.append(write_status_row(candidate, "coordinate_id_setid_failed", f"Coordinate ID set-id failed; see {annotate_log}"))
                continue

            if str(candidate["anchor_id"]) != expected_coordinate_id(candidate):
                all_rows.append(write_status_row(candidate, "coordinate_id_anchor_mismatch", "Internal anchor ID does not match expected coordinate ID."))
                continue

            index_anchor_log = OUTDIR / f"{query}_anchored_window_index.log"
            index_anchor_cmd = ["bcftools", "index", "-t", "-f", str(anchored_vcf)]
            index_anchor_proc = run_cmd(index_anchor_cmd, index_anchor_log)
            command_logs.append(str(index_anchor_log))
            if index_anchor_proc.returncode != 0:
                all_rows.append(write_status_row(candidate, "tool_error", f"Anchored window indexing failed; see {index_anchor_log}"))
                continue

            out_prefix = TMPDIR / f"{query}_plink2_coordinate_anchor_ld"
            ld_log = OUTDIR / f"{query}_plink2_coordinate_anchor_ld.log"
            ld_cmd = [
                "plink2",
                "--vcf", str(anchored_vcf),
                "--allow-extra-chr",
                "--r2-phased", "cols=+chrom,+pos,+ref,+alt",
                "--ld-snp", str(candidate["anchor_id"]),
                "--ld-window-kb", str(WINDOW_BP // 1000),
                "--ld-window", "999999",
                "--ld-window-r2", "0",
                "--out", str(out_prefix),
            ]
            ld_proc = run_cmd(ld_cmd, ld_log)
            command_logs.append(str(ld_log))

            if ld_proc.returncode != 0:
                all_rows.append(write_status_row(candidate, "tool_error", f"plink2 coordinate-anchor LD failed; see {ld_log}"))
                continue

            ld_files = sorted(TMPDIR.glob(f"{query}_plink2_coordinate_anchor_ld*.vcor*")) + sorted(
                TMPDIR.glob(f"{query}_plink2_coordinate_anchor_ld*.ld*")
            )
            if not ld_files:
                all_rows.append(write_status_row(candidate, "no_proxy_returned", "No PLINK2 LD output file detected."))
                continue

            parsed_rows, parsed_high, parsed_any = parse_ld_files(candidate, ld_files)
            all_rows.extend(parsed_rows)
            high_ld_rows.extend(parsed_high)

            if not parsed_any:
                all_rows.append(write_status_row(candidate, "parse_error", "PLINK2 output found but no LD rows parsed."))

    fields = [
        "query_variant",
        "anchor_id",
        "query_gene",
        "query_coordinate_GRCh38",
        "proxy_variant",
        "proxy_coordinate_GRCh38",
        "chromosome",
        "position",
        "alleles",
        "distance_bp",
        "D_prime",
        "r2",
        "population_panel",
        "genome_build",
        "source_route",
        "source_file",
        "query_timestamp",
        "query_status",
        "error_message",
        "evidence_boundary",
    ]

    all_path = OUTDIR / "phaseD_local_ld_proxy_coordinate_anchor_all_status_records_v1.tsv"
    high_path = OUTDIR / "phaseD_local_ld_proxy_coordinate_anchor_high_ld_r2_ge_0_8_records_v1.tsv"

    with all_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(all_rows)

    with high_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(high_ld_rows)

    payload = {
        "status": "COMPLETED_COORDINATE_ANCHOR_SCRIPT_EXECUTION" if not missing_tools else "TOOL_ERROR",
        "candidate_ld_query_performed": "YES",
        "coordinate_anchor_route": "YES",
        "candidate_proxy_evidence_generated": "YES_IF_SUCCESS_ROWS_PRESENT",
        "all_status_records": str(all_path),
        "high_ld_records": str(high_path),
        "high_ld_threshold": HIGH_LD_R2,
        "command_logs": command_logs,
        "TRAF3_local_coordinate_status_boundary": (
            "exact_candidate_coordinate_not_found_in_local_phased_BCF_is_not_absence_of_LD_evidence"
        ),
        "BIRC3_identity_caution_required": "YES",
        "TNFAIP3_paired_haplotype_context_required": "YES",
        "functional_validation_claim_allowed": "NO",
        "causality_claim_allowed": "NO",
        "selection_or_adaptation_claim_allowed": "NO",
        "ancestry_effect_claim_allowed": "NO",
        "population_specific_biology_claim_allowed": "NO",
        "tissue_specific_biological_claim_allowed": "NO",
        "absence_boundary": (
            "candidate_not_found_coordinate_not_found_tool_failure_empty_output_or_no_proxy_returned_is_not_absence_of_LD_evidence"
        ),
    }
    AUDIT_JSON.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
