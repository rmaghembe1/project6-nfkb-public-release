#!/usr/bin/env python3
"""Export Phase C Markdown manuscript assembly to DOCX.

This checkpoint uses pandoc to export the audited Phase C Markdown manuscript
assembly to DOCX. It does not perform new analysis, validate function, infer
selection, assign ancestry effects or authorise broad biological interpretation.
"""

from __future__ import annotations

import csv
import json
import subprocess
import zipfile
from pathlib import Path

SOURCE_MD = Path("manuscript/phaseC/assembly/phaseC_manuscript_assembly_master_v1.md")
OUT_DOCX = Path("manuscript/phaseC/export/docx/phaseC_manuscript_assembly_v1.docx")
AUDIT_DIR = Path("manuscript/phaseC/export/audit")
OUTJSON = Path("results/qc/phaseC_docx_export_v1.json")

for path in [OUT_DOCX.parent, AUDIT_DIR, OUTJSON.parent]:
    path.mkdir(parents=True, exist_ok=True)

failures = []

if not SOURCE_MD.exists():
    failures.append(f"Missing source Markdown: {SOURCE_MD}")

pandoc_path = subprocess.check_output(["bash", "-lc", "command -v pandoc || true"], text=True).strip()
if not pandoc_path:
    failures.append("pandoc not found on PATH")

if not failures:
    cmd = [
        "pandoc",
        str(SOURCE_MD),
        "-o",
        str(OUT_DOCX),
        "--standalone",
        "--toc",
        "--metadata",
        "title=Population-resolved NF-kB regulatory-circuitry variation: a genotype-only, evidence-gated Phase C manuscript assembly",
    ]
    proc = subprocess.run(cmd, text=True, capture_output=True)
    if proc.returncode != 0:
        failures.append(f"pandoc export failed: {proc.stderr.strip()}")

docx_exists = OUT_DOCX.exists()
docx_size = OUT_DOCX.stat().st_size if docx_exists else 0
if not docx_exists:
    failures.append(f"DOCX was not created: {OUT_DOCX}")
elif docx_size == 0:
    failures.append(f"DOCX is empty: {OUT_DOCX}")

docx_internal_checks = {
    "word/document.xml": "NO",
    "[Content_Types].xml": "NO",
}
docx_text = ""
if docx_exists and docx_size > 0:
    try:
        with zipfile.ZipFile(OUT_DOCX) as zf:
            names = set(zf.namelist())
            for item in docx_internal_checks:
                docx_internal_checks[item] = "YES" if item in names else "NO"
            if "word/document.xml" in names:
                docx_text = zf.read("word/document.xml").decode("utf-8", errors="replace").lower()
    except Exception as exc:
        failures.append(f"DOCX zip/read check failed: {exc}")

for item, status in docx_internal_checks.items():
    if status != "YES":
        failures.append(f"Missing DOCX internal component: {item}")

unsupported_phrases = [
    "demonstrates selection",
    "proves adaptation",
    "establishes ancestry effects",
    "causal variant",
    "population-specific variant",
    "functionally validated by phase c",
]
unsupported_hits = [phrase for phrase in unsupported_phrases if phrase in docx_text]
if unsupported_hits:
    failures.append(f"Unsupported phrase hits in DOCX XML: {unsupported_hits}")

manifest_path = AUDIT_DIR / "phaseC_docx_export_manifest_v1.tsv"
with manifest_path.open("w", newline="", encoding="utf-8") as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=["file", "exists", "size_bytes", "role"],
        delimiter="\t",
        lineterminator="\n",
    )
    writer.writeheader()
    writer.writerow({
        "file": str(SOURCE_MD),
        "exists": "YES" if SOURCE_MD.exists() else "NO",
        "size_bytes": SOURCE_MD.stat().st_size if SOURCE_MD.exists() else 0,
        "role": "source_markdown",
    })
    writer.writerow({
        "file": str(OUT_DOCX),
        "exists": "YES" if docx_exists else "NO",
        "size_bytes": docx_size,
        "role": "docx_export",
    })

summary_path = AUDIT_DIR / "phaseC_docx_export_summary_v1.tsv"
summary_rows = [
    {"metric": "source_markdown_exists", "value": "YES" if SOURCE_MD.exists() else "NO"},
    {"metric": "pandoc_found", "value": "YES" if pandoc_path else "NO"},
    {"metric": "docx_generated", "value": "YES" if docx_exists and docx_size > 0 else "NO"},
    {"metric": "docx_size_bytes", "value": docx_size},
    {"metric": "document_xml_present", "value": docx_internal_checks["word/document.xml"]},
    {"metric": "content_types_present", "value": docx_internal_checks["[Content_Types].xml"]},
    {"metric": "unsupported_phrase_hits", "value": len(unsupported_hits)},
    {"metric": "new_analysis_performed", "value": "NO"},
    {"metric": "new_biological_claims_added", "value": "NO"},
    {"metric": "selection_or_population_claim_authorised", "value": "NO"},
    {"metric": "broad_biological_interpretation_authorised", "value": "NO"},
]
with summary_path.open("w", newline="", encoding="utf-8") as handle:
    writer = csv.DictWriter(handle, fieldnames=["metric", "value"], delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(summary_rows)

payload = {
    "status": "PASS" if not failures else "FAIL",
    "operation": "Export Phase C Markdown manuscript assembly to DOCX.",
    "source_markdown": str(SOURCE_MD),
    "docx": str(OUT_DOCX),
    "pandoc_found": "YES" if pandoc_path else "NO",
    "docx_generated": "YES" if docx_exists and docx_size > 0 else "NO",
    "docx_size_bytes": docx_size,
    "unsupported_phrase_hits": len(unsupported_hits),
    "new_analysis_performed": "NO",
    "new_biological_claims_added": "NO",
    "selection_or_population_claim_authorised": "NO",
    "broad_biological_interpretation_authorised": "NO",
    "outputs": {
        "manifest": str(manifest_path),
        "summary": str(summary_path),
    },
    "failures": failures,
    "boundary": (
        "This checkpoint exports DOCX from the audited Markdown manuscript assembly only. "
        "It does not perform new analysis, validate function, infer selection, assign ancestry effects "
        "or authorise broad biological interpretation."
    ),
}

OUTJSON.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
print(json.dumps(payload, indent=2))

if payload["status"] != "PASS":
    raise SystemExit(1)
