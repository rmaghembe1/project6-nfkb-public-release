#!/usr/bin/env python3
"""Freeze Phase C gene-module context map for Project 6 NF-kB candidates.

This checkpoint creates a fixed gene-to-module context map before applying
candidate prioritisation scores. It does not score variants, rank candidates,
infer function, infer selection or make biological interpretation.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path

LANDSCAPE = Path(
    "11_results/phaseC/genotype_only_qc/vep_consequence_landscape/"
    "phaseC_vep_landscape_per_candidate_technical_summary_v1.tsv"
)
FRAMEWORK_JSON = Path("results/qc/phaseC_candidate_prioritisation_framework_design_v1.json")

OUTDIR = Path("11_results/phaseC/genotype_only_qc/gene_module_context_map")
OUTJSON = Path("results/qc/phaseC_gene_module_context_map_freeze_v1.json")

OUTDIR.mkdir(parents=True, exist_ok=True)
OUTJSON.parent.mkdir(parents=True, exist_ok=True)

failures: list[str] = []


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_rows(path: Path, fieldnames: list[str], rows: list[dict[str, object]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row.get(key, ".") for key in fieldnames})


if not LANDSCAPE.exists():
    failures.append(f"Missing VEP consequence landscape candidate summary: {LANDSCAPE}")
if not FRAMEWORK_JSON.exists():
    failures.append(f"Missing prioritisation framework design JSON: {FRAMEWORK_JSON}")

candidate_rows = read_tsv(LANDSCAPE) if LANDSCAPE.exists() else []
candidate_genes = sorted({row["gene"] for row in candidate_rows if row.get("gene") and row["gene"] != "."})

# Fixed Project 6 NF-kB circuit module map.
# This map is intentionally compact and predeclared before prioritisation is applied.
module_map = {
    "NFKB1": {
        "module": "NF-kB_transcription_factor_core",
        "submodule": "canonical_precursor_and_DNA_binding_subunit",
        "context_weight": 1,
        "context_note": "Core NF-kB transcription-factor family member in the predefined Project 6 panel.",
    },
    "NFKB2": {
        "module": "NF-kB_transcription_factor_core",
        "submodule": "noncanonical_precursor_and_DNA_binding_subunit",
        "context_weight": 1,
        "context_note": "Core NF-kB transcription-factor family member linked to noncanonical pathway context.",
    },
    "REL": {
        "module": "NF-kB_transcription_factor_core",
        "submodule": "Rel_family_transcriptional_subunit",
        "context_weight": 1,
        "context_note": "Core Rel-family transcriptional subunit in the predefined panel.",
    },
    "RELA": {
        "module": "NF-kB_transcription_factor_core",
        "submodule": "canonical_RelA_p65_transcriptional_subunit",
        "context_weight": 1,
        "context_note": "Core canonical NF-kB transcriptional subunit in the predefined panel.",
    },
    "RELB": {
        "module": "NF-kB_transcription_factor_core",
        "submodule": "noncanonical_RelB_transcriptional_subunit",
        "context_weight": 1,
        "context_note": "Core Rel-family subunit linked to noncanonical pathway context.",
    },
    "CHUK": {
        "module": "IKK_IkB_canonical_signalling_machinery",
        "submodule": "IKK_alpha",
        "context_weight": 1,
        "context_note": "IKK complex kinase component in canonical/noncanonical NF-kB regulation.",
    },
    "IKBKB": {
        "module": "IKK_IkB_canonical_signalling_machinery",
        "submodule": "IKK_beta",
        "context_weight": 1,
        "context_note": "IKK complex kinase component in canonical NF-kB activation.",
    },
    "IKBKG": {
        "module": "IKK_IkB_canonical_signalling_machinery",
        "submodule": "NEMO_regulatory_scaffold",
        "context_weight": 1,
        "context_note": "Regulatory scaffold of the IKK complex; chrX handled through sex-aware Phase C workflow.",
    },
    "NFKBIA": {
        "module": "IKK_IkB_canonical_signalling_machinery",
        "submodule": "IkB_inhibitor",
        "context_weight": 1,
        "context_note": "IkB inhibitory regulator in canonical NF-kB pathway.",
    },
    "NFKBIB": {
        "module": "IKK_IkB_canonical_signalling_machinery",
        "submodule": "IkB_inhibitor",
        "context_weight": 1,
        "context_note": "IkB inhibitory regulator in canonical NF-kB pathway.",
    },
    "NFKBIE": {
        "module": "IKK_IkB_canonical_signalling_machinery",
        "submodule": "IkB_inhibitor",
        "context_weight": 1,
        "context_note": "IkB inhibitory regulator in canonical NF-kB pathway.",
    },
    "MAP3K14": {
        "module": "noncanonical_NF-kB_branch",
        "submodule": "NIK_MAP3K14",
        "context_weight": 1,
        "context_note": "Central kinase in noncanonical NF-kB signalling.",
    },
    "BIRC2": {
        "module": "ubiquitin_editing_restraint_and_survival_regulation",
        "submodule": "cIAP1",
        "context_weight": 1,
        "context_note": "Ubiquitin-linked regulatory/survival component relevant to NF-kB pathway control.",
    },
    "BIRC3": {
        "module": "ubiquitin_editing_restraint_and_survival_regulation",
        "submodule": "cIAP2",
        "context_weight": 1,
        "context_note": "Ubiquitin-linked regulatory/survival component relevant to NF-kB pathway control.",
    },
    "TNFAIP3": {
        "module": "ubiquitin_editing_restraint_and_survival_regulation",
        "submodule": "A20_ubiquitin_editing_checkpoint",
        "context_weight": 1,
        "context_note": "A20/TNFAIP3 is the predefined Project 6 evolutionary anchor and negative-regulatory checkpoint.",
    },
    "TNIP1": {
        "module": "ubiquitin_editing_restraint_and_survival_regulation",
        "submodule": "A20_binding_regulatory_partner",
        "context_weight": 1,
        "context_note": "A20-interacting regulatory partner in NF-kB restraint context.",
    },
    "CYLD": {
        "module": "ubiquitin_editing_restraint_and_survival_regulation",
        "submodule": "deubiquitinase_checkpoint",
        "context_weight": 1,
        "context_note": "Deubiquitinase checkpoint component in NF-kB signal restraint.",
    },
    "OTULIN": {
        "module": "ubiquitin_editing_restraint_and_survival_regulation",
        "submodule": "linear_ubiquitin_deubiquitinase",
        "context_weight": 1,
        "context_note": "Linear ubiquitin editing component relevant to NF-kB regulation.",
    },
    "RBCK1": {
        "module": "ubiquitin_editing_restraint_and_survival_regulation",
        "submodule": "LUBAC_component_HOIL1",
        "context_weight": 1,
        "context_note": "LUBAC-associated component in linear ubiquitin/NF-kB regulatory context.",
    },
    "RNF31": {
        "module": "ubiquitin_editing_restraint_and_survival_regulation",
        "submodule": "LUBAC_component_HOIP",
        "context_weight": 1,
        "context_note": "LUBAC catalytic component in linear ubiquitin/NF-kB regulatory context.",
    },
    "SHARPIN": {
        "module": "ubiquitin_editing_restraint_and_survival_regulation",
        "submodule": "LUBAC_component_SHARPIN",
        "context_weight": 1,
        "context_note": "LUBAC-associated component in linear ubiquitin/NF-kB regulatory context.",
    },
    "TRAF2": {
        "module": "ubiquitin_editing_restraint_and_survival_regulation",
        "submodule": "TRAF_adaptor_ubiquitin_signalling",
        "context_weight": 1,
        "context_note": "TRAF adaptor in NF-kB-linked ubiquitin signalling context.",
    },
    "TRAF3": {
        "module": "ubiquitin_editing_restraint_and_survival_regulation",
        "submodule": "TRAF_adaptor_noncanonical_restraint",
        "context_weight": 1,
        "context_note": "TRAF adaptor with noncanonical NF-kB regulatory context.",
    },
    "BCL3": {
        "module": "NF-kB_transcriptional_co_regulation",
        "submodule": "atypical_IkB_transcriptional_cofactor",
        "context_weight": 1,
        "context_note": "Atypical IkB/co-regulatory factor in NF-kB transcriptional context.",
    },
}

module_rows = []
for gene in sorted(module_map):
    entry = module_map[gene]
    module_rows.append(
        {
            "gene": gene,
            "module": entry["module"],
            "submodule": entry["submodule"],
            "context_weight": entry["context_weight"],
            "context_note": entry["context_note"],
            "predeclared_before_scoring": "YES",
            "interpretation_boundary": "Context label only; not evidence of functional effect or selection.",
        }
    )

missing_from_map = sorted(set(candidate_genes) - set(module_map))
mapped_without_candidates = sorted(set(module_map) - set(candidate_genes))

if missing_from_map:
    failures.append(f"Candidate genes missing from module map: {','.join(missing_from_map)}")

module_path = OUTDIR / "phaseC_gene_module_context_map_v1.tsv"
write_rows(
    module_path,
    [
        "gene",
        "module",
        "submodule",
        "context_weight",
        "context_note",
        "predeclared_before_scoring",
        "interpretation_boundary",
    ],
    module_rows,
)

candidate_context_rows = []
for row in candidate_rows:
    gene = row["gene"]
    entry = module_map.get(gene, {})
    candidate_context_rows.append(
        {
            "allele_key": row["allele_key"],
            "gene": gene,
            "module": entry.get("module", "."),
            "submodule": entry.get("submodule", "."),
            "context_weight": entry.get("context_weight", "."),
            "context_available": "YES" if gene in module_map else "NO",
            "has_moderate_or_high_impact_canonical_like": row.get("has_moderate_or_high_impact_canonical_like", "."),
            "has_mane_or_appris_canonical_like_annotation": row.get("has_mane_or_appris_canonical_like_annotation", "."),
            "interpretation_boundary": "Context label only; not evidence of functional effect or selection.",
        }
    )

candidate_context_path = OUTDIR / "phaseC_candidate_gene_module_context_v1.tsv"
write_rows(
    candidate_context_path,
    [
        "allele_key",
        "gene",
        "module",
        "submodule",
        "context_weight",
        "context_available",
        "has_moderate_or_high_impact_canonical_like",
        "has_mane_or_appris_canonical_like_annotation",
        "interpretation_boundary",
    ],
    candidate_context_rows,
)

# Module-level candidate count summary.
module_counts: dict[tuple[str, str], int] = {}
for row in candidate_context_rows:
    key = (row["module"], row["has_moderate_or_high_impact_canonical_like"])
    module_counts[key] = module_counts.get(key, 0) + 1

module_count_rows = [
    {
        "module": module,
        "has_moderate_or_high_impact_canonical_like": flag,
        "candidate_count": count,
    }
    for (module, flag), count in sorted(module_counts.items())
]

module_count_path = OUTDIR / "phaseC_candidate_counts_by_module_and_canonical_impact_flag_v1.tsv"
write_rows(
    module_count_path,
    ["module", "has_moderate_or_high_impact_canonical_like", "candidate_count"],
    module_count_rows,
)

payload = {
    "status": "PASS" if not failures else "FAIL",
    "operation": "Freeze Project 6 gene-module context map before candidate prioritisation scoring.",
    "candidate_genes_detected": len(candidate_genes),
    "module_map_genes": len(module_map),
    "candidate_records_with_context": len(candidate_context_rows),
    "candidate_genes_missing_from_map": missing_from_map,
    "module_genes_without_current_candidates": mapped_without_candidates,
    "outputs": {
        "gene_module_context_map": str(module_path),
        "candidate_gene_module_context": str(candidate_context_path),
        "candidate_counts_by_module_and_canonical_impact_flag": str(module_count_path),
    },
    "failures": failures,
    "boundary": (
        "This checkpoint freezes gene-module context labels before prioritisation scoring. "
        "It does not apply scores, rank variants, select candidates, infer function, infer selection, "
        "make population-differentiation claims or authorise biological interpretation."
    ),
}

OUTJSON.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
print(json.dumps(payload, indent=2))

if payload["status"] != "PASS":
    raise SystemExit(1)
