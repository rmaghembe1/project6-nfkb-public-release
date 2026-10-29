#!/usr/bin/env python3
"""Generate Phase C Figure 3 annotation-supported evidence landscape.

This checkpoint generates a manuscript-facing evidence-landscape figure from
committed Phase C annotation-supported evidence source data. It does not perform
new analysis, write Results prose, validate function, infer selection, assign
ancestry effects or authorise broad biological interpretation.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

SRC = Path("manuscript/phaseC/figures/source_data/figure3_annotation_supported_evidence_landscape_source_v1.tsv")

SVG_OUT = Path("manuscript/phaseC/figures/svg/figure3_phaseC_annotation_supported_evidence_landscape_v1.svg")
PNG_OUT = Path("manuscript/phaseC/figures/png/figure3_phaseC_annotation_supported_evidence_landscape_v1.png")
CAPTION_OUT = Path("manuscript/phaseC/captions/figure3_caption_v1.md")
SUMMARY_OUT = Path("manuscript/phaseC/audit/phaseC_figure3_generation_summary_v1.tsv")
OUTJSON = Path("results/qc/phaseC_figure3_generation_v1.json")

for path in [SVG_OUT.parent, PNG_OUT.parent, CAPTION_OUT.parent, SUMMARY_OUT.parent, OUTJSON.parent]:
    path.mkdir(parents=True, exist_ok=True)

failures = []

if not SRC.exists():
    failures.append(f"Missing Figure 3 source data: {SRC}")
    rows = []
else:
    with SRC.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))

if len(rows) != 4:
    failures.append(f"Expected 4 annotation-supported candidate rows; observed {len(rows)}.")

def pick(row: dict[str, str], candidates: list[str], default: str = ".") -> str:
    lower_map = {k.lower(): k for k in row}
    for c in candidates:
        if c.lower() in lower_map:
            val = row.get(lower_map[c.lower()], "").strip()
            if val:
                return val
    for key, val in row.items():
        if any(c.lower() in key.lower() for c in candidates):
            if str(val).strip():
                return str(val).strip()
    return default

cards = []
for row in rows:
    candidate = pick(row, ["allele_key", "candidate", "variant", "locus"])
    gene = pick(row, ["gene", "gene_symbol"])
    mapped = pick(row, ["mapped_identity", "identity", "protein_change", "variant_identity"])
    grade = pick(row, ["overall_evidence_grade", "evidence_grade", "grade"])
    functional = pick(row, ["functional_validation_available", "functional_validation_status", "functional_validation"])
    selection = pick(row, ["selection_or_population_claim_supported", "selection_population_claim", "selection_claim"])
    boundary = pick(row, ["interpretation_boundary", "boundary", "manuscript_boundary"])
    cards.append(
        {
            "candidate": candidate,
            "gene": gene,
            "mapped": mapped,
            "grade": grade,
            "functional": functional,
            "selection": selection,
            "boundary": boundary,
        }
    )

# Stable display order by gene/module.
order = {"TRAF3": 0, "BIRC3": 1, "TNFAIP3": 2}
cards.sort(key=lambda x: (order.get(x["gene"], 99), x["candidate"]))

fig, ax = plt.subplots(figsize=(13.2, 7.2))
ax.set_xlim(0, 13.2)
ax.set_ylim(0, 7.2)
ax.axis("off")

ax.text(
    6.6,
    6.78,
    "Phase C annotation-supported evidence landscape",
    ha="center",
    va="center",
    fontsize=14,
    fontweight="bold",
)
ax.text(
    6.6,
    6.38,
    "Candidate-specific evidence review with explicit manuscript boundaries",
    ha="center",
    va="center",
    fontsize=10,
)

x_positions = [0.65, 3.8, 6.95, 10.1]
box_w = 2.45
box_h = 3.45
y = 2.35

for i, card in enumerate(cards):
    x = x_positions[i] if i < len(x_positions) else 0.65 + i * 3.15
    box = FancyBboxPatch(
        (x, y),
        box_w,
        box_h,
        boxstyle="round,pad=0.08,rounding_size=0.08",
        linewidth=1.35,
        facecolor="white",
        edgecolor="black",
    )
    ax.add_patch(box)

    ax.text(x + box_w / 2, y + 3.08, card["gene"], ha="center", va="center", fontsize=12.5, fontweight="bold")
    ax.text(x + box_w / 2, y + 2.68, card["candidate"], ha="center", va="center", fontsize=7.4)
    ax.text(x + box_w / 2, y + 2.23, f"Identity: {card['mapped']}", ha="center", va="center", fontsize=7.8, wrap=True)
    ax.text(x + box_w / 2, y + 1.76, f"Evidence: {card['grade']}", ha="center", va="center", fontsize=8.2, fontweight="bold")
    ax.text(x + box_w / 2, y + 1.28, f"Functional validation:\n{card['functional']}", ha="center", va="center", fontsize=7.4)
    ax.text(x + box_w / 2, y + 0.72, f"Selection/population claim:\n{card['selection']}", ha="center", va="center", fontsize=7.4)

boundary = FancyBboxPatch(
    (0.85, 0.55),
    11.5,
    1.15,
    boxstyle="round,pad=0.06,rounding_size=0.05",
    linewidth=1.2,
    facecolor="white",
    edgecolor="black",
)
ax.add_patch(boundary)
ax.text(
    6.6,
    1.16,
    "Boundary: exact-variant or local VEP support is reported only within candidate-specific evidence limits; no causal, adaptive, ancestry-effect, population-specific or independent single-variant functional-validation claims are authorised.",
    fontsize=8.8,
    ha="center",
    va="center",
    wrap=True,
)

fig.savefig(SVG_OUT, format="svg", bbox_inches="tight")
fig.savefig(PNG_OUT, format="png", dpi=1800, bbox_inches="tight")
plt.close(fig)

# Clean Matplotlib SVG trailing whitespace.
text = SVG_OUT.read_text(encoding="utf-8")
SVG_OUT.write_text("\n".join(line.rstrip(" \t") for line in text.splitlines()) + "\n", encoding="utf-8")

CAPTION_OUT.write_text(
    "Figure 3. Phase C annotation-supported evidence landscape. "
    "The figure summarises the first-pass evidence-review status for four annotation-supported candidates across TRAF3, BIRC3 and paired TNFAIP3. "
    "Candidate identity, evidence grade, functional-validation status and selection/population-claim status are displayed with explicit interpretation boundaries. "
    "The figure does not imply causality, selection, adaptation, population specificity, ancestry effects or independent single-variant functional validation.\n",
    encoding="utf-8",
)

summary_rows = [
    {"metric": "figure_id", "value": "FIGURE_3"},
    {"metric": "candidate_rows", "value": len(cards)},
    {"metric": "genes_displayed", "value": ",".join(sorted({c["gene"] for c in cards}))},
    {"metric": "svg_generated", "value": "YES" if SVG_OUT.exists() else "NO"},
    {"metric": "png_generated", "value": "YES" if PNG_OUT.exists() else "NO"},
    {"metric": "png_dpi", "value": "1800"},
    {"metric": "caption_written", "value": "YES" if CAPTION_OUT.exists() else "NO"},
    {"metric": "new_analysis_performed", "value": "NO"},
    {"metric": "results_prose_written", "value": "NO"},
    {"metric": "selection_or_population_claim_authorised", "value": "NO"},
    {"metric": "broad_biological_interpretation_authorised", "value": "NO"},
]
with SUMMARY_OUT.open("w", newline="", encoding="utf-8") as handle:
    writer = csv.DictWriter(handle, fieldnames=["metric", "value"], delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(summary_rows)

payload = {
    "status": "PASS" if not failures else "FAIL",
    "operation": "Generate Phase C Figure 3 annotation-supported evidence landscape.",
    "figure_id": "FIGURE_3",
    "candidate_rows": len(cards),
    "genes_displayed": sorted({c["gene"] for c in cards}),
    "svg_generated": "YES" if SVG_OUT.exists() else "NO",
    "png_generated": "YES" if PNG_OUT.exists() else "NO",
    "png_dpi": 1800,
    "caption_written": "YES" if CAPTION_OUT.exists() else "NO",
    "new_analysis_performed": "NO",
    "results_prose_written": "NO",
    "selection_or_population_claim_authorised": "NO",
    "broad_biological_interpretation_authorised": "NO",
    "outputs": {
        "svg": str(SVG_OUT),
        "png": str(PNG_OUT),
        "caption": str(CAPTION_OUT),
        "summary": str(SUMMARY_OUT),
    },
    "failures": failures,
    "boundary": (
        "This checkpoint generates Figure 3 as an annotation-supported evidence landscape only. "
        "It does not perform new analysis, write Results prose, validate function, infer selection, "
        "make population-differentiation claims, assign ancestry effects or authorise broad biological interpretation."
    ),
}

OUTJSON.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
print(json.dumps(payload, indent=2))

if payload["status"] != "PASS":
    raise SystemExit(1)
