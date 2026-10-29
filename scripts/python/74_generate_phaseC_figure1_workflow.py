#!/usr/bin/env python3
"""Generate Phase C Figure 1 workflow diagram.

This checkpoint generates a manuscript-facing workflow/audit-control figure
from committed Phase C manuscript source data. It does not perform new analysis,
write Results prose, validate function, infer selection, assign ancestry effects
or authorise broad biological interpretation.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

SRC = Path("manuscript/phaseC/figures/source_data/figure1_phaseC_workflow_gate_source_v1.tsv")

SVG_OUT = Path("manuscript/phaseC/figures/svg/figure1_phaseC_genotype_only_qc_workflow_v1.svg")
PNG_OUT = Path("manuscript/phaseC/figures/png/figure1_phaseC_genotype_only_qc_workflow_v1.png")
CAPTION_OUT = Path("manuscript/phaseC/captions/figure1_caption_v1.md")
SUMMARY_OUT = Path("manuscript/phaseC/audit/phaseC_figure1_generation_summary_v1.tsv")
OUTJSON = Path("results/qc/phaseC_figure1_generation_v1.json")

for path in [SVG_OUT.parent, PNG_OUT.parent, CAPTION_OUT.parent, SUMMARY_OUT.parent, OUTJSON.parent]:
    path.mkdir(parents=True, exist_ok=True)

failures = []

if not SRC.exists():
    failures.append(f"Missing Figure 1 source data: {SRC}")

rows = []
if SRC.exists():
    with SRC.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))

expected_gate_rows = 4
if len(rows) != expected_gate_rows:
    failures.append(f"Expected {expected_gate_rows} manuscript-readiness gate rows; observed {len(rows)}.")

steps = [
    ("1", "Strict genotype-only\nQC", "Allele-count and\ngroup-frequency QC"),
    ("2", "Candidate\nprioritisation", "Predeclared technical\ntriage rules"),
    ("3", "VEP annotation\nand evidence review", "TRAF3, BIRC3 and\npaired TNFAIP3 candidates"),
    ("4", "Manuscript\nreadiness gate", "Allowed statements,\nblocked claims and caveats"),
    ("5", "Results package\nand display blueprint", "Traceable files,\ntables and figures"),
]

fig, ax = plt.subplots(figsize=(13.5, 5.2))
ax.set_xlim(0, 13.5)
ax.set_ylim(0, 5.2)
ax.axis("off")

x_positions = [0.6, 3.2, 5.8, 8.4, 11.0]
box_w = 1.9
box_h = 1.45
y = 2.35

for i, (number, title, subtitle) in enumerate(steps):
    x = x_positions[i]
    box = FancyBboxPatch(
        (x, y),
        box_w,
        box_h,
        boxstyle="round,pad=0.08,rounding_size=0.08",
        linewidth=1.4,
        facecolor="white",
        edgecolor="black",
    )
    ax.add_patch(box)
    ax.text(x + 0.18, y + box_h - 0.28, number, fontsize=11, fontweight="bold", va="top")
    ax.text(x + box_w / 2, y + 0.88, title, fontsize=10.5, fontweight="bold", ha="center", va="center")
    ax.text(x + box_w / 2, y + 0.34, subtitle, fontsize=8.6, ha="center", va="center")

    if i < len(steps) - 1:
        arrow = FancyArrowPatch(
            (x + box_w + 0.08, y + box_h / 2),
            (x_positions[i + 1] - 0.08, y + box_h / 2),
            arrowstyle="-|>",
            mutation_scale=14,
            linewidth=1.3,
            color="black",
        )
        ax.add_patch(arrow)

# Boundary bar
boundary = FancyBboxPatch(
    (1.0, 0.55),
    11.5,
    0.85,
    boxstyle="round,pad=0.06,rounding_size=0.05",
    linewidth=1.2,
    facecolor="white",
    edgecolor="black",
)
ax.add_patch(boundary)
ax.text(
    6.75,
    1.08,
    "Interpretation boundary: descriptive genotype-only workflow; no causal, adaptive, ancestry-effect, population-specific or independent functional-validation claims",
    fontsize=9.0,
    ha="center",
    va="center",
)

ax.text(
    6.75,
    4.65,
    "Phase C genotype-only QC and manuscript-readiness workflow",
    fontsize=14,
    fontweight="bold",
    ha="center",
    va="center",
)

ax.text(
    6.75,
    4.25,
    "All downstream manuscript displays must preserve the committed readiness gate",
    fontsize=10,
    ha="center",
    va="center",
)

fig.savefig(SVG_OUT, format="svg", bbox_inches="tight")
fig.savefig(PNG_OUT, format="png", dpi=1800, bbox_inches="tight")
plt.close(fig)

CAPTION_OUT.write_text(
    "Figure 1. Phase C genotype-only QC and manuscript-readiness workflow. "
    "The workflow summarises the controlled progression from strict genotype-only QC through candidate prioritisation, "
    "annotation-supported evidence review, manuscript-readiness gating, results-package indexing and display blueprinting. "
    "The figure is an audit-control workflow diagram only and does not imply causality, selection, adaptation, "
    "population specificity, ancestry effects or independent functional validation.\n",
    encoding="utf-8",
)

summary_rows = [
    {"metric": "figure_id", "value": "FIGURE_1"},
    {"metric": "source_rows", "value": len(rows)},
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
    "operation": "Generate Phase C Figure 1 workflow diagram.",
    "figure_id": "FIGURE_1",
    "source_rows": len(rows),
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
        "This checkpoint generates Figure 1 as a workflow/audit-control display only. It does not perform new analysis, "
        "write Results prose, validate function, infer selection, make population-differentiation claims, assign ancestry "
        "effects or authorise broad biological interpretation."
    ),
}

OUTJSON.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
print(json.dumps(payload, indent=2))

if payload["status"] != "PASS":
    raise SystemExit(1)
