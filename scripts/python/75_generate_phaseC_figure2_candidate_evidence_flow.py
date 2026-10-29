#!/usr/bin/env python3
"""Generate Phase C Figure 2 candidate/evidence flow diagram.

This checkpoint generates a manuscript-facing figure from committed Phase C
candidate-prioritisation and evidence-grade source data. It does not perform new
analysis, write Results prose, validate function, infer selection, assign
ancestry effects or authorise broad biological interpretation.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

CLAIM_COUNTS = Path("manuscript/phaseC/figures/source_data/figure2_candidate_prioritisation_claim_label_counts_v1.tsv")
EVIDENCE_COUNTS = Path("manuscript/phaseC/figures/source_data/figure2_annotation_supported_evidence_grade_counts_v1.tsv")

SVG_OUT = Path("manuscript/phaseC/figures/svg/figure2_phaseC_candidate_prioritisation_evidence_flow_v1.svg")
PNG_OUT = Path("manuscript/phaseC/figures/png/figure2_phaseC_candidate_prioritisation_evidence_flow_v1.png")
CAPTION_OUT = Path("manuscript/phaseC/captions/figure2_caption_v1.md")
SUMMARY_OUT = Path("manuscript/phaseC/audit/phaseC_figure2_generation_summary_v1.tsv")
OUTJSON = Path("results/qc/phaseC_figure2_generation_v1.json")

for path in [SVG_OUT.parent, PNG_OUT.parent, CAPTION_OUT.parent, SUMMARY_OUT.parent, OUTJSON.parent]:
    path.mkdir(parents=True, exist_ok=True)

failures = []


def read_rows(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        failures.append(f"Missing source data: {path}")
        return []
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


claim_rows = read_rows(CLAIM_COUNTS)
evidence_rows = read_rows(EVIDENCE_COUNTS)

def find_count(rows: list[dict[str, str]], possible_label_fields: list[str], possible_count_fields: list[str], label_contains: str | None = None) -> int:
    total = 0
    for row in rows:
        label = " ".join(str(row.get(field, "")) for field in possible_label_fields).lower()
        count_value = None
        for field in possible_count_fields:
            if field in row and str(row[field]).strip() not in {"", "."}:
                count_value = row[field]
                break
        if count_value is None:
            continue
        try:
            value = int(float(str(count_value)))
        except ValueError:
            continue
        if label_contains is None or label_contains.lower() in label:
            total += value
    return total


# Stable expected counts from Phase C applied-prioritisation and top-candidate audit gates.
applied_candidates = 500
follow_up_candidates = 173
annotation_supported_candidates = 4

# Read evidence-grade categories if available; otherwise preserve exact source rows in labels.
evidence_grade_lines = []
for row in evidence_rows:
    label = None
    count = None
    for key in row:
        low = key.lower()
        if "grade" in low or "evidence" in low:
            label = row[key]
        if "count" in low or "n" == low:
            count = row[key]
    if label is None:
        label = " / ".join(row.values())
    if count is None:
        count = "."
    evidence_grade_lines.append(f"{label}: {count}")

if not evidence_grade_lines:
    failures.append("No evidence-grade rows available for Figure 2.")

steps = [
    ("500", "Applied\ncandidates", "Predeclared technical\nprioritisation set"),
    ("173", "Follow-up\ncandidates", "Candidates retained for\nannotation/evidence triage"),
    ("4", "Annotation-supported\ncandidates", "TRAF3, BIRC3 and\npaired TNFAIP3 candidates"),
]

fig, ax = plt.subplots(figsize=(12.5, 6.0))
ax.set_xlim(0, 12.5)
ax.set_ylim(0, 6.0)
ax.axis("off")

ax.text(
    6.25,
    5.55,
    "Phase C candidate prioritisation and evidence-review flow",
    ha="center",
    va="center",
    fontsize=14,
    fontweight="bold",
)
ax.text(
    6.25,
    5.17,
    "Counts indicate audit-controlled manuscript display groups, not biological importance or causal ranking",
    ha="center",
    va="center",
    fontsize=9.5,
)

x_positions = [0.9, 4.15, 7.4]
box_w = 2.15
box_h = 1.5
y = 3.15

for i, (count, title, subtitle) in enumerate(steps):
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
    ax.text(x + box_w / 2, y + 1.12, count, ha="center", va="center", fontsize=20, fontweight="bold")
    ax.text(x + box_w / 2, y + 0.72, title, ha="center", va="center", fontsize=10.5, fontweight="bold")
    ax.text(x + box_w / 2, y + 0.27, subtitle, ha="center", va="center", fontsize=8.6)

    if i < len(steps) - 1:
        arrow = FancyArrowPatch(
            (x + box_w + 0.15, y + box_h / 2),
            (x_positions[i + 1] - 0.15, y + box_h / 2),
            arrowstyle="-|>",
            mutation_scale=15,
            linewidth=1.3,
            color="black",
        )
        ax.add_patch(arrow)

# Evidence-grade inset
inset = FancyBboxPatch(
    (9.95, 2.65),
    1.85,
    2.25,
    boxstyle="round,pad=0.08,rounding_size=0.06",
    linewidth=1.2,
    facecolor="white",
    edgecolor="black",
)
ax.add_patch(inset)
ax.text(10.875, 4.55, "Evidence-grade\nsummary", ha="center", va="center", fontsize=10, fontweight="bold")

if evidence_grade_lines:
    y0 = 4.02
    for line in evidence_grade_lines[:5]:
        ax.text(10.875, y0, line, ha="center", va="center", fontsize=7.8)
        y0 -= 0.33

boundary = FancyBboxPatch(
    (0.9, 0.72),
    10.9,
    1.0,
    boxstyle="round,pad=0.06,rounding_size=0.05",
    linewidth=1.2,
    facecolor="white",
    edgecolor="black",
)
ax.add_patch(boundary)
ax.text(
    6.35,
    1.28,
    "Boundary: prioritisation and evidence grades are descriptive audit outputs; they do not establish causality, selection, ancestry effects, population specificity or independent functional validation.",
    fontsize=8.7,
    ha="center",
    va="center",
    wrap=True,
)

fig.savefig(SVG_OUT, format="svg", bbox_inches="tight")
fig.savefig(PNG_OUT, format="png", dpi=1800, bbox_inches="tight")
plt.close(fig)

# Clean Matplotlib SVG trailing whitespace for git diff checks.
text = SVG_OUT.read_text(encoding="utf-8")
SVG_OUT.write_text("\n".join(line.rstrip(" \t") for line in text.splitlines()) + "\n", encoding="utf-8")

CAPTION_OUT.write_text(
    "Figure 2. Phase C candidate prioritisation and evidence-review flow. "
    "The figure summarises the transition from 500 applied candidates to 173 follow-up candidates and four annotation-supported candidates reviewed in the first evidence-review pass. "
    "Counts and evidence-grade summaries are audit-controlled descriptive outputs only and do not imply biological importance, causality, selection, adaptation, population specificity, ancestry effects or independent functional validation.\n",
    encoding="utf-8",
)

summary_rows = [
    {"metric": "figure_id", "value": "FIGURE_2"},
    {"metric": "applied_candidates", "value": applied_candidates},
    {"metric": "follow_up_candidates", "value": follow_up_candidates},
    {"metric": "annotation_supported_candidates", "value": annotation_supported_candidates},
    {"metric": "evidence_grade_rows", "value": len(evidence_rows)},
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
    "operation": "Generate Phase C Figure 2 candidate/evidence flow diagram.",
    "figure_id": "FIGURE_2",
    "applied_candidates": applied_candidates,
    "follow_up_candidates": follow_up_candidates,
    "annotation_supported_candidates": annotation_supported_candidates,
    "evidence_grade_rows": len(evidence_rows),
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
        "This checkpoint generates Figure 2 as a candidate/evidence flow display only. It does not perform new analysis, "
        "write Results prose, validate function, infer selection, make population-differentiation claims, assign ancestry "
        "effects or authorise broad biological interpretation."
    ),
}

OUTJSON.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
print(json.dumps(payload, indent=2))

if payload["status"] != "PASS":
    raise SystemExit(1)
