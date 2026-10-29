#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$ROOT"
mkdir -p data/raw data/external data/metadata data/processed results/qc results/phaseA results/tables results/figures logs
required=(python git)
optional=(snakemake bcftools tabix bedtools plink2 curl gsutil gcloud)
{
  echo "Project 6 local validation: $(date -Is)"
  echo "Repository: $ROOT"
  echo "Required tools:"
  for x in "${required[@]}"; do
    if command -v "$x" >/dev/null 2>&1; then echo "  PASS: $x -> $(command -v "$x")"; else echo "  FAIL: missing $x"; exit 1; fi
  done
  echo "Optional/full-analysis tools:"
  for x in "${optional[@]}"; do
    if command -v "$x" >/dev/null 2>&1; then echo "  FOUND: $x -> $(command -v "$x")"; else echo "  NOT YET FOUND: $x"; fi
  done
} | tee logs/phaseA_environment_check.log
echo "Repository directories initialized and minimal tool check completed."
