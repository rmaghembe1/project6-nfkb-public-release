#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$ROOT"
OUT="logs/phaseA_remote_resource_access_check.log"
: > "$OUT"
echo "Remote resource access check: $(date -Is)" | tee -a "$OUT"
# Data URLs/objects are logged, not committed as downloaded large data files.
GNOMAD_META="gs://gcp-public-data--gnomad/release/3.1/secondary_analyses/hgdp_1kg_v2/metadata_and_qc/gnomad_meta_updated.tsv"
GNOMAD_PHASED="gs://gcp-public-data--gnomad/resources/hgdp_1kg/phased_haplotypes_v2/"
if command -v gcloud >/dev/null 2>&1; then
  gcloud storage ls "$GNOMAD_META" | tee -a "$OUT"
  gcloud storage ls "$GNOMAD_PHASED" | sed -n '1,20p' | tee -a "$OUT"
elif command -v gsutil >/dev/null 2>&1; then
  gsutil ls "$GNOMAD_META" | tee -a "$OUT"
  gsutil ls "$GNOMAD_PHASED" | sed -n '1,20p' | tee -a "$OUT"
else
  echo "SKIP: install Google Cloud CLI or gsutil to verify gnomAD bucket objects locally." | tee -a "$OUT"
fi
for url in \
  "https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/release_1.5/MANE.GRCh38.v1.5.summary.txt.gz" \
  "https://screen.encodeproject.org/" \
  "https://www.ensembl.org/vep"; do
  if command -v curl >/dev/null 2>&1; then
    code=$(curl -L -s --connect-timeout 15 --max-time 45 -o /dev/null -w "%{http_code}" "$url" || true)
    echo "$url HTTP=$code" | tee -a "$OUT"
  fi
done
echo "Review output and record exact downloadable objects/checksums before Phase B." | tee -a "$OUT"
