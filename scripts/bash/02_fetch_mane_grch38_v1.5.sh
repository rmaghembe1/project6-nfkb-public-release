#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$ROOT"

BASE="https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/release_1.5"
DEST="02_reference_data/MANE/GRCh38_v1.5"
LOCKDIR="15_logs_and_audit/reference_locks"

mkdir -p "$DEST" "$LOCKDIR"

FILES=(
  "MANE.GRCh38.v1.5.summary.txt.gz"
  "MANE.GRCh38.v1.5.refseq_genomic.gtf.gz"
  "README_versions.txt"
)

for f in "${FILES[@]}"; do
  echo "Retrieving $f"
  rm -f "$DEST/$f.part"
  curl -fL -sS \
    --retry 3 \
    --retry-delay 5 \
    --connect-timeout 30 \
    --max-time 900 \
    "$BASE/$f" \
    -o "$DEST/$f.part"
  mv "$DEST/$f.part" "$DEST/$f"
done

gzip -t "$DEST/MANE.GRCh38.v1.5.summary.txt.gz"
gzip -t "$DEST/MANE.GRCh38.v1.5.refseq_genomic.gtf.gz"

sha256sum "$DEST"/* \
  > "$LOCKDIR/MANE_GRCh38_v1.5_sha256_latest.txt"

echo "PASS: MANE GRCh38 v1.5 retrieval, gzip validation and SHA256 locking completed."
