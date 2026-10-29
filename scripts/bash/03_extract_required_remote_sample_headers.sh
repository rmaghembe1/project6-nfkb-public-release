#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$ROOT"

INVENTORY="04_gene_panel_and_regions/nfkb_required_chromosome_query_inventory_GRCh38.tsv"
OUTDIR="03_raw_genomic_resources/sample_compatibility"
LOGDIR="15_logs_and_audit/sample_compatibility"

GENOTYPE_BASE="https://storage.googleapis.com/gcp-public-data--gnomad/release/3.1.2/vcf/genomes"
PHASED_BASE="https://storage.googleapis.com/gcp-public-data--gnomad/resources/hgdp_1kg/phased_haplotypes_v2"

mkdir -p "$OUTDIR" "$LOGDIR"

LOG="$LOGDIR/required_chromosome_remote_sample_header_access_2026-05-27.log"
: > "$LOG"

{
  echo "Required-chromosome remote sample-header extraction: $(date -Is)"
  echo
  echo "bcftools version:"
  bcftools --version | head -n 2
  echo
  echo -e "chromosome\tgenotype_samples\tphased_samples\tphased_object"
} | tee -a "$LOG"

tail -n +2 "$INVENTORY" | cut -f1 | while read -r chrom; do
    genotype_url="${GENOTYPE_BASE}/gnomad.genomes.v3.1.2.hgdp_tgp.${chrom}.vcf.bgz"
    genotype_out="${OUTDIR}/gnomad_v3.1.2_hgdp_tgp_${chrom}_genotype_samples.txt"
    genotype_warn="${OUTDIR}/gnomad_v3.1.2_hgdp_tgp_${chrom}_genotype_header_warnings.log"

    if [ "$chrom" = "chrX" ]; then
        phased_label="chrX_nonPAR"
        phased_object="hgdp1kgp_chrX_non_par.full.shapeit5_rare.bcf"
    else
        phased_label="$chrom"
        phased_object="hgdp1kgp_${chrom}.filtered.SNV_INDEL.phased.shapeit5.bcf"
    fi

    phased_url="${PHASED_BASE}/${phased_object}"
    phased_out="${OUTDIR}/hgdp1kgp_${phased_label}_phased_haplotype_samples.txt"
    phased_warn="${OUTDIR}/hgdp1kgp_${phased_label}_phased_header_warnings.log"

    echo "Reading genotype header: ${chrom}" >&2
    timeout 300 bcftools query -l "$genotype_url" \
      > "$genotype_out" 2> "$genotype_warn"

    echo "Reading phased header: ${phased_label}" >&2
    timeout 300 bcftools query -l "$phased_url" \
      > "$phased_out" 2> "$phased_warn"

    genotype_count=$(wc -l < "$genotype_out")
    phased_count=$(wc -l < "$phased_out")

    printf '%s\t%s\t%s\t%s\n' \
      "$chrom" "$genotype_count" "$phased_count" "$phased_object" \
      | tee -a "$LOG"
done

echo "PASS: remote sample-header extraction completed for all required chromosomes." \
  | tee -a "$LOG"
