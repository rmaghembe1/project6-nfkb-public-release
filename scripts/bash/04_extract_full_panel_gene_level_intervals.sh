#!/usr/bin/env bash
set -uo pipefail

MANIFEST="config/phaseB/full_panel_extraction/nfkb_full_panel_gene_level_extraction_manifest_GRCh38.tsv"
RAWBASE="03_raw_genomic_resources/interval_extraction/full_panel"
OUTDIR="11_results/phaseB/interval_extraction/full_panel"
AUDITDIR="15_logs_and_audit/phaseB/interval_extraction/full_panel"
STATUS="${OUTDIR}/nfkb_full_panel_gene_level_extraction_status_v1.tsv"
LOG="${AUDITDIR}/nfkb_full_panel_gene_level_extraction_run_2026-05-28.log"

mkdir -p "${RAWBASE}/genotype" "${RAWBASE}/phased" "$OUTDIR" "$AUDITDIR"

printf 'run_started\t%s\n' "$(date -Is)" > "$LOG"
printf 'bcftools\t%s\n' "$(bcftools --version | head -n 1)" >> "$LOG"
printf 'gene_symbol\tmodule_id\tchromosome\tregion\tresource\tstatus\toutput\tindex\tsample_count\tvariant_count\tfirst_variant\tlast_variant\n' > "${STATUS}.tmp"

failures=0
completed=0

validate_and_report() {
    local gene="$1" module="$2" chrom="$3" region="$4" resource="$5" out="$6" idx="$7" status="$8"
    local samples variants first last

    samples="$(bcftools query -l "$out" 2>>"$LOG" | wc -l | tr -d ' ')"
    variants="$(bcftools view -H "$out" 2>>"$LOG" | wc -l | tr -d ' ')"
    first="$(bcftools query -f '%CHROM:%POS\n' "$out" 2>>"$LOG" | head -n 1)"
    last="$(bcftools query -f '%CHROM:%POS\n' "$out" 2>>"$LOG" | tail -n 1)"

    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
      "$gene" "$module" "$chrom" "$region" "$resource" "$status" \
      "$out" "$idx" "$samples" "$variants" "$first" "$last" \
      >> "${STATUS}.tmp"
}

extract_genotype() {
    local gene="$1" module="$2" chrom="$3" region="$4" url="$5" out="$6"
    local idx="${out}.tbi" tmp="${out}.tmp" tmpidx="${out}.tmp.tbi"

    mkdir -p "$(dirname "$out")"

    if [[ -s "$out" && -s "$idx" ]] && bcftools view -h "$out" >/dev/null 2>>"$LOG"; then
        echo "[$gene] genotype: SKIP valid existing output"
        validate_and_report "$gene" "$module" "$chrom" "$region" "genotype" "$out" "$idx" "SKIPPED_VALID_EXISTING"
        return 0
    fi

    rm -f "$tmp" "$tmpidx"
    echo "[$gene] genotype: START $region"

    if timeout 1800 bcftools view -r "$region" -Oz -o "$tmp" "$url" >>"$LOG" 2>&1 \
       && tabix -f -p vcf "$tmp" >>"$LOG" 2>&1 \
       && bcftools view -h "$tmp" >/dev/null 2>>"$LOG"; then
        mv "$tmp" "$out"
        mv "$tmpidx" "$idx"
        echo "[$gene] genotype: DONE"
        validate_and_report "$gene" "$module" "$chrom" "$region" "genotype" "$out" "$idx" "PASS"
        return 0
    fi

    rm -f "$tmp" "$tmpidx"
    echo "[$gene] genotype: FAILED"
    printf '%s\t%s\t%s\t%s\tgenotype\tFAIL\t%s\t%s\tNA\tNA\tNA\tNA\n' \
      "$gene" "$module" "$chrom" "$region" "$out" "$idx" >> "${STATUS}.tmp"
    return 1
}

extract_phased() {
    local gene="$1" module="$2" chrom="$3" region="$4" url="$5" out="$6"
    local idx="${out}.csi" tmp="${out}.tmp" tmpidx="${out}.tmp.csi"

    mkdir -p "$(dirname "$out")"

    if [[ -s "$out" && -s "$idx" ]] && bcftools view -h "$out" >/dev/null 2>>"$LOG"; then
        echo "[$gene] phased: SKIP valid existing output"
        validate_and_report "$gene" "$module" "$chrom" "$region" "phased" "$out" "$idx" "SKIPPED_VALID_EXISTING"
        return 0
    fi

    rm -f "$tmp" "$tmpidx"
    echo "[$gene] phased: START $region"

    if timeout 1800 bcftools view -r "$region" -Ob -o "$tmp" "$url" >>"$LOG" 2>&1 \
       && bcftools index -f -o "$tmpidx" "$tmp" >>"$LOG" 2>&1 \
       && bcftools view -h "$tmp" >/dev/null 2>>"$LOG"; then
        mv "$tmp" "$out"
        mv "$tmpidx" "$idx"
        echo "[$gene] phased: DONE"
        validate_and_report "$gene" "$module" "$chrom" "$region" "phased" "$out" "$idx" "PASS"
        return 0
    fi

    rm -f "$tmp" "$tmpidx"
    echo "[$gene] phased: FAILED"
    printf '%s\t%s\t%s\t%s\tphased\tFAIL\t%s\t%s\tNA\tNA\tNA\tNA\n' \
      "$gene" "$module" "$chrom" "$region" "$out" "$idx" >> "${STATUS}.tmp"
    return 1
}

while IFS=$'\t' read -r gene module chrom region geno_url phase_url geno_out phase_out; do
    [[ "$gene" == "gene_symbol" ]] && continue

    if extract_genotype "$gene" "$module" "$chrom" "$region" "$geno_url" "$geno_out"; then
        completed=$((completed + 1))
    else
        failures=$((failures + 1))
    fi

    if extract_phased "$gene" "$module" "$chrom" "$region" "$phase_url" "$phase_out"; then
        completed=$((completed + 1))
    else
        failures=$((failures + 1))
    fi
done < "$MANIFEST"

mv "${STATUS}.tmp" "$STATUS"

printf 'run_finished\t%s\n' "$(date -Is)" >> "$LOG"
printf 'successful_or_existing_outputs\t%s\n' "$completed" >> "$LOG"
printf 'failed_outputs\t%s\n' "$failures" >> "$LOG"

echo
echo "Full-panel extraction run finished."
echo "Successful/existing resource outputs: $completed / 48"
echo "Failed resource outputs: $failures"
echo "Status table: $STATUS"
echo "Run log: $LOG"

[[ "$failures" -eq 0 ]]
