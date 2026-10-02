#!/usr/bin/env bash
# Step 6: Extract SNP-only sites (drop indels) from each selected
# per-chromosome VCF (step08_select_variants.sh).
# See Methods: "Biallelic SNPs were further selected using bcftools (v1.14)."

set -euo pipefail

for vcf in Results_each_chr/*sel.vcf.gz; do
  out=$(basename "$vcf" .vcf.gz)
  bcftools view "$vcf" -V indels -O z -o "Results_each_chr/${out}.snp.vcf.gz" \
    2> "Results_each_chr/${out}.snp.vcf.gz.err" &
done
wait
