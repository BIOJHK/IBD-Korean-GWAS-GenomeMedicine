#!/usr/bin/env bash
# Step 7: Keep only biallelic SNPs from the SNP-only VCFs (step06_extract_snp.sh).
# See Methods: "Biallelic SNPs were further selected using bcftools (v1.14)
# with the parameters -m 2 -M 2."

set -euo pipefail

for vcf in Results_each_chr/*sel.snp.vcf.gz; do
  out=$(basename "$vcf" .vcf.gz)
  bcftools view "$vcf" -m 2 -M 2 -O z -o "Results_each_chr/${out}.snp.bi.vcf.gz" \
    2> "Results_each_chr/${out}.snp.bi.vcf.gz.err" &
done
wait
