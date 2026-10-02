#!/usr/bin/env bash
# Step 9: Apply minor-allele-frequency and missingness filtering to the
# biallelic SNP set (step07_extract_biallelic_snp.sh), then index the result.
# See Methods: "excluding variants based on the following criteria: ...
# missing rate: 1% (--max-missing 0.99)." Adjust MAX_MISSING/MAF below if
# you are reproducing a different filtering pass than the one used here.

set -euo pipefail

QUEUE="${QUEUE:-your_queue.q}"
NOTIFY_EMAIL="${NOTIFY_EMAIL:-your_email@example.com}"
SLOTS="${SLOTS:-1}"
MAX_MISSING="0.1"
MAF="0.01"

for vcf in Results4/BiSNP/*vcf.gz; do
  out=$(basename "$vcf" | cut -d '.' -f 1)
  outfile="Results4/MAF_MISS/${out}.maf${MAF}.miss0.9"
  cmd="vcftools --gzvcf ${vcf} \
    --max-missing ${MAX_MISSING} --maf ${MAF} \
    --recode --recode-INFO-all --stdout \
    2> ${outfile}.err \
    | bgzip > ${outfile}.vcf.gz 2> ${outfile}.vcf.gz.err \
    && tabix -p vcf ${outfile}.vcf.gz"

  echo "$cmd" | qsub -N "MM.${out}" -M "$NOTIFY_EMAIL" -m e -cwd \
    -pe pe_slots "$SLOTS" -V -q "$QUEUE"
done
