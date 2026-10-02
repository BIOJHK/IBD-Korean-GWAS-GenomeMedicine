#!/usr/bin/env bash
# Step 4: Select the final variant set (drop filtered / non-variant sites)
# with GATK SelectVariants, following the hard-filtering in step03.

set -euo pipefail

QUEUE="${QUEUE:-your_queue.q}"
NOTIFY_EMAIL="${NOTIFY_EMAIL:-your_email@example.com}"
SLOTS="${SLOTS:-4}"
REF="Ref/Homo_sapiens_assembly38_reduced.fa"

for i in {1..22} X Y M; do
  vcf="Results_each_chr/IBD_4samples.chr${i}.GT.filter.vcf.gz"
  out=$(basename "$vcf" .vcf.gz)
  cmd="gatk SelectVariants \
    -R ${REF} \
    -V ${vcf} \
    --exclude-filtered true \
    --exclude-non-variants true \
    --set-filtered-gt-to-nocall true \
    -O ${out}.sel.vcf.gz \
    2> ${out}.sel.vcf.err"

  echo "$cmd" | qsub -N "SE.${i}" -M "$NOTIFY_EMAIL" -m e -cwd \
    -pe pe_slots "$SLOTS" -V -q "$QUEUE"
done
