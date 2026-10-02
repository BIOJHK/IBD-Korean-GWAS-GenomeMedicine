#!/usr/bin/env bash
# Step 2: Joint genotyping per genomic interval with GATK GenotypeGVCFs.
# See Methods: "We merged the individual genotype to a single massive GVCF
# format for each chromosome using CombineVariant, and the variant quality
# of the joint genotypes was recalibrated using the VQSR module in GATK."

set -euo pipefail

QUEUE="${QUEUE:-your_queue.q}"
NOTIFY_EMAIL="${NOTIFY_EMAIL:-your_email@example.com}"
SLOTS="${SLOTS:-4}"
REF="Ref/Homo_sapiens_assembly38_reduced.fa"
INTERVALS="Ref/Homo_sapiens_assembly38_reduced.intervals"

while read -r interval; do
  [ -z "$interval" ] && continue
  name=$(basename "$interval" | cut -d '.' -f 2)
  cmd="gatk GenotypeGVCFs \
    -L ${interval} \
    -R ${REF} \
    -V Results_each_chr/IBD_4samples.${name}.vcf.gz \
    -O Results_each_chr/IBD_4samples.${name}.GT.vcf.gz \
    2> Results_each_chr/IBD_4samples.${name}.GT.vcf.err"

  echo "$cmd" | qsub -N "G.${name}" -M "$NOTIFY_EMAIL" -m e -cwd \
    -pe pe_slots "$SLOTS" -V -q "$QUEUE"
done < "$INTERVALS"
