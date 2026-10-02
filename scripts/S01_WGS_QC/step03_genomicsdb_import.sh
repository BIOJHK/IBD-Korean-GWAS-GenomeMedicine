#!/usr/bin/env bash
# Step 3: Consolidate per-sample GVCFs into a GenomicsDB workspace per
# genomic interval with GATK GenomicsDBImport, ahead of joint genotyping
# (step04_genotype_gvcfs.sh).

set -euo pipefail

QUEUE="${QUEUE:-your_queue.q}"
NOTIFY_EMAIL="${NOTIFY_EMAIL:-your_email@example.com}"
SLOTS="${SLOTS:-16}"
REF="Ref/Homo_sapiens_assembly38_reduced.fa"
INTERVALS="Ref/Homo_sapiens_assembly38_reduced.intervals"
GVCF_LIST="Results2/gvcf2.list1.txt"

# Build the "-V <gvcf> -V <gvcf> ..." argument list once; it's the same
# for every interval below.
gvcf_args=""
while read -r gvcf_file; do
  [ -z "$gvcf_file" ] && continue
  gvcf_args="${gvcf_args} -V ${gvcf_file}"
done < "$GVCF_LIST"

while read -r interval; do
  [ -z "$interval" ] && continue
  name=$(echo "$interval" | cut -d '.' -f 2)
  workspace="Normal3000_${name}"
  cmd="gatk GenomicsDBImport \
    --tmp-dir ./TMP \
    -L ${interval} \
    -R ${REF} \
    ${gvcf_args} \
    --batch-size 50 \
    --consolidate true \
    --genomicsdb-workspace-path GenomeDB/${workspace} \
    --genomicsdb-shared-posixfs-optimizations true \
    1> ./Results3/GenomicDBimport/${workspace}.log \
    2> ./Results3/GenomicDBimport/${workspace}.err"

  echo "$cmd" | qsub -N "DB.${name}" -M "$NOTIFY_EMAIL" -m e -cwd \
    -pe pe_slots "$SLOTS" -V -q "$QUEUE"
done < "$INTERVALS"
