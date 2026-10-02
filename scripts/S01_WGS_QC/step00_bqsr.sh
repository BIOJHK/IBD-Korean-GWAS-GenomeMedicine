#!/usr/bin/env bash
# Step 0: Compute the Base Quality Score Recalibration (BQSR) table with
# GATK BaseRecalibrator, ahead of applying it (step01_apply_bqsr.sh).
# See Methods: "... following the GATK best practices."

set -euo pipefail

QUEUE="${QUEUE:-your_queue.q}"
NOTIFY_EMAIL="${NOTIFY_EMAIL:-your_email@example.com}"
SLOTS="${SLOTS:-12}"
REF="Ref/Homo_sapiens_assembly38_reduced.fa"
KNOWN_SITES="GCF_000001405.39.Chr.vcf.gz"

for bam in Results/040*fixmate.bam; do
  sample=$(basename "$bam" .out)
  cmd="gatk BaseRecalibrator \
    --tmp-dir TMP/ \
    --input ${bam} \
    --output Results/${sample}.BQSR.tbl \
    --reference ${REF} \
    --known-sites ${KNOWN_SITES} \
    1> Results/${sample}.BQSR.out \
    2> Results/${sample}.BQSR.err"

  echo "$cmd" | qsub -N "BQS.${sample}" -M "$NOTIFY_EMAIL" -m e -cwd \
    -pe pe_slots "$SLOTS" -V -q "$QUEUE"
done
