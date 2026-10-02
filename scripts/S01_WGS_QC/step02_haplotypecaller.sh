#!/usr/bin/env bash
# Step 1: Per-sample germline variant calling with GATK HaplotypeCaller (GVCF mode).
# See Methods: "Individual genotypes were called in the GVCF format by
# HaplotypeCaller in GATK with these options: --max-alternate-alleles 6,
# --native-pair-hmm-threads 4 and --emit-ref-confidence GVCF."

set -euo pipefail

QUEUE="${QUEUE:-your_queue.q}"
NOTIFY_EMAIL="${NOTIFY_EMAIL:-your_email@example.com}"
SLOTS="${SLOTS:-12}"
REF="Ref/Homo_sapiens_assembly38_reduced.fa"
DBSNP="GCF_000001405.39.Chr.vcf.gz"

for bam in Results2/*realigned.bam; do
  sample=$(basename "$bam" .out)
  cmd="gatk HaplotypeCaller \
    --tmp-dir TMP/ \
    --input ${bam} \
    --output Results2/${sample}.g.vcf \
    --reference ${REF} \
    --dbsnp ${DBSNP} \
    --max-alternate-alleles 6 \
    --native-pair-hmm-threads 4 \
    --emit-ref-confidence GVCF \
    1> Results2/${sample}.g.vcf.log \
    2> Results2/${sample}.g.vcf.err"

  echo "$cmd" | qsub -N Hap -M "$NOTIFY_EMAIL" -m e -cwd \
    -pe pe_slots "$SLOTS" -V -q "$QUEUE"
done
