#!/usr/bin/env bash
# Step 0: Apply Base Quality Score Recalibration (BQSR) to aligned reads.
# See Methods: "the PCR duplicates were removed using Picard MarkDuplicates
# ... following the GATK best practices."
#
# For each fixmate BAM, submits one GATK ApplyBQSR job to the cluster job
# scheduler (SGE/qsub). Adjust QUEUE / NOTIFY_EMAIL / SLOTS below to match
# your own cluster environment before running.

set -euo pipefail

# --- Cluster / job-scheduler settings (adjust to your environment) ---
QUEUE="${QUEUE:-your_queue.q}"
NOTIFY_EMAIL="${NOTIFY_EMAIL:-your_email@example.com}"
SLOTS="${SLOTS:-12}"
REF="Ref/Homo_sapiens_assembly38_reduced.fa"

for bam in Results/040*fixmate.bam; do
  sample=$(basename "$bam" .out)
  cmd="gatk ApplyBQSR \
    --tmp-dir TMP/ \
    --input ${bam} \
    --bqsr-recal-file Results/${sample}.BQSR.tbl \
    --output Results/${sample}.BQSR.bam \
    --reference ${REF} \
    1> Results/${sample}.BQSRappl.out \
    2> Results/${sample}.BQSRappl.err"

  echo "$cmd" | qsub -N BQSa -M "$NOTIFY_EMAIL" -m e -cwd \
    -pe pe_slots "$SLOTS" -V -q "$QUEUE"
done
