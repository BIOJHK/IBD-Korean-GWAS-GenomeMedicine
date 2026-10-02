#!/usr/bin/env bash
# Step 3: Variant quality hard-filtering per chromosome with GATK3 VariantFiltration.
# See Methods: "Basic SNP extraction was performed using vcftools ..., excluding
# variants based on the following criteria: QD < 3.0, FS : 30.0, DP < 7, GQ < 10.0,
# MQRankSum < -2.0, ReadPosRankSum < -2.0, MQ < 30.0."
#
# Two variant-selection strategies are kept below, matching what was actually run:
#   (A) a fixed chromosome list (chr1-22, X, Y, M)
#   (B) whichever per-chromosome GT VCFs are actually present under Results_each_chr/
# Run whichever matches your intermediate files; (B) is the one used if both
# produce the same chromosome set.

set -euo pipefail

QUEUE="${QUEUE:-your_queue.q}"
NOTIFY_EMAIL="${NOTIFY_EMAIL:-your_email@example.com}"
SLOTS="${SLOTS:-4}"
REF="./Ref/Homo_sapiens_assembly38_reduced.fa"
GATK3_JAR="GenomeAnalysisTK.38.jar"

# Filter expressions must keep their literal double quotes when the command
# is later interpreted by the job's shell, otherwise "<"/">" are parsed as
# redirections instead of comparison operators.
filter_args=(
  --filterName LowReadPosRankSum --filterExpression "\"ReadPosRankSum < -2.0\""
  --filterName LowMQRankSum      --filterExpression "\"MQRankSum < -2.0\""
  --filterName LowQual           --filterExpression "\"QUAL < 30.0\""
  --filterName QD                --filterExpression "\"QD < 3.0\""
  --filterName FS                --filterExpression "\"FS > 30.0\""
  --filterName MQ                --filterExpression "\"MQ < 30.0\""
  --filterName DP0               --filterExpression "\"DP < 7\""
  --genotypeFilterName DP        --genotypeFilterExpression "\"DP < 7\""
  --genotypeFilterName GQ        --genotypeFilterExpression "\"GQ < 10.0\""
)

submit_variant_filtration() {
  local vcf="$1" job_name="$2"
  local out
  out=$(basename "$vcf" .vcf.gz)
  local cmd="java -jar ${GATK3_JAR} -T VariantFiltration -R ${REF} -V ${vcf} \
    ${filter_args[*]} \
    2> Results_each_chr/${out}.filter.err \
    | bgzip 1> Results_each_chr/${out}.filter.vcf.gz \
    2> Results_each_chr/${out}.filter.vcf.gz.err"

  echo "$cmd" | qsub -N "VF.${job_name}" -M "$NOTIFY_EMAIL" -m e -cwd \
    -pe pe_slots "$SLOTS" -V -q "$QUEUE"
}

# (A) Fixed chromosome list
for i in {1..22} X Y M; do
  vcf="Results_each_chr/IBD_4samples.chr${i}.GT.vcf.gz"
  submit_variant_filtration "$vcf" "$i"
done

# (B) Whichever per-chromosome GT VCFs are actually present
for vcf in Results_each_chr/*GT.vcf.gz; do
  out=$(basename "$vcf" .vcf.gz)
  submit_variant_filtration "$vcf" "$out"
done
