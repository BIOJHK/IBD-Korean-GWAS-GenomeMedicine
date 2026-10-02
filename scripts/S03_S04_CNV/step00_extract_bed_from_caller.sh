#!/bin/bash
# ============================================================
# Step 1: Delly CNV BCF / CNVkit .cns -> DEL/DUP BED 추출
# ============================================================
# 요구사항: bcftools, awk
#
# 결과물 (4개 파일, 모두 아래 공통 컬럼 구조):
#   Chromosome  Start  End  Type(DEL/DUP)  Length  Sample_Raw  BED_Status
#
#   delly_DEL.bed   delly_DUP.bed
#   cnvkit_DEL.bed  cnvkit_DUP.bed
#
# 이 4개 파일을 step2_reciprocal_merge_multi_threshold.sh 에 그대로 사용합니다.
# ============================================================
set -uo pipefail   # 개별 파일 오류로 전체가 중단되지 않도록 -e는 빼둠

# ------------------------------------------------------------
# CONFIG — 경로/옵션을 실제 환경에 맞게 수정하세요
# ------------------------------------------------------------
DELLY_DIR="/input/data/CNV_all/DELLY"          # 재귀적으로 *.cnv.bcf 탐색
CNVKIT_DIRS=(
    "/input/data/CNV_all/CNVKIT/IBD_cns"
    "/input/data/CNV_all/CNVKIT/KDC_cns"
)
OUT_DIR="/input/data/CNV_all/extracted_bed"

# Delly: PASS만 쓸지(권장) LowQual까지 포함할지
DELLY_PASS_ONLY=true

# CNVkit: germline 샘플 권장 threshold (CNVkit 공식 문서 기준).
#   tumor 기본값은 -0.25/+0.2 이지만, germline은 -0.4/+0.3이 더 정밀하다고 명시되어 있음.
CNVKIT_LOSS_THRESH=-0.4
CNVKIT_GAIN_THRESH=0.3

mkdir -p "$OUT_DIR"

# ============================================================
# ★ 샘플 ID 추출 규칙 — 반드시 실제 파일명 패턴에 맞게 확인/수정하세요.
#   KDC0001_SoMd.cnv.bcf              -> KDC0001
#   123456_sorted_MarkDuplicate.cnv.bcf -> 123456
#   123456.cov.cns                    -> 123456
#   KDC0001.cov.cns                   -> KDC0001
# CNVkit과 Delly의 샘플 ID가 서로 정확히 일치해야 이후 병합(같은 사람끼리
# 매칭)이 제대로 됩니다. 두 caller의 명명 규칙이 다르면 여기서 반드시
# 통일시켜 주세요 (예: 앞의 0 유무, KDC/IBD 접두어 유무 등).
# ============================================================
delly_sample_id () { echo "$1" | sed -E 's/(_SoMd)?(_sorted_MarkDuplicate)?\.cnv\.bcf$//'; }
cnvkit_sample_id () { echo "$1" | sed -E 's/\.cov\.cns$//; s/\.cns$//'; }

# ============================================================
# 1) Delly BCF -> BED
# ============================================================
echo "=> [1/2] Delly BCF에서 DEL/DUP 추출..."
DELLY_RAW="$OUT_DIR/delly_all_raw.tsv"
: > "$DELLY_RAW"

FILTER_EXPR='INFO/SVTYPE="CNV"'
if [ "$DELLY_PASS_ONLY" = true ]; then
    FILTER_EXPR='INFO/SVTYPE="CNV" && FILTER="PASS"'
fi

n_files=0
find "$DELLY_DIR" -name "*.cnv.bcf" | while read -r f; do
    n_files=$((n_files + 1))
    base=$(basename "$f")
    sid=$(delly_sample_id "$base")

    bcftools query -f '%CHROM\t%POS\t%INFO/END\t%FILTER\t[%CN]\n' \
        -i "$FILTER_EXPR" "$f" 2>>"$OUT_DIR/delly_extract_errors.log" \
      | awk -v sample="$sid" 'BEGIN{OFS="\t"}
        {
          chrom=$1; sub(/^chr/,"",chrom);
          start=$2-1; end=$3; cn=$5;
          if (cn=="." || cn=="") next;
          if (cn+0==2) next;                 # 정상 카피수는 CNV가 아니므로 제외
          type=(cn+0<2)?"DEL":"DUP";
          len=end-start;
          if (len<=0) next;
          print chrom, start, end, type, len, sample, "PASS"
        }' >> "$DELLY_RAW"
done
echo "   Delly raw 콜: $(wc -l < "$DELLY_RAW")개"
[ -s "$OUT_DIR/delly_extract_errors.log" ] && \
    echo "   [주의] 일부 파일에서 오류 발생 — delly_extract_errors.log 확인 필요"

grep -P '\tDEL\t' "$DELLY_RAW" > "$OUT_DIR/delly_DEL.bed"
grep -P '\tDUP\t' "$DELLY_RAW" > "$OUT_DIR/delly_DUP.bed"
echo "   -> delly_DEL.bed ($(wc -l < "$OUT_DIR/delly_DEL.bed")개) / delly_DUP.bed ($(wc -l < "$OUT_DIR/delly_DUP.bed")개)"

# ============================================================
# 2) CNVkit .cns -> BED
# ============================================================
echo "=> [2/2] CNVkit .cns에서 DEL/DUP 추출 (LOSS<=${CNVKIT_LOSS_THRESH}, GAIN>=${CNVKIT_GAIN_THRESH})..."
CNVKIT_RAW="$OUT_DIR/cnvkit_all_raw.tsv"
: > "$CNVKIT_RAW"

for dir in "${CNVKIT_DIRS[@]}"; do
    find "$dir" -name "*.cns" | while read -r f; do
        base=$(basename "$f")
        sid=$(cnvkit_sample_id "$base")

        awk -v sample="$sid" -v loss="$CNVKIT_LOSS_THRESH" -v gain="$CNVKIT_GAIN_THRESH" \
            'BEGIN{OFS="\t"; FS="\t"}
             NR==1 {next}                     # 헤더 줄 스킵
             {
               chrom=$1; sub(/^chr/,"",chrom);
               start=$2; end=$3; log2=$5;
               if (log2=="" || log2=="NA") next;
               log2v = log2 + 0;
               type="";
               if (log2v <= loss) type="DEL";
               else if (log2v >= gain) type="DUP";
               else next;                      # neutral 구간은 CNV 아님
               len = end-start;
               if (len<=0) next;
               print chrom, start, end, type, len, sample, "PASS"
             }' "$f" >> "$CNVKIT_RAW"
    done
done
echo "   CNVkit raw 콜: $(wc -l < "$CNVKIT_RAW")개"

grep -P '\tDEL\t' "$CNVKIT_RAW" > "$OUT_DIR/cnvkit_DEL.bed"
grep -P '\tDUP\t' "$CNVKIT_RAW" > "$OUT_DIR/cnvkit_DUP.bed"
echo "   -> cnvkit_DEL.bed ($(wc -l < "$OUT_DIR/cnvkit_DEL.bed")개) / cnvkit_DUP.bed ($(wc -l < "$OUT_DIR/cnvkit_DUP.bed")개)"

echo ""
echo "=== 완료 ==="
echo "생성된 파일:"
echo "  $OUT_DIR/delly_DEL.bed"
echo "  $OUT_DIR/delly_DUP.bed"
echo "  $OUT_DIR/cnvkit_DEL.bed"
echo "  $OUT_DIR/cnvkit_DUP.bed"
echo ""
echo "다음: step2_reciprocal_merge_multi_threshold.sh 에서 EXTRACT_DIR을"
echo "  \"$OUT_DIR\" 로 설정하고 실행하세요."
