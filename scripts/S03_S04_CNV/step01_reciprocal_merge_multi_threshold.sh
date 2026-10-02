

#!/bin/bash
# ============================================================
# Step 2: 샘플별 reciprocal overlap 병합 (여러 threshold 동시 생성)
# ============================================================
# 요구사항: bedtools, awk, sort
#
# 핵심 아이디어 ("가상 염색체" 트릭):
#   bedtools intersect는 "같은 사람의 콜끼리만 비교"라는 개념이 없습니다.
#   그래서 실제 염색체 이름 앞에 샘플ID를 붙여(예: "KDC0001__1") 마치 각
#   샘플이 서로 다른 염색체 세트를 가진 것처럼 인코딩하면, bedtools가
#   자동으로 "같은 샘플 + 같은 염색체"인 콜끼리만 겹침을 계산하게 됩니다.
#   수천 명을 일일이 for-loop로 도는 것보다 훨씬 빠릅니다.
#
# 출력: merged_consensus_02.bed / _03.bed / _04.bed / _05.bed
#   (기존 메인/진단 스크립트가 기대하는 것과 동일한 컬럼 구조)
#   Chromosome  Start  End  Type  Length  Sample_Raw  BED_Status
# ============================================================
set -euo pipefail

# ------------------------------------------------------------
# CONFIG
# ------------------------------------------------------------
EXTRACT_DIR="/input/data/CNV_all/extracted_bed"     # step1 결과 위치
OUT_DIR="/input/data/CNV_all/merged_by_threshold"

THRESHOLDS=(0.2 0.3 0.4 0.5)

mkdir -p "$OUT_DIR"

DELLY_DEL="$EXTRACT_DIR/delly_DEL.bed"
DELLY_DUP="$EXTRACT_DIR/delly_DUP.bed"
CNVKIT_DEL="$EXTRACT_DIR/cnvkit_DEL.bed"
CNVKIT_DUP="$EXTRACT_DIR/cnvkit_DUP.bed"

for f in "$DELLY_DEL" "$DELLY_DUP" "$CNVKIT_DEL" "$CNVKIT_DUP"; do
    [ -s "$f" ] || { echo "[오류] $f 가 없거나 비어있습니다. step1을 먼저 실행하세요."; exit 1; }
done

# ------------------------------------------------------------
# 가상 염색체 인코딩: Sample_Raw + "__" + 원래 염색체
# ------------------------------------------------------------
make_virtual_bed () {
    # $1 = 입력 bed (Chromosome Start End Type Length Sample_Raw BED_Status)
    # $2 = 출력 bed (vchrom Start End origChrom Type Length Sample_Raw)
    awk 'BEGIN{OFS="\t"}
         { vchrom = $6 "__" $1;
           print vchrom, $2, $3, $1, $4, $5, $6 }' "$1" \
      | sort -k1,1 -k2,2n > "$2"
}

echo "=> 가상 염색체 인코딩 (샘플ID + 염색체)..."
make_virtual_bed "$DELLY_DEL"  "$OUT_DIR/delly_DEL.vchrom.bed"
make_virtual_bed "$DELLY_DUP"  "$OUT_DIR/delly_DUP.vchrom.bed"
make_virtual_bed "$CNVKIT_DEL" "$OUT_DIR/cnvkit_DEL.vchrom.bed"
make_virtual_bed "$CNVKIT_DUP" "$OUT_DIR/cnvkit_DUP.vchrom.bed"

for T in "${THRESHOLDS[@]}"; do
    # 0.3 -> "03", 0.5 -> "05"  (기존 merged_consensus_03.bed / _05.bed 명명 규칙과 동일)
    TAG=$(awk -v t="$T" 'BEGIN{printf "%02d", t*10 + 0.5}')
    echo ""
    echo "=> threshold=${T} (merged_consensus_${TAG}.bed) 병합 중..."

    for TYPE in DEL DUP; do
        DELLY_V="$OUT_DIR/delly_${TYPE}.vchrom.bed"
        CNVKIT_V="$OUT_DIR/cnvkit_${TYPE}.vchrom.bed"

        # -a: Delly (breakpoint 정밀도가 상대적으로 높다고 알려져 있어 좌표 기준으로 우선 사용)
        # -b: CNVkit
        # -f T -r : 두 콜이 서로 T 비율 이상 reciprocal overlap일 때만 매칭
        #           (가상 염색체 덕분에 자동으로 "같은 샘플"끼리만 비교됨)
        bedtools intersect -a "$DELLY_V" -b "$CNVKIT_V" -f "$T" -r -wa -wb \
            > "$OUT_DIR/consensus_${TYPE}_${TAG}.raw.tsv"

        n=$(wc -l < "$OUT_DIR/consensus_${TYPE}_${TAG}.raw.tsv")
        echo "   ${TYPE}: ${n}개 매칭"
    done

    # ------------------------------------------------------------
    # -wa -wb 결과(14컬럼: Delly 7 + CNVkit 7)를 최종 스키마로 변환.
    # 좌표는 두 콜의 union(바깥쪽 경계)을 사용합니다.
    # ------------------------------------------------------------
    awk 'BEGIN{OFS="\t"}
         {
           # 1:vchrom 2:start 3:end 4:origChrom 5:type 6:length 7:sample  (Delly, A)
           # 8:vchrom 9:start 10:end 11:origChrom 12:type 13:length 14:sample (CNVkit, B)
           chrom  = $4;
           sample = $7;
           type   = $5;
           s = ($2 < $9)  ? $2 : $9;
           e = ($3 > $10) ? $3 : $10;
           len = e - s;
           print chrom, s, e, type, len, sample, "consensus_2caller"
         }' "$OUT_DIR/consensus_DEL_${TAG}.raw.tsv" "$OUT_DIR/consensus_DUP_${TAG}.raw.tsv" \
        > "$OUT_DIR/merged_consensus_${TAG}.bed"

    n_total=$(wc -l < "$OUT_DIR/merged_consensus_${TAG}.bed")
    echo "   -> merged_consensus_${TAG}.bed  (총 ${n_total}개)"
done

echo ""
echo "=== 완료 ==="
echo ""
echo "[주의] bedtools intersect -wa -wb는 하나의 Delly 콜이 여러 CNVkit 콜과"
echo "겹치면(혹은 그 반대) 여러 줄을 출력할 수 있습니다. 결과 파일에서 동일"
echo "(sample, chrom, type)에 대해 거의 중복된 행이 여러 개 보이면, 아래처럼"
echo "pandas에서 가장 큰 overlap 하나만 남기는 후처리를 추가하는 걸 권합니다:"
echo ""
echo "  df = df.sort_values('Length', ascending=False)"
echo "         .drop_duplicates(subset=['Chromosome','Start','End','Sample_Raw','Type'])"
echo ""
echo "cnv_diagnostics_4_threshold_sensitivity.py 의 THRESHOLD_FILES에 아래 경로를 채우세요:"
for T in "${THRESHOLDS[@]}"; do
    TAG=$(awk -v t="$T" 'BEGIN{printf "%02d", t*10 + 0.5}')
    PCT=$(awk -v t="$T" 'BEGIN{printf "%d", t*100 + 0.5}')
    echo "  \"reciprocal_${PCT}pct\": \"$OUT_DIR/merged_consensus_${TAG}.bed\","
done
