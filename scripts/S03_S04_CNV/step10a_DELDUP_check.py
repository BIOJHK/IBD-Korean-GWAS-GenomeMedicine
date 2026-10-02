import pandas as pd
import os

# ══════════════════════════════════════════════════════
# 0. 설정
# ══════════════════════════════════════════════════════
WINDOW_SIZE = 1_000_000

# 이미 생성된 cnv_frequency_by_window_{SUFFIX}.csv 경로
# (document 19/21 스크립트의 6.5단계 산출물)
FREQ_FILES = {
    "CD30": "260605_survivor_check/UC527KDC_260718CD_check_recipropal30cnv_frequency_by_window_PCadj_loo.csv",
    "CD40": "260605_survivor_check/UC527KDC_260718CD_check_recipropal40cnv_frequency_by_window_PCadj_loo.csv",
    "UC30": "260605_survivor_check/UC526KDC_260716.intcheck30cnv_frequency_by_window_PCadj_loo.csv",
    "UC40": "260605_survivor_check/UC526KDC_260716.intcheck40cnv_frequency_by_window_PCadj_loo.csv",
}

# 22개 unique candidate locus (chrom은 'Chr NN' 포맷으로 맞춤, Methods 표기와 동일)
def to_chr_label(chrom_str):
    p = str(chrom_str).replace("chr", "").strip()
    return f"Chr {int(p):02d}" if p.isdigit() else f"Chr {p}"

CANDIDATES = [
    # (analysis_key, chrom, win_mid, gene_label)
    ("CD30", "20", 14_500_000, "MACROD2"),
    ("CD30", "5",  32_500_000, "GOLPH3"),
    ("CD30", "6",  500_000,    "IRF4"),
    ("CD30", "9",  12_500_000, "LURAP1L"),
    ("CD30", "9",  11_500_000, ""),
    ("CD30", "18", 68_500_000, "CCDC102B/TMX3"),
    ("CD30", "1",  189_500_000, ""),
    ("CD30", "2",  152_500_000, "FMNL2"),
    ("CD30", "10", 66_500_000, "CTNNA3"),
    ("CD30", "14", 39_500_000, "FBXO33/GEMIN2"),
    ("CD40", "6",  164_500_000, "MEAT6"),
    ("UC30", "10", 57_500_000, "(PCDH15 인접)"),
    ("UC30", "1",  151_500_000, "1q21/SHC1"),
    ("UC30", "3",  4_500_000,  "BHLHE40"),
    ("UC30", "1",  166_500_000, "FAM78B"),
    ("UC30", "14", 42_500_000, "FSCB"),
    ("UC30", "11", 25_500_000, "LUZP2"),
    ("UC30", "9",  9_500_000,  "PTPRD"),
    ("UC30", "11", 98_500_000, "CNTN5"),
    ("UC30", "8",  5_500_000,  "CSMD1"),
    ("UC40", "4",  160_500_000, ""),
    ("UC40", "4",  161_500_000, "FSTL5"),
]

# ══════════════════════════════════════════════════════
# 1. 분석별로 freq_df 로드 (한 번만)
# ══════════════════════════════════════════════════════
freq_cache = {}
for key, path in FREQ_FILES.items():
    if not os.path.exists(path):
        print(f"[경고] 파일 없음: {path}")
        continue
    df = pd.read_csv(path)
    freq_cache[key] = df

# ══════════════════════════════════════════════════════
# 2. 각 candidate에 대해 DEL/DUP frequency 조회
# ══════════════════════════════════════════════════════
rows = []
for analysis_key, chrom, win_mid, gene_label in CANDIDATES:
    if analysis_key not in freq_cache:
        continue
    df = freq_cache[analysis_key]
    chrom_label = to_chr_label(chrom)
    window_idx = win_mid // WINDOW_SIZE

    sub = df[(df["Chromosome"] == chrom_label) & (df["Window"] == window_idx)]

    del_row = sub[sub["TypeNorm"] == "DEL"]
    dup_row = sub[sub["TypeNorm"] == "DUP"]

    del_case = del_row["Freq_case"].values[0] if len(del_row) else 0.0
    del_ctrl = del_row["Freq_ctrl"].values[0] if len(del_row) else 0.0
    dup_case = dup_row["Freq_case"].values[0] if len(dup_row) else 0.0
    dup_ctrl = dup_row["Freq_ctrl"].values[0] if len(dup_row) else 0.0

    del_dir = "case>ctrl" if del_case > del_ctrl else ("ctrl>case" if del_ctrl > del_case else "동일")
    dup_dir = "case>ctrl" if dup_case > dup_ctrl else ("ctrl>case" if dup_ctrl > dup_case else "동일")

    total_del = del_case + del_ctrl
    total_dup = dup_case + dup_ctrl
    if total_del == 0 and total_dup == 0:
        dominant = "매칭 없음"
    elif total_del >= total_dup:
        dominant = "DEL 우세"
    else:
        dominant = "DUP 우세"

    # DEL과 DUP이 둘 다 존재하면서 방향이 반대인 경우 = 신호가 상쇄될 수 있는 case
    mixed_opposite = (total_del > 0 and total_dup > 0 and del_dir != dup_dir
                       and del_dir != "동일" and dup_dir != "동일")

    rows.append({
        "analysis": analysis_key,
        "locus": f"{chrom_label.replace('Chr ', 'chr')}:{win_mid:,}",
        "gene": gene_label,
        "DEL_case_freq": round(del_case, 4), "DEL_ctrl_freq": round(del_ctrl, 4), "DEL_direction": del_dir,
        "DUP_case_freq": round(dup_case, 4), "DUP_ctrl_freq": round(dup_ctrl, 4), "DUP_direction": dup_dir,
        "dominant_type": dominant,
        "mixed_opposite_direction": mixed_opposite,
    })

out_df = pd.DataFrame(rows)
out_df.to_csv("candidate_del_dup_breakdown.tsv", sep="\t", index=False)

print(f"[저장] candidate_del_dup_breakdown.tsv\n")
print(out_df.to_string(index=False))

n_mixed = out_df["mixed_opposite_direction"].sum()
print(f"\n[요약] DEL/DUP 방향이 반대로 섞여 신호가 상쇄될 수 있는 locus: {n_mixed}개")
if n_mixed > 0:
    print("  -> 이 locus들은 DEL-only / DUP-only로 나눠서 burden test를 재실행하는 걸 권장합니다.")
