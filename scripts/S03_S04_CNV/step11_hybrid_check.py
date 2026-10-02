import subprocess
import pandas as pd
import os

# ══════════════════════════════════════════════════════
# 0. 설정
# ══════════════════════════════════════════════════════
CD_GWAS = "/input/data/CDKDC.250120.firthfall.04.log.DIS.glm.logistic.hybrid"
UC_GWAS = "/input/data/UC526KDC.03.ff.log.Pheno.glm.logistic.hybrid"

OUT_DIR = "snp_lookup"
os.makedirs(OUT_DIR, exist_ok=True)

WINDOW_HALF = 500_000  # CNV candidate 중심좌표 기준 ±500kb (AnnotSV 때와 동일 기준)

# (label, gwas_file, chrom, win_mid, CNV_OR, CNV_p, note)
CANDIDATES = [
    ("MACROD2",    CD_GWAS, "20", 14_500_000, 3.36,  0.037, "CD, DEL, risk 방향"),
    ("FMNL2",      CD_GWAS, "2",  152_500_000, 0.78,  0.666, "CD, DEL, 무신호(p=0.666)"),
    ("BHLHE40",    UC_GWAS, "3",  4_500_000,  3.19,  0.030, "UC, DEL, risk 방향"),
    ("1q21_SHC1",  UC_GWAS, "1",  151_500_000, 0.363, 0.013, "UC, DUP, protective 방향"),
    ("chr10_57.5M", UC_GWAS, "10", 57_500_000, None,  2.3e-5, "UC, genome-wide significant CNV window"),
]

# ══════════════════════════════════════════════════════
# 1. awk로 좌표 구간만 우선 추출 (전체 로드 방지)
# ══════════════════════════════════════════════════════
def extract_region(gwas_path, chrom, pos_min, pos_max, out_path):
    """대용량 genome-wide summary stats에서 지정 구간만 awk로 뽑아 out_path에 저장"""
    awk_cmd = (
        f"awk -F'\\t' 'NR==1 || ($1==\"{chrom}\" && $2>={pos_min} && $2<={pos_max})' "
        f"{gwas_path} > {out_path}"
    )
    subprocess.run(awk_cmd, shell=True, check=True)

for label, gwas_path, chrom, win_mid, cnv_or, cnv_p, note in CANDIDATES:
    if not os.path.exists(gwas_path):
        print(f"[경고] GWAS 파일 없음: {gwas_path}")
        continue
    out_path = os.path.join(OUT_DIR, f"{label}_snps.tsv")
    extract_region(gwas_path, chrom, win_mid - WINDOW_HALF, win_mid + WINDOW_HALF, out_path)
    n_lines = sum(1 for _ in open(out_path)) - 1
    print(f"[추출] {label}: {n_lines}개 SNP -> {out_path}")

# ══════════════════════════════════════════════════════
# 2. 각 locus의 lead SNP(최소 P) 조회 + 방향성 비교
# ══════════════════════════════════════════════════════
def cnv_direction(cnv_or):
    if cnv_or is None:
        return "N/A(genome-wide sig, 개별 OR 별도 확인 필요)"
    return "risk(OR>1)" if cnv_or > 1 else "protective(OR<1)"

def snp_direction(or_val):
    return "risk(OR>1)" if or_val > 1 else "protective(OR<1)"

rows = []
for label, gwas_path, chrom, win_mid, cnv_or, cnv_p, note in CANDIDATES:
    path = os.path.join(OUT_DIR, f"{label}_snps.tsv")
    if not os.path.exists(path):
        continue

    df = pd.read_csv(path, sep="\t")
    df.columns = [c.strip() for c in df.columns]

    # ERRCODE가 '.'(에러 없음)인 것만, P가 숫자로 유효한 것만 사용
    df = df[df["ERRCODE"] == "."]
    df["P"] = pd.to_numeric(df["P"], errors="coerce")
    df["OR"] = pd.to_numeric(df["OR"], errors="coerce")
    df = df.dropna(subset=["P", "OR"])

    if df.empty:
        rows.append({
            "locus": label, "n_snps_valid": 0,
            "lead_SNP": None, "lead_POS": None, "lead_OR": None, "lead_P": None,
            "SNP_direction": "N/A", "CNV_direction": cnv_direction(cnv_or),
            "concordant": "N/A(SNP 없음)", "note": note,
        })
        continue

    lead = df.loc[df["P"].idxmin()]
    snp_dir = snp_direction(lead["OR"])
    cnv_dir = cnv_direction(cnv_or)
    concordant = (
        "N/A" if cnv_or is None else
        ("일치" if (("risk" in snp_dir) == ("risk" in cnv_dir)) else "불일치")
    )

    rows.append({
        "locus": label,
        "n_snps_valid": len(df),
        "lead_SNP": lead["ID"], "lead_POS": int(lead["POS"]),
        "lead_OR": round(lead["OR"], 3), "lead_P": lead["P"],
        "SNP_direction": snp_dir, "CNV_direction": cnv_dir,
        "concordant": concordant, "note": note,
    })

result_df = pd.DataFrame(rows)
result_df.to_csv("cnv_snp_direction_comparison.tsv", sep="\t", index=False)

print(f"\n[저장] cnv_snp_direction_comparison.tsv\n")
print(result_df.to_string(index=False))
