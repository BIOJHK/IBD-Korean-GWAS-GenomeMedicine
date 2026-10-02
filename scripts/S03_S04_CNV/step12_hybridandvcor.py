import subprocess
import pandas as pd
import os
import re

# ══════════════════════════════════════════════════════
# 0. 설정
# ══════════════════════════════════════════════════════
CD_CAUSAL_DIR = "/input/data/CD_causal"
UC_CAUSAL_DIR = "/input/data/UC_causal"

CD_GWAS = "/input/data/CDKDC.250120.firthfall.04.log.DIS.glm.logistic.hybrid"
UC_GWAS = "/input/data/UC526KDC.03.ff.log.Pheno.glm.logistic.hybrid"

WINDOW_HALF = 500_000

# 22개 unique CNV candidate locus (chrom, win_mid, gene_label, disease_source)
CNV_CANDIDATES = [
    ("20", 14_500_000, "MACROD2", "CD"),
    ("5",  32_500_000, "GOLPH3", "CD/UC"),
    ("6",  500_000,    "IRF4", "CD/UC"),
    ("9",  12_500_000, "LURAP1L", "CD"),
    ("9",  11_500_000, "", "CD"),
    ("18", 68_500_000, "CCDC102B/TMX3", "CD"),
    ("1",  189_500_000, "", "CD"),
    ("2",  152_500_000, "FMNL2", "CD"),
    ("10", 66_500_000, "CTNNA3", "CD"),
    ("14", 39_500_000, "FBXO33/GEMIN2", "CD"),
    ("6",  164_500_000, "MEAT6", "CD"),
    ("10", 57_500_000, "(genome-wide sig)", "UC"),
    ("1",  151_500_000, "1q21/SHC1", "UC"),
    ("3",  4_500_000,  "BHLHE40", "UC"),
    ("1",  166_500_000, "FAM78B", "UC"),
    ("14", 42_500_000, "FSCB", "UC"),
    ("11", 25_500_000, "LUZP2", "UC"),
    ("9",  9_500_000,  "PTPRD", "UC"),
    ("11", 98_500_000, "CNTN5", "UC"),
    ("8",  5_500_000,  "CSMD1", "UC"),
    ("4",  160_500_000, "", "UC"),
    ("4",  161_500_000, "FSTL5", "UC"),
]

# 이전 SNP 대조에서 뽑았던 lead SNP들 (참고 비교용)
PREVIOUS_LEAD_SNPS = {
    "rs406348": "MACROD2 lead (CD)",
    "rs142762584": "FMNL2 lead (CD)",
    "rs117777449": "BHLHE40 lead (UC)",
    "rs78218226": "1q21_SHC1 lead (UC)",
    "rs55846923": "chr10_57.5M lead (UC)",
}

# ══════════════════════════════════════════════════════
# 1. causal SNP 폴더에서 rsID 목록 추출
# ══════════════════════════════════════════════════════
def extract_rsids_from_dir(dir_path):
    if not os.path.isdir(dir_path):
        print(f"[경고] 폴더 없음: {dir_path}")
        return []
    files = os.listdir(dir_path)
    rsids = set()
    for f in files:
        rsids.update(re.findall(r"rs\d+", f))
    return sorted(rsids)

cd_causal_rsids = extract_rsids_from_dir(CD_CAUSAL_DIR)
uc_causal_rsids = extract_rsids_from_dir(UC_CAUSAL_DIR)
print(f"[목록] CD_causal 고유 rsID: {len(cd_causal_rsids)}개")
print(f"[목록] UC_causal 고유 rsID: {len(uc_causal_rsids)}개")

# ══════════════════════════════════════════════════════
# 2. 각 rsID의 좌표를 GWAS 파일에서 grep으로 조회
# ══════════════════════════════════════════════════════
def lookup_positions(rsid_list, gwas_path, tmp_id_file):
    if not rsid_list or not os.path.exists(gwas_path):
        return pd.DataFrame(columns=["CHROM", "POS", "ID"])
    with open(tmp_id_file, "w") as f:
        f.write("\n".join(rsid_list) + "\n")

    header = subprocess.run(f"head -1 {gwas_path}", shell=True,
                             capture_output=True, text=True).stdout.strip().split("\t")
    out_path = tmp_id_file + ".hits.tsv"
    cmd = f"grep -Fwf {tmp_id_file} {gwas_path} > {out_path}"
    subprocess.run(cmd, shell=True)

    if os.path.getsize(out_path) == 0:
        return pd.DataFrame(columns=["CHROM", "POS", "ID"])

    df = pd.read_csv(out_path, sep="\t", header=None, names=header)
    df.columns = [c.strip() for c in df.columns]
    return df[["#CHROM", "POS", "ID"]].rename(columns={"#CHROM": "CHROM"})

os.makedirs("snp_lookup", exist_ok=True)
cd_pos_df = lookup_positions(cd_causal_rsids, CD_GWAS, "snp_lookup/cd_causal_ids.txt")
uc_pos_df = lookup_positions(uc_causal_rsids, UC_GWAS, "snp_lookup/uc_causal_ids.txt")

cd_pos_df["source"] = "CD_causal"
uc_pos_df["source"] = "UC_causal"
all_causal = pd.concat([cd_pos_df, uc_pos_df], ignore_index=True)
all_causal["CHROM"] = all_causal["CHROM"].astype(str)
all_causal["POS"] = pd.to_numeric(all_causal["POS"], errors="coerce")
all_causal.to_csv("snp_lookup/all_causal_snps_with_pos.tsv", sep="\t", index=False)
print(f"[저장] 좌표 확인된 causal SNP {len(all_causal)}개 -> snp_lookup/all_causal_snps_with_pos.tsv")

# ══════════════════════════════════════════════════════
# 3. 22개 CNV candidate와 좌표 겹침 확인
# ══════════════════════════════════════════════════════
overlap_rows = []
for chrom, win_mid, gene, disease in CNV_CANDIDATES:
    win_start, win_end = win_mid - WINDOW_HALF, win_mid + WINDOW_HALF
    hits = all_causal[
        (all_causal["CHROM"] == chrom) &
        (all_causal["POS"] >= win_start) & (all_causal["POS"] <= win_end)
    ]
    if len(hits) > 0:
        for _, h in hits.iterrows():
            overlap_rows.append({
                "CNV_locus": f"chr{chrom}:{win_mid:,}", "gene": gene, "CNV_disease": disease,
                "causal_SNP": h["ID"], "causal_SNP_POS": h["POS"], "causal_source": h["source"],
                "dist_from_win_mid_kb": round((h["POS"] - win_mid) / 1000, 1),
            })

overlap_df = pd.DataFrame(overlap_rows)
overlap_df.to_csv("snp_lookup/cnv_causal_snp_overlap.tsv", sep="\t", index=False)

print(f"\n[결과] 22개 CNV candidate 중 causal SNP와 겹치는 window: "
      f"{overlap_df['CNV_locus'].nunique() if len(overlap_df) else 0}개")
if len(overlap_df):
    print(overlap_df.to_string(index=False))
else:
    print("  -> 직접 겹치는 causal SNP 없음")

# ══════════════════════════════════════════════════════
# 4. 이전에 뽑은 lead SNP 5개가 causal 목록에 있는지 확인
# ══════════════════════════════════════════════════════
print(f"\n[확인] 이전 분석의 lead SNP가 causal 목록에 포함되는지:")
for rsid, note in PREVIOUS_LEAD_SNPS.items():
    in_cd = rsid in cd_causal_rsids
    in_uc = rsid in uc_causal_rsids
    status = "CD_causal에 있음" if in_cd else ("UC_causal에 있음" if in_uc else "목록에 없음")
    print(f"  {rsid} ({note}): {status}")
