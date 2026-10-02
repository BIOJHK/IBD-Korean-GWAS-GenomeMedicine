#CNV 파이프라인 진단 스크립트 (5) — 고정 샘플셋 + TREATMENT 표준 OLS 회귀
#====================================================================================
# 기존 step02_reciprocal_merge_with_CD_check.py (diag4b, case/control MWU+permutation)
# 구조를 그대로 유지하되, phenotype을 TREATMENT(순서형 3그룹: 1=24개월내 2회+ 변경,
# 2=1회 변경, 3=변경 없음)로 바꾸고, case/control 비교 대신 표준 OLS 선형회귀로
# CNV burden ~ TREATMENT 관계를 분석합니다.
#
#  유지되는 부분:
#   1) union 데이터 기준으로 outlier를 한 번만 정하고, 그 샘플 집합을
#      모든 threshold에 동일하게 고정 (fixed sample set)
#   2) threshold별 요약 테이블(csv) + 비교 플롯(png) 저장
#
#  바뀐 부분:
#   1) outcome: CNV 개수(total_cnv), CNV 총 길이(total_length) 모두 계산하고,
#      논문에서 표준적으로 쓰는 방식대로 log1p 변환하여 연속형으로 사용
#      (zero-inflated / right-skewed count라서 raw scale OLS는 가정 위반)
#   2) predictor: TREATMENT (1/2/3, 순서형 -> 선형 dose-response로 취급)
#   3) 검정: 표준 OLS + HC3 robust SE (등분산 가정 깨짐에 대응, 논문에서 표준)
#      + 참고용으로 label-permutation 기반 slope p-value도 병행 산출
#      (원래 diag4b가 MWU와 permutation을 나란히 비교했던 철학을 그대로 유지)


import warnings
warnings.filterwarnings("ignore")

import os
import pandas as pd
import numpy as np
import scipy.stats as stats
import statsmodels.api as sm
from intervaltree import IntervalTree
import matplotlib.pyplot as plt

# ============================================================
# CONFIG
# ============================================================
PHENOTYPE_PATH  = "/input/data/UC_PHENO_GROUP3.txt"   # IID / TREATMENT 2컬럼
CENTROMERE_PATH = "260605_survivor_check/centromeres.txt.gz"
SEGDUP_PATH     = "260605_survivor_check/genomicSuperDups.txt.gz"
CYTOBAND_PATH   = "260605_survivor_check/cytoBand.txt.gz"
OUT_DIR         = "260605_survivor_check/diagnostics/UCTREAT"

PERICENTRO_BUFFER = 5_000_000
SAMPLE_SD_CUT     = 3            # outlier 판정은 union 기준으로 딱 한 번만 적용
CNV_SIZE_MIN      = 1_000
CNV_SIZE_MAX      = 3_000_000
ACROCENTRIC       = {"Chr 13", "Chr 14", "Chr 15", "Chr 21", "Chr 22"}

N_PERM       = 5000
RANDOM_SEED  = 42

# TREATMENT 컬럼: 1=24개월내 2회이상 변경(중증) / 2=1회 변경 / 3=변경없음(양호)
# 값이 클수록 "CNV burden이 낮을 것"이라는 방향으로 해석하고 싶으면 REVERSE_TREATMENT=True
# (그러면 severity = 4 - TREATMENT 로 재코딩되어, beta>0 == burden 높을수록 중증)
REVERSE_TREATMENT = False

UNION_LABEL = "union(min_supp=1)"
THRESHOLD_FILES = {
    UNION_LABEL:         "260605_survivor_check/merged_10000.v1.bed",
    "reciprocal_20pct":  "/input/data/CNV_all/merged_by_threshold/merged_consensus_02.bed",
    "reciprocal_30pct":  "/input/data/CNV_all/merged_by_threshold/merged_consensus_03.bed",
    "reciprocal_40pct":  "/input/data/CNV_all/merged_by_threshold/merged_consensus_04.bed",
    "reciprocal_50pct":  "/input/data/CNV_all/merged_by_threshold/merged_consensus_05.bed",
}


def norm(c):
    p = str(c).replace('chr', '').strip()
    return f"Chr {int(p):02d}" if p.isdigit() else f"Chr {p}"


# ============================================================
# 블랙리스트 (모든 threshold 공통) — 원본과 동일
# ============================================================
def build_blacklist_tree():
    intervals = []
    cyto = pd.read_csv(CYTOBAND_PATH, sep="\t", header=None,
                        names=["chrom", "start", "end", "band", "stain"])
    cyto['chrom'] = cyto['chrom'].apply(norm)
    band_sub = cyto[cyto['stain'].isin(["acen", "stalk", "gvar"])]
    is_acen = (band_sub['stain'] == 'acen').values
    s_vals = band_sub['start'].astype(int).values.copy()
    e_vals = band_sub['end'].astype(int).values.copy()
    s_vals[is_acen] = np.maximum(0, s_vals[is_acen] - PERICENTRO_BUFFER)
    e_vals[is_acen] = e_vals[is_acen] + PERICENTRO_BUFFER
    intervals.extend(zip(band_sub['chrom'].values, s_vals, e_vals))

    for c in ACROCENTRIC:
        parm = cyto[(cyto.chrom == c) & (cyto.band.astype(str).str.startswith("p"))]
        if len(parm):
            intervals.append((c, 0, int(parm.end.max())))

    cen = pd.read_csv(CENTROMERE_PATH, sep="\t", header=None, usecols=[1, 2, 3],
                       names=["c", "s", "e"])
    cen['c'] = cen['c'].apply(norm)
    s_buf = (cen['s'].astype(int) - PERICENTRO_BUFFER).clip(lower=0)
    e_buf = cen['e'].astype(int) + PERICENTRO_BUFFER
    intervals.extend(zip(cen['c'].values, s_buf.values, e_buf.values))

    sd = pd.read_csv(SEGDUP_PATH, sep="\t", header=None, usecols=[1, 2, 3],
                      names=["c", "s", "e"])
    sd['c'] = sd['c'].apply(norm)
    intervals.extend(zip(sd['c'].values, sd['s'].astype(int).values, sd['e'].astype(int).values))

    trees = {}
    for c, s, e in intervals:
        if s < e:
            trees.setdefault(c, IntervalTree()).addi(int(s), int(e))
    return trees


# ============================================================
# ★ phenotype 로드 — TREATMENT 버전
#    IID에 04395251 같은 leading-zero ID가 섞여있어 dtype=str 필수
# ============================================================
def load_phenotype():
    pheno = pd.read_csv(PHENOTYPE_PATH, sep=r"\s+", dtype={"IID": str}).set_index("IID")
    pheno['TREATMENT'] = pheno['TREATMENT'].astype(float)
    if REVERSE_TREATMENT:
        pheno['treatment_x'] = 4 - pheno['TREATMENT']
    else:
        pheno['treatment_x'] = pheno['TREATMENT']
    return pheno


def load_masked_filtered_cnv(bed_path, blacklist_trees):
    cnv = pd.read_csv(bed_path, sep="\t", header=None,
                       names=["Chromosome", "Start", "End", "Type", "Length", "Sample_Raw", "BED_Status"])
    cnv['Chromosome'] = cnv['Chromosome'].apply(norm)
    cnv['IID'] = (cnv['Sample_Raw'].str.replace('.survivor.chr.vcf', '', regex=False)
                                   .str.replace('_filtered.vcf', '', regex=False))
    cnv = cnv[~cnv['Chromosome'].isin(['Chr X', 'Chr Y'])].reset_index(drop=True)
    n_raw = len(cnv)

    artifact = [bool(c in blacklist_trees and blacklist_trees[c].overlap(int(s), int(e)))
                for c, s, e in zip(cnv.Chromosome, cnv.Start, cnv.End)]
    cnv = cnv[~pd.Series(artifact)].reset_index(drop=True)
    n_masked = len(cnv)

    cnv['size'] = cnv['End'] - cnv['Start']
    cnv = cnv[(cnv['size'] >= CNV_SIZE_MIN) & (cnv['size'] <= CNV_SIZE_MAX)].reset_index(drop=True)
    n_final = len(cnv)

    return cnv, n_raw, n_masked, n_final


# ============================================================
# ★ 고정 샘플셋 결정 (union 기준, 한 번만) — y 컬럼 의존성 제거
# ============================================================
def determine_fixed_sample_set(pheno, union_cnv):
    per_sample = union_cnv.groupby('IID').size().rename('total_cnv')
    tmp = pheno[[]].join(per_sample, how='left').fillna({'total_cnv': 0})
    mu, sdv = tmp['total_cnv'].mean(), tmp['total_cnv'].std()
    cutoff = mu + SAMPLE_SD_CUT * sdv
    kept = tmp[tmp['total_cnv'] <= cutoff]
    excluded = tmp[tmp['total_cnv'] > cutoff]
    print(f"[고정 샘플셋] union 기준 outlier cutoff={cutoff:.1f}개 -> "
          f"{len(kept)}명 유지 / {len(excluded)}명 제외")
    if len(excluded):
        print(f"   제외된 샘플 예시(최대 5개): {list(excluded.index[:5])}")
    return kept.index, tmp


# ============================================================
# 샘플별 burden 계산: CNV 개수 + 총 길이(bp), log1p 변환까지
# ============================================================
def compute_burden(cnv, fixed_index):
    per_count  = cnv.groupby('IID').size().rename('total_cnv')
    per_length = cnv.groupby('IID')['size'].sum().rename('total_length')

    burden = pd.concat([per_count, per_length], axis=1).reindex(fixed_index).fillna(0)
    burden['log_total_cnv']    = np.log1p(burden['total_cnv'])
    burden['log_total_length'] = np.log1p(burden['total_length'])
    return burden


# ============================================================
# permutation test (OLS slope 기반, outcome 고정 / predictor 셔플)
# 원본 diag4b의 "analytic vs permutation p-value 비교" 철학을 그대로 유지
# ============================================================
def permutation_slope_pvalue(x, y, n_perm=N_PERM, seed=RANDOM_SEED):
    rng = np.random.default_rng(seed)
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    obs_slope = np.polyfit(x, y, 1)[0]

    n = len(y)
    perm_slopes = np.empty(n_perm)
    for i in range(n_perm):
        y_shuf = y[rng.permutation(n)]
        perm_slopes[i] = np.polyfit(x, y_shuf, 1)[0]

    p = (np.abs(perm_slopes) >= np.abs(obs_slope)).mean()
    p = max(p, 1.0 / (n_perm + 1))
    return p, obs_slope


# ============================================================
# 표준 OLS: outcome ~ TREATMENT, HC3 robust SE (논문 표준 방식)
# ============================================================
def standard_ols(x, y):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)

    X = sm.add_constant(x)
    model = sm.OLS(y, X).fit(cov_type="HC3")

    beta = model.params[1]
    se   = model.bse[1]
    analytic_p = model.pvalues[1]
    ci_lo, ci_hi = model.conf_int(alpha=0.05)[1]
    r2 = model.rsquared

    perm_p, _ = permutation_slope_pvalue(x, y)

    return {
        "beta": beta, "SE": se,
        "CI_low": ci_lo, "CI_high": ci_hi,
        "analytic_p": analytic_p,
        "perm_p": perm_p,
        "r_squared": r2,
    }


# ============================================================
# 하나의 threshold를 "고정 샘플셋" 기준으로 OLS 요약
# ============================================================
def summarize_fixed(label, bed_path, blacklist_trees, fixed_index, treatment_x):
    print(f"\n--- {label} ---")
    if not os.path.exists(bed_path):
        print("   [건너뜀] 파일 없음")
        return None

    cnv, n_raw, n_masked, n_final = load_masked_filtered_cnv(bed_path, blacklist_trees)
    burden = compute_burden(cnv, fixed_index)

    trt = treatment_x.reindex(fixed_index)
    valid = trt.notna()
    x = trt[valid].values

    row = {
        "threshold": label,
        "n_raw": n_raw, "n_after_masking": n_masked, "n_final": n_final,
        "n_used": int(valid.sum()),
    }

    for outcome_col, prefix in [("log_total_cnv", "cnvcount"), ("log_total_length", "cnvlen")]:
        y = burden.loc[valid, outcome_col].values
        res = standard_ols(x, y)
        for k, v in res.items():
            row[f"{prefix}_{k}"] = v

    print(f"   n_final={n_final} n_used={row['n_used']} | "
          f"[log CNV count] beta={row['cnvcount_beta']:.4f} p={row['cnvcount_analytic_p']:.3g} "
          f"perm_p={row['cnvcount_perm_p']:.3g} R2={row['cnvcount_r_squared']:.3f} | "
          f"[log CNV length] beta={row['cnvlen_beta']:.4f} p={row['cnvlen_analytic_p']:.3g} "
          f"perm_p={row['cnvlen_perm_p']:.3g} R2={row['cnvlen_r_squared']:.3f}")
    return row


def plot_comparison(summary_df, out_dir):
    df = summary_df.dropna(subset=["n_final"]).reset_index(drop=True)
    x = np.arange(len(df))
    labels = df['threshold'].tolist()

    fig, axes = plt.subplots(4, 1, figsize=(9, 14), sharex=True)

    axes[0].bar(x, df['n_final'], color='#2c3e7b')
    axes[0].set_ylabel('최종 CNV 개수')
    axes[0].set_title('Threshold별 최종 CNV 개수')

    width = 0.35
    axes[1].bar(x - width/2, df['cnvcount_beta'], width, label='log(CNV count) beta', color='#c0392b')
    axes[1].bar(x + width/2, df['cnvlen_beta'], width, label='log(CNV length) beta', color='#2c7bc0')
    axes[1].axhline(0, color='grey', lw=1)
    axes[1].set_ylabel('OLS beta (vs TREATMENT)')
    axes[1].set_title('Threshold별 CNV burden ~ TREATMENT 회귀계수 (고정 샘플셋)')
    axes[1].legend()

    axes[2].plot(x, df['cnvcount_r_squared'], marker='o', label='log(CNV count) R²', color='#c0392b')
    axes[2].plot(x, df['cnvlen_r_squared'], marker='o', label='log(CNV length) R²', color='#2c7bc0')
    axes[2].set_ylabel('R²')
    axes[2].set_title('Threshold별 모델 설명력(R²)')
    axes[2].legend()

    axes[3].plot(x, -np.log10(df['cnvcount_analytic_p'].clip(lower=1e-300)),
                 marker='o', label='count: OLS(HC3) -log10(p)', color='#7b2c8c')
    axes[3].plot(x, -np.log10(df['cnvcount_perm_p'].clip(lower=1e-300)),
                 marker='s', label='count: Permutation -log10(p)', color='#2c8c5a', ls='--')
    axes[3].plot(x, -np.log10(df['cnvlen_analytic_p'].clip(lower=1e-300)),
                 marker='o', label='length: OLS(HC3) -log10(p)', color='#c08c2c')
    axes[3].plot(x, -np.log10(df['cnvlen_perm_p'].clip(lower=1e-300)),
                 marker='s', label='length: Permutation -log10(p)', color='#2c5a8c', ls='--')
    axes[3].axhline(-np.log10(0.05), color='grey', ls=':', lw=1, label='p=0.05')
    axes[3].set_ylabel('-log10(p-value)')
    axes[3].set_title('Threshold별 CNV burden ~ TREATMENT 연관성 유의성\n(OLS analytic vs Permutation 비교)')
    axes[3].legend(fontsize=8)

    plt.xticks(x, labels, rotation=30, ha='right')
    plt.tight_layout()
    plt.savefig(out_dir + "diag5_treatment_ols_comparison.png", dpi=200)
    plt.close()
    print(f"\n저장: diag5_treatment_ols_comparison.png")


if __name__ == "__main__":
    os.makedirs(OUT_DIR, exist_ok=True)
    pheno = load_phenotype()
    blacklist_trees = build_blacklist_tree()

    # 1) union 데이터로 고정 샘플셋을 먼저 확정
    print("=> 고정 샘플셋 산출을 위해 union 데이터 로드...")
    union_cnv, _, _, _ = load_masked_filtered_cnv(THRESHOLD_FILES[UNION_LABEL], blacklist_trees)
    fixed_index, _ = determine_fixed_sample_set(pheno, union_cnv)
    treatment_x = pheno['treatment_x']

    # 2) 모든 threshold를 동일한 fixed_index로 비교
    rows = []
    for label, path in THRESHOLD_FILES.items():
        r = summarize_fixed(label, path, blacklist_trees, fixed_index, treatment_x)
        if r is not None:
            rows.append(r)

    summary_df = pd.DataFrame(rows)
    print("\n===== 고정 샘플셋 기준 Threshold별 TREATMENT OLS 회귀 결과 =====")
    print(summary_df.to_string(index=False, float_format=lambda x: f"{x:.4g}"))
    summary_df.to_csv(OUT_DIR + "diag5_treatment_ols_summary.csv", index=False)

    plot_comparison(summary_df, OUT_DIR)

    print("\n[해석 가이드]")
    print("- 모든 threshold가 같은 사람들을 비교하므로, n_used는 threshold와")
    print("  무관하게 동일해야 합니다 (다르면 fixed_index 적용 버그).")
    print("- beta의 부호: REVERSE_TREATMENT=False(기본값)이면 TREATMENT가 클수록")
    print("  '약 변경 없음(양호)'이므로, beta<0 은 'CNV burden이 높을수록 약 변경이")
    print("  잦다(중증)'는 방향과 일치합니다. REVERSE_TREATMENT=True로 바꾸면 부호가")
    print("  반대로 나오니 원하는 해석 방향에 맞게 설정하세요.")
    print("- analytic_p(OLS+HC3)와 perm_p가 크게 다르면, 정규성/등분산 가정이")
    print("  깨졌다는 뜻이니 perm_p를 더 신뢰하세요.")
    print("- log(CNV count) 결과와 log(CNV length) 결과가 방향은 같은데 유의성이")
    print("  다르면, '개수는 늘었지만 크기는 그대로(또는 반대)'라는 의미이니")
    print("  두 결과를 함께 리포트하는 것을 권장합니다.")
    print("\n=== 진단 5 (TREATMENT OLS) 완료 ===")
