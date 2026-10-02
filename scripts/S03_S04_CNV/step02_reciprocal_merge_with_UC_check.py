
#CNV 파이프라인 진단 스크립트 (4b) — 고정 샘플셋 + mean/pct_zero + permutation test
#====================================================================================
#이전 diag4의 한계 3가지를 보완합니다.

#  1) threshold마다 샘플 outlier 제거 기준(SAMPLE_SD_CUT)이 그 threshold 자체의
#     total_cnv 분포에 의존 -> 비교 대상 인원이 threshold마다 미묘하게 달랐음
#     => union 데이터 기준으로 outlier를 한 번만 정하고, 그 샘플 집합을
#        모든 threshold에 동일하게 고정합니다.

#  2) median만 보면 0/1/2 사이를 오가는 이산적 상황에서 정보 손실이 큼
#     => mean, median, "0개 비율(pct_zero)"을 case/control 각각 다 계산합니다.

#  3) Mann-Whitney U는 tie(동점)가 많고 zero-inflated인 데이터에서 근사가
#     불안정할 수 있음
#     => label permutation(순열검정)으로 mean 차이에 대한 경험적 p-value를
#        별도로 계산해서 MWU와 나란히 비교합니다.

#실행하면 threshold별 요약 테이블(csv)과 비교 플롯(png)이 저장됩니다.


import warnings
warnings.filterwarnings("ignore")

import os
import pandas as pd
import numpy as np
import scipy.stats as stats
from intervaltree import IntervalTree
import matplotlib.pyplot as plt

# ============================================================
# CONFIG
# ============================================================
PHENOTYPE_PATH  = "260605_survivor_check/UC526KDC.pheno.txt"
CENTROMERE_PATH = "260605_survivor_check/centromeres.txt.gz"
SEGDUP_PATH     = "260605_survivor_check/genomicSuperDups.txt.gz"
CYTOBAND_PATH   = "260605_survivor_check/cytoBand.txt.gz"
OUT_DIR         = "260605_survivor_check/diagnostics/"

PERICENTRO_BUFFER = 5_000_000
SAMPLE_SD_CUT     = 3            # outlier 판정은 union 기준으로 딱 한 번만 적용
CNV_SIZE_MIN      = 1_000
CNV_SIZE_MAX      = 3_000_000
ACROCENTRIC       = {"Chr 13", "Chr 14", "Chr 15", "Chr 21", "Chr 22"}

N_PERM       = 5000
RANDOM_SEED  = 42

# threshold(라벨) : 병합 BED 파일 경로. UNION_LABEL로 지정한 항목을
# "고정 샘플셋을 정하는 기준"으로 사용합니다 (가장 정보량이 많아서
# outlier를 가장 잘 드러내는 데이터라고 보기 때문입니다).
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
# 블랙리스트 (모든 threshold 공통)
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


def load_phenotype():
    pheno = pd.read_csv(PHENOTYPE_PATH, sep="\t").set_index("IID")
    if pheno['DIS'].dtype == object:
        pheno['y'] = np.where(pheno['DIS'].str.lower() == 'case', 1, 0)
    else:
        pheno['y'] = pheno['DIS'] - 1 if pheno['DIS'].max() == 2 else pheno['DIS']
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
# ★ 고정 샘플셋 결정 (union 기준, 한 번만)
# ============================================================
def determine_fixed_sample_set(pheno, union_cnv):
    per_sample = union_cnv.groupby('IID').size().rename('total_cnv')
    tmp = pheno[['y']].join(per_sample, how='left').fillna({'total_cnv': 0})
    mu, sdv = tmp['total_cnv'].mean(), tmp['total_cnv'].std()
    cutoff = mu + SAMPLE_SD_CUT * sdv
    kept = tmp[tmp['total_cnv'] <= cutoff]
    excluded = tmp[tmp['total_cnv'] > cutoff]
    print(f"[고정 샘플셋] union 기준 outlier cutoff={cutoff:.1f}개 -> "
          f"{len(kept)}명 유지 / {len(excluded)}명 제외")
    if len(excluded):
        print(f"   제외된 샘플 예시(최대 5개): {list(excluded.index[:5])}")
    return kept.index, tmp  # 고정 index(IID 리스트)


# ============================================================
# permutation test (mean 차이 기반, label shuffling)
# ============================================================
def permutation_pvalue(case_v, ctrl_v, n_perm=N_PERM, seed=RANDOM_SEED):
    rng = np.random.default_rng(seed)
    case_v = np.asarray(case_v, dtype=float)
    ctrl_v = np.asarray(ctrl_v, dtype=float)
    obs_diff = case_v.mean() - ctrl_v.mean()

    pooled = np.concatenate([case_v, ctrl_v])
    n_case = len(case_v)
    n_total = len(pooled)

    diffs = np.empty(n_perm)
    for i in range(n_perm):
        idx = rng.permutation(n_total)
        shuffled = pooled[idx]
        diffs[i] = shuffled[:n_case].mean() - shuffled[n_case:].mean()

    # 양측검정: 관측된 절대차이보다 극단적인 permutation 비율
    p = (np.abs(diffs) >= np.abs(obs_diff)).mean()
    # 0이 나오면 "1/n_perm보다 작다"는 의미이므로 최소값으로 대체 (0으로 보고하지 않기 위함)
    p = max(p, 1.0 / (n_perm + 1))
    return p, obs_diff


# ============================================================
# 하나의 threshold를 "고정 샘플셋" 기준으로 요약
# ============================================================
def summarize_fixed(label, bed_path, blacklist_trees, fixed_index, y_series):
    print(f"\n--- {label} ---")
    if not os.path.exists(bed_path):
        print("   [건너뜀] 파일 없음")
        return None

    cnv, n_raw, n_masked, n_final = load_masked_filtered_cnv(bed_path, blacklist_trees)

    per_sample = cnv.groupby('IID').size().rename('total_cnv')
    # ★ 핵심: 항상 동일한 fixed_index로 reindex -> 모든 threshold에서 분모가 완전히 동일
    total_cnv = per_sample.reindex(fixed_index).fillna(0)

    y = y_series.reindex(fixed_index)
    case_v = total_cnv[y == 1]
    ctrl_v = total_cnv[y == 0]

    mwu_p = stats.mannwhitneyu(case_v, ctrl_v)[1] if len(case_v) and len(ctrl_v) else np.nan
    perm_p, obs_diff = permutation_pvalue(case_v.values, ctrl_v.values)

    pct_case_zero = (case_v == 0).mean() * 100
    pct_ctrl_zero = (ctrl_v == 0).mean() * 100

    row = {
        "threshold": label,
        "n_raw": n_raw, "n_after_masking": n_masked, "n_final": n_final,
        "n_case": len(case_v), "n_ctrl": len(ctrl_v),
        "case_mean": case_v.mean(), "ctrl_mean": ctrl_v.mean(),
        "case_median": case_v.median(), "ctrl_median": ctrl_v.median(),
        "pct_case_zero": pct_case_zero, "pct_ctrl_zero": pct_ctrl_zero,
        "obs_mean_diff(case-ctrl)": obs_diff,
        "MWU_pvalue": mwu_p,
        "perm_pvalue": perm_p,
    }
    print(f"   n_final={n_final} | case_mean={row['case_mean']:.3f}(0%:{pct_case_zero:.0f}%) "
          f"ctrl_mean={row['ctrl_mean']:.3f}(0%:{pct_ctrl_zero:.0f}%) | "
          f"MWU_p={mwu_p:.3g}  perm_p={perm_p:.3g}")
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
    axes[1].bar(x - width/2, df['case_mean'], width, label='case mean', color='#c0392b')
    axes[1].bar(x + width/2, df['ctrl_mean'], width, label='control mean', color='#2c7bc0')
    axes[1].set_ylabel('샘플당 CNV 평균')
    axes[1].set_title('Threshold별 case/control burden 평균 (고정 샘플셋)')
    axes[1].legend()

    axes[2].plot(x, df['pct_case_zero'], marker='o', label='case', color='#c0392b')
    axes[2].plot(x, df['pct_ctrl_zero'], marker='o', label='control', color='#2c7bc0')
    axes[2].set_ylabel('CNV 0개인 샘플 비율(%)')
    axes[2].set_title('Threshold별 "정보 없음(0개)" 샘플 비율')
    axes[2].legend()

    axes[3].plot(x, -np.log10(df['MWU_pvalue'].clip(lower=1e-300)),
                 marker='o', label='MWU -log10(p)', color='#7b2c8c')
    axes[3].plot(x, -np.log10(df['perm_pvalue'].clip(lower=1e-300)),
                 marker='s', label='Permutation -log10(p)', color='#2c8c5a', ls='--')
    axes[3].axhline(-np.log10(0.05), color='grey', ls=':', lw=1, label='p=0.05')
    axes[3].set_ylabel('-log10(p-value)')
    axes[3].set_title('Threshold별 case-control burden 불균형 유의성\n(MWU vs Permutation 비교)')
    axes[3].legend()

    plt.xticks(x, labels, rotation=30, ha='right')
    plt.tight_layout()
    plt.savefig(out_dir + "diag4b_fixed_sampleset_comparison.png", dpi=200)
    plt.close()
    print(f"\n저장: diag4b_fixed_sampleset_comparison.png")


if __name__ == "__main__":
    os.makedirs(OUT_DIR, exist_ok=True)
    pheno = load_phenotype()
    blacklist_trees = build_blacklist_tree()

    # 1) union 데이터로 고정 샘플셋을 먼저 확정
    print("=> 고정 샘플셋 산출을 위해 union 데이터 로드...")
    union_cnv, _, _, _ = load_masked_filtered_cnv(THRESHOLD_FILES[UNION_LABEL], blacklist_trees)
    fixed_index, _ = determine_fixed_sample_set(pheno, union_cnv)
    y_series = pheno['y']

    # 2) 모든 threshold를 동일한 fixed_index로 비교
    rows = []
    for label, path in THRESHOLD_FILES.items():
        r = summarize_fixed(label, path, blacklist_trees, fixed_index, y_series)
        if r is not None:
            rows.append(r)

    summary_df = pd.DataFrame(rows)
    print("\n===== 고정 샘플셋 기준 Threshold 비교 =====")
    print(summary_df.to_string(index=False, float_format=lambda x: f"{x:.4g}"))
    summary_df.to_csv(OUT_DIR + "diag4b_fixed_sampleset_summary.csv", index=False)

    plot_comparison(summary_df, OUT_DIR)

    print("\n[해석 가이드]")
    print("- 이제 모든 threshold가 정확히 같은 사람들을 비교하므로, n_case/n_ctrl은")
    print("  threshold와 무관하게 동일해야 합니다 (동일하지 않다면 fixed_index 적용에")
    print("  버그가 있다는 뜻이니 코드를 다시 확인하세요).")
    print("- MWU_p와 perm_p가 서로 크게 다르다면(예: 한쪽만 유의), tie/zero-inflation")
    print("  때문에 MWU 근사가 불안정했다는 뜻이니 perm_p를 더 신뢰하세요.")
    print("- pct_case_zero / pct_ctrl_zero 가 threshold가 올라갈수록 둘 다 80~90%")
    print("  이상으로 수렴한다면, 그 지점부터는 '깨끗해진 것'이 아니라 '검정력이")
    print("  거의 없어진 것'일 가능성이 큽니다.")
    print("\n=== 진단 4b 완료 ===")
