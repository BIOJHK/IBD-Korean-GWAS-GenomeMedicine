#CNV 파이프라인 진단 스크립트 (4c) — depth-adjusted burden 비교
#================================================================
#depth_check_case_control.py 에서 확인된 case/control 시퀀싱 depth 차이
#(case가 약 12% 높음, p=4.7e-30)를 covariate로 넣고, 각 threshold에서
#case-control burden 차이가 depth 보정 후에도 유의하게 남는지 확인합니다.
#
#total_cnv는 카운트 데이터(0,1,2,3...)이므로 정규분포를 가정하는 t-test/
#선형회귀 대신 Poisson GLM을 사용합니다:
#
#    total_cnv ~ group(case/ctrl) + depth_z
#
#group의 계수(및 p-value)가 depth를 통제한 뒤에도 유의하다면, 그건 depth로
#설명되지 않는 잔여 신호이므로 좀 더 신뢰할 수 있습니다.


import warnings
warnings.filterwarnings("ignore")

import os
import pandas as pd
import numpy as np
import statsmodels.api as sm
from intervaltree import IntervalTree

# ============================================================
# CONFIG
# ============================================================
PHENOTYPE_PATH  = "260605_survivor_check/CDKDC.250120.03.pca.pheno.txt"
CENTROMERE_PATH = "260605_survivor_check/centromeres.txt.gz"
SEGDUP_PATH     = "260605_survivor_check/genomicSuperDups.txt.gz"
CYTOBAND_PATH   = "260605_survivor_check/cytoBand.txt.gz"
OUT_DIR         = "260605_survivor_check/diagnostics/"

# depth_check_case_control.py 실행 결과 (IID, y, mean_depth, n_segments 포함)
DEPTH_CSV_PATH  = OUT_DIR + "depth_check_per_sample.csv"

PERICENTRO_BUFFER = 5_000_000
SAMPLE_SD_CUT     = 3
CNV_SIZE_MIN      = 1_000
CNV_SIZE_MAX      = 3_000_000
ACROCENTRIC       = {"Chr 13", "Chr 14", "Chr 15", "Chr 21", "Chr 22"}

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


def load_masked_filtered_cnv(bed_path, blacklist_trees):
    cnv = pd.read_csv(bed_path, sep="\t", header=None,
                       names=["Chromosome", "Start", "End", "Type", "Length", "Sample_Raw", "BED_Status"])
    cnv['Chromosome'] = cnv['Chromosome'].apply(norm)
    cnv['IID'] = (cnv['Sample_Raw'].str.replace('.survivor.chr.vcf', '', regex=False)
                                   .str.replace('_filtered.vcf', '', regex=False))
    cnv = cnv[~cnv['Chromosome'].isin(['Chr X', 'Chr Y'])].reset_index(drop=True)

    artifact = [bool(c in blacklist_trees and blacklist_trees[c].overlap(int(s), int(e)))
                for c, s, e in zip(cnv.Chromosome, cnv.Start, cnv.End)]
    cnv = cnv[~pd.Series(artifact)].reset_index(drop=True)

    cnv['size'] = cnv['End'] - cnv['Start']
    cnv = cnv[(cnv['size'] >= CNV_SIZE_MIN) & (cnv['size'] <= CNV_SIZE_MAX)].reset_index(drop=True)
    return cnv


def determine_fixed_sample_set(depth_df):
    tmp = depth_df.copy()
    tmp['total_cnv'] = tmp.get('total_cnv', 0)  # placeholder, union에서 채움
    return tmp.index


def poisson_group_test(total_cnv, y, depth_z):
    """total_cnv ~ group + depth_z 를 Poisson GLM으로 적합, group 계수/p 반환."""
    X = np.column_stack([np.ones(len(y)), y, depth_z])
    try:
        model = sm.GLM(total_cnv, X, family=sm.families.Poisson())
        res = model.fit()
        coef, se, p = res.params[1], res.bse[1], res.pvalues[1]
        # 과산포(overdispersion) 체크: Pearson chi2 / df_resid 가 1보다 많이 크면
        # Poisson 대신 Negative Binomial을 고려해야 함
        dispersion = res.pearson_chi2 / res.df_resid
        return coef, se, p, dispersion
    except Exception as e:
        return np.nan, np.nan, np.nan, np.nan


def poisson_group_only_test(total_cnv, y):
    """depth 보정 없는 버전 (비교용, MWU 대신 count 모델로 일관성 있게)."""
    X = np.column_stack([np.ones(len(y)), y])
    try:
        model = sm.GLM(total_cnv, X, family=sm.families.Poisson())
        res = model.fit()
        return res.params[1], res.pvalues[1]
    except Exception:
        return np.nan, np.nan


if __name__ == "__main__":
    os.makedirs(OUT_DIR, exist_ok=True)

    depth_df = pd.read_csv(DEPTH_CSV_PATH).set_index("IID")
    print(f"depth 데이터 로드: {len(depth_df)}명")

    blacklist_trees = build_blacklist_tree()

    # union 기준 고정 샘플셋 (depth 정보가 있는 사람만)
    union_cnv = load_masked_filtered_cnv(THRESHOLD_FILES[UNION_LABEL], blacklist_trees)
    per_sample_union = union_cnv.groupby('IID').size().rename('total_cnv')
    tmp = depth_df.join(per_sample_union, how='left').fillna({'total_cnv': 0})
    mu, sdv = tmp['total_cnv'].mean(), tmp['total_cnv'].std()
    cutoff = mu + SAMPLE_SD_CUT * sdv
    fixed_index = tmp[tmp['total_cnv'] <= cutoff].index
    print(f"고정 샘플셋: {len(fixed_index)}명 (union 기준 outlier {(tmp['total_cnv'] > cutoff).sum()}명 제외)")

    y_all = depth_df.loc[fixed_index, 'y'].values.astype(float)
    depth_all = depth_df.loc[fixed_index, 'mean_depth'].values.astype(float)
    depth_z = (depth_all - depth_all.mean()) / depth_all.std()

    rows = []
    for label, path in THRESHOLD_FILES.items():
        cnv = load_masked_filtered_cnv(path, blacklist_trees)
        per_sample = cnv.groupby('IID').size().rename('total_cnv')
        total_cnv = per_sample.reindex(fixed_index).fillna(0).values.astype(float)

        coef_raw, p_raw = poisson_group_only_test(total_cnv, y_all)
        coef_adj, se_adj, p_adj, disp = poisson_group_test(total_cnv, y_all, depth_z)

        case_mean = total_cnv[y_all == 1].mean()
        ctrl_mean = total_cnv[y_all == 0].mean()

        row = {
            "threshold": label,
            "case_mean": case_mean, "ctrl_mean": ctrl_mean,
            "Poisson_coef_unadjusted": coef_raw, "Poisson_p_unadjusted": p_raw,
            "Poisson_coef_depth_adjusted": coef_adj, "Poisson_p_depth_adjusted": p_adj,
            "overdispersion(pearson_chi2/df)": disp,
        }
        rows.append(row)
        print(f"\n--- {label} ---")
        print(f"   case_mean={case_mean:.3f}  ctrl_mean={ctrl_mean:.3f}")
        print(f"   depth 보정 전: coef={coef_raw:+.3f}  p={p_raw:.3g}")
        print(f"   depth 보정 후: coef={coef_adj:+.3f}  p={p_adj:.3g}  "
              f"(overdispersion={disp:.2f})")
        if disp > 2:
            print(f"   [참고] overdispersion이 {disp:.1f}로 높습니다 — Poisson 대신 "
                  f"Negative Binomial 모델을 쓰면 더 정확할 수 있습니다.")

    summary_df = pd.DataFrame(rows)
    print("\n===== depth 보정 전/후 비교 =====")
    print(summary_df.to_string(index=False, float_format=lambda x: f"{x:.4g}"))
    summary_df.to_csv(OUT_DIR + "diag4c_depth_adjusted_summary.csv", index=False)

    print("\n[해석 가이드]")
    print("- Poisson_p_unadjusted 와 Poisson_p_depth_adjusted 를 나란히 비교하세요.")
    print("  depth 보정 후 p-value가 크게 완화(예: p<0.001 -> p>0.05)된 threshold는")
    print("  그 신호의 상당 부분이 depth 차이로 설명된다는 뜻입니다.")
    print("- 반대로 depth 보정 후에도 유의성이 거의 그대로 유지된다면, depth로")
    print("  설명되지 않는 잔여 신호가 있다는 뜻이라 더 신뢰할 수 있습니다.")
    print("- Poisson_coef의 부호가 보정 전후로 바뀐다면(예: + -> -), 이건 특히")
    print("  주의 깊게 보셔야 합니다 — depth가 신호의 방향 자체를 왜곡하고")
    print("  있었다는 뜻입니다.")
