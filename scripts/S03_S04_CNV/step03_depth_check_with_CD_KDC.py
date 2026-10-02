#
#Case/Control 시퀀싱 depth 간이 점검
#=====================================
#새 데이터를 만들 필요 없이, 이미 갖고 계신 CNVkit .cns 파일의 'depth'
#컬럼만 읽어서 샘플별 평균 depth를 구하고 case/control로 나눠 비교합니다.
#
#이게 유의하게 다르면(특히 case가 전반적으로 낮다면), 지금까지 본
#union의 case 쪽 콜 인플레이션과 30~40%에서의 방향 반전을
#상당 부분 설명할 수 있는 유력한 근거가 됩니다.
#

import warnings
warnings.filterwarnings("ignore")

import os
import glob
import pandas as pd
import numpy as np
import scipy.stats as stats
import matplotlib.pyplot as plt

# ============================================================
# CONFIG
# ============================================================
CNVKIT_DIRS = [
    "/input/data/CNV_all/CNVKIT/IBD_cns",
    "/input/data/CNV_all/CNVKIT/KDC_cns",
]
PHENOTYPE_PATH = "260605_survivor_check/CDKDC.250120.03.pca.pheno.txt"
OUT_DIR = "260605_survivor_check/diagnostics/"


def cnvkit_sample_id(basename):
    return basename.replace(".cov.cns", "").replace(".cns", "")


def main():
    os.makedirs(OUT_DIR, exist_ok=True)

    pheno = pd.read_csv(PHENOTYPE_PATH, sep=" ").set_index("IID")
    if pheno['DIS'].dtype == object:
        pheno['y'] = np.where(pheno['DIS'].str.lower() == 'case', 1, 0)
    else:
        pheno['y'] = pheno['DIS'] - 1 if pheno['DIS'].max() == 2 else pheno['DIS']

    print("=> .cns 파일에서 샘플별 평균 depth 계산 중...")
    rows = []
    files = []
    for d in CNVKIT_DIRS:
        files.extend(glob.glob(os.path.join(d, "*.cns")))
    print(f"   대상 파일 {len(files)}개")

    for i, f in enumerate(files):
        sid = cnvkit_sample_id(os.path.basename(f))
        try:
            df = pd.read_csv(f, sep="\t", usecols=["depth", "probes"])
        except Exception as e:
            print(f"   [스킵] {f}: {e}")
            continue
        # probes(구간에 포함된 probe 수)로 가중평균 -> 큰 구간이 평균에 더 반영되도록
        if df['probes'].sum() > 0:
            mean_depth = np.average(df['depth'], weights=df['probes'])
        else:
            mean_depth = df['depth'].mean()
        rows.append({"IID": sid, "mean_depth": mean_depth, "n_segments": len(df)})
        if (i + 1) % 500 == 0:
            print(f"   진행률: {i+1}/{len(files)}")

    depth_df = pd.DataFrame(rows).set_index("IID")
    merged = pheno[['y']].join(depth_df, how='inner')
    n_missing = len(pheno) - len(merged)
    print(f"   phenotype-.cns 매칭: {len(merged)}명 (매칭 안 된 {n_missing}명 제외)")

    case_d = merged.loc[merged.y == 1, 'mean_depth'].dropna()
    ctrl_d = merged.loc[merged.y == 0, 'mean_depth'].dropna()

    print("\n===== Case vs Control 평균 depth =====")
    print(f"Case   (n={len(case_d)}): mean={case_d.mean():.2f}  median={case_d.median():.2f}  "
          f"sd={case_d.std():.2f}")
    print(f"Control(n={len(ctrl_d)}): mean={ctrl_d.mean():.2f}  median={ctrl_d.median():.2f}  "
          f"sd={ctrl_d.std():.2f}")

    mwu_stat, mwu_p = stats.mannwhitneyu(case_d, ctrl_d)
    t_stat, t_p = stats.ttest_ind(case_d, ctrl_d, equal_var=False)
    print(f"\nMWU p-value    : {mwu_p:.3e}")
    print(f"Welch t-test p : {t_p:.3e}")
    fold = case_d.mean() / ctrl_d.mean()
    print(f"Case/Control depth 배율: {fold:.3f}")

    merged.to_csv(OUT_DIR + "depth_check_per_sample.csv")

    plt.figure(figsize=(7, 5))
    plt.hist(case_d, bins=40, alpha=0.5, label=f'case (n={len(case_d)})', density=True)
    plt.hist(ctrl_d, bins=40, alpha=0.5, label=f'control (n={len(ctrl_d)})', density=True)
    plt.xlabel('샘플 평균 depth (probe-weighted)')
    plt.ylabel('density')
    plt.title(f'Case vs Control 평균 depth 분포 (MWU p={mwu_p:.2e})')
    plt.legend()
    plt.tight_layout()
    plt.savefig(OUT_DIR + "depth_check_distribution.png", dpi=200)
    plt.close()

    print(f"\n저장: depth_check_per_sample.csv / depth_check_distribution.png")

    print("\n[해석 가이드]")
    if mwu_p < 0.001 and abs(fold - 1) > 0.1:
        print("   [경고] Case/Control 간 depth가 통계적으로 유의하고 10% 이상 차이납니다.")
        print("   지금까지 관찰된 union의 case 쪽 콜 인플레이션, 그리고 30~40%에서의")
        print("   방향 반전이 이 depth 차이만으로 상당 부분 설명될 가능성이 높습니다.")
        print("   -> 가능하다면 depth를 별도 covariate로 넣거나, depth-matched")
        print("      subsampling을 고려하시길 권합니다.")
    elif mwu_p < 0.05:
        print("   [주의] 약한 수준의 depth 차이가 있습니다. 다른 요인(배치, 캡처 키트)과")
        print("   함께 고려해보시되, 이것만으로 전체 반전을 설명하긴 어려울 수 있습니다.")
    else:
        print("   [정상] Case/Control 간 depth 차이가 통계적으로 뚜렷하지 않습니다.")
        print("   depth 자체는 이번 반전의 주된 원인이 아닐 가능성이 높습니다 —")
        print("   다른 배치 요인(캡처 키트, reference panel, 시퀀싱 시점)을 의심해보세요.")


if __name__ == "__main__":
    main()
