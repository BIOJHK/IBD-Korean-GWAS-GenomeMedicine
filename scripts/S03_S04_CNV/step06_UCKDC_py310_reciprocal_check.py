import warnings
warnings.filterwarnings("ignore")

import os
import time
import pandas as pd
import numpy as np
import scipy.stats as stats
from scipy.stats import chi2, fisher_exact
from intervaltree import IntervalTree
import statsmodels.api as sm
import matplotlib.pyplot as plt
from firthlogist import FirthLogisticRegression
from concurrent.futures import ProcessPoolExecutor

# ============================================================
# 1. 경로
# ============================================================
PHENOTYPE_PATH  = "260605_survivor_check/UC526KDC.pheno.txt"
PC_PATH         = "260605_survivor_check/UC526KDC.02.pca30.eigenvec"
CNV_DATA_PATH   = "/input/data/CNV_all/merged_by_threshold/merged_consensus_03.bed"
CENTROMERE_PATH = "260605_survivor_check/centromeres.txt.gz"
SEGDUP_PATH     = "260605_survivor_check/genomicSuperDups.txt.gz"
CYTOBAND_PATH   = "260605_survivor_check/cytoBand.txt.gz"
OUT_DIR         = "260605_survivor_check/UC526KDC_260716.intcheck30"

# ============================================================
# 2. 옵션
# ============================================================
WINDOW_SIZE       = 1_000_000
MIN_CARRIERS      = 10
FREQ_MAX          = 0.05
PERICENTRO_BUFFER = 5_000_000
SAMPLE_SD_CUT     = 3
CNV_SIZE_MIN      = 1_000
CNV_SIZE_MAX      = 3_000_000
N_PC              = 5
ACROCENTRIC       = {"Chr 13","Chr 14","Chr 15","Chr 21","Chr 22"}
N_JOBS            = max(1, (os.cpu_count() or 4) - 2)

BURDEN_MODE = "loo"
SUFFIX = {"none": "PConly", "global": "PCadj_global", "loo": "PCadj_loo"}[BURDEN_MODE]

def norm(c):
    p = str(c).replace('chr','').strip()
    return f"Chr {int(p):02d}" if p.isdigit() else f"Chr {p}"

def genomic_lambda(pvals):
    p = np.asarray(pvals, float); p = p[(p > 0) & (p <= 1)]
    return np.median(chi2.isf(p, df=1)) / chi2.ppf(0.5, df=1)

def load_pcs(path, n_pc):
    raw = pd.read_csv(path, sep=r"\s+", header=None)
    try:
        float(raw.iloc[0, 2]); has_header = False
    except (ValueError, TypeError):
        has_header = True
    if has_header:
        raw = pd.read_csv(path, sep=r"\s+")
        cols = [str(c).replace('#','').strip() for c in raw.columns]; raw.columns = cols
        iid_col = "IID" if "IID" in cols else cols[1]
        pc_cols = [c for c in cols if c.upper().startswith("PC")][:n_pc]
        out = raw.rename(columns={iid_col: "IID"}).set_index("IID")[pc_cols]
    else:
        raw.columns = ["FID","IID"] + [f"PC{i}" for i in range(1, raw.shape[1]-1)]
        out = raw.set_index("IID")[[f"PC{i}" for i in range(1, n_pc+1)]]
    out.columns = [f"PC{i+1}" for i in range(out.shape[1])]
    return out.apply(pd.to_numeric, errors="coerce")

def _zscore(x):
    sd = x.std()
    if sd < 1e-8:
        return np.zeros_like(x)
    return (x - x.mean()) / sd

# ============================================================
# ★ 병렬 처리용 worker 함수/초기화자
# ============================================================
_WD = {}

def _init_worker(mode, sample_idx_list, id_to_pos, y_arr, n_total, pc_mat,
                  fixed_cov_mat, fixed_llf0,
                  del_n_total, dup_n_total, size_sum_total, n_cnv_total,
                  win_burden_lookup,
                  min_carriers, freq_max, window_size):
    _WD['mode'] = mode
    _WD['sample_idx_list'] = sample_idx_list
    _WD['id_to_pos'] = id_to_pos
    _WD['y_arr'] = y_arr
    _WD['n_total'] = n_total
    _WD['pc_mat'] = pc_mat
    _WD['fixed_cov_mat'] = fixed_cov_mat
    _WD['fixed_llf0'] = fixed_llf0
    _WD['del_n_total'] = del_n_total
    _WD['dup_n_total'] = dup_n_total
    _WD['size_sum_total'] = size_sum_total
    _WD['n_cnv_total'] = n_cnv_total
    _WD['win_burden_lookup'] = win_burden_lookup
    _WD['MIN_CARRIERS'] = min_carriers
    _WD['FREQ_MAX'] = freq_max
    _WD['WINDOW_SIZE'] = window_size

def _build_cov_and_llf0(chrom, ws):
    d = _WD
    window_idx = ws // d['WINDOW_SIZE']
    del_loo   = d['del_n_total'].copy()
    dup_loo   = d['dup_n_total'].copy()
    size_loo  = d['size_sum_total'].copy()
    ncnv_loo  = d['n_cnv_total'].copy()

    contrib = d['win_burden_lookup'].get((chrom, window_idx))
    if contrib:
        for iid, (dn, dpn, nw, ssum) in contrib.items():
            pos = d['id_to_pos'].get(iid)
            if pos is not None:
                del_loo[pos]  -= dn
                dup_loo[pos]  -= dpn
                ncnv_loo[pos] -= nw
                size_loo[pos] -= ssum

    size_mean_loo = np.divide(size_loo, ncnv_loo,
                               out=np.zeros_like(size_loo, dtype=float),
                               where=ncnv_loo > 0)

    cov_list = []
    for arr in (del_loo, dup_loo, size_mean_loo):
        z = _zscore(arr.astype(float))
        if np.std(z) > 1e-8:
            cov_list.append(z)
    cov_list.append(d['pc_mat'])
    cov_mat_w = np.column_stack(cov_list) if cov_list else d['pc_mat']

    X0 = np.column_stack([np.ones(d['n_total']), cov_mat_w])
    try:
        res0 = sm.Logit(d['y_arr'], X0).fit(disp=0, maxiter=100)
        if not res0.mle_retvals.get('converged', False):
            return None, None
        return cov_mat_w, res0.llf
    except Exception:
        return None, None

def _test_one_window(task):
    chrom, ws, carriers = task
    d = _WD
    sample_idx_list = d['sample_idx_list']; y_arr = d['y_arr']; n_total = d['n_total']
    MIN_CARRIERS = d['MIN_CARRIERS']; FREQ_MAX = d['FREQ_MAX']; WINDOW_SIZE = d['WINDOW_SIZE']

    carriers_set = set(carriers)
    is_carrier = np.array([iid in carriers_set for iid in sample_idx_list], dtype=float)
    nc = int(is_carrier.sum())

    if nc < MIN_CARRIERS or nc / n_total > FREQ_MAX:
        return None, "freq"
    if nc == n_total or np.std(is_carrier) == 0:
        return None, "collinear"

    n_carrier_case = int(is_carrier[y_arr == 1].sum())
    n_carrier_ctrl = int(is_carrier[y_arr == 0].sum())
    if n_carrier_case == 0 or n_carrier_ctrl == 0:
        return None, "no_info"

    a = n_carrier_case
    b = int((y_arr == 1).sum()) - n_carrier_case
    c = n_carrier_ctrl
    dd = int((y_arr == 0).sum()) - n_carrier_ctrl
    or_crude, p_crude = fisher_exact([[a, b], [c, dd]])

    if d['mode'] == 'loo':
        cov_mat_w, llf0_w = _build_cov_and_llf0(chrom, ws)
        if cov_mat_w is None:
            return None, "not_converged"
    else:
        cov_mat_w = d['fixed_cov_mat']
        llf0_w = d['fixed_llf0']

    X = np.column_stack([np.ones(n_total), is_carrier, cov_mat_w])
    use_firth = (n_carrier_case < 3 or n_carrier_ctrl < 3)

    try:
        if not use_firth:
            res = sm.Logit(y_arr, X).fit(disp=0, maxiter=100)
            if not res.mle_retvals.get('converged', False):
                return None, "not_converged"
            beta = res.params[1]
            if not np.isfinite(beta) or abs(beta) > 8:
                return None, "extreme_beta"
            lrt = 2.0 * (res.llf - llf0_w)
            pv = float(chi2.sf(lrt, df=1)) if lrt > 0 else 1.0
            if not np.isfinite(pv):
                return None, "extreme_beta"
        else:
            Xf = np.column_stack([is_carrier, cov_mat_w])
            fl = FirthLogisticRegression(max_iter=100)
            fl.fit(Xf, y_arr)
            beta = fl.coef_[0]
            pv = fl.pvals_[0]
            if not np.isfinite(beta) or not np.isfinite(pv) or abs(beta) > 8:
                return None, "extreme_beta"
    except Exception:
        return None, "exception"

    result = {"Chromosome": chrom, "Win_Mid": ws + WINDOW_SIZE // 2,
              "N_carrier": nc, "N_case": n_carrier_case, "N_ctrl": n_carrier_ctrl,
              "Method": "Firth" if use_firth else "Logit",
              "OR": float(np.exp(beta)), "P_value": pv,
              "OR_crude": float(or_crude), "P_crude": float(p_crude)}
    return result, ("firth_used" if use_firth else None)


def main():
    _t0 = time.time()
    def tick(label):
        nonlocal _t0
        now = time.time()
        print(f"   [TIME] {label}: {now-_t0:.1f}초")
        _t0 = now

    # ============================================================
    # 3. 표현형 + CNV 로드
    # ============================================================
    print("=> [1/7] 데이터 로드...")
    pheno = pd.read_csv(PHENOTYPE_PATH, sep="\t").set_index("IID")
    if pheno['DIS'].dtype == object:
        pheno['y'] = np.where(pheno['DIS'].str.lower()=='case', 1, 0)
    else:
        pheno['y'] = pheno['DIS'] - 1 if pheno['DIS'].max()==2 else pheno['DIS']

    cnv = pd.read_csv(CNV_DATA_PATH, sep="\t", header=None,
                      names=["Chromosome","Start","End","Type","Length","Sample_Raw","BED_Status"])
    cnv['Chromosome'] = cnv['Chromosome'].apply(norm)
    cnv['IID'] = (cnv['Sample_Raw'].str.replace('.survivor.chr.vcf','',regex=False)
                                   .str.replace('_filtered.vcf','',regex=False))
    cnv = cnv[~cnv['Chromosome'].isin(['Chr X', 'Chr Y'])].reset_index(drop=True)
    print(f"   chrX/Y 제외 후 CNV: {len(cnv)}개")
    tick("1단계 데이터 로드")

    # ============================================================
    # 4. 통합 블랙리스트 (벡터화)
    # ============================================================
    print("=> [2/7] cytoBand/centromere/superDup 통합 마스킹...")
    intervals = []
    cyto = pd.read_csv(CYTOBAND_PATH, sep="\t", header=None,
                       names=["chrom","start","end","band","stain"])
    cyto['chrom'] = cyto['chrom'].apply(norm)
    band_sub = cyto[cyto['stain'].isin(["acen","stalk","gvar"])]
    is_acen = (band_sub['stain'] == 'acen').values
    s_vals = band_sub['start'].astype(int).values.copy()
    e_vals = band_sub['end'].astype(int).values.copy()
    s_vals[is_acen] = np.maximum(0, s_vals[is_acen] - PERICENTRO_BUFFER)
    e_vals[is_acen] = e_vals[is_acen] + PERICENTRO_BUFFER
    intervals.extend(zip(band_sub['chrom'].values, s_vals, e_vals))

    for c in ACROCENTRIC:
        parm = cyto[(cyto.chrom==c) & (cyto.band.astype(str).str.startswith("p"))]
        if len(parm):
            intervals.append((c, 0, int(parm.end.max())))

    cen = pd.read_csv(CENTROMERE_PATH, sep="\t", header=None, usecols=[1,2,3],
                      names=["c","s","e"]); cen['c']=cen['c'].apply(norm)
    s_buf = (cen['s'].astype(int) - PERICENTRO_BUFFER).clip(lower=0)
    e_buf = cen['e'].astype(int) + PERICENTRO_BUFFER
    intervals.extend(zip(cen['c'].values, s_buf.values, e_buf.values))

    sd = pd.read_csv(SEGDUP_PATH, sep="\t", header=None, usecols=[1,2,3],
                     names=["c","s","e"]); sd['c']=sd['c'].apply(norm)
    intervals.extend(zip(sd['c'].values, sd['s'].astype(int).values, sd['e'].astype(int).values))

    trees = {}
    for c, s, e in intervals:
        if s < e:
            trees.setdefault(c, IntervalTree()).addi(int(s), int(e))

    artifact = [bool(c in trees and trees[c].overlap(int(s), int(e)))
                for c, s, e in zip(cnv.Chromosome, cnv.Start, cnv.End)]
    cnv = cnv[~pd.Series(artifact)].reset_index(drop=True)
    print(f"   통합 마스킹 후 CNV: {len(cnv)}개 (제거 {sum(artifact)}개)")
    tick("2단계 마스킹")

    # ============================================================
    # 5. CNV 크기 필터
    # ============================================================
    cnv['size'] = cnv['End'] - cnv['Start']
    cnv = cnv[(cnv['size'] >= CNV_SIZE_MIN) & (cnv['size'] <= CNV_SIZE_MAX)].reset_index(drop=True)
    print(f"   크기 필터 후 CNV: {len(cnv)}개")

    # ============================================================
    # 6. 샘플 burden QC + PC 병합
    # ============================================================
    print("=> [3/7] 샘플 QC 및 PC 병합...")
    per_sample = cnv.groupby('IID').size().rename('total_cnv')
    master = pheno[['y']].join(per_sample, how='left').fillna({'total_cnv':0})

    mu, sdv = master['total_cnv'].mean(), master['total_cnv'].std()
    master = master[master['total_cnv'] <= mu + SAMPLE_SD_CUT*sdv]

    cb, kb = master.loc[master.y==1,'total_cnv'], master.loc[master.y==0,'total_cnv']
    print(f"   burden  case중앙값={cb.median():.1f} / ctrl={kb.median():.1f} / "
          f"MWU p={stats.mannwhitneyu(cb,kb)[1]:.2e}")

    pcs = load_pcs(PC_PATH, N_PC)
    master = master.join(pcs, how='left')
    pc_cols = [f"PC{i}" for i in range(1, N_PC+1)]
    n_before = len(master)
    master = master.dropna(subset=pc_cols)
    print(f"   PC 매칭: {len(master)}/{n_before}명 (PC 결측 {n_before-len(master)}명 제외)")

    cnv = cnv[cnv.IID.isin(master.index)]
    assert len(master) > 0, "샘플 0명 — IID 매칭 확인"

    fixed_use_cov = []
    if BURDEN_MODE == "global":
        typ = cnv.assign(_t=cnv['Type'].astype(str).str.upper())
        cand = pd.DataFrame(index=master.index)
        cand['del_n']     = typ[typ._t.str.contains('DEL')].groupby('IID').size()
        cand['dup_n']     = typ[typ._t.str.contains('DUP')].groupby('IID').size()
        cand['size_mean'] = cnv.groupby('IID')['size'].mean()
        cand = cand.fillna(0)
        for c in ['del_n','dup_n','size_mean']:
            if cand[c].std(ddof=0) > 0:
                master[c] = (cand[c]-cand[c].mean())/cand[c].std()
                fixed_use_cov.append(c)

    print(f"   BURDEN_MODE={BURDEN_MODE} / 최종 분석 샘플 {len(master)}명")
    tick("3단계 샘플 QC")

    # ============================================================
    # ★ 진단: 크기별 burden 비교 (필요시 유지, 결과 확인용)
    # ============================================================
    size_bins = [1000, 2000, 5000, 10000, 20000, 50000, 100000, 300000, 3_000_000]
    cnv['size_bin'] = pd.cut(cnv['size'], bins=size_bins)

    case_iid = master.index[master['y']==1]
    ctrl_iid = master.index[master['y']==0]

    case_dist = cnv[cnv.IID.isin(case_iid)].groupby('size_bin', observed=True).size() / len(case_iid)
    ctrl_dist = cnv[cnv.IID.isin(ctrl_iid)].groupby('size_bin', observed=True).size() / len(ctrl_iid)

    comp = pd.DataFrame({'case_per_sample': case_dist, 'ctrl_per_sample': ctrl_dist})
    comp['fold'] = comp['case_per_sample'] / comp['ctrl_per_sample']
    print(comp)

    print("Case  DEL:DUP =")
    print(cnv[cnv.IID.isin(case_iid)]['Type'].str.upper().str.contains('DEL')
          .map({True: 'DEL', False: 'DUP'}).value_counts(normalize=True))
    print("Ctrl  DEL:DUP =")
    print(cnv[cnv.IID.isin(ctrl_iid)]['Type'].str.upper().str.contains('DEL')
          .map({True: 'DEL', False: 'DUP'}).value_counts(normalize=True))

    plt.figure(figsize=(8,5))
    plt.hist(master.loc[case_iid,'total_cnv'], bins=50, alpha=0.5, label='case', density=True)
    plt.hist(master.loc[ctrl_iid,'total_cnv'], bins=50, alpha=0.5, label='control', density=True)
    plt.legend(); plt.xlabel('total CNV per sample'); plt.ylabel('density')
    plt.title('Sample-level CNV burden distribution')
    plt.tight_layout()
    plt.savefig(OUT_DIR + f"burden_distribution_check_{SUFFIX}.png", dpi=200)
    plt.close()   # ★ show() 대신 close() — 터미널 실행 시 멈춤 방지
    cnv = cnv.drop(columns=['size_bin'])
    tick("진단: 크기별 burden 비교")

    # ============================================================
    # 6.5 CNV frequency profiling (TypeNorm 여기서 생성됨)
    # ============================================================
    print("=> [3.5/7] CNV frequency 계산 (vectorized)...")
    case_idx = master.index[master['y'] == 1]
    ctrl_idx = master.index[master['y'] == 0]
    n_case, n_ctrl = len(case_idx), len(ctrl_idx)

    cnv['win_start'] = cnv['Start'] // WINDOW_SIZE
    cnv['win_end']   = (cnv['End'] - 1) // WINDOW_SIZE
    cnv['TypeNorm']  = np.where(cnv['Type'].str.upper().str.contains('DEL'), 'DEL', 'DUP')

    n_win   = (cnv['win_end'] - cnv['win_start'] + 1).values
    rep_idx = np.repeat(np.arange(len(cnv)), n_win)
    cum     = np.cumsum(n_win)
    offset  = np.arange(cum[-1]) - np.repeat(cum - n_win, n_win)
    win_col = np.repeat(cnv['win_start'].values, n_win) + offset

    exploded = pd.DataFrame({
        'Chromosome': cnv['Chromosome'].values[rep_idx],
        'Window':     win_col,
        'IID':        cnv['IID'].values[rep_idx],
        'TypeNorm':   cnv['TypeNorm'].values[rep_idx],
    })
    exploded['is_case'] = exploded['IID'].isin(set(case_idx))

    freq_df = (exploded.groupby(['Chromosome','Window','TypeNorm','is_case'])['IID']
               .nunique().unstack('is_case', fill_value=0)
               .rename(columns={True: 'n_case', False: 'n_ctrl'}).reset_index())
    freq_df['Freq_case'] = freq_df['n_case'] / n_case
    freq_df['Freq_ctrl'] = freq_df['n_ctrl'] / n_ctrl
    freq_df['Win_Mid']   = freq_df['Window'] * WINDOW_SIZE + WINDOW_SIZE // 2
    freq_df.to_csv(OUT_DIR + f"cnv_frequency_by_window_{SUFFIX}.csv", index=False)
    print(f"   frequency window {len(freq_df)}개 계산 완료")
    tick("3.5단계 frequency 계산")

    # ============================================================
    # 6.6 LOO 모드 전용: raw burden totals + per-window contribution lookup
    # ============================================================
    sample_idx = master.index
    sample_idx_list = list(sample_idx)
    id_to_pos = {iid: i for i, iid in enumerate(sample_idx_list)}
    n_total = len(master)

    del_n_total = dup_n_total = size_sum_total = n_cnv_total = None
    win_burden_lookup = {}

    if BURDEN_MODE == "loo":
        print("=> [3.6/7] LOO burden 준비...")
        is_del_all = (cnv['TypeNorm'] == 'DEL').values

        del_n_total   = cnv[is_del_all].groupby('IID').size().reindex(sample_idx_list).fillna(0).values.astype(float)
        dup_n_total   = cnv[~is_del_all].groupby('IID').size().reindex(sample_idx_list).fillna(0).values.astype(float)
        size_sum_total = cnv.groupby('IID')['size'].sum().reindex(sample_idx_list).fillna(0).values.astype(float)
        n_cnv_total    = cnv.groupby('IID').size().reindex(sample_idx_list).fillna(0).values.astype(float)

        burden_exp = pd.DataFrame({
            'Chromosome': cnv['Chromosome'].values[rep_idx],
            'Window':     win_col,
            'IID':        cnv['IID'].values[rep_idx],
            'is_del':     is_del_all[rep_idx],
            'size':       cnv['size'].values[rep_idx],
        })
        bw = (burden_exp.groupby(['Chromosome','Window','IID'])
              .agg(del_n_win=('is_del','sum'), n_win=('is_del','size'), size_sum_win=('size','sum'))
              .reset_index())
        bw['dup_n_win'] = bw['n_win'] - bw['del_n_win']

        for (chrom_k, window_k), g in bw.groupby(['Chromosome','Window']):
            win_burden_lookup[(chrom_k, int(window_k))] = {
                row.IID: (row.del_n_win, row.dup_n_win, row.n_win, row.size_sum_win)
                for row in g.itertuples()
            }
        print(f"   LOO burden lookup: {len(win_burden_lookup)}개 (chrom,window) 조합")
        tick("3.6단계 LOO burden 준비")

    # ============================================================
    # 7. 윈도우별 logistic (병렬 처리)
    # ============================================================
    print(f"=> [4/7] 윈도우별 logistic 검정 (병렬, N_JOBS={N_JOBS}, BURDEN_MODE={BURDEN_MODE})...")
    y_arr = master['y'].values.astype(float)
    pc_mat = master[pc_cols].values.astype(float)

    fixed_cov_mat = pc_mat
    if BURDEN_MODE == "global" and fixed_use_cov:
        fixed_cov_mat = np.column_stack([master[c].values.astype(float) for c in fixed_use_cov] + [pc_mat])

    X0_fixed = np.column_stack([np.ones(n_total), fixed_cov_mat])
    fixed_llf0 = sm.Logit(y_arr, X0_fixed).fit(disp=0, maxiter=100).llf if BURDEN_MODE != "loo" else None

    # 7-1. 윈도우 후보 수집
    tasks = []
    for chrom, grp in cnv.groupby('Chromosome'):
        ctree = IntervalTree()
        for s, e, iid in zip(grp.Start.values, grp.End.values, grp.IID.values):
            if s < e:
                ctree.addi(int(s), int(e), iid)
        cmax = int(grp.End.max())
        for ws in range(0, cmax, WINDOW_SIZE):
            hits = ctree.overlap(ws, ws + WINDOW_SIZE)
            if hits:
                carriers = frozenset(iv.data for iv in hits)
                tasks.append((chrom, ws, carriers))

    # ★★★ 여기가 핵심 수정 — Bonferroni 분모를 여기서 고정 ★★★
    n_candidate_windows = len(tasks)   # carrier가 1명이라도 있어서 "검정을 시도한" 전체 후보 윈도우 수
    print(f"   병렬 처리 대상 윈도우 (Bonferroni 분모로 고정): {n_candidate_windows}개")
    tick("4단계 윈도우 후보 수집")

    # 7-2. 병렬 검정
    results = []
    fail_reasons = {"freq":0, "collinear":0, "no_info":0, "not_converged":0,
                     "extreme_beta":0, "exception":0, "firth_used":0}

    with ProcessPoolExecutor(
        max_workers=N_JOBS, initializer=_init_worker,
        initargs=(BURDEN_MODE, sample_idx_list, id_to_pos, y_arr, n_total, pc_mat,
                  fixed_cov_mat, fixed_llf0,
                  del_n_total, dup_n_total, size_sum_total, n_cnv_total,
                  win_burden_lookup,
                  MIN_CARRIERS, FREQ_MAX, WINDOW_SIZE)
    ) as executor:
        for i, (res, reason) in enumerate(executor.map(_test_one_window, tasks, chunksize=20)):
            if res is not None:
                results.append(res)
                if reason == "firth_used":
                    fail_reasons["firth_used"] += 1
            elif reason in fail_reasons:
                fail_reasons[reason] += 1
            if (i+1) % 500 == 0:
                print(f"   진행률: {i+1}/{len(tasks)}")

    fails = sum(v for k,v in fail_reasons.items() if k != "firth_used")
    print(f"   검정 {len(results)}개 완료 / 제외 {fails}개  (사유별: {fail_reasons})")
    tick("4단계 병렬 logistic 검정")

    win = pd.DataFrame(results)
    lam = genomic_lambda(win['P_value'])
    chisq = chi2.isf(win['P_value'].clip(lower=1e-300), df=1)
    lam_for_correction = max(lam, 1.0)
    win['P_GC'] = chi2.sf(chisq / lam_for_correction, df=1)

    # ★ 고정된 분모를 CSV에도 명시적으로 기록 (재현성/투명성 목적)
    win['N_candidate_windows_total'] = n_candidate_windows
    win['N_tested_windows'] = len(win)

    win.sort_values("P_value").to_csv(OUT_DIR+f"cnv_window_logit_{SUFFIX}.UC.csv", index=False)
    lam = lam_for_correction

    # ============================================================
    # 8. λ + QQ
    # ============================================================
    print("=> [5/7] inflation(λ) 및 QQ...")
    print(f"   genomic inflation λ = {lam:.3f}  (목표 ≲ 1.10)")
    p_sorted = np.sort(win['P_value'].clip(lower=1e-300).values)
    n = len(p_sorted); expq = (np.arange(1,n+1)-0.5)/n
    plt.figure(figsize=(6,6))
    lim = max(-np.log10(expq).max(), -np.log10(p_sorted).max())
    plt.plot([0,lim],[0,lim], color='crimson', lw=1.5)
    plt.scatter(-np.log10(expq), -np.log10(p_sorted), s=10, color='#2c3e7b')
    plt.title(f"QQ ({SUFFIX}, lambda={lam:.2f})")
    plt.xlabel("Expected -log10(P)"); plt.ylabel("Observed -log10(P)")
    plt.tight_layout(); plt.savefig(OUT_DIR+f"qq_{SUFFIX}.UC.png", dpi=300); plt.close()

    # ============================================================
    # 9. Manhattan plot
    # ============================================================
    print("=> [6/7] Manhattan plot...")
    chrom_order = sorted(win['Chromosome'].unique())
    offset, offsets, xticks = 0, {}, []
    GAP = 20_000_000
    for c in chrom_order:
        offsets[c] = offset
        cmax = win.loc[win.Chromosome==c,'Win_Mid'].max()
        xticks.append((c.replace("Chr ",""), offset + cmax/2))
        offset += cmax + GAP
    win['cum_pos']   = win['Win_Mid'] + win['Chromosome'].map(offsets)
    win['neg_log_p'] = -np.log10(win['P_value'].clip(lower=1e-300))

    # ★★★ 핵심 수정 — len(win) 대신 고정된 n_candidate_windows 사용 ★★★
    bonf = 0.05 / n_candidate_windows
    print(f"   Bonferroni 임계값 = 0.05 / {n_candidate_windows}(후보 윈도우) = {bonf:.3e}"
          f"   [참고: 검정 성공 윈도우는 {len(win)}개]")

    plt.figure(figsize=(16,6))
    colors = ['#2c3e7b','#7b9acc']
    for i,c in enumerate(chrom_order):
        s = win[win.Chromosome==c]
        plt.scatter(s['cum_pos'], s['neg_log_p'], s=14, color=colors[i%2], zorder=3)
    plt.axhline(-np.log10(bonf), color='crimson', ls='--', lw=1.5,
                label=f'Bonferroni (p={bonf:.1e}, N={n_candidate_windows})')
    plt.xticks([x for _,x in xticks], [c for c,_ in xticks], rotation=90, fontsize=8)
    plt.xlabel('Chromosome'); plt.ylabel('-log10(P-value)')
    plt.title(f"CNV Association ({SUFFIX}, lambda={lam:.2f})", fontsize=14, fontweight='bold')
    plt.legend(loc='upper right'); plt.tight_layout()
    plt.savefig(OUT_DIR+f"manhattan_{SUFFIX}.UC.png", dpi=300); plt.close()

    # ============================================================
    # 10. 상위 결과
    # ============================================================
    print("=> [7/7] 완료. 상위 10개 윈도우:")
    print(win.sort_values("P_value").head(10)[
        ["Chromosome","Win_Mid","N_carrier","N_case","N_ctrl","OR","P_value","OR_crude","P_crude"]
    ].to_string(index=False, float_format=lambda x: f"{x:.4g}"))
    print(f"\n저장: cnv_window_logit_{SUFFIX}.UC.csv / qq_{SUFFIX}.UC.png / manhattan_{SUFFIX}.UC.png")

    # ★ 재현성을 위한 메타데이터 기록 (논문 Methods 작성 시 그대로 인용 가능)
    with open(OUT_DIR + f"analysis_metadata_{SUFFIX}.txt", "w", encoding="utf-8") as f:
        f.write(f"CNV_DATA_PATH: {CNV_DATA_PATH}\n")
        f.write(f"BURDEN_MODE: {BURDEN_MODE}\n")
        f.write(f"WINDOW_SIZE: {WINDOW_SIZE}\n")
        f.write(f"MIN_CARRIERS: {MIN_CARRIERS}\n")
        f.write(f"FREQ_MAX: {FREQ_MAX}\n")
        f.write(f"N_candidate_windows (Bonferroni 분모, 고정): {n_candidate_windows}\n")
        f.write(f"N_tested_windows (실제 검정 성공): {len(win)}\n")
        f.write(f"Bonferroni_threshold: {bonf:.6e}\n")
        f.write(f"Genomic_inflation_lambda: {lam:.4f}\n")
        f.write(f"N_samples_final: {n_total}\n")
    print(f"   메타데이터 저장: analysis_metadata_{SUFFIX}.txt")

    # ============================================================
    # CNV frequency plot
    # ============================================================
    print("=> CNV frequency plot 생성...")
    chrom_order = sorted(freq_df['Chromosome'].unique())
    offset, offsets, xticks = 0, {}, []
    for c in chrom_order:
        offsets[c] = offset
        cmax = freq_df.loc[freq_df.Chromosome == c, 'Win_Mid'].max()
        xticks.append((c.replace("Chr ", ""), offset + cmax / 2))
        offset += cmax + GAP
    freq_df['cum_pos'] = freq_df['Win_Mid'] + freq_df['Chromosome'].map(offsets)

    fig, axes = plt.subplots(2, 1, figsize=(16, 8), sharex=True)
    for grp_name, ax in [('case', axes[0]), ('ctrl', axes[1])]:
        for sv_type, color in [('DUP', '#c0392b'), ('DEL', '#2c3e7b')]:
            sub = freq_df[freq_df.TypeNorm == sv_type]
            y = sub[f'Freq_{grp_name}'] * (1 if sv_type == 'DUP' else -1)
            ax.fill_between(sub['cum_pos'], y, 0, color=color, step='mid', alpha=0.8,
                             label=sv_type if grp_name == 'case' else None)
        ax.set_ylabel(f"{grp_name}\nfrequency")
        ax.axhline(0, color='black', lw=0.5)
    axes[0].legend(loc='upper right')
    axes[1].set_xticks([x for _, x in xticks])
    axes[1].set_xticklabels([c for c, _ in xticks], rotation=90, fontsize=8)
    plt.xlabel('Chromosome')
    plt.suptitle('CNV Frequency Profile (Gain up / Loss down): Case vs Control', fontweight='bold')
    plt.tight_layout()
    plt.savefig(OUT_DIR + f"cnv_frequency_profile_{SUFFIX}.png", dpi=300)
    plt.close()

    tick("5~7단계 플롯 및 저장")
    print(f"\n=== 완료 (BURDEN_MODE={BURDEN_MODE}) ===")


if __name__ == "__main__":
    main()
