import pandas as pd
import glob
import os
import re

# ══════════════════════════════════════════════════════
# 0. 공통 설정
# ══════════════════════════════════════════════════════
ANNOT_BASE_DIR = "/input/data/CNV_all"                       # split1, split2 ... 상위 폴더
SURVIVOR_DIR   = "/input/data/CNV_all/260605_survivor_check" # P_GC csv들이 있는 폴더
WINDOW         = 1_000_000   # AnnotSV 매칭 시 후보 좌표 기준 총 윈도우 크기
HALF           = WINDOW // 2
TOP_N          = 10          # 파일별로 뽑을 top N (P_GC 오름차순 = 가장 유의한 순)

def normalize_chrom_value(val) -> str:
    """
    'Chr 20' -> '20', 'chr05' -> '5', ' 09 ' -> '9', 'ChrX' -> 'X'
    CSV의 'Chr 20'(공백+0패딩)과 AnnotSV 쪽 'chr20' 표기를 동일한 값으로 맞추기 위함
    """
    s = str(val).strip()
    s = re.sub(r"(?i)^chr\s*", "", s)   # 'chr'/'Chr ' 접두사 제거
    s = s.strip()
    if s.isdigit():
        s = str(int(s))                 # 앞자리 0 제거: '05' -> '5'
    return s

# split* 폴더/파일 명명 규칙을 실행 초반에 한 번 확인 (진단용)
_split_dirs = sorted(glob.glob(os.path.join(ANNOT_BASE_DIR, "split*")))
if not _split_dirs:
    print(f"[경고] {ANNOT_BASE_DIR} 아래에 split* 폴더를 찾지 못했습니다. ANNOT_BASE_DIR 경로를 확인하세요.")
else:
    _sample = glob.glob(os.path.join(_split_dirs[0], "*.tsv"))[:3]
    print(f"[진단] split* 폴더 {len(_split_dirs)}개 발견. 예시 파일명: {_sample}")

# UC527KDC_..._recipropal{30,40}..., UC526KDC_..._intcheck{30,40}... 총 4개 파일
BASE_PATTERNS = [
    "UC527KDC_260718CD_check_recipropal{}cnv_window_logit_PCadj_loo.UC.csv",
    "UC526KDC_260716.intcheck{}cnv_window_logit_PCadj_loo.UC.csv",
]
WIN_SIZES = [30, 40]

csv_files = [
    os.path.join(SURVIVOR_DIR, pat.format(w))
    for pat in BASE_PATTERNS
    for w in WIN_SIZES
]

# ══════════════════════════════════════════════════════
# 1. 파일별 P_GC top10 -> candidate 목록 생성
# ══════════════════════════════════════════════════════
top_frames = []

for f in csv_files:
    if not os.path.exists(f):
        print(f"[SKIP] 파일 없음: {f}")
        continue

    df = pd.read_csv(f)
    df["Chromosome"] = df["Chromosome"].apply(normalize_chrom_value)  # 'Chr 20' -> '20'
    # P_GC가 NaN인 행 제외, 오름차순 정렬(가장 유의한=작은 값이 위)
    df_valid = df.dropna(subset=["P_GC"]).sort_values("P_GC", ascending=True)
    top = df_valid.head(TOP_N).copy()
    top["source_file"] = os.path.basename(f)
    top_frames.append(top)

if not top_frames:
    raise SystemExit("survivor_check CSV를 하나도 읽지 못했습니다. 경로/파일명을 확인하세요.")

top_all = pd.concat(top_frames, ignore_index=True)

top_out_path = "top_PGC_candidates.tsv"
top_all.to_csv(top_out_path, sep="\t", index=False)
print(f"[1단계] 파일 {len(top_frames)}개에서 top{TOP_N}씩 추출 -> 총 {len(top_all)}개 candidate ({top_out_path})")
print(f"[진단] 정규화된 Chromosome 값들: {sorted(top_all['Chromosome'].unique(), key=lambda x: (len(x), x))}")

# candidates: (chrom, pos) 리스트로 변환 — Chromosome/Win_Mid를 그대로 사용
# (같은 chrom이 여러 번 나와도 상관없도록 list-of-tuple 유지)
candidates = list(
    zip(top_all["Chromosome"].astype(str), top_all["Win_Mid"].astype(int))
)

# ══════════════════════════════════════════════════════
# 2. AnnotSV 매칭 (이전 스크립트와 동일 로직)
# ══════════════════════════════════════════════════════
USECOLS_WANTED = [
    "AnnotSV_ID", "SV_chrom", "SV_start", "SV_end", "SV_length", "SV_type",
    "Samples_ID", "Annotation_mode", "Gene_name",
    "B_gain_AFmax", "B_loss_AFmax", "ACMG_class",
]

def normalize_chrom(series: pd.Series) -> pd.Series:
    """SV_chrom 컬럼도 동일 규칙으로 정규화 (Chromosome 컬럼과 동일 함수 재사용)"""
    return series.apply(normalize_chrom_value)

def find_chr_files(chrom: str) -> list[str]:
    """chrom='5' 든 실제 파일이 chr05로 되어있든 둘 다 잡히도록 0패딩 버전도 같이 탐색"""
    patterns = [os.path.join(ANNOT_BASE_DIR, "split*", f"*chr{chrom}.merged.annotated.tsv")]
    if chrom.isdigit() and len(chrom) == 1:
        patterns.append(os.path.join(ANNOT_BASE_DIR, "split*", f"*chr0{chrom}.merged.annotated.tsv"))
    files = []
    for p in patterns:
        files.extend(glob.glob(p))
    return sorted(set(files))

results = []
summary_rows = []

for idx, (chrom, pos) in enumerate(candidates):
    row = top_all.iloc[idx]
    base_info = dict(
        query_chrom=chrom,
        query_pos=pos,
        query_P_GC=row["P_GC"],
        query_OR=row["OR"],
        query_source_file=row["source_file"],
    )

    files = find_chr_files(chrom)
    if not files:
        summary_rows.append({
            **base_info, "status": "annot_파일_없음", "n_matches": 0,
            "matched_genes": "", "matched_SV_types": "", "matched_ACMG_class": "",
        })
        print(f"[chr{chrom}:{pos:,}] (P_GC={row['P_GC']:.3g}) -> annot 파일 없음")
        continue

    win_start, win_end = pos - HALF, pos + HALF
    candidate_hits = []

    for f in files:
        header_cols = pd.read_csv(f, sep="\t", nrows=0).columns
        usecols = [c for c in USECOLS_WANTED if c in header_cols]

        adf = pd.read_csv(f, sep="\t", usecols=usecols, low_memory=False)
        adf["SV_chrom_norm"] = normalize_chrom(adf["SV_chrom"])

        hit = adf[
            (adf["SV_chrom_norm"] == chrom)
            & (adf["SV_start"] <= win_end)
            & (adf["SV_end"] >= win_start)
        ].copy()

        if not hit.empty:
            hit["annot_source_file"] = os.path.basename(f)
            candidate_hits.append(hit)

    if not candidate_hits:
        summary_rows.append({
            **base_info, "status": "no_match", "n_matches": 0,
            "matched_genes": "", "matched_SV_types": "", "matched_ACMG_class": "",
        })
        print(f"[chr{chrom}:{pos:,}] (P_GC={row['P_GC']:.3g}) -> ±{HALF/1e6:.1f}Mb 내 매칭 SV 없음")
        continue

    hit_all = pd.concat(candidate_hits, ignore_index=True)
    for k, v in base_info.items():
        hit_all[k] = v
    results.append(hit_all)

    genes = set()
    for g in hit_all["Gene_name"].dropna():
        genes.update(str(g).split(";"))
    genes.discard("")
    sv_types = sorted(hit_all["SV_type"].dropna().unique().tolist())
    acmg = (
        sorted(hit_all["ACMG_class"].dropna().astype(str).unique().tolist())
        if "ACMG_class" in hit_all.columns else []
    )

    summary_rows.append({
        **base_info,
        "status": "matched",
        "n_matches": len(hit_all),
        "matched_genes": ";".join(sorted(genes)),
        "matched_SV_types": ";".join(sv_types),
        "matched_ACMG_class": ";".join(acmg),
    })
    gene_preview = ", ".join(sorted(genes)[:5]) + ("..." if len(genes) > 5 else "")
    print(f"[chr{chrom}:{pos:,}] (P_GC={row['P_GC']:.3g}) -> {len(hit_all)}개 SV 매칭 | gene: {gene_preview}")

# ══════════════════════════════════════════════════════
# 3. 저장 — ① 상세 매칭 결과, ② 한눈에 보는 요약 테이블
# ══════════════════════════════════════════════════════
# ① 상세 결과 (matched candidate에 한해 SV 단위로 여러 줄)
if results:
    out = pd.concat(results, ignore_index=True)
    front_cols = ["query_chrom", "query_pos", "query_P_GC", "query_OR", "query_source_file", "annot_source_file"]
    other_cols = [c for c in out.columns if c not in front_cols + ["SV_chrom_norm"]]
    out = out[front_cols + other_cols].sort_values("query_P_GC")
    out.to_csv("matched_candidates_from_PGC.tsv", sep="\t", index=False)

# ② candidate 40개 전체를 한 줄씩 — 매칭 여부 status 포함
summary_df = pd.DataFrame(summary_rows).sort_values("query_P_GC").reset_index(drop=True)
summary_df.to_csv("candidates_summary.tsv", sep="\t", index=False)

n_matched = (summary_df["status"] == "matched").sum()
n_no_match = (summary_df["status"] == "no_match").sum()
n_no_file = (summary_df["status"] == "annot_파일_없음").sum()

print(f"\n[요약] 총 {len(summary_df)}개 candidate 중 matched={n_matched}, no_match={n_no_match}, annot_파일_없음={n_no_file}")
print(f"[요약] 전체 결과 -> candidates_summary.tsv (한눈에 보기용)")
print(f"[상세] SV 단위 상세 결과 -> matched_candidates_from_PGC.tsv (matched만)")

# 콘솔에서 바로 한눈에 볼 수 있도록 핵심 컬럼만 출력
print("\n" + summary_df[
    ["query_chrom", "query_pos", "query_P_GC", "status", "n_matches", "matched_genes"]
].to_string(index=False))
