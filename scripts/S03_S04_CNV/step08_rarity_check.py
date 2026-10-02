import pandas as pd
import glob
import os
import re

# ══════════════════════════════════════════════════════
# 0. 공통 설정
# ══════════════════════════════════════════════════════
ANNOT_BASE_DIR = "/input/data/CNV_all"
SURVIVOR_DIR   = "/input/data/CNV_all/260605_survivor_check"
WINDOW         = 1_000_000
HALF           = WINDOW // 2
TOP_N          = 10
RARE_AF_THRESHOLD = 0.01   # gnomAD-SV 등에서 이 값 미만이면 'RARE'로 판정 (필요시 조정)

BASE_PATTERNS = [
    "UC527KDC_260718CD_check_recipropal{}cnv_window_logit_PCadj_loo.UC.csv",
    "UC526KDC_260716.intcheck{}cnv_window_logit_PCadj_loo.UC.csv",
]
WIN_SIZES = [30, 40]

# ── 문헌 검색으로 확인된 known CNV hotspot / benign 영역 ──────────
# 좌표는 hg38 기준 근사치입니다. 정확한 판정이 필요하면 UCSC/ClinGen에서
# 재확인하세요 (특히 1q21.1은 assembly에 따라 흔들릴 수 있음).
KNOWN_REGIONS = [
    {
        "chrom": "6", "start": 100_000, "end": 900_000,
        "note": "6p25.3 DUSP22/IRF4/EXOC2/HUS1B - ISCA에서 benign CNV로 분류된 영역",
    },
    {
        "chrom": "1", "start": 146_000_000, "end": 153_000_000,
        "note": "1q21.1 recurrent CNV hotspot (segmental duplication/NAHR 매개, 좌표 근사치 - 검증 필요)",
    },
]

def normalize_chrom_value(val) -> str:
    s = str(val).strip()
    s = re.sub(r"(?i)^chr\s*", "", s)
    s = s.strip()
    if s.isdigit():
        s = str(int(s))
    return s

def find_chr_files(chrom: str) -> list[str]:
    patterns = [os.path.join(ANNOT_BASE_DIR, "split*", f"*chr{chrom}.merged.annotated.tsv")]
    if chrom.isdigit() and len(chrom) == 1:
        patterns.append(os.path.join(ANNOT_BASE_DIR, "split*", f"*chr0{chrom}.merged.annotated.tsv"))
    files = []
    for p in patterns:
        files.extend(glob.glob(p))
    return sorted(set(files))

def extract_max_af(value):
    """세미콜론으로 나열된 AFmax 문자열에서 최댓값을 뽑아냄. 값이 없으면 None.

    주의: float("nan"), float("inf") 는 파이썬에서 예외 없이 파싱되기 때문에,
    소스 파일에 결측치가 문자열 'nan'/'NA'/'.' 등으로 저장돼 있으면 그대로
    NaN이 섞여 들어가 max() 및 비교 연산을 오염시킨다. 반드시 명시적으로 걸러낸다.
    """
    if pd.isna(value) or value == "":
        return None
    vals = []
    for p in str(value).split(";"):
        p = p.strip()
        if not p or p.lower() in ("nan", "na", "n/a", ".", "-"):
            continue
        try:
            f = float(p)
        except ValueError:
            continue
        if f != f or f in (float("inf"), float("-inf")):  # f != f  <=>  isnan(f)
            continue
        vals.append(f)
    return max(vals) if vals else None

def check_known_region(chrom: str, win_start: int, win_end: int) -> str:
    hits = []
    for region in KNOWN_REGIONS:
        if region["chrom"] == chrom and win_start <= region["end"] and win_end >= region["start"]:
            hits.append(region["note"])
    return " | ".join(hits) if hits else ""

# ══════════════════════════════════════════════════════
# 1. 파일별 P_GC top10 -> candidate 목록
# ══════════════════════════════════════════════════════
csv_files = [
    os.path.join(SURVIVOR_DIR, pat.format(w))
    for pat in BASE_PATTERNS
    for w in WIN_SIZES
]

top_frames = []
for f in csv_files:
    if not os.path.exists(f):
        print(f"[SKIP] 파일 없음: {f}")
        continue
    df = pd.read_csv(f)
    df["Chromosome"] = df["Chromosome"].apply(normalize_chrom_value)
    df_valid = df.dropna(subset=["P_GC"]).sort_values("P_GC", ascending=True)
    top = df_valid.head(TOP_N).copy()
    top["source_file"] = os.path.basename(f)
    top_frames.append(top)

if not top_frames:
    raise SystemExit("survivor_check CSV를 하나도 읽지 못했습니다. 경로를 확인하세요.")

top_all = pd.concat(top_frames, ignore_index=True)
candidates = list(zip(top_all["Chromosome"].astype(str), top_all["Win_Mid"].astype(int)))

# ══════════════════════════════════════════════════════
# 2. AnnotSV 매칭 + AFmax 기반 rarity 판정 + known region 체크
# ══════════════════════════════════════════════════════
USECOLS_WANTED = [
    "AnnotSV_ID", "SV_chrom", "SV_start", "SV_end", "SV_length", "SV_type",
    "Samples_ID", "Annotation_mode", "Gene_name",
    "B_gain_AFmax", "B_loss_AFmax", "ACMG_class",
]

summary_rows = []

for idx, (chrom, pos) in enumerate(candidates):
    row = top_all.iloc[idx]
    base_info = dict(
        query_chrom=chrom, query_pos=pos,
        query_P_GC=row["P_GC"], query_OR=row["OR"],
        query_source_file=row["source_file"],
    )

    win_start, win_end = pos - HALF, pos + HALF
    known_note = check_known_region(chrom, win_start, win_end)

    files = find_chr_files(chrom)
    if not files:
        summary_rows.append({
            **base_info, "status": "annot_파일_없음", "n_matches": 0,
            "matched_genes": "", "max_AF_in_window": None, "rarity": "N/A",
            "known_region_note": known_note,
        })
        continue

    candidate_hits = []
    for f in files:
        header_cols = pd.read_csv(f, sep="\t", nrows=0).columns
        usecols = [c for c in USECOLS_WANTED if c in header_cols]
        adf = pd.read_csv(f, sep="\t", usecols=usecols, low_memory=False)
        adf["SV_chrom_norm"] = adf["SV_chrom"].apply(normalize_chrom_value)

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
            "matched_genes": "", "max_AF_in_window": None, "rarity": "N/A",
            "known_region_note": known_note,
        })
        continue

    hit_all = pd.concat(candidate_hits, ignore_index=True)

    # gnomAD-SV 등 benign SV DB의 AFmax 중 이 윈도우 내 최댓값 추출
    # 주의: 한 컬럼 안에 None과 float가 섞이면 pandas apply() 결과 dtype이
    # float64로 강제 변환되면서 None이 NaN으로 바뀐다. 그래서 `is not None`이
    # 아니라 pd.notna()로 걸러야 한다 (nan is not None 은 True라서 안 걸러짐).
    af_values = []
    for col in ["B_gain_AFmax", "B_loss_AFmax"]:
        if col in hit_all.columns:
            af_values += [v for v in hit_all[col].apply(extract_max_af).tolist() if pd.notna(v)]
    max_af = max(af_values) if af_values else None

    if max_af is None:
        rarity = "AF_정보없음"
    elif max_af < RARE_AF_THRESHOLD:
        rarity = f"RARE (max_AF={max_af:.4f})"
    else:
        rarity = f"COMMON (max_AF={max_af:.4f}) - benign 가능성"

    genes = set()
    for g in hit_all["Gene_name"].dropna():
        genes.update(str(g).split(";"))
    genes.discard("")

    summary_rows.append({
        **base_info,
        "status": "matched",
        "n_matches": len(hit_all),
        "matched_genes": ";".join(sorted(genes)),
        "max_AF_in_window": max_af,
        "rarity": rarity,
        "known_region_note": known_note,
    })

# ══════════════════════════════════════════════════════
# 3. 저장 + 우선순위 정렬
# ══════════════════════════════════════════════════════
summary_df = pd.DataFrame(summary_rows)

# 우선순위: matched & RARE & known_region_note 없음 & P_GC 낮은 순
def priority_key(r):
    is_rare = str(r["rarity"]).startswith("RARE")
    has_flag = bool(r["known_region_note"])
    return (0 if (r["status"] == "matched" and is_rare and not has_flag) else 1,
            0 if r["status"] == "matched" else 1,
            r["query_P_GC"])

summary_df["_sort_key"] = summary_df.apply(priority_key, axis=1)
summary_df = summary_df.sort_values("_sort_key").drop(columns="_sort_key").reset_index(drop=True)

out_path = "candidates_rarity_summary.tsv"
summary_df.to_csv(out_path, sep="\t", index=False)

n_rare = summary_df["rarity"].astype(str).str.startswith("RARE").sum()
n_common = summary_df["rarity"].astype(str).str.startswith("COMMON").sum()
n_flagged = (summary_df["known_region_note"] != "").sum()

print(f"[요약] 총 {len(summary_df)}개 candidate 중 RARE={n_rare}, COMMON(benign 가능성)={n_common}, "
      f"known hotspot/benign 겹침={n_flagged}")
print(f"[저장] {out_path}\n")

print(summary_df[
    ["query_chrom", "query_pos", "query_P_GC", "status", "rarity", "known_region_note", "matched_genes"]
].to_string(index=False))
