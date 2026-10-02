import pandas as pd
import glob
import os
import re

# ══════════════════════════════════════════════════════
# 0. 공통 설정 (기존 스크립트와 동일)
# ══════════════════════════════════════════════════════
ANNOT_BASE_DIR = "/input/data/CNV_all"
SURVIVOR_DIR   = "/input/data/CNV_all/260605_survivor_check"
WINDOW         = 1_000_000
HALF           = WINDOW // 2
TOP_N          = 10

BASE_PATTERNS = [
    "UC527KDC_260718CD_check_recipropal{}cnv_window_logit_PCadj_loo.UC.csv",
    "UC526KDC_260716.intcheck{}cnv_window_logit_PCadj_loo.UC.csv",
]
WIN_SIZES = [30, 40]

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

def parse_af_list(value):
    """세미콜론 리스트를 float로 파싱. 실패/nan/inf 항목은 None으로 남겨서
    source 리스트와 인덱스를 그대로 맞출 수 있게 한다."""
    if pd.isna(value) or value == "":
        return []
    out = []
    for p in str(value).split(";"):
        p = p.strip()
        if not p or p.lower() in ("nan", "na", "n/a", ".", "-"):
            out.append(None)
            continue
        try:
            f = float(p)
        except ValueError:
            out.append(None)
            continue
        if f != f or f in (float("inf"), float("-inf")):
            out.append(None)
            continue
        out.append(f)
    return out

def parse_source_list(value):
    if pd.isna(value) or value == "":
        return []
    return [p.strip() for p in str(value).split(";")]

def best_source_for_max_af(af_str, source_str):
    """AFmax 리스트와 source 리스트를 같은 인덱스로 짝지어서,
    최댓값과 그 값을 낸 source를 함께 반환. 값이 없으면 (None, None)."""
    afs = parse_af_list(af_str)
    srcs = parse_source_list(source_str)
    best_af, best_src = None, None
    for i, af in enumerate(afs):
        if af is None:
            continue
        src = srcs[i] if i < len(srcs) else "(source 불명)"
        if best_af is None or af > best_af:
            best_af, best_src = af, src
    return best_af, best_src

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
    raise SystemExit("survivor_check CSV를 하나도 읽지 못했습니다.")

top_all = pd.concat(top_frames, ignore_index=True)
candidates = list(zip(top_all["Chromosome"].astype(str), top_all["Win_Mid"].astype(int)))

# ══════════════════════════════════════════════════════
# 2. AnnotSV 매칭 + 최댓값의 source 추적
# ══════════════════════════════════════════════════════
USECOLS_WANTED = [
    "SV_chrom", "SV_start", "SV_end", "SV_type", "Gene_name",
    "B_gain_AFmax", "B_gain_source", "B_loss_AFmax", "B_loss_source",
]

rows_out = []

for idx, (chrom, pos) in enumerate(candidates):
    row = top_all.iloc[idx]
    win_start, win_end = pos - HALF, pos + HALF

    files = find_chr_files(chrom)
    if not files:
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
            candidate_hits.append(hit)

    if not candidate_hits:
        continue

    hit_all = pd.concat(candidate_hits, ignore_index=True)

    # 이 candidate 윈도우 안의 모든 SV 행에서 (max_af, source) 후보를 다 모은다
    overall_best_af, overall_best_src, overall_best_type = None, None, None
    for col_af, col_src, tag in [
        ("B_gain_AFmax", "B_gain_source", "gain"),
        ("B_loss_AFmax", "B_loss_source", "loss"),
    ]:
        if col_af not in hit_all.columns or col_src not in hit_all.columns:
            continue
        for af_str, src_str in zip(hit_all[col_af], hit_all[col_src]):
            af, src = best_source_for_max_af(af_str, src_str)
            if af is None:
                continue
            if overall_best_af is None or af > overall_best_af:
                overall_best_af, overall_best_src, overall_best_type = af, src, tag

    rows_out.append({
        "query_chrom": chrom,
        "query_pos": pos,
        "query_P_GC": row["P_GC"],
        "max_AF": overall_best_af,
        "max_AF_type": overall_best_type,   # gain(중복) 또는 loss(결실) 유래
        "max_AF_source": overall_best_src,  # 어느 DB/ID에서 온 값인지
    })

# ══════════════════════════════════════════════════════
# 3. 저장 + 출력
# ══════════════════════════════════════════════════════
out_df = pd.DataFrame(rows_out).sort_values("query_P_GC").reset_index(drop=True)
out_df.to_csv("candidates_af_source.tsv", sep="\t", index=False)

print(f"[저장] candidates_af_source.tsv\n")
print(out_df.to_string(index=False))

print("\n[요약] max_AF_source 값별 등장 빈도:")
print(out_df["max_AF_source"].value_counts().to_string())
