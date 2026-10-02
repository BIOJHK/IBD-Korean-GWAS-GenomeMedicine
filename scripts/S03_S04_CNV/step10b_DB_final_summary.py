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
RARE_AF_THRESHOLD = 0.01

BASE_PATTERNS = [
    "UC527KDC_260718CD_check_recipropal{}cnv_window_logit_PCadj_loo.UC.csv",
    "UC526KDC_260716.intcheck{}cnv_window_logit_PCadj_loo.UC.csv",
]
WIN_SIZES = [30, 40]

KNOWN_REGIONS = [
    {"chrom": "6", "start": 100_000, "end": 900_000,
     "note": "6p25.3 DUSP22/IRF4/EXOC2/HUS1B - ISCA benign CNV"},
    {"chrom": "1", "start": 146_000_000, "end": 153_000_000,
     "note": "1q21.1 recurrent CNV hotspot (좌표 근사치)"},
]

# evidence tier 우선순위: 낮을수록 "아직 살아있는 candidate" (follow-up 우선)
TIER_PRIORITY = {
    "RARE": 0,
    "NO_DATA": 1,
    "UNRELIABLE_dbVar_placeholder": 1,
    "CAUTION_DDD": 2,
    "RELIABLE_COMMON": 3,
}

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
    if pd.isna(value) or value == "":
        return []
    out = []
    for p in str(value).split(";"):
        p = p.strip()
        if not p or p.lower() in ("nan", "na", "n/a", ".", "-"):
            out.append(None); continue
        try:
            f = float(p)
        except ValueError:
            out.append(None); continue
        if f != f or f in (float("inf"), float("-inf")):
            out.append(None); continue
        out.append(f)
    return out

def parse_source_list(value):
    if pd.isna(value) or value == "":
        return []
    return [p.strip() for p in str(value).split(";")]

def best_source_for_max_af(af_str, source_str):
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

def check_known_region(chrom, win_start, win_end):
    hits = [r["note"] for r in KNOWN_REGIONS
            if r["chrom"] == chrom and win_start <= r["end"] and win_end >= r["start"]]
    return " | ".join(hits) if hits else ""

def classify_tier(max_af, source):
    if max_af is None or source is None:
        return "NO_DATA"
    src_lower = str(source).lower()
    if src_lower.startswith("ddd"):
        return "CAUTION_DDD"
    if src_lower.startswith("dbvar") and abs(max_af - RARE_AF_THRESHOLD) < 1e-6:
        return "UNRELIABLE_dbVar_placeholder"
    return "RARE" if max_af < RARE_AF_THRESHOLD else "RELIABLE_COMMON"

# ══════════════════════════════════════════════════════
# 1. 파일별 P_GC top10 -> candidate 목록
# ══════════════════════════════════════════════════════
csv_files = [os.path.join(SURVIVOR_DIR, pat.format(w)) for pat in BASE_PATTERNS for w in WIN_SIZES]

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
# 2. AnnotSV 매칭 + tier 분류
# ══════════════════════════════════════════════════════
USECOLS_WANTED = [
    "SV_chrom", "SV_start", "SV_end", "SV_type", "Gene_name",
    "B_gain_AFmax", "B_gain_source", "B_loss_AFmax", "B_loss_source",
]

rows = []

for idx, (chrom, pos) in enumerate(candidates):
    row = top_all.iloc[idx]
    win_start, win_end = pos - HALF, pos + HALF
    known_note = check_known_region(chrom, win_start, win_end)

    base = dict(
        query_chrom=chrom, query_pos=pos,
        query_P_GC=row["P_GC"], query_OR=row["OR"],
        query_source_file=row["source_file"],
    )

    files = find_chr_files(chrom)
    if not files:
        rows.append({**base, "status": "annot_파일_없음", "matched_genes": "",
                     "max_AF": None, "max_AF_source": None,
                     "evidence_tier": "NO_DATA", "known_region_note": known_note})
        continue

    candidate_hits = []
    for f in files:
        header_cols = pd.read_csv(f, sep="\t", nrows=0).columns
        usecols = [c for c in USECOLS_WANTED if c in header_cols]
        adf = pd.read_csv(f, sep="\t", usecols=usecols, low_memory=False)
        adf["SV_chrom_norm"] = adf["SV_chrom"].apply(normalize_chrom_value)
        hit = adf[(adf["SV_chrom_norm"] == chrom) & (adf["SV_start"] <= win_end) & (adf["SV_end"] >= win_start)].copy()
        if not hit.empty:
            candidate_hits.append(hit)

    if not candidate_hits:
        rows.append({**base, "status": "no_match", "matched_genes": "",
                     "max_AF": None, "max_AF_source": None,
                     "evidence_tier": "NO_DATA", "known_region_note": known_note})
        continue

    hit_all = pd.concat(candidate_hits, ignore_index=True)

    overall_best_af, overall_best_src = None, None
    for col_af, col_src in [("B_gain_AFmax", "B_gain_source"), ("B_loss_AFmax", "B_loss_source")]:
        if col_af not in hit_all.columns or col_src not in hit_all.columns:
            continue
        for af_str, src_str in zip(hit_all[col_af], hit_all[col_src]):
            af, src = best_source_for_max_af(af_str, src_str)
            if af is not None and (overall_best_af is None or af > overall_best_af):
                overall_best_af, overall_best_src = af, src

    genes = set()
    for g in hit_all["Gene_name"].dropna():
        genes.update(str(g).split(";"))
    genes.discard("")

    tier = classify_tier(overall_best_af, overall_best_src)

    rows.append({
        **base, "status": "matched",
        "matched_genes": ";".join(sorted(genes)),
        "max_AF": overall_best_af, "max_AF_source": overall_best_src,
        "evidence_tier": tier, "known_region_note": known_note,
    })

# ══════════════════════════════════════════════════════
# 3. 최종 정렬 + 저장
# ══════════════════════════════════════════════════════
full_df = pd.DataFrame(rows)
full_df["_tier_rank"] = full_df["evidence_tier"].map(TIER_PRIORITY)
full_df = full_df.sort_values(["_tier_rank", "query_P_GC"]).drop(columns="_tier_rank").reset_index(drop=True)
full_df.to_csv("final_candidates_full.tsv", sep="\t", index=False)

# locus(chrom,pos) 단위로 중복 제거 — 여러 그룹(CD30/CD40/UC30/UC40)에 걸쳐 나온 candidate는
# 가장 낮은 P_GC(가장 유의한) 값만 대표로 남긴다.
unique_df = (
    full_df.sort_values("query_P_GC")
    .drop_duplicates(subset=["query_chrom", "query_pos"], keep="first")
)
unique_df["_tier_rank"] = unique_df["evidence_tier"].map(TIER_PRIORITY)
unique_df = unique_df.sort_values(["_tier_rank", "query_P_GC"]).drop(columns="_tier_rank").reset_index(drop=True)
unique_df.to_csv("final_candidates_unique_loci.tsv", sep="\t", index=False)

print(f"[저장] final_candidates_full.tsv (40행, 그룹별 중복 포함)")
print(f"[저장] final_candidates_unique_loci.tsv ({len(unique_df)}행, locus 단위 중복 제거)\n")

print("[Tier별 개수 — unique locus 기준]")
print(unique_df["evidence_tier"].value_counts().to_string())

print("\n[최종 우선순위 테이블 — unique locus, RARE/NO_DATA/DDD 순]")
print(unique_df[
    ["query_chrom", "query_pos", "query_P_GC", "evidence_tier",
     "max_AF", "max_AF_source", "known_region_note", "matched_genes"]
].to_string(index=False))
