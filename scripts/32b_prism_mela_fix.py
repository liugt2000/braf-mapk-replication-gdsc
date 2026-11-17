"""32b 补丁: PRISM 黑色素瘤分类修正 (PRISM 组织标签用 skin 而非 melanoma) + PRISM 解离检验
从 results/32_prism_per_line.csv 直接重算, 不重解析 RDS。"""
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(r"D:\AI_work\workbuddy\2026-11-16-16-13-17")
OUT = ROOT / "results"

per = pd.read_csv(OUT / "32_prism_per_line.csv")
per["is_mela"] = per["lab"].str.contains("melanoma|skin", na=False)

pv = per[per["status"] == "V600"].pivot_table(
    index=["cellid", "key", "is_crc", "is_mela"], columns="drugid",
    values="aac_recomputed").reset_index()
pv = pv.dropna(subset=["Dabrafenib", "Trametinib"])
pv["daac"] = pv["Trametinib"] - pv["Dabrafenib"]

n_crc = int(pv["is_crc"].sum()); n_mel = int(pv["is_mela"].sum())
print(f"[PRISM dab+tra 修正后] V600 双药可配 = {len(pv)} (CRC={n_crc}, 黑色素瘤={n_mel})")
print("  CRC 系:", pv.loc[pv.is_crc, "cellid"].tolist())
print("  黑色素瘤系:", pv.loc[pv.is_mela, "cellid"].tolist())
print("  组织标签核对:", per[per.cellid.isin(pv.loc[pv.is_mela, "cellid"])]
      [["cellid", "lab"]].drop_duplicates().to_string(index=False))


def mwu(a, b):
    a = np.asarray(a, float); a = a[~np.isnan(a)]
    b = np.asarray(b, float); b = b[~np.isnan(b)]
    if len(a) < 2 or len(b) < 2:
        return dict(n1=len(a), n2=len(b), med1=np.nan, med2=np.nan, U=np.nan, p=np.nan)
    U, p = stats.mannwhitneyu(a, b, alternative="two-sided")
    return dict(n1=len(a), n2=len(b), med1=float(np.median(a)), med2=float(np.median(b)),
                U=float(U), p=float(p))


new_rows = []
if n_mel >= 2:
    for label, sel in [("AAC_Dabrafenib", "Dabrafenib"),
                       ("AAC_Trametinib", "Trametinib"),
                       ("dAAC_Trametinib_minus_Dabrafenib", "daac")]:
        r = mwu(pv.loc[pv.is_crc, sel], pv.loc[pv.is_mela, sel])
        new_rows.append(dict(dataset="PRISM", test=f"V600_CRC_vs_melanoma_{label}",
                             drug="Dabrafenib+Trametinib", metric="aac_recomputed", **r))
    r = mwu(pv.loc[pv.is_crc, "daac"], pv.loc[~pv.is_crc, "daac"])
    new_rows.append(dict(dataset="PRISM", test="V600_CRC_vs_allothers_dAAC",
                         drug="Dabrafenib+Trametinib", metric="aac_recomputed", **r))

# 追加到结果 CSV (替换旧 PRISM 解离行, 防重)
res = pd.read_csv(OUT / "32_third_party_replication.csv")
res = res[~((res.dataset == "PRISM") & (res.test.str.startswith("V600_")))]
res = pd.concat([res, pd.DataFrame(new_rows)], ignore_index=True)
res.to_csv(OUT / "32_third_party_replication.csv", index=False)
pv.to_csv(OUT / "32_prism_pairs_dab_tra.csv", index=False)

# 更新 summary JSON
with open(OUT / "32_summary.json", encoding="utf-8") as f:
    summary = json.load(f)
summary["prism"]["dissociation_note"] = (
    "黑色素瘤分类修正为 lab 含 melanoma|skin (PRISM 标签用 skin); "
    "CRC n=%d, melanoma n=%d" % (n_crc, n_mel))
summary["prism"]["tests"] = [r for r in res.to_dict("records")
                             if r["dataset"] == "PRISM"]
with open(OUT / "32_summary.json", "w", encoding="utf-8") as f:
    json.dump(summary, f, indent=2, ensure_ascii=False, default=float)

print("\n新增 PRISM 解离行:")
print(pd.DataFrame(new_rows).to_string(index=False))
