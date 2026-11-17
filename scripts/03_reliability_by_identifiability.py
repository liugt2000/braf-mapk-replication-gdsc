"""
03_reliability_by_identifiability.py
把"外推"与"真噪声"分开
-------------------------------------------------
问题：GDSC1/GDSC2 的 ICC 低，到底是因为 IC50 被外推到测试浓度之外（设计问题，
可以用筛选规则修好），还是因为测量本身不可重复（本质问题，修不好）？

做法：按"两个活动中 IC50 是否都落在各自实测浓度区间内"分层，重算可靠性。
      若限定可识别后 ICC 大幅回升 -> 问题在设计，有救。
      若仍然很低         -> 问题在测量本质，模型上限被钉死。

输出：results/03_reliability_stratified.csv, results/03_reliability_stratified.json
"""
import json
import numpy as np
import pandas as pd


def load(path, camp):
    d = pd.read_csv(path)
    d = d[d.TCGA_DESC == "COREAD"].copy()
    d["campaign"] = camp
    d["ic50"] = np.exp(d.LN_IC50)
    d["identifiable"] = (d.ic50 >= d.MIN_CONC) & (d.ic50 <= d.MAX_CONC)
    # 活动内按 (系,药) 聚合：保留均值与"是否全部可识别"
    agg = d.groupby(["CELL_LINE_NAME", "DRUG_NAME"]).agg(
        ln=("LN_IC50", "mean"),
        ident_all=("identifiable", "all"),
        n=("LN_IC50", "size"),
    ).reset_index()
    return agg.rename(columns={"ln": f"ln_{camp}", "ident_all": f"ident_{camp}", "n": f"n_{camp}"})


a = load("data/raw/GDSC1_fitted_dose_response_24Jul22.csv", "g1")
b = load("data/raw/GDSC2_fitted_dose_response_24Jul22.csv", "g2")
m = a.merge(b, on=["CELL_LINE_NAME", "DRUG_NAME"], how="inner")


def icc(x, y):
    x = np.asarray(x, float); y = np.asarray(y, float); n = len(x)
    if n < 5:
        return np.nan, np.nan
    grand = np.mean(np.concatenate([x, y]))
    ms_b = 2 * (((x + y) / 2 - grand) ** 2).sum() / (n - 1)
    ms_w = ((x - y) ** 2).sum() / n
    r = np.corrcoef(x, y)[0, 1]
    return (ms_b - ms_w) / (ms_b + ms_w), r


CHEMO = ["5-Fluorouracil", "Oxaliplatin", "Irinotecan", "SN-38"]
m["axis"] = np.where(m.DRUG_NAME.isin(CHEMO), "chemo", "other")
m["both_identifiable"] = m.ident_g1 & m.ident_g2

rows = []
for axis, sub in [("ALL", m), ("chemo", m[m.axis == "chemo"])]:
    for lab, s in [("all pairs", sub), ("both identifiable", sub[sub.both_identifiable])]:
        if len(s) < 5:
            continue
        i, r = icc(s.ln_g1, s.ln_g2)
        d = np.abs(s.ln_g1 - s.ln_g2)
        rows.append({
            "subset": axis, "stratum": lab, "n_pairs": len(s),
            "pct_of_subset": round(100 * len(s) / len(sub), 1),
            "icc_1_1": round(float(i), 3),
            "pearson_r": round(float(r), 3),
            "r2_ceiling": round(float(r ** 2), 3),
            "median_fold_diff": round(float(np.exp(np.median(d))), 2),
            "pct_fold_gt_2": round(100 * float((d > np.log(2)).mean()), 1),
        })

# 逐药
for drug in CHEMO + ["Trametinib", "Selumetinib", "Dabrafenib", "Fulvestrant"]:
    s = m[m.DRUG_NAME == drug]
    for lab, ss in [("all pairs", s), ("both identifiable", s[s.both_identifiable])]:
        if len(ss) < 5:
            continue
        i, r = icc(ss.ln_g1, ss.ln_g2)
        d = np.abs(ss.ln_g1 - ss.ln_g2)
        rows.append({
            "subset": f"drug:{drug}", "stratum": lab, "n_pairs": len(ss),
            "pct_of_subset": round(100 * len(ss) / len(s), 1),
            "icc_1_1": round(float(i), 3), "pearson_r": round(float(r), 3),
            "r2_ceiling": round(float(r ** 2), 3),
            "median_fold_diff": round(float(np.exp(np.median(d))), 2),
            "pct_fold_gt_2": round(100 * float((d > np.log(2)).mean()), 1),
        })

tab = pd.DataFrame(rows)
tab.to_csv("results/03_reliability_stratified.csv", index=False)
print("=== 可靠性按可识别性分层 ===")
print(tab.to_string(index=False))

json.dump({
    "n_overlap_pairs": int(len(m)),
    "pct_pairs_both_identifiable": round(100 * float(m.both_identifiable.mean()), 1),
    "table": tab.to_dict("records"),
}, open("results/03_reliability_stratified.json", "w", encoding="utf-8"),
    indent=2, ensure_ascii=False)
m.to_csv("results/coread_pairs_with_identifiability.csv", index=False)
