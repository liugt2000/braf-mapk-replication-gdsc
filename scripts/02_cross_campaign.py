"""
02_cross_campaign.py
跨筛选活动（GDSC1 vs GDSC2）的标签可靠性
-------------------------------------------------
GDSC1 与 GDSC2 是两个独立进行的筛选活动（不同化合物批次、不同浓度梯度设计、
不同时间）。同一细胞系 × 同一药物在两个活动中各测一次，构成一个天然的
"重测信度"实验，且这两个活动都不是为本次分析做的。

输出：results/02_cross_campaign.json, results/table3_cross_campaign.csv
"""
import json
import numpy as np
import pandas as pd

g1 = pd.read_csv("data/raw/GDSC1_fitted_dose_response_24Jul22.csv")
g2 = pd.read_csv("data/raw/GDSC2_fitted_dose_response_24Jul22.csv")

c1 = g1[g1.TCGA_DESC == "COREAD"]
c2 = g2[g2.TCGA_DESC == "COREAD"]
print(f"GDSC1 COREAD: {len(c1)} curves, {c1.CELL_LINE_NAME.nunique()} lines, {c1.DRUG_NAME.nunique()} drugs")
print(f"GDSC2 COREAD: {len(c2)} curves, {c2.CELL_LINE_NAME.nunique()} lines, {c2.DRUG_NAME.nunique()} drugs")

# 活动内先按 (系, 药) 取均值，避免一对多
a = c1.groupby(["CELL_LINE_NAME", "DRUG_NAME"], as_index=False)["LN_IC50"].mean()
a.columns = ["line", "drug", "ln1"]
b = c2.groupby(["CELL_LINE_NAME", "DRUG_NAME"], as_index=False)["LN_IC50"].mean()
b.columns = ["line", "drug", "ln2"]

m = a.merge(b, on=["line", "drug"], how="inner")
print(f"\n两活动共同覆盖的 (细胞系, 药物) 对: {len(m)}  "
      f"({m.line.nunique()} 株系 x {m.drug.nunique()} 个药)")


def stats(x, y):
    x = np.asarray(x, float); y = np.asarray(y, float)
    n = len(x)
    r = np.corrcoef(x, y)[0, 1]
    d = y - x
    grand = np.mean(np.concatenate([x, y]))
    ms_b = 2 * (((x + y) / 2 - grand) ** 2).sum() / (n - 1)
    ms_w = ((x - y) ** 2).sum() / n
    icc = (ms_b - ms_w) / (ms_b + ms_w)
    return {
        "n_pairs": int(n),
        "pearson_r": round(float(r), 3),
        "r_squared": round(float(r ** 2), 3),
        "spearman_rho": round(float(pd.Series(x).corr(pd.Series(y), method="spearman")), 3),
        "icc_1_1": round(float(icc), 3),
        "median_abs_diff_ln": round(float(np.median(np.abs(d))), 3),
        "median_fold_diff": round(float(np.exp(np.median(np.abs(d)))), 2),
        "pct_pairs_fold_gt_2": round(100 * float((np.abs(d) > np.log(2)).mean()), 1),
        "pct_pairs_fold_gt_4": round(100 * float((np.abs(d) > np.log(4)).mean()), 1),
    }


rows = [{"subset": "ALL COREAD drug pairs", "n_drugs": m.drug.nunique(),
         **stats(m.ln1, m.ln2)}]

CHEMO = ["5-Fluorouracil", "Oxaliplatin", "Irinotecan", "SN-38"]
mc = m[m.drug.isin(CHEMO)]
if len(mc) > 3:
    rows.append({"subset": "CRC chemo backbone", "n_drugs": mc.drug.nunique(),
                 **stats(mc.ln1, mc.ln2)})

for drug in CHEMO + ["Trametinib", "Dabrafenib", "Selumetinib", "Fulvestrant"]:
    d = m[m.drug == drug]
    if len(d) > 5:
        rows.append({"subset": f"drug: {drug}", "n_drugs": 1, **stats(d.ln1, d.ln2)})

tab = pd.DataFrame(rows)
tab.to_csv("results/table3_cross_campaign.csv", index=False)
print("\n=== TABLE 3: GDSC1 vs GDSC2 cross-campaign reliability ===")
print(tab.to_string(index=False))

out = {
    "coread_lines_gdsc1": int(c1.CELL_LINE_NAME.nunique()),
    "coread_lines_gdsc2": int(c2.CELL_LINE_NAME.nunique()),
    "overlap_pairs": int(len(m)),
    "overlap_lines": int(m.line.nunique()),
    "overlap_drugs": int(m.drug.nunique()),
    "overall": stats(m.ln1, m.ln2),
}
if len(mc) > 3:
    out["chemo_backbone"] = stats(mc.ln1, mc.ln2)
with open("results/02_cross_campaign.json", "w", encoding="utf-8") as f:
    json.dump(out, f, indent=2, ensure_ascii=False)

# 供后续使用的宽表
m.to_csv("results/coread_ic50_two_campaigns.csv", index=False)
