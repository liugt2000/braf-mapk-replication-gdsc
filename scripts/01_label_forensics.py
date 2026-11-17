"""
01_label_forensics.py
GDSC2 / COREAD 药物反应标签的取证审计
-------------------------------------------------
目的：在建模之前先回答一个问题——"要预测的这个标签，本身有多可靠？"
输出：results/01_label_forensics.json  +  results/table1_drug_panel.csv
      results/table2_replicate_reliability.csv
"""
import json, os
import numpy as np
import pandas as pd

RAW = "data/raw/GDSC2_fitted_dose_response_24Jul22.csv"
OUT = "results"
os.makedirs(OUT, exist_ok=True)

df = pd.read_csv(RAW)
df["ln_ic50"] = df["LN_IC50"]
df["ic50"] = np.exp(df["ln_ic50"])

# ---------- COREAD 子集 ----------
crc = df[df.TCGA_DESC == "COREAD"].copy()
crc["below_range"] = crc.ic50 < crc.MIN_CONC
crc["above_range"] = crc.ic50 > crc.MAX_CONC
crc["identifiable"] = ~(crc.below_range | crc.above_range)

# ---------- 药物面板：临床用药 + 通路匹配 + 阴性对照 ----------
PANEL = {
    "5-Fluorouracil": "A_chemo",
    "Oxaliplatin":    "A_chemo",
    "Irinotecan":     "A_chemo",
    "SN-38":          "A_chemo",
    "Dabrafenib":     "B_pathway",
    "Trametinib":     "B_pathway",
    "Selumetinib":    "B_pathway",
    "Fulvestrant":    "C_negative",   # ER 拮抗剂，与 CRC 无关
}

rows = []
for drug, axis in PANEL.items():
    d = crc[crc.DRUG_NAME == drug]
    if d.empty:
        continue
    rows.append({
        "drug": drug,
        "axis": axis,
        "n_curves": len(d),
        "n_lines": d.CELL_LINE_NAME.nunique(),
        "pct_extrapolated_high": round(100 * d.above_range.mean(), 1),
        "pct_identifiable": round(100 * d.identifiable.mean(), 1),
        "median_RMSE": round(d.RMSE.median(), 3),
        "median_AUC": round(d.AUC.median(), 3),
        "ln_ic50_range": f"{d.ln_ic50.min():.2f} to {d.ln_ic50.max():.2f}",
        "ln_ic50_sd": round(d.ln_ic50.std(), 3),
    })
panel = pd.DataFrame(rows)
panel.to_csv(f"{OUT}/table1_drug_panel.csv", index=False)

# ---------- 复孔可靠性：同系同药的独立重复筛选 ----------
def reliability(frame, key_cols, label):
    """成对重复的一致性：pooled within-pair SD + ICC(1)"""
    g = frame.groupby(key_cols)["ln_ic50"]
    agg = g.agg(["size", "mean", "std"]).reset_index()
    rep = agg[agg["size"] > 1]
    if rep.empty:
        return None
    # 组内合并标准差（等组内方差假设）
    ss_within = ((rep["std"] ** 2) * (rep["size"] - 1)).sum()
    df_within = (rep["size"] - 1).sum()
    pooled_sd = np.sqrt(ss_within / df_within)
    # ICC(1)：组间方差 / (组间方差 + 组内方差)
    k = rep["size"].mean()
    ms_between = ((rep["mean"] - rep["mean"].mean()) ** 2).sum() / (len(rep) - 1)
    var_between = max(ms_between - pooled_sd ** 2, 1e-9)
    icc = var_between / (var_between + pooled_sd ** 2)
    # 典型成对差异倍数
    typical_pair_ratio = float(np.exp(pooled_sd * np.sqrt(2)))
    return {
        "scheme": label,
        "key": " + ".join(key_cols),
        "n_comparisons": int(len(rep)),
        "n_curves_involved": int(rep["size"].sum()),
        "pooled_within_pair_sd_ln": round(float(pooled_sd), 3),
        "median_pair_sd_ln": round(float(rep["std"].median()), 3),
        "pct_pairs_sd_gt_1": round(100 * float((rep["std"] > 1).mean()), 1),
        "pct_pairs_ratio_gt_7.4x": round(100 * float((rep["std"] * np.sqrt(2) > 2).mean()), 1),
        "typical_pair_fold_ratio": round(typical_pair_ratio, 2),
        "icc_1": round(float(icc), 3),
    }

rel_rows = [
    reliability(crc, ["CELL_LINE_NAME", "DRUG_ID"], "strict: line + DRUG_ID (同一化合物批次)"),
    reliability(crc, ["CELL_LINE_NAME", "DRUG_NAME"], "loose: line + DRUG_NAME (可能跨批次)"),
    reliability(crc[crc.DRUG_NAME.isin(["5-Fluorouracil", "Oxaliplatin", "Irinotecan", "SN-38"])],
                ["CELL_LINE_NAME", "DRUG_NAME"], "chemo panel only: line + DRUG_NAME"),
]
rel = pd.DataFrame([r for r in rel_rows if r])
rel.to_csv(f"{OUT}/table2_replicate_reliability.csv", index=False)

# ---------- 汇总 ----------
summary = {
    "dataset": "GDSC2 (24Jul22)",
    "n_total_curves": int(len(df)),
    "n_coread_curves": int(len(crc)),
    "n_coread_lines": int(crc.CELL_LINE_NAME.nunique()),
    "n_coread_drugs": int(crc.DRUG_NAME.nunique()),
    "coread_pct_ic50_extrapolated_above_max_conc": round(100 * crc.above_range.mean(), 1),
    "coread_pct_ic50_extrapolated_below_min_conc": round(100 * crc.below_range.mean(), 1),
    "coread_pct_ic50_identifiable": round(100 * crc.identifiable.mean(), 1),
    "coread_median_RMSE": round(float(crc.RMSE.median()), 3),
    "coread_pct_RMSE_gt_0.15": round(100 * float((crc.RMSE > 0.15).mean()), 1),
    "no_ci_provided_in_public_file": True,
    "drugs_absent_that_matter_clinically": [
        "Cetuximab", "Panitumumab", "Capecitabine", "Trifluridine", "Regorafenib", "Encorafenib"
    ],
}
with open(f"{OUT}/01_label_forensics.json", "w", encoding="utf-8") as f:
    json.dump(summary, f, indent=2, ensure_ascii=False)

print(json.dumps(summary, indent=2, ensure_ascii=False))
print("\n=== TABLE 1: drug panel ===")
print(panel.to_string(index=False))
print("\n=== TABLE 2: replicate reliability ===")
print(rel.to_string(index=False))
