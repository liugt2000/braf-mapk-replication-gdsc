"""
11_well_level_noise.py
在孔一级直接测量噪声地板
-------------------------------------------------
前面所有"标签不可靠"的说法都依赖拟合结果的间接推断。
这里回到 189 万个原始孔，直接测：
 (a) 同一块板上、同一药物、同一浓度的复孔变异（技术噪声）
 (b) 同一浓度下不同细胞系之间的变异（生物学信号）
信噪比 = (b) / (a)。若接近 1，则模型不可能学到东西。

输出：results/11_well_noise.json, results/table7_noise_floor.csv
"""
import json
import numpy as np
import pandas as pd

raw = pd.read_parquet("data/processed/coread_raw_wells.parquet")
raw = raw[raw.TAG.astype(str) != "UN-USED"]
print(f"原始孔 {len(raw):,}  板 {raw.BARCODE.nunique()}  系 {raw.CELL_LINE_NAME.nunique()}")

# ---------- (a) 复孔变异 ----------
trt = raw[raw.DRUG_ID.notna() & raw.CONC.notna()]
rep = (trt.groupby(["BARCODE", "DRUG_ID", "CONC", "CELL_LINE_NAME"])["INTENSITY"]
       .agg(n="size", mean="mean", std="std").reset_index())
rep = rep[(rep.n >= 2) & (rep["mean"] > 0)]
rep["cv"] = rep["std"] / rep["mean"]
print(f"\n复孔组 {len(rep):,}（每组 >=2 孔，中位孔数 {int(rep.n.median())}）")
print(f"  复孔 CV 中位数 = {rep.cv.median():.3f}   75 分位 = {rep.cv.quantile(.75):.3f}"
      f"   90 分位 = {rep.cv.quantile(.9):.3f}")

dmso = (raw[raw.TAG.astype(str) == "DMSO"]
        .groupby("BARCODE")["INTENSITY"].agg(n="size", mean="mean", std="std"))
dmso = dmso[(dmso.n >= 3) & (dmso["mean"] > 0)]
dmso["cv"] = dmso["std"] / dmso["mean"]
print(f"  DMSO 对照板内 CV 中位数 = {dmso.cv.median():.3f}（{len(dmso)} 块板）")

# ---------- (b) 同一浓度下不同细胞系之间的变异 ----------
# 用归一化存活率：每孔 / 该板 DMSO 中位数
ctrl = (raw[raw.TAG.astype(str) == "DMSO"]
        .groupby("BARCODE")["INTENSITY"].median().rename("ctrl"))
trt2 = trt.merge(ctrl, on="BARCODE", how="inner")
trt2 = trt2[trt2.ctrl > 0]
trt2["viab"] = trt2.INTENSITY / trt2.ctrl

# 选每块板最常见的药物浓度布局：按 (DRUG_ID, CONC) 聚合到细胞系均值
line_conc = (trt2.groupby(["DRUG_ID", "CONC", "CELL_LINE_NAME"])["viab"]
             .mean().reset_index())
g = line_conc.groupby(["DRUG_ID", "CONC"]).agg(
    n_lines=("CELL_LINE_NAME", "nunique"),
    between_sd=("viab", "std"),
    between_mean=("viab", "mean")).reset_index()
g = g[g.n_lines >= 20]

# 同 (DRUG_ID, CONC) 的复孔 CV 作为技术噪声
tech = (trt2.groupby(["DRUG_ID", "CONC", "CELL_LINE_NAME"])["viab"]
        .agg(n="size", mean="mean", std="std").reset_index())
tech = tech[(tech.n >= 2) & (tech["mean"] > 0)]
tech["cv"] = tech["std"] / tech["mean"]
tech_g = tech.groupby(["DRUG_ID", "CONC"])["cv"].median().rename("tech_cv").reset_index()

m = g.merge(tech_g, on=["DRUG_ID", "CONC"], how="inner")
m["snr"] = m.between_sd / m.tech_cv

print("\n=== 信噪比：细胞系间 SD / 复孔 CV ===")
print(f"  参与比较的 (药, 浓度) 组合: {len(m)}")
print(f"  细胞系间 SD 中位数   : {m.between_sd.median():.3f}")
print(f"  复孔 CV 中位数       : {m.tech_cv.median():.3f}")
print(f"  信噪比 中位          : {m.snr.median():.2f}")
print(f"  信噪比 < 1 的比例     : {100*(m.snr<1).mean():.1f}%")
print(f"  信噪比 < 2 的比例     : {100*(m.snr<2).mean():.1f}%")

# ---------- 逐药 ----------
drugs = {1003: "Camptothecin"}
# 从拟合表取药名
fit = pd.read_csv("data/raw/GDSC2_fitted_dose_response_24Jul22.csv")
dmap = fit.drop_duplicates("DRUG_ID").set_index("DRUG_ID")["DRUG_NAME"].to_dict()
m["drug"] = m.DRUG_ID.map(dmap)
sel = ["5-Fluorouracil", "Oxaliplatin", "Irinotecan", "SN-38",
       "Trametinib", "Selumetinib", "Dabrafenib", "Fulvestrant"]
sub = m[m.drug.isin(sel)]
per = (sub.groupby("drug").agg(n_conc=("snr", "size"),
                               between_sd=("between_sd", "median"),
                               tech_cv=("tech_cv", "median"),
                               snr_median=("snr", "median"),
                               pct_snr_lt1=("snr", lambda s: round(100 * (s < 1).mean(), 1)))
       .reset_index())
per.to_csv("results/table7_noise_floor.csv", index=False)
print("\n=== TABLE 7: 逐药噪声地板 ===")
print(per.to_string(index=False))
if len(sub):
    ex = sub[sub.drug == "5-Fluorouracil"].sort_values("CONC")
    print(f"\n5-FU 逐浓度信噪比（前 8 个）:")
    print(ex[["CONC", "n_lines", "between_sd", "tech_cv", "snr"]].head(8).to_string(index=False))

json.dump({
    "n_wells_used": int(len(raw)),
    "n_replicate_groups": int(len(rep)),
    "replicate_cv_median": round(float(rep.cv.median()), 3),
    "replicate_cv_p90": round(float(rep.cv.quantile(.9)), 3),
    "dmso_plate_cv_median": round(float(dmso.cv.median()), 3),
    "n_drug_conc_comparisons": int(len(m)),
    "between_line_sd_median": round(float(m.between_sd.median()), 3),
    "tech_cv_median": round(float(m.tech_cv.median()), 3),
    "snr_median": round(float(m.snr.median()), 2),
    "pct_snr_lt1": round(100 * float((m.snr < 1).mean()), 1),
    "pct_snr_lt2": round(100 * float((m.snr < 2).mean()), 1),
    "per_drug": per.to_dict("records"),
}, open("results/11_well_noise.json", "w", encoding="utf-8"), indent=2, ensure_ascii=False)
