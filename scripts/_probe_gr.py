"""快速探针: GR 重推导的前置数据核查(不写正式结果)"""
import re
import pandas as pd

BASE = r"D:\AI_work\workbuddy\2026-11-16-16-13-17\data\raw"
RES = r"D:\AI_work\workbuddy\2026-11-16-16-13-17\results"


def norm(s):
    return re.sub(r"[^A-Z0-9]", "", str(s).upper())


# 1. 目标系: GDSC2 v600 有双药可测的系
A = pd.read_csv(rf"{RES}\22_crc_exception.csv")
g2 = A[(A.campaign == "GDSC2") & A.v600]
print("GDSC2 v600 lines:", len(g2), " COREAD:", (g2.TCGA_DESC == "COREAD").sum())
print(g2[["nm", "TCGA_DESC", "Dabrafenib", "Trametinib"]]
      .sort_values("TCGA_DESC").to_string(index=False))

# 2. 倍增时间覆盖
gr = pd.read_csv(rf"{BASE}\growth_rate_20220907.csv")
gr["nm"] = gr.model_name.map(norm)
gr["dt"] = pd.to_numeric(gr.doubling_time_hours, errors="coerce")
cov = g2[["nm", "TCGA_DESC"]].merge(gr[["nm", "dt"]], on="nm", how="left")
print("\nDT coverage:", cov.dt.notna().sum(), "/", len(cov))
print(cov.to_string(index=False))
print("\nDT missing lines:", cov[cov.dt.isna()].nm.tolist())

# 3. 药物 ID
fit_cols = ["DRUG_NAME", "DRUG_ID"]
f = pd.read_csv(rf"{BASE}\GDSC2_fitted_dose_response_24Jul22.csv",
                usecols=["DRUG_NAME", "DRUG_ID", "MIN_CONC", "MAX_CONC"])
for drug in ["Dabrafenib", "Trametinib"]:
    sub = f[f.DRUG_NAME == drug]
    print(f"\n{drug}: DRUG_ID={sub.DRUG_ID.unique()}, "
          f"MIN={sub.MIN_CONC.median()}, MAX={sub.MAX_CONC.median()}")

# 4. raw 文件头部结构: DURATION / TAG / DMSO 标记
raw = pd.read_csv(rf"{BASE}\GDSC2_public_raw_data_24Jul22.csv",
                  nrows=200000, low_memory=False)
print("\nDURATION values:", raw.DURATION.unique())
print("DMSO rows in first 200k:", raw.DRUG_ID.isna().sum())
print("TAG samples (drug wells):",
      raw[raw.DRUG_ID.notna()].TAG.dropna().unique()[:5])
print("TAG samples (DMSO wells):",
      raw[raw.DRUG_ID.isna()].TAG.dropna().unique()[:5])
print("CONC samples:", raw.CONC.dropna().unique()[:5])
print("BARCODE dup example:", raw.BARCODE.nunique(), "plates in first 200k rows")
