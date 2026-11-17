"""
04_extract_raw_coread.py
从 2.1 GB 原始板读数中抽取 COREAD 子集
-------------------------------------------------
原始文件是长格式：每一行 = 一个孔的原始荧光强度。
字段 TAG 区分 DMSO 对照 / 药物处理 / 阳性对照，INTENSITY 是未归一化的原始读数。
这正是"回到原始材料"所要求的那一层。

输出：data/processed/coread_raw_wells.parquet  results/04_raw_extract.json
"""
import json, os
import numpy as np
import pandas as pd

RAW = "data/raw/GDSC2_public_raw_data_24Jul22.csv"
os.makedirs("data/processed", exist_ok=True)

fit = pd.read_csv("data/raw/GDSC2_fitted_dose_response_24Jul22.csv")
coread = fit[fit.TCGA_DESC == "COREAD"]
lines = sorted(coread.CELL_LINE_NAME.unique())
sidm = set(coread.SANGER_MODEL_ID.dropna().unique())
print(f"COREAD 细胞系 {len(lines)} 株, SIDM {len(sidm)} 个")

usecols = ["BARCODE", "SCAN_ID", "SCAN_DATE", "CELL_LINE_NAME", "SANGER_MODEL_ID",
           "SEEDING_DENSITY", "DRUGSET_ID", "ASSAY", "DURATION", "POSITION",
           "TAG", "DRUG_ID", "CONC", "INTENSITY"]
dtypes = {"BARCODE": "category", "CELL_LINE_NAME": "object", "SANGER_MODEL_ID": "object",
          "DRUGSET_ID": "category", "ASSAY": "category", "DURATION": "category",
          "TAG": "category", "DRUG_ID": "Int32", "CONC": "float64",
          "INTENSITY": "float64", "SEEDING_DENSITY": "float64", "POSITION": "category"}

keep, tag_counts, n_rows = [], {}, 0
for chunk in pd.read_csv(RAW, usecols=usecols, dtype=dtypes,
                         chunksize=2_000_000, low_memory=False):
    n_rows += len(chunk)
    for t, c in chunk["TAG"].value_counts().items():
        tag_counts[str(t)] = tag_counts.get(str(t), 0) + int(c)
    sub = chunk[chunk.SANGER_MODEL_ID.isin(sidm)]
    if len(sub):
        keep.append(sub)

raw = pd.concat(keep, ignore_index=True)
print(f"扫描 {n_rows:,} 行 -> COREAD 保留 {len(raw):,} 孔")
print("TAG 分布(全库抽样统计):", dict(sorted(tag_counts.items(), key=lambda x: -x[1])[:8]))
print(f"\nCOREAD 子集: {raw.BARCODE.nunique()} 块板, {raw.CELL_LINE_NAME.nunique()} 株系, "
      f"{raw.DRUG_ID.nunique()} 个药, 浓度点 {raw.CONC.nunique()} 种")
print("TAG 分布(COREAD):", raw.TAG.value_counts().to_dict())

raw.to_parquet("data/processed/coread_raw_wells.parquet", index=False, compression="zstd")

json.dump({
    "total_rows_scanned": int(n_rows),
    "coread_wells": int(len(raw)),
    "coread_plates": int(raw.BARCODE.nunique()),
    "coread_lines": int(raw.CELL_LINE_NAME.nunique()),
    "coread_drugs": int(raw.DRUG_ID.nunique()),
    "tag_counts_global": tag_counts,
    "tag_counts_coread": {str(k): int(v) for k, v in raw.TAG.value_counts().items()},
    "n_positive_control_wells": int((raw.TAG == "positive control").sum()),
    "n_dmso_wells": int((raw.TAG == "DMSO").sum()),
}, open("results/04_raw_extract.json", "w", encoding="utf-8"), indent=2, ensure_ascii=False)
print("\n已写出 data/processed/coread_raw_wells.parquet")
