"""
10_auc_vs_ic50.py
IC50 vs AUC：换一个目标量能不能救回来？
-------------------------------------------------
IC50 在 73% 的情况下不存在（被外推）。GDSC 同时提供 AUC ——
它是在实测浓度区间上算的曲线下面积，对每条曲线都有定义、有界、不需要外推。

检验三件事：
1. 跨活动重测信度：AUC 是否比 LN_IC50 稳定？
2. 可识别性：AUC 是否存在"不存在"的问题？（不存在——它总有定义）
3. 模型可预测性：以 AUC 为目标，是否更多药物能打败置换零？

输出：results/10_auc_vs_ic50.json, results/table6_auc_vs_ic50.csv
"""
import gzip, json, re, warnings
import numpy as np
import pandas as pd
from sklearn.linear_model import ElasticNetCV
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import KFold, cross_val_predict
from sklearn.metrics import r2_score
warnings.filterwarnings("ignore")
rng = np.random.default_rng(3)


def norm(s):
    return re.sub(r"[^A-Z0-9]", "", str(s).upper())


# ---------- 1. 跨活动信度：AUC vs IC50 ----------
g1 = pd.read_csv("data/raw/GDSC1_fitted_dose_response_24Jul22.csv")
g2 = pd.read_csv("data/raw/GDSC2_fitted_dose_response_24Jul22.csv")
c1, c2 = g1[g1.TCGA_DESC == "COREAD"], g2[g2.TCGA_DESC == "COREAD"]

def agg(df, v, tag):
    a = df.groupby(["CELL_LINE_NAME", "DRUG_NAME"], as_index=False)[v].mean()
    return a.rename(columns={v: f"{v}_{tag}"})

rows = []
for metric in ["LN_IC50", "AUC"]:
    m = agg(c1, metric, "1").merge(agg(c2, metric, "2"), on=["CELL_LINE_NAME", "DRUG_NAME"])
    x, y = m[f"{metric}_1"].values, m[f"{metric}_2"].values
    r = np.corrcoef(x, y)[0, 1]
    d = np.abs(y - x)
    rows.append({"metric": metric, "n_pairs": len(m),
                 "pearson_r": round(float(r), 3),
                 "r_squared": round(float(r ** 2), 3),
                 "median_abs_diff": round(float(np.median(d)), 4)})
    # 逐药
    for drug in ["5-Fluorouracil", "Oxaliplatin", "Irinotecan", "SN-38",
                 "Trametinib", "Selumetinib", "Dabrafenib", "Fulvestrant"]:
        s = m[m.DRUG_NAME == drug]
        if len(s) < 10:
            continue
        rr = float(np.corrcoef(s[f"{metric}_1"], s[f"{metric}_2"])[0, 1])
        rows.append({"metric": metric, "n_pairs": len(s), "drug": drug,
                     "pearson_r": round(rr, 3), "r_squared": round(rr ** 2, 3)})
t = pd.DataFrame(rows)
print("=== 跨活动重测信度：IC50 vs AUC ===")
print(t[t.drug.isna()].to_string(index=False))
print()
piv = t.dropna(subset=["drug"]).pivot_table(index="drug", columns="metric",
                                            values="pearson_r")
piv["AUC_minus_IC50"] = (piv["AUC"] - piv["LN_IC50"]).round(3)
print(piv.round(3).to_string())

# ---------- 2. 模型：AUC 为目标 ----------
with gzip.open("data/raw/CCLE_RNAseq_genes_rpkm_20180929.gct.gz", "rt") as f:
    f.readline(); f.readline()
    S = f.readline().rstrip("\n").split("\t")[2:]
lk = {}
for s in S:
    lk.setdefault(norm(s.split("_", 1)[0]), []).append(s)
lk = {k: v[0] for k, v in lk.items()}

crc = c2.copy()
crc["ccle"] = crc.CELL_LINE_NAME.map(lambda n: lk.get(norm(n)))
need = sorted(crc.ccle.dropna().unique())
raw = pd.read_csv("data/raw/CCLE_RNAseq_genes_rpkm_20180929.gct.gz", sep="\t",
                  skiprows=2, usecols=["Name", "Description"] + need,
                  dtype={c: "float32" for c in need}, compression="gzip")
raw = raw.dropna(subset=["Description"])
expr = np.log2(raw.set_index("Description")[need].groupby(level=0).max() + 1)
expr = expr.loc[expr.var(axis=1).sort_values(ascending=False).index[:2000]]


def model():
    return Pipeline([("sc", StandardScaler()),
                     ("en", ElasticNetCV(l1_ratio=[0.5], cv=4,
                                         alphas=np.logspace(-2, 1.5, 15),
                                         max_iter=1000, random_state=0, n_jobs=1))])


DRUGS = ["5-Fluorouracil", "Oxaliplatin", "Irinotecan", "SN-38",
         "Trametinib", "Selumetinib", "Dabrafenib", "Fulvestrant"]
res = []
for drug in DRUGS:
    d = crc[crc.DRUG_NAME == drug].dropna(subset=["ccle"]).drop_duplicates("ccle")
    if len(d) < 20:
        continue
    X = expr[d.ccle.values].T.values.astype(np.float32)
    cv = KFold(5, shuffle=True, random_state=0)
    row = {"drug": drug, "n": len(d)}
    for target in ["LN_IC50", "AUC"]:
        y = d[target].values.astype(float)
        pred = cross_val_predict(model(), X, y, cv=cv)
        r2 = r2_score(y, pred)
        null = np.array([r2_score(yp, cross_val_predict(model(), X, yp, cv=cv))
                         for yp in (rng.permutation(y) for _ in range(12))])
        row[f"{target}_r2"] = round(float(r2), 3)
        row[f"{target}_null95"] = round(float(np.percentile(null, 95)), 3)
        row[f"{target}_beats"] = bool(r2 > np.percentile(null, 95))
    res.append(row)
    print(f"  {drug:16s} IC50 R2={row['LN_IC50_r2']:+.3f} ({'beat' if row['LN_IC50_beats'] else 'FAIL'})"
          f"   AUC R2={row['AUC_r2']:+.3f} ({'beat' if row['AUC_beats'] else 'FAIL'})")

rt = pd.DataFrame(res)
rt.to_csv("results/table6_auc_vs_ic50.csv", index=False)
print("\n=== TABLE 6: 目标量对模型可预测性的影响 ===")
print(rt.to_string(index=False))
print(f"\n打败随机的药物数：IC50 {int(rt.LN_IC50_beats.sum())}/{len(rt)}  "
      f"AUC {int(rt.AUC_beats.sum())}/{len(rt)}")

json.dump({"reliability": t.to_dict("records"),
           "model": rt.to_dict("records")},
          open("results/10_auc_vs_ic50.json", "w", encoding="utf-8"),
          indent=2, ensure_ascii=False)
