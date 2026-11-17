"""
06_decisive_test.py
5-FU 模型的预测力，是否由"外推的 IC50"撑起来？
-------------------------------------------------
5-FU 是唯一同时满足两个条件的药：模型显著优于随机（R2=0.364），
且标签有可观的重测信度（rho=0.667）。

但如果这 0.364 主要来自那些 IC50 被外推到最高测试浓度之外的细胞系
（即"根本没被抑制"的那些），那模型学到的是"这条曲线没做出来"，
而不是"这个细胞耐药"。

做法：按 IC50 是否可识别（落在实测浓度区间内）分层，比较模型 R2，
      并给出 bootstrap 置信区间。
输出：results/06_decisive_test.json
"""
import gzip, json, re, warnings, gc
import numpy as np
import pandas as pd
from sklearn.linear_model import ElasticNetCV
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import KFold, cross_val_predict
from sklearn.metrics import r2_score
warnings.filterwarnings("ignore")
rng = np.random.default_rng(7)


def norm(s):
    return re.sub(r"[^A-Z0-9]", "", str(s).upper())


with gzip.open("data/raw/CCLE_RNAseq_genes_rpkm_20180929.gct.gz", "rt") as f:
    f.readline(); f.readline()
    samples = f.readline().rstrip("\n").split("\t")[2:]
lk = {}
for s in samples:
    lk.setdefault(norm(s.split("_", 1)[0]), []).append(s)
lk = {k: v[0] for k, v in lk.items()}

fit = pd.read_csv("data/raw/GDSC2_fitted_dose_response_24Jul22.csv")
crc = fit[fit.TCGA_DESC == "COREAD"].copy()
crc["ic50"] = np.exp(crc.LN_IC50)
crc["ident"] = (crc.ic50 >= crc.MIN_CONC) & (crc.ic50 <= crc.MAX_CONC)
crc["ccle"] = crc.CELL_LINE_NAME.map(lambda n: lk.get(norm(n)))

need = sorted(crc.ccle.dropna().unique())
raw = pd.read_csv("data/raw/CCLE_RNAseq_genes_rpkm_20180929.gct.gz", sep="\t",
                  skiprows=2, usecols=["Name", "Description"] + need,
                  dtype={c: "float32" for c in need}, compression="gzip")
raw = raw.dropna(subset=["Description"])
expr = np.log2(raw.set_index("Description")[need].groupby(level=0).max() + 1)
var = expr.var(axis=1).sort_values(ascending=False)
expr = expr.loc[var.index[:2000]]
del raw; gc.collect()


def model():
    return Pipeline([("sc", StandardScaler()),
                     ("en", ElasticNetCV(l1_ratio=[0.5], cv=4,
                                         alphas=np.logspace(-2, 1.5, 15),
                                         max_iter=1000, random_state=0, n_jobs=1))])


def evaluate(sub, label, n_boot=300):
    d = sub.dropna(subset=["ccle"]).drop_duplicates("ccle")
    if len(d) < 12:
        return {"stratum": label, "n_lines": len(d), "note": "n<12, 不评估"}
    X = expr[d.ccle.values].T.values.astype(np.float32)
    y = d.LN_IC50.values.astype(np.float64)
    cv = KFold(5, shuffle=True, random_state=0)
    pred = cross_val_predict(model(), X, y, cv=cv)
    r2 = r2_score(y, pred)
    sp = float(pd.Series(pred).corr(pd.Series(y), method="spearman"))
    # 置换零分布
    null = np.array([r2_score(yp, cross_val_predict(model(), X, yp, cv=cv))
                     for yp in (rng.permutation(y) for _ in range(15))])
    # bootstrap CI（对样本重采样，在固定预测上算，成本低）
    boots = []
    for _ in range(n_boot):
        idx = rng.integers(0, len(y), len(y))
        if len(np.unique(y[idx])) < 3:
            continue
        try:
            boots.append(r2_score(y[idx], pred[idx]))
        except Exception:
            pass
    lo, hi = (np.percentile(boots, [2.5, 97.5]) if boots else (np.nan, np.nan))
    return {"stratum": label, "n_lines": len(d),
            "model_r2": round(float(r2), 3),
            "r2_ci95": [round(float(lo), 3), round(float(hi), 3)],
            "spearman": round(sp, 3),
            "null_r2_p95": round(float(np.percentile(null, 95)), 3),
            "beats_null": bool(r2 > np.percentile(null, 95)),
            "pct_identifiable_in_subset": round(100 * float(d.ident.mean()), 1)}


rows = []
for drug in ["5-Fluorouracil", "Oxaliplatin", "Irinotecan"]:
    d = crc[crc.DRUG_NAME == drug].groupby(["CELL_LINE_NAME"], as_index=False).agg(
        LN_IC50=("LN_IC50", "mean"), ident=("ident", "all"), ccle=("ccle", "first"))
    rows.append(evaluate(d, f"{drug} | ALL lines"))
    rows.append(evaluate(d[d.ident], f"{drug} | 仅 IC50 可识别"))
    print(f"{drug}: all n={len(d)}, identifiable n={int(d.ident.sum())}")

tab = pd.DataFrame(rows)
tab.to_csv("results/06_decisive_test.csv", index=False)
print("\n=== 决定性检验：模型预测力 vs 标签可识别性 ===")
print(tab.to_string(index=False))
json.dump(tab.to_dict("records"),
          open("results/06_decisive_test.json", "w", encoding="utf-8"),
          indent=2, ensure_ascii=False)
