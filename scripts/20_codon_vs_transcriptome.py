"""
20_codon_vs_transcriptome.py
----------------------------
论文的方法学核心: 同一个药, 同一个样本集,
  (a) 只用 BRAF 状态(一个二值变量) 能解释多少 IC50 方差?
  (b) 用两万个转录本训练弹性网(嵌套CV) 能解释多少?
  (c) 转录本模型能不能自己把 BRAF 找出来?
  (d) 在 BRAF 已突变者内部, 转录组能否解释残余差异?

这一步决定了"基因 vs 机器学习"的定量对照, 也是本专刊主题的落点。
"""
import gzip
import re
import warnings
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.linear_model import ElasticNetCV
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import KFold, cross_val_predict

warnings.filterwarnings("ignore")
BASE = "data/raw"
DRUGS = ["Dabrafenib", "Selumetinib", "Trametinib", "PLX-4720",
         "Refametinib", "PD0325901", "SCH772984", "Trametinib"]


def norm(s):
    return re.sub(r"[^A-Z0-9]", "", str(s).upper())


def braf_set():
    m = pd.read_csv(f"{BASE}/mutations_latest.csv.gz")
    m = m[(m.gene_symbol == "BRAF") &
          (m.cancer_driver.astype(str).str.lower() == "true")]
    pm = m.protein_mutation.fillna("").astype(str)
    return {norm(x) for x in m[pm.str.contains("V600")].model_name.unique()}


def main():
    v600 = braf_set()
    # 表达谱表头与列映射
    with gzip.open(f"{BASE}/CCLE_RNAseq_genes_rpkm_20180929.gct.gz", "rt") as f:
        f.readline(); f.readline()
        samples = f.readline().rstrip("\n").split("\t")[2:]
    lk = {}
    for s in samples:
        lk.setdefault(norm(s.split("_", 1)[0]), []).append(s)
    lk = {k: v[0] for k, v in lk.items()}

    d = pd.read_csv(f"{BASE}/GDSC2_fitted_dose_response_24Jul22.csv")
    d["ic50"] = np.exp(d.LN_IC50)
    d["ident"] = (d.ic50 >= d.MIN_CONC) & (d.ic50 <= d.MAX_CONC)
    d = d[d.ident].copy()
    d["nm"] = d.CELL_LINE_NAME.map(norm)
    d["ccle"] = d.nm.map(lk)
    allcells = sorted(d.ccle.dropna().unique())
    print(f"可匹配表达谱的细胞系 {len(allcells)} 株; 读取表达矩阵 ...")

    raw = pd.read_csv(f"{BASE}/CCLE_RNAseq_genes_rpkm_20180929.gct.gz", sep="\t",
                      skiprows=2, usecols=["Name", "Description"] + allcells,
                      dtype={x: "float32" for x in allcells}, compression="gzip")
    raw = raw.dropna(subset=["Description"])
    expr = np.log2(raw.set_index("Description")[allcells]
                   .groupby(level=0).max() + 1)
    print(f"  表达矩阵 {expr.shape}")
    # 无监督高变基因预筛(不接触标签)
    expr = expr.loc[expr.var(axis=1).sort_values(ascending=False).index[:4000]]
    print(f"  高变基因筛至 {expr.shape[0]}")

    rows = []
    for drug in sorted(set(DRUGS)):
        s = d[d.DRUG_NAME == drug].dropna(subset=["ccle"]).drop_duplicates("ccle")
        if len(s) < 40:
            continue
        cells = [c for c in s.ccle if c in expr.columns]
        if len(cells) < 40:
            continue
        sub = s[s.ccle.isin(cells)].set_index("ccle").loc[cells]
        y = sub.LN_IC50.values.astype(float)
        isv = np.array([c.split("_")[0] in v600 or norm(c.split("_")[0]) in v600
                        for c in cells])

        # (a) BRAF 单独
        if isv.sum() >= 5 and (~isv).sum() >= 20:
            r2_braf = stats.pearsonr(isv.astype(float), y)[0] ** 2
            d_braf = (y[isv].mean() - y[~isv].mean()) / np.sqrt(
                ((isv.sum() - 1) * y[isv].var(ddof=1)
                 + ((~isv).sum() - 1) * y[~isv].var(ddof=1))
                / (len(y) - 2))
        else:
            r2_braf, d_braf = np.nan, np.nan

        # (b) 转录组弹性网, 嵌套CV
        X = expr[cells].T.values.astype(np.float32)
        pipe = Pipeline([("sc", StandardScaler()),
                         ("en", ElasticNetCV(l1_ratio=[0.5],
                                             alphas=np.logspace(-2, 1.5, 15),
                                             cv=4, max_iter=500, n_jobs=1,
                                             random_state=0))])
        pred = cross_val_predict(pipe, X, y, cv=KFold(5, shuffle=True,
                                                      random_state=0))
        r2_expr = 1 - ((y - pred) ** 2).sum() / ((y - y.mean()) ** 2).sum()
        r_expr = stats.pearsonr(pred, y)[0]

        # (c) 转录组能否预测 BRAF 状态
        auc_braf = np.nan
        if isv.sum() >= 5 and (~isv).sum() >= 20:
            from sklearn.linear_model import LogisticRegressionCV
            from sklearn.model_selection import StratifiedKFold
            from sklearn.metrics import roc_auc_score
            lp = Pipeline([("sc", StandardScaler()),
                           ("lr", LogisticRegressionCV(Cs=5, cv=3, max_iter=500,
                                                       solver="liblinear"))])
            pp = cross_val_predict(lp, X, isv.astype(int),
                                   cv=StratifiedKFold(5, shuffle=True,
                                                      random_state=0),
                                   method="predict_proba")[:, 1]
            auc_braf = roc_auc_score(isv.astype(int), pp)

        # (d) BRAF 突变者内部的残余可解释性
        r2_within = np.nan
        if isv.sum() >= 15:
            Xw, yw = X[isv], y[isv]
            pw = cross_val_predict(pipe, Xw, yw,
                                   cv=KFold(5, shuffle=True, random_state=0))
            r2_within = 1 - ((yw - pw) ** 2).sum() / ((yw - yw.mean()) ** 2).sum()

        rows.append(dict(drug=drug, n=len(y), n_v600=int(isv.sum()),
                         r2_braf=r2_braf, d_braf=d_braf,
                         r2_expr=r2_expr, r_expr=r_expr,
                         auc_expr_for_braf=auc_braf, r2_expr_within_v600=r2_within))
        print(f"  {drug:12s} n={len(y):4d} v600={isv.sum():3d} | "
              f"BRAF R2={r2_braf:.3f} | 转录组 R2={r2_expr:+.3f} | "
              f"转录组判BRAF AUC={auc_braf if np.isnan(auc_braf) else round(auc_braf,3)}")

    r = pd.DataFrame(rows).sort_values("r2_braf", ascending=False)
    r.to_csv("results/20_codon_vs_transcriptome.csv", index=False)
    print("\n" + "=" * 78)
    print("一个密码子 vs 两万个转录本")
    print("=" * 78)
    print(r.to_string(index=False, float_format=lambda v: f"{v:.3f}"))
    m = r.dropna(subset=["r2_braf"])
    if len(m):
        print(f"\n  BRAF 状态解释方差 中位 = {m.r2_braf.median():.3f}")
        print(f"  转录组模型解释方差 中位 = {m.r2_expr.median():+.3f}")
        print(f"  比值 = {m.r2_braf.median()/max(m.r2_expr.median(),1e-6):.1f} 倍")
        print(f"  转录组能判出 BRAF 状态的 AUC 中位 = {m.auc_expr_for_braf.median():.3f}")


if __name__ == "__main__":
    main()
