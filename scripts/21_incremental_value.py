"""
21_incremental_value.py
-----------------------
转录组能解释的方差略高于 BRAF 状态, 但这不意味着它带来新信息 ——
因为 BRAF 突变会在转录组上留下很强的 MAPK 输出印记(判 BRAF 的 AUC 约 0.91)。
所以真正的问题是: 条件于 BRAF 状态, 转录组还能加多少?

做法(严格无损): 先在训练折内用 BRAF 状态拟合 y 得到残差, 再用转录组在
**内层 CV** 里预测残差; 报告外层 CV 的增量 R2。若增量 ≈ 0 或为负,
说明转录组只是突变状态的替代变量。
"""
import gzip
import re
import warnings
import numpy as np
import pandas as pd
from sklearn.linear_model import ElasticNetCV
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import KFold

warnings.filterwarnings("ignore")
BASE = "data/raw"
DRUGS = ["Dabrafenib", "Selumetinib", "Trametinib", "Refametinib",
         "PD0325901", "SCH772984"]


def norm(s):
    return re.sub(r"[^A-Z0-9]", "", str(s).upper())


def braf_set():
    m = pd.read_csv(f"{BASE}/mutations_latest.csv.gz")
    m = m[(m.gene_symbol == "BRAF") &
          (m.cancer_driver.astype(str).str.lower() == "true")]
    pm = m.protein_mutation.fillna("").astype(str)
    return {norm(x) for x in m[pm.str.contains("V600")].model_name.unique()}


def r2(y, p):
    return 1 - ((y - p) ** 2).sum() / ((y - y.mean()) ** 2).sum()


def main():
    v600 = braf_set()
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
    cells_all = sorted(d.ccle.dropna().unique())

    raw = pd.read_csv(f"{BASE}/CCLE_RNAseq_genes_rpkm_20180929.gct.gz", sep="\t",
                      skiprows=2, usecols=["Name", "Description"] + cells_all,
                      dtype={x: "float32" for x in cells_all}, compression="gzip")
    raw = raw.dropna(subset=["Description"])
    expr = np.log2(raw.set_index("Description")[cells_all]
                   .groupby(level=0).max() + 1)
    expr = expr.loc[expr.var(axis=1).sort_values(ascending=False).index[:4000]]
    print(f"表达矩阵 {expr.shape}")

    rows = []
    for drug in DRUGS:
        s = d[d.DRUG_NAME == drug].dropna(subset=["ccle"]).drop_duplicates("ccle")
        cells = [c for c in s.ccle if c in expr.columns]
        if len(cells) < 40:
            continue
        sub = s[s.ccle.isin(cells)].set_index("ccle").loc[cells]
        y = sub.LN_IC50.values.astype(float)
        braf = np.array([1.0 if norm(c.split("_")[0]) in v600 else 0.0
                         for c in cells])
        X = expr[cells].T.values.astype(np.float32)

        outer = KFold(5, shuffle=True, random_state=0)
        pred_braf = np.zeros_like(y)
        pred_both = np.zeros_like(y)
        for tr, te in outer.split(X):
            # --- 仅 BRAF ---
            b0 = braf[tr].mean()
            m1 = np.polyfit(braf[tr], y[tr], 1)
            pred_braf[te] = np.polyval(m1, braf[te])

            # --- BRAF + 转录组(残差法) ---
            res_tr = y[tr] - np.polyval(m1, braf[tr])
            pipe = Pipeline([("sc", StandardScaler()),
                             ("en", ElasticNetCV(l1_ratio=[0.5],
                                                 alphas=np.logspace(-2, 1, 12),
                                                 cv=4, max_iter=400, n_jobs=1,
                                                 random_state=0))])
            pipe.fit(X[tr], res_tr)
            pred_both[te] = np.polyval(m1, braf[te]) + pipe.predict(X[te])

        r2_b = r2(y, pred_braf)
        r2_both = r2(y, pred_both)
        rows.append(dict(drug=drug, n=len(y), n_v600=int(braf.sum()),
                         r2_braf=r2_b, r2_braf_plus_expr=r2_both,
                         delta_r2=r2_both - r2_b))
        print(f"  {drug:12s} n={len(y):4d} | BRAF R2={r2_b:.3f} | "
              f"BRAF+转录组 R2={r2_both:.3f} | 增量={r2_both-r2_b:+.3f}")

    r = pd.DataFrame(rows).sort_values("delta_r2", ascending=False)
    r.to_csv("results/21_incremental_value.csv", index=False)
    print("\n" + "=" * 70)
    print("条件于 BRAF 状态后, 转录组的增量价值")
    print("=" * 70)
    print(r.to_string(index=False, float_format=lambda v: f"{v:.3f}"))
    print(f"\n  增量 R2 中位 = {r.delta_r2.median():+.3f}   "
          f"为正的药物数 = {(r.delta_r2>0).sum()}/{len(r)}")
    print("  解读: 增量≈0 或为负 => 转录组模型只是突变状态的替代变量")


if __name__ == "__main__":
    main()
