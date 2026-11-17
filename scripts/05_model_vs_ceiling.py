"""
05_model_vs_ceiling.py
模型性能 vs 标签可靠性天花板
-------------------------------------------------
核心检验：以 GDSC2 IC50 为标签、在 CRC 内部做交叉验证的模型，其可解释方差 R^2
在数学上不可能超过该标签的重测信度 rho（GDSC1 与 GDSC2 对同一系同一药的相关系数）：
      X1 = T + e1,  X2 = T + e2,  rho = corr(X1,X2) = var(T)/var(X)
      完美模型 M = T 时 corr(M, X1) = sqrt(rho)  =>  R^2 <= rho
若模型 R^2 系统性接近或超过 rho，说明评估存在泄漏，或它只是在拟合同一次测量的噪声。

同时给出：置换零分布（模型能否打败随机）、增殖基线（模型是否只是在预测增殖率）。
输出：results/05_model_vs_ceiling.csv / .json
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
rng = np.random.default_rng(42)


def norm(s):
    return re.sub(r"[^A-Z0-9]", "", str(s).upper())


# ---------- 细胞系名匹配 ----------
with gzip.open("data/raw/CCLE_RNAseq_genes_rpkm_20180929.gct.gz", "rt") as f:
    f.readline(); f.readline()
    samples = f.readline().rstrip("\n").split("\t")[2:]
ccle_lookup = {}
for s in samples:
    ccle_lookup.setdefault(norm(s.split("_", 1)[0]), []).append(s)
ccle_lookup = {k: v[0] for k, v in ccle_lookup.items()}

fit = pd.read_csv("data/raw/GDSC2_fitted_dose_response_24Jul22.csv")
crc = fit[fit.TCGA_DESC == "COREAD"]
need = sorted({ccle_lookup.get(norm(n)) for n in crc.CELL_LINE_NAME.unique()} - {None})
print(f"需要 {len(need)} 个 CCLE 细胞系列")

# ---------- 表达谱：只读需要的列，用基因名做索引 ----------
print("读取表达谱 ...")
raw = pd.read_csv("data/raw/CCLE_RNAseq_genes_rpkm_20180929.gct.gz",
                  sep="\t", skiprows=2, usecols=["Name", "Description"] + need,
                  dtype={c: "float32" for c in need}, compression="gzip")
raw = raw.dropna(subset=["Description"])
expr = raw.set_index("Description")[need]
expr = expr.groupby(level=0).max()          # 同名基因取最大
expr = np.log2(expr + 1)
print("  表达矩阵:", expr.shape)

expr_full = expr.copy()   # 保留全基因矩阵，供增殖基线使用
# 无监督高变基因预筛（不接触标签，不构成泄漏）
var = expr.var(axis=1).sort_values(ascending=False)
TOP = 2000
genes = var.index[:TOP].tolist()
expr = expr.loc[genes]
print(f"  高变基因筛至 {TOP}")

# ---------- 标签 ----------
lab = (crc.groupby(["CELL_LINE_NAME", "DRUG_NAME"], as_index=False)["LN_IC50"]
       .mean().rename(columns={"LN_IC50": "ln_ic50"}))
lab["ccle"] = lab.CELL_LINE_NAME.map(lambda n: ccle_lookup.get(norm(n)))

# ---------- 天花板 ----------
pairs = pd.read_csv("results/coread_pairs_with_identifiability.csv")
ceiling = {}
for drug, g in pairs.groupby("DRUG_NAME"):
    if len(g) >= 8:
        ceiling[drug] = float(np.corrcoef(g.ln_g1, g.ln_g2)[0, 1])

PROLIF = [g for g in ["MKI67", "PCNA", "CCNB1", "CDK1", "TOP2A", "BUB1", "AURKA",
                      "PLK1", "CDC20", "CCNA2", "RRM2", "TYMS", "MCM2", "MCM6"]
          if g in expr_full.index]
print("增殖基线基因:", PROLIF)


def model():
    return Pipeline([
        ("sc", StandardScaler()),
        ("en", ElasticNetCV(l1_ratio=[0.5], cv=4, alphas=np.logspace(-2, 1.5, 15),
                            max_iter=1000, random_state=0, n_jobs=1)),
    ])


DRUGS = ["5-Fluorouracil", "Oxaliplatin", "Irinotecan", "SN-38",
         "Trametinib", "Selumetinib", "Dabrafenib", "Fulvestrant"]

rows = []
for drug in DRUGS:
    d = lab[(lab.DRUG_NAME == drug) & lab.ccle.notna()].drop_duplicates("ccle")
    if len(d) < 20:
        print(f"跳过 {drug}: 匹配株系仅 {len(d)}")
        continue
    X = expr[d.ccle.values].T.values.astype(np.float32)
    y = d.ln_ic50.values.astype(np.float64)
    cv = KFold(5, shuffle=True, random_state=0)

    pred = cross_val_predict(model(), X, y, cv=cv)
    r2 = r2_score(y, pred)
    rho_s = float(pd.Series(pred).corr(pd.Series(y), method="spearman"))

    null = []
    for _ in range(10):
        yp = rng.permutation(y)
        null.append(r2_score(yp, cross_val_predict(model(), X, yp, cv=cv)))
    null = np.array(null)

    r2_prolif = None
    if PROLIF:
        Xp = expr_full.loc[PROLIF, d.ccle.values].T.values.astype(np.float32).mean(1, keepdims=True)
        r2_prolif = float(r2_score(y, cross_val_predict(model(), Xp, y, cv=cv)))

    gc.collect()
    ceil = ceiling.get(drug, np.nan)
    rows.append({
        "drug": drug, "n_lines": len(d),
        "model_cv_r2": round(float(r2), 3),
        "model_spearman": round(rho_s, 3),
        "null_r2_p95": round(float(np.percentile(null, 95)), 3),
        "null_r2_median": round(float(np.median(null)), 3),
        "beats_null": bool(r2 > np.percentile(null, 95)),
        "prolif_baseline_r2": round(r2_prolif, 3) if r2_prolif is not None else None,
        "ceiling_rho": round(ceil, 3) if np.isfinite(ceil) else None,
        "r2_minus_ceiling": round(float(r2 - ceil), 3) if np.isfinite(ceil) else None,
    })
    pd.DataFrame(rows).to_csv("results/05_model_vs_ceiling.csv", index=False)
    print(f"  {drug:16s} n={len(d):2d} R2={r2:+.3f} null95={np.percentile(null,95):+.3f} "
          f"ceil={ceil if np.isfinite(ceil) else float('nan'):+.3f} prolif={r2_prolif if r2_prolif is not None else float('nan'):+.3f}")

tab = pd.DataFrame(rows)
tab.to_csv("results/05_model_vs_ceiling.csv", index=False)
print("\n=== TABLE 5 ===")
print(tab.to_string(index=False))
json.dump({"n_matched_lines": len(need), "n_genes_used": TOP,
           "table": tab.to_dict("records")},
          open("results/05_model_vs_ceiling.json", "w", encoding="utf-8"),
          indent=2, ensure_ascii=False)
