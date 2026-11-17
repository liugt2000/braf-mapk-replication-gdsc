"""
16_rtk_score_vs_brafi.py
------------------------
论文的核心断言待检:
  在 BRAF V600 模型内部, 对 BRAF 抑制剂的敏感性由基线 RTK/MAPK 反馈张力决定,
  而不是由 BRAF 本身决定 —— 结直肠恰好是这个张力的高位端。

检验:
  1. RTK/feedback 转录评分 与 达拉非尼 ln(IC50) 的关系(全体 BRAF V600)
  2. 校正瘤种后是否仍成立(否则只是"学到的瘤种标签")
  3. 增殖评分作为竞争解释(排除"只是增殖快慢")
  4. 结直肠在这个评分上的位置
"""
import gzip
import re
import numpy as np
import pandas as pd
from scipy import stats

BASE = "data/raw"
RTK = ["EGFR", "ERBB2", "ERBB3", "AREG", "EREG", "TGFA", "EPGN",
       "DUSP4", "DUSP6", "SPRY2", "SPRY4", "ETV4", "ETV5", "PHLDA1",
       "EPHA2", "FGFR1", "MET", "IGF1R", "NRG1", "HGF"]
PROLIF = ["MKI67", "PCNA", "CCNB1", "CDK1", "TOP2A", "BUB1", "AURKA",
          "PLK1", "CDC20", "CCNA2", "RRM2", "MCM2", "MCM6"]


def norm(s):
    return re.sub(r"[^A-Z0-9]", "", str(s).upper())


def braf_v600():
    m = pd.read_csv(f"{BASE}/mutations_latest.csv.gz")
    m = m[(m.gene_symbol == "BRAF") &
          (m.cancer_driver.astype(str).str.lower() == "true")]
    pm = m.protein_mutation.fillna("").astype(str)
    return {norm(x) for x in m[pm.str.contains("V600")].model_name.unique()}


def main():
    v600 = braf_v600()
    d = pd.read_csv(f"{BASE}/GDSC2_fitted_dose_response_24Jul22.csv")
    d = d[d.DRUG_NAME == "Dabrafenib"].copy()
    d["ic50"] = np.exp(d.LN_IC50)
    d["ident"] = (d.ic50 >= d.MIN_CONC) & (d.ic50 <= d.MAX_CONC)
    d = d[d.ident]
    d["v600"] = d.CELL_LINE_NAME.map(lambda x: norm(x) in v600)
    print(f"达拉非尼可识别样本 {len(d)} 中, BRAF V600 {d.v600.sum()} 株")

    # ---- 表达谱 ----
    with gzip.open(f"{BASE}/CCLE_RNAseq_genes_rpkm_20180929.gct.gz", "rt") as f:
        f.readline(); f.readline()
        hdr = f.readline().rstrip("\n").split("\t")
    samples = hdr[2:]
    lk = {}
    for s in samples:
        lk.setdefault(norm(s.split("_", 1)[0]), []).append(s)
    lk = {k: v[0] for k, v in lk.items()}
    d["ccle"] = d.CELL_LINE_NAME.map(lambda x: lk.get(norm(x)))
    need = sorted(d.ccle.dropna().unique())
    print(f"  其中可匹配表达谱 {len(need)} 株")

    raw = pd.read_csv(f"{BASE}/CCLE_RNAseq_genes_rpkm_20180929.gct.gz", sep="\t",
                      skiprows=2, usecols=["Name", "Description"] + need,
                      dtype={x: "float32" for x in need}, compression="gzip")
    raw = raw.dropna(subset=["Description"])
    expr = np.log2(raw.set_index("Description")[need]
                   .groupby(level=0).max() + 1)
    print(f"  表达矩阵 {expr.shape}")

    def score(genes):
        g = [x for x in genes if x in expr.index]
        z = expr.loc[g].sub(expr.loc[g].mean(axis=1), axis=0).div(
            expr.loc[g].std(axis=1).replace(0, np.nan), axis=0)
        return z.mean(axis=0), g

    rtk, rtk_g = score(RTK)
    prol, prol_g = score(PROLIF)
    print(f"  RTK 评分基因 {len(rtk_g)}/{len(RTK)}: 缺失 {set(RTK)-set(rtk_g)}")
    print(f"  增殖评分基因 {len(prol_g)}/{len(PROLIF)}")

    w = d.dropna(subset=["ccle"]).drop_duplicates("ccle").copy()
    w["rtk"] = w.ccle.map(rtk)
    w["prol"] = w.ccle.map(prol)
    w = w.dropna(subset=["rtk", "prol"])
    w.to_csv("results/16_brafi_rtk_score.csv", index=False)

    v = w[w.v600]
    print("\n" + "=" * 72)
    print("1. BRAF V600 内部: RTK 评分 vs 达拉非尼敏感性")
    print("=" * 72)
    print(f"  n = {len(v)}")
    r, p = stats.pearsonr(v.rtk, v.LN_IC50)
    rs, ps = stats.spearmanr(v.rtk, v.LN_IC50)
    print(f"  Pearson r = {r:+.3f} (p={p:.4f})   Spearman rho = {rs:+.3f} (p={ps:.4f})")
    r2_, p2_ = stats.pearsonr(v.prol, v.LN_IC50)
    print(f"  增殖评分   r = {r2_:+.3f} (p={p2_:.4f})   [竞争解释]")

    # 校正瘤种
    if len(v) >= 10:
        vv = v.copy()
        for g in ["rtk", "prol"]:
            vv[g] = (vv[g] - vv[g].mean()) / vv[g].std()
        y = (vv.LN_IC50 - vv.LN_IC50.mean()) / vv.LN_IC50.std()
        X = np.column_stack([np.ones(len(vv)), vv.rtk, vv.prol])
        beta, *_ = np.linalg.lstsq(X, y, rcond=None)
        resid = y - X @ beta
        dof = len(y) - 3
        s2 = resid @ resid / dof
        cov = s2 * np.linalg.pinv(X.T @ X)
        se = np.sqrt(np.diag(cov))
        tval = beta / se
        pval = 2 * stats.t.sf(np.abs(tval), dof)
        print(f"\n  多变量(标准化): RTK beta={beta[1]:+.3f} (p={pval[1]:.4f})"
              f" | 增殖 beta={beta[2]:+.3f} (p={pval[2]:.4f})")

    print("\n" + "=" * 72)
    print("2. 各瘤种的 RTK 评分位置与达拉非尼敏感性")
    print("=" * 72)
    tab = v.groupby("TCGA_DESC").agg(
        n=("rtk", "size"), rtk_mean=("rtk", "mean"),
        ln_ic50_med=("LN_IC50", "median")).sort_values("rtk_mean", ascending=False)
    print(tab.to_string(float_format=lambda x: f"{x:.2f}"))

    print("\n" + "=" * 72)
    print("3. 全 panel 视角: 结直肠在 RTK 张力谱上的位置(BRAF 突变与否分层)")
    print("=" * 72)
    allw = w.copy()
    print(allw.groupby(["v600"]).agg(n=("rtk", "size"),
                                     rtk_mean=("rtk", "mean")).to_string(
        float_format=lambda x: f"{x:.2f}"))
    crc = w[w.TCGA_DESC == "COREAD"]
    oth = w[w.TCGA_DESC != "COREAD"]
    u = stats.mannwhitneyu(crc.rtk, oth.rtk)
    print(f"\n  COREAD(n={len(crc)}) RTK 均值={crc.rtk.mean():.3f} vs "
          f"其他(n={len(oth)}) {oth.rtk.mean():.3f}  p={u.pvalue:.4f}")


if __name__ == "__main__":
    main()
