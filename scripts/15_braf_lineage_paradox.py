"""
15_braf_lineage_paradox.py
--------------------------
核心问题: BRAF 突变 -> BRAF/MEK 抑制剂敏感, 是 GDSC 全资源里唯一通过
"双活动同向 + 大效应量 + 双 p<0.001" 过滤的基因-药物轴。

但它是在泛癌上成立的。本脚本检验:
  1. 按瘤种拆开, BRAF 突变对 MAPK 抑制剂的效果在哪一类肿瘤里最强、哪一类里失效
  2. 结直肠癌是否是这个轴的例外(即"CRC BRAF 悖论"能否在公共数据中复现)
  3. 若失效, 能否用一个可测量的转录程序(EGFR/RTK 反馈)解释它

所有统计仅在 IC50 可识别(落在实测浓度区间内)的样本上进行。
"""
import re
import numpy as np
import pandas as pd
from scipy import stats

BASE = "data/raw"
MAPKi = ["Dabrafenib", "Trametinib", "Selumetinib", "Refametinib",
         "PLX-4720", "PD0325901", "CI-1040"]


def norm(s):
    return re.sub(r"[^A-Z0-9]", "", str(s).upper())


def braf_v600():
    m = pd.read_csv(f"{BASE}/mutations_latest.csv.gz")
    m = m[(m.gene_symbol == "BRAF") & (m.cancer_driver.astype(str).str.lower() == "true")]
    pm = m.protein_mutation.fillna("").astype(str)
    v600 = m[pm.str.contains("V600", na=False)]
    return ({norm(x) for x in v600.model_name.unique()},
            {norm(x) for x in m.model_name.unique()})


def prep(fn):
    d = pd.read_csv(f"{BASE}/{fn}")
    d = d[d.DRUG_NAME.isin(MAPKi)].copy()
    d["ic50"] = np.exp(d.LN_IC50)
    d["ident"] = (d.ic50 >= d.MIN_CONC) & (d.ic50 <= d.MAX_CONC)
    return d[d.ident]


def main():
    v600, anybraf = braf_v600()
    print(f"BRAF 任意 driver 突变 {len(anybraf)} 株; 其中 V600 {len(v600)} 株")

    rows = []
    for fn, tag in [("GDSC2_fitted_dose_response_24Jul22.csv", "GDSC2"),
                    ("GDSC1_fitted_dose_response_24Jul22.csv", "GDSC1")]:
        d = prep(fn)
        d["mut"] = d.CELL_LINE_NAME.map(lambda x: "V600" if norm(x) in v600
                                        else ("BRAFother" if norm(x) in anybraf else "WT"))
        for drug in MAPKi:
            s = d[d.DRUG_NAME == drug]
            if s.empty:
                continue
            for lin, g in s.groupby("TCGA_DESC"):
                a = g[g.mut == "V600"].LN_IC50
                w = g[g.mut == "WT"].LN_IC50
                if len(a) >= 3 and len(w) >= 5:
                    sp = np.sqrt(((len(a) - 1) * a.var(ddof=1) + (len(w) - 1) * w.var(ddof=1))
                                 / (len(a) + len(w) - 2))
                    rows.append(dict(campaign=tag, drug=drug, lineage=lin,
                                     n_v600=len(a), n_wt=len(w),
                                     med_v600=a.median(), med_wt=w.median(),
                                     d=(a.mean() - w.mean()) / sp if sp > 0 else np.nan,
                                     p=stats.mannwhitneyu(a, w).pvalue))
    r = pd.DataFrame(rows)
    r.to_csv("results/15_braf_by_lineage.csv", index=False)

    print("\n" + "=" * 78)
    print("BRAF V600 vs WT 的 ln(IC50) 差异, 按瘤种拆分 (仅可识别样本)")
    print("=" * 78)
    for drug in ["Dabrafenib", "Selumetinib", "Trametinib"]:
        t = r[(r.drug == drug) & (r.campaign == "GDSC2")].sort_values("d")
        if t.empty:
            continue
        print(f"\n--- {drug} (GDSC2) ---")
        print(t[["lineage", "n_v600", "n_wt", "med_v600", "med_wt", "d", "p"]]
              .to_string(index=False, float_format=lambda v: f"{v:.2f}"))
        c = r[(r.drug == drug) & (r.campaign == "GDSC1") & (r.lineage == "COREAD")]
        if not c.empty:
            x = c.iloc[0]
            print(f"  GDSC1 复制 (COREAD): n_v600={x.n_v600:.0f} n_wt={x.n_wt:.0f} "
                  f"d={x.d:.2f} p={x.p:.3f}")

    print("\n" + "=" * 78)
    print("结直肠 vs 其他瘤种: BRAF 突变带来的敏感性获益(以倍数表示)")
    print("=" * 78)
    for drug in ["Dabrafenib", "Selumetinib"]:
        t = r[(r.drug == drug) & (r.campaign == "GDSC2")].copy()
        if t.empty:
            continue
        t["fold"] = np.exp(t.med_wt - t.med_v600)
        t = t.sort_values("fold", ascending=False)
        print(f"\n{drug}:  倍数 = exp(WT中位 - V600中位), 越大说明 BRAF 突变带来的敏感性越强")
        print(t[["lineage", "n_v600", "fold", "p"]].to_string(
            index=False, float_format=lambda v: f"{v:.2f}"))

    # ---- 剂量反应曲线形态: V600 CRC vs V600 其他 ----
    print("\n" + "=" * 78)
    print("BRAF V600 内部: CRC 与其他瘤种的达拉非尼敏感性直接对比")
    print("=" * 78)
    d = prep("GDSC2_fitted_dose_response_24Jul22.csv")
    d = d[d.DRUG_NAME == "Dabrafenib"].copy()
    d["mut"] = d.CELL_LINE_NAME.map(lambda x: norm(x) in v600)
    v = d[d.mut]
    for lin, g in v.groupby("TCGA_DESC"):
        if len(g) >= 3:
            print(f"  {lin:14s} n={len(g):3d}  中位 ln(IC50)={g.LN_IC50.median():6.2f} "
                  f" 中位 IC50={np.exp(g.LN_IC50.median()):8.1f} nM")
    crc = v[v.TCGA_DESC == "COREAD"].LN_IC50
    oth = v[v.TCGA_DESC != "COREAD"].LN_IC50
    if len(crc) >= 3 and len(oth) >= 5:
        u = stats.mannwhitneyu(crc, oth)
        print(f"\n  COREAD(n={len(crc)}) vs 其他瘤种(n={len(oth)}): "
              f"p={u.pvalue:.4f}  倍数={np.exp(oth.median()-crc.median()):.2f}x")
        print(f"  若 >1, 表示 BRAF 突变的 CRC 比 BRAF 突变的他瘤更耐药")


if __name__ == "__main__":
    main()
