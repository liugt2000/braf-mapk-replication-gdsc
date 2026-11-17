"""
14_census_verify_and_pancancer.py
---------------------------------
两件事:
A. 验证 13 号脚本的细胞系名匹配是否成立(排除"普查无效"这一失败模式)
B. 把同一套效能普查放到泛癌 panel 上跑 —— 判断"基因->药物"这类主题到底是
   "CRC 这 46 株不够", 还是"公共数据整体做不到"
"""
import re
import numpy as np
import pandas as pd
from scipy import stats

BASE = "data/raw"
DRIVERS = ["KRAS", "NRAS", "BRAF", "MAP2K1", "NF1", "PIK3CA", "PTEN", "AKT1",
           "TP53", "CDKN2A", "RB1", "ATM", "APC", "CTNNB1", "RNF43", "SMAD2",
           "SMAD4", "TGFBR2", "MLH1", "MSH2", "MSH6", "PMS2", "BRCA1", "BRCA2",
           "ERCC1", "ERCC2", "ERCC5", "ERBB2", "EGFR", "MET", "ARID1A",
           "FBXW7", "POLE"]


def norm(s):
    return re.sub(r"[^A-Z0-9]", "", str(s).upper())


def load_mut():
    m = pd.read_csv(f"{BASE}/mutations_latest.csv.gz")
    pm = m.protein_mutation.fillna("-").astype(str)
    trunc = pm.str.contains(r"\*|fs", regex=True)
    m["damaging"] = m.cancer_driver.astype(str).str.lower().eq("true") | trunc
    m = m[m.damaging & m.gene_symbol.isin(DRIVERS)]
    return m[["model_name", "gene_symbol"]].drop_duplicates()


def main():
    mut = load_mut()
    g2 = pd.read_csv(f"{BASE}/GDSC2_fitted_dose_response_24Jul22.csv")
    g1 = pd.read_csv(f"{BASE}/GDSC1_fitted_dose_response_24Jul22.csv")
    for d in (g2, g1):
        d["ic50"] = np.exp(d.LN_IC50)
        d["ident"] = (d.ic50 >= d.MIN_CONC) & (d.ic50 <= d.MAX_CONC)

    # ---------- A. 匹配验证 ----------
    print("=" * 70)
    print("A. 细胞系名匹配验证")
    print("=" * 70)
    lines = set(g2.CELL_LINE_NAME.unique())
    mut_models = set(mut.model_name.unique())
    print(f"  GDSC2 全部细胞系 {len(lines)}   突变文件 model_name {mut_models.__len__()}")
    inter = {x for x in lines if norm(x) in {norm(y) for y in mut_models}}
    print(f"  GDSC2 与突变文件可匹配细胞系: {len(inter)} ({100*len(inter)/len(lines):.1f}%)")

    crc = set(g2.loc[g2.TCGA_DESC == "COREAD", "CELL_LINE_NAME"].unique())
    crc_hit = {x for x in crc if norm(x) in {norm(y) for y in mut_models}}
    print(f"  COREAD {len(crc)} 株, 匹配上突变 {len(crc_hit)} 株")
    print(f"  未匹配: {sorted(crc - crc_hit)}")

    print("\n  关键驱动基因在 COREAD 中的突变株数(仅可匹配者):")
    for g in ["KRAS", "BRAF", "TP53", "APC", "SMAD4", "PIK3CA", "ERBB2", "CTNNB1"]:
        s = {norm(x) for x in mut.loc[mut.gene_symbol == g, "model_name"]}
        n = len({x for x in crc_hit if norm(x) in s})
        print(f"    {g:8s} {n:3d} / {len(crc_hit)}  ({100*n/len(crc_hit):.0f}%)")

    # ---------- B. 泛癌普查 ----------
    print("\n" + "=" * 70)
    print("B. 泛癌 panel 上的同一套普查")
    print("=" * 70)

    def run(d, label):
        d = d[d.ident].copy()
        P = d.pivot_table(index="CELL_LINE_NAME", columns="DRUG_NAME",
                          values="LN_IC50", aggfunc="mean")
        print(f"  {label}: {P.shape[0]} 系 x {P.shape[1]} 药 (可识别)")
        alt = {g: {norm(x) for x in mut.loc[mut.gene_symbol == g, "model_name"]}
               for g in DRIVERS}
        rows = []
        for drug in P.columns:
            s = P[drug].dropna()
            if len(s) < 40:
                continue
            for g in DRIVERS:
                idx = [norm(x) for x in s.index]
                isalt = np.array([i in alt[g] for i in idx])
                if isalt.sum() < 8 or (~isalt).sum() < 20:
                    continue
                va = s.values[isalt]
                vw = s.values[~isalt]
                sp = np.sqrt(((isalt.sum() - 1) * va.var(ddof=1)
                              + ((~isalt).sum() - 1) * vw.var(ddof=1))
                             / (len(va) + len(vw) - 2))
                if sp == 0:
                    continue
                rows.append(dict(drug=drug, gene=g, n_alt=int(isalt.sum()),
                                 n_wt=int((~isalt).sum()),
                                 d=(va.mean() - vw.mean()) / sp,
                                 p=stats.mannwhitneyu(va, vw).pvalue))
        r = pd.DataFrame(rows)
        return r

    R2 = run(g2, "GDSC2 泛癌")
    R1 = run(g1, "GDSC1 泛癌")
    print(f"\n  可检验组合: GDSC2 {len(R2)}  GDSC1 {len(R1)}")

    mg = R2.merge(R1, on=["drug", "gene"], suffixes=("_2", "_1"))
    print(f"  两活动共同可检验: {len(mg)}")
    if len(mg):
        same = np.sign(mg.d_2) == np.sign(mg.d_1)
        print(f"  符号一致率: {100*same.mean():.1f}%   (随机期望 50%)")
        ok = mg[(mg.d_2.abs() >= 0.8) & (mg.d_1.abs() >= 0.5) & same]
        print(f"  效应量>=0.8(GDSC2) 且 >=0.5(GDSC1) 且同向: {len(ok)} 个")
        ok = ok.reindex(ok.d_2.abs().sort_values(ascending=False).index)
        print()
        print(ok[["drug", "gene", "n_alt_2", "n_wt_2", "d_2", "p_2",
                  "d_1", "p_1"]].head(30).to_string(
            index=False, float_format=lambda v: f"{v:.3f}"))
        mg.to_csv("results/14_pancancer_census.csv", index=False)

    # 对照: COREAD 内部的符号一致率
    print("\n" + "=" * 70)
    print("C. 同一口径下 COREAD 单独的符号一致率")
    print("=" * 70)
    prev = pd.read_csv("results/13_gene_power_census.csv")
    prev = prev.dropna(subset=["d_gdsc1"])
    if len(prev):
        s = np.sign(prev.d_gdsc2) == np.sign(prev.d_gdsc1)
        print(f"  可配对 {len(prev)} 个 (基因,药物) 组合, 符号一致率 {100*s.mean():.1f}%")
        print(f"  其中 GDSC1 侧 p<0.05 的比例: {100*(prev.p_gdsc1<0.05).mean():.1f}%")
        print(f"  GDSC2 侧 p<0.05 的比例:      {100*(prev.p_gdsc2<0.05).mean():.1f}%")


if __name__ == "__main__":
    main()
