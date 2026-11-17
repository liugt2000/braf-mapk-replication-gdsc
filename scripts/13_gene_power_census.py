"""
13_gene_power_census.py
-----------------------
在"某个基因/通路导致耐药"这一主题下, 先用已测得的删失率与跨活动可靠性做一次效能普查,
而不是直接跳进建模。

对每个 (候选驱动基因, 药物) 组合输出:
  n_altered / n_wt   —— 仅统计 IC50 落在实测浓度区间内(可识别)的细胞系
  cohen_d            —— 效应量(改变组 - 野生组, 以 ln(IC50) 为单位)
  p_mwu              —— Mann-Whitney
  d_gdsc1, p_gdsc1   —— 在 GDSC1 独立筛选中同向复制的程度(最关键的一道过滤)

结论口径: 只有"可识别样本量足够 + 效应量大 + 在 GDSC1 同向复制"的基因-药物对,
才值得作为论文主角。其余一律判为不可做。
"""
import numpy as np
import pandas as pd
from scipy import stats

BASE = "data/raw"
DRIVERS = [
    # MAPK
    "KRAS", "NRAS", "BRAF", "MAP2K1", "NF1", "RASA1",
    # PI3K
    "PIK3CA", "PIK3R1", "PTEN", "AKT1", "MTOR", "TSC1", "TSC2",
    # p53 / cell cycle
    "TP53", "CDKN2A", "RB1", "ATM", "CHEK2",
    # Wnt / TGF-beta (CRC 主轴)
    "APC", "CTNNB1", "RNF43", "AXIN1", "AXIN2", "SMAD2", "SMAD3", "SMAD4",
    "TGFBR1", "TGFBR2", "BMPR1A",
    # DDR / MMR
    "MLH1", "MSH2", "MSH6", "PMS2", "POLE", "POLD1", "BRCA1", "BRCA2",
    "ERCC1", "ERCC2", "ERCC4", "ERCC5", "XRCC1",
    # RTK / 其他
    "ERBB2", "ERBB3", "EGFR", "MET", "FGFR2", "IGF1R", "ALK", "RET",
    "ARID1A", "ARID2", "FBXW7", "SOX9", "TCF7L2", "GPC5", "MYC",
]
# 机制上已知的耐药轴, 用于给结果贴标签
AXIS = {
    "KRAS": "MAPK", "NRAS": "MAPK", "BRAF": "MAPK", "MAP2K1": "MAPK",
    "NF1": "MAPK", "RASA1": "MAPK",
    "PIK3CA": "PI3K", "PIK3R1": "PI3K", "PTEN": "PI3K", "AKT1": "PI3K",
    "MTOR": "PI3K", "TSC1": "PI3K", "TSC2": "PI3K",
    "TP53": "p53", "CDKN2A": "cellcycle", "RB1": "cellcycle",
    "ATM": "DDR", "CHEK2": "DDR",
    "APC": "Wnt", "CTNNB1": "Wnt", "RNF43": "Wnt", "AXIN1": "Wnt", "AXIN2": "Wnt",
    "SMAD2": "TGFb", "SMAD3": "TGFb", "SMAD4": "TGFb", "TGFBR1": "TGFb",
    "TGFBR2": "TGFb", "BMPR1A": "TGFb",
    "MLH1": "MMR", "MSH2": "MMR", "MSH6": "MMR", "PMS2": "MMR",
    "POLE": "DDR", "POLD1": "DDR", "BRCA1": "DDR", "BRCA2": "DDR",
    "ERCC1": "NER", "ERCC2": "NER", "ERCC4": "NER", "ERCC5": "NER", "XRCC1": "BER",
    "ERBB2": "RTK", "ERBB3": "RTK", "EGFR": "RTK", "MET": "RTK",
    "FGFR2": "RTK", "IGF1R": "RTK", "ALK": "RTK", "RET": "RTK",
}


def norm(s):
    import re
    return re.sub(r"[^A-Z0-9]", "", str(s).upper())


def load_mut():
    m = pd.read_csv(f"{BASE}/mutations_latest.csv.gz")
    m = m[m.gene_symbol.isin(DRIVERS)].copy()
    # 只保留能破坏功能的变异: 蛋白改变且被标记为 driver; 或明确的截断/移码
    pm = m.protein_mutation.fillna("-").astype(str)
    trunc = pm.str.contains(r"\*|fs|del|Ter", regex=True)
    m["damaging"] = m.cancer_driver.astype(str).str.lower().eq("true") | trunc
    m = m[m.damaging]
    return m[["model_name", "gene_symbol"]].drop_duplicates()


def load_drug(fn, tag):
    d = pd.read_csv(f"{BASE}/{fn}")
    d = d[d.TCGA_DESC == "COREAD"].copy()
    d["ic50"] = np.exp(d.LN_IC50)
    d["ident"] = (d.ic50 >= d.MIN_CONC) & (d.ic50 <= d.MAX_CONC)
    d["src"] = tag
    return d


def main():
    print("载入突变 ...")
    mut = load_mut()
    print(f"  damaging 变异记录 {len(mut):,}  涉及基因 {mut.gene_symbol.nunique()}")

    g2 = load_drug("GDSC2_fitted_dose_response_24Jul22.csv", "GDSC2")
    g1 = load_drug("GDSC1_fitted_dose_response_24Jul22.csv", "GDSC1")

    def pivot(d, ident_only):
        x = d[d.ident] if ident_only else d
        return x.pivot_table(index="CELL_LINE_NAME", columns="DRUG_NAME",
                             values="LN_IC50", aggfunc="mean")

    P2 = pivot(g2, True)
    P1 = pivot(g1, True)
    print(f"  GDSC2 COREAD 可识别矩阵: {P2.shape[0]} 系 x {P2.shape[1]} 药")
    print(f"  GDSC1 COREAD 可识别矩阵: {P1.shape[0]} 系 x {P1.shape[1]} 药")

    # 建 (系, 基因) -> 1
    alt = {}
    for g in DRIVERS:
        alt[g] = set(mut.loc[mut.gene_symbol == g, "model_name"].map(norm))

    rows = []
    for drug in P2.columns:
        s2 = P2[drug].dropna()
        if len(s2) < 12:
            continue
        for g in DRIVERS:
            if g not in AXIS:
                continue
            a = set(s2.index) & alt[g]
            w = set(s2.index) - alt[g]
            if len(a) < 4 or len(w) < 6:
                continue
            va, vw = s2[list(a)].values, s2[list(w)].values
            d_eff = (va.mean() - vw.mean()) / np.sqrt(
                ((len(va) - 1) * va.var(ddof=1) + (len(vw) - 1) * vw.var(ddof=1))
                / max(len(va) + len(vw) - 2, 1))
            p = stats.mannwhitneyu(va, vw, alternative="two-sided").pvalue

            # GDSC1 独立复制
            d1 = np.nan
            p1 = np.nan
            n1a = n1w = 0
            if drug in P1.columns:
                t1 = P1[drug].dropna()
                a1 = [x for x in t1.index if x in alt[g]]
                w1 = [x for x in t1.index if x not in alt[g]]
                n1a, n1w = len(a1), len(w1)
                if n1a >= 3 and n1w >= 4:
                    va1, vw1 = t1[a1].values, t1[w1].values
                    d1 = (va1.mean() - vw1.mean()) / np.sqrt(
                        ((len(va1) - 1) * va1.var(ddof=1) + (len(vw1) - 1) * vw1.var(ddof=1))
                        / max(len(va1) + len(vw1) - 2, 1))
                    p1 = stats.mannwhitneyu(va1, vw1, alternative="two-sided").pvalue

            rows.append(dict(
                drug=drug, gene=g, axis=AXIS[g],
                n_alt=len(a), n_wt=len(w), d_gdsc2=d_eff, p_gdsc2=p,
                n_alt1=n1a, n_wt1=n1w, d_gdsc1=d1, p_gdsc1=p1,
                replicated=(not np.isnan(d1)) and np.sign(d1) == np.sign(d_eff),
            ))

    r = pd.DataFrame(rows)
    r.to_csv("results/13_gene_power_census.csv", index=False)
    print(f"\n共检验 {len(r)} 个 (基因, 药物) 组合 -> results/13_gene_power_census.csv")

    # 头号候选: 效应量大 + 两活动同向 + GDSC1 也显著
    key = r[(r.replicated) & (r.d_gdsc2.abs() >= 0.8)].copy()
    key["score"] = key.d_gdsc2.abs() * np.minimum(key.n_alt, 12) / 12
    key = key.sort_values("score", ascending=False)
    print(f"\n=== 满足'效应量>=0.8 且 GDSC1 同向复制'的组合: {len(key)} 个 ===")
    cols = ["drug", "gene", "axis", "n_alt", "n_wt", "d_gdsc2", "p_gdsc2",
            "n_alt1", "d_gdsc1", "p_gdsc1"]
    print(key[cols].head(25).to_string(index=False,
                                       float_format=lambda v: f"{v:.3f}"))

    # 按 (基因, 轴) 汇总: 同一个基因在多少个药上同向复制
    print("\n=== 通路级证据: 每个驱动基因在多少药物上产生同向效应 ===")
    agg = r.groupby(["gene", "axis"]).agg(
        n_pairs=("drug", "size"),
        n_rep=("replicated", "sum"),
        max_d=("d_gdsc2", lambda s: s.abs().max()),
        median_d=("d_gdsc2", lambda s: s.abs().median()),
    ).reset_index()
    agg["rep_frac"] = agg.n_rep / agg.n_pairs
    agg = agg[agg.n_pairs >= 8].sort_values(["rep_frac", "max_d"], ascending=False)
    print(agg.head(20).to_string(index=False, float_format=lambda v: f"{v:.3f}"))


if __name__ == "__main__":
    main()
