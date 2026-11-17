"""
19_mutational_burden_confounder.py
----------------------------------
发现: 单活动 FDR(q<0.05, |d|>=0.5) 给出 61 个"突变->耐药"关联, 且**方向 100% 为正**。
100% 同向在纯偶然下概率约 2^-61, 因此必然存在系统性混杂。

假设: 混杂是"细胞系突变负荷"。高突变负荷的细胞系(MSI/POLE/MMR 缺陷)若同时
      在任意驱动基因上更常被标注为突变, 就会在每一个基因上制造出"突变->耐药"的假象。

检验分四步:
 1. 有无线级别的全局混杂: 突变负荷 vs 跨药平均 IC50
 2. 61 个 FDR 命中的效应量能否被负荷解释(偏相关 / 分层)
 3. BRAF 是否豁免于该混杂
 4. 若控制负荷, 存活者是否仍只有 BRAF
"""
import re
import numpy as np
import pandas as pd
from scipy import stats
from scipy.stats import false_discovery_control

BASE = "data/raw"
DRIVERS = ["KRAS", "NRAS", "BRAF", "MAP2K1", "NF1", "PIK3CA", "PTEN", "AKT1",
           "TP53", "CDKN2A", "RB1", "ATM", "APC", "CTNNB1", "RNF43", "SMAD2",
           "SMAD4", "TGFBR2", "MLH1", "MSH2", "MSH6", "PMS2", "BRCA1", "BRCA2",
           "ERCC1", "ERCC2", "ERCC5", "ERBB2", "EGFR", "MET", "ARID1A",
           "FBXW7", "POLE"]


def norm(s):
    return re.sub(r"[^A-Z0-9]", "", str(s).upper())


def main():
    m = pd.read_csv(f"{BASE}/mutations_latest.csv.gz")
    pm = m.protein_mutation.fillna("").astype(str)
    m["damaging"] = (m.cancer_driver.astype(str).str.lower() == "true") | \
                    pm.str.contains(r"\*|fs", regex=True)
    md = m[m.damaging]
    burden = (md.drop_duplicates(["model_name", "gene_symbol"])
              .groupby("model_name").size().rename("n_dmg_genes"))
    burden.index = [norm(x) for x in burden.index]
    burden = burden[~burden.index.duplicated()]

    d = pd.read_csv(f"{BASE}/GDSC2_fitted_dose_response_24Jul22.csv")
    d["ic50"] = np.exp(d.LN_IC50)
    d["ident"] = (d.ic50 >= d.MIN_CONC) & (d.ic50 <= d.MAX_CONC)
    d = d[d.ident].copy()
    d["nm"] = d.CELL_LINE_NAME.map(norm)
    d["burden"] = d.nm.map(burden)

    # ---- 1. 线级别全局混杂 ----
    print("=" * 74)
    print("1. 细胞系层面的全局混杂: 突变负荷 vs 跨药平均敏感性")
    print("=" * 74)
    P = d.pivot_table(index="nm", columns="DRUG_NAME", values="LN_IC50",
                      aggfunc="mean")
    P = P[(P.notna().sum(axis=0) >= 40).index[P.notna().sum(axis=0) >= 40]]
    R = P.rank(axis=0)                       # 每药内排序, 消除药物尺度差异
    mean_rank = R.mean(axis=1)
    b = pd.Series({k: burden.get(k, np.nan) for k in mean_rank.index})
    ok = b.notna() & mean_rank.notna()
    r, p = stats.pearsonr(b[ok], mean_rank[ok])
    print(f"  细胞系 n={int(ok.sum())}")
    print(f"  突变负荷 vs 跨药平均 IC50 排序: r={r:+.3f} (p={p:.2e})")
    q = pd.qcut(b[ok], 4, labels=["Q1 最低", "Q2", "Q3", "Q4 最高"])
    print("  按负荷四分位看平均 IC50 排序(排序值越大 = 越耐药):")
    print(mean_rank[ok].groupby(q, observed=True).agg(["mean", "size"])
          .to_string(float_format=lambda v: f"{v:.3f}"))

    # ---- 2 & 3. 基因层面: 负荷是否解释效应 ----
    print("\n" + "=" * 74)
    print("2. 基因层面: 每个基因的突变组是否富含高负荷细胞系")
    print("=" * 74)
    alt = {g: {norm(x) for x in md.loc[md.gene_symbol == g, "model_name"]}
           for g in DRIVERS}
    rows = []
    for g in DRIVERS:
        isalt = np.array([k in alt[g] for k in mean_rank.index])
        if isalt.sum() < 8 or (~isalt).sum() < 20:
            continue
        bb = b.values
        ok2 = pd.notna(bb) & mean_rank.notna().values
        am = isalt[ok2]
        d_b = (bb[ok2][am].mean() - bb[ok2][~am].mean()) / (
            np.sqrt(((am.sum() - 1) * bb[ok2][am].var(ddof=1)
                     + ((~am).sum() - 1) * bb[ok2][~am].var(ddof=1))
                    / (am.sum() + (~am).sum() - 2)) + 1e-9)
        rows.append(dict(gene=g, n_alt=int(isalt.sum()), d_burden=d_b))
    B = pd.DataFrame(rows).sort_values("d_burden", ascending=False)
    print(B.to_string(index=False, float_format=lambda v: f"{v:+.2f}"))
    print(f"\n  突变组更富高负荷的基因占 {(B.d_burden>0).mean()*100:.0f}%  "
          f"(中位 d={B.d_burden.median():+.2f})")

    # ---- 4. 负荷作为协变量后, BRAF 效应是否仍在 ----
    print("\n" + "=" * 74)
    print("3. 负荷校正后重跑关键药物: BRAF 效应是否仍然成立")
    print("=" * 74)
    for drug in ["Dabrafenib", "Selumetinib", "Trametinib"]:
        s = d[d.DRUG_NAME == drug].dropna(subset=["burden"])
        s = s.drop_duplicates("nm")
        if len(s) < 40:
            continue
        s["is_braf"] = s.nm.isin(alt["BRAF"]).astype(int)
        x = s[["LN_IC50", "is_braf", "burden"]].dropna()
        # 未校正
        a = x[x.is_braf == 1].LN_IC50
        w = x[x.is_braf == 0].LN_IC50
        raw_d = (a.mean() - w.mean()) / np.sqrt(
            ((len(a) - 1) * a.var(ddof=1) + (len(w) - 1) * w.var(ddof=1))
            / (len(a) + len(w) - 2))
        # 校正负荷(OLS)
        X = np.column_stack([np.ones(len(x)), x.is_braf, x.burden])
        beta, *_ = np.linalg.lstsq(X, x.LN_IC50.values, rcond=None)
        resid = x.LN_IC50.values - X @ beta
        s2 = resid @ resid / (len(x) - 3)
        se = np.sqrt(np.diag(s2 * np.linalg.pinv(X.T @ X)))
        t = beta / se
        pv = 2 * stats.t.sf(np.abs(t), len(x) - 3)
        print(f"  {drug:12s} n={len(x):4d}  未校正 d={raw_d:+.2f}  "
              f"负荷校正后 BRAF beta={beta[1]:+.3f} (p={pv[1]:.2e})  "
              f"负荷 beta={beta[2]:+.3f} (p={pv[2]:.2e})")

    # ---- 5. 全表: 负荷校正后的基因-药物扫描 ----
    print("\n" + "=" * 74)
    print("4. 对所有 (基因,药物) 做负荷校正, 再看还剩哪些显著")
    print("=" * 74)
    drugs = [c for c in P.columns if P[c].notna().sum() >= 40]
    out = []
    for drug in drugs:
        s = d[d.DRUG_NAME == drug].dropna(subset=["burden"]).drop_duplicates("nm")
        if len(s) < 40:
            continue
        for g in DRIVERS:
            fa = s.nm.isin(alt[g]).astype(int).values
            if fa.sum() < 8 or (1 - fa).sum() < 20:
                continue
            X = np.column_stack([np.ones(len(s)), fa, s.burden.values])
            y = s.LN_IC50.values
            beta, *_ = np.linalg.lstsq(X, y, rcond=None)
            resid = y - X @ beta
            s2 = resid @ resid / (len(s) - 3)
            se = np.sqrt(np.diag(s2 * np.linalg.pinv(X.T @ X)))
            t = beta[1] / se[1]
            out.append(dict(drug=drug, gene=g, n_alt=int(fa.sum()),
                            beta_adj=beta[1],
                            p_adj=2 * stats.t.sf(abs(t), len(s) - 3)))
    O = pd.DataFrame(out)
    O["q_adj"] = false_discovery_control(O.p_adj.values, method="bh")
    O.to_csv("results/19_burden_adjusted_scan.csv", index=False)
    sig = O[(O.q_adj < 0.05) & (O.beta_adj.abs() >= 0.3)]
    print(f"  检验 {len(O)} 个组合; 负荷校正后 q<0.05 且 |beta|>=0.3: {len(sig)} 个")
    print(sig.sort_values("beta_adj").head(15).to_string(
        index=False, float_format=lambda v: f"{v:.3g}"))
    print(f"\n  这些显著组合涉及的基因: {sorted(sig.gene.unique())}")


if __name__ == "__main__":
    main()
