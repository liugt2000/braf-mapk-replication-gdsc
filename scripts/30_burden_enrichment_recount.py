"""
30_burden_enrichment_recount.py
-------------------------------
目的: 精确重算稿件 Results S2 引用的 "96%" 富集比例, 给出真实分子/分母,
并核对 19 号脚本全部 stdout 数字(环境重建后的复现审计).

背景: v3 稿件写 "96% of the 34 driver genes". 事实核查:
  - 驱动基因面板 = 33 个(脚本 19 DRIVERS 列表与 18_gene_frequencies.csv 均为 33 行);
  - "96%" 的分母 = 通过 n_alt>=8 且 n_wt>=20 门槛的**可检验**基因数
    (既不是 33 个面板基因, 也不是图 1C 硬编码显示的 20 个);
  - 本脚本按 19 号脚本的确定性逻辑逐行复刻重算, 输出精确分子/分母.

输出:
  results/30_burden_enrichment.csv              每基因 n_alt / n_wt / d_burden (33 基因全表)
  results/30_burden_adjusted_scan_recomputed.csv 负荷校正全扫描独立复算(不覆盖 19 号产物)
  results/30_recount_summary.json               全部校验数字
"""
import json
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
    out = {}
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

    # ---- 1. 线级别全局混杂(复刻 19 号脚本第 1 节) ----
    P = d.pivot_table(index="nm", columns="DRUG_NAME", values="LN_IC50",
                      aggfunc="mean")
    P = P[(P.notna().sum(axis=0) >= 40).index[P.notna().sum(axis=0) >= 40]]
    R = P.rank(axis=0)
    mean_rank = R.mean(axis=1)
    b = pd.Series({k: burden.get(k, np.nan) for k in mean_rank.index})
    ok = b.notna() & mean_rank.notna()
    r, p = stats.pearsonr(b[ok], mean_rank[ok])
    out["line_level"] = {"n": int(ok.sum()), "r": round(float(r), 4),
                         "p": float(p)}

    # ---- 2. 基因层面富集(复刻 19 号脚本第 2 节, 全部 33 基因) ----
    alt = {g: {norm(x) for x in md.loc[md.gene_symbol == g, "model_name"]}
           for g in DRIVERS}
    rows = []
    for g in DRIVERS:
        isalt = np.array([k in alt[g] for k in mean_rank.index])
        n_alt = int(isalt.sum())
        n_wt = int((~isalt).sum())
        if n_alt < 8 or n_wt < 20:
            rows.append(dict(gene=g, n_alt=n_alt, n_wt=n_wt,
                             d_burden=np.nan, testable=False))
            continue
        bb = b.values
        ok2 = pd.notna(bb) & mean_rank.notna().values
        am = isalt[ok2]
        d_b = (bb[ok2][am].mean() - bb[ok2][~am].mean()) / (
            np.sqrt(((am.sum() - 1) * bb[ok2][am].var(ddof=1)
                     + ((~am).sum() - 1) * bb[ok2][~am].var(ddof=1))
                    / (am.sum() + (~am).sum() - 2)) + 1e-9)
        rows.append(dict(gene=g, n_alt=n_alt, n_wt=n_wt,
                         d_burden=float(d_b), testable=True))
    B = pd.DataFrame(rows)
    Bt = B[B.testable].copy()
    n_pos = int((Bt.d_burden > 0).sum())
    n_test = int(Bt.shape[0])
    out["gene_level"] = {
        "n_genes_panel": len(DRIVERS),
        "n_testable": n_test,
        "n_positive": n_pos,
        "positive_pct_exact": round(n_pos / n_test * 100, 2),
        "printed_pct_19style": f"{n_pos / n_test * 100:.0f}%",
        "median_d": round(float(Bt.d_burden.median()), 4),
        "positive_genes": sorted(Bt[Bt.d_burden > 0].gene.tolist()),
        "nonpositive_genes": Bt[Bt.d_burden <= 0].gene.tolist(),
        "nonpositive_details": Bt[Bt.d_burden <= 0][
            ["gene", "n_alt", "d_burden"]].to_dict("records"),
    }
    B.sort_values("d_burden", ascending=False, na_position="last").to_csv(
        "results/30_burden_enrichment.csv", index=False)

    # ---- 3. 负荷校正全扫描独立复算(不覆盖 19 号产物) ----
    drugs = [c for c in P.columns if P[c].notna().sum() >= 40]
    scan_rows = []
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
            scan_rows.append(dict(drug=drug, gene=g, n_alt=int(fa.sum()),
                                  beta_adj=float(beta[1]),
                                  p_adj=float(2 * stats.t.sf(abs(t), len(s) - 3))))
    O = pd.DataFrame(scan_rows)
    O["q_adj"] = false_discovery_control(O.p_adj.values, method="bh")
    O.to_csv("results/30_burden_adjusted_scan_recomputed.csv", index=False)
    sig = O[(O.q_adj < 0.05) & (O.beta_adj.abs() >= 0.3)]
    top8 = sig.reindex(sig.beta_adj.abs().sort_values(ascending=False).index).head(8)
    out["burden_adjusted"] = {
        "n_combos": int(len(O)),
        "n_significant": int(len(sig)),
        "genes_involved": sorted(sig.gene.unique().tolist()),
        "top8_by_abs_beta": [
            {"drug": rr.drug, "gene": rr.gene,
             "beta_adj": round(float(rr.beta_adj), 4)}
            for rr in top8.itertuples()],
    }
    try:
        O19 = pd.read_csv("results/19_burden_adjusted_scan.csv")
        merged = O19.merge(O, on=["drug", "gene"], suffixes=("_19", "_30"))
        out["reproduce_19_scan"] = {
            "rows_19": int(len(O19)),
            "rows_30": int(len(O)),
            "max_beta_abs_diff": float(
                (merged.beta_adj_19 - merged.beta_adj_30).abs().max()),
            "max_p_abs_diff": float(
                (merged.p_adj_19 - merged.p_adj_30).abs().max()),
        }
    except FileNotFoundError:
        out["reproduce_19_scan"] = "19_burden_adjusted_scan.csv not found"

    # ---- 4. 稿件引用值逐项核对 ----
    checks = {
        "r_expect_-0.03": f"{r:+.2f}" == "-0.03",
        "p_expect_0.37": f"{p:.2f}" == "0.37",
        "n_lines_expect_962": int(ok.sum()) == 962,
        "median_d_expect_+1.40": f"{Bt.d_burden.median():+.2f}" == "+1.40",
        "n_combos_expect_2891": len(O) == 2891,
        "n_sig_expect_48": len(sig) == 48,
    }
    for g, expect in [("MLH1", 3.18), ("BRCA2", 3.14), ("MSH6", 2.62),
                      ("FBXW7", 2.30), ("BRAF", 0.20)]:
        got = Bt[Bt.gene == g].d_burden
        checks[f"{g}_expect_{expect:+.2f}"] = ((not got.empty)
                                               and f"{got.iloc[0]:+.2f}" == f"{expect:+.2f}")
    out["manuscript_value_checks"] = checks

    print(json.dumps(out, indent=2, ensure_ascii=False, default=str))
    with open("results/30_recount_summary.json", "w", encoding="utf-8") as f:
        json.dump(out, f, indent=2, ensure_ascii=False, default=str)
    print("\n写入 results/30_burden_enrichment.csv / 30_burden_adjusted_scan_recomputed.csv"
          " / 30_recount_summary.json")


if __name__ == "__main__":
    main()
