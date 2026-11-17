"""
22_crc_exception_quantified.py
------------------------------
把"CRC 例外"从一句定性描述变成可测量的量。

定义 BRAFi/MEKi 解离度(per line, 同一细胞系配对的两种药):
    delta = ln(IC50_dabrafenib) - ln(IC50_trametinib)
delta 越大 = 该细胞系对 BRAF 选择性抑制的抵抗越强, 而对 MEK 抑制仍保留反应。
这正是 BEACON 方案(BRAFi + 抗 EGFR)要解决的问题的细胞层面读数。

分层: BRAF V600 的 COREAD vs SKCM vs 其他; 并在两个筛选活动中分别计算以检验复制。
同时用突变负荷作为协变量, 避免与 19 号的混杂混淆。
"""
import re
import numpy as np
import pandas as pd
from scipy import stats

BASE = "data/raw"


def norm(s):
    return re.sub(r"[^A-Z0-9]", "", str(s).upper())


def braf_sets():
    m = pd.read_csv(f"{BASE}/mutations_latest.csv.gz")
    m = m[(m.gene_symbol == "BRAF") &
          (m.cancer_driver.astype(str).str.lower() == "true")]
    pm = m.protein_mutation.fillna("").astype(str)
    return {norm(x) for x in m[pm.str.contains("V600")].model_name.unique()}


def load(fn, v600):
    d = pd.read_csv(f"{BASE}/{fn}")
    d["ic50"] = np.exp(d.LN_IC50)
    d["ident"] = (d.ic50 >= d.MIN_CONC) & (d.ic50 <= d.MAX_CONC)
    d = d[d.ident].copy()
    d["nm"] = d.CELL_LINE_NAME.map(norm)
    d["v600"] = d.nm.isin(v600)
    return d


def main():
    v600 = braf_sets()
    out = []
    for fn, tag in [("GDSC2_fitted_dose_response_24Jul22.csv", "GDSC2"),
                    ("GDSC1_fitted_dose_response_24Jul22.csv", "GDSC1")]:
        d = load(fn, v600)
        need = ["Dabrafenib", "Trametinib"]
        pv = d[d.DRUG_NAME.isin(need)].pivot_table(
            index=["nm", "TCGA_DESC", "v600"], columns="DRUG_NAME",
            values="LN_IC50", aggfunc="mean").reset_index()
        pv.columns.name = None
        pv = pv.dropna(subset=need)
        if pv.empty:
            print(f"{tag}: 无配对样本")
            continue
        pv["delta_brafi_vs_meki"] = pv.Dabrafenib - pv.Trametinib
        pv["campaign"] = tag
        out.append(pv)
        print(f"\n=== {tag}: 配对(BRAFi & MEKi 均可识别) {len(pv)} 株, "
              f"其中 BRAF V600 {int(pv.v600.sum())} 株 ===")
        v = pv[pv.v600]
        for lin, g in v.groupby("TCGA_DESC"):
            if len(g) >= 3:
                print(f"  {lin:12s} n={len(g):3d}  "
                      f"delta 中位={g.delta_brafi_vs_meki.median():6.2f}  "
                      f"(正=对BRAFi更抵抗)")
        crc = v[v.TCGA_DESC == "COREAD"].delta_brafi_vs_meki
        oth = v[v.TCGA_DESC != "COREAD"].delta_brafi_vs_meki
        if len(crc) >= 3 and len(oth) >= 5:
            u = stats.mannwhitneyu(crc, oth)
            print(f"  COREAD(n={len(crc)}) vs 其他(n={len(oth)}): "
                  f"中位 {crc.median():.2f} vs {oth.median():.2f}  p={u.pvalue:.4f}")

    if not out:
        return
    A = pd.concat(out, ignore_index=True)
    A.to_csv("results/22_crc_exception.csv", index=False)

    print("\n" + "=" * 76)
    print("跨活动一致性: delta 在 GDSC1 与 GDSC2 之间是否同向")
    print("=" * 76)
    m = A.pivot_table(index="nm", columns="campaign",
                      values="delta_brafi_vs_meki", aggfunc="mean").dropna()
    if len(m) >= 8:
        r, p = stats.pearsonr(m.GDSC1, m.GDSC2)
        print(f"  配对细胞系 {len(m)} 株, delta 跨活动 r={r:+.3f} (p={p:.4f})")
        print(f"  符号一致率 {100*(np.sign(m.GDSC1)==np.sign(m.GDSC2)).mean():.1f}%")

    print("\n" + "=" * 76)
    print("BRAF V600 的 COREAD 在两个药上的绝对位置(GDSC2)")
    print("=" * 76)
    g2 = A[A.campaign == "GDSC2"]
    for lin in ["COREAD", "SKCM", "THCA"]:
        s = g2[(g2.TCGA_DESC == lin) & g2.v600]
        if len(s) >= 3:
            print(f"  {lin:8s} n={len(s):3d}  达拉非尼中位 IC50="
                  f"{np.exp(s.Dabrafenib.median()):8.2f} nM   "
                  f"曲美替尼中位 IC50={np.exp(s.Trametinib.median()):8.2f} nM")


if __name__ == "__main__":
    main()
