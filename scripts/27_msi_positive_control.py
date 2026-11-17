"""
27. 阳性对照：结直肠癌 panel 能否复现 MSI-H 与氟尿嘧啶辅助化疗无效这一临床定论（Ribic 2003）。
只用可识别 IC50；报告两个活动各自的方向、效应量、p 值、最小可检出效应。
"""
import re

import numpy as np
import pandas as pd
from scipy import stats

BASE = "data/raw"


def norm(s):
    return re.sub(r"[^A-Z0-9]", "", str(s).upper())


def load(camp, fname):
    d = pd.read_csv(f"{BASE}/{fname}")
    d["ic50"] = np.exp(d.LN_IC50)
    d["ident"] = (d.ic50 >= d.MIN_CONC) & (d.ic50 <= d.MAX_CONC)
    return d[(d.TCGA_DESC == "COREAD") & d.ident].copy()


def mde_d(n1, n2, alpha=0.05, power=0.8):
    """两样本 t 检验的最小可检出标准化效应量（迭代）"""
    from scipy.stats import nct, t as tdist
    df = n1 + n2 - 2
    tc = tdist.ppf(1 - alpha / 2, df)
    lo, hi = 0.01, 10.0
    for _ in range(100):
        mid = (lo + hi) / 2
        ncp = mid * np.sqrt(n1 * n2 / (n1 + n2))
        p = nct.sf(tc, df, ncp) + nct.cdf(-tc, df, ncp)
        if p < power:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2


ml = pd.read_csv(f"{BASE}/model_list_latest.csv.gz")
msi = {norm(r.model_name): r.msi_status for r in ml.itertuples()}

rows = []
for camp, fname in [("GDSC2", "GDSC2_fitted_dose_response_24Jul22.csv"),
                    ("GDSC1", "GDSC1_fitted_dose_response_24Jul22.csv")]:
    c = load(camp, fname)
    c["msi"] = c.CELL_LINE_NAME.map(lambda x: msi.get(norm(x)))
    for drug in ["5-Fluorouracil", "SN-38", "Irinotecan", "Oxaliplatin"]:
        q = (c[c.DRUG_NAME == drug].dropna(subset=["msi"])
             .drop_duplicates("CELL_LINE_NAME"))
        a = q[q.msi == "MSI"].LN_IC50.values
        b = q[q.msi == "MSS"].LN_IC50.values
        if len(a) < 4 or len(b) < 4:
            rows.append(dict(campaign=camp, drug=drug, n_msi=len(a), n_mss=len(b),
                             fold=np.nan, p_welch=np.nan, p_mwu=np.nan,
                             d=np.nan, mde_d=np.nan, note="样本不足"))
            continue
        sp = np.sqrt(((len(a) - 1) * a.var(ddof=1) + (len(b) - 1) * b.var(ddof=1))
                     / (len(a) + len(b) - 2))
        rows.append(dict(
            campaign=camp, drug=drug, n_msi=len(a), n_mss=len(b),
            fold=float(np.exp(np.median(a) - np.median(b))),
            p_welch=float(stats.ttest_ind(a, b).pvalue),
            p_mwu=float(stats.mannwhitneyu(a, b).pvalue),
            d=float((a.mean() - b.mean()) / sp),
            mde_d=float(mde_d(len(a), len(b))),
            mde_ln=float(mde_d(len(a), len(b)) * sp),
            note=""))

R = pd.DataFrame(rows)
R.to_csv("results/27_msi_positive_control.csv", index=False)

pd.set_option("display.width", 160)
print("=== 阳性对照：MSI vs MSS，可识别 IC50 ===\n")
print(R.to_string(index=False, float_format=lambda v: f"{v:.3f}"))

print("\n=== 判读 ===")
f = R[(R.drug == "5-Fluorouracil") & R.fold.notna()]
if len(f) == 2:
    dirs = np.sign(np.log(f.fold.values))
    print(f"5-FU 中位倍数: GDSC2 {f.fold.iloc[0]:.2f}x, GDSC1 {f.fold.iloc[1]:.2f}x")
    print(f"方向一致: {'是' if dirs[0] == dirs[1] else '否（符号翻转）'}")
    print(f"任一活动 p<0.05: {'是' if (f.p_welch < 0.05).any() else '否'}")
    mde_ln = float(f.mde_ln.mean())
    print(f"最小可检出效应 (80% 功效, alpha=0.05): d≈{f.mde_d.mean():.2f} 标准化单位 "
          f"= {mde_ln:.2f} ln 单位 ≈ {np.exp(mde_ln):.1f} 倍 IC50")
print(f"\nCOREAD 模型 MSI/MSS: {sum(1 for k in msi if msi[k]=='MSI')} / "
      f"{sum(1 for k in msi if msi[k]=='MSS')}（全 panel）")
