"""
09_interaction_test.py
病人层交互检验（主队列：coad_silu_2022 / AC-ICAM, n=348）
-------------------------------------------------
主估计量：score × 辅助化疗 的交互项，不是 score 主效应。
只有当"高耐药分数患者在化疗组更差、未化疗组无差异"时，才算抓到耐药而非预后。

为什么不用 TCGA：PHARMACEUTICAL_TX_ADJUVANT 仅 36/636 有值（20 YES / 16 NO），
无法做交互检验——这一点本身写进论文的局限性。

输出：results/table4_interaction.csv, results/09_interaction.json
"""
import json, warnings
import numpy as np
import pandas as pd
import requests
import statsmodels.api as sm
from statsmodels.duration.hazard_regression import PHReg
warnings.filterwarnings("ignore")
API = "https://www.cbioportal.org/api"


def clinical(study):
    r = requests.get(f"{API}/studies/{study}/clinical-data",
                     params={"clinicalDataType": "PATIENT"}, timeout=120)
    c = pd.DataFrame(r.json())
    return c.pivot_table(index="patientId", columns="clinicalAttributeId",
                         values="value", aggfunc="first")


def map_chemo(v):
    s = str(v).strip()
    if s.lower() in ("no", "none", "nan", ""):
        return 0
    if "Nno" in s:          # 数据字典中的 "Chemotherapy Nno" = 未化疗
        return 0
    if "chemo" in s.lower():
        return 1
    return np.nan


cl = clinical("coad_silu_2022")
d = pd.DataFrame(index=cl.index)
d["time"] = pd.to_numeric(cl.get("OS_MONTHS"), errors="coerce")
d["event"] = cl.get("OS_STATUS").astype(str).str.startswith("1").astype(int)
d["chemo"] = cl.get("ADJUVANT_TREATMENT").map(map_chemo)
d["stage"] = cl.get("AJCC_PATH_STAGE").astype(str).str.extract(r"(\d)")[0]
d["recur"] = (cl.get("DISTANT_RECURRENCE_STATUS").astype(str).str.lower() == "yes").astype(int)
d["age"] = pd.to_numeric(cl.get("AGE_AT_DIAGNOSIS", cl.get("AGE")), errors="coerce")
sex = cl.get("SEX", pd.Series(index=cl.index, dtype=object)).astype(str).str.upper()
d["male"] = sex.str.startswith("M").astype(float)

print("化疗映射:", d.chemo.value_counts(dropna=False).to_dict())
print("分期:", d.stage.value_counts(dropna=False).to_dict())

sc = pd.read_csv("results/patient_scores.csv")
s = sc[sc.study == "coad_silu_2022"].set_index("sampleId")
s.index = [i.rsplit("-PT", 1)[0] for i in s.index]
s = s[~s.index.duplicated()]
d = d.join(s, how="inner")
print(f"合并后 n={len(d)}  死亡={int(d.event.sum())}  化疗={int(d.chemo.sum())}")

d23 = d[d.stage.isin(["2", "3"])].copy()
print(f"II/III 期 n={len(d23)}  死亡={int(d23.event.sum())}  化疗={int(d23.chemo.sum())}"
      f"  远处复发={int(d23.recur.sum())}")


def cox(df, score_col, covars, label):
    x = df.dropna(subset=["time", "event", "chemo", score_col]).copy()
    x = x[x.time > 0]
    if len(x) < 40 or x.event.sum() < 15:
        return None
    x["s"] = (x[score_col] - x[score_col].mean()) / x[score_col].std()
    x["int"] = x["s"] * x["chemo"]
    cols = ["s", "chemo", "int"] + covars
    X = sm.add_constant(x[cols].astype(float))
    try:
        m = PHReg(x["time"].values, X.values,
                  status=x["event"].values.astype(int), ties="efron")
        r = m.fit(disp=0)
    except Exception as e:
        print("  Cox 失败:", e); return None
    nm = list(X.columns); ci = r.conf_int()
    o = {}
    for t in ["s", "chemo", "int"]:
        i = nm.index(t)
        o[t] = dict(HR=round(float(np.exp(r.params[i])), 3),
                    lo=round(float(np.exp(ci[i][0])), 3),
                    hi=round(float(np.exp(ci[i][1])), 3),
                    p=round(float(r.pvalues[i]), 4))
    o.update(n=int(len(x)), events=int(x.event.sum()), chemo_n=int(x.chemo.sum()),
             label=label, score=score_col)
    return o


rows, detail = [], {}
for tag, sub in [("stage II-III", d23), ("all stages", d)]:
    covars = [c for c in ["age", "male", "stage"] if sub[c].notna().sum() > 0.7 * len(sub)]
    for scol in ["score_5FU", "score_random", "score_prolif"]:
        r = cox(sub, scol, covars, tag)
        if not r:
            continue
        rows.append({"subset": tag, "score": scol, "n": r["n"], "events": r["events"],
                     "n_chemo": r["chemo_n"],
                     "HR_score": r["s"]["HR"], "p_score": r["s"]["p"],
                     "HR_chemo": r["chemo"]["HR"], "p_chemo": r["chemo"]["p"],
                     "HR_int": r["int"]["HR"], "int_lo": r["int"]["lo"],
                     "int_hi": r["int"]["hi"], "p_int": r["int"]["p"]})
        detail[f"{tag}|{scol}"] = r
        print(f"  {scol:12s} {tag:13s} n={r['n']:3d} ev={r['events']:3d} "
              f"| HR_int={r['int']['HR']:.2f} [{r['int']['lo']:.2f}-{r['int']['hi']:.2f}] "
              f"p={r['int']['p']:.3f} | HR_score={r['s']['HR']:.2f} p={r['s']['p']:.3f}")

tab = pd.DataFrame(rows)
tab.to_csv("results/table4_interaction.csv", index=False)
print("\n=== TABLE 4 ===")
print(tab.to_string(index=False))

# 可解释展示：分数三分位 × 化疗 的死亡率
d23["tert"] = pd.qcut(d23.score_5FU, 3, labels=["low", "mid", "high"])
xt = d23.groupby(["tert", "chemo"], observed=True)["event"].agg(["size", "sum", "mean"]).round(3)
print("\n=== 分数三分位 × 化疗 的死亡比例 ===")
print(xt.to_string())
xt.to_csv("results/table5_stratified_event_rates.csv")

json.dump(detail, open("results/09_interaction.json", "w", encoding="utf-8"),
          indent=2, ensure_ascii=False)
