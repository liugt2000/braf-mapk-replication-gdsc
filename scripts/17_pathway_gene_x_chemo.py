"""
17_pathway_gene_x_chemo.py
--------------------------
以"某个信号通路基因导致耐药"为主角, 在唯一具备真实辅助化疗标注的公开队列
(AC-ICAM / coad_silu_2022, n=348)中做预设基因 x 化疗 的交互检验。

预设候选(不做全基因组搜索, 避免数据挖掘):
  TGF-beta 轴: SMAD4, SMAD2, TGFBR2
  MAPK 轴:     KRAS, BRAF, NRAS
  PI3K 轴:     PIK3CA, PTEN
  p53 轴:      TP53
  Wnt 轴:      APC, CTNNB1, RNF43
  MMR:         MLH1, MSH2, MSH6, PMS2

判据(事先写死): 只有 interaction HR 显著且方向为"突变组化疗获益更小"时,
才可宣称该基因导致化疗耐药。主效应显著不算。
"""
import numpy as np
import pandas as pd
import requests

API = "https://www.cbioportal.org/api"
STUDY = "coad_silu_2022"
CAND = ["SMAD4", "SMAD2", "TGFBR2", "KRAS", "BRAF", "NRAS", "PIK3CA",
        "PTEN", "TP53", "APC", "CTNNB1", "RNF43", "MLH1", "MSH2", "MSH6", "PMS2"]


def get(path, **kw):
    r = requests.get(f"{API}{path}", params=kw, timeout=180)
    r.raise_for_status()
    return r.json()


def main():
    # ---- 临床 ----
    cd = pd.DataFrame(get(f"/studies/{STUDY}/clinical-data",
                          clinicalDataType="PATIENT"))
    cw = cd.pivot_table(index="patientId", columns="clinicalAttributeId",
                        values="value", aggfunc="first")
    print(f"队列 {STUDY}: {len(cw)} 例")
    for c in ["ADJUVANT_TREATMENT", "OS_MONTHS", "OS_STATUS",
              "AJCC_PATH_STAGE", "SEX", "AGE"]:
        if c in cw:
            print(f"  {c}: {cw[c].value_counts(dropna=False).head(5).to_dict()}")

    # ---- 突变 ----
    prof = [p for p in get(f"/studies/{STUDY}/molecular-profiles")
            if p.get("molecularAlterationType") == "MUTATION_EXTENDED"]
    pid = prof[0]["molecularProfileId"]
    print(f"\n突变谱: {pid}")
    rm = requests.post(
        f"{API}/molecular-profiles/{pid}/mutations/fetch",
        json={"sampleListId": f"{STUDY}_all"}, timeout=300)
    rm.raise_for_status()
    muts = pd.DataFrame(rm.json())
    print(f"  突变记录 {len(muts):,}")
    # entrez -> hugo 映射
    gmap = {}
    for g in CAND:
        try:
            info = get(f"/genes/{g}")
            gmap[int(info["entrezGeneId"])] = info["hugoGeneSymbol"]
        except Exception:
            pass
    muts["entrezGeneId"] = muts.entrezGeneId.astype(int)
    muts = muts[muts.entrezGeneId.isin(gmap)]
    muts["gene"] = muts.entrezGeneId.map(gmap)
    muts = muts.drop_duplicates(["patientId", "gene"])
    muts["patient"] = muts.patientId
    print("  预设基因突变例数:")
    alt = {}
    for g in CAND:
        s = set(muts.loc[muts.gene == g, "patient"])
        alt[g] = s
        print(f"    {g:8s} {len(s):3d}")

    # ---- 生存 ----
    d = pd.DataFrame(index=cw.index)
    d["arm"] = cw.get("ADJUVANT_TREATMENT")
    d["os_m"] = pd.to_numeric(cw.get("OS_MONTHS"), errors="coerce")
    st = cw.get("OS_STATUS").astype(str).str.upper()
    d["event"] = st.str.startswith("1").astype(int)
    d["stage"] = cw.get("AJCC_PATH_STAGE")
    d["age"] = pd.to_numeric(cw.get("AGE"), errors="coerce")
    d = d.dropna(subset=["arm", "os_m"])
    d = d[d.os_m > 0]
    print(f"\n有化疗标注 + 生存: {len(d)} 例  死亡 {int(d.event.sum())}")
    d["chemo"] = (d.arm.astype(str).str.lower() != "no").astype(int)
    print(f"  接受化疗: {int(d.chemo.sum())}  未化疗: {int((1-d.chemo).sum())}")

    # ---- Cox 交互 ----
    try:
        from lifelines import CoxPHFitter
        from lifelines.statistics import logrank_test
    except ImportError:
        import subprocess
        subprocess.run(["python", "-m", "pip", "install", "--quiet", "lifelines"])
        from lifelines import CoxPHFitter
        from lifelines.statistics import logrank_test

    rows = []
    for g in CAND:
        d["alt"] = d.index.isin(alt[g]).astype(int)
        if d.alt.sum() < 12:
            continue
        x = d.copy()
        x["inter"] = x.alt * x.chemo
        x = x[["os_m", "event", "alt", "chemo", "inter"]].dropna()
        try:
            cph = CoxPHFitter()
            cph.fit(x, duration_col="os_m", event_col="event")
            s = cph.summary
            rows.append(dict(gene=g, n_alt=int(x.alt.sum()), n=len(x),
                             hr_alt=s.loc["alt", "exp(coef)"],
                             p_alt=s.loc["alt", "p"],
                             hr_int=s.loc["inter", "exp(coef)"],
                             p_int=s.loc["inter", "p"],
                             ci_lo=s.loc["inter", "exp(coef) lower 95%"],
                             ci_hi=s.loc["inter", "exp(coef) upper 95%"]))
        except Exception as e:
            print(f"  {g} 拟合失败: {e}")

    r = pd.DataFrame(rows).sort_values("p_int")
    r.to_csv("results/17_pathway_gene_x_chemo.csv", index=False)
    print("\n" + "=" * 78)
    print("基因 x 化疗 交互 (Cox, 终点 OS)   HR_int<1 = 突变组从化疗获益更多;")
    print("                                    HR_int>1 = 突变组化疗获益更少(即耐药)")
    print("=" * 78)
    print(r.to_string(index=False, float_format=lambda v: f"{v:.3f}"))

    if len(r):
        hit = r[(r.p_int < 0.05)]
        print(f"\n交互 p<0.05 的基因: {len(hit)} 个" +
              (f" -> {hit.gene.tolist()}" if len(hit) else " (无)"))


if __name__ == "__main__":
    main()
