"""
23_patient_layer_braf.py
------------------------
把细胞系的结论落到病人: BRAF V600 的结直肠癌患者在真实队列里是什么样。

队列: AC-ICAM / coad_silu_2022 (n=348, 有真实辅助化疗标注与 OS)
终点: 总生存; 并检验 BRAF 状态 x 辅助化疗 的交互(该不该用化疗替代靶向)
补充: TCGA COADREAD 作为第二个队列的一致性检查
"""
import numpy as np
import pandas as pd
import requests
from scipy import stats
from lifelines import CoxPHFitter
from lifelines.statistics import logrank_test

API = "https://www.cbioportal.org/api"


def get(path, **kw):
    r = requests.get(f"{API}{path}", params=kw, timeout=180)
    r.raise_for_status()
    return r.json()


def clinical(study):
    cd = pd.DataFrame(get(f"/studies/{study}/clinical-data",
                          clinicalDataType="PATIENT"))
    return cd.pivot_table(index="patientId", columns="clinicalAttributeId",
                          values="value", aggfunc="first")


def mutations(study, genes):
    prof = [p for p in get(f"/studies/{study}/molecular-profiles")
            if p.get("molecularAlterationType") == "MUTATION_EXTENDED"]
    if not prof:
        return None
    r = requests.post(f"{API}/molecular-profiles/{prof[0]['molecularProfileId']}"
                             f"/mutations/fetch",
                      json={"sampleListId": f"{study}_all"}, timeout=300)
    r.raise_for_status()
    m = pd.DataFrame(r.json())
    gmap = {}
    for g in genes:
        try:
            info = get(f"/genes/{g}")
            gmap[int(info["entrezGeneId"])] = info["hugoGeneSymbol"]
        except Exception:
            pass
    m["entrezGeneId"] = m.entrezGeneId.astype(int)
    m = m[m.entrezGeneId.isin(gmap)].copy()
    m["gene"] = m.entrezGeneId.map(gmap)
    return m


def main():
    study = "coad_silu_2022"
    cw = clinical(study)
    m = mutations(study, ["BRAF"])
    m["v600"] = m.proteinChange.fillna("").str.contains("V600")
    braf_v600 = set(m.loc[m.v600, "patientId"])
    braf_any = set(m.loc[m.gene == "BRAF", "patientId"])
    print(f"队列 {study}: {len(cw)} 例")
    print(f"  BRAF 任意突变 {len(braf_any)} 例; 其中 V600 {len(braf_v600)} 例")

    d = pd.DataFrame(index=cw.index)
    d["os_m"] = pd.to_numeric(cw.get("OS_MONTHS"), errors="coerce")
    d["event"] = cw.get("OS_STATUS").astype(str).str.upper().str.startswith("1").astype(int)
    d["arm"] = cw.get("ADJUVANT_TREATMENT")
    d["stage"] = pd.to_numeric(cw.get("AJCC_PATH_STAGE"), errors="coerce")
    agecol = "AGE_AT_DX" if "AGE_AT_DX" in cw else "AGE"
    d["age"] = pd.to_numeric(cw.get(agecol), errors="coerce")
    d["braf"] = d.index.isin(braf_v600).astype(int)
    d = d.dropna(subset=["os_m"])
    d = d[d.os_m > 0]
    d["chemo"] = (d.arm.astype(str).str.lower() != "no").astype(int)
    # 只在 I-III 期里看辅助化疗, IV 期不适用
    dz = d[(d.stage <= 3) & d.stage.notna()].copy()
    print(f"\n可分析 {len(d)} 例 (死亡 {int(d.event.sum())}); "
          f"I-III 期 {len(dz)} 例 (死亡 {int(dz.event.sum())})")
    print(f"  BRAF V600 在 I-III 期中 {int(dz.braf.sum())} 例")

    print("\n" + "=" * 72)
    print("1. BRAF V600 的预后效应 (I-III 期, 校正分期与年龄)")
    print("=" * 72)
    for tag, s in [("全套", dz), ("仅化疗组", dz[dz.chemo == 1]),
                   ("仅未化疗组", dz[dz.chemo == 0])]:
        if s.braf.sum() < 5:
            print(f"  {tag}: BRAF V600 仅 {int(s.braf.sum())} 例, 跳过")
            continue
        x = s[["os_m", "event", "braf", "stage", "age"]].dropna()
        if len(x) < 15 or x.braf.sum() < 5 or x.braf.nunique() < 2:
            print(f"  {tag}: 有效样本不足 (n={len(x)}), 跳过")
            continue
        try:
            cph = CoxPHFitter()
            cph.fit(x, duration_col="os_m", event_col="event")
            r = cph.summary.loc["braf"]
        except Exception as e:
            print(f"  {tag}: 拟合失败 {type(e).__name__}, 跳过")
            continue
        print(f"  {tag:12s} n={len(x):3d} 突变{int(x.braf.sum()):2d} 例  "
              f"HR={r['exp(coef)']:.2f} ({r['exp(coef) lower 95%']:.2f}-"
              f"{r['exp(coef) upper 95%']:.2f})  p={r['p']:.4f}")

    print("\n" + "=" * 72)
    print("2. BRAF 状态 x 辅助化疗 交互 (BRAF 突变者该不该用化疗)")
    print("=" * 72)
    x = dz[["os_m", "event", "braf", "chemo", "stage", "age"]].dropna().copy()
    x["inter"] = x.braf * x.chemo
    print(f"  交互模型样本 {len(x)}, BRAF 突变 {int(x.braf.sum())}")
    cph = CoxPHFitter()
    cph.fit(x, duration_col="os_m", event_col="event")
    print(cph.summary[["exp(coef)", "exp(coef) lower 95%",
                       "exp(coef) upper 95%", "p"]].to_string(
        float_format=lambda v: f"{v:.3f}"))

    print("\n  分层死亡率(BRAF V600 组内):")
    b = dz[dz.braf == 1]
    for tag, s in [("化疗", b[b.chemo == 1]), ("未化疗", b[b.chemo == 0])]:
        if len(s):
            print(f"    {tag:6s} n={len(s):2d}  死亡 {int(s.event.sum()):2d} "
                  f"({100*s.event.mean():.0f}%)  中位随访 OS "
                  f"{s.os_m.median():.1f} 月")
    nb = dz[dz.braf == 0]
    for tag, s in [("化疗", nb[nb.chemo == 1]), ("未化疗", nb[nb.chemo == 0])]:
        print(f"    BRAF-WT {tag:6s} n={len(s):3d}  死亡 {int(s.event.sum()):3d} "
              f"({100*s.event.mean():.0f}%)")

    print("\n" + "=" * 72)
    print("3. 第二队列一致性: TCGA COADREAD")
    print("=" * 72)
    try:
        tw = clinical("coadread_tcga")
        tm = mutations("coadread_tcga", ["BRAF"])
        tv = set(tm.loc[tm.proteinChange.fillna("").str.contains("V600"),
                        "patientId"])
        tany = set(tm.loc[tm.gene == "BRAF", "patientId"])
        print(f"  TCGA COADREAD: {len(tw)} 例, BRAF 突变 {len(tany)}, V600 {len(tv)}")
        t = pd.DataFrame(index=tw.index)
        t["os_m"] = pd.to_numeric(tw.get("OS_MONTHS"), errors="coerce")
        t["event"] = tw.get("OS_STATUS").astype(str).str.upper().str.startswith("1").astype(int)
        t["braf"] = t.index.isin(tv).astype(int)
        t = t.dropna(subset=["os_m"])
        t = t[t.os_m > 0]
        if t.braf.sum() >= 5:
            lr = logrank_test(t[t.braf == 1].os_m, t[t.braf == 0].os_m,
                              t[t.braf == 1].event, t[t.braf == 0].event)
            print(f"  BRAF V600(n={int(t.braf.sum())}) vs WT: "
                  f"log-rank p={lr.p_value:.2e}")
            print(f"    中位 OS: {t[t.braf==1].os_m.median():.1f} vs "
                  f"{t[t.braf==0].os_m.median():.1f} 月")
            cph = CoxPHFitter()
            xx = t[["os_m", "event", "braf"]]
            cph.fit(xx, duration_col="os_m", event_col="event")
            r = cph.summary.loc["braf"]
            print(f"    未校正 HR={r['exp(coef)']:.2f} "
                  f"({r['exp(coef) lower 95%']:.2f}-"
                  f"{r['exp(coef) upper 95%']:.2f}) p={r['p']:.2e}")
    except Exception as e:
        print(f"  TCGA 检查失败: {e}")

    x.to_csv("results/23_patient_braf_analysis.csv", index=False)

    print("\n" + "=" * 72)
    print("4. BRAF V600 与分子分型 CMS 的关系(若测到, 可解释为何 CRC 例外)")
    print("=" * 72)
    if "CMS" in cw:
        cm = cw.CMS.reindex(dz.index)
        print(pd.crosstab(dz.braf, cm, dropna=False).to_string())
    if "TUMOR_ANATOMIC_LOCATION" in cw:
        loc = cw.TUMOR_ANATOMIC_LOCATION.reindex(dz.index)
        print("\n  原发部位(右半 vs 左半)与 BRAF:")
        print(pd.crosstab(dz.braf, loc, dropna=False).to_string())


if __name__ == "__main__":
    main()
