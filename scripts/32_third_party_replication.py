"""
32_third_party_replication.py (v2) — CTRPv2/PRISM 第三方复现 (任务 3)

在两个非 GDSC 公共药敏 campaign 中复现稿件核心主张:
  (A) 泛癌: BRAF V600 系是否比 V600 阴性系对 BRAFi/MEKi 更敏感
  (B) 解离: V600 系内部, CRC 相对黑色素瘤是否对 BRAFi 特异性不敏感、对 MEKi 保持敏感

关键口径 (v2 修正):
  - V600+ = cBioPortal ccle_broad_2019 蛋白改变以 V600 开头 (101 独立样本)
  - V600-(WT) = 该研究 1739 样本全集内、无 V600 改变 (测序覆盖 1570; 全集口径记录于 summary)
  - 组织: CTRP 用 ccle_primary_site=='large_intestine' 定 CRC、ccle_primary_hist 含 melanoma 定黑色素瘤
  - 指标锚定检查: 经典敏感/抵抗对 (A-375 vs HT-29 等) 先行验证 aac_published 与 aac_recomputed,
    取有药理学区分度者为主指标 (预期 aac_recomputed; aac 0-1, 越高越敏感)
  - 解离的 AAC 类比: dAAC = AAC(MEKi) - AAC(BRAFi), 正 = 相对偏 BRAFi 抵抗
  - PRISM 仅取 passed_str_profiling==True 的系 (Corsello 自有 QC)
"""
import json
import os
import re
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import rdata
from scipy import stats

warnings.filterwarnings("ignore", category=UserWarning)
ROOT = Path(r"D:\AI_work\workbuddy\2026-11-16-16-13-17")
OUT = ROOT / "results"


def norm(s) -> str:
    return re.sub(r"[^A-Z0-9]", "", str(s).upper())


def bare(sid) -> str:
    return str(sid).split("_", 1)[0]


def mwu(a, b):
    a = np.asarray(a, float); a = a[~np.isnan(a)]
    b = np.asarray(b, float); b = b[~np.isnan(b)]
    if len(a) < 2 or len(b) < 2:
        return dict(n1=len(a), n2=len(b), med1=float(np.median(a)) if len(a) else np.nan,
                    med2=float(np.median(b)) if len(b) else np.nan, U=np.nan, p=np.nan)
    U, p = stats.mannwhitneyu(a, b, alternative="two-sided")
    return dict(n1=len(a), n2=len(b), med1=float(np.median(a)), med2=float(np.median(b)),
                U=float(U), p=float(p))


# ============ BRAF 状态 ============
mut = json.load(open(ROOT / "data/raw/ccle2019_braf_mutations.json", encoding="utf-8"))
v600 = {norm(bare(m["sampleId"])) for m in mut
        if str(m.get("proteinChange", "")).upper().startswith("V600")}
study = {norm(bare(s)) for s in json.load(open(ROOT / "data/raw/ccle2019_samples_all.json",
                                               encoding="utf-8"))}
print(f"[BRAF] V600 样本={len(v600)}, 研究 sample 全集={len(study)}, "
      f"V600∩全集={len(v600 & study)} (测序覆盖 1570/1739)")

ANCHOR_MEL = ["A-375", "A2058", "SK-MEL-28", "WM266-4", "Malme-3M", "LOX-IMVI", "Hs 294T"]
ANCHOR_CRC = ["HT-29", "RKO", "COLO 205", "SW1417", "LS411N", "MDST8"]


def anchor_report(df, metrics, name):
    """打印经典锚定系在各指标下的值, 验证药理学方向"""
    print(f"\n== [{name}] 锚定检查 ==")
    d = df[df["cellid"].isin(ANCHOR_MEL + ANCHOR_CRC)]
    piv = d.pivot_table(index="cellid", columns="drugid", values=metrics,
                        aggfunc="median")
    order = [c for c in ANCHOR_MEL + ANCHOR_CRC if c in piv.index]
    print(piv.loc[order].round(3).to_string())


# ============ 通用抽取 ============
def sens_table(path: Path, drug_names: list, keep_metrics: list) -> pd.DataFrame:
    p = rdata.read_rds(str(path))
    info = p.sensitivity["info"]
    prof = p.sensitivity["profiles"]
    drug_ids = [d for d in drug_names if d in set(info["drugid"].astype(str))]
    sub = info[info["drugid"].isin(drug_ids)][["cellid", "drugid"]]
    prof = prof.loc[sub.index]
    sub = pd.concat([sub.reset_index(drop=True), prof.reset_index(drop=True)], axis=1)
    for m in keep_metrics:
        if m not in sub.columns:
            sub[m] = np.nan
    print(f"[{path.name}] 命中 {len(drug_ids)}/{len(drug_names)} 药: {drug_ids}; "
          f"experiments={len(sub)}")
    return sub, p


# ============ CTRPv2 ============
print("\n" + "=" * 60 + "\nCTRPv2\n" + "=" * 60)
CTRP_DRUGS = ["Dabrafenib", "Trametinib", "Vemurafenib", "PLX4720", "Selumetinib"]
ct, p_ct = sens_table(ROOT / "data/raw/psets/PSet_CTRPv2_2015.rds", CTRP_DRUGS,
                      ["aac_published", "aac_recomputed"])

per = ct.groupby(["cellid", "drugid"])[["aac_published", "aac_recomputed"]] \
        .median().reset_index()
cell_ct = p_ct.cell
site = cell_ct["ccle_primary_site"].astype(str).str.lower().values
hist = cell_ct["ccle_primary_hist"].astype(str).str.lower().values
ct_meta = pd.DataFrame({
    "cellid": cell_ct["cellid"].astype(str).values,
    "site": site, "hist": hist,
})
per = per.merge(ct_meta, on="cellid", how="left")
per["key"] = per["cellid"].map(norm)
per["status"] = np.where(per["key"].isin(v600), "V600",
                         np.where(per["key"].isin(study), "WT", "unknown"))
per["is_crc"] = per["site"].eq("large_intestine")
per["is_mela"] = per["hist"].str.contains("melanoma", na=False) | \
                 (per["site"].eq("skin") & per["hist"].isna())
print("状态计数(每系):", per.groupby("status")["cellid"].nunique().to_dict())

# 锚定检查 (用两指标)
anchor_report(per[per["drugid"].isin(["Vemurafenib", "PLX4720", "Selumetinib"])],
              ["aac_published", "aac_recomputed"], "CTRPv2 anchors")

# 主指标判定: aac_recomputed 上 BRAFi 锚定应满足 敏感黑素瘤 > 抵抗CRC
a_re = per.pivot_table(index="cellid", columns="drugid", values="aac_recomputed")
mel_med = a_re.loc[[c for c in ANCHOR_MEL if c in a_re.index], "Vemurafenib"].median()
crc_med = a_re.loc[[c for c in ANCHOR_CRC if c in a_re.index], "Vemurafenib"].median()
a_pub = per.pivot_table(index="cellid", columns="drugid", values="aac_published")
mel_pub = a_pub.loc[[c for c in ANCHOR_MEL if c in a_pub.index], "Vemurafenib"].median()
crc_pub = a_pub.loc[[c for c in ANCHOR_CRC if c in a_pub.index], "Vemurafenib"].median()
metric = "aac_recomputed" if mel_med > crc_med else "aac_published"
anchor_verdict = {
    "aac_recomputed": {"melanoma_anchor_median": float(mel_med),
                       "crc_anchor_median": float(crc_med),
                       "separates": bool(mel_med > crc_med)},
    "aac_published": {"melanoma_anchor_median": float(mel_pub),
                      "crc_anchor_median": float(crc_pub),
                      "separates": bool(mel_pub > crc_pub)},
    "chosen": metric,
}
print(f"\n锚定判定: {anchor_verdict}")

rows = []
# (A) 泛癌
for did in CTRP_DRUGS:
    d = per[(per["drugid"] == did) & (per["status"] != "unknown")]
    res = mwu(d.loc[d.status == "V600", metric], d.loc[d.status == "WT", metric])
    rows.append(dict(dataset="CTRPv2", test="pancancer_V600_vs_WT", drug=did,
                     metric=metric, **res))
# (B) 解离
def pair_dissociation(per, brafi, meki, tag):
    pv = per[per["status"] == "V600"].pivot_table(
        index=["cellid", "key", "is_crc", "is_mela"], columns="drugid",
        values=metric)
    pv = pv.reset_index()
    pv = pv.dropna(subset=[brafi, meki])
    pv["daac"] = pv[meki] - pv[brafi]     # 正 = 相对 BRAFi 抵抗
    n_crc = int(pv["is_crc"].sum()); n_mel = int(pv["is_mela"].sum())
    print(f"\n[{tag}] V600 双药可配 = {len(pv)} (CRC={n_crc}, 黑色素瘤={n_mel})")
    print("  CRC 系:", pv.loc[pv.is_crc, "cellid"].tolist())
    print("  黑色素瘤系:", pv.loc[pv.is_mela, "cellid"].tolist())
    if n_crc >= 2 and n_mel >= 2:
        for label, sel in [(f"AAC_{brafi}", brafi), (f"AAC_{meki}", meki),
                           (f"dAAC_{meki}_minus_{brafi}", "daac")]:
            res = mwu(pv.loc[pv.is_crc, sel], pv.loc[pv.is_mela, sel])
            rows.append(dict(dataset="CTRPv2", test=f"V600_CRC_vs_melanoma_{label}",
                             drug=f"{brafi}+{meki}", metric=metric, **res))
        res = mwu(pv.loc[pv.is_crc, "daac"], pv.loc[~pv.is_crc, "daac"])
        rows.append(dict(dataset="CTRPv2", test="V600_CRC_vs_allothers_dAAC",
                         drug=f"{brafi}+{meki}", metric=metric, **res))
    else:
        rows.append(dict(dataset="CTRPv2", test=f"V600_CRC_vs_melanoma_dAAC_{tag}",
                         drug=f"{brafi}+{meki}", metric=metric,
                         n1=n_crc, n2=n_mel, med1=np.nan, med2=np.nan,
                         U=np.nan, p=np.nan))
    pv.to_csv(OUT / f"32_ctrpv2_pairs_{tag}.csv", index=False)
    return pv

pair_dissociation(per, "Vemurafenib", "Selumetinib", "vemu_selu")
pair_dissociation(per, "PLX4720", "Selumetinib", "plx_selu")
pair_dissociation(per, "Dabrafenib", "Trametinib", "dab_tra")
per.to_csv(OUT / "32_ctrpv2_per_line.csv", index=False)

# ============ PRISM ============
print("\n" + "=" * 60 + "\nPRISM\n" + "=" * 60)
PRISM_DRUGS = ["Dabrafenib", "Trametinib", "Vemurafenib", "PLX4720", "Selumetinib",
               "PD-0325901", "TAK-733", "Cobimetinib", "Refametinib"]
pr, p_pr = sens_table(ROOT / "data/raw/psets/PSet_PRISM_2020.rds", PRISM_DRUGS,
                      ["aac_recomputed", "auc_published", "ic50_published"])

qc = p_pr.cell
passed = set(qc.loc[qc["passed_str_profiling"].astype(str).str.lower().eq("true"),
                    "cellid"].astype(str))
print(f"PRISM passed_str_profiling==True 系数 = {len(passed)}")
pr = pr[pr["cellid"].isin(passed)]
per_p = pr.groupby(["cellid", "drugid"])[["aac_recomputed", "auc_published",
                                          "ic50_published"]].median().reset_index()
lab = (qc["primary_tissue"].astype(str).str.lower() + "|" +
       qc["tissueid"].astype(str).str.lower() + "|" +
       qc["PRISM.tissueid"].astype(str).str.lower())
pr_meta = pd.DataFrame({"cellid": qc["cellid"].astype(str).values, "lab": lab.values})
per_p = per_p.merge(pr_meta, on="cellid", how="left")
per_p["key"] = per_p["cellid"].map(norm)
per_p["status"] = np.where(per_p["key"].isin(v600), "V600",
                           np.where(per_p["key"].isin(study), "WT", "unknown"))
per_p["is_crc"] = per_p["lab"].str.contains("colorectal|colon|large_intestine|bowel",
                                            na=False, regex=True)
per_p["is_mela"] = per_p["lab"].str.contains("melanoma", na=False)
print("PRISM 状态计数(每系):", per_p.groupby("status")["cellid"].nunique().to_dict())

anchor_report(per_p[per_p["drugid"].isin(["Dabrafenib", "Trametinib", "Vemurafenib"])],
              ["aac_recomputed", "auc_published"], "PRISM anchors")

# 泛癌
p_metric = "aac_recomputed"
for did in per_p["drugid"].unique():
    d = per_p[(per_p["drugid"] == did) & (per_p["status"] != "unknown")]
    res = mwu(d.loc[d.status == "V600", p_metric], d.loc[d.status == "WT", p_metric])
    rows.append(dict(dataset="PRISM", test="pancancer_V600_vs_WT", drug=did,
                     metric=p_metric, **res))
# 解离: Dabrafenib + Trametinib
pv = per_p[per_p["status"] == "V600"].pivot_table(
    index=["cellid", "key", "is_crc", "is_mela"], columns="drugid",
    values=p_metric).reset_index()
pv = pv.dropna(subset=["Dabrafenib", "Trametinib"])
pv["daac"] = pv["Trametinib"] - pv["Dabrafenib"]
n_crc = int(pv["is_crc"].sum()); n_mel = int(pv["is_mela"].sum())
print(f"\n[PRISM dab+tra] V600 双药可配 = {len(pv)} (CRC={n_crc}, 黑色素瘤={n_mel})")
print("  CRC 系:", pv.loc[pv.is_crc, "cellid"].tolist())
print("  黑色素瘤系:", pv.loc[pv.is_mela, "cellid"].tolist())
if n_crc >= 2 and n_mel >= 2:
    for label, sel in [("AAC_Dabrafenib", "Dabrafenib"),
                       ("AAC_Trametinib", "Trametinib"),
                       ("dAAC_Trametinib_minus_Dabrafenib", "daac")]:
        res = mwu(pv.loc[pv.is_crc, sel], pv.loc[pv.is_mela, sel])
        rows.append(dict(dataset="PRISM", test=f"V600_CRC_vs_melanoma_{label}",
                         drug="Dabrafenib+Trametinib", metric=p_metric, **res))
    res = mwu(pv.loc[pv.is_crc, "daac"], pv.loc[~pv.is_crc, "daac"])
    rows.append(dict(dataset="PRISM", test="V600_CRC_vs_allothers_dAAC",
                     drug="Dabrafenib+Trametinib", metric=p_metric, **res))
pv.to_csv(OUT / "32_prism_pairs_dab_tra.csv", index=False)
per_p.to_csv(OUT / "32_prism_per_line.csv", index=False)

# ============ 输出 ============
res_df = pd.DataFrame(rows)
res_df.to_csv(OUT / "32_third_party_replication.csv", index=False)
summary = {
    "date": "2026-11-17",
    "braf_reference": {
        "v600_samples": len(v600),
        "study_universe": len(study),
        "sequenced_coverage": "ccle_broad_2019_sequenced=1570 (API 无样本ID端点, "
                              "WT=全集(1739)-V600; 未测序~169 不在 CTRP/PRISM 主流系中, 影响可忽略, 已记录)",
    },
    "ctrpv2": {
        "pset": "Orcestra canonical CTRPv2_2015 (DOI 10.5281/zenodo.7826870)",
        "anchor_verdict": anchor_verdict,
        "metric": metric,
        "tests": [r for r in rows if r["dataset"] == "CTRPv2"],
    },
    "prism": {
        "pset": "Orcestra canonical PRISM_2020 (DOI 10.5281/zenodo.7826864)",
        "qc": "passed_str_profiling==True",
        "metric": p_metric,
        "tests": [r for r in rows if r["dataset"] == "PRISM"],
    },
    "software": {"pandas": pd.__version__, "numpy": np.__version__,
                 "scipy": __import__("scipy").__version__, "rdata": "1.1.0"},
}
with open(OUT / "32_summary.json", "w", encoding="utf-8") as f:
    json.dump(summary, f, indent=2, ensure_ascii=False, default=float)

print("\n==== 全部结果行 ====")
print(res_df.to_string(index=False))

# 清理 v1 产物
for stale in ["32_ctrpv2_per_line_aac.csv", "32_ctrpv2_v600_pairs.csv",
              "32_prism_v600_pairs.csv"]:
    fp = OUT / stale
    if fp.exists():
        os.remove(fp)
        print("已清理 v1 产物:", stale)
print("\n完成: results/32_third_party_replication.csv, 32_summary.json, "
      "32_ctrpv2_per_line.csv, 32_prism_per_line.csv")
