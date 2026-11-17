"""诊断: CTRPv2 AAC 分布 + SW1417 组织 + PRISM 药物表结构"""
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import rdata

warnings.filterwarnings("ignore", category=UserWarning)
ROOT = Path(r"D:\AI_work\workbuddy\2026-11-16-16-13-17")

p = rdata.read_rds(str(ROOT / "data/raw/psets/PSet_CTRPv2_2015.rds"))
info = p.sensitivity["info"]
prof = p.sensitivity["profiles"]

# 关注药
drugs = ["Dabrafenib", "Trametinib", "Vemurafenib", "PLX4720", "Selumetinib"]
sub = info[info["drugid"].isin(drugs)][["cellid", "drugid"]]
prof = prof.loc[sub.index]
sub = pd.concat([sub.reset_index(drop=True), prof.reset_index(drop=True)], axis=1)

print("== AAC 分位数(按药) ==")
for d in drugs:
    a = sub.loc[sub["drugid"] == d, "aac_published"].dropna()
    print(f"{d}: n={len(a)} p25={a.quantile(.25):.2f} med={a.median():.2f} "
          f"p75={a.quantile(.75):.2f} p90={a.quantile(.90):.2f} max={a.max():.2f}")

print("\n== 经典系原始 AAC ==")
for cell in ["A-375", "A375", "A2058", "SK-MEL-28", "COLO 205", "HT-29", "RKO",
             "SW1417", "LS411N", "MDST8", "WM266-4", "Malme-3M", "LOX-IMVI"]:
    rows = sub[sub["cellid"].astype(str).str.upper() == cell.upper()]
    if len(rows):
        vals = {r["drugid"]: round(r["aac_published"], 2)
                if pd.notna(r["aac_published"]) else None
                for _, r in rows.iterrows()}
        print(f"{cell}: {vals}")

print("\n== CTRP tissueid 里的肠道/黑色素瘤相关取值 ==")
cell = p.cell
tis = cell[["cellid", "tissueid", "ccle_primary_site", "ccle_primary_hist",
            "ccle_hist_subtype_1"]].astype(str)
interesting = tis[tis.apply(
    lambda r: any(k in " ".join(r.values).upper() for k in
                  ["INTESTINE", "COLON", "COLORECTAL", "MELANOMA", "SKIN"]), axis=1)]
print(interesting["tissueid"].value_counts().to_string())

print("\n== V600 相关 CRC 系的 tissueid ==")
for c in ["SW1417", "HT-29", "RKO", "COLO 205", "COLO 741", "LS411N", "MDST8",
          "SNU-C5", "COLO 201"]:
    m = tis[tis["cellid"] == c]
    if len(m):
        print(f"{c}: tissueid={m.iloc[0]['tissueid']!r} site={m.iloc[0]['ccle_primary_site']!r}")

print("\n== PRISM 药物/细胞表结构 ==")
q = rdata.read_rds(str(ROOT / "data/raw/psets/PSet_PRISM_2020.rds"))
print("drug cols:", list(q.drug.columns))
print("cell cols:", list(q.cell.columns))
print(q.drug.head(3).to_string())
print("\nPRISM drug 表中 RAF/MEK 药:")
namecol = "cpd_name" if "cpd_name" in q.drug.columns else (
    "name" if "name" in q.drug.columns else q.drug.columns[0])
pat = "dabrafenib|trametinib|vemurafenib|PLX-4720|PLX4720|selumetinib|AZD6244|encorafenib|binimetinib|GDC-0973|cobimetinib|refametinib|TAK-733|PD-0325901"
hits = q.drug[q.drug[namecol].astype(str).str.contains(pat, case=False, na=False)]
print(hits[[namecol]].to_string() if len(hits) else "  无")
print("\nPRISM cell 表样例:")
print(q.cell.head(3).to_string())
print("PRISM cellid 样例:", q.sensitivity["info"]["cellid"].unique()[:6].tolist())
print("PRISM tissueid 取值数:", q.cell["tissueid"].nunique() if "tissueid" in q.cell.columns else "n/a")
