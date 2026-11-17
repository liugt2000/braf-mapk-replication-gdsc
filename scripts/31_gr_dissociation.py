"""
31_gr_dissociation.py
---------------------
两件事:

[0] 单位审计 —— scripts/15 与 22 的 stdout 把 exp(LN_IC50) 标成了 "nM",
    但 GDSC 的 LN_IC50 是 ln(IC50, µM)。这里用 22 号结果表复算那些中位数,
    给出正确单位(µM)与 nM 换算, 供稿件修正。锚定验证: A375 trametinib
    exp(-4.87)=0.0076 µM = 7.6 nM, 与文献 2–10 nM 一致 → LN_IC50 必为 ln(µM)。

[1] GR 校正重推导 (Hafner et al. 2016, Nat Methods, PMC5099923)
    质疑: CRC 的 BRAFi/MEKi 解离(Δ=5.28 vs 2.66)会不会只是生长率伪影?
    做法: 从 GDSC2 raw 逐孔数据出发, 每板用 DMSO(NC)对照孔归一化得
    相对活力 Vc/Vctrl, 结合 CMP doubling time (DT) 计算
        GR(c) = 2·(Vc/Vctrl)^(DT/72) − 1
    GR50 ⇔ Vc/Vctrl = 0.75^(72/DT) (比 0.5 的经典阈值更符合"净增殖减半")
    Δ_GR = ln GR50(dabrafenib) − ln GR50(trametinib), CRC vs 其他瘤种。
    附: 用同一归一化管线重算经典 IC50(thr=0.5)与拟合表对照(Spearman),
    验证 raw→归一化这条链本身可信。
"""
import json
import re
import numpy as np
import pandas as pd
from scipy import stats

BASE = r"D:\AI_work\workbuddy\2026-11-16-16-13-17\data\raw"
RES = r"D:\AI_work\workbuddy\2026-11-16-16-13-17\results"
DAB_ID, TRA_ID = 1373, 1372


def norm(s):
    return re.sub(r"[^A-Z0-9]", "", str(s).upper())


# ---------------------------------------------------------------- [0] 单位审计
def unit_audit():
    A = pd.read_csv(rf"{RES}\22_crc_exception.csv")
    g2 = A[(A.campaign == "GDSC2") & A.v600]
    out = {}
    print("=" * 76)
    print("[0] 单位审计: 22 号表(gexp exp(median ln)) 的正确单位是 µM, 不是 nM")
    print("=" * 76)
    for lin in ["COREAD", "SKCM"]:
        s = g2[g2.TCGA_DESC == lin]
        for drug in ["Dabrafenib", "Trametinib"]:
            med_ln = s[drug].median()
            med_um = float(np.exp(med_ln))
            out[f"{lin}_{drug}_median_uM"] = round(med_um, 4)
            out[f"{lin}_{drug}_n"] = int(len(s))
            print(f"  {lin:7s} {drug:12s} n={len(s):2d}  exp(median ln)="
                  f"{med_um:8.4f} µM  (= {med_um*1000:8.2f} nM)")
    ratio = out["COREAD_Dabrafenib_median_uM"] / out["SKCM_Dabrafenib_median_uM"]
    out["dabrafenib_CRC_vs_SKCM_fold"] = round(ratio, 1)
    print(f"  → 达拉非尼 CRC/SKCM 倍数 = {ratio:.1f}× (稿件 30-fold 不变, 仅单位改 µM)")
    a375 = g2[g2.nm == "A375"]
    if not a375.empty:
        print(f"  锚定: A375 trametinib exp(ln)="
              f"{float(np.exp(a375.Trametinib.iloc[0])):.4f} µM "
              f"(文献 2–10 nM → µM 口径成立)")
    return out


# ---------------------------------------------------------------- 数据准备
def load_targets_and_dt():
    A = pd.read_csv(rf"{RES}\22_crc_exception.csv")
    g2 = A[(A.campaign == "GDSC2") & A.v600][
        ["nm", "TCGA_DESC", "Dabrafenib", "Trametinib"]].copy()
    gr = pd.read_csv(rf"{BASE}\growth_rate_20220907.csv")
    gr["nm"] = gr.model_name.map(norm)
    gr["dt"] = pd.to_numeric(gr.doubling_time_hours, errors="coerce")
    dtm = gr.groupby("nm").dt.median()
    g2["DT"] = g2.nm.map(dtm)
    g2["thr"] = 0.75 ** (72.0 / g2.DT)
    n_dup = gr.groupby("nm").size()
    print(f"\n[1a] 目标系 {len(g2)} 株(GDSC2 V600 双药可识别), DT 覆盖 "
          f"{int(g2.DT.notna().sum())}/{len(g2)}; 重复估计的系: "
          f"{sorted(n_dup[n_dup > 1].index.tolist())}")
    return g2, gr


def scan_raw(nms):
    keep = pd.read_csv(
        rf"{BASE}\GDSC2_public_raw_data_24Jul22.csv",
        usecols=["CELL_LINE_NAME", "BARCODE", "TAG", "DRUG_ID",
                 "CONC", "INTENSITY", "DURATION"],
        chunksize=1_000_000, low_memory=False)
    chunks = []
    for i, c in enumerate(keep):
        nmv = c.CELL_LINE_NAME.map(norm)
        drug = pd.to_numeric(c.DRUG_ID, errors="coerce")
        is_dmso = drug.isna() & c.TAG.fillna("").str.startswith("NC")
        is_drug = drug.isin([DAB_ID, TRA_ID])
        sel = nmv.isin(nms) & (is_dmso | is_drug)
        if sel.any():
            d = c[sel].copy()
            d["nm"] = nmv[sel]
            d["is_dmso"] = is_dmso[sel]
            chunks.append(d)
        if (i + 1) % 10 == 0:
            print(f"    chunk {i+1}, 累计保留 {sum(len(x) for x in chunks):,} 孔")
    raw = pd.concat(chunks, ignore_index=True)
    print(f"[1b] raw 过滤完成: {len(raw):,} 孔 "
          f"(DMSO {int(raw.is_dmso.sum()):,} / 药孔 {int((~raw.is_dmso).sum()):,}), "
          f"DURATION={sorted(raw.DURATION.unique().tolist())}")
    return raw


def gr50(logc, rel, thr):
    """沿浓度升序找第一个 rel<=thr 的交叉点, log-conc 线性插值; 返回(µM, censor)"""
    order = np.argsort(logc)
    lc, rv = np.asarray(logc)[order], np.asarray(rel)[order]
    below = rv <= thr
    if not below.any():
        return np.nan, "right"      # 最高浓度仍未压到阈值
    j = int(np.argmax(below))
    if j == 0:
        return 10 ** lc[0], "left"   # 最低浓度已低于阈值
    x0, x1, y0, y1 = lc[j - 1], lc[j], rv[j - 1], rv[j]
    if y1 == y0:
        return 10 ** x1, "ok"
    x = x0 + (thr - y0) * (x1 - x0) / (y1 - y0)
    return 10 ** x, "ok"


def main():
    audit = unit_audit()
    g2, gr = load_targets_and_dt()
    raw = scan_raw(set(g2.nm))

    ctrl = (raw[raw.is_dmso]
            .groupby(["BARCODE"]).INTENSITY.median().rename("ctrl"))
    drug = raw[~raw.is_dmso].merge(ctrl, on="BARCODE", how="left")
    n_no_ctrl = int(drug.ctrl.isna().sum())
    drug = drug[drug.ctrl.notna() & (drug.ctrl > 0)]
    drug["rel"] = drug.INTENSITY / drug.ctrl
    drug["drug"] = np.where(drug.DRUG_ID == DAB_ID, "dab", "tra")
    print(f"[1c] 有对照板的药孔 {len(drug):,} (无对照弃 {n_no_ctrl}); "
          f"板级 DMSO 中位数范围 {drug.ctrl.min():.0f}–{drug.ctrl.max():.0f}")

    med = (drug.groupby(["nm", "drug", "CONC"]).rel.median()
           .reset_index())

    rows, anchor_raw, anchor_fit = [], [], []
    for _, t in g2.iterrows():
        r = dict(nm=t.nm, tissue=t.TCGA_DESC, DT=t.DT,
                 fitted_ln_dab=t.Dabrafenib, fitted_ln_tra=t.Trametinib)
        for key, did, thr in [("dab", DAB_ID, 0.5), ("tra", TRA_ID, 0.5)]:
            s = med[(med.nm == t.nm) & (med.drug == key)]
            if len(s) < 5 or pd.isna(t.DT):
                r[f"n_conc_{key}"] = len(s)
                r[f"IC50raw_{key}"] = np.nan
                r[f"GR50_{key}"] = np.nan
                continue
            lc = np.log10(s.CONC.values)
            r[f"n_conc_{key}"] = len(s)
            r[f"IC50raw_{key}"], _ = gr50(lc, s.rel.values, thr)
            r[f"GR50_{key}"], r[f"censor_{key}"] = gr50(
                lc, s.rel.values, float(t.thr))
            anchor_raw.append(r[f"IC50raw_{key}"])
            anchor_fit.append(r[f"fitted_ln_{key}"])
        rows.append(r)
    R = pd.DataFrame(rows)

    ok = R.dropna(subset=["GR50_dab", "GR50_tra"])
    ok = ok[(ok.censor_dab == "ok") & (ok.censor_tra == "ok")]
    ok["delta_GR"] = np.log(ok.GR50_dab) - np.log(ok.GR50_tra)
    ok["delta_lnIC50"] = ok.fitted_ln_dab - ok.fitted_ln_tra
    crc, oth = ok[ok.tissue == "COREAD"], ok[ok.tissue != "COREAD"]
    p_gr = (stats.mannwhitneyu(crc.delta_GR, oth.delta_GR).pvalue
            if len(crc) >= 3 and len(oth) >= 5 else np.nan)

    # 验证锚: raw 经典 IC50 vs 拟合表
    v = R.dropna(subset=["IC50raw_dab"])
    rho_dab = (stats.spearmanr(np.log(v.IC50raw_dab), v.fitted_ln_dab)
               .statistic if len(v) >= 8 else np.nan)
    v2 = R.dropna(subset=["IC50raw_tra"])
    rho_tra = (stats.spearmanr(np.log(v2.IC50raw_tra), v2.fitted_ln_tra)
               .statistic if len(v2) >= 8 else np.nan)

    print("\n" + "=" * 76)
    print("[2] GR 校正后的解离 (Δ_GR = ln GR50_dab − ln GR50_tra)")
    print("=" * 76)
    print(f"  双药 GR50 均可识别: {len(ok)}/{len(g2)} 株 "
          f"(CRC {len(crc)}, 非CRC {len(oth)})")
    for lab, s in [("COREAD", crc), ("非COREAD", oth)]:
        if len(s):
            print(f"  {lab:8s} n={len(s):2d}  Δ_GR 中位={s.delta_GR.median():6.2f}"
                  f"  GR50 dab 中位={s.GR50_dab.median():.3f} µM  "
                  f"GR50 tra 中位={s.GR50_tra.median():.3f} µM")
    print(f"  Mann-Whitney p = {p_gr:.2e}" if pd.notna(p_gr) else "  p=n/a")
    print(f"  [同批线 Δ_lnIC50(未校正) 中位: CRC {crc.delta_lnIC50.median():.2f}"
          f" vs 非CRC {oth.delta_lnIC50.median():.2f}]")

    R.to_csv(rf"{RES}\31_gr_dissociation.csv", index=False)
    summary = dict(
        unit_audit=audit,
        n_lines=int(len(g2)), n_dt=int(g2.DT.notna().sum()),
        n_gr_identifiable=int(len(ok)),
        crc_n=int(len(crc)), oth_n=int(len(oth)),
        delta_GR_crc=round(float(crc.delta_GR.median()), 3) if len(crc) else None,
        delta_GR_oth=round(float(oth.delta_GR.median()), 3) if len(oth) else None,
        delta_GR_p=float(p_gr) if pd.notna(p_gr) else None,
        gr50_dab_crc_uM=round(float(crc.GR50_dab.median()), 4) if len(crc) else None,
        gr50_dab_oth_uM=round(float(oth.GR50_dab.median()), 4) if len(oth) else None,
        gr50_tra_crc_uM=round(float(crc.GR50_tra.median()), 4) if len(crc) else None,
        gr50_tra_oth_uM=round(float(oth.GR50_tra.median()), 4) if len(oth) else None,
        rawIC50_vs_fitted_spearman_dab=round(float(rho_dab), 3) if pd.notna(rho_dab) else None,
        rawIC50_vs_fitted_spearman_tra=round(float(rho_tra), 3) if pd.notna(rho_tra) else None,
    )
    with open(rf"{RES}\31_gr_summary.json", "w") as f:
        json.dump(summary, f, indent=2)
    print(f"\n写入 {RES}\\31_gr_dissociation.csv 与 31_gr_summary.json")
    print(ok[["nm", "tissue", "DT", "GR50_dab", "GR50_tra", "delta_GR",
              "delta_lnIC50"]].to_string(index=False))


if __name__ == "__main__":
    main()
