#!/usr/bin/env python
"""
preflight_audit.py — 药物反应标签的预检审计
=================================================
在任何建模之前，先回答三个问题：
  1. 我的目标量有多少是可测量的（而不是外推的）？
  2. 这个目标量在两次独立测定之间能复现吗？
  3. 单浓度下，细胞系之间的差异是否大于复孔差异？

用法
----
  python preflight_audit.py --fitted GDSC2_fitted_dose_response.csv \
                            [--fitted2 GDSC1_fitted_dose_response.csv] \
                            [--disease COREAD] \
                            [--drugs "5-Fluorouracil,Oxaliplatin"] \
                            [--out audit_report]

输入要求（GDSC 格式的拟合表，最少需要这些列）：
  CELL_LINE_NAME, DRUG_NAME, TCGA_DESC, MIN_CONC, MAX_CONC, LN_IC50
可选：AUC

输出
----
  <out>.json   机器可读
  <out>.txt    人可读报告
  <out>_by_drug.csv

判定阈值（可在代码顶部修改）
--------------------------
  censoring_fraction  > 0.30   → 目标量不可作为连续回归目标
  cross_campaign_r    < 0.70   → 跨测定不可复现，单队列 CV 会高估
  identifiable_pairs  < 0.50   → 大部分配对无法用于信度估计
"""
from __future__ import annotations

import argparse
import json
import sys

import numpy as np
import pandas as pd

THRESHOLDS = {
    "censoring_fraction": 0.30,
    "cross_campaign_r": 0.70,
    "identifiable_pairs": 0.50,
}


def load(path: str, disease: str | None) -> pd.DataFrame:
    df = pd.read_csv(path)
    need = {"CELL_LINE_NAME", "DRUG_NAME", "MIN_CONC", "MAX_CONC", "LN_IC50"}
    missing = need - set(df.columns)
    if missing:
        sys.exit(f"[错误] {path} 缺少必需列: {sorted(missing)}")
    if disease and "TCGA_DESC" in df.columns:
        df = df[df.TCGA_DESC == disease].copy()
        if df.empty:
            sys.exit(f"[错误] TCGA_DESC=='{disease}' 无数据；可选值: "
                     f"{sorted(pd.read_csv(path, usecols=['TCGA_DESC']).TCGA_DESC.dropna().unique())}")
    df["ic50"] = np.exp(df["LN_IC50"])
    df["identifiable"] = (df.ic50 >= df.MIN_CONC) & (df.ic50 <= df.MAX_CONC)
    return df


def censoring(df: pd.DataFrame) -> dict:
    hi = float((df.ic50 > df.MAX_CONC).mean())
    lo = float((df.ic50 < df.MIN_CONC).mean())
    return {"n_curves": int(len(df)),
            "n_lines": int(df.CELL_LINE_NAME.nunique()),
            "n_drugs": int(df.DRUG_NAME.nunique()),
            "frac_extrapolated_high": round(hi, 4),
            "frac_extrapolated_low": round(lo, 4),
            "frac_identifiable": round(1 - hi - lo, 4)}


def per_drug(df: pd.DataFrame) -> pd.DataFrame:
    g = df.groupby("DRUG_NAME")
    out = g.apply(lambda s: pd.Series({
        "n_curves": len(s),
        "n_lines": s.CELL_LINE_NAME.nunique(),
        "frac_extrapolated": float((~s.identifiable).mean()),
        "ln_ic50_sd": float(s.LN_IC50.std()),
        "fold_span": float(np.exp(s.LN_IC50.max() - s.LN_IC50.min())),
    }), include_groups=False).reset_index()
    return out.sort_values("frac_extrapolated", ascending=False)


def reliability(a: pd.DataFrame, b: pd.DataFrame) -> dict:
    ka = a.groupby(["CELL_LINE_NAME", "DRUG_NAME"], as_index=False).agg(
        v=("LN_IC50", "mean"), ida=("identifiable", "all"))
    kb = b.groupby(["CELL_LINE_NAME", "DRUG_NAME"], as_index=False).agg(
        v2=("LN_IC50", "mean"), idb=("identifiable", "all"))
    m = ka.merge(kb, on=["CELL_LINE_NAME", "DRUG_NAME"])
    if len(m) < 10:
        return {"error": "共同配对不足 10 个，无法估计信度", "n_pairs": int(len(m))}
    x, y = m.v.values, m.v2.values
    d = np.abs(y - x)
    n = len(m)
    grand = np.mean(np.concatenate([x, y]))
    ms_b = 2 * (((x + y) / 2 - grand) ** 2).sum() / (n - 1)
    ms_w = ((x - y) ** 2).sum() / n
    both = m.ida & m.idb
    out = {
        "n_pairs": int(n),
        "pearson_r": round(float(np.corrcoef(x, y)[0, 1]), 4),
        "icc_1_1": round(float((ms_b - ms_w) / (ms_b + ms_w)), 4),
        "median_fold_diff": round(float(np.exp(np.median(d))), 3),
        "frac_pairs_gt_2fold": round(float((d > np.log(2)).mean()), 4),
        "frac_both_identifiable": round(float(both.mean()), 4),
    }
    if both.sum() >= 8:
        out["pearson_r_both_identifiable"] = round(
            float(np.corrcoef(m.v[both], m.v2[both])[0, 1]), 4)
        out["n_both_identifiable"] = int(both.sum())
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description="药物反应标签预检审计")
    ap.add_argument("--fitted", required=True, help="拟合表 CSV（如 GDSC2）")
    ap.add_argument("--fitted2", help="第二个独立测定的拟合表（如 GDSC1）")
    ap.add_argument("--disease", default="COREAD", help="TCGA_DESC 取值，默认 COREAD")
    ap.add_argument("--drugs", default="", help="逗号分隔，仅关注这些药")
    ap.add_argument("--out", default="audit_report", help="输出前缀")
    a = ap.parse_args()

    for f in (a.fitted, a.fitted2):
        if f and not __import__("os").path.exists(f):
            sys.exit(f"[错误] 文件不存在: {f}")

    d1 = load(a.fitted, a.disease)
    cen = censoring(d1)
    pdrug = per_drug(d1)

    if a.drugs:
        keep = [s.strip() for s in a.drugs.split(",") if s.strip()]
        pdrug = pdrug[pdrug.DRUG_NAME.isin(keep)]

    rel = None
    if a.fitted2:
        d2 = load(a.fitted2, a.disease)
        rel = reliability(d1, d2)

    # ---- 判定 ----
    verdict = []
    if cen["frac_identifiable"] < (1 - THRESHOLDS["censoring_fraction"]):
        verdict.append(
            f"✗ 删失过重：仅 {cen['frac_identifiable']:.1%} 的 IC50 落在实测浓度区间内"
            f"（阈值 ≥ {1-THRESHOLDS['censoring_fraction']:.0%}）。"
            "IC50 不可作为连续回归目标；考虑改用删失/区间回归，或改用可测量的替代终点。")
    else:
        verdict.append(f"✓ 可识别比例 {cen['frac_identifiable']:.1%} 达到阈值。")

    if rel and "pearson_r" in rel:
        if rel["pearson_r"] < THRESHOLDS["cross_campaign_r"]:
            verdict.append(
                f"✗ 跨测定信度不足：r = {rel['pearson_r']:.2f}"
                f"（阈值 ≥ {THRESHOLDS['cross_campaign_r']}），"
                f"中位差 {rel['median_fold_diff']:.1f} 倍。"
                "单队列交叉验证会系统性高估性能。")
        else:
            verdict.append(f"✓ 跨测定信度 r = {rel['pearson_r']:.2f} 达到阈值。")
        if rel.get("frac_both_identifiable", 1) < THRESHOLDS["identifiable_pairs"]:
            verdict.append(
                f"✗ 仅 {rel['frac_both_identifiable']:.1%} 的配对在两次测定中都可识别，"
                "表观一致性可能来自共同外推。")

    report = {"file": a.fitted, "second_file": a.fitted2, "disease": a.disease,
              "censoring": cen, "reliability": rel, "verdict": verdict,
              "thresholds": THRESHOLDS}

    with open(f"{a.out}.json", "w", encoding="utf-8") as f:
        json.dump(report, f, indent=2, ensure_ascii=False)
    pdrug.to_csv(f"{a.out}_by_drug.csv", index=False)

    lines = ["=" * 72,
             "药物反应标签 · 预检审计报告",
             "=" * 72,
             f"数据: {a.fitted}" + (f"  vs  {a.fitted2}" if a.fitted2 else ""),
             f"病种: {a.disease}",
             "",
             f"[1] 规模      {cen['n_curves']:,} 条曲线 / {cen['n_lines']} 株系 / {cen['n_drugs']} 个药",
             f"[2] 删失      外推(高) {cen['frac_extrapolated_high']:.1%}"
             f"  外推(低) {cen['frac_extrapolated_low']:.1%}"
             f"  可识别 {cen['frac_identifiable']:.1%}"]
    if rel and "pearson_r" in rel:
        lines += ["",
                  f"[3] 信度      配对 {rel['n_pairs']:,}",
                  f"              Pearson r   {rel['pearson_r']:.3f}",
                  f"              ICC(1,1)    {rel['icc_1_1']:.3f}",
                  f"              中位倍数差  {rel['median_fold_diff']:.2f}x",
                  f"              >2倍的配对   {rel['frac_pairs_gt_2fold']:.1%}",
                  f"              两活动中都可识别 {rel['frac_both_identifiable']:.1%}"]
        if "pearson_r_both_identifiable" in rel:
            lines.append(f"              仅可识别配对 r  {rel['pearson_r_both_identifiable']:.3f}"
                         f" (n={rel['n_both_identifiable']})")
    lines += ["", "逐药结果（按删失比例降序）:", pdrug.head(25).to_string(index=False), "",
              "-" * 72, "判定:", *[f"  {v}" for v in verdict], "=" * 72]
    txt = "\n".join(lines)
    with open(f"{a.out}.txt", "w", encoding="utf-8") as f:
        f.write(txt)
    print(txt)


if __name__ == "__main__":
    main()
