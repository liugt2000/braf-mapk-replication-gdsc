"""
33_figures_jtm.py — JTM/BMC 规格投稿图 (任务 4)

与 24_figures_v2.py 的差异 (JTM 规格清单):
  - Arial, 全部文字 ≥7pt (BMC 图规范), mathtext 同步 Arial
  - 色盲安全 Okabe-Ito 配色 (蓝 #0072B2 / 朱红 #D55E00 / 橙 #E69F00 / 灰 #999999)
  - 面板字母大写加粗
  - 图 1C: 26 个可检验基因 (results/30_burden_enrichment.csv 重算口径), 弃用 20 基因硬编码;
           注记 25/26 (96%)、ERCC5 +5.99、BRAF +0.20 (中位 +1.40 的 1/7)、NRAS −0.04 唯一非正向
  - 图 1A: 48 为已核验值, 保留
  - 图 1D: 抖动改为固定种子 (20261116)
  - 图 2A: ylabel "IC50, µM" (GDSC LN_IC50 = ln µM; v2 的 nM 是错的), med 标注即 µM 值
  - 图 2B: CMS 数字 28/0/7/2/10 与 10/67/47/58/58 已核验, 保留
  - 输出 PNG 400dpi + PDF 矢量双格式, 文件名 *_jtm 以保留 v2 原图作审计
"""
import json

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

plt.rcParams.update({
    "font.family": "Arial",
    "mathtext.fontset": "custom",
    "mathtext.rm": "Arial",
    "mathtext.it": "Arial:italic",
    "mathtext.bf": "Arial:bold",
    "font.size": 7,
    "axes.linewidth": 0.6,
    "xtick.major.width": 0.6,
    "ytick.major.width": 0.6,
    "xtick.labelsize": 7,
    "ytick.labelsize": 7,
    "axes.labelsize": 7.5,
    "axes.titlesize": 8,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "pdf.fonttype": 42,
    "legend.fontsize": 7,
})
BLUE, VERM, ORANGE, GREY, BLACK = ("#0072B2", "#D55E00", "#E69F00",
                                   "#999999", "#000000")
SEED = 20261116


def fig1():
    summ = json.load(open("results/18_summary.json"))
    null = pd.read_csv("results/18_null_distribution.csv")

    fig, axes = plt.subplots(
        2, 2, figsize=(7.2, 6.6),
        gridspec_kw={"height_ratios": [1, 1.22]})
    rng = np.random.default_rng(SEED)

    # ---- A 漏斗 ----
    ax = axes[0, 0]
    stages = ["Testable\npairs", "FDR q<0.05\n(one campaign)",
              "Burden-\nadjusted", "Replicated\n(both campaigns)"]
    vals = [summ["n_testable"], summ["fdr_hits"], 48, summ["n_pass"]]
    cols = [GREY, GREY, GREY, BLUE]
    y = np.arange(len(stages))[::-1]
    ax.barh(y, vals, color=cols, height=0.55, edgecolor="none")
    for yy, v in zip(y, vals):
        ax.text(v * 1.35, yy, f"{v:,}", va="center", fontsize=7.5,
                weight="bold" if v == summ["n_pass"] else "normal")
    ax.set_yticks(y)
    ax.set_yticklabels(stages, fontsize=7)
    ax.set_xscale("log")
    ax.set_xlim(1, 4000)
    ax.set_xlabel("number of gene–drug pairs")
    ax.set_title("A  Three filters converge on 5 pairs", loc="left",
                 weight="bold")
    ax.text(1.2, -0.72, "all 5 are BRAF → BRAF/MEK inhibitor",
            fontsize=7, color=BLUE, va="center")
    ax.set_ylim(-1.0, len(stages) - 0.4)

    # ---- B 置换零分布 ----
    ax = axes[0, 1]
    vc = null.n_pass.value_counts().sort_index()
    ax.bar(vc.index, vc.values, color=GREY, width=0.6)
    ax.axvline(summ["n_pass"], color=BLUE, lw=1.2)
    ax.annotate("observed = 5\n($p$ = 0.005)",
                xy=(summ["n_pass"], vc.max() * 0.55),
                xytext=(summ["n_pass"] * 0.24, vc.max() * 0.78),
                fontsize=7, color=BLUE,
                arrowprops=dict(arrowstyle="->", color=BLUE, lw=0.7))
    ax.set_xlim(-0.6, 6)
    ax.set_xlabel("pairs passing the filter\nunder permuted gene labels")
    ax.set_ylabel("permutations (n = 200)")
    ax.set_title("B  Not explainable by multiple testing", loc="left",
                 weight="bold")

    # ---- C 负荷富集 (26 可检验基因, 重算口径) ----
    ax = axes[1, 0]
    b = pd.read_csv("results/30_burden_enrichment.csv")
    b = b[b.testable].dropna(subset=["d_burden"])
    b = b.sort_values("d_burden", ascending=True).reset_index(drop=True)
    colors = []
    for g in b.gene:
        colors.append(BLUE if g == "BRAF" else
                      (VERM if g == "NRAS" else GREY))
    ax.barh(np.arange(len(b)), b.d_burden, color=colors, height=0.66,
            edgecolor="none")
    ax.axvline(0, color="black", lw=0.6)
    med = float(b.d_burden.median())
    ax.axvline(med, color=BLUE, lw=0.7, ls=(0, (4, 3)), alpha=0.7)
    ax.set_yticks(np.arange(len(b)))
    ax.set_yticklabels(b.gene, fontsize=7)
    ax.set_xlabel("enrichment of high-burden lines in mutation carriers\n"
                  "(Cohen's $d$; testable: $\\geq$8 altered, $\\geq$20 wild-type)")
    ax.set_title("C  25 of 26 testable genes track burden;\n"
                 "     BRAF barely does, NRAS does not", loc="left",
                 weight="bold")
    ax.text(med * 1.55, 3.2,
            f"25/26 positive (96%)\nmedian $d$ = +{med:.2f}\n"
            "ERCC5 +5.99 (largest)\n"
            "BRAF +0.20 (1/7 of median)\nNRAS −0.04 (only negative)",
            fontsize=7, va="center", ha="left",
            bbox=dict(facecolor="white", edgecolor="none", alpha=0.0))
    ax.set_xlim(min(b.d_burden.min() - 0.15, -0.2),
                b.d_burden.max() * 1.06)

    # ---- D 解离 ----
    ax = axes[1, 1]
    A = pd.read_csv("results/22_crc_exception.csv")
    g2 = A[A.campaign == "GDSC2"]
    groups = ["COREAD", "THCA", "SKCM"]
    labs = ["Colorectal", "Thyroid", "Melanoma"]
    data = [g2[(g2.TCGA_DESC == g) & g2.v600].delta_brafi_vs_meki.dropna().values
            for g in groups]
    pos = [0, 1, 2]
    bp = ax.boxplot(data, positions=pos, widths=0.5, patch_artist=True,
                    showfliers=False,
                    medianprops=dict(color="black", lw=1.0))
    for i, patch in enumerate(bp["boxes"]):
        patch.set_facecolor(BLUE if i == 0 else "#dfe6ee")
        patch.set_edgecolor("black")
        patch.set_linewidth(0.6)
    for i, arr in enumerate(data):
        ax.scatter(rng.normal(pos[i], 0.055, len(arr)), arr, s=10,
                   color="black", alpha=0.55, linewidth=0, zorder=3)
    ax.set_xticks(pos)
    ax.set_xticklabels([f"{l}\n(n={len(d_)})" for l, d_ in zip(labs, data)],
                       fontsize=7)
    ax.set_ylabel("Δ = ln IC$_{50}$ dabrafenib − ln IC$_{50}$ trametinib")
    ax.set_title("D  Colorectal dissociation of BRAFi vs MEKi",
                 loc="left", weight="bold")
    ax.text(0.03, 0.99, "$p$ = 1.0×10$^{-4}$ (colorectal vs others)\n"
            "replicated in GDSC1: $p$ = 0.005\n"
            "cross-campaign $r$ = 0.73, sign 100%",
            transform=ax.transAxes, fontsize=7, va="top", ha="left")
    ax.set_ylim(min(d_.min() for d_ in data) - 0.8,
                max(d_.max() for d_ in data) + 3.1)

    plt.tight_layout()
    for ext in ("png", "pdf"):
        plt.savefig(f"results/fig1_evidence_chain_jtm.{ext}", dpi=400,
                    bbox_inches="tight")
    plt.close()
    print("写入 results/fig1_evidence_chain_jtm.png / .pdf")


def fig2():
    fig, axes = plt.subplots(1, 3, figsize=(7.2, 2.6),
                             gridspec_kw={"width_ratios": [1.05, 1, 1]})
    # ---- A 绝对效力 (µM!) ----
    ax = axes[0]
    C = pd.read_csv("results/22_crc_exception.csv")
    C = C[C.campaign == "GDSC2"]
    for i, (drug, col) in enumerate([("Dabrafenib", BLUE), ("Trametinib", ORANGE)]):
        for j, lin in enumerate(["COREAD", "SKCM"]):
            s = C[(C.TCGA_DESC == lin) & C.v600][drug].dropna()
            med = np.exp(s.median())
            lo, hi = np.exp(s.quantile(.25)), np.exp(s.quantile(.75))
            x = j * 2.6 + i * 0.75
            ax.errorbar(x, med, yerr=[[med - lo], [hi - med]], fmt="o",
                        color=col, ms=5, capsize=2.5, lw=1.0)
            ax.text(x + 0.16, med, f"{med:.2f}", ha="left", va="center",
                    fontsize=7, color=col, clip_on=False)
    ax.set_yscale("log")
    ax.set_xticks([0.38, 1.12, 2.98, 3.72])
    ax.set_xticklabels(["Colo", "Mela", "Colo", "Mela"], fontsize=7)
    ax.set_ylabel("IC$_{50}$, µM")
    ax.set_ylim(0.0025, 14)
    ax.set_xlim(-0.55, 4.55)
    ax.set_title("A  BRAF V600: MEKi equipotent,\n     BRAFi 30-fold weaker in CRC",
                 loc="left", fontsize=8, weight="bold")
    ax.legend(handles=[Line2D([], [], color=BLUE, marker="o", ls="", ms=4,
                              label="dabrafenib"),
                       Line2D([], [], color=ORANGE, marker="o", ls="", ms=4,
                              label="trametinib")],
              frameon=False, loc="upper right")

    # ---- B CMS ----
    ax = axes[1]
    cms = {"CMS1": 28, "CMS2": 0, "CMS3": 7, "CMS4": 2, "mixed": 10}
    cms_wt = {"CMS1": 10, "CMS2": 67, "CMS3": 47, "CMS4": 58, "mixed": 58}
    cats = ["CMS1", "CMS2", "CMS3", "CMS4", "mixed"]
    x = np.arange(len(cats))
    v1 = np.array([cms[c] / 47 * 100 for c in cats])
    v2 = np.array([cms_wt[c] / 240 * 100 for c in cats])
    ax.bar(x - 0.19, v1, width=0.36, color=BLUE, label="BRAF V600")
    ax.bar(x + 0.19, v2, width=0.36, color="#c9c9c9", label="BRAF WT")
    ax.set_xticks(x)
    ax.set_xticklabels(cats, fontsize=7, rotation=30, ha="right")
    ax.set_ylabel("% of patients")
    ax.set_title("B  BRAF V600 concentrates\n     in CMS1 (n=47 vs 240)",
                 loc="left", fontsize=8, weight="bold")
    ax.legend(frameon=False, loc="upper right")

    # ---- C OS ----
    ax = axes[2]
    df = pd.read_csv("results/23_patient_braf_analysis.csv")
    for flag, col, lab in [(0, GREY, "BRAF WT"), (1, BLUE, "BRAF V600")]:
        t = df[df.braf == flag]
        ts = np.sort(t.os_m.values)
        ev = t.event.values[np.argsort(t.os_m.values)]
        n = len(ts)
        at_risk, s = n, 1.0
        xs, ys = [0], [1.0]
        for i, ti in enumerate(ts):
            if ev[i] == 1:
                s *= (at_risk - 1) / at_risk if at_risk > 1 else 0
                xs.append(ti); ys.append(s)
            at_risk -= 1
        ax.step(xs, ys, where="post", color=col, lw=1.1,
                label=f"{lab} (n={n})")
    ax.set_xlabel("months")
    ax.set_ylabel("overall survival")
    ax.set_xlim(0, 120)
    ax.set_ylim(0, 1.02)
    ax.set_title("C  Overall survival, stage I–III\n     (adjusted HR 0.96, $p$=0.86)",
                 loc="left", fontsize=8, weight="bold")
    ax.legend(frameon=False, loc="lower left")

    plt.tight_layout()
    for ext in ("png", "pdf"):
        plt.savefig(f"results/fig2_translational_jtm.{ext}", dpi=400,
                    bbox_inches="tight")
    plt.close()
    print("写入 results/fig2_translational_jtm.png / .pdf")


if __name__ == "__main__":
    fig1()
    fig2()
