"""
24_figures_v2.py
----------------
Fig 2 (方法学主图, 4 面板):
  A 复制漏斗: 1338 可检验 -> 61 FDR -> 48 负荷校正 -> 5 跨活动复制(全部 BRAF)
  B 置换零分布 vs 观测值 5
  C 突变负荷混杂: 每个驱动基因"突变组在本系中更富高负荷"的效应量
  D BRAFi/MEKi 解离度按瘤种(两活动)

Fig 3 (转化主图, 3 面板):
  E BRAF V600 与野生型在达拉非尼/曲美替尼上的绝对 IC50
  F 病人队列 CMS 分布
  G 队列 OS 分层
"""
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 7,
    "axes.linewidth": 0.6, "xtick.major.width": 0.6,
    "ytick.major.width": 0.6, "axes.spines.top": False,
    "axes.spines.right": False, "pdf.fonttype": 42,
})
BLUE, RED, GREY, GOLD = "#2b6cb0", "#c53030", "#8a8a8a", "#b7791f"


def fig2():
    summ = json.load(open("results/18_summary.json"))
    null = pd.read_csv("results/18_null_distribution.csv")

    fig, axes = plt.subplots(2, 2, figsize=(7.2, 5.4))
    ax = axes[0, 0]

    # ---- A 漏斗 ----
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
    ax.set_yticklabels(stages, fontsize=6.5)
    ax.set_xscale("log")
    ax.set_xlim(1, 4000)
    ax.set_xlabel("number of gene–drug pairs")
    ax.set_title("A  Three filters converge on 5 pairs", loc="left",
                 fontsize=7.8, weight="bold")
    ax.text(1.2, -0.62, "all 5 are BRAF → BRAF/MEK inhibitor",
            fontsize=6.4, color=BLUE, va="center")
    ax.set_ylim(-1.0, len(stages) - 0.4)

    # ---- B 置换零 ----
    ax = axes[0, 1]
    vc = null.n_pass.value_counts().sort_index()
    ax.bar(vc.index, vc.values, color=GREY, width=0.6)
    ax.axvline(summ["n_pass"], color=BLUE, lw=1.2)
    ax.annotate(f"observed = {summ['n_pass']}\n(empirical p = {summ['empirical_p']:.3f})",
                xy=(summ["n_pass"], vc.max() * 0.55),
                xytext=(summ["n_pass"] * 0.28, vc.max() * 0.75),
                fontsize=6.4, color=BLUE,
                arrowprops=dict(arrowstyle="->", color=BLUE, lw=0.7))
    ax.set_xlim(-0.6, 6)
    ax.set_xlabel("pairs passing the filter under permuted gene labels")
    ax.set_ylabel("permutations (n=200)")
    ax.set_title("B  Not explainable by multiple testing", loc="left",
                 fontsize=7.8, weight="bold")

    # ---- C 突变负荷混杂 ----
    ax = axes[1, 0]
    b = pd.read_csv("results/19_burden_adjusted_scan.csv")
    # 用 19 号脚本的口径重算每基因的负荷效应
    import re
    from scipy import stats as st
    mm = pd.read_csv("data/raw/mutations_latest.csv.gz")
    pm = mm.protein_mutation.fillna("").astype(str)
    mm["dmg"] = (mm.cancer_driver.astype(str).str.lower() == "true") | \
                pm.str.contains(r"\*|fs", regex=True)
    mm = mm[mm.dmg]
    burden = (mm.drop_duplicates(["model_name", "gene_symbol"])
              .groupby("model_name").size())
    burden.index = [re.sub(r"[^A-Z0-9]", "", str(x).upper()) for x in burden.index]
    burden = burden[~burden.index.duplicated()]
    d = pd.read_csv("data/raw/GDSC2_fitted_dose_response_24Jul22.csv")
    d["ic50"] = np.exp(d.LN_IC50)
    d = d[(d.ic50 >= d.MIN_CONC) & (d.ic50 <= d.MAX_CONC)].copy()
    d["nm"] = d.CELL_LINE_NAME.map(lambda x: re.sub(r"[^A-Z0-9]", "", str(x).upper()))
    P = d.pivot_table(index="nm", columns="DRUG_NAME", values="LN_IC50", aggfunc="mean")
    keep = P.notna().sum(axis=0) >= 40
    P = P[keep.index[keep]]
    mr = P.rank(axis=0).mean(axis=1)
    genes = ["KRAS", "BRAF", "TP53", "APC", "PIK3CA", "SMAD4", "ERBB2",
             "MLH1", "MSH6", "BRCA2", "ARID1A", "FBXW7", "ATM", "RNF43",
             "PTEN", "NF1", "CTNNB1", "CDKN2A", "RB1", "NRAS"]
    vals = []
    for g in genes:
        s = {re.sub(r"[^A-Z0-9]", "", str(x).upper())
             for x in mm.loc[mm.gene_symbol == g, "model_name"]}
        isalt = np.array([k in s for k in mr.index])
        if isalt.sum() < 8 or (~isalt).sum() < 20:
            continue
        bb = burden.reindex(mr.index).values
        ok = ~np.isnan(bb)
        a, w = bb[ok & isalt], bb[ok & ~isalt]
        sp = np.sqrt(((len(a) - 1) * a.var(ddof=1) + (len(w) - 1) * w.var(ddof=1))
                     / (len(a) + len(w) - 2))
        vals.append((g, (a.mean() - w.mean()) / sp))
    V = pd.DataFrame(vals, columns=["gene", "d"]).sort_values("d")
    cols = [BLUE if g == "BRAF" else GREY for g in V.gene]
    ax.barh(np.arange(len(V)), V.d, color=cols, height=0.62)
    ax.axvline(0, color="black", lw=0.6)
    ax.set_yticks(np.arange(len(V)))
    ax.set_yticklabels(V.gene, fontsize=5.6)
    ax.set_xlabel("enrichment of high-burden lines\nin mutation-carriers (Cohen's d)")
    ax.set_title("C  Most driver genes track mutational burden;\n    BRAF does not",
                 loc="left", fontsize=7.8, weight="bold")
    bi = list(V.gene).index("BRAF")
    ax.text(V.d.iloc[bi] + 0.10, bi, "BRAF", fontsize=6.4, color=BLUE,
            va="center", ha="left", weight="bold")

    # ---- D BRAFi/MEKi 解离度 ----
    ax = axes[1, 1]
    A = pd.read_csv("results/22_crc_exception.csv")
    g2 = A[A.campaign == "GDSC2"]
    groups = ["COREAD", "THCA", "SKCM"]
    labs = ["Colorectal", "Thyroid", "Melanoma"]
    data = [g2[(g2.TCGA_DESC == g) & g2.v600].delta_brafi_vs_meki.dropna().values
            for g in groups]
    pos = [0, 1, 2]
    bp = ax.boxplot(data, positions=pos, widths=0.5, patch_artist=True,
                    showfliers=False, medianprops=dict(color="black", lw=1.0))
    for i, patch in enumerate(bp["boxes"]):
        patch.set_facecolor(BLUE if i == 0 else "#dfe6ee")
        patch.set_edgecolor("black")
        patch.set_linewidth(0.6)
    for i, arr in enumerate(data):
        ax.scatter(np.random.normal(pos[i], 0.055, len(arr)), arr, s=9,
                   color="black", alpha=0.55, linewidth=0, zorder=3)
    ax.set_xticks(pos)
    ax.set_xticklabels([f"{l}\n(n={len(d_)})" for l, d_ in zip(labs, data)],
                       fontsize=6.4)
    ax.set_ylabel("ln IC50 dabrafenib\n− ln IC50 trametinib")
    ax.set_title("D  Colorectal dissociation of BRAFi vs MEKi", loc="left",
                 fontsize=7.8, weight="bold")
    ax.text(0.03, 0.99, "p = 1.0×10⁻⁴ (colorectal vs others)\n"
            "replicated in GDSC1: p = 0.005\n"
            "cross-campaign r = 0.73, sign 100%",
            transform=ax.transAxes, fontsize=5.8, va="top", ha="left")
    ax.set_ylim(min(d_.min() for d_ in data) - 0.8,
                max(d_.max() for d_ in data) + 2.6)

    plt.tight_layout()
    plt.savefig("results/fig1_evidence_chain.png", dpi=400, bbox_inches="tight")
    plt.close()
    print("写入 results/fig1_evidence_chain.png")


def fig3():
    fig, axes = plt.subplots(1, 3, figsize=(7.2, 2.5),
                             gridspec_kw={"width_ratios": [1.05, 1, 1]})
    # ---- E 绝对 IC50 ----
    ax = axes[0]
    C = pd.read_csv("results/22_crc_exception.csv")
    C = C[C.campaign == "GDSC2"]
    for i, (drug, col) in enumerate([("Dabrafenib", BLUE), ("Trametinib", GOLD)]):
        for j, lin in enumerate(["COREAD", "SKCM"]):
            s = C[(C.TCGA_DESC == lin) & C.v600][drug].dropna()
            med = np.exp(s.median())
            lo, hi = np.exp(s.quantile(.25)), np.exp(s.quantile(.75))
            x = j * 2.6 + i * 0.75
            ax.errorbar(x, med, yerr=[[med - lo], [hi - med]], fmt="o",
                        color=col, ms=5, capsize=2.5, lw=1.0)
            ax.text(x + 0.16, med, f"{med:.2f}", ha="left", va="center",
                    fontsize=5.6, color=col, clip_on=False)
    ax.set_yscale("log")
    ax.set_xticks([0.38, 1.12, 2.98, 3.72])
    ax.set_xticklabels(["Colo", "Mela", "Colo", "Mela"], fontsize=6)
    ax.set_ylabel("IC$_{50}$, nM")
    ax.set_ylim(0.0025, 14)
    ax.set_xlim(-0.55, 4.55)
    ax.set_title("A  BRAF V600: MEKi equipotent,\n    BRAFi 30-fold weaker in CRC",
                 loc="left", fontsize=7.4, weight="bold")
    ax.legend(handles=[Line2D([], [], color=BLUE, marker="o", ls="", ms=4,
                              label="dabrafenib"),
                       Line2D([], [], color=GOLD, marker="o", ls="", ms=4,
                              label="trametinib")],
              frameon=False, fontsize=5.8, loc="upper right")

    # ---- F CMS ----
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
    ax.set_xticklabels(cats, fontsize=5.8, rotation=30, ha="right")
    ax.set_ylabel("% of patients")
    ax.set_title("B  BRAF V600 concentrates\n    in CMS1 (n=47 vs 240)",
                 loc="left", fontsize=7.4, weight="bold")
    ax.legend(frameon=False, fontsize=5.6, loc="upper right")

    # ---- G OS ----
    ax = axes[2]
    df = pd.read_csv("results/23_patient_braf_analysis.csv")
    for flag, col, lab in [(0, GREY, "BRAF WT"), (1, BLUE, "BRAF V600")]:
        t = df[df.braf == flag]
        # KM by hand
        ts = np.sort(t.os_m.values)
        ev = t.event.values[np.argsort(t.os_m.values)]
        n = len(ts)
        surv, at_risk, s = [], n, 1.0
        xs, ys = [0], [1.0]
        for i, ti in enumerate(ts):
            if ev[i] == 1:
                s *= (at_risk - 1) / at_risk if at_risk > 1 else 0
                xs.append(ti); ys.append(s)
            at_risk -= 1
        ax.step(xs, ys, where="post", color=col, lw=1.1, label=f"{lab} (n={n})")
    ax.set_xlabel("months")
    ax.set_ylabel("overall survival")
    ax.set_xlim(0, 120)
    ax.set_ylim(0, 1.02)
    ax.set_title("C  OS, stage I–III\n    (adjusted HR 0.96, p=0.86)",
                 loc="left", fontsize=7.4, weight="bold")
    ax.legend(frameon=False, fontsize=5.6, loc="lower left")

    plt.tight_layout()
    plt.savefig("results/fig2_translational.png", dpi=400, bbox_inches="tight")
    plt.close()
    print("写入 results/fig2_translational.png")


if __name__ == "__main__":
    fig2()
    fig3()
