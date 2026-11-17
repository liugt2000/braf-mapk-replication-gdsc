"""07_figure1.py  主图：标签取证 -> 可靠性 -> 信号分解"""
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 8,
                     "axes.linewidth": 0.7, "axes.spines.top": False,
                     "axes.spines.right": False, "figure.dpi": 300})
C1, C2, C3 = "#C0392B", "#2C6E9B", "#7F8C8D"

fit = pd.read_csv("data/raw/GDSC2_fitted_dose_response_24Jul22.csv")
crc = fit[fit.TCGA_DESC == "COREAD"].copy()
crc["ic50"] = np.exp(crc.LN_IC50)
crc["ident"] = (crc.ic50 >= crc.MIN_CONC) & (crc.ic50 <= crc.MAX_CONC)

fig, axes = plt.subplots(2, 2, figsize=(7.2, 5.6))

# ---- A: 5-FU IC50 vs 测试上限 ----
ax = axes[0, 0]
d = crc[crc.DRUG_NAME == "5-Fluorouracil"]
d = d.sort_values("ic50")
x = np.arange(len(d))
ax.axhline(32, color=C1, lw=1.0, ls="--")
ax.text(25, 34, "max tested conc. 32 $\\mu$M", color=C1, fontsize=6.2, va="bottom")
ax.scatter(x[d.ident.values], d.ic50.values[d.ident.values], s=16,
           color=C2, label=f"measured (n={int(d.ident.sum())})", zorder=3)
ax.scatter(x[~d.ident.values], d.ic50.values[~d.ident.values], s=16,
           color=C1, marker="^", label=f"extrapolated (n={int((~d.ident).sum())})", zorder=3)
ax.set_yscale("log")
ax.set_xlabel("CRC cell line (GDSC2)")
ax.set_ylabel("5-FU IC$_{50}$ ($\\mu$M)")
ax.set_title("A  Only 37% of 5-FU IC$_{50}$ values are measurements", loc="left", fontsize=8.5, weight="bold")
ax.legend(frameon=False, fontsize=6.5, loc="lower right")
ax.set_ylim(1, 3e4)

# ---- B: 各药外推比例 ----
ax = axes[0, 1]
t1 = pd.read_csv("results/table1_drug_panel.csv")
t1 = t1.sort_values("pct_extrapolated_high")
cols = [C1 if v > 50 else C2 for v in t1.pct_extrapolated_high]
ax.barh(t1.drug, t1.pct_extrapolated_high, color=cols, height=0.65)
ax.axvline(50, color=C3, lw=0.8, ls=":")
for y, v in enumerate(t1.pct_extrapolated_high):
    ax.text(v + 1.5, y, f"{v:.0f}%", va="center", fontsize=6.5)
ax.set_xlim(0, 118)
ax.set_xlabel("% of cell lines with IC$_{50}$ beyond tested range")
ax.set_title("B  Censoring is the rule, not the exception", loc="left", fontsize=8.5, weight="bold")
ax.tick_params(labelsize=6.5)

# ---- C: 跨活动重测信度 ----
ax = axes[1, 0]
t3 = pd.read_csv("results/table3_cross_campaign.csv")
t3 = t3[t3.subset.str.startswith("drug:")].copy()
t3["name"] = t3.subset.str.replace("drug: ", "", regex=False)
t3 = t3.sort_values("pearson_r")
cols = [C1 if r < 0.7 else C2 for r in t3.pearson_r]
ax.barh(t3.name, t3.pearson_r, color=cols, height=0.6)
for y, (r, f) in enumerate(zip(t3.pearson_r, t3.median_fold_diff)):
    ax.text(r + 0.02, y, f"r={r:.2f}, {f:.1f}x", va="center", fontsize=6.3)
ax.set_xlim(0, 1.25)
ax.set_xlabel("Pearson r between GDSC1 and GDSC2 (same line, same drug)")
ax.set_title("C  Two independent screens disagree by 2-5 fold", loc="left", fontsize=8.5, weight="bold")
ax.tick_params(labelsize=6.5)

# ---- D: 信号分解 ----
ax = axes[1, 1]
labels = ["Overall\nR$^2$", "Sensitive vs\nresistant\n(AUC)", "Rank within\nresponsive\n(Spearman)", "Rank within\nresistant\n(Spearman)"]
vals = [0.286, 0.778, 0.329, 0.222]
alt = [-0.358, None, -0.358, None]
x = np.arange(len(vals))
cols = [C2, C2, C1, C1]
ax.bar(x, vals, color=cols, width=0.6)
ax.scatter([0, 2], [-0.358, -0.358], marker="x", s=40, color="black", zorder=4,
           label="alternative specification")
ax.axhline(0, color="black", lw=0.7)
for i, v in enumerate(vals):
    ax.text(i, v + 0.03, f"{v:.2f}", ha="center", fontsize=6.8)
ax.set_xticks(x); ax.set_xticklabels(labels, fontsize=6.2)
ax.set_ylim(-0.6, 1.0)
ax.set_ylabel("performance")
ax.set_title("D  The signal is a censoring boundary, not potency", loc="left", fontsize=8.5, weight="bold")
ax.legend(frameon=False, fontsize=6.2, loc="upper right")

plt.tight_layout()
plt.savefig("results/figure1_forensics.png", bbox_inches="tight")
print("saved results/figure1_forensics.png")
