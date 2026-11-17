"""12_figures.py  图 2 误差预算 / 图 3 可靠性-可预测性 / 图 4 病人层"""
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 8,
                     "axes.linewidth": 0.7, "axes.spines.top": False,
                     "axes.spines.right": False, "figure.dpi": 300})
RED, BLUE, GREY, GOLD = "#C0392B", "#2C6E9B", "#7F8C8D", "#D68910"

# ============ 图 2：误差预算 ============
w = json.load(open("results/11_well_noise.json"))
t2 = pd.read_csv("results/table2_replicate_reliability.csv")
t3 = pd.read_csv("results/table3_cross_campaign.csv")
t1 = pd.read_csv("results/table1_drug_panel.csv")

well_cv = w["replicate_cv_median"]
snr = w["snr_median"]
rep_fold = float(t2[t2.scheme.str.startswith("loose")].typical_pair_fold_ratio.iloc[0])
cross_fold = float(t3[t3.subset == "drug: 5-Fluorouracil"].median_fold_diff.iloc[0])
fu_row = t1[t1.drug == "5-Fluorouracil"].iloc[0]
cens_fold = 13924 / 35.4     # 5-FU 外推范围 / 外推下限

items = [
    ("Well-to-well replicate\n(CV 9.4%)", 1 / (1 - well_cv), BLUE),
    ("Single-dose ranking\n(SNR 1.5)", snr, BLUE),
    ("Same line, same drug,\ndifferent screen (6.0x)", rep_fold, GOLD),
    ("GDSC1 vs GDSC2, same line\nand drug (3.5x)", cross_fold, GOLD),
    ("Extrapolated IC$_{50}$\nrange (400x)", cens_fold, RED),
]

fig, ax = plt.subplots(figsize=(4.6, 2.9))
ys = np.arange(len(items))[::-1]
for y, (lab, v, c) in zip(ys, items):
    ax.barh(y, v, color=c, height=0.62)
    ax.text(v * 1.12, y, f"{v:.1f}x" if v < 100 else f"{v:.0f}x",
            va="center", fontsize=7, weight="bold")
ax.set_yticks(ys)
ax.set_yticklabels([i[0] for i in items], fontsize=6.4)
ax.set_xscale("log")
ax.set_xlim(1, 2000)
ax.set_xlabel("fold-uncertainty contributed (log scale)")
ax.axvline(2, color="black", lw=0.7, ls=":")
ax.text(2.15, -0.55, "2-fold clinical decision threshold", fontsize=5.6, color="black")
ax.set_ylim(-0.9, len(items) - 0.4)
ax.set_title("Error budget: batch effects and censoring dominate, pipetting does not",
             loc="left", fontsize=7.4, weight="bold", pad=10)
plt.tight_layout()
plt.savefig("results/figure2_error_budget.png", bbox_inches="tight")
plt.close()

# ============ 图 3：可靠性 × 可预测性 ============
m5 = pd.read_csv("results/05_model_vs_ceiling.csv")
rel = t3[t3.subset.str.startswith("drug:")].copy()
rel["drug"] = rel.subset.str.replace("drug: ", "", regex=False)
mg = m5.merge(rel[["drug", "pearson_r", "median_fold_diff"]], on="drug", how="inner")

fig, ax = plt.subplots(figsize=(4.4, 3.6))
ax.axhline(0, color="black", lw=0.7)
ax.axvline(0.6, color=GREY, lw=0.7, ls="--")
OFF = {"5-Fluorouracil": (0.012, 0.012, "left"),
       "Trametinib": (0.012, -0.030, "left"),
       "Dabrafenib": (-0.012, 0.012, "right"),
       "Selumetinib": (0.012, -0.014, "left"),
       "SN-38": (0.012, 0.010, "left"),
       "Fulvestrant": (0.012, 0.010, "left")}
for _, r in mg.iterrows():
    ax.scatter(r.pearson_r, r.model_cv_r2, s=48,
               color=BLUE if r.beats_null else RED, zorder=3,
               edgecolor="white", linewidth=0.6)
    dx, dy, ha = OFF.get(r.drug, (0.012, 0.010, "left"))
    ax.annotate(r.drug, (r.pearson_r, r.model_cv_r2), fontsize=5.8,
                xytext=(r.pearson_r + dx, r.model_cv_r2 + dy), ha=ha)
ax.set_xlim(0.30, 0.82)
ax.set_ylim(-0.30, 0.42)
ax.text(0.605, 0.395, "r = 0.60", fontsize=5.6, color=GREY)
ax.set_xlabel("label reliability\n(Pearson r, GDSC1 vs GDSC2)")
ax.set_ylabel("model predictability\n(cross-validated R$^2$)")
ax.set_title("Reliability and predictability are dissociated", loc="left",
             fontsize=7.6, weight="bold", pad=8)
ax.scatter([], [], color=BLUE, s=30, label="beats permutation null")
ax.scatter([], [], color=RED, s=30, label="fails null")
ax.legend(frameon=False, fontsize=5.8, loc="lower left",
          handletextpad=0.3, borderpad=0.2)
plt.tight_layout()
plt.savefig("results/figure3_reliability_vs_predictability.png", bbox_inches="tight")
plt.close()

# ============ 图 4：病人层 ============
t4 = pd.read_csv("results/table4_interaction.csv")
st = pd.read_csv("results/table5_stratified_event_rates.csv", index_col=[0, 1])

fig, axes = plt.subplots(1, 2, figsize=(7.2, 2.8),
                         gridspec_kw={"width_ratios": [1.25, 1]})
ax = axes[0]
sub = t4[t4.subset == "stage II-III"].reset_index(drop=True)
ys = np.arange(len(sub))[::-1]
for y, r in zip(ys, sub.itertuples()):
    c = BLUE if r.score == "score_5FU" else GREY
    ax.errorbar(r.HR_int, y, xerr=[[r.HR_int - r.int_lo], [r.int_hi - r.HR_int]],
                fmt="o", color=c, capsize=2.5, ms=4, lw=1.1)
    ax.text(2.35, y, f"{r.HR_int:.2f} (p={r.p_int:.2f})", fontsize=5.8, va="center")
ax.axvline(1, color="black", lw=0.7, ls="--")
ax.set_yticks(ys)
ax.set_yticklabels([s.replace("score_", "") for s in sub.score], fontsize=6.5)
ax.set_xlim(0.4, 2.4)
ax.set_xlabel("HR for score $\\times$ adjuvant chemotherapy (stage II-III)")
ax.set_title("A  No predictive interaction in patients", loc="left",
             fontsize=7.6, weight="bold")
ax.text(0.45, len(sub) - 0.55, "n = 232, 90 deaths", fontsize=5.8, color=GREY)

ax = axes[1]
st2 = st.reset_index()
w_ = 0.36
for i, (ch, lab, c) in enumerate([(0, "no chemo", GREY), (1, "chemo", RED)]):
    s = st2[st2.chemo == ch]
    ax.bar(np.arange(len(s)) + (i - 0.5) * w_, s["mean"] * 100, width=w_,
           color=c, label=lab)
    for j, v in enumerate(s["mean"] * 100):
        ax.text(j + (i - 0.5) * w_, v + 1.2, f"{v:.0f}", ha="center", fontsize=5.6)
ax.set_xticks(np.arange(3))
ax.set_xticklabels(["low", "mid", "high"], fontsize=6.5)
ax.set_xlabel("predicted 5-FU resistance score (tertile)")
ax.set_ylabel("death rate (%)")
ax.set_ylim(0, 62)
ax.legend(frameon=False, fontsize=5.8)
ax.set_title("B  Chemo helps in every stratum", loc="left",
             fontsize=7.6, weight="bold")
plt.tight_layout()
plt.savefig("results/figure4_patient_layer.png", bbox_inches="tight")
print("saved figure2/3/4")
