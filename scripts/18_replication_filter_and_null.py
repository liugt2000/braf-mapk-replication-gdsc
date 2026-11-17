"""
18_replication_filter_and_null.py
---------------------------------
论文的统计骨干。

主张: 在 GDSC 两个独立筛选活动之间, 按预设的复制标准过滤后, 1133 个可检验的
(驱动基因, 药物) 组合里只有极少数能存活, 且存活者全部是 BRAF -> MAPK 抑制剂。

本脚本必须回答三个审稿人一定会问的问题:
  Q1 严格过滤后到底剩几个? 分别是什么?
  Q2 这个数量是否只是多重检验的偶然产物?  -> 置换零分布(保留改变频率)
  Q3 通过 FDR 校正后是否仍显著?

置换零的构造(关键): 对每个基因, 随机抽取同等数量的细胞系作为"改变组",
两个筛选活动共用同一份随机分组(从而保留跨活动的相关结构)。
这样零假设是"任意一组细胞系都会产生同样大的双活动同向效应"。
"""
import re
import numpy as np
import pandas as pd

BASE = "data/raw"
DRIVERS = ["KRAS", "NRAS", "BRAF", "MAP2K1", "NF1", "PIK3CA", "PTEN", "AKT1",
           "TP53", "CDKN2A", "RB1", "ATM", "APC", "CTNNB1", "RNF43", "SMAD2",
           "SMAD4", "TGFBR2", "MLH1", "MSH2", "MSH6", "PMS2", "BRCA1", "BRCA2",
           "ERCC1", "ERCC2", "ERCC5", "ERBB2", "EGFR", "MET", "ARID1A",
           "FBXW7", "POLE"]
MIN_N_ALT, MIN_N_WT, MIN_TOTAL = 8, 20, 40
N_PERM = 200
RNG = np.random.default_rng(20261116)


def norm(s):
    return re.sub(r"[^A-Z0-9]", "", str(s).upper())


def load_mut():
    m = pd.read_csv(f"{BASE}/mutations_latest.csv.gz")
    pm = m.protein_mutation.fillna("").astype(str)
    m["damaging"] = (m.cancer_driver.astype(str).str.lower() == "true") | \
                    pm.str.contains(r"\*|fs", regex=True)
    m = m[m.damaging & m.gene_symbol.isin(DRIVERS)]
    return m[["model_name", "gene_symbol"]].drop_duplicates()


def main():
    mut = load_mut()
    alt = {g: {norm(x) for x in mut.loc[mut.gene_symbol == g, "model_name"]}
           for g in DRIVERS}

    frames, mats = {}, {}
    for fn, tag in [("GDSC2_fitted_dose_response_24Jul22.csv", "G2"),
                    ("GDSC1_fitted_dose_response_24Jul22.csv", "G1")]:
        d = pd.read_csv(f"{BASE}/{fn}")
        d["ic50"] = np.exp(d.LN_IC50)
        d["ident"] = (d.ic50 >= d.MIN_CONC) & (d.ic50 <= d.MAX_CONC)
        d = d[d.ident]
        P = d.pivot_table(index="CELL_LINE_NAME", columns="DRUG_NAME",
                          values="LN_IC50", aggfunc="mean")
        n = P.notna().sum(axis=0)
        P = P[n[n >= MIN_TOTAL].index]
        frames[tag] = P
        print(f"{tag}: {P.shape[0]} 系 x {P.shape[1]} 药 (每药 >= {MIN_TOTAL} 可识别样本)")

    # 统一行索引(细胞系集合 = 两活动并集)
    all_lines = sorted(set(frames["G2"].index) | set(frames["G1"].index))
    idx = {l: i for i, l in enumerate(all_lines)}
    mask_alt = np.zeros((len(DRIVERS), len(all_lines)), dtype=bool)
    for i, g in enumerate(DRIVERS):
        for l in all_lines:
            if norm(l) in alt[g]:
                mask_alt[i, idx[l]] = True
    print(f"细胞系全集 {len(all_lines)} 株")

    def to_mat(P):
        M = np.full((P.shape[1], len(all_lines)), np.nan)
        for j, drug in enumerate(P.columns):
            col = P[drug]
            ok = col.notna()
            M[j, [idx[l] for l in col.index[ok]]] = col[ok].values
        return M, list(P.columns)

    M2, drugs2 = to_mat(frames["G2"])
    M1, drugs1 = to_mat(frames["G1"])
    common_drugs = [d for d in drugs2 if d in set(drugs1)]
    # 把两个矩阵都对齐到共同药物(顺序一致)
    M2 = M2[[drugs2.index(d) for d in common_drugs], :]
    M1 = M1[[drugs1.index(d) for d in common_drugs], :]
    print(f"两活动共同的药物 {len(common_drugs)}")

    def stats_for(M_mat, drugs, masks):
        """masks: (n_gene, n_line) bool。返回 d 与 p 矩阵 (n_gene, n_drug)。"""
        ng, nd = masks.shape[0], len(drugs)
        D = np.full((ng, nd), np.nan)
        T = np.full((ng, nd), np.nan)
        valid = ~np.isnan(M_mat)
        for gi in range(ng):
            mk = masks[gi]
            a_sel = valid & mk[None, :]
            w_sel = valid & (~mk)[None, :]
            na = a_sel.sum(axis=1).astype(float)
            nw = w_sel.sum(axis=1).astype(float)
            Ma = np.where(a_sel, M_mat, np.nan)
            Mw = np.where(w_sel, M_mat, np.nan)
            with np.errstate(invalid="ignore"):
                ma = np.nanmean(Ma, axis=1)
                mw = np.nanmean(Mw, axis=1)
                va = np.nanvar(Ma, axis=1, ddof=1)
                vw = np.nanvar(Mw, axis=1, ddof=1)
            sp = np.sqrt(((na - 1) * va + (nw - 1) * vw)
                         / np.maximum(na + nw - 2, 1))
            with np.errstate(divide="ignore", invalid="ignore"):
                D[gi] = (ma - mw) / sp
            # Welch t -> p (正态近似足够: n>=40)
            se = np.sqrt(va / np.maximum(na, 1) + vw / np.maximum(nw, 1))
            with np.errstate(divide="ignore", invalid="ignore"):
                t = (ma - mw) / se
            T[gi] = 2 * (1 - _norm_cdf(np.abs(t)))
            bad = (na < MIN_N_ALT) | (nw < MIN_N_WT)
            D[gi, bad] = np.nan
            T[gi, bad] = np.nan
        return D, T

    def _norm_cdf(z):
        from scipy.special import erf
        return 0.5 * (1 + erf(np.asarray(z) / np.sqrt(2)))

    # ---- 观测 ----
    print("\n计算观测统计 ...")
    D2, T2 = stats_for(M2, common_drugs, mask_alt)
    D1, T1 = stats_for(M1, common_drugs, mask_alt)
    same = np.sign(D2) == np.sign(D1)
    pass_mask = ((np.abs(D2) >= 0.8) & (np.abs(D1) >= 0.5) & same
                 & (T2 < 1e-3) & (T1 < 1e-3))
    n_pass = int(np.nansum(pass_mask))
    print(f"  可检验组合(两活动均有足够样本): {int(np.isfinite(D2).sum() & np.isfinite(D1).sum())}")
    print(f"  观测通过数: {n_pass}")

    rec = []
    for gi, g in enumerate(DRIVERS):
        for di, drug in enumerate(common_drugs):
            if pass_mask[gi, di]:
                rec.append(dict(gene=g, drug=drug, d_gdsc2=D2[gi, di],
                                d_gdsc1=D1[gi, di], p_gdsc2=T2[gi, di],
                                p_gdsc1=T1[gi, di]))
    obs = pd.DataFrame(rec)
    print("\n通过严格过滤的组合:")
    print(obs.to_string(index=False, float_format=lambda v: f"{v:.3g}")
          if len(obs) else "  无")

    # ---- FDR ----
    from scipy.stats import false_discovery_control
    flat = T2[np.isfinite(T2)]
    q = false_discovery_control(flat, method="bh")
    print(f"\nBH-FDR: {len(flat)} 个检验, q<0.05 的 {int((q<0.05).sum())} 个")
    t2q = np.full(T2.shape, np.nan)
    t2q[np.isfinite(T2)] = q
    sig_fdr = np.isfinite(t2q) & (t2q < 0.05)
    hits_fdr = [(DRIVERS[gi], common_drugs[di], t2q[gi, di])
                for gi in range(len(DRIVERS)) for di in range(len(common_drugs))
                if sig_fdr[gi, di]]
    print(f"  q<0.05 且 |d|>=0.5 的组合 (仅 GDSC2 侧):")
    hh = [(g, d, qq, D2[DRIVERS.index(g), common_drugs.index(d)])
          for g, d, qq in hits_fdr
          if abs(D2[DRIVERS.index(g), common_drugs.index(d)]) >= 0.5]
    for g, d, qq, dd in sorted(hh, key=lambda x: x[3]):
        print(f"    {g:8s} {d:22s} q={qq:.2e} d={dd:+.2f}")
    print(f"  合计 {len(hh)} 个")

    # ---- 置换零 ----
    print(f"\n置换零分布 ({N_PERM} 次) ...")
    n_alt_obs = mask_alt.sum(axis=1)
    null_counts = []
    for b in range(N_PERM):
        rnd = np.zeros_like(mask_alt)
        for gi in range(len(DRIVERS)):
            k = int(n_alt_obs[gi])
            pick = RNG.choice(len(all_lines), size=k, replace=False)
            rnd[gi, pick] = True
        d2b, t2b = stats_for(M2, common_drugs, rnd)
        d1b, t1b = stats_for(M1, common_drugs, rnd)
        sm = np.sign(d2b) == np.sign(d1b)
        pm = ((np.abs(d2b) >= 0.8) & (np.abs(d1b) >= 0.5) & sm
              & (t2b < 1e-3) & (t1b < 1e-3))
        null_counts.append(int(np.nansum(pm)))
        if (b + 1) % 25 == 0:
            print(f"    {b+1}/{N_PERM}  当前零分布均值={np.mean(null_counts):.2f}")

    nc = np.array(null_counts)
    p_emp = (1 + (nc >= n_pass).sum()) / (N_PERM + 1)
    print(f"\n  零分布通过数: 均值={nc.mean():.2f} 中位={np.median(nc):.0f} "
          f"95%分位={np.quantile(nc,.95):.0f} 最大={nc.max()}")
    print(f"  观测 {n_pass} vs 零分布: 经验 p = {p_emp:.4f}")

    pd.DataFrame(dict(perm=np.arange(N_PERM), n_pass=nc)).to_csv(
        "results/18_null_distribution.csv", index=False)
    obs.to_csv("results/18_surviving_pairs.csv", index=False)
    pd.DataFrame(dict(gene=DRIVERS,
                      n_altered=n_alt_obs)).to_csv(
        "results/18_gene_frequencies.csv", index=False)

    # 汇总
    summary = dict(n_testable=int(np.isfinite(D2).sum()),
                   n_pass=n_pass,
                   null_mean=float(nc.mean()),
                   null_q95=float(np.quantile(nc, 0.95)),
                   empirical_p=float(p_emp),
                   fdr_hits=len(hh),
                   only_genes=sorted(obs.gene.unique()) if len(obs) else [])
    pd.Series(summary).to_json("results/18_summary.json", indent=2)
    print("\n写入 results/18_summary.json")


if __name__ == "__main__":
    main()
