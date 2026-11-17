"""
08_build_patient_score.py
把 CCLE 上训练的 5-FU 敏感性模型迁移到病人队列
-------------------------------------------------
1. 在 CCLE(39 株 CRC) 上拟合弹性网，取非零系数基因作为 signature
2. 把基因名映射为 Entrez ID
3. 从 cBioPortal 拉取病人队列表达值，队列内逐基因 z 标准化
4. score = sum(beta_i * z_i) / sqrt(sum(beta_i^2))   ——平台无关的线性预测子
5. 同时计算增殖分数与随机基因集分数（阴性对照）

输出：results/patient_scores.csv, results/08_signature.json
"""
import gzip, json, re, time, warnings
import numpy as np
import pandas as pd
import requests
from sklearn.linear_model import ElasticNetCV
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
warnings.filterwarnings("ignore")

API = "https://www.cbioportal.org/api"
rng = np.random.default_rng(11)


def norm(s):
    return re.sub(r"[^A-Z0-9]", "", str(s).upper())


def api_get(path, **kw):
    for i in range(4):
        try:
            r = requests.get(f"{API}{path}", timeout=90, **kw)
            if r.status_code == 200:
                return r.json()
        except Exception:
            pass
        time.sleep(2)
    raise RuntimeError(f"GET failed: {path}")


def api_post(path, payload):
    for i in range(4):
        try:
            r = requests.post(f"{API}{path}", json=payload, timeout=180)
            if r.status_code == 200:
                return r.json()
        except Exception:
            pass
        time.sleep(2)
    raise RuntimeError(f"POST failed: {path}")


# ---------- 1. CCLE 上拟合 5-FU 模型 ----------
with gzip.open("data/raw/CCLE_RNAseq_genes_rpkm_20180929.gct.gz", "rt") as f:
    f.readline(); f.readline()
    samples = f.readline().rstrip("\n").split("\t")[2:]
lk = {}
for s in samples:
    lk.setdefault(norm(s.split("_", 1)[0]), []).append(s)
lk = {k: v[0] for k, v in lk.items()}

fit = pd.read_csv("data/raw/GDSC2_fitted_dose_response_24Jul22.csv")
crc = fit[fit.TCGA_DESC == "COREAD"].copy()
crc["ccle"] = crc.CELL_LINE_NAME.map(lambda n: lk.get(norm(n)))
d = (crc[crc.DRUG_NAME == "5-Fluorouracil"].dropna(subset=["ccle"])
     .drop_duplicates("ccle"))
need = sorted(crc.ccle.dropna().unique())

raw = pd.read_csv("data/raw/CCLE_RNAseq_genes_rpkm_20180929.gct.gz", sep="\t",
                  skiprows=2, usecols=["Name", "Description"] + need,
                  dtype={c: "float32" for c in need}, compression="gzip")
raw = raw.dropna(subset=["Description"])
expr = np.log2(raw.set_index("Description")[need].groupby(level=0).max() + 1)
expr = expr.loc[expr.var(axis=1).sort_values(ascending=False).index[:2000]]

X = expr[d.ccle.values].T.values.astype(np.float64)
y = d.LN_IC50.values
pipe = Pipeline([("sc", StandardScaler()),
                 ("en", ElasticNetCV(l1_ratio=[0.5], cv=4,
                                     alphas=np.logspace(-2, 1.5, 15),
                                     max_iter=2000, random_state=0, n_jobs=1))])
pipe.fit(X, y)
beta = pipe.named_steps["en"].coef_
genes = expr.index.values
sig = pd.Series(beta, index=genes)
sig = sig[sig != 0].sort_values(key=abs, ascending=False)
print(f"5-FU signature: {len(sig)} 个非零基因 (原始 2000)")
print("  上调(耐药) top10:", sig.head(10).index.tolist())
print("  下调(敏感) top10:", sig.tail(10).index.tolist())

# ---------- 2. Entrez 映射 ----------
print("拉取基因注释 ...")
allg = api_get("/genes")
sym2ent = {g["hugoGeneSymbol"]: g["entrezGeneId"] for g in allg}
sig_ent = {g: sym2ent[g] for g in sig.index if g in sym2ent}
print(f"  映射成功 {len(sig_ent)}/{len(sig)}")

# 阴性对照：同数量随机基因（表达量匹配）
pool = [g for g in expr.index if g in sym2ent and g not in sig.index]
ctrl_genes = list(rng.choice(pool, size=len(sig_ent), replace=False))
ctrl_ent = {g: sym2ent[g] for g in ctrl_genes}

PROLIF = [g for g in ["MKI67", "PCNA", "CCNB1", "CDK1", "TOP2A", "BUB1", "AURKA",
                      "PLK1", "CDC20", "CCNA2", "RRM2", "MCM2", "MCM6"]
          if g in sym2ent]
prolif_ent = {g: sym2ent[g] for g in PROLIF}

# ---------- 3. 病人队列 ----------
COHORTS = {
    "coadread_tcga": "coadread_tcga_rna_seq_v2_mrna",
    "coad_silu_2022": "coad_silu_2022_rna_seq_mrna",
}
out = []
for study, prof in COHORTS.items():
    print(f"\n=== {study} ({prof}) ===")
    ent = list(set(list(sig_ent.values()) + list(ctrl_ent.values()) + list(prolif_ent.values())))
    rows = api_post(f"/molecular-profiles/{prof}/molecular-data/fetch",
                    {"sampleListId": f"{study}_all", "entrezGeneIds": ent})
    md = pd.DataFrame(rows)
    print(f"  {md.sampleId.nunique()} 样本 x {md.entrezGeneId.nunique()} 基因")
    mat = md.pivot_table(index="sampleId", columns="entrezGeneId", values="value")
    mat = np.log2(mat + 1)
    z = (mat - mat.mean()) / mat.std(ddof=0)      # 队列内逐基因标准化

    def score(ent_map):
        ids = [v for v in ent_map.values() if v in z.columns]
        b = np.array([sig[sig.index[i]] for i, g in enumerate(ent_map) if ent_map[g] in z.columns]
                     ) if ent_map is sig_ent else np.ones(len(ids))
        vals = z[ids].values
        return (vals @ b) / np.sqrt((b ** 2).sum())

    s = pd.Series(score(sig_ent), index=z.index, name="score_5FU")
    sc = pd.Series(z[[v for v in ctrl_ent.values() if v in z.columns]].values.mean(1),
                   index=z.index, name="score_random")
    sp = pd.Series(z[[v for v in prolif_ent.values() if v in z.columns]].values.mean(1),
                   index=z.index, name="score_prolif")
    out.append(pd.concat([s, sc, sp], axis=1).assign(study=study))
    print(f"  score_5FU  range [{s.min():.2f}, {s.max():.2f}]  sd={s.std():.2f}")

res = pd.concat(out)
res.index.name = "sampleId"
res.to_csv("results/patient_scores.csv")
print("\n已写出 results/patient_scores.csv", res.shape)

json.dump({"n_signature_genes": int(len(sig)),
           "n_mapped": int(len(sig_ent)),
           "up_genes": sig.head(20).index.tolist(),
           "down_genes": sig.tail(20).index.tolist(),
           "n_samples": {k: int((res.study == k).sum()) for k in COHORTS}},
          open("results/08_signature.json", "w", encoding="utf-8"),
          indent=2, ensure_ascii=False)
