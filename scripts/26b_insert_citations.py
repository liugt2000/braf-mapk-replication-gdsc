"""在正文插入引用占位符。每处替换断言命中次数==1。"""

SRC = "manuscript_JTM_v2_BRAF.md"
t = open(SRC, encoding="utf-8").read()

EDITS = [
    # Background: 资源清单
    (
        "Public pharmacogenomic resources — GDSC, the Cancer Cell Line Encyclopedia (CCLE), the NCI-60, PRISM — were built",
        "Public pharmacogenomic resources — GDSC {{GDSC_2013,IORIO_2016}}, the Cancer Cell Line Encyclopedia (CCLE) "
        "{{GHANDI_2019,BARRETINA_2012}}, the NCI-60 {{SHOEMAKER_2006}}, PRISM {{CORSELLO_2020}} — were built",
    ),
    # Background: 跨研究不一致
    (
        "Cross-study discordance in large pharmacogenomic datasets is well documented at the level of drug response itself, "
        "but the consequence for *gene-level* claims is usually left implicit.",
        "Cross-study discordance in large pharmacogenomic datasets is well documented at the level of drug response itself "
        "{{HAIBAKAINS_2013,SAFIKHANI_2016}}, but the consequence for *gene-level* claims is usually left implicit.",
    ),
    # Methods: BH 校正
    (
        "and the adjusted *β* for alteration status was reported with Benjamini–Hochberg correction across all pairs.",
        "and the adjusted *β* for alteration status was reported with Benjamini–Hochberg correction "
        "{{BENJAMINI_1995}} across all pairs.",
    ),
    # Methods: 报告规范
    (
        "Reporting follows REMARK; bias assessment follows PROBAST.",
        "Reporting follows REMARK {{SAUERBREI_2018}}; bias assessment follows PROBAST {{WOLFF_2019}}, and the "
        "expression models are reported against TRIPOD+AI {{COLLINS_2024}}.",
    ),
    # Results: BRAF 豁免负荷富集的机制说明（改为可引用的表述）
    (
        "BRAF is the exception (*d* = +0.20): BRAF V600E arises in non-hypermutated, mismatch-repair-proficient "
        "colorectal tumours, so BRAF-mutant lines are drawn from the same burden distribution as wild-type lines.",
        "BRAF is the exception (*d* = +0.20): BRAF V600E is enriched in right-sided, microsatellite-stable colorectal "
        "tumours {{YAEGER_2018}}, so BRAF-mutant lines are not drawn from the hypermutated tail.",
    ),
    # Results: BEACON
    (
        "which is precisely the pattern that motivated combined BRAF and EGFR blockade in the BEACON trial.",
        "which is precisely the pattern that motivated combined BRAF and EGFR blockade in the BEACON trial "
        "{{KOPETZ_2019}}.",
    ),
    # Results: CMS
    (
        "These patients were concentrated in CMS1 (28/47, 60%) versus 10/240 (4%) in wild-type tumours,",
        "These patients were concentrated in CMS1 {{GUINNEY_2015}} (28/47, 60%) versus 10/240 (4%) in wild-type tumours,",
    ),
    # Discussion: 计算对应物 + 机制
    (
        "This is the cellular counterpart of the clinical observation that BRAF-selective monotherapy fails in "
        "colorectal cancer while MEK inhibition retains activity, and it is the reason the effective regimen adds "
        "EGFR blockade.",
        "This is the cellular counterpart of the clinical observation that BRAF-selective monotherapy fails in "
        "colorectal cancer while MEK inhibition retains activity {{CORCORAN_2012,PRAHALLAD_2012}}, and it is the "
        "reason the effective regimen adds EGFR blockade {{KOPETZ_2019}}.",
    ),
    # Limitations: 增加第七、第八条
    (
        "Sixth, we tested 34 pre-specified genes and did not search the genome; genes outside that set are untested, "
        "not excluded.",
        "Sixth, we tested 34 pre-specified genes and did not search the genome; genes outside that set are untested, "
        "not excluded. Seventh, we did not integrate CRISPR-based fitness data for the same models "
        "{{TSHERNIAK_2017}}, which is the natural way to ask whether the colorectal dissociation reflects a "
        "difference in BRAF dependency rather than in drug pharmacology. Eighth, the dissociation rests on "
        "small-molecule screens only; patient-derived organoid panels {{VLACHOGIANNIS_2018}} would be the "
        "appropriate next test.",
    ),
]

for old, new in EDITS:
    n = t.count(old)
    assert n == 1, f"命中 {n} 次（应为 1）: {old[:70]}"
    t = t.replace(old, new)

# Discussion 新增一段：用 Ribic 说明"不可检出 ≠ 不存在"
anchor = ("That asymmetry is what makes the replication filter more than a robustness check.")
assert t.count(anchor) == 1
new_para = (
    "Absence from our list also does not mean absence of a real interaction. The benefit of fluorouracil-based "
    "adjuvant chemotherapy is absent in mismatch-repair-deficient colon cancer {{RIBIC_2003}} — a gene–drug "
    "interaction that has held up clinically for two decades — yet it could not have emerged here: 5-FU IC₅₀ is "
    "censored in 63% of colorectal lines, and mismatch-repair status was not in the pre-specified gene set. The "
    "correct reading of \"only one pair replicated\" is therefore *only one pair replicated among those this "
    "resource can actually measure*, and the resource's measurement range is far narrower than its catalogue.\n\n"
)
t = t.replace(anchor, new_para + anchor)

open(SRC, "w", encoding="utf-8").write(t)
print(f"插入完成，剩余占位符 {t.count('{{')} 处待编号")
