"""
把正文的引用占位符 {{KEY,KEY}} 替换为首次出现顺序的 Vancouver 编号，
并重排参考文献表。带断言，任何一处未命中即中止。
"""
import re
import sys

SRC = "manuscript_JTM_v2_BRAF.md"

# 已逐条用 Crossref / Europe PMC 核对过的书目（键 -> 格式化字符串，不含编号）
REFS = {
    "GDSC_2013": "Yang W, Soares J, Greninger P, et al. Genomics of Drug Sensitivity in Cancer (GDSC): a resource for therapeutic biomarker discovery in cancer cells. *Nucleic Acids Res.* 2013;41:D955–61.",
    "IORIO_2016": "Iorio F, Knijnenburg TA, Vis DJ, et al. A landscape of pharmacogenomic interactions in cancer. *Cell.* 2016;166:740–54.",
    "GHANDI_2019": "Ghandi M, Huang FW, Jané-Valbuena J, et al. Next-generation characterization of the Cancer Cell Line Encyclopedia. *Nature.* 2019;569:503–8.",
    "BARRETINA_2012": "Barretina J, Caponigro G, Stransky N, et al. The Cancer Cell Line Encyclopedia enables predictive modelling of anticancer drug sensitivity. *Nature.* 2012;483:603–7.",
    "SHOEMAKER_2006": "Shoemaker RH. The NCI60 human tumour cell line anticancer drug screen. *Nat Rev Cancer.* 2006;6:813–23.",
    "CORSELLO_2020": "Corsello SM, Nagari RT, Spangler RD, et al. Discovering the anti-cancer potential of non-oncology drugs by systematic viability profiling. *Nat Cancer.* 2020;1:235–48.",
    "HAIBAKAINS_2013": "Haibe-Kains B, El-Hachem N, Birkbak NJ, et al. Inconsistency in large pharmacogenomic studies. *Nature.* 2013;504:389–93.",
    "SAFIKHANI_2016": "Safikhani Z, Smirnov P, Freeman M, et al. Revisiting inconsistency in large pharmacogenomic studies. *F1000Res.* 2016;5:2333.",
    "BENJAMINI_1995": "Benjamini Y, Hochberg Y. Controlling the false discovery rate: a practical and powerful approach to multiple testing. *J R Stat Soc Series B.* 1995;57:289–300.",
    "SAUERBREI_2018": "Sauerbrei W, Taube SE, McShane LM, et al. Reporting recommendations for tumor marker prognostic studies (REMARK). *J Natl Cancer Inst.* 2018;110:803–11.",
    "WOLFF_2019": "Wolff RF, Moons KGM, Riley RD, et al. PROBAST: a tool to assess the risk of bias and applicability of prediction model studies. *Ann Intern Med.* 2019;170:51–8.",
    "COLLINS_2024": "Collins GS, Moons KGM, Dhiman P, et al. TRIPOD+AI statement: updated guidance for reporting clinical prediction models that use regression or machine learning methods. *BMJ.* 2024;385:e078378.",
    "YAEGER_2018": "Yaeger R, Chatila WK, Lipsyc MD, et al. Clinical sequencing defines the genomic landscape of metastatic colorectal cancer. *Cancer Cell.* 2018;33:125–136.e3.",
    "KOPETZ_2019": "Kopetz S, Grothey A, Yaeger R, et al. Encorafenib, binimetinib, and cetuximab in BRAF V600E-mutated colorectal cancer. *N Engl J Med.* 2019;381:1632–43.",
    "CORCORAN_2012": "Corcoran RB, Ebi H, Turke AB, et al. EGFR-mediated re-activation of MAPK signaling contributes to insensitivity of BRAF mutant colorectal cancers to RAF inhibition with vemurafenib. *Cancer Discov.* 2012;2:227–35.",
    "PRAHALLAD_2012": "Prahallad A, Sun C, Huang S, et al. Unresponsiveness of colon cancer to BRAF(V600E) inhibition through feedback activation of EGFR. *Nature.* 2012;483:100–3.",
    "GUINNEY_2015": "Guinney J, Dienstmann R, Wang X, et al. The consensus molecular subtypes of colorectal cancer. *Nat Med.* 2015;21:1350–6.",
    "RIBIC_2003": "Ribic CM, Sargent DJ, Moore MJ, et al. Tumor microsatellite-instability status as a predictor of benefit from fluorouracil-based adjuvant chemotherapy for colon cancer. *N Engl J Med.* 2003;349:247–57.",
    "TSHERNIAK_2017": "Tsherniak A, Vazquez F, Montgomery PG, et al. Defining a Cancer Dependency Map. *Cell.* 2017;170:564–576.e16.",
    "VLACHOGIANNIS_2018": "Vlachogiannis G, Hedayat S, Vatsiou A, et al. Patient-derived organoids model treatment response of metastatic gastrointestinal cancers. *Science.* 2018;359:920–6.",
}

text = open(SRC, encoding="utf-8").read()

# ---- 1. 编号：按占位符首次出现顺序 ----
order = []
for m in re.finditer(r"\{\{([A-Z0-9_,]+)\}\}", text):
    for k in m.group(1).split(","):
        k = k.strip()
        if k not in REFS:
            sys.exit(f"未知引用键: {k}")
        if k not in order:
            order.append(k)
num = {k: i + 1 for i, k in enumerate(order)}


def fmt(nums):
    nums = sorted(set(nums))
    out, i = [], 0
    while i < len(nums):
        j = i
        while j + 1 < len(nums) and nums[j + 1] == nums[j] + 1:
            j += 1
        if j - i >= 2:
            out.append(f"{nums[i]}–{nums[j]}")
        else:
            out.extend(str(n) for n in nums[i:j + 1])
        i = j + 1
    return "[" + ",".join(out) + "]"


text = re.sub(r"\{\{([A-Z0-9_,]+)\}\}",
              lambda m: fmt([num[k.strip()] for k in m.group(1).split(",")]), text)

# ---- 2. 重建参考文献表 ----
missing = [k for k in REFS if k not in order]
if missing:
    sys.exit("以下文献在正文中未被引用: " + ", ".join(missing))

ref_block = "## References\n\n" + "\n".join(
    f"{num[k]}. {REFS[k]}" for k in order) + "\n"
text = re.sub(r"## References\n\n(?:[0-9]+\. .*\n|\n)*$", ref_block, text)

open(SRC, "w", encoding="utf-8").write(text)
print(f"共 {len(order)} 条文献，编号按首次出现顺序：")
for k in order:
    print(f"  [{num[k]:2d}] {k}")
