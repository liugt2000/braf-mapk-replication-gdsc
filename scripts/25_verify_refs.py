"""
核对稿件参考文献：逐条向 Europe PMC 检索，比对标题/期刊/年/卷/页/作者首位。
输出 results/ref_check.csv 并打印需人工修正的条目。
"""
import json
import re
import time
import urllib.parse
import urllib.request

API = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"

REFS = [
    dict(n=1, title="Genomics of Drug Sensitivity in Cancer (GDSC): a resource for therapeutic biomarker discovery in cancer cells",
         journal="Nucleic Acids Res", year=2013, vol="41", pages="D955-61", first="Yang W"),
    dict(n=2, title="A landscape of pharmacogenomic interactions in cancer",
         journal="Cell", year=2016, vol="166", pages="740-54", first="Iorio F"),
    dict(n=3, title="Next-generation characterization of the Cancer Cell Line Encyclopedia",
         journal="Nature", year=2019, vol="569", pages="503-8", first="Ghandi M"),
    dict(n=4, title="The Cancer Cell Line Encyclopedia enables predictive modelling of anticancer drug sensitivity",
         journal="Nature", year=2012, vol="483", pages="603-7", first="Barretina J"),
    dict(n=5, title="Inconsistency in large pharmacogenomic studies",
         journal="Nature", year=2013, vol="504", pages="389-93", first="Haibe-Kains B"),
    dict(n=6, title="Revisiting inconsistency in large pharmacogenomic studies",
         journal="F1000Res", year=2016, vol="5", pages="2333", first="Safikhani Z"),
    dict(n=7, title="Discovering the anti-cancer potential of non-oncology drugs by systematic viability profiling",
         journal="Nat Cancer", year=2020, vol="1", pages="235-48", first="Corsello SM"),
    dict(n=8, title="Defining a Cancer Dependency Map",
         journal="Cell", year=2017, vol="170", pages="564-76", first="Tsherniak A"),
    dict(n=9, title="EGFR-mediated re-activation of MAPK signaling contributes to insensitivity of BRAF mutant colorectal cancers to RAF inhibition with vemurafenib",
         journal="Cancer Discov", year=2012, vol="2", pages="227-35", first="Corcoran RB"),
    dict(n=10, title="Unresponsiveness of colon cancer to BRAF(V600E) inhibition through feedback activation of EGFR",
         journal="Nature", year=2012, vol="483", pages="100-3", first="Prahallad A"),
    dict(n=11, title="Encorafenib, binimetinib, and cetuximab in BRAF V600E-mutated colorectal cancer",
         journal="N Engl J Med", year=2019, vol="381", pages="1632-43", first="Kopetz S"),
    dict(n=12, title="The consensus molecular subtypes of colorectal cancer",
         journal="Nat Med", year=2015, vol="21", pages="1350-6", first="Guinney J"),
    dict(n=13, title="Tumor microsatellite-instability status as a predictor of benefit from fluorouracil-based adjuvant chemotherapy for colon cancer",
         journal="N Engl J Med", year=2003, vol="349", pages="247-57", first="Ribic CM"),
    dict(n=14, title="Patient-derived organoids model treatment response of metastatic gastrointestinal cancers",
         journal="Science", year=2018, vol="359", pages="920-6", first="Vlachogiannis G"),
    dict(n=15, title="Controlling the false discovery rate: a practical and powerful approach to multiple testing",
         journal="J R Stat Soc Series B", year=1995, vol="57", pages="289-300", first="Benjamini Y"),
    dict(n=16, title="TRIPOD+AI statement: updated guidance for reporting clinical prediction models that use regression or machine learning methods",
         journal="BMJ", year=2024, vol="385", pages="e078378", first="Collins GS"),
    dict(n=17, title="PROBAST: a tool to assess the risk of bias and applicability of prediction model studies",
         journal="Ann Intern Med", year=2019, vol="170", pages="51-8", first="Wolff RF"),
    dict(n=18, title="Reporting recommendations for tumor marker prognostic studies (REMARK)",
         journal="J Natl Cancer Inst", year=2018, vol="110", pages="803-11", first="Sauerbrei W"),
]


def get(url, tries=3):
    for i in range(tries):
        try:
            with urllib.request.urlopen(url, timeout=45) as r:
                return json.loads(r.read().decode("utf-8"))
        except Exception as e:
            if i == tries - 1:
                raise
            time.sleep(2)


def norm_pages(p):
    return re.sub(r"[^0-9]", "", str(p or ""))


def search(title):
    q = 'TITLE:"%s"' % title.replace('"', "")
    url = API + "?" + urllib.parse.urlencode({
        "query": q, "format": "json", "resultType": "core", "pageSize": 5})
    return get(url).get("resultList", {}).get("result", [])


rows = []
for r in REFS:
    try:
        hits = search(r["title"])
    except Exception as e:
        print(f"[{r['n']:2d}] 检索失败: {e}")
        continue
    if not hits:
        print(f"[{r['n']:2d}] ✗ 无命中  {r['title'][:70]}")
        rows.append(dict(n=r["n"], status="NOT_FOUND", title=r["title"]))
        continue
    # 优先取标题最接近的
    best = None
    for h in hits:
        t = (h.get("title") or "").rstrip(".")
        if t.lower()[:40] == r["title"].lower()[:40]:
            best = h
            break
    best = best or hits[0]
    ji = best.get("journalInfo", {}) or {}
    jt = (ji.get("journal") or {}).get("title") or best.get("bookOrReportDetails", {}).get("publisher", "")
    got_year = best.get("pubYear")
    got_vol = ji.get("volume")
    got_pages = best.get("pageInfo")
    got_first = (best.get("authorString") or "").split(",")[0].strip()

    issues = []
    if str(got_year) != str(r["year"]):
        issues.append(f"年 {r['year']}→{got_year}")
    if str(got_vol or "") != str(r["vol"]):
        issues.append(f"卷 {r['vol']}→{got_vol}")
    if norm_pages(got_pages)[:4] != norm_pages(r["pages"])[:4]:
        issues.append(f"页 {r['pages']}→{got_pages}")
    if got_first.replace(".", "").split()[0].lower() != r["first"].replace(".", "").split()[0].lower():
        issues.append(f"首作者 {r['first']}→{got_first}")
    if r["journal"].lower().replace(".", "").split()[0] not in jt.lower().replace(".", ""):
        issues.append(f"刊 {r['journal']}→{jt}")

    status = "OK" if not issues else "CHECK"
    mark = "✓" if status == "OK" else "△"
    print(f"[{r['n']:2d}] {mark} {jt} {got_year};{got_vol}:{got_pages}  {got_first}"
          + ("" if not issues else "   << " + "; ".join(issues)))
    rows.append(dict(n=r["n"], status=status, title=(best.get("title") or "")[:90],
                     journal=jt, year=got_year, vol=got_vol, pages=got_pages,
                     first_author=got_first, doi=best.get("doi", ""),
                     pmid=best.get("pmid", ""), issues="; ".join(issues)))
    time.sleep(0.4)

import csv
with open("results/ref_check.csv", "w", newline="", encoding="utf-8-sig") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
    w.writeheader()
    w.writerows(rows)

bad = [r for r in rows if r["status"] != "OK"]
print(f"\n合计 {len(rows)} 条，完全一致 {len(rows)-len(bad)} 条，需修 {len(bad)} 条")
