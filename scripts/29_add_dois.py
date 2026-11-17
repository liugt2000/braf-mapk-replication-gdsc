"""
29. 给参考文献补 DOI，并向 Crossref 反查每条 DOI 的标题做交叉验证。
标题相似度低于阈值即标记需人工复核。
"""
import difflib
import json
import re
import time
import urllib.request

SRC = "manuscript_JTM_v2_BRAF.md"

DOIS = {
    1: "10.1093/nar/gks1111",
    2: "10.1016/j.cell.2016.06.017",
    3: "10.1038/s41586-019-1186-3",
    4: "10.1038/nature11003",
    5: "10.1038/nrc1951",
    6: "10.1038/s43018-019-0018-6",
    7: "10.1038/nature12831",
    8: "10.12688/f1000research.9611.3",
    9: "10.1111/j.2517-6161.1995.tb02031.x",
    10: "10.1093/jnci/djy088",
    11: "10.7326/m18-1376",
    12: "10.1136/bmj-2023-078378",
    13: "10.1016/j.ccell.2017.12.004",
    14: "10.1056/NEJMoa1908075",
    15: "10.1038/nm.3967",
    16: "10.1056/NEJMoa022289",
    17: "10.1158/2159-8290.CD-11-0341",
    18: "10.1038/nature10868",
    19: "10.1016/j.cell.2017.06.010",
    20: "10.1126/science.aao2774",
}


def cr(doi, tries=3):
    for i in range(tries):
        try:
            req = urllib.request.Request(
                f"https://api.crossref.org/works/{doi}",
                headers={"User-Agent": "refaudit/1.0 (mailto:noreply@example.com)"})
            with urllib.request.urlopen(req, timeout=40) as r:
                return json.loads(r.read().decode())["message"]
        except Exception:
            if i == tries - 1:
                return None
            time.sleep(2)


t = open(SRC, encoding="utf-8").read()
body, refs = t.split("## References", 1)
lines = refs.splitlines()

print(f"{'#':>3} {'DOI':42s} {'标题相似度':>8}  判定")
print("-" * 100)

for n, doi in DOIS.items():
    idx = next((i for i, l in enumerate(lines) if re.match(rf"^{n}\.\s", l)), None)
    if idx is None:
        print(f"{n:3d} 文献表缺该条")
        continue
    line = lines[idx]
    # 抽出文献标题（第一个句号前的部分，去掉作者段）
    m = re.match(rf"^{n}\.\s(.+?)\.\s\*", line)
    local_title = m.group(1) if m else ""
    # 去掉作者前缀：以 "et al." 或 "XX Y." 之后为标题起点
    local_title = re.sub(r"^.*?(?:et al\.|[A-Z]{1,3}\.)\s*", "", local_title)

    meta = cr(doi)
    if meta is None:
        print(f"{n:3d} {doi:42s} {'—':>8}  DOI 无法解析")
        lines[idx] = line.rstrip() + f" doi:{doi}"
        continue
    remote = (meta.get("title") or [""])[0]
    remote = re.sub(r"<[^>]+>", "", remote).strip().rstrip(".")
    sim = difflib.SequenceMatcher(
        None, local_title.lower(), remote.lower()).ratio()
    verdict = "一致" if sim >= 0.85 else "需复核"
    print(f"{n:3d} {doi:42s} {sim:8.2f}  {verdict}")
    if verdict == "需复核":
        print(f"      本地: {local_title[:90]}")
        print(f"      远端: {remote[:90]}")
    lines[idx] = line.rstrip() + f" doi:{doi}"
    time.sleep(0.3)

open(SRC, "w", encoding="utf-8").write(body + "## References" + "\n".join(lines))
print("\nDOI 已写入文献表")
