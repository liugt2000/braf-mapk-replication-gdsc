"""
28. 参考文献完整性审计：
  1) 正文编号是否在 1..N 内连续、无跳号
  2) 首次出现顺序是否严格递增（Vancouver 要求）
  3) 文献表条目数是否与 N 一致
  4) 每条文献表条目是否都在正文被引用
"""
import re

SRC = "manuscript_JTM_v2_BRAF.md"
t = open(SRC, encoding="utf-8").read()

body, refs = t.split("## References", 1)
entries = [l for l in refs.splitlines() if re.match(r"^\d+\.\s", l)]
listed = [int(l.split(".")[0]) for l in entries]

cited_seq = []
for m in re.finditer(r"\[([0-9,–\-]+)\]", body):
    for part in m.group(1).split(","):
        if "–" in part or "-" in part:
            a, b = re.split(r"[–\-]", part)
            cited_seq.extend(range(int(a), int(b) + 1))
        else:
            cited_seq.append(int(part))

first_order = []
for n in cited_seq:
    if n not in first_order:
        first_order.append(n)

N = len(listed)
ok = True


def check(cond, msg):
    global ok
    print(("  ✓ " if cond else "  ✗ ") + msg)
    if not cond:
        ok = False


print(f"文献表条目: {N}")
print(f"正文引用次数: {len(cited_seq)}  去重: {len(set(cited_seq))}")
print()
print("【1】编号范围与连续性")
check(listed == list(range(1, N + 1)), f"文献表编号 1..{N} 连续且升序")
check(set(cited_seq) <= set(range(1, N + 1)),
      "正文编号均在 1..N 内" if set(cited_seq) <= set(range(1, N + 1))
      else f"越界编号: {sorted(set(cited_seq) - set(range(1, N + 1)))}")
missing = set(range(1, N + 1)) - set(cited_seq)
check(not missing, "无未被引用的文献" if not missing else f"未被引用: {sorted(missing)}")

print("\n【2】首次出现顺序严格递增")
check(first_order == list(range(1, N + 1)),
      "首次出现顺序 = 1,2,...,N" if first_order == list(range(1, N + 1))
      else f"实际顺序: {first_order}")

print("\n【3】占位符与格式残留")
check("{{" not in t, "无未替换占位符")
check(not re.search(r"\[\s*\]", t), "无空引用括号")
bad_style = re.findall(r"\[,|,\s*\]|\[\s+\d|,\s*,|\[\d+\s+\d+\]", t)
check(not bad_style, f"无格式异常括号 {bad_style}" if bad_style else "无格式异常括号")

print("\n【4】DOI 覆盖率")
withdoi = [l for l in entries if "doi:" in l.lower()]
check(len(withdoi) == N, f"{len(withdoi)}/{N} 条含 DOI")

print("\n" + ("全部通过" if ok else "存在未通过项"))
