# -*- coding: utf-8 -*-
"""
paper/ 안에서 두 문서가 같은 라벨로 진술하는 것들을 뽑아 나란히 보여준다.
병합 전에 무엇이 갈라져 있는지 눈으로 확인하기 위한 일회용 도구.

    python gate/_diff_shared.py
"""

import difflib
import io
import os
import re
import sys

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
PAPER = os.path.join(ROOT, "paper")
KINDS = "theorem|proposition|lemma|corollary|conjecture"


def grab(path):
    src = io.open(path, encoding="utf-8", newline="").read()
    out = {}
    pat = re.compile(
        r"\\begin\{(" + KINDS + r")\}(\[[^\]]*\])?\\label\{([^}]*)\}"
        r"(.*?)\\end\{\1\}", re.S)
    for m in pat.finditer(src):
        out[m.group(3)] = (m.group(1), m.group(2) or "", m.group(4).strip())
    return out


def norm(t):
    return [l.rstrip() for l in t.splitlines() if l.strip()]


def main():
    files = sorted(f for f in os.listdir(PAPER) if f.endswith(".tex"))
    docs = {f: grab(os.path.join(PAPER, f)) for f in files}
    labels = {}
    for f, d in docs.items():
        for lab in d:
            labels.setdefault(lab, []).append(f)
    shared = {k: v for k, v in labels.items() if len(v) > 1}

    print(f"두 문서에 함께 있는 진술: {len(shared)}")
    print("=" * 70)
    identical = []
    for lab, fs in shared.items():
        a, b = docs[fs[0]][lab], docs[fs[1]][lab]
        same = norm(a[2]) == norm(b[2])
        if same:
            identical.append(lab)
        print(f"\n### {lab}   [{fs[0]} / {fs[1]}]   "
              f"{'동일' if same else '다름'}")
        if a[1] != b[1]:
            print(f"    optional arg: {a[1]!r}  vs  {b[1]!r}")
        if same:
            continue
        for line in difflib.unified_diff(
                norm(a[2]), norm(b[2]),
                fromfile=fs[0], tofile=fs[1], lineterm="", n=1):
            print("   ", line)
    print()
    print(f"본문이 완전히 동일한 것: {len(identical)} / {len(shared)}")
    if identical:
        print("   ", ", ".join(identical))
    return 0


if __name__ == "__main__":
    sys.exit(main())
