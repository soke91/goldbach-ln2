# -*- coding: utf-8 -*-
"""Map every source anchor to its projection counterpart.

The projection renames as it regrades: rem:X -> note:X / meas:X / obs:X,
prop:X -> obs:X / meas:X.  So label identity alone (labels.py) understates
survival and misses regrades.  Here we match on the STEM after the colon,
and on title text, and report the grade move.
"""
import io
import json
import os
import re
import sys
from difflib import SequenceMatcher

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

R = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                 "..", "results"))
load = lambda n: json.load(io.open(os.path.join(R, n), encoding="utf-8"))

src = load("inventory_source.json")
prj = load("inventory_proj.json")

GRADE = {"Theorem": 5, "Corollary": 4, "Proposition": 4, "Lemma": 3,
         "Conjecture": 2, "Definition": 1, "Remark": 0,
         "theorem": 5, "corollary": 4, "proposition": 4, "lemma": 3,
         "conjecture": 2, "definition": 1,
         "note": 0, "observation": 1, "measurement": 1, "remark": 0}

stem = lambda k: k.split(":", 1)[1].lower() if ":" in k else k.lower()
norm = lambda s: re.sub(r"[^a-z0-9 ]", " ",
                        re.sub(r"\$[^$]*\$", " ", (s or "").lower()))
norm = lambda s, _n=norm: re.sub(r"\s+", " ", _n(s)).strip()


def sim(a, b):
    return SequenceMatcher(None, norm(a), norm(b)).ratio()


pstem = {}
for x in prj:
    if x["key"]:
        pstem.setdefault(stem(x["key"]), []).append(x)

rows = []
for x in src:
    k = x["key"]
    best, how = None, ""
    if k:
        cand = pstem.get(stem(k), [])
        if cand:
            best, how = cand[0], "stem"
    if best is None and x["title"]:
        sc = [(sim(x["title"], p["title"]), p) for p in prj if p["title"]]
        sc.sort(key=lambda t: -t[0])
        if sc and sc[0][0] >= 0.62:
            best, how = sc[0][1], "title %.2f" % sc[0][0]
    rows.append((x, best, how))

print("STATISTIC: survival and grade-move of every numbered source object "
      "in the projection, matched on label stem then on title text")
print("FIELD:     289 source objects (251 Remark + 38 numbered) vs 77 "
      "projection environments")
print("CONSTANTS: title-similarity acceptance threshold 0.62 "
      "(SequenceMatcher on lowercased, math-stripped titles)")
print("NULL:      none -- a correspondence table, not a test; the "
      "threshold is fixed before reading the output")
print("DENOM:     289 source objects")
print()

nr = [r for r in rows if r[0]["kind"] == "Remark"]
nn = [r for r in rows if r[0]["kind"] != "Remark"]

print("=" * 74)
print("A. NUMBERED (non-Remark) SOURCE STATEMENTS -- grade moves")
print("=" * 74)
for x, p, how in nn:
    if p is None:
        print("  MISSING  %-18s %-12s %s"
              % (x["key"] or "(none)", x["kind"], x["title"][:44]))
        continue
    gs, gp = GRADE[x["kind"]], GRADE[p["env"]]
    if gs != gp:
        print("  %-8s %-18s %-12s -> %-12s %-28s"
              % ("DOWN" if gp < gs else "UP", x["key"], x["kind"],
                 p["env"], os.path.basename(p["file"])[:28]))
        print("           key %s -> %s   [%s]" % (x["key"], p["key"], how))
        print("           %s" % x["title"][:62])
print()
print("  numbered statements: %d ; missing %d ; regraded %d"
      % (len(nn), sum(1 for _, p, _ in nn if p is None),
         sum(1 for x, p, _ in nn
             if p is not None and GRADE[x["kind"]] != GRADE[p["env"]])))

print()
print("=" * 74)
print("B. REMARKS THAT SURVIVE AS A PROJECTION ENVIRONMENT")
print("=" * 74)
surv = [(x, p, how) for x, p, how in nr if p is not None]
for x, p, how in sorted(surv, key=lambda t: t[1]["file"]):
    print("  %-26s -> %-12s %-28s %-4s  %s"
          % (x["key"] or "(no anchor)", p["env"],
             os.path.basename(p["file"])[:28], p["line"], how))
print("  surviving as an environment: %d of %d" % (len(surv), len(nr)))

print()
print("=" * 74)
print("C. REMARKS WITH NO ENVIRONMENT IN THE PROJECTION (%d)"
      % (len(nr) - len(surv)))
print("=" * 74)
for x, p, how in nr:
    if p is None:
        print("  %-9s %-5d %-28s %s"
              % (os.path.basename(x["file"])[:9], x["line"],
                 x["key"] or "(no anchor)", x["title"][:52]))

with io.open(os.path.join(R, "corresp.json"), "w",
             encoding="utf-8", newline="\n") as f:
    json.dump([{"src_key": x["key"], "src_kind": x["kind"],
                "src_file": x["file"], "src_line": x["line"],
                "src_title": x["title"],
                "proj_key": p["key"] if p else None,
                "proj_env": p["env"] if p else None,
                "proj_file": p["file"] if p else None,
                "proj_line": p["line"] if p else None,
                "how": how} for x, p, how in rows],
              f, ensure_ascii=False, indent=1)
