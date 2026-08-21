# -*- coding: utf-8 -*-
"""Cross-check the numbered statements: source anchors vs projection labels.

Reads the inventories written by extract.py.  Reports, for every numbered
non-remark statement in the source, whether the projection carries the same
label, in which paper, and under which environment (i.e. whether the grade
moved).  Also reports projection labels with no source anchor.
"""
import io
import json
import os
import sys

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

R = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                 "..", "results"))


def load(n):
    return json.load(io.open(os.path.join(R, n), encoding="utf-8"))


src = load("inventory_source.json")
prj = load("inventory_proj.json")
plab = load("labels_proj.json")
slab = load("labels_source.json")

GRADE = {"Theorem": 5, "Corollary": 4, "Proposition": 4, "Lemma": 3,
         "Conjecture": 2, "Definition": 1, "Remark": 0,
         "theorem": 5, "corollary": 4, "proposition": 4, "lemma": 3,
         "conjecture": 2, "definition": 1,
         "note": 0, "observation": 1, "measurement": 1, "remark": 0}

pby = {}
for x in prj:
    if x["key"]:
        pby[x["key"]] = x

print("STATISTIC: presence and grade of every numbered source statement "
      "in the projection")
print("FIELD:     the 38 non-Remark numbered statements of "
      "v2/paper/*.md, and all 77 numbered environments of "
      "deploy/papers/*.tex")
print("CONSTANTS: grade ladder Theorem=5 Corollary=Proposition=4 Lemma=3 "
      "Conjecture=2 Definition=Observation=Measurement=1 Remark=Note=0")
print("NULL:      a faithful projection carries every source label at the "
      "same grade; deviations are the signal")
print("DENOM:     38 source statements; 77 projection environments")
print()

named = [x for x in src if x["kind"] != "Remark"]
print("=" * 72)
print("A. SOURCE NUMBERED STATEMENTS (%d)" % len(named))
print("=" * 72)
miss, regrade, ok = [], [], []
for x in sorted(named, key=lambda y: (y["file"], y["line"])):
    k = x["key"]
    if not k:
        print("  %-9s %-12s %-4d  (NO ANCHOR IN SOURCE) %s"
              % (os.path.basename(x["file"])[:9], x["kind"], x["line"],
                 x["title"][:40]))
        continue
    p = pby.get(k)
    if p is None:
        where = "label used elsewhere in proj" if k in plab else "ABSENT"
        miss.append((k, x, where))
        print("  %-9s %-12s %-4d  %-22s -> %s"
              % (os.path.basename(x["file"])[:9], x["kind"], x["line"],
                 k, where))
    else:
        gs, gp = GRADE[x["kind"]], GRADE[p["env"]]
        tag = "same" if gs == gp else ("DOWN" if gp < gs else "UP")
        if tag != "same":
            regrade.append((k, x, p, tag))
        else:
            ok.append(k)
        print("  %-9s %-12s %-4d  %-22s -> %-11s %-28s %s"
              % (os.path.basename(x["file"])[:9], x["kind"], x["line"],
                 k, p["env"], os.path.basename(p["file"])[:28], tag))

print()
print("  present at same grade : %d" % len(ok))
print("  regraded              : %d" % len(regrade))
print("  absent from projection: %d" % len(miss))

print()
print("=" * 72)
print("B. REGRADES (detail)")
print("=" * 72)
for k, x, p, tag in regrade:
    print("  %-22s %s(%s) -> %s(%s)   [%s]"
          % (k, x["kind"], os.path.basename(x["file"]), p["env"],
             os.path.basename(p["file"]), tag))
    print("      source title: %s" % x["title"][:66])
    print("      proj   title: %s" % p["title"][:66])

print()
print("=" * 72)
print("C. PROJECTION ENVIRONMENTS WITH NO SOURCE ANCHOR OF THAT NAME")
print("=" * 72)
n = 0
for x in sorted(prj, key=lambda y: (y["file"], y["line"])):
    if x["key"] and x["key"] not in slab:
        n += 1
        print("  %-30s %-5d %-12s %-24s %s"
              % (os.path.basename(x["file"])[:30], x["line"], x["env"],
                 x["key"], x["title"][:34]))
print("  total: %d" % n)

print()
print("=" * 72)
print("D. SOURCE REMARK ANCHORS PRESENT IN THE PROJECTION AS LABELS")
print("=" * 72)
rem = [x for x in src if x["kind"] == "Remark" and x["key"]]
hit = [x for x in rem if x["key"] in plab]
print("  source remarks carrying an anchor: %d of %d"
      % (len(rem), sum(1 for x in src if x["kind"] == "Remark")))
print("  of those, anchor survives as a projection label: %d" % len(hit))
for x in sorted(hit, key=lambda y: y["key"]):
    p = pby.get(x["key"])
    print("    %-26s -> %-12s %s"
          % (x["key"], p["env"] if p else "(ref only)",
             os.path.basename(p["file"])[:28] if p else
             ",".join(os.path.basename(f) for f in plab[x["key"]])[:28]))
