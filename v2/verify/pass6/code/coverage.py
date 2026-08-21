# -*- coding: utf-8 -*-
"""Which dropped remarks could the projection possibly depend on?

PREREG.md sec.B fixes the line: a remark whose subject the projection
never raises is a legitimate exclusion (L3, not counted); one whose
subject the projection does raise, minus the remark's qualification, is
L1 or L2.

The mechanical proxy for "the projection raises the subject" used here
is the remark's own evidence pointer: each source remark carries
<!-- evidence: SCRIPT.py -->.  deploy/code/ holds 34 scripts, a subset
of v2/code/'s 235.  A remark whose evidence script was NOT deployed
belongs to a thread the projection did not carry at all; a remark whose
evidence script WAS deployed sits on measurements the projection is
quoting, and its qualification is then load-bearing.

This is a screen, not a verdict: the shortlist it produces is read by
hand.
"""
import io
import json
import os
import re
import sys

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

ROOT = r"z:\업무\goldbach-ln2-real"
R = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                 "..", "results"))
load = lambda n: json.load(io.open(os.path.join(R, n), encoding="utf-8"))

src = load("inventory_source.json")
corr = {c["src_key"]: c for c in load("corresp.json") if c["src_key"]}
cls = {c["key"]: c for c in load("classify.json")}

deployed = set(os.listdir(os.path.join(ROOT, "deploy", "code")))
v2code = set(os.listdir(os.path.join(ROOT, "v2", "code")))

EV = re.compile(r"<!--\s*evidence:\s*(.*?)\s*-->")

rows = []
for x in src:
    if x["kind"] != "Remark":
        continue
    evs = EV.findall(x["body"])
    scripts = [e.strip() for e in evs if e.strip().endswith(".py")]
    scripts = [os.path.basename(s) for s in scripts]
    dep = [s for s in scripts if s in deployed]
    kept = bool(corr.get(x["key"], {}).get("proj_key"))
    nums = cls.get(x["key"] or "(none)", {}).get("nums", [])
    intex = cls.get(x["key"] or "(none)", {}).get("in_tex", [])
    rows.append({"key": x["key"] or "(none)",
                 "file": os.path.basename(x["file"]), "line": x["line"],
                 "title": x["title"], "scripts": scripts, "dep": dep,
                 "kept": kept, "nums": nums, "intex": intex,
                 "analytic": (not scripts)})

print("STATISTIC: for each source remark, whether its evidence script "
      "was carried into deploy/code/, whether the remark itself "
      "survived as a projection environment, and whether any of its "
      "numbers reach the projection text")
print("FIELD:     251 source remarks; deploy/code/ (34 scripts) against "
      "v2/code/ (235)")
print("CONSTANTS: none; the evidence pointer is the remark's own "
      "<!-- evidence: --> comment")
print("NULL:      a remark standing on a deployed script, dropped, and "
      "with no number reaching the projection, is the shortlist this "
      "pass reads by hand; the rest are threads the projection never "
      "raised")
print("DENOM:     251 remarks")
print()

n_dep = [r for r in rows if r["dep"]]
n_nodep = [r for r in rows if r["scripts"] and not r["dep"]]
n_anal = [r for r in rows if r["analytic"]]

print("=" * 74)
print("A. PARTITION")
print("=" * 74)
print("  remarks total                                   : %d" % len(rows))
print("  standing on a script that WAS deployed          : %d" % len(n_dep))
print("  standing on a script that was NOT deployed      : %d"
      % len(n_nodep))
print("  carrying no evidence pointer (analytic/prose)   : %d"
      % len(n_anal))
print()
print("  of the deployed-script remarks:")
print("     survive as a projection environment          : %d"
      % sum(1 for r in n_dep if r["kept"]))
print("     dropped, but some number reaches the tex     : %d"
      % sum(1 for r in n_dep if not r["kept"] and r["intex"]))
print("     dropped, and NO number reaches the tex       : %d"
      % sum(1 for r in n_dep if not r["kept"] and not r["intex"]))

print()
print("=" * 74)
print("B. SHORTLIST -- deployed script, remark dropped, no number in the "
      "projection")
print("=" * 74)
short = [r for r in n_dep if not r["kept"] and not r["intex"]]
short.sort(key=lambda r: (r["file"], r["line"]))
for r in short:
    print("  %-9s %-5d %-24s %s" % (r["file"][:9], r["line"], r["key"],
                                    r["title"][:46]))
    print("        evidence: %s" % ", ".join(r["dep"]))
print("  count: %d" % len(short))

print()
print("=" * 74)
print("C. DEPLOYED-SCRIPT REMARKS THAT DID REACH THE PROJECTION")
print("=" * 74)
for r in n_dep:
    if r["kept"] or r["intex"]:
        print("  %-24s %-8s %s"
              % (r["key"], "env" if r["kept"] else "nums",
                 r["title"][:44]))

print()
print("=" * 74)
print("D. SCRIPTS CITED BY SOURCE REMARKS BUT ABSENT FROM v2/code/")
print("=" * 74)
missing = sorted({s for r in rows for s in r["scripts"]
                  if s not in v2code})
for s in missing:
    who = [r["key"] for r in rows if s in r["scripts"]]
    print("  %-40s cited by %s" % (s, ", ".join(who[:4])))
print("  count: %d" % len(missing))
