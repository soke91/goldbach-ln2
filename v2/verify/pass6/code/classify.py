# -*- coding: utf-8 -*-
"""Triage all 251 source remarks into H / S / M per the pre-registered rule.

The rule (PREREG.md sec.A) is applied sentence-by-sentence.  This script
does the mechanical half: it tags each remark with the pre-registered
cue phrases, extracts its distinctive numeric literals, and reports how
many of those literals appear anywhere in the projection.  The
adjudication of the borderline cases is done by hand against this
output; the script never decides alone, it ranks what must be read.

A remark is a MISSING-S/M CANDIDATE when it carries S or M cues and
none of its distinctive numbers reaches the projection.
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
PROJ = ["P1-mobius-fixed-class.tex", "P2-no-go-divisor-switch.tex",
        "P3-wall-second-moment.tex", "P4-coherent-cell-floor.tex",
        "P5-negative-map.tex"]

# --- pre-registered cue phrases (PREREG.md sec.A) ---------------------
H_CUE = ["an earlier version", "earlier version", "we previously",
         "previously", "withdrawn", "withdraw", "retract", "was wrong",
         "were wrong", "the wrong one", "pre-registered", "preregistered",
         "our rule", "rule failed", "fails there", "the audit's rule",
         "declined null", "this program reported", "version 3",
         "v1 said", "v2 said", "the code was", "a bug", "the bug",
         "superseded", "my own", "the objection was mine",
         "the control was owed", "nobody had", "artifact",
         "could not have seen", "misspecified", "the error was"]
S_CUE = ["only", "not a constant", "does not", "cannot", "no longer",
         "is not", "must state", "must be stated", "the range",
         "at every scale reached", "one family", "is not closed",
         "not closed", "the escape", "survives", "does not survive",
         "convention", "we write", "one symbol", "two meanings",
         "scope", "is dead", "closed", "no-go", "below its own",
         "does not grow", "not reproduced", "is not the", "caveat",
         "expiry", "out of sample", "and it fails", "the price",
         "buys nothing", "is free", "not a proof", "measurement and not",
         "not thereby"]
M_CUE = [r"\breads\b", r"\bmeasured\b", r"\bfitted\b", r"\bslope\b",
         r"\bcorrelation\b", r"\bmedian\b", r"\bsd\b",
         r"\bstandard deviation\b", r"\bdraws\b", r"\boctave", r"\brung",
         r"\bat \$N", r"\bover \$N", r"\bfactor \$", r"\bratio\b",
         r"\bpercentile\b", r"\bband\b", r"\bnull\b"]

NUM = re.compile(r"(?<![\w.])(\d+\.\d{3,})")
EV = re.compile(r"<!--\s*evidence:\s*(.*?)\s*-->")


def rd(p):
    return io.open(p, encoding="utf-8", errors="replace").read()


src = json.load(io.open(os.path.join(R, "inventory_source.json"),
                        encoding="utf-8"))
corresp = {c["src_key"]: c for c in
           json.load(io.open(os.path.join(R, "corresp.json"),
                             encoding="utf-8")) if c["src_key"]}
projtxt = "\n".join(rd(os.path.join(ROOT, "deploy", "papers", f))
                    for f in PROJ)
projres = []
for d, _, fs in os.walk(os.path.join(ROOT, "deploy", "results")):
    for f in sorted(fs):
        projres.append(rd(os.path.join(d, f)))
projres = "\n".join(projres)

rows = []
for x in src:
    if x["kind"] != "Remark":
        continue
    body = x["body"]
    low = body.lower()
    h = sorted({c for c in H_CUE if c in low})
    s = sorted({c for c in S_CUE if c in low})
    m = sorted({c.strip("\\b") for c in M_CUE if re.search(c, low)})
    nums = sorted(set(NUM.findall(body)))
    intex = [n for n in nums if n in projtxt]
    inres = [n for n in nums if n in projres]
    ev = EV.search(body)
    cls = ("H" if h else "") + ("S" if s else "") + ("M" if (m or nums)
                                                     else "")
    rows.append({
        "key": x["key"] or "(none)", "file": os.path.basename(x["file"]),
        "line": x["line"], "title": x["title"],
        "cls": cls or "?", "h": h, "s": s, "m": m,
        "nums": nums, "in_tex": intex, "in_res": inres,
        "evidence": ev.group(1) if ev else "",
        "proj": corresp.get(x["key"], {}).get("proj_key"),
        "len": len(body),
    })

print("STATISTIC: pre-registered H/S/M cue tags per source remark, its "
      "distinctive numeric literals, and how many of those reach the "
      "projection text or deploy/results/")
print("FIELD:     all 251 remarks of v2/paper/*.md")
print("CONSTANTS: cue-phrase lists fixed in PREREG.md sec.A and printed "
      "below; numeric literal = >=3 decimal places")
print("NULL:      none -- a triage that ranks remarks for hand reading; "
      "no verdict is taken from the tags alone")
print("DENOM:     251 remarks")
print()
print("H cues (%d): %s" % (len(H_CUE), ", ".join(H_CUE)))
print("S cues (%d): %s" % (len(S_CUE), ", ".join(S_CUE)))
print("M cues (%d): %s" % (len(M_CUE), ", ".join(M_CUE)))
print()

from collections import Counter
print("=" * 74)
print("A. TAG DISTRIBUTION (mechanical, before hand adjudication)")
print("=" * 74)
for k, v in sorted(Counter(r["cls"] for r in rows).items()):
    print("   %-6s %d" % (k, v))
print("   total  %d" % len(rows))
print()
print("   carries H cues            : %d" % sum(1 for r in rows if r["h"]))
print("   carries S cues            : %d" % sum(1 for r in rows if r["s"]))
print("   carries M cues or numbers : %d"
      % sum(1 for r in rows if r["m"] or r["nums"]))
print("   H and S both (mixed)      : %d"
      % sum(1 for r in rows if r["h"] and r["s"]))
print("   survives as a projection environment: %d"
      % sum(1 for r in rows if r["proj"]))

print()
print("=" * 74)
print("B. RANKED CANDIDATES -- S or M cues, no number reaching the "
      "projection, no surviving environment")
print("=" * 74)
cand = [r for r in rows
        if (r["s"] or r["nums"]) and not r["in_tex"] and not r["proj"]]
cand.sort(key=lambda r: (-len(r["nums"]), -r["len"]))
print("  %d of %d remarks" % (len(cand), len(rows)))
print()
for r in cand:
    print("  %-9s %-5d %-26s %-4s n=%-3d res=%-3d  %s"
          % (r["file"][:9], r["line"], r["key"], r["cls"],
             len(r["nums"]), len(r["in_res"]), r["title"][:44]))
    if r["evidence"]:
        print("        evidence: %s" % r["evidence"])

print()
print("=" * 74)
print("C. REMARKS WHOSE NUMBERS DO REACH THE PROJECTION TEXT")
print("=" * 74)
reach = [r for r in rows if r["in_tex"]]
for r in sorted(reach, key=lambda r: -len(r["in_tex"])):
    print("  %-9s %-5d %-26s %2d/%-2d nums in tex  %s"
          % (r["file"][:9], r["line"], r["key"], len(r["in_tex"]),
             len(r["nums"]), r["title"][:40]))
print("  total: %d" % len(reach))

print()
print("=" * 74)
print("D. REMARKS WITH NO NUMBERS AND NO S CUES (pre-registered -> H)")
print("=" * 74)
pureh = [r for r in rows if not r["nums"] and not r["s"]]
for r in pureh:
    print("  %-9s %-5d %-26s %s"
          % (r["file"][:9], r["line"], r["key"], r["title"][:48]))
print("  total: %d" % len(pureh))

json.dump(rows, io.open(os.path.join(R, "classify.json"), "w",
                        encoding="utf-8", newline="\n"),
          ensure_ascii=False, indent=1)
