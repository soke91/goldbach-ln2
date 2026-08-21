# -*- coding: utf-8 -*-
r"""BRIEF sec.2.2 item 2 / PREREG sec.C(c): is any statement text different?

For each numbered statement carried into the projection, normalise both
sides to a comparable token stream -- strip markup, unify the notation
the projection deliberately renamed (A -> \AAA, E_mu -> \Emu), drop
whitespace -- and report the token-level difference.

The normalisation is a fixed list, printed below.  It cannot hide a
dropped hypothesis or a flipped inequality: those are token changes it
does not touch.  Everything the diff reports is then read by hand.
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

ROOT = r"z:\업무\goldbach-ln2-real"
R = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                 "..", "results"))
load = lambda n: json.load(io.open(os.path.join(R, n), encoding="utf-8"))

src = {x["key"]: x for x in load("inventory_source.json") if x["key"]}
prj = {x["key"]: x for x in load("inventory_proj.json") if x["key"]}
corr = {c["src_key"]: c for c in load("corresp.json") if c["src_key"]}

# fixed normalisation list, printed with the output
SUBS = [
    (r"\\begin\{(equation|align|equation\*)\}", " "),
    (r"\\end\{(equation|align|equation\*)\}", " "),
    (r"\\begin\{[a-z]+\}(\[[^\]]*\])?", " "),
    (r"\\end\{[a-z]+\}", " "),
    (r"\\label\{[^}]*\}", " "),
    (r"\\(ref|eqref|texorpdfstring|emph|textbf|text)\b", " "),
    (r"\\cite\[[^\]]*\]\{[^}]*\}", " "),
    (r"\\cite\{[^}]*\}", " "),
    (r"<!--.*?-->", " "),
    (r"\*\*", " "), (r"\*", " "),
    (r"\[(thm|cor|prop|lem|conj|eq|rem|sec|note|meas|obs):[A-Za-z0-9]*\]",
     " REF "),
    (r"\\(AAA|Emu|SS)\b", lambda m: {"AAA": "A", "Emu": "EMU",
                                     "SS": "S"}[m.group(1)]),
    (r"\bA\(N\)", "A(N)"),
    (r"\\!+", " "), (r"\\,", " "), (r"\\;", " "), (r"\\ ", " "),
    (r"\\bigl|\\bigr|\\Bigl|\\Bigr|\\left|\\right", " "),
    (r"[{}$~]", " "),
    (r"—|--|---", "-"),
    (r"\s+", " "),
]


def norm(s):
    for pat, rep in SUBS:
        s = re.sub(pat, rep, s, flags=re.S)
    return s.strip()


print("STATISTIC: token-level similarity between each source statement "
      "and its projection counterpart, after a fixed normalisation")
print("FIELD:     the 38 numbered source statements and their "
      "projection counterparts")
print("CONSTANTS: the %d normalisation rules printed below; report "
      "threshold: every statement is printed with its ratio, and every "
      "one below 0.92 is dumped in full" % len(SUBS))
print("NULL:      a faithful projection restates each result with the "
      "same hypotheses, the same quantifiers and the same inequality "
      "directions; the normalisation touches none of those")
print("DENOM:     38 statements")
print()
print("normalisation rules (pattern -> replacement):")
for pat, rep in SUBS:
    print("   %-64s -> %r" % (pat, rep if isinstance(rep, str) else "<fn>"))
print()

rows = []
for k, x in src.items():
    if x["kind"] == "Remark":
        continue
    c = corr.get(k)
    if not c or not c["proj_key"]:
        continue
    p = prj.get(c["proj_key"])
    if not p:
        continue
    a, b = norm(x["body"]), norm(p["body"])
    # compare only up to the proof, on both sides
    a = re.split(r"\bProof\b|\bproof\b", a)[0]
    b = re.split(r"\bProof\b|\bproof\b", b)[0]
    rows.append((SequenceMatcher(None, a, b).ratio(), k, c["proj_key"],
                 x["kind"], p["env"], a, b))

rows.sort()
print("=" * 74)
print("A. SIMILARITY, WORST FIRST")
print("=" * 74)
for r, k, pk, kind, env, a, b in rows:
    print("  %.3f  %-18s -> %-18s %s/%s" % (r, k, pk, kind, env))

print()
print("=" * 74)
print("B. FULL TEXT OF EVERY STATEMENT BELOW 0.92")
print("=" * 74)
for r, k, pk, kind, env, a, b in rows:
    if r >= 0.92:
        continue
    print()
    print("-" * 70)
    print("  %s  (%.3f)   %s -> %s" % (k, r, kind, env))
    print("-" * 70)
    print("  SOURCE    : %s" % a[:1400])
    print()
    print("  PROJECTION: %s" % b[:1400])
