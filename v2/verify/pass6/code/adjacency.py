# -*- coding: utf-8 -*-
"""BRIEF sec.2.2 item 4: did a surviving statement lose its qualifying remark?

For every source object that the projection carries, list the source
remarks that sit between it and the next numbered statement -- the
remarks that were narrowing its scope in situ -- and report which of
those reached the projection.

A surviving statement whose in-situ remarks all vanished is where a
claim can silently widen.  Those are printed first and read by hand.
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
load = lambda n: json.load(io.open(os.path.join(R, n), encoding="utf-8"))

src = load("inventory_source.json")
corr = {c["src_key"]: c for c in load("corresp.json") if c["src_key"]}

by = {}
for x in src:
    by.setdefault(x["file"], []).append(x)
for v in by.values():
    v.sort(key=lambda x: x["line"])

print("STATISTIC: for each source statement carried into the projection, "
      "the source remarks attached to it in situ and whether each "
      "reached the projection")
print("FIELD:     v2/paper/*.md ordered by line; 'attached' = every "
      "Remark between the statement and the next non-Remark statement")
print("CONSTANTS: none")
print("NULL:      a faithful projection carries a statement together "
      "with the remarks that narrow it; a survivor with all attached "
      "remarks dropped is the signal")
print("DENOM:     the 38 numbered source statements")
print()

flag = []
for f, xs in by.items():
    for i, x in enumerate(xs):
        if x["kind"] == "Remark":
            continue
        att = []
        for y in xs[i + 1:]:
            if y["kind"] != "Remark":
                break
            att.append(y)
        c = corr.get(x["key"])
        if not c or not c["proj_key"]:
            continue
        kept = [y for y in att
                if corr.get(y["key"], {}).get("proj_key")]
        flag.append((x, c, att, kept))

flag.sort(key=lambda t: (len(t[3]), -len(t[2])))

print("=" * 74)
print("A. SURVIVING STATEMENTS WITH ATTACHED REMARKS, NONE OF WHICH "
      "SURVIVED")
print("=" * 74)
n = 0
for x, c, att, kept in flag:
    if att and not kept:
        n += 1
        print("  %-16s %-12s %s:%d  ->  %s in %s"
              % (x["key"], x["kind"], os.path.basename(x["file"]),
                 x["line"], c["proj_key"],
                 os.path.basename(c["proj_file"])))
        for y in att:
            print("        dropped  %-24s %s"
                  % (y["key"] or "(no anchor)", y["title"][:44]))
print("  count: %d" % n)

print()
print("=" * 74)
print("B. SURVIVING STATEMENTS, SOME ATTACHED REMARKS KEPT")
print("=" * 74)
for x, c, att, kept in flag:
    if att and kept:
        print("  %-16s -> %-16s  kept %d of %d"
              % (x["key"], c["proj_key"], len(kept), len(att)))
        for y in att:
            mark = "kept   " if corr.get(y["key"], {}).get("proj_key") \
                else "dropped"
            print("        %s  %-24s %s"
                  % (mark, y["key"] or "(no anchor)", y["title"][:42]))

print()
print("=" * 74)
print("C. SURVIVING STATEMENTS WITH NO ATTACHED REMARKS IN THE SOURCE")
print("=" * 74)
for x, c, att, kept in flag:
    if not att:
        print("  %-16s %-12s -> %s" % (x["key"], x["kind"], c["proj_key"]))
