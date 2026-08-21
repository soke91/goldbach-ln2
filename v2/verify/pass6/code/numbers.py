# -*- coding: utf-8 -*-
"""Two-directional numeric comparison between projection and source.

Direction 2 (contamination) is the one that matters: every number in
deploy/papers/*.tex with >=3 decimal places is looked up in the source
and in v2/results/ + deploy/results/.  A number that appears nowhere in
the source, or appears only inside a sentence that withdraws it, is a
candidate defect.

"Withdrawal context" is detected by keyword in the source sentence that
carries the number; the keyword list is fixed here and printed, so the
screen is auditable and can be re-run.
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

SOURCE = ["v2/paper/theorem_A.md", "v2/paper/wall_v3.md"]
PROJ = ["deploy/papers/P1-mobius-fixed-class.tex",
        "deploy/papers/P2-no-go-divisor-switch.tex",
        "deploy/papers/P3-wall-second-moment.tex",
        "deploy/papers/P4-coherent-cell-floor.tex",
        "deploy/papers/P5-negative-map.tex"]

WITHDRAW = ["withdraw", "withdrawn", "retract", "was wrong", "were wrong",
            "superseded", "no longer", "earlier version", "an earlier",
            "previously", "we previously", "the bug", "buggy", "a bug",
            "incorrect", "in error", "mistaken", "misspecified",
            "this was an artifact", "artifact", "corrected below",
            "since been", "obsolete", "does not survive", "not reproduced",
            "declined", "stale", "replaced by", "supersedes"]

NUM = re.compile(r"(?<![\w.])(\d+\.\d{3,})")


def rd(rel):
    with io.open(os.path.join(ROOT, rel.replace("/", os.sep)),
                 encoding="utf-8", errors="replace") as f:
        return f.read()


def walk(rel):
    base = os.path.join(ROOT, rel.replace("/", os.sep))
    for d, _, fs in os.walk(base):
        if "__pycache__" in d:
            continue
        for f in sorted(fs):
            yield os.path.join(d, f)


srctxt = {rel: rd(rel) for rel in SOURCE}
srcall = "\n".join(srctxt.values())

resblob = {}
for tree in ["v2/results", "deploy/results", "v2/code", "deploy/code"]:
    buf = []
    for p in walk(tree):
        try:
            buf.append(io.open(p, encoding="utf-8", errors="replace").read())
        except Exception:
            pass
    resblob[tree] = "\n".join(buf)


def sentences_with(txt, tok):
    out = []
    for m in re.finditer(re.escape(tok), txt):
        a = max(0, txt.rfind("\n\n", 0, m.start()))
        b = txt.find("\n\n", m.end())
        b = len(txt) if b < 0 else b
        out.append(txt[a:b].strip())
    return out


print("STATISTIC: for each numeric literal with >=3 decimals in the "
      "projection, its presence in the source and in the evidence trees, "
      "and whether the source paragraph carrying it withdraws it")
print("FIELD:     deploy/papers/*.tex (5 files) against v2/paper/*.md, "
      "v2/results/, deploy/results/, v2/code/, deploy/code/")
print("CONSTANTS: decimal-place threshold 3; withdrawal keyword list of "
      "%d phrases, fixed before the run and printed below" % len(WITHDRAW))
print("NULL:      a faithful projection quotes only numbers that the "
      "source states as current and that the evidence trees reproduce; "
      "the null is 0 unsourced and 0 withdrawn-but-quoted")
print("DENOM:     every distinct such literal in the projection")
print()
print("withdrawal keywords: " + ", ".join(WITHDRAW))
print()

print("=" * 74)
print("DIRECTION 2 -- CONTAMINATION (projection number -> source)")
print("=" * 74)
tot = 0
nosrc = []
withdrawn = []
seen = {}
for rel in PROJ:
    txt = rd(rel)
    for m in NUM.finditer(txt):
        tok = m.group(1)
        ln = txt[:m.start()].count("\n") + 1
        seen.setdefault(tok, []).append((os.path.basename(rel), ln))

for tok in sorted(seen, key=lambda t: (-len(seen[t]), t)):
    tot += 1
    where = seen[tok]
    ins = tok in srcall
    inres = [t for t in resblob if tok in resblob[t]]
    flag = []
    if not ins:
        flag.append("NOT-IN-SOURCE")
    if not inres:
        flag.append("NOT-IN-EVIDENCE")
    ctx = sentences_with(srcall, tok) if ins else []
    wd = []
    for c in ctx:
        low = c.lower()
        hits = [w for w in WITHDRAW if w in low]
        if hits:
            wd.append(hits)
    if wd and len(wd) == len(ctx):
        flag.append("ALL-SOURCE-CONTEXTS-FLAGGED")
    elif wd:
        flag.append("SOME-CONTEXT-FLAGGED(%d/%d)" % (len(wd), len(ctx)))
    line = "  %-14s x%-2d %-34s src=%-3s res=%-22s %s" % (
        tok, len(where),
        ",".join("%s:%d" % w for w in where[:2])[:34],
        "Y" if ins else "N",
        ",".join(t.split("/")[-1] for t in inres)[:22] or "-",
        " ".join(flag))
    print(line)
    if not ins:
        nosrc.append((tok, where, inres))
    if wd:
        withdrawn.append((tok, where, wd, len(ctx)))

print()
print("  distinct literals in projection : %d" % tot)
print("  absent from source              : %d" % len(nosrc))
print("  source context carries a "
      "withdrawal keyword: %d" % len(withdrawn))

print()
print("=" * 74)
print("DETAIL -- literals absent from the source")
print("=" * 74)
for tok, where, inres in nosrc:
    print("  %s  at %s" % (tok, ", ".join("%s:%d" % w for w in where)))
    print("      evidence trees containing it: %s" % (inres or "NONE"))
    for rel in PROJ:
        t = rd(rel)
        for m in re.finditer(re.escape(tok), t):
            a = t.rfind("\n\n", 0, m.start())
            b = t.find("\n\n", m.end())
            print("      | " + re.sub(
                r"\s+", " ", t[max(0, a):b if b > 0 else len(t)])[:300])
            break

print()
print("=" * 74)
print("DETAIL -- literals whose source paragraph carries a withdrawal word")
print("=" * 74)
for tok, where, wd, n in withdrawn:
    print("  %s (%d of %d source paragraphs flagged) at %s"
          % (tok, len(wd), n, ", ".join("%s:%d" % w for w in where)))
    print("      keywords: %s" % "; ".join(",".join(h) for h in wd[:4]))

json.dump({"total": tot,
           "absent_from_source": [t for t, _, _ in nosrc],
           "withdrawal_flagged": [t for t, _, _, _ in withdrawn]},
          io.open(os.path.join(R, "numbers.json"), "w",
                  encoding="utf-8", newline="\n"),
          ensure_ascii=False, indent=1)
