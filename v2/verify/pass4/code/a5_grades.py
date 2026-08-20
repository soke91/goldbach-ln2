# -*- coding: utf-8 -*-
"""
a5_grades.py  --  pass4, blind mathematical re-verification.

PRE-REGISTRATION (fixed before the run).

WHAT IS MEASURED.  Two purely textual audits of the five .tex sources.

 (G1) The grade audit of BRIEF A section 3.1.  For every numbered
      environment in each paper, whether a \\begin{proof} follows its
      \\end{...} before the next numbered environment starts.  A
      theorem / corollary / proposition / lemma with no proof is a
      grade violation; an observation / measurement / conjecture
      without one is not.
      FALSIFIED (no finding) if every plain-style statement has a
      proof block.

 (G2) The symbol audit of BRIEF A section 3.6, restricted to the one
      letter the papers themselves single out.  P1 declares "one
      symbol must not carry two meanings" when it renames Huang-Li's
      A to \\AAA.  This lists every distinct definitional use of the
      letter B in P1 and P2.
      FALSIFIED (no finding) if B is defined exactly once per paper.

NULL.  None applies: both are exact textual enumerations of a fixed
finite corpus with no sampling and no sign input.
"""
import io
import re
import os

# The papers under audit.  Overridable so the pass can be re-run against
# a copy held anywhere; it defaults to the deployed set in this
# repository, which is where those files live once a pass is over.
PAPERS = os.environ.get(
    "AUDIT_PAPERS_DIR",
    os.path.join(os.path.dirname(os.path.abspath(__file__)),
                 "..", "..", "..", "..", "deploy", "papers"),
)
FILES = sorted(f for f in os.listdir(PAPERS) if f.endswith(".tex"))
PLAIN = {"theorem", "corollary", "proposition", "lemma"}
LOOSE = {"observation", "measurement", "conjecture"}
ALL = PLAIN | LOOSE

print(__doc__.strip())
print()
print("G1  grade audit: every numbered statement, and whether a proof follows")
print()
for fn in FILES:
    src = io.open(os.path.join(PAPERS, fn), encoding="utf-8").read()
    toks = []
    for m in re.finditer(r"\\(begin|end)\{(%s|proof)\}" % "|".join(ALL), src):
        toks.append((m.start(), m.group(1), m.group(2)))
    print("  == %s" % fn)
    i = 0
    viol = 0
    n = 0
    while i < len(toks):
        pos, be, env = toks[i]
        if be == "begin" and env in ALL:
            n += 1
            # find its end
            j = i + 1
            while j < len(toks) and not (toks[j][1] == "end" and toks[j][2] == env):
                j += 1
            lab = ""
            mm = re.search(r"\\label\{([^}]*)\}", src[pos:toks[j][0]] if j < len(toks) else src[pos:])
            if mm:
                lab = mm.group(1)
            nxt = toks[j + 1] if j + 1 < len(toks) else None
            has = bool(nxt and nxt[1] == "begin" and nxt[2] == "proof")
            line = src[:pos].count("\n") + 1
            flag = ""
            if env in PLAIN and not has:
                flag = "   <-- NO PROOF BLOCK"
                viol += 1
            print("     %-12s %-22s line %-5d proof=%-5s%s"
                  % (env, lab, line, has, flag))
            i = j + 1
        else:
            i += 1
    print("     statements: %d   plain-style without a proof block: %d" % (n, viol))
    print()

print("G2  the letter B, definitional uses")
for fn in ("P1-mobius-fixed-class.tex", "P2-no-go-divisor-switch.tex"):
    src = io.open(os.path.join(PAPERS, fn), encoding="utf-8").read()
    print("  == %s" % fn)
    pats = [r"B\(K\)\s*=", r"B\(N\)\s*:?=", r"B_\{?\\?log\}?\(K\)\s*:?=",
            r"B_w\s*:?=", r"B_H\(N\)\s*=", r"\\tau_3\(q\)\^\{?B",
            r"with\s+\$?B=3", r"B\(s\)\s*=", r"for all \$A,B>0",
            r"B\(x\)|W_j"]
    for p in pats:
        for m in re.finditer(p, src):
            a = max(0, m.start() - 70)
            ctx = " ".join(src[a:m.end() + 70].split())
            print("     %-22s line %-5d ...%s..."
                  % (p[:22], src[:m.start()].count("\n") + 1, ctx))
    print()
