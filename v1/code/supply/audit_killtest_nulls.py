# -*- coding: utf-8 -*-
"""
What null does each kill-test actually use? (increment 310)

WHY. `OPEN_QUESTIONS.md` Register B, written at increment 308, puts
eleven closures in one row:

    Forge K1-K4, R1-R4; Construction C1, C2, C4
    "all measured no signal against permutation nulls, which hazard 7
     invalidates"

That is a claim about eleven files, and I wrote it without opening
them. Opening two was enough to show it is wrong: `e1_constr_c1.py`
says in its own header "Null: 8 draws of random signs on the real
support" -- a coin control, the very thing hazard 7 asks for -- and
`e1_forge_r1.py` compares against random-frequency templates, not a
permutation of anything.

So before any of these closures can be re-audited, the register has to
be corrected, and the way to correct it is to read the nulls
mechanically rather than to assert them again.

WHAT THIS CLASSIFIES. For each kill-test source, which of these
appears:

  COIN        random signs on the real support -- rng.choice([-1,1]),
              integers(0,2)*2-1, and similar. Varies one input and
              leaves the rest byte-identical. Hazard 7 satisfied.
  PERMUTE     rng.permutation / rng.shuffle of a data array. Destroys
              every structure at once. Hazard 7's target.
  SURROGATE   phase randomisation, random frequencies, random
              templates -- a null built from a different object rather
              than from a rearrangement of this one.
  ANALYTIC    a z-score against a stated distribution with no
              resampling at all.
  NONE        no null detected.

A file may carry more than one. The classification is by AST and by
regex over the source, and it is deliberately crude: the point is to
replace an unchecked assertion about eleven files with something
checkable, not to adjudicate the statistics.

WHAT THIS DOES NOT DO. It does not say whether a null is the RIGHT
null. A coin control on the wrong statistic is still wrong, and an
analytic z against a correct distribution is fine. The verdict per
closure is the next increment's work; this one fixes the register.

SELF-TEST runs first on synthetic sources.
"""
import io
import os
import re
import sys

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

HERE = os.path.dirname(os.path.abspath(__file__))

FILES = [
    ("Forge K1", "e1_forge_kt1.py"),
    ("Forge K2", "e1_forge_kt2.py"),
    ("Forge K3", "e1_forge_kt3.py"),
    ("Forge K4", "e1_forge_kt4.py"),
    ("Forge R1", "e1_forge_r1.py"),
    ("Forge R2", "e1_forge_r2.py"),
    ("Forge R2b", "e1_forge_r2b.py"),
    ("Forge R4", "e1_forge_r4.py"),
    ("Forge R4b", "e1_forge_r4b.py"),
    ("Constr C1", "e1_constr_c1.py"),
    ("Constr C2", "e1_constr_c2.py"),
    ("Constr C2b", "e1_constr_c2b.py"),
    ("Constr C4", "e1_constr_c4.py"),
]

PATS = [
    ("COIN", re.compile(
        r"choice\(\s*\[\s*-1(\.0)?\s*,\s*1(\.0)?\s*\]|"
        r"integers\(\s*0\s*,\s*2\s*\)\s*\*\s*2|"
        r"random_sign|random signs|randint\(\s*0\s*,\s*2\s*\)\s*\*\s*2")),
    ("PERMUTE", re.compile(r"\.permutation\(|\.shuffle\(|permuted")),
    ("SURROGATE", re.compile(
        r"phase[_ ]random|random[_ ]freq|uniform\(\s*\d+\s*,\s*\d+\s*,"
        r"\s*size|random template|surrogate")),
    ("ANALYTIC", re.compile(r"\bz\s*=\s*\(|zscore|z_score|sqrt\(n\)")),
]


def classify(src):
    found = []
    for name, rx in PATS:
        if rx.search(src):
            found.append(name)
    return found or ["NONE"]


SELFTEST = [
    ("coin", "s = rng.choice([-1.0, 1.0], size=t.shape)", "COIN"),
    ("permute", "y = rng.permutation(x)", "PERMUTE"),
    ("surrogate", "freqs = rng.uniform(10, 105, size=len(Z))",
     "SURROGATE"),
    ("none", "print('hello')", "NONE"),
]


def selftest():
    print("SELF-TEST")
    ok = True
    for name, src, want in SELFTEST:
        got = classify(src)
        good = want in got
        ok &= good
        print(f"    {name:<12} -> {'+'.join(got):<22} "
              f"{'as expected' if good else 'WRONG'}")
    print(f"    {'SELF-TEST OK' if ok else 'SELF-TEST FAILED'}\n")
    return ok


def main():
    if not selftest():
        print("DONE (self-test failed)")
        sys.exit(1)
    print("(1) the null each kill-test actually uses")
    print(f"{'closure':<12} {'file':<22} {'null':<28} evidence")
    tally = {}
    missing = []
    for label, fn in FILES:
        p = os.path.join(HERE, fn)
        if not os.path.exists(p):
            missing.append((label, fn))
            print(f"{label:<12} {fn:<22} {'FILE NOT FOUND':<28}")
            continue
        src = io.open(p, encoding="utf-8").read()
        kinds = classify(src)
        for k in kinds:
            tally[k] = tally.get(k, 0) + 1
        ev = ""
        for name, rx in PATS:
            m = rx.search(src)
            if m and name in kinds:
                ln = src[:m.start()].count("\n") + 1
                ev = f"line {ln}: {m.group(0)[:34]}"
                break
        print(f"{label:<12} {fn:<22} {'+'.join(kinds):<28} {ev}")

    print(f"\n(2) tally over {len(FILES)} files")
    for k in ("COIN", "PERMUTE", "SURROGATE", "ANALYTIC", "NONE"):
        if k in tally:
            print(f"    {k:<10} {tally[k]}")

    nperm = tally.get("PERMUTE", 0)
    ncoin = tally.get("COIN", 0)
    print(f"\n    Register B (increment 308) asserted that all of these")
    print(f"    'measured no signal against permutation nulls'.")
    print(f"    Files actually containing a permutation null: {nperm}")
    print(f"    Files actually containing a coin null:        {ncoin}")
    if missing:
        print(f"    Files named in the register that do not exist: "
              f"{len(missing)}")
    if nperm == len(FILES):
        v = ("the register's claim holds as written: every kill-test "
             "uses a permutation null")
    elif nperm == 0:
        v = ("the register's claim is FALSE for every file: not one of "
             "these kill-tests uses a permutation null, so hazard 7's "
             "objection does not apply to any of them as stated")
    else:
        v = (f"the register's claim is false as written: {nperm} of "
             f"{len(FILES)} use a permutation null and {ncoin} already "
             f"use the coin control hazard 7 asks for. The row must be "
             f"split before any of these closures is judged")
    print(f"    {v}")
    print("DONE")


if __name__ == "__main__":
    main()
