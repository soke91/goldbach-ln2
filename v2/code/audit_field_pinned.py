# -*- coding: utf-8 -*-
r"""
Does this repository's own FIELD block pin its own field?

WHAT IS AT STAKE

verify/pass2 and verify/pass3 examined two of Conjecture conj:L's four
stamps and found the same defect in both: the number is right as far
as anything could tell, and cannot be reconstructed from what was
printed.  One was missing its grid -- 100 grids consistent with the
two facts it published give its count.  The other was missing its cell
index and its prime range -- twenty-one configurations were tried and
none reproduced its interval.

Recording that about another tree and not asking it here would be
incoherent.  The gate has G4, which requires every result file to open
with STATISTIC: and FIELD:, and pass3 found that G4 was silently
skipping the verify tree entirely -- so the one check that exists for
this had a hole in it, which is not evidence that the check works.
**And G4 only looks for the presence of the line.**  A FIELD: line
that names no numbers pins nothing.

THE TEST

A script fixes its field in module-level constants.  If a constant
that enters the computation has a value that appears nowhere in the
result file, then a reader holding only that file cannot rebuild what
was run -- exactly the position verify/pass2 and pass3 were in.  So:
for every script in code/ and in verify/*/code/, take its module-level
numeric constants, and ask whether each value appears in its own
result file.

Two kinds of constant are set aside before the run, by name, because
they cannot change a number: BLOCK and CHUNK (memory tiling) and
WIDTH, DEC, DIGITS, PREC (printing).  **SEED is deliberately NOT set
aside** -- a null whose seed is unrecorded is unreproducible, which is
the defect under test rather than an exception to it.

BACKS: Remark {#rem:fieldpinned} in paper/wall_v3.md.

PRE-REGISTERED PREDICTIONS (written before this script was run)

  F1  At least one result file in this repository has a field-bearing
      constant whose value appears nowhere in it.  The defect is not
      only the other tree's.
  F2  But most files are clean: the share with no absence is above
      0.80.  If it is not, this repository has no standing to have
      reported the defect elsewhere.
  F3  **No absence is a field bound.**  Among the absent constants,
      none has a name in the bound family -- N, K, P, MAX, MIN, LO,
      HI, CAP, RANGE, BAND, SCAN, DEPTH, SEED -- so what is missing is
      cosmetic and the grid, the range and the seed are always
      printed.  If this fails, the repository has the same defect it
      recorded in the other tree, and the failure names the files.
  F4  The five files verify/pass2 and verify/pass3 wrote are clean,
      the last commit having rewritten their headers for exactly this
      reason.

REFUTATION RULE (fixed before the run)

  F1  REFUTED if every file is clean.  That would be a better outcome
      than the prediction and would say the discipline is tighter than
      expected; it is written as a prediction so the opposite result
      cannot be presented as a surprise.
  F2  REFUTED below 0.80.  Then the defect is the rule here too, and
      the honest reading is that verify/pass2 and pass3 documented in
      another tree a fault this one shares.
  F3  REFUTED by a single absent constant in the bound family.  The
      name list is fixed above and is a heuristic: a bound can be
      named something else, and a constant in the family can be
      cosmetic in a particular script.  So a failure of F3 is a
      pointer to files to read, not a verdict on them, and the files
      are named in the output for that reason.
  F4  REFUTED by any absence in those five.  This one gates nothing
      but it is the check on the fix the last commit claimed.

  **The test is one-sided and that is stated before the run.**  A
  value that appears in the file may appear by coincidence -- 2 or 100
  will match something in almost any output -- so "clean" here means
  "no absence detected", never "reconstructible".  The failures are
  informative and the passes are weak, and no count below is presented
  as a measure of how reconstructible anything is.

  NO NULL IS RUN and none applies: this is a census of a deterministic
  property of files on disk, not a signal against a background.
"""

import io
import os
import re
import sys

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
CODE = os.path.join(ROOT, "code")
RESULTS = os.path.join(ROOT, "results")
VERIFY = os.path.join(ROOT, "verify")
OUT = os.path.join(RESULTS, "audit_field_pinned.txt")

CONST = re.compile(r"^([A-Z][A-Z0-9_]*)\s*=\s*(.+?)\s*(?:#.*)?$", re.M)
NUMTOK = re.compile(r"-?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?")
SETASIDE = ("BLOCK", "CHUNK", "WIDTH", "DEC", "DIGITS", "PREC")
BOUNDS = ("N", "K", "P", "MAX", "MIN", "LO", "HI", "CAP",
          "RANGE", "BAND", "SCAN", "DEPTH", "SEED")
CLEAN_SHARE = 0.80


def pairs():
    """(script, result) for everything the repository publishes"""
    out = []
    for f in sorted(os.listdir(CODE)):
        if f.endswith(".py"):
            out.append((os.path.join(CODE, f),
                        os.path.join(RESULTS, f[:-3] + ".txt")))
    if os.path.isdir(VERIFY):
        for base, _, fs in os.walk(VERIFY):
            if os.path.basename(base) != "code":
                continue
            for f in sorted(fs):
                if not f.endswith(".py"):
                    continue
                out.append((os.path.join(base, f),
                            os.path.join(os.path.dirname(base),
                                         "results", f[:-3] + ".txt")))
    return out


def constants(src):
    """module-level numeric constants, name -> list of values"""
    out = {}
    for m in CONST.finditer(src):
        name, rhs = m.group(1), m.group(2)
        if name in ("OUT", "ROOT", "RES", "RESULTS", "CODE", "VERIFY"):
            continue
        if any(w in name for w in SETASIDE):
            continue
        if rhs.startswith(("os.", "r\"", "r'", "\"", "'", "re.")):
            continue
        if "(" in rhs and not rhs.startswith(("(", "[")):
            continue                      # a call, not a literal
        toks = NUMTOK.findall(rhs)
        if not toks or len(toks) > 24:
            continue
        if re.search(r"[A-Za-z_]", re.sub(r"[eE]", "", rhs)):
            continue                      # names on the right: not a literal
        out[name] = [float(t) for t in toks]
    return out


def present(v, nums):
    """is v among the numerals of the result file?"""
    if v != v:
        return False
    s = ("%r" % v).rstrip("0").rstrip(".")
    dec = len(s.split(".")[1]) if "." in s else 0
    tol = 0.5 * 10.0 ** (-dec)
    return any(abs(v - x) <= tol for x in nums)


def main():
    lines = []

    def say(t=""):
        print(t)
        sys.stdout.flush()
        lines.append(t)

    rows, absent_rows = [], []
    for script, result in pairs():
        if not os.path.exists(result):
            continue
        src = io.open(script, encoding="utf-8", errors="replace").read()
        body = src.split('"""', 2)[-1] if src.count('"""') >= 2 else src
        cons = constants(body)
        if not cons:
            continue
        txt = io.open(result, encoding="utf-8", errors="replace").read()
        nums = [float(t) for t in NUMTOK.findall(txt)]
        miss = []
        for name, vals in sorted(cons.items()):
            gone = [v for v in vals if not present(v, nums)]
            if gone:
                miss.append((name, gone))
        rel = os.path.relpath(script, ROOT).replace(os.sep, "/")
        rows.append((rel, len(cons), len(miss)))
        for name, gone in miss:
            absent_rows.append((rel, name, gone))

    say("a census of %d script/result pairs that carry constants"
        % len(rows))
    say("  set aside by name before the run: %s"
        % ", ".join(SETASIDE))
    say("  the test is one-sided: an absence is informative, a match "
        "may be")
    say("  coincidence, and \"clean\" means no absence was detected")
    dec = 4
    say("PRINTBOUND audit_field_pinned %d %.8f"
        % (dec, 0.5 * 10.0 ** (-dec)))

    clean = sum(1 for _, _, m in rows if m == 0)
    share = clean / len(rows)

    # -------------------------------------------------------------- F1
    say()
    say("F1  is there any absence at all?")
    say("  %d of %d pairs have at least one absent constant, %d "
        "absences in all" % (len(rows) - clean, len(rows),
                             len(absent_rows)))
    f1 = len(absent_rows) > 0
    say("  F1 %s   (cap: at least one)" % ("hold" if f1 else "REFUTED"))
    say("COUNT fieldpinned_absences %d" % len(absent_rows))

    # -------------------------------------------------------------- F2
    say()
    say("F2  are most files clean?")
    say("  clean %d of %d, share %.4f against the cap %.2f"
        % (clean, len(rows), share, CLEAN_SHARE))
    f2 = share > CLEAN_SHARE
    say("  F2 %s   (cap: above %.2f)"
        % ("hold" if f2 else "REFUTED", CLEAN_SHARE))
    say("SHARE fieldpinned_clean %.4f" % share)

    # -------------------------------------------------------------- F3
    say()
    say("F3  is any absence a field bound?")
    bad = [(f, n, g) for f, n, g in absent_rows
           if any(w in n for w in BOUNDS)]
    say("  absences whose name is in the bound family %s:"
        % ", ".join(BOUNDS))
    if not bad:
        say("    none")
    for f, n, g in bad:
        say("    %-52s %-12s %s"
            % (f, n, ", ".join("%g" % v for v in g)))
    f3 = len(bad) == 0
    say("  F3 %s   (cap: a single one)" % ("hold" if f3 else "REFUTED"))
    say("COUNT fieldpinned_bound_absences %d" % len(bad))

    # -------------------------------------------------------------- F4
    say()
    say("F4  are the five files the last commit rewrote clean?")
    five = [(f, c, m) for f, c, m in rows if f.startswith("verify/")]
    for f, c, m in five:
        say("  %-58s %2d constants, %d absent" % (f, c, m))
    f4 = len(five) > 0 and all(m == 0 for _, _, m in five)
    say("  F4 %s   (cap: any absence)" % ("hold" if f4 else "REFUTED"))

    say()
    say("  every absence, so the files can be read rather than "
        "trusted")
    say("      %-52s %-12s %s" % ("script", "constant", "values"))
    for f, n, g in absent_rows:
        say("      %-52s %-12s %s"
            % (f, n, ", ".join("%g" % v for v in g)))
    say("SCALES 1")

    say()
    say("=" * 70)
    say("F1 %s  F2 %s  F3 %s  F4 %s"
        % tuple("hold" if v else "REFUTED"
                for v in (f1, f2, f3, f4)))
    say()
    if f2 and f3:
        say("what this repository recorded about the other tree does "
            "not describe")
        say("it: the grid, the range and the seed are printed, and "
            "what goes")
        say("unprinted is cosmetic. That is a weak statement by "
            "construction --")
        say("the test cannot certify reconstructibility, only fail to "
            "find a")
        say("hole -- and it is the statement the test supports.")
    elif not f3:
        say("this repository has the defect it recorded elsewhere. "
            "The files")
        say("are named above and the bound family is a heuristic, so "
            "each one")
        say("is a pointer to read rather than a verdict.")

    head = [
        "STATISTIC: per script, the module-level numeric constants",
        "           whose value appears nowhere in that script's own",
        "           result file -- a one-sided test for whether a",
        "           reader holding only the result can rebuild the",
        "           field. BLOCK, CHUNK, WIDTH, DEC, DIGITS and PREC",
        "           are set aside by name before the run; SEED is",
        "           deliberately not.",
        "FIELD: every script in v2/code/ and v2/verify/*/code/ that",
        "       has a result file beside it and declares at least one",
        "       module-level numeric constant, at the commit this was",
        "       run on. A constant counts as present when some",
        "       numeral in the result file matches it to half a unit",
        "       in its own last printed place.",
        "",
    ]
    io.open(OUT, "w", encoding="utf-8", newline="\n").write(
        "\n".join(head + lines) + "\n")
    print("\nwrote %s" % os.path.normpath(OUT))
    return 0


if __name__ == "__main__":
    sys.exit(main())
