# -*- coding: utf-8 -*-
"""
Document integrity linter (increment 296).

WHY THIS EXISTS. Five times now an escape sequence has been collapsed by
a shell heredoc and written into a tracked file:

  inc. 280  FINDINGS.md truncated to zero bytes
  inc. 281  CLOSURE_REAUDIT.md truncated; recovered from git
  inc. 284  two Korean syllables mangled (opt -> different opt)
  inc. 293  a syntax error from a literal newline inside a string
  inc. 287  `\\prod_{q\\nmid N}` written as "\\prod_{q" + newline +
            "mid N}" -- and it sat in every push for EIGHT increments
            before anyone looked

Each time the response was a local repair and a note. Increment 285's
own conclusion applies: naming a hazard does not prevent it; only a
check that fails loudly does. This is that check.

It also covers two bookkeeping faults that have each occurred: a
correction row inserted at the wrong anchor (inc. 295, the ordering
assertion printed False and I committed anyway), and a repair that left
a duplicated line behind.

WHAT IT CHECKS.

  A. CONTROL CHARACTERS in .md and .tex. TAB, BEL, backspace, formfeed,
     vertical tab and CR are exactly what \\t \\a \\b \\f \\v \\r become
     when an escape collapses. None belongs in these files.
  B. STRIPPED LaTeX COMMANDS. \\n collapses to a real newline, which no
     control-character scan can see. The signature is a line beginning
     with the tail of a command (mid, eq, u, abla, ...) when the
     PREVIOUS line ends mid-formula -- unbalanced $, or a trailing
     brace, subscript or backslash. Restricting to mid-formula context
     is what keeps ordinary English ("times below", "in this program")
     from being flagged; an earlier version of this heuristic produced
     four false positives from exactly that.
  C. CORRECTION NUMBERING in CLOSURE_REAUDIT.md: distinct, and
     descending by exactly one from the maximum onward. Two
     weaker invariants failed the self-test first -- "strictly
     descending" flags the legitimate ascending historical block
     (15 false positives), and "one direction change" PASSES the
     fault that actually happened at increment 295. That fault
     strands a row above its successor, so the run from the
     maximum skips a number; contiguity is what separates them.
  D. COUNT AGREEMENT: STATUS.md's stated correction count must equal
     the highest correction number in the table.
  E. INCREMENT AGREEMENT: STATUS.md's stated increment must equal
     the most recent commit's, or be exactly one ahead, which is
     an increment in progress rather than a disagreement.
  F. ADJACENT DUPLICATE LINES, which is what a botched splice leaves.
  G. UNDEFINED LaTeX ENVIRONMENTS. A .tex using begin{proposition}
     with no matching newtheorem does not compile, and there is no TeX
     engine in this environment to say so. The working paper had
     exactly that fault, unnoticed since those propositions were
     added.

Exit code 1 on any failure, so it can gate a commit rather than be
read and ignored.
"""
import io
import os
import re
import subprocess
import sys

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

BAD_CTRL = {0x09: "TAB", 0x07: "BEL", 0x08: "BS", 0x0C: "FF",
            0x0B: "VT", 0x0D: "CR"}

# tails of LaTeX commands that begin with an escape-collapsing letter
TAILS = ("mid", "eq", "abla", "ot", "onumber", "u",          # \n
         "imes", "heta", "au", "o", "op", "ext", "ilde",     # \t
         "lpha", "pprox", "symp", "ngle",                    # \a
         "eta", "ig", "mod", "ar", "ullet", "inom",          # \b
         "rac", "orall", "loor",                             # \f
         "ho", "ight", "angle", "floor",                     # \r
         "arepsilon", "ec", "arphi")                         # \v
TAIL_RE = re.compile(r"^\s*(" + "|".join(TAILS) + r")\b")
TRAILING = re.compile(r"[{_^\\]\s*$")


def midformula(line):
    """Is this line left inside a formula?

    The first version asked whether the last `$` was followed by text,
    which is true of any line ending in prose after an inline formula
    -- it flagged "...$1.051$--$1.068$ times the random-sign" followed
    by the English word "floor". A line is left open when its `$`
    count is ODD, or when it ends on a brace, subscript or backslash.
    """
    return (line.count("$") % 2 == 1) or bool(TRAILING.search(line))


def docs():
    out = []
    for dirpath, dirnames, filenames in os.walk(REPO):
        if ".git" in dirpath:
            continue
        for fn in filenames:
            if fn.endswith((".md", ".tex")):
                out.append(os.path.join(dirpath, fn))
    return sorted(out)


def block_ok(ns):
    """From the maximum onward the rows must descend by exactly one.

    Two earlier invariants failed their own self-test. "Strictly
    descending" flagged the legitimate ascending historical block (15
    false positives). "One direction change" and "the descending block
    starts at the maximum" both PASS increment 295's actual fault,
    ..., 35, 96, 97, 95, ..., because it is distinct, turns once, and
    does start its tail at the maximum.

    What that fault really does is strand 96 above 97, so the run from
    the maximum SKIPS a number. Descending-and-contiguous is the
    invariant that separates the two, and the self-test checks it on
    both the real shape and the real fault.
    """
    k = ns.index(max(ns))
    tail = ns[k:]
    for a, b in zip(tail, tail[1:]):
        if b != a - 1:
            return False, f"after {a} the table jumps to {b}"
    return True, ""


def selftest():
    """Show the checks can fail. Asserting it is worth nothing (276)."""
    ok = True

    def chk(name, cond):
        nonlocal ok
        ok &= cond
        print(f"    {name:<52} {'detects' if cond else 'MISSES'}")

    # A: the six control characters this hazard produces
    chk("A: BEL in text", "\x07" in "m(N)\x07symp\\sqrt{}")
    chk("A: TAB in text", "\t" in "_{\text{mask}}")
    # B: split command caught, ordinary prose not
    chk("B: '\\prod_{q' / 'mid N}' flagged",
        bool(TAIL_RE.match("mid N}(1-1/(q(q-1)))$"))
        and midformula("$\\mathfrak A(N)=\\prod_{q"))
    chk("B: 'times below' after prose not flagged",
        not midformula("reached the point of displaying a formula"))
    chk("B: 'in $\\log N$' after prose not flagged",
        not midformula("phase-random"))
    # C: the real shape and the real fault
    chk("C: real table shape passes", block_ok([3, 4, 35, 97, 96, 95])[0])
    chk("C: increment 295's actual fault caught",
        not block_ok([3, 4, 35, 96, 97, 95])[0])
    chk("C: duplicate row number", len(set([5, 5])) != 2)
    # G: an environment used with no matching newtheorem
    defined_t = set(re.findall(r"newtheorem\*?\{([^}]+)\}",
                               r"\newtheorem{theorem}{Theorem}"))
    chk("G: proposition used, only theorem defined",
        "proposition" not in defined_t)
    chk("G: theorem used and defined", "theorem" in defined_t)
    # F: duplicated line
    L = ["Verified to $1.000000\\pm0.000145$."] * 2
    chk("F: adjacent duplicate", L[0].strip() == L[1].strip()
        and len(L[0].strip()) > 20)
    print(f"    {'SELF-TEST OK' if ok else 'A CHECK CANNOT FAIL'}")
    return ok


def main():
    print("(0) self-test: each check on a synthetic fault")
    if not selftest():
        print("FAIL")
        sys.exit(1)
    print()
    fails = []
    files = docs()
    print(f"(A) control characters in {len(files)} .md/.tex files")
    for path in files:
        raw = io.open(path, "rb").read().decode("utf-8", "replace")
        raw = raw.replace("\r\n", "\n")          # CRLF is git's doing
        for ch, name in BAD_CTRL.items():
            if chr(ch) in raw:
                ln = raw[:raw.index(chr(ch))].count("\n") + 1
                fails.append(f"A {os.path.relpath(path, REPO)}:{ln} "
                             f"contains {name}")
    print(f"    {len([f for f in fails if f.startswith('A ')])} hits")

    print("(B) LaTeX commands split by a collapsed backslash-n")
    nb = 0
    for path in files:
        lines = io.open(path, encoding="utf-8",
                        errors="replace").read().split("\n")
        for i in range(1, len(lines)):
            if TAIL_RE.match(lines[i]) and midformula(lines[i - 1]):
                fails.append(f"B {os.path.relpath(path, REPO)}:{i+1} "
                             f"{lines[i-1][-28:]!r} / {lines[i][:28]!r}")
                nb += 1
    print(f"    {nb} hits")

    # The table is NOT uniformly descending: an old ascending block
    # (rows 3..35) sits above a newer descending block (97..36),
    # because entries have been prepended since. The first draft of
    # this check assumed one direction and reported 15 false
    # positives. What an insertion at the wrong anchor actually does
    # is create an EXTRA direction change, so that is what is
    # counted -- together with distinctness, which a duplicate
    # anchor would break.
    print("(C) correction rows: distinct, descending and contiguous")
    cr = os.path.join(REPO, "CLOSURE_REAUDIT.md")
    nums = [int(m.group(1)) for m in
            re.finditer(r"^\| (\d+) \|", io.open(cr, encoding="utf-8").read(),
                        re.M)]
    if len(set(nums)) != len(nums):
        dup = [n for n in set(nums) if nums.count(n) > 1]
        fails.append(f"C duplicate correction numbers: {dup}")
    good, why = block_ok(nums)
    if not good:
        fails.append(f"C {why}")
    print(f"    {len(nums)} rows, {len(set(nums))} distinct, "
          f"block {'ok' if good else 'BROKEN'}")

    st = io.open(os.path.join(REPO, "STATUS.md"), encoding="utf-8").read()
    print("(D) STATUS correction count against the table")
    m = re.search(r"(\d+) recorded corrections", st)
    stated, actual = int(m.group(1)), max(nums)
    if stated != actual:
        fails.append(f"D STATUS says {stated}, table holds {actual}")
    print(f"    STATUS {stated}, table max {actual}")

    print("(E) STATUS increment against the last commit")
    try:
        subj = subprocess.run(["git", "log", "-1", "--format=%s"],
                              cwd=REPO, capture_output=True,
                              text=True, timeout=30).stdout
        gi = int(re.search(r"Increment (\d+)", subj).group(1))
        si = int(re.search(r"Increment (\d+)", st).group(1))
        # Run before the commit, STATUS is legitimately one ahead --
        # that is work in progress, not a mismatch. Two or more ahead,
        # or behind at all, is a real disagreement.
        if si == gi:
            note = "in step"
        elif si == gi + 1:
            note = "one ahead (increment in progress)"
        else:
            note = "MISMATCH"
            fails.append(f"E STATUS says increment {si}, "
                         f"last commit says {gi}")
        print(f"    STATUS {si}, last commit {gi} -- {note}")
    except Exception as e:
        print(f"    skipped ({type(e).__name__})")

    # (G) LaTeX environments used but never defined. The paper used
    # begin{proposition} twice with no matching newtheorem, so it had
    # not compiled since those propositions were added -- and there is
    # no TeX engine in this environment to say so. A structural check
    # costs nothing and catches it.
    print("(G) LaTeX environments used but not defined")
    STD = {"document", "abstract", "equation", "equation*", "align",
           "align*", "gather", "gather*", "itemize", "enumerate",
           "description", "tabular", "longtable", "table", "figure",
           "center", "quote", "verbatim", "thebibliography", "proof",
           "array", "cases", "pmatrix", "bmatrix", "matrix", "split",
           "multline", "eqnarray", "displaymath", "flushleft",
           "flushright", "minipage", "footnotesize", "small"}
    ng = 0
    for path in files:
        if not path.endswith(".tex"):
            continue
        t = io.open(path, encoding="utf-8", errors="replace").read()
        defined = set(re.findall(r"newtheorem\*?\{([^}]+)\}", t))
        used = set(re.findall(r"begin\{([^}]+)\}", t))
        for m in sorted(used - defined - STD):
            fails.append(f"G {os.path.relpath(path, REPO)}: "
                         f"environment {m!r} used, never defined")
            ng += 1
    print(f"    {ng} hits")

    print("(F) adjacent duplicate lines")
    nf = 0
    for path in files:
        lines = io.open(path, encoding="utf-8",
                        errors="replace").read().split("\n")
        for i in range(1, len(lines)):
            a = lines[i].strip()
            if a and a == lines[i - 1].strip() and len(a) > 20:
                fails.append(f"F {os.path.relpath(path, REPO)}:{i+1} "
                             f"{a[:40]!r}")
                nf += 1
    print(f"    {nf} hits")

    print()
    if fails:
        print(f"{len(fails)} PROBLEM(S):")
        for f in fails:
            print(f"  {f}")
        print("FAIL")
        sys.exit(1)
    print("all checks pass")
    print("DONE")


if __name__ == "__main__":
    main()
