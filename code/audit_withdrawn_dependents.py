# -*- coding: utf-8 -*-
"""
What still leans on a withdrawn claim? (increment 310)

WHY THIS IS THE MISSING REGISTER. `OPEN_QUESTIONS.md` says a withdrawn
POSITIVE claim "costs nothing but a claim", and that is wrong. A
positive claim that was used to derive something else takes the
derivation with it, and the derived statement does not get withdrawn
automatically -- it just stops having support. Nothing in this repo
tracks that.

Two registers already exist and neither covers it:

  Register A   holes with nothing in them -- the question the withdrawn
               claim answered.
  Register B   closures whose premise a correction moved.
  MISSING      **live claims that were DERIVED from a withdrawn one.**

`code/audit_withdrawn_forms.py` catches the textual case -- a document
still printing a withdrawn number. This is the logical case: a document
stating something else, correctly, on the strength of a claim that no
longer stands. No regex finds that. What a machine can do is narrow the
reading list, which is what this does.

WHAT IT PRODUCES. For every correction that withdrew a positive claim,
the list of places in live documents that cite it by number. Each
citation is one of:

  RECORD    "history in CLOSURE_REAUDIT #69" -- a pointer, harmless.
  SUPPORT   the citing sentence rests on the withdrawn claim. This is
            what has to be found and re-derived.

The classification between those two is a judgement and is NOT
automated here; the tool prints the citing line so the judgement can be
made and recorded. Pretending to automate it would be the same fault as
increment 308's Register B, which asserted a property of eleven files
without opening them and was wrong about twelve of thirteen (#134).

SELF-TEST runs first.
"""
import io
import os
import re
import sys

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RECORD = "CLOSURE_REAUDIT.md"

# Corrections that withdrew a POSITIVE claim (as opposed to repairing a
# tool, a document, or a test). Hand-made and stated as such.
WITHDRAWN_POSITIVES = {
    36: "the wall's variance law Var C = k*S*N*log N",
    47: "the cancellation ratio's decay exponent",
    67: "the mask's share of the variance",
    68: "the wall's variance exponent in log N",
    69: "the mask's scaling N^{1/4}",
    84: "the direction of rho's trend",
    86: "the quantitative half of the demasking claim",
    94: "the wall carries the zeta zeros",
    96: "the wall's spectral share 0.39%",
    99: "the rate constant b = 2.68",
    110: "the spectral lines are the wall's",
    112: "the mask's measured scaling",
    118: "the mask's amplitude is unresolved at large N",
    119: "a single exponent for the mask's decay",
}

CITE = re.compile(r"#(\d{1,3})\b")
POINTER = re.compile(
    r"history|History|CLOSURE_REAUDIT|withdrawn|superseded|철회|"
    r"corrections?\s*#|recorded (?:at|in)")


def scan(text):
    """Returns [(lineno, correction, line, looks_like_pointer)]."""
    out = []
    for i, line in enumerate(text.splitlines()):
        for m in CITE.finditer(line):
            n = int(m.group(1))
            if n in WITHDRAWN_POSITIVES:
                out.append((i + 1, n, line.strip(),
                            bool(POINTER.search(line))))
    return out


SELFTEST = [
    ("pointer", "History: CLOSURE_REAUDIT.md #36, #67.", True),
    ("support", "Because the mask scales as #69 gives, the tail is safe.",
     False),
    ("not a withdrawn one", "See #1 and #2 for the frame.", None),
]


def selftest():
    print("SELF-TEST")
    ok = True
    for name, src, want in SELFTEST:
        hits = scan(src)
        if want is None:
            good = len(hits) == 0
            got = "no hit"
        else:
            good = len(hits) > 0 and hits[0][3] == want
            got = ("pointer" if hits and hits[0][3] else "support") \
                if hits else "no hit"
        ok &= good
        print(f"    {name:<22} -> {got:<10} "
              f"{'as expected' if good else 'WRONG'}")
    print(f"    {'SELF-TEST OK' if ok else 'SELF-TEST FAILED'}\n")
    return ok


def main():
    if not selftest():
        print("DONE (self-test failed)")
        sys.exit(1)
    targets = []
    for base, _d, files in os.walk(ROOT):
        if os.sep + ".git" in base:
            continue
        for f in files:
            if f.endswith((".md", ".tex")):
                targets.append(os.path.join(base, f))
    targets.sort()

    npoint, nread = 0, 0
    print(f"(1) citations of a withdrawn positive claim, in live "
          f"documents")
    for p in targets:
        rel = os.path.relpath(p, ROOT).replace(os.sep, "/")
        if rel == RECORD or rel == "OPEN_QUESTIONS.md":
            continue
        try:
            text = io.open(p, encoding="utf-8").read()
        except Exception:
            continue
        hits = scan(text)
        if not hits:
            continue
        shown = [h for h in hits if not h[3]]
        npoint += len(hits) - len(shown)
        nread += len(shown)
        if not shown:
            continue
        print(f"\n  {rel}")
        for lineno, n, line, _ in shown:
            print(f"    line {lineno:>4}  #{n} — {WITHDRAWN_POSITIVES[n]}")
            print(f"        {line[:104]}")

    print(f"\n(2) tally")
    print(f"    citations that look like pointers to the record: {npoint}")
    print(f"    citations that need reading:                     {nread}")
    print(f"    withdrawn positive claims tracked:               "
          f"{len(WITHDRAWN_POSITIVES)}")
    print(f"\n    Every line above is a place where a live document")
    print(f"    mentions a withdrawn claim OUTSIDE a history pointer.")
    print(f"    Whether it merely mentions it or RESTS on it is a")
    print(f"    judgement, and this tool does not make it -- the")
    print(f"    reading list is the deliverable. Automating the")
    print(f"    judgement is what produced #134.")
    print("DONE")


if __name__ == "__main__":
    main()
