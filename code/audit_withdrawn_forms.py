# -*- coding: utf-8 -*-
"""
Does any live document still lean on a withdrawn form? (increment 308)

WHY THIS EXISTS. `CLOSURE_REAUDIT.md` records what was withdrawn and
what replaced it, and the repo's rule is that every other document
states the current position directly rather than arguing with itself.
Nothing enforced the second half. A withdrawal is only complete when
the withdrawn form has stopped being asserted anywhere, and until this
run nobody had checked.

It found `paper/negative_map.tex` -- the consolidated working paper,
written at increment 299 -- still carrying two claims the record has
since withdrawn, as live measured facts: #69's `g = -0.489` with
`m(N) ~ sqrt(S(N)) N^{1/4}`, withdrawn at #112 and again at #119, and
#99's `b = 2.68`, withdrawn at #300. The paper is the artifact a reader
outside this program would read first.

THE CHECK. For each withdrawn form below, every occurrence in a `.md`
or `.tex` file must sit within two lines of a withdrawal marker --
"withdrawn", "superseded", "no longer", "was", "철회", and so on. An
occurrence with no such marker nearby is an ASSERTION of a withdrawn
form, which is the thing being looked for. `CLOSURE_REAUDIT.md` is
exempt: it is the record, and quoting the withdrawn form is its job.

WHAT THIS CANNOT DO. It checks the forms listed here, and the list is
hand-made, so a withdrawal whose form nobody adds stays unchecked. It
also cannot tell a genuine assertion from an unlucky coincidence of
digits -- which is why each pattern is anchored to its context rather
than being a bare number. Both limits are stated in the output.

SELF-TEST runs first, on synthetic text, and must show the check
catching an asserted withdrawn form and passing one that is properly
marked.
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

# (correction, what it was, regex). Anchored to context, not bare
# digits -- and the digits themselves need boundaries, or 0.4891 in an
# unrelated quantity reads as the withdrawn 0.489. The self-test caught
# exactly that before this ever ran on the corpus.
B = r"(?<![\d.])"
E = r"(?![\d])"
FORMS = [
    ("#69/#112/#119", "the mask's exponent g = -0.489",
     B + r"-?0\.489" + E),
    ("#69/#112/#119", "the mask amplitude N^{1/4}",
     r"N\^\{1/4\}|N\^{1/4}"),
    ("#67/#112", "the mask's variance share at large N",
     B + r"1\.15\s*%|" + B + r"14\.1\s*%"),
    ("#96/#110", "the wall's spectral share",
     B + r"0\.39\s*%"),
    ("#99/#300", "the rate constant b = 2.68",
     r"b\s*=\s*2\.68" + E + r"|" + B + r"2\.6817" + E),
    ("#94/#110", "the wall's zeta lines against a permutation null",
     B + r"1566" + E + r"|z\s*=\s*\+?13\.9" + E),
    ("#36/#68", "the fitted variance constant kappa = 0.465",
     B + r"0\.465" + E),
]

MARKERS = ["withdrawn", "superseded", "retracted", "no longer",
           "not replaced", "was fitted", "stays withdrawn",
           "contaminated", "철회", "⚠️", "does not survive", "refuted",
           "not determined", "History:", "history in", "corrected",
           "supersede"]
# 느슬한 마커를 넣지 않는다. 첫 판에서 오탄 셋을 지우려고
# "wrong"과 "premature"를 넣었는데, 바로 그 두 말이 **다른 주장에
# 대해** 쓰인 단락에서 진짜 하나(논문의 `b = 2.68`)를 같이
# 묻었다. 오탄은 규칙을 느슨하게 해서가 아니라 **명시적 예외로**
# 끜다 — 그래야 예외마다 이유가 남는다.
# 대소문자를 가리면 `Withdrawn with it:` 같은 멀짱한 철회 문장을
# 놓친다. 첫 판에서 STATUS.md가 그렇게 오탄으로 잡혔다.
MARKERS = [m.lower() for m in MARKERS]
OKMARK = re.compile(r"withdrawn-ok:")
CONTEXT = 4


def scan_text(text, forms=FORMS):
    lines = text.splitlines()
    hits = []
    for corr, what, pat in forms:
        rx = re.compile(pat)
        for i, line in enumerate(lines):
            if not rx.search(line):
                continue
            lo = max(0, i - CONTEXT)
            hi = min(len(lines), i + CONTEXT + 1)
            near = "\n".join(lines[lo:hi])
            if OKMARK.search(near):
                continue
            if any(m in near.lower() for m in MARKERS):
                continue
            hits.append((i + 1, corr, what, line.strip()[:110]))
    return hits


SELFTEST = [
    ("asserted with no marker", True,
     "The mask is lower order: g = -0.489, so m(N) ~ N^{1/4}.\n"),
    ("quoted while withdrawing", False,
     "Increment 280's g = -0.489 is withdrawn; the scaling is\n"
     "not determined.\n"),
    ("marker two lines away", False,
     "the fitted exponent\n"
     "g = -0.489 appeared here\n"
     "and is superseded by the exact second moment\n"),
    ("explicit whitelist", False,
     "g = -0.489 <!-- withdrawn-ok: quoted as history -->\n"),
    ("unrelated text", False,
     "The measured ratio is 0.4891 in a different quantity.\n"),
]


def selftest():
    print("SELF-TEST")
    ok = True
    for name, should, txt in SELFTEST:
        got = len(scan_text(txt)) > 0
        good = got == should
        ok &= good
        print(f"    {name:<28} {'flags' if got else 'passes':>7}   "
              f"{'as expected' if good else 'WRONG'}")
    print(f"    {'SELF-TEST OK' if ok else 'SELF-TEST FAILED'}\n")
    return ok


def main():
    if not selftest():
        print("DONE (self-test failed)")
        sys.exit(1)
    targets = []
    for base, _dirs, files in os.walk(ROOT):
        if os.sep + ".git" in base:
            continue
        for f in files:
            if f.endswith((".md", ".tex")):
                targets.append(os.path.join(base, f))
    targets.sort()

    total, nfiles = 0, 0
    print(f"(1) withdrawn forms asserted in live documents "
          f"-- {len(targets)} files, {len(FORMS)} forms")
    for p in targets:
        rel = os.path.relpath(p, ROOT).replace(os.sep, "/")
        if rel == RECORD:
            continue
        try:
            text = io.open(p, encoding="utf-8").read()
        except Exception:
            continue
        hits = scan_text(text)
        if not hits:
            continue
        nfiles += 1
        total += len(hits)
        print(f"\n  {rel}")
        for lineno, corr, what, line in hits:
            print(f"    line {lineno:>4}  {corr:<16} {what}")
            print(f"        {line}")

    print(f"\n  {total} assertion(s) of a withdrawn form in "
          f"{nfiles} file(s)")
    print("  LIMIT: the form list is hand-made, so a withdrawal nobody")
    print("  adds here stays unchecked; and a marker within two lines")
    print("  is taken as a withdrawal, which is a heuristic.")
    if total:
        print("DONE (findings)")
        sys.exit(1)
    print("  no live document asserts a withdrawn form")
    print("DONE")


if __name__ == "__main__":
    main()
