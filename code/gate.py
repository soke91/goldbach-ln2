# -*- coding: utf-8 -*-
"""
One gate, one exit code (increment 323)

WHY. #161: every commit this session ran the checkers as

    python code/lint_docs.py 2>&1 | tail -4 && git commit ...

and a pipeline's exit status is the LAST command's, so `tail` succeeding
masked the linter exiting 1. Increment 321's follow-up went out with two
live failures. The invocation was fixed by hand; this makes the mistake
unavailable. There is now one thing to run, it prints one line per
checker, and it exits nonzero if any of them does.

It also runs a RETROACTIVE audit that #161 left open: the bookkeeping
invariant (STATUS's stated correction count equals the highest row in
CLOSURE_REAUDIT) is computable from any commit without a working tree,
so it can be replayed over the session's history. If the gate never
fired, the question of whether anything shipped broken is answerable
rather than assumed.

CHECKERS RUN
  lint_docs.py               escape collapse, numbering, counts, LaTeX
  lint_gates.py              a verdict announced with no failure path
  lint_verdicts.py           a conclusion composed rather than computed
  audit_withdrawn_forms.py   a withdrawn form still asserted
  audit_quoted_numbers.py    a figure that drifted from results/

`verify_all.py` and `verify_deep.py` are NOT run here: they are minutes
long and are the reproduction stamp, not the commit gate. Run them
before a release, not before every commit.
"""
import io
import os
import re
import subprocess
import sys

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)

CHECKERS = [
    "lint_docs.py",
    "lint_gates.py",
    "lint_verdicts.py",
    "audit_withdrawn_forms.py",
    "audit_quoted_numbers.py",
]

CNT = re.compile(r"(\d+)\s+recorded corrections")
INC = re.compile(r"Increment\s+(\d+)")
ROW = re.compile(r"^\|\s*(\d{1,3})\s*\|", re.M)


def git(*args):
    p = subprocess.run(["git"] + list(args), cwd=REPO,
                       capture_output=True)
    return p.returncode, p.stdout.decode("utf-8", "replace")


def retro(n=30):
    """Replay the count invariant over the last n commits."""
    rc, out = git("log", "--format=%h %s", "-n", str(n))
    if rc:
        return None
    bad = []
    seen = 0
    for line in out.strip().split("\n"):
        if not line.strip():
            continue
        h = line.split()[0]
        rc1, st = git("show", f"{h}:STATUS.md")
        rc2, cr = git("show", f"{h}:CLOSURE_REAUDIT.md")
        if rc1 or rc2:
            continue
        mc = CNT.search(st)
        mi = INC.search(st)
        rows = [int(x) for x in ROW.findall(cr)]
        if not mc or not rows:
            continue
        seen += 1
        stated, top = int(mc.group(1)), max(rows)
        if stated != top:
            bad.append((h, mi.group(1) if mi else "?", stated, top,
                        line.split(" ", 1)[1][:56]))
    return seen, bad


def main():
    print("(1) checkers, by exit code")
    fails = []
    for c in CHECKERS:
        p = subprocess.run([sys.executable, os.path.join(HERE, c)],
                           capture_output=True, cwd=REPO)
        ok = p.returncode == 0
        if not ok:
            fails.append(c)
        print(f"    {c:<28} {'PASS' if ok else 'FAIL':>4}   "
              f"exit {p.returncode}")
        if not ok:
            tail = p.stdout.decode("utf-8", "replace").strip().split("\n")
            for ln in tail[-6:]:
                print(f"        {ln}")

    print(f"\n(2) retroactive: the count invariant over the last "
          f"30 commits")
    r = retro()
    if r is None:
        print("    not a git repository, or git unavailable")
    else:
        seen, bad = r
        print(f"    {seen} commits carried both files")
        if bad:
            for h, inc, stated, top, subj in bad:
                print(f"    {h}  inc {inc}: STATUS says {stated}, "
                      f"table max {top}   {subj}")
            print(f"    ⚠️  {len(bad)} commit(s) shipped with the "
                  f"count invariant broken")
        else:
            print("    none broken -- the gate's bypass never let a "
                  "bookkeeping fault through, which is luck rather "
                  "than design")

    if fails:
        print(f"\n{len(fails)} checker(s) failed: "
              f"{', '.join(fails)}")
        print("DONE (blocked)")
        sys.exit(1)
    print("\nall checkers pass")
    print("DONE")


if __name__ == "__main__":
    main()
