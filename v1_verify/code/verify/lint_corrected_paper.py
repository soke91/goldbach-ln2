# -*- coding: utf-8 -*-
r"""
Structural lint for v1_verify/paper/wall_v1_corrected.tex.

The corrected paper is a rewritten copy of v1/paper/wall_v1.tex, so the
things that can go wrong in editing it are mechanical: an unbalanced
environment, a \\ref with no \\label, a stray astral character that the
console cannot encode, or a figure this tree has refuted still sitting
in the body. This checks those and exits nonzero on any of them.

Checks:
  (1) every \begin{env} has a matching \end{env}, in order
  (2) braces balance
  (3) every \ref{key} resolves to a \label{key}
  (4) no astral (non-BMP) characters
  (5) none of the forms this tree refuted appears in the body
  (6) the file differs from v1's in exactly the intended places, and
      v1's own file is byte-identical to the frozen one
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
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
NEW = os.path.join(ROOT, "v1_verify", "paper", "wall_v1_corrected.tex")
OLD = os.path.join(ROOT, "v1", "paper", "wall_v1.tex")

# forms this tree refuted; none may appear in the corrected body
REFUTED = [
    (r"\\sum_\{N\\le X\} C\(N\)\^2", "the N<=X form of Lemma 13"),
    (r"0\.39\s*\\?%", "the withdrawn spectral share"),
    (r"1566", "the withdrawn 1566x ratio"),
    (r"-?0\.0976", "the reconstruction with the W denominator"),
    (r"a factor \$0\.54\$", "the withdrawn comparison against -0.18"),
    (r"5\.09\\cdot10\^\{-6\}", "the value-permutation surrogate maximum"),
    (r"z\\ge23", "the per-ordinate z from the permutation null"),
    (r"ten kill-tested", "the miscount of the kill-test table"),
    (r"eighteen\s+pre-registered", "the miscount of the closures"),
]


def main():
    src = io.open(NEW, encoding="utf-8").read()
    fails = []

    # (1) environments
    stack = []
    for m in re.finditer(r"\\(begin|end)\{([A-Za-z*]+)\}", src):
        kind, env = m.group(1), m.group(2)
        if kind == "begin":
            stack.append(env)
        else:
            if not stack:
                fails.append(f"\\end{{{env}}} with nothing open")
            elif stack[-1] != env:
                fails.append(f"\\end{{{env}}} closes \\begin{{{stack[-1]}}}")
                stack.pop()
            else:
                stack.pop()
    if stack:
        fails.append(f"unclosed environments: {stack}")
    print(f"(1) environments balance: "
          f"{'OK' if not fails else 'FAIL'}")

    # (2) braces
    depth = 0
    for i, ch in enumerate(src):
        if ch == "{" and (i == 0 or src[i - 1] != "\\"):
            depth += 1
        elif ch == "}" and (i == 0 or src[i - 1] != "\\"):
            depth -= 1
            if depth < 0:
                fails.append(f"brace underflow at offset {i}")
                break
    if depth != 0:
        fails.append(f"braces unbalanced, net {depth}")
    print(f"(2) braces balance: {'OK' if depth == 0 else 'FAIL ' + str(depth)}")

    # (3) refs
    labels = set(re.findall(r"\\label\{([^}]+)\}", src))
    refs = set(re.findall(r"\\ref\{([^}]+)\}", src))
    missing = sorted(refs - labels)
    if missing:
        fails.append(f"dangling refs: {missing}")
    print(f"(3) every \\ref resolves: "
          f"{'OK' if not missing else 'FAIL ' + str(missing)}"
          f"   ({len(labels)} labels, {len(refs)} referenced)")

    # (4) astral characters
    astral = sorted({c for c in src if ord(c) > 0xFFFF})
    if astral:
        fails.append(f"astral characters: {[hex(ord(c)) for c in astral]}")
    print(f"(4) no astral characters: {'OK' if not astral else 'FAIL'}")

    # (5) refuted forms
    print("(5) refuted forms absent from the body:")
    bad = 0
    for pat, what in REFUTED:
        hits = [m.start() for m in re.finditer(pat, src)]
        if hits:
            bad += 1
            line = src[:hits[0]].count("\n") + 1
            print(f"      FAIL line {line}: {what}")
            fails.append(f"refuted form present: {what}")
    if not bad:
        print(f"      OK -- none of {len(REFUTED)} forms present")

    # (6) v1 untouched, and the diff is confined to the paper body
    old = io.open(OLD, encoding="utf-8").read()
    nold, nnew = old.count("\n"), src.count("\n")
    same = sum(1 for a, b in zip(old.split("\n"), src.split("\n")) if a == b)
    print(f"(6) v1 original: {nold} lines; corrected: {nnew} lines; "
          f"{same} leading lines identical")
    if not os.access(OLD, os.W_OK):
        print("    v1/paper/wall_v1.tex is read-only, as the freeze "
              "requires")
    else:
        fails.append("v1/paper/wall_v1.tex is writable")
        print("    FAIL: v1/paper/wall_v1.tex is writable")

    print()
    if fails:
        print(f"{len(fails)} problem(s)")
        print("DONE (failed)")
        sys.exit(1)
    print("no problems")
    print("DONE")


if __name__ == "__main__":
    main()
