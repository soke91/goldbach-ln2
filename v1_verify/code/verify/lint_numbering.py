r"""
Resolve every numbered statement in v1/paper/wall_v1.tex to its number,
and check that the numbers cited elsewhere in the repository agree.

The paper declares one shared counter --
    \newtheorem{theorem}{Theorem}
    \newtheorem{conjecture}[theorem]{Conjecture}   ... etc, all [theorem]
-- so every theorem, conjecture, corollary, proposition, lemma AND
remark advances the same count, in source order. Documents that cite
"Proposition 15" or "Lemma 18" are citing that count, and a remark
inserted or removed anywhere above shifts everything below it.

This prints the resolved table and then greps the live documents for
"<Kind> <n> (`label`)" citations and flags any that disagree.
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
PAPERS = [os.path.join(ROOT, "v1", "paper", "wall_v1.tex"),
          os.path.join(ROOT, "v1", "paper", "theorem_A.tex")]

SHARED = {"theorem", "conjecture", "corollary", "proposition", "lemma",
          "remark"}
KIND = {"theorem": "Theorem", "conjecture": "Conjecture",
        "corollary": "Corollary", "proposition": "Proposition",
        "lemma": "Lemma", "remark": "Remark"}


def resolve(path):
    src = io.open(path, encoding="utf-8").read()
    body = src.split(r"\begin{document}", 1)[1]
    out, n = [], 0
    for m in re.finditer(r"\\begin\{(" + "|".join(SHARED) + r")\}"
                         r"(\[[^\]]*\])?(\s*\\label\{([^}]+)\})?", body):
        env = m.group(1)
        n += 1
        lab = m.group(4)
        if lab is None:
            tail = body[m.end():m.end() + 260]
            lm = re.search(r"\\label\{([^}]+)\}", tail)
            lab = lm.group(1) if lm else None
        out.append((n, KIND[env], lab))
    return out


def main():
    bylabel = {}
    for path in PAPERS:
        table = resolve(path)
        rel = os.path.relpath(path, ROOT).replace(os.sep, "/")
        print(f"numbered statements in {rel} "
              f"(one shared counter, {len(table)} entries)")
        print()
        for n, k, lab in table:
            print(f"  {k:>11} {n:>2}   {lab if lab else '(unlabelled)'}")
            if lab and lab not in bylabel:
                bylabel[lab] = (n, k)
        print()

    # citations of the form "Kind N (`label`)" or "Kind N (label)"
    docs = []
    for base, _d, fs in os.walk(ROOT):
        if any(s in base for s in (".git", "__pycache__", "archive")):
            continue
        for f in fs:
            if f.endswith((".md", ".tex")):
                docs.append(os.path.join(base, f))
    pat = re.compile(
        r"\b(Theorem|Conjecture|Corollary|Proposition|Lemma|Remark)s?\s+"
        r"(\d+)\s*\(`?\\?r?e?f?\{?([a-zA-Z:]+)\}?`?\)")
    print()
    print("citations of the form  <Kind> <n> (`label`)  in live documents")
    bad = 0
    seen = 0
    for p in sorted(docs):
        rel = os.path.relpath(p, ROOT).replace(os.sep, "/")
        text = io.open(p, encoding="utf-8", errors="replace").read()
        for m in pat.finditer(text):
            kind, num, lab = m.group(1), int(m.group(2)), m.group(3)
            if lab not in bylabel:
                continue
            seen += 1
            n, k = bylabel[lab]
            if n != num or k != kind:
                bad += 1
                line = text[:m.start()].count("\n") + 1
                print(f"  MISMATCH {rel}:{line}  cites {kind} {num} "
                      f"for `{lab}`, which is {k} {n}")
    print(f"  {seen} resolvable citations, {bad} mismatched")
    print()
    if bad:
        print("DONE (failed)")
        sys.exit(1)
    print("DONE")


if __name__ == "__main__":
    main()
