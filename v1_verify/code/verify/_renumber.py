r"""
One-shot: fix the statement numbers this tree inherited from
v1/PROVENANCE.md, which is off by one from `conj:L` onward.

Only v1_verify's own documents are touched. v1 is frozen.
"""
import io
import os
import sys

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))

FIX = [
    ("Conjecture 11 (`conj:L`)", "Conjecture 10 (`conj:L`)"),
    ("Conjecture 11 (`conj:L`)**", "Conjecture 10 (`conj:L`)**"),
    ("Proposition 12 (`prop:V`)", "Proposition 11 (`prop:V`)"),
    ("Proposition 12 of the same paper", "Proposition 11 of the same paper"),
    ("Proposition 12 defines", "Proposition 11 defines"),
    ("Proposition 12 (the `A`", "Proposition 11 (the `A`"),
    ("Proposition 21 (`prop:coh`)", "Proposition 20 (`prop:coh`)"),
    ("Lemma 20 gives the error bar", "Lemma 19 gives the error bar"),
    ("Lemma 20's three terms", "Lemma 19's three terms"),
    ("Lemma 18 (`lem:coin`)", "Lemma 17 (`lem:coin`)"),
    ("Lemma 18, applied to", "Lemma 17, applied to"),
    ("Lemma 18 says", "Lemma 17 says"),
    ("of Lemma 18", "of Lemma 17"),
    ("Lemma 18's control", "Lemma 17's control"),
    ("Proposition 12 of the paper", "Proposition 11 of the paper"),
]

TARGETS = [
    "v1_verify/README.md",
    "v1_verify/paper/ADVERSARIAL_FINDINGS.md",
]


def main():
    total = 0
    for rel in TARGETS:
        p = os.path.join(ROOT, rel)
        if not os.path.exists(p):
            continue
        t = io.open(p, encoding="utf-8").read()
        n = 0
        for a, b in FIX:
            c = t.count(a)
            if c:
                t = t.replace(a, b)
                n += c
        if n:
            io.open(p, "w", encoding="utf-8").write(t)
        print(f"  {rel}: {n} replacement(s)")
        total += n
    print(f"{total} replacement(s) total")
    print("DONE")


if __name__ == "__main__":
    main()
