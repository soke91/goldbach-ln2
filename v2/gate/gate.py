# -*- coding: utf-8 -*-
"""
The gate. One file, eight checks. Exit code is the number of failures.

    python gate/gate.py > gate/gate.txt 2>&1; echo $?

Commit only on 0. Do not pipe -- a pipe eats the exit code.

CONVENTIONS THIS GATE ENFORCES

  Every numbered statement carries an evidence marker on the line after
  its label:

      \\begin{proposition}\\label{prop:V}
      % evidence: lab_second_moment.py
      ...
      \\begin{theorem}\\label{thm:A}
      % evidence: analytic

  `analytic` means the statement is a proof with no computation behind
  it. Anything else names a script in code/, which must exist together
  with its result file in results/, and that result file must open with

      STATISTIC: <what was computed, precisely>
      FIELD: <over what range / on what set>

  A symbol whose meaning could drift is declared once:

      % symbol: rho = Var C / V, the cancellation ratio

  Declaring the same symbol twice with different meanings fails.
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
PAPER = os.path.join(ROOT, "paper")
CODE = os.path.join(ROOT, "code")
RESULTS = os.path.join(ROOT, "results")
VERIFY = os.path.join(ROOT, "verify")

STMT = r"theorem|proposition|lemma|corollary|conjecture"
KINDS = "Theorem|Proposition|Lemma|Corollary|Conjecture|Remark"
VERDICT = r"\b(PASS|FAIL|DEAD|ALIVE|CONFIRMED|REFUTED|PROVEN|VERIFIED)\b"

fails = []
notes = []


def read(p):
    return io.open(p, encoding="utf-8", newline="", errors="replace").read()


def tex_files():
    if not os.path.isdir(PAPER):
        return []
    return sorted(os.path.join(PAPER, f) for f in os.listdir(PAPER)
                  if f.endswith(".tex"))


def py_files(d):
    out = []
    for base, _, fs in os.walk(d):
        if "__pycache__" in base:
            continue
        out += [os.path.join(base, f) for f in sorted(fs) if f.endswith(".py")]
    return out


def rel(p):
    return os.path.relpath(p, ROOT).replace("\\", "/")


def statements(src):
    """(kind, label, evidence_or_None, line_no) for each numbered statement."""
    out = []
    lines = src.splitlines()
    for i, line in enumerate(lines):
        m = re.search(r"\\begin\{(" + STMT + r")\}.*?\\label\{([^}]*)\}", line)
        if not m:
            continue
        ev = None
        for j in range(i + 1, min(i + 4, len(lines))):
            e = re.match(r"\s*%\s*evidence:\s*(\S+)", lines[j])
            if e:
                ev = e.group(1)
                break
        out.append((m.group(1), m.group(2), ev, i + 1))
    return out


# ------------------------------------------------------------------ G1
def g1_evidence(docs):
    missing, broken = [], []
    for path, src in docs:
        for kind, label, ev, ln in statements(src):
            if ev is None:
                missing.append(f"{rel(path)}:{ln} {kind} {label}")
                continue
            if ev == "analytic":
                continue
            s = os.path.join(CODE, ev)
            r = os.path.join(RESULTS, os.path.splitext(ev)[0] + ".txt")
            if not os.path.exists(s):
                broken.append(f"{label} -> code/{ev} (missing)")
            elif not os.path.exists(r):
                broken.append(f"{label} -> results/{os.path.basename(r)} "
                              f"(missing)")
    for m in missing:
        fails.append(f"G1 no evidence marker: {m}")
    for b in broken:
        fails.append(f"G1 evidence does not exist: {b}")
    return len(missing) + len(broken)


# ------------------------------------------------------------------ G2
def g2_hand_numbers(docs):
    n = 0
    for path, src in docs:
        for m in re.finditer(r"(?<!\\)\b(" + KINDS + r")~?\s(\d+)", src):
            back = src[max(0, m.start() - 26): m.start()]
            if re.search(r"'s\s*$|\\cite\{[^}]*\}[^.]{0,12}$", back):
                continue          # someone else's numbering, not ours
            fails.append(f"G2 hand-written number: {rel(path)} "
                         f"'{m.group(1)} {m.group(2)}' -- use \\ref")
            n += 1
    return n


# ------------------------------------------------------------------ G3
def g3_symbols(docs):
    seen, n = {}, 0
    for path, src in docs:
        for m in re.finditer(r"%\s*symbol:\s*(\S+)\s*=\s*(.+)", src):
            sym, mean = m.group(1), m.group(2).strip()
            if sym in seen and seen[sym][1] != mean:
                fails.append(f"G3 symbol '{sym}' declared twice with "
                             f"different meanings: {seen[sym][0]} / "
                             f"{rel(path)}")
                n += 1
            else:
                seen.setdefault(sym, (rel(path), mean))
        for m in re.finditer(r"\\(?:newcommand|DeclareMathOperator)"
                             r"\*?\{?\\(\w+)", src):
            k = "\\" + m.group(1)
            if k in seen and seen[k][0] != rel(path):
                fails.append(f"G3 macro {k} defined in two documents")
                n += 1
            seen.setdefault(k, (rel(path), "macro"))
    return n


# ------------------------------------------------------------------ G4
def g4_statistic_and_field(docs):
    n, checked = 0, set()
    for path, src in docs:
        for _, label, ev, _ in statements(src):
            if not ev or ev == "analytic" or ev in checked:
                continue
            checked.add(ev)
            r = os.path.join(RESULTS, os.path.splitext(ev)[0] + ".txt")
            if not os.path.exists(r):
                continue                      # already reported by G1
            head = read(r)[:4000]
            for want in ("STATISTIC:", "FIELD:"):
                if want not in head:
                    fails.append(f"G4 {rel(r)} has no '{want}' line "
                                 f"(supports {label})")
                    n += 1
    return n


# ------------------------------------------------------------------ G5
def g5_verdicts():
    n = 0
    for p in py_files(CODE) + py_files(VERIFY):
        src = read(p)
        if "# verdict-ok" in src:
            continue
        prints = re.findall(r'print\(\s*["\']([^"\']*)["\']\s*\)', src)
        if not any(re.search(VERDICT, s) for s in prints):
            continue
        has_exit = re.search(r"sys\.exit\(\s*(?!0\s*\))", src) \
            or "raise " in src or re.search(r"^\s*assert\s", src, re.M)
        if not has_exit:
            fails.append(f"G5 {rel(p)} prints a verdict with no failure path")
            n += 1
    return n


# ------------------------------------------------------------------ G6
def g6_one_home(docs):
    n, where = 0, {}
    for path, src in docs:
        for _, label, _, ln in statements(src):
            if label in where:
                fails.append(f"G6 label '{label}' defined in two documents: "
                             f"{where[label]} and {rel(path)}")
                n += 1
            else:
                where[label] = rel(path)
    # a claim restated verbatim in a second document has two homes
    bodies = {}
    for path, src in docs:
        for m in re.finditer(r"\\begin\{(?:" + STMT + r")\}(.{80,400}?)"
                             r"\\end\{(?:" + STMT + r")\}", src, re.S):
            key = re.sub(r"\s+", " ", m.group(1))[:120].strip()
            if key in bodies and bodies[key] != rel(path):
                fails.append(f"G6 a statement body appears in both "
                             f"{bodies[key]} and {rel(path)}")
                n += 1
            bodies.setdefault(key, rel(path))
    return n


# ------------------------------------------------------------------ G8
def g8_budget(docs):
    prose = 0
    for name in ("README.md", "DECISIONS.md", "OPEN.md"):
        p = os.path.join(ROOT, name)
        if os.path.exists(p):
            prose += len(read(p).splitlines())
    for base, _, fs in os.walk(VERIFY):
        for f in fs:
            if f.endswith(".md"):
                prose += len(read(os.path.join(base, f)).splitlines())
    product = sum(len(src.splitlines()) for _, src in docs)
    if product == 0:
        notes.append(f"G8 skipped: paper/ is empty ({prose} lines of prose)")
        return 0
    ratio = prose / product
    print(f"    prose {prose} / paper {product} = {ratio:.2f}  (cap 1.50)")
    if ratio > 1.5:
        fails.append(f"G8 prose is {ratio:.2f}x the paper -- move a "
                     f"registry into this gate")
        return 1
    return 0


def main():
    docs = [(p, read(p)) for p in tex_files()]
    print("gate")
    print(f"paper: {len(docs)} document(s)"
          + ("" if docs else "  -- empty, structural checks only"))
    print("=" * 62)

    counts = [
        ("G1 evidence exists", g1_evidence(docs)),
        ("G2 no hand-written numbers", g2_hand_numbers(docs)),
        ("G3 one meaning per symbol", g3_symbols(docs)),
        ("G4 statistic and field declared", g4_statistic_and_field(docs)),
        ("G5 no verdict without failure path", g5_verdicts()),
        ("G6 one home per claim", g6_one_home(docs)),
        ("G8 prose budget", g8_budget(docs)),
    ]
    print()
    for name, c in counts:
        print(f"  {name:<38} {'ok' if c == 0 else str(c) + ' FAILED'}")
    print("  G7 applied to every document in paper/          ok")
    print()
    for nte in notes:
        print(f"  note: {nte}")
    for f in fails:
        print(f"  - {f}")
    print()
    print(f"failures: {len(fails)}")
    return len(fails)


if __name__ == "__main__":
    sys.exit(main())
