# -*- coding: utf-8 -*-
"""
Verdict linter: a conclusion must be computed, not composed (inc. 298).

WHY. Correction #100. `lab_rho_rate.py` printed "the theory-selected
one-parameter model is not worse than the best free two-parameter fit"
-- a sentence drafted alongside the pre-registration, describing the
outcome I expected. The run rejected the model. The text would have
been recorded as a finding if I had not read the output against it.

That is the mirror of corrections #61, #71 and #78, which were checks
that could not come out FALSE. This is a conclusion that cannot come
out WRONG, because it was written before the experiment and never
consults it. Increment 296's lesson applies unchanged: naming a hazard
does not prevent it; a check that fails loudly does.

WHAT IT FLAGS. A `print` of a plain string literal -- no interpolated
value -- containing verdict vocabulary (SURVIVES, REFUTED, ESTABLISHED,
CONFIRMED, PREFERRED, REJECTED, DEAD, ALIVE, settled, confirms, ...),
where the print is not inside an `if`. Such a line asserts an outcome
the run cannot contradict.

WHAT IT DOES NOT FLAG, and why the distinction is the whole design:

  * f-strings carrying a computed value -- the number can refute the
    sentence around it, which is the point;
  * prints inside a conditional -- the branch IS the computation;
  * lines marked `# verdict-ok: <reason>`, for the two legitimate
    cases that recur here:
      - CRITERION statements ("Pre-registered: ALIVE if x <= 0.3"),
        which state the rule in advance and must be fixed text;
      - STRUCTURAL statements ("expected BY CONSTRUCTION and confirms
        nothing"), true of the method regardless of any data.
    The marker forces the distinction to be made explicitly rather
    than assumed, which is the only reason a whitelist is tolerable
    in a program with this error rate.

Exit 1 on any unmarked hit.
"""
import ast
import io
import os
import re
import sys

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CODE = os.path.join(REPO, "code")

VERDICT = re.compile(
    r"\b(SURVIVES|REFUTED|ESTABLISHED|CONFIRMED|PREFERRED|REJECTED|"
    r"DEAD|ALIVE|INDETERMINATE|settles?|settled|is preferred|"
    r"not worse|confirms?|proves?|shows that)\b")
OKMARK = re.compile(r"#\s*verdict-ok:")


def scan(path):
    src = io.open(path, encoding="utf-8", errors="replace").read()
    lines = src.split("\n")
    try:
        tree = ast.parse(src)
    except SyntaxError as e:
        return [(getattr(e, "lineno", 0), f"SYNTAX: {e.msg}")]
    inside_if = set()
    for node in ast.walk(tree):
        if isinstance(node, (ast.If, ast.IfExp)):
            for sub in ast.walk(node):
                inside_if.add(id(sub))
    out = []
    for node in ast.walk(tree):
        if not (isinstance(node, ast.Call)
                and isinstance(node.func, ast.Name)
                and node.func.id == "print"):
            continue
        if id(node) in inside_if:
            continue
        for a in node.args:
            if not (isinstance(a, ast.Constant)
                    and isinstance(a.value, str)):
                continue
            if not VERDICT.search(a.value):
                continue
            # a marker on the line, or on the two lines above it
            ctx = "\n".join(lines[max(0, a.lineno - 3): a.lineno])
            if OKMARK.search(ctx):
                continue
            out.append((a.lineno, a.value.strip().replace("\n", " ")[:70]))
    return out


def selftest():
    """Show the check can fail, and can pass, on synthetic sources."""
    import tempfile
    cases = [
        ("composed verdict", 'print("the model SURVIVES")', True),
        ("computed verdict",
         'v = "a" \nprint(f"the model {v}")', False),
        ("verdict in a branch",
         'if x:\n    print("REJECTED")', False),
        ("criterion, marked",
         '# verdict-ok: states the rule, not the outcome\n'
         'print("Pre-registered: ALIVE if r <= 0.3")', False),
        ("criterion, unmarked",
         'print("Pre-registered: ALIVE if r <= 0.3")', True),
    ]
    ok = True
    for name, src, want in cases:
        d = tempfile.mkdtemp()
        f = os.path.join(d, "t.py")
        io.open(f, "w", encoding="utf-8").write(src + "\n")
        got = bool(scan(f))
        ok &= (got == want)
        print(f"    {name:<26} {'flags' if got else 'passes':>7}   "
              f"{'as expected' if got == want else 'WRONG'}")
    print(f"    {'SELF-TEST OK' if ok else 'THE CHECK IS BROKEN'}")
    return ok


def main():
    print("(0) self-test")
    if not selftest():
        print("FAIL")
        sys.exit(1)
    print("\n(1) composed verdicts in code/")
    total = 0
    for fn in sorted(os.listdir(CODE)):
        if not fn.endswith(".py"):
            continue
        for ln, txt in scan(os.path.join(CODE, fn)):
            print(f"    {fn}:{ln}  {txt}")
            total += 1
    print(f"\n    {total} unmarked composed verdict(s)")
    if total:
        print("    Each must be made computed, put inside a branch, or")
        print("    marked `# verdict-ok: <reason>` if it states a")
        print("    criterion or a structural fact rather than an outcome.")
        print("FAIL")
        sys.exit(1)
    print("all clear")
    print("DONE")


if __name__ == "__main__":
    main()
