# -*- coding: utf-8 -*-
"""
A verifier that cannot come out false (lint_gates.py, increment 307).

WHY THIS EXISTS. Increment 285 found that `verify_all.py` -- the file
STATUS calls the CI stamp -- had no assertions and no failure path: its
criteria were characters inside output strings, `f"평균 r=... (기준
0.80)"`, and 0.80 appeared nowhere in the code. It printed a summary and
exited 0 whatever came out. Increment 272 had found the same shape
elsewhere and named it hazard 6's third form.

It was fixed in that one file. Increment 307 looked at the siblings:

  verify_deep.py          same defect, untouched -- criteria in prose,
                          no comparison, no failure path. It also still
                          carries the sample defect 285 fixed in the
                          same pass: `// 6 * 6 + 2` forces
                          N = 2 (mod 6) and so excludes every N
                          divisible by 3. A file called verify_DEEP
                          structurally excludes the deep N.
  verify_propositions.py  prints "ALL PASS" or "SOMETHING FAILED", and
                          "SENSITIVITY OK" or "A CHECK CANNOT FAIL",
                          and exits 0 either way. This is the file
                          CLOSURE_REAUDIT #61 holds up as the ANSWER to
                          "a check that cannot fail", and STATUS cites
                          it as a reproduction command.

So the repair was applied once and never swept. The answer that has
worked twice -- lint_docs.py for escape collapse, lint_verdicts.py for
composed conclusions -- is to make the check mechanical.

THE CHECK, and why it is this one. A first draft flagged criterion
numbers that appear only inside printed strings. It missed
verify_deep.py, its own motivating example, because that file assembles
its criteria into `rows.append(...)` and prints the rows later; and it
fired 22 times, mostly on reference values printed for the reader in
exploratory scripts. Precision matters more than reach here: a linter
that cries wolf is ignored within a week. So:

  (A) THE GATE CHECK, which fails the build. A file whose output
      ANNOUNCES a verification -- a summary banner, ALL PASS,
      SOMETHING FAILED, 전체완료, 종합검증 -- must have a failure path:
      a nonzero `sys.exit`, a `raise`, or an `assert`. Announcing a
      verdict while being unable to return the bad one is exactly the
      defect, and it is decidable from the source.

  (B) AN INFORMATIONAL LISTING, which does not fail the build.
      Criterion numbers that appear in a non-docstring string and
      nowhere in code. These are usually harmless context in a lab
      script, and occasionally the real thing; they are printed so a
      reader can judge, not counted against anyone.

Marked exemptions: `# gate-ok: <reason>` on the line.

SELF-TEST. Runs first, on synthetic sources, and must show each check
catching a planted fault and passing clean code. A verifier nobody has
run against a known-false input is not known to work (#61) -- which is
the whole subject of this file.
"""
import ast
import io
import os
import re
import sys

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

HERE = os.path.dirname(os.path.abspath(__file__))

# 배너는 **판정**을 알리는 말만이다. `전체완료`와 `DONE`은 이 코퍼스에서
# 실행이 끝났다는 표시일 뿐이라 뺀다 — 첫 판에 넣었더니 236파일 중
# 27개가 걸렸고, 전부 판정이 아니라 완료 표시였다. 정밀도가 이 검사의
# 전부다.
BANNERS = ["종합검증", "전체통과", "검증 요약", "ALL PASS",
           "SOMETHING FAILED", "ALL CHECKS PASS", "VERIFIED",
           "A CHECK CANNOT FAIL"]
WORDS = ["기준", "이어야", "반정규", "criterion", "threshold",
         "should be", "must be", "expected", "cutoff"]
NUM = re.compile(r"(?<![\w.])(\d+\.\d+|\d+)(?![\w.])")
SKIP = re.compile(r"#\s*gate-ok:")


def docstring_nodes(tree):
    out = set()
    for node in ast.walk(tree):
        if isinstance(node, (ast.Module, ast.FunctionDef,
                             ast.AsyncFunctionDef, ast.ClassDef)):
            body = getattr(node, "body", None)
            if body and isinstance(body[0], ast.Expr) and isinstance(
                    body[0].value, ast.Constant) and isinstance(
                    body[0].value.value, str):
                out.add(id(body[0].value))
    return out


def live_strings(tree):
    """String constants that are not docstrings, with their line."""
    docs = docstring_nodes(tree)
    out = []
    for node in ast.walk(tree):
        if isinstance(node, ast.Constant) and isinstance(node.value, str) \
                and id(node) not in docs:
            out.append((getattr(node, "lineno", 0), node.value))
    return out


def has_failure_path(tree):
    for node in ast.walk(tree):
        if isinstance(node, ast.Assert):
            return True
        if isinstance(node, ast.Raise):
            return True
        if isinstance(node, ast.Call) and isinstance(
                node.func, ast.Attribute) and node.func.attr == "exit":
            # 인자 없는 sys.exit()과 sys.exit(0)은 실패 경로가 아니다.
            # 첫 판에서 이 두 줄이 빠져 sys.exit(0)을 게이트로 셌고,
            # 자기검정이 코퍼스에 돌기 전에 잡았다.
            if not node.args:
                continue
            if all(isinstance(a, ast.Constant) and a.value in (0, None)
                   for a in node.args):
                continue
            return True
        if isinstance(node, ast.Call) and isinstance(
                node.func, ast.Name) and node.func.id == "SystemExit":
            return True
    return False


def code_numbers(tree):
    out = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Constant) and isinstance(
                node.value, (int, float)) and not isinstance(
                node.value, bool):
            out.add(float(node.value))
    return out


def scan_source(src):
    """Returns (announces, has_exit, prose_hits)."""
    tree = ast.parse(src)
    lines = src.splitlines()
    strs = live_strings(tree)
    announces = None
    for lineno, text in strs:
        if 1 <= lineno <= len(lines) and SKIP.search(lines[lineno - 1]):
            continue
        for b in BANNERS:
            if b in text:
                announces = (lineno, b)
                break
        if announces:
            break
    nums = code_numbers(tree)
    prose = []
    for lineno, text in strs:
        if 1 <= lineno <= len(lines) and SKIP.search(lines[lineno - 1]):
            continue
        if not any(w in text for w in WORDS):
            continue
        for m in NUM.finditer(text):
            val = float(m.group(1))
            if val not in nums:
                prose.append((lineno, text.strip()[:66], val))
    return announces, has_failure_path(tree), prose


SELFTEST = [
    ("announces, cannot fail", True, (
        "rows = [('r', '평균 0.7 (기준 0.80)')]\n"
        "print('===== 종합검증 요약 =====')\n"
        "for a, b in rows:\n"
        "    print(a, b)\n"
        "print('전체완료')\n")),
    ("announces, can fail", False, (
        "import sys\n"
        "ok = 0.7 >= 0.80\n"
        "print('===== 종합검증 요약 =====', ok)\n"
        "if not ok:\n"
        "    sys.exit(1)\n"
        "print('전체완료')\n")),
    ("announces, fails by assert", False, (
        "ok = 0.7 >= 0.80\n"
        "print('ALL PASS')\n"
        "assert ok\n")),
    ("no banner, no gate needed", False, (
        "print('mean r = 0.7, 반정규 기준 0.798 for reference')\n")),
    ("banner whitelisted", False, (
        "print('종합검증 요약')  # gate-ok: exploratory, not a gate\n")),
    ("완료 marker is not a verdict", False, (
        "print('measured r = 0.7')\nprint('전체완료')\n")),
    ("exit 0 is not a failure path", True, (
        "import sys\n"
        "print('ALL PASS')\n"
        "sys.exit(0)\n")),
]


def selftest():
    print("SELF-TEST  (A) announces a verdict but cannot return the bad one")
    ok = True
    for name, should_flag, src in SELFTEST:
        ann, hasx, _ = scan_source(src)
        got = bool(ann) and not hasx
        good = (got == should_flag)
        ok &= good
        print(f"    {name:<32} {'flags' if got else 'passes':>7}   "
              f"{'as expected' if good else 'WRONG'}")
    doc = "'''기준 0.80 in a docstring'''\nprint('hi')\n"
    _, _, prose = scan_source(doc)
    good = len(prose) == 0
    ok &= good
    print(f"    {'(B) docstring is not prose-in-code':<32} "
          f"{'passes':>7}   {'as expected' if good else 'WRONG'}")
    live = "x = 1\nprint(f'r={x} (기준 0.80)')\n"
    _, _, prose = scan_source(live)
    good = len(prose) == 1
    ok &= good
    print(f"    {'(B) live string is':<32} {'flags':>7}   "
          f"{'as expected' if good else 'WRONG'}")
    print(f"    {'SELF-TEST OK' if ok else 'SELF-TEST FAILED'}\n")
    return ok


def main():
    if not selftest():
        print("DONE (self-test failed)")
        sys.exit(1)
    files = sorted(f for f in os.listdir(HERE) if f.endswith(".py"))
    bad, prose_all = [], []
    for f in files:
        p = os.path.join(HERE, f)
        try:
            src = io.open(p, encoding="utf-8").read()
            ann, hasx, prose = scan_source(src)
        except Exception as e:
            print(f"  {f}: skipped ({type(e).__name__})")
            continue
        if ann and not hasx:
            bad.append((f, ann))
        if prose:
            prose_all.append((f, prose))

    print(f"(A) files that announce a verification but cannot fail "
          f"-- over {len(files)} files")
    if bad:
        for f, (lineno, banner) in bad:
            print(f"    {f:<28} line {lineno:>4}  prints "
                  f"“{banner}”, no failure path")
    else:
        print("    none")

    print(f"\n(B) criterion numbers that appear in no code, informational")
    n = sum(len(h) for _, h in prose_all)
    for f, hits in prose_all:
        seen = sorted({v for _, _, v in hits})
        print(f"    {f:<28} {', '.join(f'{v:g}' for v in seen)}")
    print(f"    {n} number(s) in {len(prose_all)} file(s) -- usually "
          f"context for a reader, occasionally the real thing")

    if bad:
        print(f"\n  {len(bad)} file(s) announce a verdict while being "
              f"unable to return the bad one.")
        print("DONE (findings)")
        sys.exit(1)
    print("\n  no verifier announces a verdict it cannot fail")
    print("DONE")


if __name__ == "__main__":
    main()
