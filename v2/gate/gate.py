# -*- coding: utf-8 -*-
"""
The gate. One file, fourteen checks. Exit code is the number of failures.

    python gate/gate.py > gate/gate.txt 2>&1; echo $?

Commit only on 0. Do not pipe -- a pipe eats the exit code.

CONVENTIONS THIS GATE ENFORCES

  The paper is Markdown. A numbered statement is a level-4 heading
  carrying its label, and its evidence marker is the comment under it:

      #### Theorem (w_k = 1) {#thm:A}
      <!-- evidence: analytic -->

      #### Proposition {#prop:V}
      <!-- evidence: lab_second_moment.py -->

  `analytic` means the statement is a proof with no computation behind
  it. Anything else names a script in code/, which must exist together
  with its result file in results/, and that result file must open with

      STATISTIC: <what was computed, precisely>
      FIELD: <over what range / on what set>

  A symbol whose meaning could drift is declared once:

      <!-- symbol: rho = Var C / V, the cancellation ratio -->

  Declaring the same symbol twice with different meanings fails.
"""

import hashlib
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


def evidence_paths(ev):
    """어디에 그 증거가 있고 그 결과는 어디에 있는가.

    `code/` 밑이 기본이다. 재검증 트리(`verify/passN/code/`)도 증거가 될
    수 있고 -- pass2 가 처음으로 논문 remark 를 떠받쳤다 -- 그 결과는 같은
    패스의 `results/` 에 있다. G1 과 G70 이 code/ 만 보던 동안 그 인용은
    실재하는데도 없다고 보고됐다. 규칙("인용된 증거는 실재해야 한다")은
    그대로이고 닿는 범위만 넓힌다.
    """
    stem = os.path.splitext(ev)[0] + ".txt"
    if ev.replace(os.sep, "/").startswith("verify/"):
        return (os.path.join(ROOT, ev),
                os.path.join(ROOT, stem.replace("/code/", "/results/")))
    return os.path.join(CODE, ev), os.path.join(RESULTS, stem)

STMT = r"Theorem|Proposition|Lemma|Corollary|Conjecture"
KINDS = "Theorem|Proposition|Lemma|Corollary|Conjecture|Remark"
VERDICT = r"\b(PASS|FAIL|DEAD|ALIVE|CONFIRMED|REFUTED|PROVEN|VERIFIED)\b"

fails = []
notes = []


def read(p):
    return io.open(p, encoding="utf-8", newline="", errors="replace").read()


def paper_files():
    if not os.path.isdir(PAPER):
        return []
    return sorted(os.path.join(PAPER, f) for f in os.listdir(PAPER)
                  if f.endswith(".md"))


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
    """(kind, label, evidence_or_None, line_no) for each numbered statement.

    Markdown form:  #### Theorem (optional) {#label}
                    <!-- evidence: foo.py -->
    """
    out = []
    lines = src.splitlines()
    for i, line in enumerate(lines):
        m = re.match(r"\s*#{3,4}\s+(" + STMT + r")\b[^{]*\{#([^}]*)\}",
                     line, re.I)
        if not m:
            continue
        ev = None
        for j in range(i + 1, min(i + 5, len(lines))):
            e = re.match(r"\s*<!--\s*evidence:\s*(\S+)\s*-->", lines[j])
            if e:
                ev = e.group(1)
                break
        out.append((m.group(1).lower(), m.group(2), ev, i + 1))
    return out


def cited(src):
    """(label, script, line) for EVERY evidence marker, not just the ones
    under a numbered statement.

    statements() only matches Theorem|Proposition|Lemma|Corollary|
    Conjecture, because only those are required to carry evidence. But a
    marker under a Remark is still a citation, and most of this work's
    findings are recorded in Remarks -- so G10, G12 and G14 were blind to
    exactly the places that needed them. The label is the nearest
    preceding {#...} on a heading.
    """
    out = []
    lines = src.splitlines()
    label = "?"
    for i, line in enumerate(lines):
        h = re.match(r"\s*#{1,4}\s+.*\{#([^}]*)\}", line)
        if h:
            label = h.group(1)
            continue
        e = re.match(r"\s*<!--\s*evidence:\s*(\S+)\s*-->", line)
        if e and e.group(1) != "analytic":
            out.append((label, e.group(1), i + 1))
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
            s, r = evidence_paths(ev)
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
        for m in re.finditer(r"(?<!\\)\b(" + KINDS + r")\s(\d+)", src):
            # 남의 논문 번호는 우리 카운터가 아니다. 귀속 표지를 앞에서
            # 찾되, 줄바꿈을 공백으로 눌러서 본다 -- 마크다운 변환이
            # 문장 중간에 개행을 넣기 때문에 좁은 창은 놓친다.
            back = re.sub(r"\s+", " ", src[max(0, m.start() - 44): m.start()])
            if re.search(r"(?:their|its|his|her)\s$|'s\s$|\]\s?$", back, re.I):
                continue
            fails.append(f"G2 hand-written number: {rel(path)} "
                         f"'{m.group(1)} {m.group(2)}' -- use \\ref")
            n += 1
    return n


# ------------------------------------------------------------------ G3
def g3_symbols(docs):
    seen, n = {}, 0
    for path, src in docs:
        for m in re.finditer(r"<!--\s*symbol:\s*(\S+)\s*=\s*([^>]+?)\s*-->", src):
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
    # cited(), not statements(): a result file reached only through a
    # Remark had never been checked for its header. Same blind spot as
    # G12 and G14 had.
    n, checked = 0, set()
    for path, src in docs:
        for label, ev, _ in cited(src):
            if ev in checked:
                continue
            checked.add(ev)
            _, r = evidence_paths(ev)
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
        for m in re.finditer(r"^#{3,4}[ \t]+(?:" + STMT + r")\b[^\n]*\n"
                             r"(.{80,600}?)(?=\n#{1,4}[ \t])",
                             src, re.S | re.M | re.I):
            body = re.sub(r"<!--.*?-->", "", m.group(1), flags=re.S)
            key = re.sub(r"\s+", " ", body)[:120].strip()
            if len(key) < 80:
                continue
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


# ----------------------------------------------------------------- G14
def g14_evidence_relevant(docs):
    """증거로 지목된 스크립트는 그 진술을 이름으로 알아야 한다.

    G1은 파일이 있는지, G4는 헤더가 있는지만 본다. 둘 다 통과하면서
    스크립트가 그 진술을 아예 계산하지 않을 수 있다 -- prop:scaleinv가
    lab_cell_floor.py를 가리키고 있었는데 그 스크립트는 S_2를 만들지도
    않았다. 관련성은 린트가 읽을 수 없지만 **지목**은 읽을 수 있다:
    증거로 불린 스크립트나 그 결과 파일 안에 그 라벨이 있어야 한다.
    """
    n = 0
    for path, src in docs:
        for label, ev, ln in cited(src):
            hay = ""
            for p in (os.path.join(CODE, ev),
                      os.path.join(RESULTS,
                                   os.path.splitext(ev)[0] + ".txt")):
                if os.path.exists(p):
                    hay += read(p)
            if not hay:
                continue                      # G1 reports it
            if label not in hay:
                fails.append(f"G14 {rel(path)}:{ln} names code/{ev} as "
                             f"evidence for '{label}', but neither the "
                             f"script nor its result mentions {label}")
                n += 1
    return n


# ----------------------------------------------------------------- G13
def g13_lab_nulls():
    """새 측정(`lab_`)은 자기 널을 결과 파일에 적어야 한다.

    M2("임계를 정하기 전에 널을 계산한다")와 M4("대조군을 먼저 돌린다")는
    이 리포의 방법 규칙인데 강제가 없었다. 그 틈으로 레벨 측정 하나가
    널 없이 논문에 실렸고, 뒤늦게 돌린 동전이 그걸 죽였다 -- 동전 쪽이
    오히려 더 높은 레벨까지 살아남았다. 침묵하는 고장이라 검사로 세운다.

    `audit_`는 이미 있는 주장을 대조하는 것이므로 대상이 아니다. 널이
    필요 없는 측정이라면 왜 필요 없는지를 NULL: 줄에 적으면 된다 --
    없다고 말하는 것과 안 적는 것은 다르다.
    """
    n = 0
    if not os.path.isdir(RESULTS):
        return 0
    for f in sorted(os.listdir(RESULTS)):
        if not (f.startswith("lab_") and f.endswith(".txt")):
            continue
        if not re.search(r"^NULL:", read(os.path.join(RESULTS, f)),
                         re.M):
            fails.append(f"G13 results/{f} has no 'NULL:' line -- a new "
                         f"measurement must state its control, or state "
                         f"why it needs none")
            n += 1
    return n


# ----------------------------------------------------------------- G12
VERDICT_TAIL = re.compile(
    r"^\s*(\b[A-Z]{1,2}\d+\b\s+(hold|REFUTED)\s*)*$")


def g12_refutations_disclosed(docs):
    """근거 파일 안에서 죽은 사전등록 검사는 논문이 이름으로 불러야 한다.

    사전등록한 규칙이 깨졌는데 논문이 그 스크립트를 아무 일 없다는 듯
    인용하면, 독자에게는 통과한 근거로 보인다. 스크립트를 못 대게 하는
    건 과하다 -- 깨진 검사도 정보다. 그러니 **공개**를 강제한다: 결과
    파일에서 REFUTED가 붙은 태그(U3, T6 …)는 그 스크립트를 인용한
    문서 안에 그 태그가 있어야 한다.
    """
    n = 0
    # (VERDICT_TAIL 은 아래 루프에서 판정줄과 사전등록 규칙을 가른다)
    for path, src in docs:
        named = set()
        for m in re.finditer(r"`([A-Za-z0-9_\\/]+\.py)`", src):
            named.add(os.path.basename(m.group(1).replace("\\", "")))
        for _, ev, _ in cited(src):
            named.add(os.path.basename(ev))
        for base in sorted(named):
            r = os.path.join(RESULTS, os.path.splitext(base)[0] + ".txt")
            if not os.path.exists(r):
                continue                      # G10 reports it
            dead = set()
            for line in read(r).splitlines():
                for m in re.finditer(r"REFUTED", line):
                    # 사전등록 절은 "X1 REFUTED if ..." 로 쓴다. 그것은 규칙이지
                    # 판정이 아니다. 판정줄은 REFUTED 뒤에 아무것도 없거나 다른
                    # 태그의 판정만 잇는다 -- 그 형태일 때만 죽은 것으로 센다.
                    if not VERDICT_TAIL.match(line[m.end():]):
                        continue
                    tags = re.findall(r"\b[A-Z]{1,2}\d+\b", line[:m.start()])
                    if tags:
                        dead.add(tags[-1])
            # 태그 이름공간이 논문 자체의 킬테스트 라벨(K1..K5, R1..R5)과
            # 충돌한다. 문서 어딘가에 R4가 있다는 것만으로는 공개가 아니다 --
            # 실제로 그 스크립트를 인용한 자리 근처에 있어야 한다.
            flat = src.replace("\\", "")
            spots = [m.start() for m in re.finditer(re.escape(base), flat)]
            for tag in sorted(dead):
                near = False
                for t in re.finditer(r"\b" + tag + r"\b", flat):
                    if any(abs(t.start() - s) <= 6000 for s in spots):
                        near = True
                        break
                if not near:
                    fails.append(f"G12 {rel(path)} cites code/{base}, whose "
                                 f"check {tag} is REFUTED, without naming "
                                 f"{tag} anywhere near that citation")
                    n += 1
    return n


# ----------------------------------------------------------------- G11
DEC = re.compile(r"(?<![\d.])(\d+\.\d{3,})(?![\d])")


def g11_numbers_backed(docs):
    """논문이 인쇄한 소수는 어느 결과 파일엔가 있어야 한다.

    G10은 이름만 본다. 이름이 실재해도 표의 숫자가 그 스크립트에서 나온
    숫자가 아닐 수 있다 -- 옮겨 적다 어긋나거나, 다른 파라미터로 돌린
    것이거나, 애초에 재현되지 않은 표이거나. 이번에 잡힌 것은 세 번째와
    네 번째의 혼합이었다: 표의 색인(K/d_0)이 항목을 결정하지 않아서
    독자가 확인할 수 없었다.

    규칙: 소수점 아래 세 자리 이상인 리터럴은, 어느 results/*.txt의
    숫자를 그 자릿수로 반올림한 값과 같아야 한다. 자릿수를 논문이
    정하므로 반올림은 논문 쪽 정밀도로 판정한다.
    """
    pool = []
    for top in (RESULTS, VERIFY):
        if not os.path.isdir(top):
            continue
        for base, _, fs in os.walk(top):
            if top is VERIFY and os.path.basename(base) != "results":
                continue
            for f in sorted(fs):
                if not f.endswith(".txt"):
                    continue
                for m in re.finditer(r"\d+\.\d+(?:[eE][-+]?\d+)?",
                                     read(os.path.join(base, f))):
                    try:
                        pool.append(float(m.group(0)))
                    except ValueError:
                        pass
    n = 0
    for path, src in docs:
        miss = []
        for m in DEC.finditer(src):
            lit = m.group(1)
            back = src[max(0, m.start() - 60):m.start()]
            # 1.5128\cdot10^{3} 은 가수만 떼면 결과 파일의 1.5128e+03과
            # 절대 만나지 않는다. 지수를 붙여서 하나의 값으로 본다.
            ex = re.match(r"\s*\\cdot\s*10\^\{?(-?\d+)\}?", src[m.end():])
            if ex:
                k = int(ex.group(1))
                d = len(lit.split(".")[1])
                v = float(lit) * 10.0 ** k
                tol = (0.5 * 10 ** (-d) * (1 + 1e-9) + 1e-12) * 10.0 ** k
                if not any(abs(r - v) <= tol for r in pool):
                    miss.append(lit + "e%d" % k)
                continue
            # 측정값이 아닌 것 둘. arXiv 식별자는 숫자가 아니라 이름이고,
            # 남의 논문에서 인용한 분수의 소수 전개는 우리가 잰 것이 아니다.
            if re.search(r"arXiv:\s*$", back):
                continue
            fr = re.search(r"(\d+)\s*/\s*(\d+)\s*(?:\\approx|≈)\s*$", back)
            d = len(lit.split(".")[1])
            v = float(lit)
            # 부동소수 경계에서 정확히 반 단위인 경우를 떨어뜨리지 않는다.
            tol = 0.5 * 10 ** (-d) * (1 + 1e-9) + 1e-12
            if fr and abs(int(fr.group(1)) / int(fr.group(2)) - v) <= tol:
                continue
            if not any(abs(r - v) <= tol for r in pool):
                miss.append(lit)
        if miss:
            uniq = sorted(set(miss))
            fails.append(f"G11 {rel(path)}: {len(uniq)} printed decimals "
                         f"with no result file behind them, e.g. "
                         f"{', '.join(uniq[:6])}")
            n += len(uniq)
    return n


# ----------------------------------------------------------------- G10
def g10_named_scripts(docs):
    """본문이 이름을 부른 스크립트는 실재해야 한다.

    G1은 진술에 달린 evidence 마커만 본다. 그런데 수치는 진술이 아니라
    그 옆 산문과 표에 있고, 산문은 `code/foo.py`라고 이름만 부르면 된다.
    이름이 실재하지 않아도 아무 검사도 울지 않는다 -- 재현 불가능한 표가
    조용히 남는 경로가 그것이다. 그래서 부른 이름을 전부 검사한다.
    """
    n = 0
    for path, src in docs:
        for m in re.finditer(r"`([A-Za-z0-9_\\/]+\.py)`", src):
            name = m.group(1).replace("\\", "")
            base = os.path.basename(name)
            # 이름이 재검증 트리를 가리키면 거기서 찾는다. evidence_paths 는
            # 이미 그렇게 하는데 G10 만 code/ 를 봤다 -- G1 과 G70 이 고친
            # 것과 같은 맹점이고, 실재하는 인용이 없다고 보고됐다.
            if name.replace(os.sep, "/").startswith("verify/"):
                s, r = evidence_paths(name)
                where = name
            else:
                s = os.path.join(CODE, base)
                r = os.path.join(RESULTS, os.path.splitext(base)[0] + ".txt")
                where = "code/" + base
            if not os.path.exists(s):
                fails.append(f"G10 {rel(path)} names {where}, "
                             f"which does not exist")
                n += 1
            elif not os.path.exists(r):
                fails.append(f"G10 {rel(path)} names {where} but "
                             f"{rel(r)} is missing")
                n += 1
    return n


# ------------------------------------------------------------------ G9
def g9_line_endings():
    """홀CR 금지.

    변환기가 CRLF를 안고 온 뒤 텍스트 모드로 쓰면 \\r\\r\\n이 되고, 파일
    안에 홀CR이 남는다. 그러면 splitlines()가 그것도 줄바꿈으로 세므로
    이 게이트가 보고하는 줄 번호가 편집기의 줄 번호와 어긋나고, 근거를
    가리키는 좌표가 전부 틀린다. 침묵하는 고장이라 검사로 세운다.
    """
    n = 0
    targets = paper_files()
    for d in (RESULTS, VERIFY, CODE):
        for base, _, fs in os.walk(d):
            targets += [os.path.join(base, f) for f in sorted(fs)
                        if f.endswith((".txt", ".md", ".py"))]
    for p in targets:
        raw = io.open(p, "rb").read()
        stray = raw.count(b"\r") - raw.count(b"\r\n")
        if stray:
            fails.append(f"G9 {rel(p)} has {stray} stray CR "
                         f"(line numbers will not match an editor)")
            n += 1
    return n


# ----------------------------------------------------------------- G15
def g15_refs_resolve(docs):
    """상호참조는 실재하는 앵커를 가리켜야 한다.

    두 건이 실제로 있었다. 하나는 존재하지 않는 절을 가리키는
    `\\S[sec:design]`이었고, 하나는 굵게 표시가 대괄호 안으로 밀려들어간
    `[thm:A.**.]`이었다. 둘 다 어느 검사도 울리지 않았고 손으로 grep을
    해서야 나왔다. 끊긴 참조는 조용한 고장이다 -- 렌더러는 대괄호를
    그대로 뱉고, 읽는 사람은 근거가 있다고 믿는다.

    앵커는 제목의 {#label}과 수식의 \\label{...} 둘 다에서 모은다.
    참조는 `label:...` 꼴만 본다: [HL], [GY] 같은 서지 약호는 앵커가
    아니고, 서지의 arXiv 식별자도 아니다.
    """
    anchors = set()
    for _, src in docs:
        anchors |= set(re.findall(r"\{#([^}\s]+)\}", src))
        anchors |= set(re.findall(r"\\label\{([^}]+)\}", src))
    n = 0
    for path, src in docs:
        body = re.sub(r"\{#[^}]*\}", "", src)
        for m in re.finditer(r"\[([A-Za-z][A-Za-z0-9]*:[^\]\s]+)\]", body):
            key = m.group(1)
            if key.startswith("arXiv:") or key in anchors:
                continue
            line = body.count("\n", 0, m.start()) + 1
            fails.append(f"G15 {rel(path)}:{line} refers to [{key}], "
                         f"which no heading or equation defines")
            n += 1
    return n


# ----------------------------------------------------------------- G16
def _defbody(src, name):
    """소스에서 최상위 def 한 덩어리를 뽑아 공백만 정규화한다."""
    m = re.search(r"^def %s\(.*?$" % re.escape(name), src, re.M)
    if not m:
        return None
    lines = src[m.start():].splitlines()
    out = [lines[0]]
    for ln in lines[1:]:
        if ln.strip() and not ln.startswith((" ", "\t")):
            break
        out.append(ln)
    return "\n".join(l.rstrip() for l in out).strip()


def g16_sieve_manifest():
    """체 구현이 갈라지면 감사는 감사가 아니다.

    code/의 스크립트 스무 개가 각자 자기 `sieves`를 들고 있다. 그것을
    독립 구현과 대조하는 audit_sieve.py는 자기가 실제로 돌려 본 구현만
    보증할 수 있는데, 실제로 이 저장소에는 서로 다른 구현이 넷 있었고
    감사는 하나만 덮고 있었다. 나머지는 검증된 적 없이 논문의 수치를
    만들고 있었고, 어떤 검사도 그 사실을 말하지 않았다 -- 스크립트끼리는
    서로 일관되기 때문이다.

    그래서 audit_sieve.py가 SIEVE_HASHES로 자기가 본 구현 목록을
    선언하게 하고, 디스크의 실제 목록과 대조한다. 새 변종이 하나라도
    나타나면 감사를 확장하기 전까지 게이트가 통과하지 않는다.
    """
    aud = os.path.join(CODE, "audit_sieve.py")
    if not os.path.exists(aud):
        fails.append("G16 code/audit_sieve.py is missing, so no sieve "
                     "implementation in code/ is audited")
        return 1
    src = read(aud)
    m = re.search(r"SIEVE_HASHES\s*=\s*\{(.*?)\n\}", src, re.S)
    if not m:
        fails.append("G16 code/audit_sieve.py declares no SIEVE_HASHES "
                     "manifest")
        return 1
    declared = set(re.findall(r'"([0-9a-f]{10})"', m.group(1)))
    found = {}
    for f in sorted(os.listdir(CODE)):
        if not f.endswith(".py"):
            continue
        b = _defbody(read(os.path.join(CODE, f)), "sieves")
        if b is None:
            continue
        h = hashlib.sha1(b.encode("utf-8")).hexdigest()[:10]
        found.setdefault(h, []).append(f)
    n = 0
    for h in sorted(set(found) - declared):
        fails.append(f"G16 code/{found[h][0]} has a sieve implementation "
                     f"({h}) that audit_sieve.py does not compare "
                     f"({len(found[h])} script(s) use it)")
        n += 1
    for h in sorted(declared - set(found)):
        fails.append(f"G16 audit_sieve.py declares sieve {h}, which no "
                     f"script in code/ has -- the manifest is stale")
        n += 1
    return n


# ----------------------------------------------------------------- G17
# N에 의존하는 양인데 상수처럼 적히던 것들. 값과, 그 의존성을 실제로
# 잰 스크립트를 짝지어 둔다. 그 스크립트 안에서는 리터럴이 허용된다 --
# 거기서는 발표된 값이 귀무값으로 쓰이기 때문이다.
NDEPENDENT = [
    (r"0\.3745(?!\d)", "audit_threshold_arithmetic.py",
     "the Goldbach threshold S(N)(1-A(N))"),
]


def g17_no_typed_thresholds():
    """N에 의존하는 임계값을 리터럴로 적지 않는다.

    lab_extend_range.py는 `THR = 0.3745`를 여덟 개 N 전체에 상수로
    적용하고 있었다. 그런데 그 값은 S(N)(1-A(N))이고 N의 소인수에
    의존한다 -- audit_threshold_arithmetic.py가 같은 크기의 N에서
    0.073312부터 0.374487까지 다섯 배로 움직이는 것을 쟀다.

    스윕의 여덟 N이 전부 홀수 근기 5를 가져서 값이 실제로 같았기
    때문에 숫자 자체는 맞았고, 그래서 어떤 검사도 울지 않았다. 맞는
    숫자가 우연히 맞는 것을 게이트가 구별하지 못하면 다음 N에서
    조용히 틀린다. 그러니 계산하게 하고, 타이핑을 금지한다.
    """
    n = 0
    for f in sorted(os.listdir(CODE)):
        if not f.endswith(".py"):
            continue
        src = read(os.path.join(CODE, f))
        for pat, owner, what in NDEPENDENT:
            if f == owner:
                continue
            m = re.search(pat, src)
            if m:
                line = src.count("\n", 0, m.start()) + 1
                fails.append(
                    f"G17 code/{f}:{line} writes {what} as a literal; "
                    f"it depends on N (see code/{owner}) and must be "
                    f"computed")
                n += 1
    return n


# ----------------------------------------------------------------- G18
def g18_results_not_stale():
    """결과 파일은 그것을 만든 소스보다 새것이어야 한다.

    이번에 실제로 당했다. lab_direct_route.py가 서식 오류로 중간에
    죽었는데, 이전 실행이 남긴 results 파일이 그대로 디스크에 있었고
    최신처럼 보였다. G10은 파일이 있는지만 보고, G11은 논문의 숫자가
    거기 있는지만 보고, G14는 라벨이 불리는지만 본다 -- 셋 다 통과한다.
    낡은 숫자를 읽고 그것을 방금 잰 값으로 보고하는 경로가 그것이다.

    스크립트를 고친 뒤 다시 돌리지 않는 것도 같은 고장이고, 이쪽이 더
    흔하다. 어느 쪽이든 저장소는 '이 소스가 만든 결과'와 '그 뒤에 고친
    소스'를 구별하지 못한다. mtime으로 구별한다.
    """
    n = 0
    for f in sorted(os.listdir(CODE)):
        if not f.endswith(".py") or f.startswith("_"):
            continue
        src = os.path.join(CODE, f)
        res = os.path.join(RESULTS, f[:-3] + ".txt")
        if not os.path.exists(res):
            continue          # G10이 볼 몫이다
        d = os.path.getmtime(src) - os.path.getmtime(res)
        if d > 1.0:
            fails.append(f"G18 results/{f[:-3]}.txt is {d/60:.0f} min "
                         f"older than code/{f}; rerun it or the numbers "
                         f"are not the ones that source produces")
            n += 1
    return n


# ----------------------------------------------------------------- G19
SAYLIT = re.compile(r'say\(\s*("(?:[^"\\]|\\.)*")', re.S)
SAYDEC = re.compile(r"(?<![\d.%])\d+\.\d{3,}")
# 리터럴 옆에 이 말이 있으면 그 숫자는 잰 값이 아니라 참조값이다.
REFWORD = re.compile(r"\bpub\b|published|printed|\btol\b|\bcap\b|target",
                     re.I)


def g19_no_typed_measurements():
    """잰 값을 출력 산문에 손으로 적지 않는다.

    G2는 논문에만 걸린다. 그런데 논문이 인용하는 것은 results 파일이고,
    스크립트는 자기 출력 산문에 아무 숫자나 적을 수 있다. 이번에
    lab_direct_level.py의 진단이 그랬고 -- 계산하지 않은 '4분의 1 이내,
    5% 이내'가 실제 숫자(2.18, 1.70)와 어긋난 채로 결과 파일에 들어갔다
    -- 게이트 열여덟 검사가 전부 통과했다. 전수 조사하니 다른 데도
    있었다: lab_extend_range.py의 '1.027, 1.068, 1.109'와
    lab_cell_singular.py의 표본오차·백분율이 전부 손으로 적힌 파생값이다.
    다시 돌리면 숫자가 바뀌는데 산문은 안 바뀐다.

    참조값은 예외다. audit_ 스크립트의 일이 바로 '발표된 0.114와
    대조'하는 것이고, 허용오차·상한도 규칙이지 측정이 아니다. 그래서
    같은 리터럴 안에 published/printed/pub/tol/cap/target 중 하나가
    있으면 통과시킨다.
    """
    n = 0
    for f in sorted(os.listdir(CODE)):
        if not f.endswith(".py"):
            continue
        src = read(os.path.join(CODE, f))
        body = src.split('"""', 2)[-1] if src.count('"""') >= 2 else src
        for m in SAYLIT.finditer(body):
            lit = m.group(1)
            if REFWORD.search(lit):
                continue
            d = SAYDEC.search(lit)
            if d:
                line = body.count("\n", 0, m.start()) + 1
                fails.append(
                    f"G19 code/{f}: a say() literal types {d.group(0)} "
                    f"with no published/tol/cap label, so it is a "
                    f"measurement written by hand: "
                    f"{' '.join(lit.split())[:60]}")
                n += 1
    return n


# ----------------------------------------------------------------- G20
EULER = re.compile(r"\*=\s*1\.0 - 1\.0 / \(p - 1\.0\) \*\* 2"
                   r"|\*=\s*1\.0 - 1\.0 / \(p \* \(p - 1\.0\)\)")
CBOUND = 4_000_000


def _resolve_int(src, name):
    """모듈 상수든 함수 기본인자든, 이름이 가리키는 정수 리터럴."""
    # 꼬리 주석을 허용한다 -- 상수 줄에 왜 그 값인지 적어 두는 것이
    # 이 저장소의 관례라, 이것을 놓치면 검사가 자기 관례에 걸린다.
    m = re.search(r"^\s*%s\s*=\s*([\d_]+)\s*(?:#.*)?$"
                  % re.escape(name), src, re.M)
    if m:
        return int(m.group(1).replace("_", ""))
    m = re.search(r"\b%s\s*=\s*([\d_]+)\s*[,)]" % re.escape(name), src)
    if m:
        return int(m.group(1).replace("_", ""))
    return None


def g20_euler_bound():
    """특이급수 상수는 측정 범위가 아니라 고정 경계에서 만든다.

    audit_constants.py가 잰 것: 곱을 1e5에서 끊으면 골드바흐 임계값이
    0.374486, 2e5 이상이면 0.374487이다. 즉 인쇄되는 마지막 자리가
    그 스크립트가 마침 체질한 범위에 달려 있다. 아홉 구현 중 일곱이
    상수를 자기 측정용 소수 목록 위에서 만들고 있었으므로, 다음에
    작은 N을 재는 스크립트가 하나 생기면 논문 안에서 같은 임계값이
    두 값으로 인쇄된다 -- G11은 둘 다 통과시킨다. 어느 결과 파일엔가
    있기만 하면 되기 때문이다.

    그래서 오일러 곱의 반복 대상을 고정 상수로 못박는다. 경계 4e6에서
    꼬리는 1.6e-8이고, 여섯째 자리까지 예순 배 여유가 있다.
    """
    n = 0
    for f in sorted(os.listdir(CODE)):
        if not f.endswith(".py") or f == "audit_constants.py":
            continue
        src = read(os.path.join(CODE, f))
        lines = src.splitlines()
        for i, ln in enumerate(lines):
            if not EULER.search(ln):
                continue
            it, at = None, i
            for j in range(i, max(-1, i - 8), -1):
                m = re.match(r"\s*for\s+p\s+in\s+(.*):\s*$", lines[j])
                if m:
                    it, at = m.group(1).strip(), j
                    break
            if it is None:
                continue
            it = re.sub(r"\[[^\]]*\]$", "", it).strip()
            m = re.match(r"primes_upto\(\s*(\w+)\s*\)$", it)
            if not m:
                # 변수라면 루프 바로 위의 배정을 따라간다 -- 파일 첫
                # 배정이 아니라. sieves(n)의 본문에 있는 pr = primes_upto(n)
                # 을 잡으면 바운드가 함수 인자로 읽혀 해석이 어긋난다.
                pat = re.compile(
                    r"^\s*(?:\w+\s*,\s*)*%s(?:\s*,\s*\w+)*\s*=\s*"
                    r"(?:primes_upto|sieves)\(\s*(\w+)\s*\)"
                    % re.escape(it))
                for j in range(at, -1, -1):
                    a = pat.match(lines[j])
                    if a:
                        m = a
                        break
            if not m:
                fails.append(f"G20 code/{f}:{i+1} builds an Euler product "
                             f"over '{it}', which is not a fixed bound")
                n += 1
                break
            v = _resolve_int(src, m.group(1))
            if v is None or v < CBOUND:
                fails.append(f"G20 code/{f}:{i+1} builds an Euler product "
                             f"at bound {m.group(1)}={v}, under the "
                             f"{CBOUND} the printed precision needs")
                n += 1
            break
    return n


# ----------------------------------------------------------------- G21
AGREE = re.compile(r"^AGREE\s+(\S+)\s+(\S+)\s+(-?[\d.eE+-]+)\s+"
                   r"([\d.eE+-]+)\s*$", re.M)


def g21_cross_check():
    """두 스크립트가 같은 양을 재면 값이 맞아야 한다.

    lab_signed_level.py는 A(N;k)를 담아 놓고 그것을 H(N;k)로 썼다.
    [eq:dilate]가 H = mu(k)A이므로 mu(k) 하나가 통째로 빠진 것이고,
    그 결과가 논문에 실렸다. 게이트 스무 검사가 전부 통과했다 --
    스크립트 하나 안에서는 아무것도 모순되지 않기 때문이다. 그것을
    잡은 것은 lab_direct_identity.py가 [eq:untrunc]에서 같은 양에
    독립적으로 도달한 것뿐이었다.

    그래서 그 대조를 손이 아니라 게이트가 하게 한다. 결과 파일이

        AGREE <라벨> <키> <값> <상대허용오차>

    를 찍으면, 같은 (라벨, 키)를 찍은 모든 파일의 값이 그중 가장 좁은
    허용오차 안에 있어야 한다. 한 파일에만 있는 라벨은 실패로 센다 --
    대조하는 것이 없는 대조는 확신만 주고 아무것도 검사하지 않는다.
    """
    seen = {}
    for f in sorted(os.listdir(RESULTS)):
        if not f.endswith(".txt"):
            continue
        for m in AGREE.finditer(read(os.path.join(RESULTS, f))):
            lab, key, val, tol = m.groups()
            try:
                seen.setdefault((lab, key), []).append(
                    (f, float(val), float(tol)))
            except ValueError:
                fails.append(f"G21 results/{f} has an unparsable AGREE "
                             f"line for {lab} {key}")
    n = 0
    # 라벨 단위로 본다. 두 스크립트의 N 범위가 달라 한쪽에만 있는 키가
    # 생기는 것은 정상이고, 그 키는 비교하지 않으면 그만이다. 보증은
    # 라벨마다 실제로 비교된 키가 하나 이상 있어야 한다는 것으로 지킨다.
    compared = {}
    for (lab, key), rows in sorted(seen.items()):
        compared.setdefault(lab, 0)
        if len(rows) < 2:
            continue
        compared[lab] += 1
        tol = min(r[2] for r in rows)
        lo = min(r[1] for r in rows)
        hi = max(r[1] for r in rows)
        d = abs(hi - lo) / max(abs(hi), abs(lo), 1e-300)
        if d > tol:
            who = ", ".join("%s=%.6g" % (r[0], r[1]) for r in rows)
            fails.append(f"G21 {lab} {key} disagrees by {d:.3g} "
                         f"(tol {tol:g}): {who}")
            n += 1
    for lab, c in sorted(compared.items()):
        if c == 0:
            who = sorted({r[0] for (l, k), rs in seen.items()
                          if l == lab for r in rs})
            fails.append(f"G21 {lab} has no key that two result files "
                         f"both report ({', '.join(who)}), so it "
                         f"cross-checks nothing")
            n += 1
    return n


# ----------------------------------------------------------------- G22
def g22_consumer_order():
    """남의 결과 파일을 읽는 스크립트는 그 파일보다 나중에 돌아야 한다.

    G18은 스크립트와 자기 결과 파일만 본다. 그런데 결과 파일을 읽는
    스크립트가 생기면 의존이 하나 더 늘고, 그건 보이지 않는다. 실제로
    걸렸다: audit_sieve.py는 lab_extend_range.txt의 발표표와 자기가
    독립 계산한 B(N)/N을 대조하는데, 그 표가 42분 뒤에 다시 생성됐다.
    즉 W4는 지금 디스크에 없는 판본과 비교한 결과였고, G18을 포함해
    스물한 검사가 전부 통과했다.

    의존은 소스 안의 "<이름>.txt" 리터럴로 잡는다. 그 파일이 results/에
    있고 자기 것이 아니면 소비다.
    """
    n = 0
    for f in sorted(os.listdir(CODE)):
        if not f.endswith(".py"):
            continue
        src = read(os.path.join(CODE, f))
        own = os.path.join(RESULTS, f[:-3] + ".txt")
        if not os.path.exists(own):
            continue
        for name in sorted(set(re.findall(r'"([A-Za-z0-9_]+)\.txt"', src))):
            if name == f[:-3]:
                continue
            dep = os.path.join(RESULTS, name + ".txt")
            if not os.path.exists(dep):
                continue
            d = os.path.getmtime(dep) - os.path.getmtime(own)
            if d > 1.0:
                fails.append(
                    f"G22 code/{f} reads results/{name}.txt, which is "
                    f"{d/60:.0f} min newer than results/{f[:-3]}.txt; "
                    f"rerun it or the comparison is against a version "
                    f"that no longer exists")
                n += 1
    return n


# ----------------------------------------------------------------- G23
def g23_sources_compile():
    """code/의 스크립트는 비어 있지 않고 파싱돼야 한다.

    이번에 실제로 당했다. 패치 스크립트가 io.open(p, "w", ...)로 열면서
    잘못된 newline 인자를 줬는데, 파이썬은 인자를 검증하기 전에 파일을
    truncate한다. lab_layer_tail.py가 0바이트가 됐고, 스물두 검사 중
    울린 것은 G14 하나뿐이었다 -- 그것도 파일이 비었기 때문이 아니라
    라벨이 사라졌기 때문에, 우연히.

    빈 소스나 깨진 소스는 자기 결과 파일을 만들 수 없으므로, 그 결과를
    인용하는 논문은 그 순간 근거를 잃는다. 조용한 고장이라 검사로
    세운다. 파싱은 컴파일까지만 한다 -- 실행하지 않는다.
    """
    import ast
    n = 0
    for base, _, fs in os.walk(CODE):
        if "__pycache__" in base:
            continue
        for f in sorted(fs):
            if not f.endswith(".py"):
                continue
            p = os.path.join(base, f)
            src = read(p)
            if not src.strip():
                fails.append(f"G23 {rel(p)} is empty")
                n += 1
                continue
            try:
                ast.parse(src.replace("\r\n", "\n"), filename=f)
            except SyntaxError as e:
                fails.append(f"G23 {rel(p)} does not parse: line "
                             f"{e.lineno}: {e.msg}")
                n += 1
    return n


# ----------------------------------------------------------------- G24
FITEXP = re.compile(r"N\^\{-?\d+\.\d{3,}\}")
SWEPT = re.compile(r"^SWEPT\s+(\S+)\s+(\S+)\s+([\d.eE+-]+)\s*$", re.M)


def g24_exponents_swept():
    """적합 지수를 인쇄하는 결과 파일은 그 지수의 견고성도 계산해야 한다.

    이 저장소는 지수 둘을 이미 철회했다. Q*와 M*는 손으로 고른 허용오차로
    정의된 통계의 기울기였고, audit_truncation_exponent.py가 그 허용오차를
    쓸자 지수가 단위구간을 가로질렀다. 인쇄된 하나의 기울기는 그 자체로는
    아무 정보도 주지 않는다 -- 자유 매개변수를 움직이기 전까지는.

    그래서 규칙: `~ N^{...}` 꼴의 적합 지수를 찍는 결과 파일은

        SWEPT <이름> <움직인 매개변수> <퍼짐>

    을 하나 이상 함께 찍어야 한다. 퍼짐은 계산된 값이어야 하며, G19가
    say() 리터럴에 잰 값을 손으로 적는 것을 금지하므로 표지만 붙여
    통과시킬 수 없다.
    """
    n = 0
    for f in sorted(os.listdir(RESULTS)):
        if not f.endswith(".txt"):
            continue
        src = read(os.path.join(RESULTS, f))
        if not FITEXP.search(src):
            continue
        sw = SWEPT.findall(src)
        if not sw:
            fails.append(f"G24 results/{f} prints a fitted exponent but "
                         f"no SWEPT line: nothing was varied to test it")
            n += 1
            continue
        for name, param, val in sw:
            try:
                float(val)
            except ValueError:
                fails.append(f"G24 results/{f} has an unparsable SWEPT "
                             f"spread for {name}")
                n += 1
    return n


# ----------------------------------------------------------------- G25
FROMPY = re.compile(r"#[^\n]*\bfrom\s+([A-Za-z0-9_]+)\.py\b")


def g25_no_copied_values():
    """다른 스크립트에서 온 값은 읽는다. 손으로 옮겨 적지 않는다.

    lab_cell_singular.py는 정확한 바닥 se_c를 여섯 개 숫자로 적어 두고
    주석에 출처를 밝히고 있었다 -- `from lab_mask_placebo.py`. 그건
    어느 검사도 볼 수 없는 의존이다: G18은 스크립트와 자기 결과를
    비교하고 G22는 스크립트와 그것이 읽는 파일을 비교하는데, **타이핑된
    숫자는 읽지 않으므로 둘 다 보지 못한다**. 출처 스크립트를 다시
    돌려 값이 바뀌어도 복사본은 조용히 남는다.

    그래서 출처를 밝힌 주석이 있으면 그 파일을 실제로 읽으라고
    요구한다. 읽으면 의존이 드러나고 G22가 순서를 지킨다.
    """
    n = 0
    for f in sorted(os.listdir(CODE)):
        if not f.endswith(".py"):
            continue
        src = read(os.path.join(CODE, f))
        body = src.split('"""', 2)[-1] if src.count('"""') >= 2 else src
        for m in FROMPY.finditer(body):
            who = m.group(1)
            if who == f[:-3]:
                continue
            if not os.path.exists(os.path.join(CODE, who + ".py")):
                continue
            if who + ".txt" in src:
                continue          # it reads the file: G22 takes over
            line = body.count("\n", 0, m.start()) + 1
            # Only VALUES. A copied function body is a different thing,
            # and G16's manifest already covers the one case of it, so
            # require a numeric literal to be assigned just below the
            # comment before calling this a copied value.
            after = body.splitlines()[line:line + 6]
            if not any(re.match(r"\s*[A-Za-z_]\w*\s*=\s*[\[\(]?\s*-?\d",
                                a) for a in after):
                continue
            fails.append(f"G25 code/{f}:{line} says its value comes "
                         f"from {who}.py but never reads "
                         f"results/{who}.txt; a copied number is a "
                         f"dependency no check can see")
            n += 1
    return n


# ----------------------------------------------------------------- G26
DECLINE = re.compile(r"none is run|none applies|none would mean|"
                     r"no null|not applicable|none is needed", re.I)


def g26_declines_are_audited():
    """널을 사양한 lab_ 결과에는 실제로 널을 돌린 짝이 있어야 한다.

    G13은 `NULL:` 줄이 있는지만 본다. 그런데 그 줄이 "해당 없음"이라고
    적으면 자기 선언으로 통과한다. 저장소에는 그렇게 사양한 lab_ 결과가
    일곱 있었고, 전부 감사했다: 넷은 다른 스크립트를 가리키는 포인터였고
    **넷 다 다른 통계를 가리키고 있었다**; 하나는 조건 논변이었고 네 규칙
    중 둘에서 틀렸다; 하나는 크기 주장을 무통제로 남겼다가 동전이 아무것도
    흡수하지 못한다는 것으로 뒤늦게 세워졌고; 하나는 값을 손으로 복사하고
    있었다.

    그래서 사양은 허용하되 혼자 두지 않는다: 사양한 결과 파일의 이름이
    다른 결과 파일 본문에 나와야 한다 -- 그 파일이 대조군을 공급한 것이다.
    자기 안에서 널을 돌리면 사양이 아니게 되므로 이 검사에 걸리지 않는다.
    """
    if not os.path.isdir(RESULTS):
        return 0
    files = {f: read(os.path.join(RESULTS, f))
             for f in sorted(os.listdir(RESULTS)) if f.endswith(".txt")}
    n = 0
    for f, src in files.items():
        if not f.startswith("lab_"):
            continue
        m = re.search(r"^NULL:(.*?)^[A-Z]{4,}:", src, re.S | re.M)
        if not m or not DECLINE.search(m.group(1)):
            continue
        stem = f[:-4]
        if any(g != f and stem in s for g, s in files.items()):
            continue
        fails.append(f"G26 results/{f} declines a null and no other "
                     f"result file supplies one for it; every such "
                     f"decline this repository has audited was wrong "
                     f"about something")
        n += 1
    return n


# ----------------------------------------------------------------- G27
PREBLK = re.compile(r"PRE-REGISTERED PREDICTIONS(.*?)REFUTATION RULE",
                    re.S)
PRELAB = re.compile(r"^\s{2}([A-Z]\d)\s{2}", re.M)
VERDICT = re.compile(r"\bhold\b|\bholds\b|\bheld\b|REFUTED|\bfails\b",
                     re.I)


def g27_predictions_adjudicated():
    """사전등록한 예측은 전부 판정이 나와야 한다.

    이 저장소의 규율은 돌리기 전에 예측과 반증 규칙을 적어두는 것이다.
    그런데 적어두기만 하고 판정을 출력하지 않으면, 남은 것들 중 통과한
    것만 눈에 띄고 나머지는 조용히 사라진다 -- 사전등록의 요점이 정확히
    거기서 무너진다.

    전수 조사에서 예순 남짓한 스크립트 중 하나가 그랬다:
    audit_extraction_tradeoff.py 의 T5는 "엄격한 j < K/d0 규약에서는
    발표된 자릿수가 재현되지 않는다"고 예측하고 반증 규칙까지 적어놓고,
    결과 파일에는 판정을 내지 않았다. 실제로는 여섯 열이 전부 그 규약에서
    재현되어 T5는 반증돼 있었고, FINDING 산문은 자기가 센 0과 어긋나는
    문장을 적고 있었다.

    검사는 순수하게 구조적이다: docstring 의 PRE-REGISTERED 블록에서
    라벨을 읽고, 같은 이름의 results 파일에서 그 라벨이 언급된 자리부터
    다음 라벨이 언급되는 자리까지 안에 판정어가 있는지 본다. 임계값도
    어휘 조율도 없다. 여러 라벨을 한 제목에 묶어 쓴 "V2/V3" 같은 경우는
    한 블록으로 본다.
    """
    if not os.path.isdir(RESULTS):
        return 0
    n = 0
    for f in sorted(os.listdir(CODE)):
        if not f.endswith(".py"):
            continue
        m = PREBLK.search(read(os.path.join(CODE, f)))
        if not m:
            continue
        labs = sorted(set(PRELAB.findall(m.group(1))))
        if not labs:
            continue
        rp = os.path.join(RESULTS, f[:-3] + ".txt")
        if not os.path.exists(rp):
            continue          # G10 and G18 own the missing-results case
        out = read(rp)
        body = out[out.index("FIELD:"):] if "FIELD:" in out else out
        raw = sorted(p.start() for L in labs
                     for p in re.finditer(r"\b" + L + r"\b", body))
        marks = [q for i, q in enumerate(raw)
                 if i == 0 or q - raw[i - 1] > 6]
        for L in labs:
            ok = False
            for p in re.finditer(r"\b" + L + r"\b", body):
                nxt = [q for q in marks if q > p.start()]
                end = nxt[0] if nxt else len(body)
                if VERDICT.search(body[p.start():end]):
                    ok = True
                    break
            if not ok:
                fails.append(
                    f"G27 code/{f}: {L} is pre-registered with a "
                    f"refutation rule but results/{f[:-3]}.txt never "
                    f"says whether it holds or is refuted")
                n += 1
    return n


# ----------------------------------------------------------------- G28
FORECAST = re.compile(r"crosses 1 at N =|extrapolation to N =|"
                      r"first at N =|reaching [\d.]+: N =", re.I)
BRACKET = re.compile(r"^\s*BRACKET\s+\S+\s+[\d.eE+-]+\s+[\d.eE+-]+"
                     r"\s+[\d.eE+-]+\s*$", re.M)


def g28_forecasts_carry_a_bracket():
    """계산 범위 밖으로 나가는 예보에는 구간이 붙어야 한다.

    `rem:forecast`는 `N = 2.077e8` 하나를 인용했다. 그 모형은 보정을
    적합하는 모형이고(`rem:modeltransfer`), 상수의 5% 표류가 측정된
    교차를 10% 옮긴다. 그 10%를 예보로 전파해 보니 답은 factor 2가
    아니라 **10진수 3분의 2**였다 -- 교차에서 `√K`가 지수를 두 배로
    만들고, 로그가 감쇠가 아니라 증폭을 하기 때문이다. 점추정 하나만
    적힌 예보는 그래서 독자를 오도한다.

    그러니 규약을 강제한다: 결과 파일이 계산 범위 밖의 N을 예보하면
    같은 파일이나 그 파일을 이름으로 부르는 다른 결과 파일에
    `BRACKET <이름> <점추정> <하한> <상한>` 줄이 있어야 한다. G26과
    같은 구조다 -- 자기가 못 대면 남이 대 준다.
    """
    if not os.path.isdir(RESULTS):
        return 0
    files = {f: read(os.path.join(RESULTS, f))
             for f in sorted(os.listdir(RESULTS)) if f.endswith(".txt")}
    n = 0
    for f, src in files.items():
        if not FORECAST.search(src):
            continue
        if BRACKET.search(src):
            continue
        stem = f[:-4]
        if any(g != f and stem in s and BRACKET.search(s)
               for g, s in files.items()):
            continue
        fails.append(f"G28 results/{f} forecasts an N beyond the "
                     f"computed range with no BRACKET line, in this "
                     f"file or in one that names it; the one forecast "
                     f"this repository has bracketed moved by "
                     f"two-thirds of a decade")
        n += 1
    return n


# ----------------------------------------------------------------- G29
OCTFIT = re.compile(r"^SWEPT \S+ octave-range ", re.M)
UNBOUNDED = re.compile(r"\[\s*[\d.e+]+\s*,\s*inf\s*\)")


def g29_octave_fits_are_bounded():
    """옥타브 적합에 끝이 열린 구간을 두지 않는다.

    lab_elementary_reach.py 가 `|P|`의 지수를 `H`가 닿는 범위보다 멀리
    밀어 보고 실패했는데, 실패의 이유가 이 검사의 근거다. 안쪽 길이 `L`에
    떨어지는 `(N,k)` 쌍의 수는 `1/L`로 줄어서, 네 배 멀리 갈 때마다 표본이
    네 배 사라진다. 풀링한 옥타브 개수가 `598, 668, 583, 190, 48, 8`이었고
    맨 위 두 칸이 적합을 끌어 지수를 `0.3674`로 만들었다 -- 같은 데이터의
    잘 조건화된 세 옥타브만 쓰면 `0.5178`이고 동전 대역이 `0.0251` 폭의
    `[0.4844, 0.5095]`다. 얇은 끝이 부호를 바꾼다.

    끝이 열린 구간은 그 얇은 끝을 **한 점으로 접어 숨긴다.** 가로좌표가
    무엇인지도 정해지지 않고(`rem:residue`가 이미 한 번 물린 자리다),
    범위의 맨 위라 점이 가장 적다. 그래서 `SWEPT ... octave-range`를 내는
    결과 파일에는 `[a, inf)` 꼴의 구간이 없어야 한다.
    """
    if not os.path.isdir(RESULTS):
        return 0
    n = 0
    for f in sorted(os.listdir(RESULTS)):
        if not f.endswith(".txt"):
            continue
        src = read(os.path.join(RESULTS, f))
        if not OCTFIT.search(src):
            continue
        m = UNBOUNDED.search(src)
        if not m:
            continue
        fails.append(f"G29 results/{f} fits an exponent over octaves "
                     f"with an unbounded top bin {m.group(0)}; that bin "
                     f"has no abscissa and holds the fewest points, and "
                     f"the thin end of such a fit changed 0.5178 into "
                     f"0.3674 in lab_elementary_reach")
        n += 1
    return n


# ----------------------------------------------------------------- G30
SWEPTOCT = re.compile(r"^SWEPT (\S+) octave-range ", re.M)
POP = re.compile(r"^POP (\S+) (\d+)\s*$", re.M)


def g30_octave_fits_declare_population():
    """옥타브 적합은 가장 얇은 칸의 개체수를 선언해야 한다.

    G29는 끝이 열린 칸을 막았다. 그런데 칸을 닫는 것만으로는 부족하다는
    것이 두 번 연속 드러났다: `lab_elementary_size`는 닫자마자 최상단
    칸이 `k` 한 개짜리가 되어 지수를 `0.3249`, 상관을 `0.82468`로
    끌었고, `lab_residue_size`도 같은 자리에서 같은 식으로 틀려 있었다.
    `lab_elementary_reach`는 같은 실패를 더 긴 레버에서 재어 `0.3674` 대
    잘 조건화된 `0.5178`을 얻었다. 세 번 다 원인은 하나 -- **적합에
    들어간 가장 얇은 칸이 결과 파일 어디에도 보이지 않았다.**

    그래서 `SWEPT <이름> octave-range`를 내는 파일은 같은 이름으로
    `POP <이름> <최소개체수>` 줄도 내야 한다. AGREE·BRACKET과 같은
    규약이고, 게이트가 얇은 끝을 직접 볼 수 있게 만든다.
    """
    if not os.path.isdir(RESULTS):
        return 0
    n = 0
    for f in sorted(os.listdir(RESULTS)):
        if not f.endswith(".txt"):
            continue
        src = read(os.path.join(RESULTS, f))
        have = {m.group(1) for m in POP.finditer(src)}
        want = [m.group(1) for m in SWEPTOCT.finditer(src)]
        miss = [w for w in want if w not in have]
        if not miss:
            continue
        fails.append(
            f"G30 results/{f} fits {len(want)} exponent(s) over octaves "
            f"and declares no POP line for {len(miss)} of them "
            f"({', '.join(miss[:3])}{' ...' if len(miss) > 3 else ''}); "
            f"a single-k bin has moved a fit of this kind by 0.15 in "
            f"the exponent twice in this repository")
        n += 1
    return n


# ----------------------------------------------------------------- G31
CORR = re.compile(r"^CORR (\S+) ([\d.]+)\s*$", re.M)


def g31_octave_fits_declare_correlation():
    """옥타브 적합은 자기 상관계수도 선언해야 한다.

    lab_residue_cancellation.py 의 V4가 옥타브 적합의 지수를
    -4.5321, -3.5975, -3.4549, -1.6274, +0.3328 로 내놓았다. 상관은
    -0.89576, -0.82128, -0.75874, -0.56012, +0.98266 -- 부호가 왔다
    갔다 하는, 아무 의미 없는 적합이다. 그런데 게이트의 어느 검사도
    이의를 제기하지 않았다. G29는 끝이 열린 칸을, G30은 얇은 칸을
    막지만, 둘 다 통과한 적합이 그냥 나쁠 수 있다.

    그래서 SWEPT <이름> octave-range 를 내는 파일은 같은 이름으로
    CORR <이름> <상관계수 절댓값> 도 내야 한다. POP·BRACKET·AGREE와
    같은 규약이고, 나쁜 적합을 게이트가 직접 본다.
    """
    if not os.path.isdir(RESULTS):
        return 0
    n = 0
    for f in sorted(os.listdir(RESULTS)):
        if not f.endswith(".txt"):
            continue
        src = read(os.path.join(RESULTS, f))
        have = {m.group(1) for m in CORR.finditer(src)}
        want = [m.group(1) for m in SWEPTOCT.finditer(src)]
        miss = [w for w in want if w not in have]
        if not miss:
            continue
        fails.append(
            f"G31 results/{f} fits {len(want)} exponent(s) over octaves "
            f"and declares no CORR line for {len(miss)} of them "
            f"({', '.join(miss[:3])}{' ...' if len(miss) > 3 else ''}); "
            f"a fit with correlation -0.56 has been reported as a "
            f"measurement in this repository with nothing objecting")
        n += 1
    return n


# ----------------------------------------------------------------- G32
KEXP = re.compile(r"log K\*\S*/log N")
BUDGET = re.compile(r"^BUDGET \S+ [\d.eE+-]+\s*$", re.M)


def g32_levels_declare_their_budget():
    """K* 지수를 내는 파일은 어느 예산에서 건넜는지 선언해야 한다.

    lab_split_budget.py 는 각 반쪽을 S(N)N 에서 건너게 하고 지수를
    0.7477 ... 0.7382 로 인쇄한다. 그런데 경로가 요구하는 예산은
    prop:nolog 의 S(N)(1-A(N))N 이고 4.7009 배 작다. 같은 이름
    K*_R 에 같은 모양의 지수인데 값이 0.74 대 0.56 -- 지수에서 0.18
    차이이고, theta' > 1/2 를 여유롭게 넘느냐 겨우 넘느냐가 갈린다.
    audit_residue_level.py 가 그걸 계산하기 전까지 두 수를 구분해 주는
    것은 산문뿐이었고, 결과 파일 열여섯 개가 K* 를 보고하면서 여섯은
    어느 예산인지 표시조차 없다.

    그래서 log K*/log N 을 인쇄하는 파일은 BUDGET <이름> <값> 줄로
    자기가 건넌 문턱의 상수를 선언해야 한다. AGREE·BRACKET·POP·CORR와
    같은 규약이고, 두 지수를 비교하기 전에 같은 예산인지 게이트가
    본다.
    """
    if not os.path.isdir(RESULTS):
        return 0
    n = 0
    for f in sorted(os.listdir(RESULTS)):
        if not f.endswith(".txt"):
            continue
        src = read(os.path.join(RESULTS, f))
        if not KEXP.search(src):
            continue
        if BUDGET.search(src):
            continue
        fails.append(
            f"G32 results/{f} prints a K* exponent but declares no "
            f"BUDGET line; the same K*_R reads 0.74 against S(N)N and "
            f"0.56 against S(N)(1-A(N))N, and only the second is what "
            f"prop:nolog asks for")
        n += 1
    return n


# ----------------------------------------------------------------- G33
DRIFT = re.compile(r"^DRIFT \S+ [\d.eE+-]+\s*$", re.M)


def g33_forecasts_declare_drift():
    """예보의 구간은 자기가 외삽한 상수의 표류를 선언해야 한다.

    G28은 계산 범위 밖 예보에 구간을 강제했다. 그런데 구간은 대개
    상수를 +-몇 퍼센트 흔들어 만든다 -- 그 상수가 실제로 얼마나
    표류하는지 재지 않은 채로. 흔든 폭보다 실제 표류가 크면 구간은
    좁고, 좁은 구간은 없는 구간보다 나쁘다: 정밀해 보이기 때문이다.

    audit_residue_constant.py 가 그 경계에서 멈췄다. 잔여의 상수를
    재니 c_R/sqrt(log N) 이 0.4558 ... 0.4958 로 평균의 14.36% 를
    오갔다 -- H 의 2.95% 에 대해. 사전등록이 예보를 Q2 통과에 걸어
    두었으므로 스크립트가 스스로 예보를 거부했다. 만약 걸어 두지
    않았다면 +-10% 로 만든 구간이 실제 표류의 절반짜리였을 것이다.

    그래서 BRACKET 을 내는 파일은 DRIFT <이름> <상대 퍼짐> 도 내야
    한다. 그 둘을 나란히 놓으면 구간이 자기가 외삽한 것의 흔들림을
    감당하는지 독자가 바로 본다.
    """
    if not os.path.isdir(RESULTS):
        return 0
    n = 0
    for f in sorted(os.listdir(RESULTS)):
        if not f.endswith(".txt"):
            continue
        src = read(os.path.join(RESULTS, f))
        if not BRACKET.search(src):
            continue
        if DRIFT.search(src):
            continue
        fails.append(
            f"G33 results/{f} publishes a BRACKET but declares no "
            f"DRIFT line for the constant it extrapolated; a bracket "
            f"built on a wobble smaller than the constant's own drift "
            f"is narrower than the truth and reads as precision")
        n += 1
    return n


# ----------------------------------------------------------------- G34
RADICALS = re.compile(r"^RADICALS \d+\s*$", re.M)


def g34_levels_declare_their_arithmetic():
    """K* 지수를 내는 파일은 자기 스윕의 홀수 근기 수를 선언해야 한다.

    G32는 어느 문턱에서 건넜는지를 강제한다. 그런데 문턱 자체가
    N의 산술에 따라 다섯 배 움직인다 -- audit_threshold_arithmetic.py
    가 S(N)(1-A(N))을 0.073312에서 0.374487까지 재놓았다. 그러면
    지수도 따라 움직인다.

    R에 대한 측정이 전부 N = 2e5*2^j 위에서 이뤄졌고, 그 다섯은 모두
    2^a 5^b라 홀수 근기가 하나다. lab_elementary_provable.py가 우연히
    그걸 드러냈다 -- 밀도 인자 d_L이 모든 N에서 소수 넷째 자리까지
    같았다. 허용 k-집합이 아예 안 바뀌기 때문이다. 그 위에서 잰
    `rem:residuelevel`의 지수 0.5654 ... 0.5799는 1/2을 0.06 넘겼는데,
    산술형을 바꾸니 원시근사곱 N 둘에서 0.4808과 0.4747로 **1/2 아래로
    떨어진다**. 조건부 환원이 가장 어려운 자리에서 깨지는 것이다.

    그래서 log K*/log N 을 내는 파일은 RADICALS <개수> 도 내야 한다.
    BUDGET이 어느 문턱인지를 말한다면 이건 어느 산술 위에서인지를
    말한다. 지수는 둘 다 없으면 아무 뜻이 없다.
    """
    if not os.path.isdir(RESULTS):
        return 0
    n = 0
    for f in sorted(os.listdir(RESULTS)):
        if not f.endswith(".txt"):
            continue
        src = read(os.path.join(RESULTS, f))
        if not KEXP.search(src):
            continue
        if RADICALS.search(src):
            continue
        fails.append(
            f"G34 results/{f} prints a K* exponent but declares no "
            f"RADICALS line; the same measurement reads 0.5675 over "
            f"one odd radical and 0.4747 over another, on either side "
            f"of the square-root barrier")
        n += 1
    return n


# ----------------------------------------------------------------- G35
DRIFTSLOPE = re.compile(r"^DRIFT (\S*slope\S*) [\d.eE+-]+\s*$", re.M)
SCATTER = re.compile(r"^SCATTER (\S+) [\d.eE+-]+\s*$", re.M)


def g35_slope_forecasts_declare_scatter():
    """기울기로 만든 예보는 적합의 잔차 산포도 선언해야 한다.

    G33은 외삽하는 상수의 표류를 강제한다. 그런데 그것으로 부족하다는
    것이 이 저장소 최초의 '시험받고 깨진 구간'으로 드러났다.

    lab_primorial_ladder.py 가 사다리 지수의 기울기 leave-one-out
    퍼짐으로 구간을 만들어 `1/2` 도달을 10^7.10 [10^7.07, 10^7.36]로
    예보했다. audit_primorial_reach.py 가 세 가로대를 더 올려 시험하니
    **기울기는 옳았다** -- 일곱에서 +0.006623, 열에서 +0.006643, 상관
    0.92378에서 0.95981로 좋아지고 LOO 퍼짐은 0.001827에서 0.000195로
    줄었다. 그런데 지수는 0.4941에서 멈췄고 교차는 사다리 안에 없었다.

    빠진 것은 가로대들이 자기 직선에서 얼마나 떨어져 있는가다. 잔차의
    r.m.s.가 0.0039이고 추세가 한 번 두 배에 0.0046이니, 산포가 한
    가로대분 추세의 여덟 할이다. 수준은 직선이 말하는 곳보다 몇 가로대
    앞이나 뒤에서 건너지고, 기울기가 아무리 정밀해도 그건 모른다.

    그래서 DRIFT 라벨에 slope 가 들어가면 SCATTER 도 내야 한다.
    """
    if not os.path.isdir(RESULTS):
        return 0
    n = 0
    for f in sorted(os.listdir(RESULTS)):
        if not f.endswith(".txt"):
            continue
        src = read(os.path.join(RESULTS, f))
        want = [m.group(1) for m in DRIFTSLOPE.finditer(src)]
        if not want:
            continue
        have = {m.group(1) for m in SCATTER.finditer(src)}
        miss = [w for w in want if w not in have]
        if not miss:
            continue
        fails.append(
            f"G35 results/{f} brackets a forecast by a slope's drift "
            f"({', '.join(miss)}) and declares no SCATTER; the one "
            f"bracket built that way and then tested failed, with the "
            f"slope right to 0.000195 and the rungs 0.0039 off their "
            f"own line")
        n += 1
    return n


# ----------------------------------------------------------------- G36
SHAPES = re.compile(r"^SHAPES \d+\s*$", re.M)


def g36_slope_forecasts_declare_shapes():
    """기울기로 만든 예보는 몇 가지 모양을 비교했는지 선언해야 한다.

    G33은 외삽하는 상수의 표류를, G35는 적합 주위의 산포를 강제한다.
    둘 다 **모양이 옳다는 전제** 위에서의 불확실성이다. 그 전제가
    답을 정한다는 것이 사다리에서 드러났다.

    audit_primorial_rung10.py 가 theta' = 0.56 을 log10 N = 11.2680 에
    놓았는데, 그건 log N 에 대한 직선의 답이다. audit_ladder_shape.py 가
    같은 열한 점에 다섯 모양을 맞춰 보니 r.m.s. 자체의 표준오차가 23.6%
    라 셋이 배제되지 않고, 그 셋이 1/2 에는 0.2148 진수로 합의하면서
    0.56 에는 **8.3508 진수**로 갈린다. 1.15 표준오차로 겨우 배제된
    넷째는 0.56 을 10^82.58 에 놓고 지수를 0.5663 에서 영원히 멈춘다.

    그래서 기울기 표류로 구간을 내는 파일은 SHAPES <개수> 도 내야
    한다. SHAPES 1 은 결함이 아니라 사실의 선언이다 -- RADICALS 1 이
    그랬듯이, 그 한 줄이 그 구간이 무엇에 대한 구간인지 말한다.
    """
    if not os.path.isdir(RESULTS):
        return 0
    n = 0
    for f in sorted(os.listdir(RESULTS)):
        if not f.endswith(".txt"):
            continue
        src = read(os.path.join(RESULTS, f))
        if not DRIFTSLOPE.search(src):
            continue
        if SHAPES.search(src):
            continue
        fails.append(
            f"G36 results/{f} brackets a forecast by a fitted slope "
            f"and declares no SHAPES count; three shapes that fit "
            f"these data equally well put theta' = 0.56 over 8.35 "
            f"decades apart, and a fourth never reaches it")
        n += 1
    return n


# ----------------------------------------------------------------- G37
SPANCLAIM = re.compile(r"(?:span|spread) [\d.]+ +\(floor [\d.]+\)")
FLOOR = re.compile(r"^FLOOR \S+ [\d.eE+-]+\s*$", re.M)


def g37_spans_declare_their_noise_floor():
    """작은 표본의 span 을 의존성이라 부르려면 잡음 바닥을 선언해야 한다.

    audit_residue_kexponent.py 가 일곱 산술형에 걸친 |R| 의 k-지수
    span 0.0727 을 재고 L2 의 문턱 0.05 를 넘겼다. 그런데 같은 통계량이
    **산술을 고정한 채** 얼마나 흔들리는지는 사다리가 말해 준다 -- 한
    근기의 열한 가로대가 자기 추세 주위로 r.m.s. 0.0191 이다. 그 폭에서
    일곱을 뽑으면 기대 범위가 0.0519 이고, 잰 span 은 그 1.40 배다.
    게다가 같은 근기가 두 파일에 거의 같은 N 으로 나타나는데 그 둘이
    0.0550 벌어져 있다 -- 일곱 근기 전체 span 의 4분의 3 이 한 근기
    안에서 나온다.

    그 대조 없이 L3 과 L4 의 부호를 읽었다면 산술 의존성을 발표했을
    것이다. README 의 M6 이 말하는 것("개수는 오차 막대가 아니다")이
    그건데, 강제되고 있지 않았다.

    그래서 `span x (floor y)` 꼴로 판정하는 결과 파일은
    FLOOR <이름> <값> 도 내야 한다: 같은 통계량이 조건을 고정했을 때
    내는 폭이다.
    """
    if not os.path.isdir(RESULTS):
        return 0
    n = 0
    for f in sorted(os.listdir(RESULTS)):
        if not f.endswith(".txt"):
            continue
        src = read(os.path.join(RESULTS, f))
        if not SPANCLAIM.search(src):
            continue
        if FLOOR.search(src):
            continue
        fails.append(
            f"G37 results/{f} judges a span against a floor and "
            f"declares no FLOOR line; the one span this repository "
            f"checked was 1.40 times what the same statistic gives "
            f"with the arithmetic held fixed")
        n += 1
    return n


SLOPECLAIM = re.compile(r"least-squares slope[^=\n]*= *[+-][\d.]+")
TSTAT = re.compile(r"^TSTAT (\S+) ([\d.eE+-]+)\s*$", re.M)
UNRESOLVED = re.compile(r"^UNRESOLVED SIGN (\S+)\s*$", re.M)


def g38_published_slopes_declare_significance():
    """발표한 기울기는 자기 표준오차에 대해 어디 서 있는지 말해야 한다.

    audit_residue_level.py 의 U4 는 "여유가 닫히지 않는다"를 다섯 점을
    지나는 기울기 +0.004692 로 판정하고, 근거로 leave-one-out 재적합
    셋이 모두 음이 아니라는 것을 들었다. 그건 유의성 검정이 아니다 --
    다섯 중 하나를 빼는 것은 잡음이 얼마든 최소제곱 적합의 유계 섭동이라
    세 값이 같은 부호인 건 데이터가 아니라 산술의 성질이다. 실제로
    재어 보니 그 기울기는 1.60 표준오차, 2σ 구간이 0 을 포함하고,
    leave-one-out 세 값은 그 구간 **안**에 들어 있었다. 같은 잣대로
    사다리의 수준 기울기는 11.99, k-지수는 5.03 이다.

    G35 는 예보에 쓰인 기울기의 산포를 요구하지만, 예보에 쓰이지 않고
    부호만 읽히는 기울기는 그물을 빠져나갔다. 그래서 결과 파일이
    `least-squares slope ... = <값>` 꼴로 값을 발표하면 어딘가에서
    TSTAT slope_<파일이름> <비율> 이 나와야 하고, 그 비율이 2 미만이면
    UNRESOLVED SIGN slope_<파일이름> 도 함께 나와야 한다.
    """
    if not os.path.isdir(RESULTS):
        return 0
    tstats, unresolved = {}, set()
    for f in sorted(os.listdir(RESULTS)):
        if not f.endswith(".txt"):
            continue
        src = read(os.path.join(RESULTS, f))
        for lab, val in TSTAT.findall(src):
            try:
                tstats[lab] = float(val)
            except ValueError:
                pass
        unresolved.update(UNRESOLVED.findall(src))
    n = 0
    for f in sorted(os.listdir(RESULTS)):
        if not f.endswith(".txt"):
            continue
        src = read(os.path.join(RESULTS, f))
        if not SLOPECLAIM.search(src):
            continue
        lab = "slope_" + f[:-4]
        if lab not in tstats:
            fails.append(
                f"G38 results/{f} publishes a least-squares slope and "
                f"no file declares TSTAT {lab}; the one slope this "
                f"repository checked stood at 1.60 standard errors "
                f"with a two-sigma interval containing zero")
            n += 1
        elif tstats[lab] < 2.0 and lab not in unresolved:
            fails.append(
                f"G38 results/{f} publishes a least-squares slope at "
                f"{tstats[lab]:.2f} standard errors and no file "
                f"declares UNRESOLVED SIGN {lab}; below 2 the sign is "
                f"not resolved and may not be read")
            n += 1
    return n


SPREAD = re.compile(r"^SPREAD (\S+) [\d.eE+-]+\s*$", re.M)


def g39_tstats_declare_their_range():
    """t-비율은 그것을 잰 구간 없이는 비교될 수 없다.

    G38 이 붙자 audit_slope_significance.txt 에 t 가 1.60, 11.99,
    5.03, 2.87 로 나란히 실렸고, 가족의 1.60 만 미해결로 표시됐다.
    그 표에서 빠진 건 각 계열이 몇 배의 N 에 걸쳐 있는가다 -- 가족은
    log N 폭 2.7726, 사다리는 6.9315 였다. 부호가 안 정해진 이유가
    통계량이 시끄러워서인지 쓸기가 짧아서인지는 t 만 봐서는 알 수
    없고, 둘째라면 더 계산하면 된다.

    실제로 그랬다: audit_level_slope_reach.py 가 같은 가족을 두 옥타브
    더 밀어 폭을 4.1589 로 넓히자 기울기가 +0.005112, 3.71 표준오차로
    **부호가 정해졌다** -- 새 논증이 아니라 새 구간이 답을 줬다.

    그래서 TSTAT <이름> <값> 을 내는 파일은 SPREAD <이름> <폭> 도
    내야 한다.
    """
    if not os.path.isdir(RESULTS):
        return 0
    n = 0
    for f in sorted(os.listdir(RESULTS)):
        if not f.endswith(".txt"):
            continue
        src = read(os.path.join(RESULTS, f))
        have = set(SPREAD.findall(src))
        for lab, _val in TSTAT.findall(src):
            if lab in have:
                continue
            fails.append(
                f"G39 results/{f} declares TSTAT {lab} and no SPREAD "
                f"{lab}; the one unresolved sign this repository had "
                f"was unresolved because the sweep was short, and two "
                f"more octaves settled it")
            n += 1
    return n


CROSSCLAIM = re.compile(r"AT OR (?:ABOVE|BELOW) [\d.]+")
MARGIN = re.compile(r"^MARGIN (\S+) ([\d.eE+-]+) ([\d.eE+-]+)\s*$", re.M)
INSIDEFLOOR = re.compile(r"^INSIDE FLOOR (\S+)\s*$", re.M)


def g40_crossings_are_wider_than_their_floor():
    """한 점으로 장벽을 넘었다고 선언하면 그 여유가 바닥보다 넓어야 한다.

    audit_primorial_rung10.py 는 원시근사 사다리의 수준 지수가
    N = 30750720 에서 0.5023 이라고 재고 "AT OR ABOVE 0.5" 로
    제곱근 장벽이 **관측으로** 넘겼다고 적었다. 그런데 같은 사다리가
    자기 적합선 주위로 내는 산포가 r.m.s. 0.0037 이다. 여유 0.0023 은
    그 바닥보다 **좁다** -- 한 가로대만으로는 넘은 것과 흔들린 것을
    가를 수 없었다.

    G37 은 span 에만 바닥을 요구했고, 고정 문턱에 대한 **여유**는
    그물 밖이었다. 다음 가로대 N = 61501440 을 실제로 재 보니 지수가
    0.5099, 여유 0.0099 로 바닥의 2.7 배이고 P2·P3 이 통과했다 --
    결론은 살았지만 rung 10 하나로는 정당화되지 않았었다.

    그래서 `AT OR ABOVE/BELOW <문턱>` 로 판정하는 결과 파일은
    MARGIN <파일이름> <여유> <바닥> 으로 덮여 있어야 하고, 여유가
    바닥 이하면 INSIDE FLOOR <파일이름> 도 있어야 한다.
    """
    if not os.path.isdir(RESULTS):
        return 0
    margins, inside = {}, set()
    for f in sorted(os.listdir(RESULTS)):
        if not f.endswith(".txt"):
            continue
        src = read(os.path.join(RESULTS, f))
        for lab, mg, fl in MARGIN.findall(src):
            try:
                margins[lab] = (float(mg), float(fl))
            except ValueError:
                pass
        inside.update(INSIDEFLOOR.findall(src))
    n = 0
    for f in sorted(os.listdir(RESULTS)):
        if not f.endswith(".txt"):
            continue
        src = read(os.path.join(RESULTS, f))
        if not CROSSCLAIM.search(src):
            continue
        lab = f[:-4]
        if lab not in margins:
            fails.append(
                f"G40 results/{f} declares a barrier crossed at a "
                f"single point and no file declares MARGIN {lab}; the "
                f"one such crossing this repository made was 0.0023 "
                f"against a floor of 0.0037")
            n += 1
            continue
        mg, fl = margins[lab]
        if mg <= fl and lab not in inside:
            fails.append(
                f"G40 results/{f} declares a barrier crossed by "
                f"{mg:g} against a floor of {fl:g} and no file "
                f"declares INSIDE FLOOR {lab}; a margin inside the "
                f"floor is not an observation")
            n += 1
    return n


CENSOR = re.compile(r"no crossing below k = \d+|^\s*\d+\s+>[\d.]+", re.M)
CENSORED = re.compile(r"^CENSORED (\S+) (\d+)\s*$", re.M)
UNCENSORED = re.compile(r"^UNCENSORED (\S+) (\d+)\s*$", re.M)
TRUNCBIAS = re.compile(r"^TRUNCATION BIAS (\S+)\s*$", re.M)


def g41_censored_crossings_are_declared():
    """잘린 관측은 세어서 선언해야 한다 — 잘리는 건 큰 값 쪽이다.

    lab_residue_signed.txt 는 다섯 N 중 넷의 부호 보행 교차를 싣고
    다섯째를 "none (no crossing below k = 100000)" 으로, 배수를
    ">16.9" 로 적었다. 논문은 보이는 네 값이 단조로 떨어진다고 읽고
    그걸 사슬 최대 손실의 성질로 썼다.

    그런데 k-상한이 감추는 건 정확히 **가장 큰** 교차다. 잘림이
    무작위가 아니므로 살아남은 값들로 맞춘 추세는 작은 쪽이 하는
    일 쪽으로 치우친다. beta 는 k < 100000 에서 적합하되 보행만
    k < 400000 까지 이어 다섯째를 실제로 재 보니 교차가 155333,
    이득이 +0.2181 로 **네 번째보다 올라간다** -- 단조 하락은
    잘림의 산물이었다. 기울기 자체는 -0.027526, 4.20 표준오차로
    살아남았지만 읽기는 고쳐야 했다.

    그래서 잘린 교차를 싣는 결과 파일은 CENSORED <파일이름> <개수>
    로 세어 선언해야 하고, 개수가 0 이 아니면 같은 수의
    UNCENSORED <파일이름> <개수> 로 해소되었거나
    TRUNCATION BIAS <파일이름> 로 편향이 남아 있음을 적어야 한다.
    """
    if not os.path.isdir(RESULTS):
        return 0
    counts, fixed, owned = {}, {}, set()
    for f in sorted(os.listdir(RESULTS)):
        if not f.endswith(".txt"):
            continue
        src = read(os.path.join(RESULTS, f))
        for lab, k in CENSORED.findall(src):
            counts[lab] = int(k)
        for lab, k in UNCENSORED.findall(src):
            fixed[lab] = int(k)
        owned.update(TRUNCBIAS.findall(src))
    n = 0
    for f in sorted(os.listdir(RESULTS)):
        if not f.endswith(".txt"):
            continue
        src = read(os.path.join(RESULTS, f))
        if not CENSOR.search(src):
            continue
        lab = f[:-4]
        if lab not in counts:
            fails.append(
                f"G41 results/{f} reports a crossing that was cut off "
                f"by a cap and no file declares CENSORED {lab}; the "
                f"cap hides the largest crossings, so a trend fitted "
                f"on the survivors is biased")
            n += 1
            continue
        c = counts[lab]
        if c and fixed.get(lab, 0) < c and lab not in owned:
            fails.append(
                f"G41 results/{f} has {c} censored crossing(s), "
                f"{fixed.get(lab, 0)} resolved, and no file declares "
                f"TRUNCATION BIAS {lab}; the reading that survives "
                f"censoring here reversed when the fifth point was "
                f"measured")
            n += 1
    return n


PERN = re.compile(r"^PERN (\S+) (\d+) ([\d.eE+-]+) ([\d.eE+-]+)\s*$",
                  re.M)
RATIO = re.compile(r"^RATIO (\S+) (\S+) ([\d.eE+-]+) ([\d.eE+-]+)\s*$",
                   re.M)


def g42_two_ranges_require_a_paired_ratio():
    """두 계열의 범위를 나란히 실으면 끝점끼리 나누게 된다.

    rem:residuesigned 는 "분해 전체는 θ' 에서 약 0.06, 이 한 걸음은
    0.21–0.29" 로 끝난다 -- 한쪽은 점, 한쪽은 범위다. 실제로는 둘 다
    계열이고, 같은 N 에서 짝지으면 비가 2.98–3.28 로 평평한데
    (기울기 0.06 표준오차) 끝점끼리 나누면 4.41 이 나온다. 큰 쪽
    끝과 작은 쪽 끝은 서로 다른 N 의 값이다.

    그래서 결과 파일이 PERN <이름> <점수> <최소> <최대> 를 둘 이상
    실으면 RATIO <A> <B> <최소> <최대> 도 내야 하고, 그 비는 같은 N
    에서 만들어졌어야 한다: 점수가 같아야 하고, 비의 범위가 끝점
    나눗셈이 주는 [loA/hiB, hiA/loB] 안에 **더 좁게** 들어가야 한다.
    끝점끼리 나눈 값을 그대로 RATIO 로 적으면 이 검사가 잡는다.
    """
    if not os.path.isdir(RESULTS):
        return 0
    n = 0
    for f in sorted(os.listdir(RESULTS)):
        if not f.endswith(".txt"):
            continue
        src = read(os.path.join(RESULTS, f))
        pern = {}
        for lab, cnt, lo, hi in PERN.findall(src):
            try:
                pern[lab] = (int(cnt), float(lo), float(hi))
            except ValueError:
                pass
        if len(pern) < 2:
            continue
        rats = RATIO.findall(src)
        if not rats:
            fails.append(
                f"G42 results/{f} declares {len(pern)} per-N ranges "
                f"and no RATIO line; dividing the extremes of two "
                f"ranges takes the two sides at different N, which "
                f"inflated one comparison here from 3.28 to 4.41")
            n += 1
            continue
        for A, B, lo, hi in rats:
            if A not in pern or B not in pern:
                fails.append(
                    f"G42 results/{f} declares RATIO {A} {B} and one "
                    f"of them has no PERN range")
                n += 1
                continue
            (nA, loA, hiA), (nB, loB, hiB) = pern[A], pern[B]
            try:
                lo, hi = float(lo), float(hi)
            except ValueError:
                continue
            if nA != nB:
                fails.append(
                    f"G42 results/{f} pairs {A} ({nA} points) with "
                    f"{B} ({nB} points); a per-N ratio needs the same "
                    f"N on both sides")
                n += 1
                continue
            wlo, whi = loA / hiB, hiA / loB
            # tol: a matched ratio must be materially narrower than
            # the extreme-division window, not merely inside it
            if not (lo >= wlo * 0.999 and hi <= whi * 1.001
                    and (hi - lo) <= 0.90 * (whi - wlo)):
                fails.append(
                    f"G42 results/{f} declares RATIO {A} {B} "
                    f"[{lo:g}, {hi:g}], which is not materially "
                    f"narrower (cap 0.90 of the width) than the "
                    f"[{wlo:.4f}, {whi:.4f}] that dividing the "
                    f"extremes gives; it was not formed at matching N")
                n += 1
    return n


SCALESM = re.compile(r"^SCALES (\S+) (\d+)\s*$", re.M)
ONESCALE = re.compile(r"^ONE SCALE (\S+)\s*$", re.M)


def g43_spans_declare_how_many_scales():
    """한 규모에서 잰 의존성은 서 있는 것과 유한 N 인 것을 못 가른다.

    audit_residue_arithmetic.py 는 `N ≈ 1.6·10^6` 하나에서 일곱
    산술형의 수준 지수 span 0.0928 을 재고 바닥 0.0134 의 6.9 배라며
    산술 의존성을 발표했다. 그 판정은 옳다 -- 그런데 그게 `N` 이
    커지면 닫히는 것인지 서 있는 것인지는 한 규모로는 알 수 없고,
    논문의 모든 수준 결과에 붙는 "근기 하나" 단서가 유한 `N` 산물인지
    방법의 성질인지가 거기 걸려 있었다.

    같은 일곱을 `2N, 4N, 8N` 까지 따라가 보니(2 를 곱해도 홀근기는
    그대로다) span 이 0.0928 → 0.1063 → 0.1011 → 0.0981 로 **닫히지
    않는다**. 기울기 +0.001535, 0.36 표준오차 -- 평평하다. 단서는
    계산을 더 해서 없앨 수 있는 것이 아니었다.

    그래서 `span x (floor y)` 로 판정하는 결과 파일은
    SCALES <파일이름> <규모 수> 도 내야 하고, 규모가 하나뿐이면
    ONE SCALE <파일이름> 로 그 한계를 적어야 한다.
    """
    if not os.path.isdir(RESULTS):
        return 0
    scales, single = {}, set()
    for f in sorted(os.listdir(RESULTS)):
        if not f.endswith(".txt"):
            continue
        src = read(os.path.join(RESULTS, f))
        for lab, k in SCALESM.findall(src):
            scales[lab] = int(k)
        single.update(ONESCALE.findall(src))
    n = 0
    for f in sorted(os.listdir(RESULTS)):
        if not f.endswith(".txt"):
            continue
        src = read(os.path.join(RESULTS, f))
        if not SPANCLAIM.search(src):
            continue
        lab = f[:-4]
        if lab not in scales:
            fails.append(
                f"G43 results/{f} judges a span against a floor and "
                f"no file declares SCALES {lab}; the one such span "
                f"this repository followed up the scale did not "
                f"close over a factor 8 in N")
            n += 1
        elif scales[lab] < 2 and lab not in single:
            fails.append(
                f"G43 results/{f} judges a span measured at one "
                f"scale of N and no file declares ONE SCALE {lab}; "
                f"one scale cannot separate a standing dependence "
                f"from a finite-N one")
            n += 1
    return n


MEANREAD = re.compile(r"mean (\w+) \(")
SERIESM = re.compile(r"^SERIES (\S+) (\d+) ([\d.eE+-]+) ([\d.eE+-]+)"
                     r"\s*$", re.M)
FLATM = re.compile(r"^FLAT (\S+) ([\d.eE+-]+)\s*$", re.M)
DRIFTSM = re.compile(r"^DRIFTS (\S+) ([\d.eE+-]+) ([\d.eE+-]+)\s*$",
                     re.M)


def g44_quoted_means_declare_their_series():
    """평균을 상수처럼 인용하려면 그 계열이 흐르지 않음을 보여야 한다.

    audit_model_transfer.txt 는 "mean gap 0.1677 -- 예산 인자 4.7009 가
    지수에서 그만큼 든다"로 끝나고, audit_residue_level.py 가 그 값을
    읽어 가며 논문은 그걸 방법의 상수로 쓴다. 뒤에 있는 계열은
    0.1824, 0.1688, 0.1706, 0.1635, 0.1532 -- 자기 값의 10분의 1 넘게
    떨어진다. 실제로 재 보니 기울기 -0.009165, 5.04 표준오차로
    **상수가 아니고**, 인용된 0.1677 은 쓸기의 어느 N 에서도 나오지
    않는 값이다. audit_residue_arithmetic.py 는 그 평균을 모형 예측으로
    써서 +0.1084 를 비교 대상으로 삼았는데 자기 N 에서의 값은 0.1057
    이었다.

    그래서 결과 파일이 `mean <이름> <값>` 으로 요약값을 발표하면
    SERIES <이름> <점수> <최소> <최대> 로 계열을 내야 하고,
    FLAT <이름> <t> 또는 DRIFTS <이름> <기울기> <t> 중 정확히 하나로
    그 계열이 흐르는지 판정해야 하며, 인용된 평균이 계열 범위 안에
    있어야 한다.
    """
    if not os.path.isdir(RESULTS):
        return 0
    series, flat, drifts = {}, set(), set()
    for f in sorted(os.listdir(RESULTS)):
        if not f.endswith(".txt"):
            continue
        src = read(os.path.join(RESULTS, f))
        for lab, cnt, lo, hi in SERIESM.findall(src):
            try:
                series[lab] = (int(cnt), float(lo), float(hi))
            except ValueError:
                pass
        flat.update(lab for lab, _ in FLATM.findall(src))
        drifts.update(lab for lab, _s, _t in DRIFTSM.findall(src))
    n = 0
    for path in sorted(py_files(CODE)):
        src = read(path)
        base = os.path.basename(path)
        for lab in sorted(set(MEANREAD.findall(src))):
            if lab not in series:
                fails.append(
                    f"G44 code/{base} reads a mean {lab} out of "
                    f"another results file and no file declares "
                    f"SERIES {lab}; the one such mean this repository "
                    f"checked came from a series that falls at 5.04 "
                    f"standard errors")
                n += 1
                continue
            if (lab in flat) == (lab in drifts):
                fails.append(
                    f"G44 code/{base} consumes the mean {lab} as a "
                    f"constant and its series is judged neither FLAT "
                    f"nor DRIFTS, or both; that needs one verdict")
                n += 1
    return n


ACROSS = re.compile(r"^ACROSS (\S+) ([\d.eE+-]+) ([\d.eE+-]+)\s*$",
                    re.M)
SENSITIVITY = re.compile(r"^SENSITIVITY (\S+) ([\d.eE+-]+) "
                         r"([\d.eE+-]+)\s*$", re.M)
UNDERSTATED = re.compile(r"^SENSITIVITY UNDERSTATED (\S+)\s*$", re.M)


def g45_cross_radical_spans_are_judged_against_the_declared_sweep():
    """예보의 민감도 쓸기가 무엇을 덮는지 실제로 대조해야 한다.

    lab_elementary_provable.py 는 초등 절반의 예보 10^5474.8 에 두 개의
    폭을 단다 -- 해석적 지수 c 를 쓴 [10^2092.7, 10^13093.3] 과
    A·dL 을 2.1071 배씩 흔든 [10^4838.5, 10^6139.9]. 그리고 자기 본문에
    "dL 은 여기서 전혀 흐르지 않는다 -- 이 쓸기의 모든 N 이 같은
    홀근기라서다, **이 쓸기에 대한 사실이지 일반적인 사실이 아니다**"
    라고 적어 두었다. 지목만 하고 시험하지 않은 노출이다.

    일곱 산술형에서 세 입력을 다시 재니 dL 이 0.1730 에서 0.4726 까지
    (상대 퍼짐 0.9309) 움직이고 예보가 두 규약에 걸쳐 4784.6 부터
    6184.3 까지 간다 -- 선언된 A·dL 폭 밖이다. c-쓸기 안에는 들어가므로
    결론은 안 바뀌지만, **선언된 민감도가 실제 노출을 과소하게
    적고 있었다.**

    그래서 어떤 예보의 산술형 간 폭을 재어 ACROSS <이름> <최소> <최대>
    로 발표하면 그 이름의 BRACKET 이 어딘가에 있어야 하고,
    SENSITIVITY <이름> <최소> <최대> 로 선언된 쓸기를 함께 적어야 하며,
    ACROSS 가 그 안에 들어가지 않으면
    SENSITIVITY UNDERSTATED <이름> 을 붙여야 한다.
    """
    if not os.path.isdir(RESULTS):
        return 0
    brackets, sens, under = set(), {}, set()
    for f in sorted(os.listdir(RESULTS)):
        if not f.endswith(".txt"):
            continue
        src = read(os.path.join(RESULTS, f))
        brackets.update(m.split()[1] for m in
                        re.findall(r"^BRACKET \S+ .*$", src, re.M))
        for lab, lo, hi in SENSITIVITY.findall(src):
            try:
                sens[lab] = (float(lo), float(hi))
            except ValueError:
                pass
        under.update(UNDERSTATED.findall(src))
    n = 0
    for f in sorted(os.listdir(RESULTS)):
        if not f.endswith(".txt"):
            continue
        src = read(os.path.join(RESULTS, f))
        for lab, lo, hi in ACROSS.findall(src):
            try:
                lo, hi = float(lo), float(hi)
            except ValueError:
                continue
            if lab not in brackets:
                fails.append(
                    f"G45 results/{f} declares ACROSS {lab} and no "
                    f"file publishes a BRACKET {lab}; a cross-radical "
                    f"span has nothing to be judged against")
                n += 1
                continue
            if lab not in sens:
                fails.append(
                    f"G45 results/{f} declares ACROSS {lab} and no "
                    f"SENSITIVITY {lab}; the one forecast this "
                    f"repository swept across radicals ran outside "
                    f"its own declared sweep")
                n += 1
                continue
            slo, shi = sens[lab]
            covered = slo <= lo and hi <= shi
            if covered == (lab in under):
                fails.append(
                    f"G45 results/{f} declares ACROSS {lab} "
                    f"[{lo:g}, {hi:g}] against SENSITIVITY "
                    f"[{slo:g}, {shi:g}] and the "
                    f"SENSITIVITY UNDERSTATED {lab} marker "
                    f"{'is present anyway' if covered else 'is missing'}")
                n += 1
    return n


SIGMACOUNT = re.compile(r"two standard deviations? [\d.]+")
EXCHANGE = re.compile(r"^EXCHANGE (\S+) (\d+) (\d+) (\d+)\s*$", re.M)


def g46_counts_are_judged_against_an_exchangeable_null():
    """개수를 √n 로 재려면 단위가 독립이어야 한다 -- 대개 아니다.

    audit_residue_coin_rank.py 의 V4 는 mu 의 옥타브 평균이 동전
    중앙값 아래로 가는 횟수를 세어 30 중 21, 기대 15.0, 2σ 5.5 로
    **반증**을 냈다. 그런데 그 널은 서른 옥타브가 독립이라고 가정하는데
    아니다 -- 같은 열여섯 부호 벡터를 모든 옥타브와 모든 N 에서 쓰고,
    다른 N 의 옥타브들이 겹치는 m 위를 달린다.

    교환가능 버전은 공짜다: 같은 개수를 동전 하나하나에 대해 세면
    22, 24 가 나온다. 즉 21 은 17 중 3 등 사건이지 2.2σ 가 아니다.
    크기를 쓰는 짝 검정도 mu 를 17 중 4 등에 놓아 같은 말을 한다.
    V4 가 잡은 건 데이터가 아니라 자기 널이었다.

    그래서 결과 파일이 개수를 `two standard deviations <값>` 으로
    판정하면 EXCHANGE <파일이름> <개수> <순위> <총 뽑기수> 로 같은
    통계량을 교환가능한 뽑기들에 대해 센 결과를 함께 내야 한다.
    """
    if not os.path.isdir(RESULTS):
        return 0
    seen = {}
    for f in sorted(os.listdir(RESULTS)):
        if not f.endswith(".txt"):
            continue
        src = read(os.path.join(RESULTS, f))
        for lab, cnt, rank, tot in EXCHANGE.findall(src):
            seen[lab] = (int(cnt), int(rank), int(tot))
    n = 0
    for f in sorted(os.listdir(RESULTS)):
        if not f.endswith(".txt"):
            continue
        src = read(os.path.join(RESULTS, f))
        if not SIGMACOUNT.search(src):
            continue
        lab = f[:-4]
        if lab not in seen:
            fails.append(
                f"G46 results/{f} judges a count against two standard "
                f"deviations and no file declares EXCHANGE {lab}; the "
                f"one such count this repository checked was 2.2 "
                f"sigma on a binomial null and 3-in-17 on an "
                f"exchangeable one")
            n += 1
            continue
        _cnt, rank, tot = seen[lab]
        if not (1 <= rank <= tot):
            fails.append(
                f"G46 results/{f} declares EXCHANGE {lab} with rank "
                f"{rank} of {tot}, which is not a rank")
            n += 1
    return n


SQRTREF = re.compile(r"sqrt ?\(?#|sqrt of the\s+number of")
REFERENCE = re.compile(r"^REFERENCE (\S+) (\d+) ([\d.eE+-]+) "
                       r"([\d.eE+-]+)\s*$", re.M)


def g47_count_references_declare_the_magnitude_one():
    """`독립 부호면 √n` 은 크기가 같을 때의 값이다 -- 재서 써야 한다.

    lab_positive_weights.py 는 확대에 걸친 상쇄 이득 G = 1.834 … 2.789
    를 sqrt(#k) = 17.7 … 38.5 에 대고 "독립 부호면 이만큼", 즉
    n_eff = G^2 = 3.4 … 7.8 인데 항은 수백~천 수백이라고 읽는다.
    그런데 sqrt(#k) 는 크기가 **모두 같을 때**의 l1/l2 다. 여기 크기는
    같지 않고(자기 규칙 T4 가 최상위 십분위 몫을 0.3486 … 0.3587 로
    쟀다), 실제 l1/l2 는 sqrt(#k) 의 0.6622 … 0.6854 배다.

    그래서 결손이 갈라진다: sqrt(#k)/G = 9.65 … 13.82 가
    크기 집중 1.459 … 1.510 (평평) 곱하기 부호 상관 6.52 … 9.40
    (오름)이다. 상관은 실재하고 결론은 서지만, 발표된 비교는 그
    크기를 집중분만큼 부풀리고 있었다.

    그래서 `sqrt(#k)` 꼴의 개수 기준을 인쇄하는 결과 파일은
    REFERENCE <파일이름> <점수> <최소비> <최대비> 로 실제 크기 기준
    l1/l2 가 그 개수 기준의 몇 배인지도 내야 한다. 그 비는 코시-슈바르츠로
    1 을 넘을 수 없으므로 넘으면 계산이 틀린 것이다.
    """
    if not os.path.isdir(RESULTS):
        return 0
    seen = {}
    for f in sorted(os.listdir(RESULTS)):
        if not f.endswith(".txt"):
            continue
        src = read(os.path.join(RESULTS, f))
        for lab, cnt, lo, hi in REFERENCE.findall(src):
            try:
                seen[lab] = (int(cnt), float(lo), float(hi))
            except ValueError:
                pass
    n = 0
    for f in sorted(os.listdir(RESULTS)):
        if not f.endswith(".txt"):
            continue
        src = read(os.path.join(RESULTS, f))
        if not SQRTREF.search(src):
            continue
        lab = f[:-4]
        if lab not in seen:
            fails.append(
                f"G47 results/{f} prints a count-based reference and "
                f"no file declares REFERENCE {lab}; the one such "
                f"reference this repository checked was 1.46 to 1.51 "
                f"times the magnitude-based one")
            n += 1
            continue
        cnt, lo, hi = seen[lab]
        if not (cnt >= 1 and 0.0 < lo <= hi <= 1.0 + 1e-9):
            fails.append(
                f"G47 results/{f} declares REFERENCE {lab} with ratio "
                f"range [{lo:g}, {hi:g}] over {cnt} points; l1/l2 "
                f"cannot exceed sqrt(n), so a ratio above 1 is an "
                f"arithmetic error")
            n += 1
    return n


REFLEVEL = re.compile(r"reference level")
FITTEDDECAY = re.compile(r"~ N\^\{-|~ \(log N\)\^\{-")
FLOORTREND = re.compile(r"^FLOORTREND (\S+) ([-+][\d.eE+-]+) "
                        r"([\d.eE+-]+)\s*$", re.M)
FLOORDELEG = re.compile(r"^FLOOR DELEGATED (\S+) (\S+)\s*$", re.M)


def g48_reference_levels_declare_their_own_trend():
    """기준선을 중심만 보고 추세를 안 보면 바닥의 움직임을 μ 것으로 읽는다.

    lab_lean_decay.py 는 질량가중 부호 기울기 |0.5-f| 가
    N^{-0.1673} 로 준다고 재고 "기울기는 μ 의 구조적 사실이 아니라
    thm:C 의 유한 N 오차이고 사라진다"고 읽는다. 곁에 실은 동전 팔은
    N 당 두 뽑기이고 "동전은 내내 1/2 에 앉는다, 그래야 하듯이"로만
    읽혔다. 1/2 에 앉느냐가 물음이 아니다 -- 같은 크기 위 무작위 부호
    마당도 자기 기울기 l2/(2 l1) 를 갖고, 그 **바닥이 N 과 함께
    움직인다**.

    실제로 재니 바닥의 기울기가 -0.315933 (17.68 표준오차)로 μ 의
    -0.167257 보다 **가파르다** -- 차이가 5.77 표준오차다. 바닥으로
    나눈 기울기는 9.86 에서 17.14 로 **오른다** (기울기 +0.148676,
    5.81 표준오차). 즉 줄어든 건 기울기가 아니라 바닥이었다.

    그래서 어떤 양의 감쇠를 `~ N^{-...}` 꼴로 적합하면서 널 팔을
    reference level 로 부르는 결과 파일은
    FLOORTREND <파일이름> <기울기> <t> 로 그 기준선 자신의 추세를
    내야 하고, 다른 파일에 위임했다면
    FLOOR DELEGATED <파일이름> <위임처> 를 적되 위임처가 FLOORTREND
    를 갖고 있어야 한다.
    """
    if not os.path.isdir(RESULTS):
        return 0
    trends, deleg = {}, {}
    for f in sorted(os.listdir(RESULTS)):
        if not f.endswith(".txt"):
            continue
        src = read(os.path.join(RESULTS, f))
        for lab, sl, t in FLOORTREND.findall(src):
            trends[lab] = (sl, t)
        for lab, tgt in FLOORDELEG.findall(src):
            deleg[lab] = tgt
    n = 0
    for f in sorted(os.listdir(RESULTS)):
        if not f.endswith(".txt"):
            continue
        src = read(os.path.join(RESULTS, f))
        if not (REFLEVEL.search(src) and FITTEDDECAY.search(src)):
            continue
        lab = f[:-4]
        if lab in trends:
            continue
        tgt = deleg.get(lab)
        if tgt and tgt in trends:
            continue
        fails.append(
            f"G48 results/{f} fits a decay against N and calls a null "
            f"arm a reference level, and no file declares FLOORTREND "
            f"{lab} or a FLOOR DELEGATED pointing at one; the one "
            f"such floor this repository measured moved faster than "
            f"the quantity it was the floor for")
        n += 1
    return n


FLOORRANGE = re.compile(r"^FLOORRANGE (\S+) (\d+) (\d+)\s*$", re.M)


def g49_delegated_floors_cover_the_range_they_are_used_on():
    """빌려 온 바닥은 그것이 측정된 구간에서만 바닥이다.

    lab_extend_range.py 는 N = 2e5 에서 2.56e7 까지 f 를 싣고 널을
    돌리지 않으면서 "f 의 동전 기준선은 1/2 이고 거기서 쟀다"고
    lab_lean_decay.py 를 가리킨다. 그런데 그 파일은 6.4e6 에서
    멈춘다 -- 가장 큰 두 N 은 바닥을 가진 적이 없었다.

    실제로 여덟 N 전부에 대해 재니 바닥의 기울기가 -0.313205
    (24.21 표준오차)로 μ 의 -0.153911 보다 9.29 표준오차만큼 가파르고,
    바닥으로 나눈 기울기가 8.49 에서 21.36 으로 +0.159294 (9.58
    표준오차) 오른다. 여섯 점과 여덟 점의 기울기 차는 0.22 표준오차라
    {#rem:leanfloor} 의 뒤집기는 짧은 쓸기의 산물이 아니었다.

    그래서 FLOORTREND <이름> 은 FLOORRANGE <이름> <최소N> <최대N> 을
    달아야 하고, FLOOR DELEGATED <a> <b> 는 b 의 구간이 a 의 구간을
    **덮을 때만** 유효하다.
    """
    if not os.path.isdir(RESULTS):
        return 0
    trends, ranges, deleg = set(), {}, {}
    for f in sorted(os.listdir(RESULTS)):
        if not f.endswith(".txt"):
            continue
        src = read(os.path.join(RESULTS, f))
        trends.update(lab for lab, _s, _t in FLOORTREND.findall(src))
        for lab, lo, hi in FLOORRANGE.findall(src):
            ranges[lab] = (int(lo), int(hi))
        for a, b in FLOORDELEG.findall(src):
            deleg[a] = b
    n = 0
    for lab in sorted(trends):
        if lab not in ranges:
            fails.append(
                f"G49 FLOORTREND {lab} carries no FLOORRANGE; a floor "
                f"is a floor only over the N it was measured on, and "
                f"the one delegation this repository checked ran a "
                f"factor 4 past its own floor")
            n += 1
    for a, b in sorted(deleg.items()):
        if a not in ranges or b not in ranges:
            fails.append(
                f"G49 FLOOR DELEGATED {a} {b} and one of them "
                f"declares no FLOORRANGE, so the coverage cannot be "
                f"checked")
            n += 1
            continue
        (alo, ahi), (blo, bhi) = ranges[a], ranges[b]
        if not (blo <= alo and ahi <= bhi):
            fails.append(
                f"G49 FLOOR DELEGATED {a} {b}: {b} was measured over "
                f"[{blo}, {bhi}] and {a} is used over [{alo}, {ahi}], "
                f"which it does not cover")
            n += 1
    return n


WINDOW = re.compile(r"^WINDOW (\S+) (\d+) (\d+)\s*$", re.M)
EXPLAINS = re.compile(r"^EXPLAINS (\S+) (\S+)\s*$", re.M)
WINDISJ = re.compile(r"^WINDOWS DISJOINT (\S+) (\S+)\s*$", re.M)


def g50_mechanisms_declare_the_window_they_were_shown_on():
    """기전은 그것이 보여진 창에서만 기전이다.

    lab_lean_oddmertens.py 는 sign H(N;k) 가 sign Modd(floor(N/k)) 와
    0.7657 로 맞는다는 것을 **2 <= N/k <= 1000** 에서 보이고, 그걸로
    lab_lean_decay.py 의 기울기 f 의 감쇠를 설명했다고 적었다. 그런데
    f 는 k < N^0.56 위에서 재고, 그 안쪽 길이는 N/k = 215 에서
    2133333 까지다 -- 보여진 창의 위끝을 2133 배 넘어선다.

    그 창에서 실제로 재 보니 일치도가 0.7657 이 아니라 0.5201 … 0.6161
    이고, Modd 예측자를 μ 자신의 크기에 얹으면 기울기를 1.8 … 3.0 배
    과대예측하며, 기울기가 -0.004298 (0.94 표준오차)로 **평평**한데
    μ 는 -0.167257 이다 -- 8.51 표준오차 차이. 기전은 통계량이 있는
    곳에서 통계량을 나르지 않는다.

    그래서 두 창을 나란히 선언하면(WINDOW <이름> <lo> <hi>) 그 관계를
    판정해야 한다: 포함되면 EXPLAINS <a> <b>, 아니면
    WINDOWS DISJOINT <a> <b>. 선언과 수가 어긋나면 실패한다.
    """
    if not os.path.isdir(RESULTS):
        return 0
    n = 0
    for f in sorted(os.listdir(RESULTS)):
        if not f.endswith(".txt"):
            continue
        src = read(os.path.join(RESULTS, f))
        win = {}
        for lab, lo, hi in WINDOW.findall(src):
            win[lab] = (int(lo), int(hi))
        if len(win) < 2:
            continue
        pairs = EXPLAINS.findall(src) + WINDISJ.findall(src)
        if not pairs:
            fails.append(
                f"G50 results/{f} declares {len(win)} windows and "
                f"adjudicates none; the one mechanism this repository "
                f"checked was shown on a window its target overran by "
                f"a factor 2133")
            n += 1
            continue
        for a, b in EXPLAINS.findall(src):
            if a not in win or b not in win:
                fails.append(
                    f"G50 results/{f} declares EXPLAINS {a} {b} and "
                    f"one of them has no WINDOW")
                n += 1
                continue
            (alo, ahi), (blo, bhi) = win[a], win[b]
            if not (alo <= blo and bhi <= ahi):
                fails.append(
                    f"G50 results/{f} declares EXPLAINS {a} {b} but "
                    f"{b}'s window [{blo}, {bhi}] is not inside {a}'s "
                    f"[{alo}, {ahi}]")
                n += 1
        for a, b in WINDISJ.findall(src):
            if a not in win or b not in win:
                fails.append(
                    f"G50 results/{f} declares WINDOWS DISJOINT {a} "
                    f"{b} and one of them has no WINDOW")
                n += 1
                continue
            (alo, ahi), (blo, bhi) = win[a], win[b]
            if alo <= blo and bhi <= ahi:
                fails.append(
                    f"G50 results/{f} declares WINDOWS DISJOINT {a} "
                    f"{b} but {b}'s window [{blo}, {bhi}] does sit "
                    f"inside {a}'s [{alo}, {ahi}]")
                n += 1
    return n


CARRIES = re.compile(r"^CARRIES (\S+) (\S+) ([\d.eE+-]+) "
                     r"([\d.eE+-]+)\s*$", re.M)


def g51_disjoint_windows_are_adjudicated_on_the_target():
    """창이 어긋난다고 적는 것으로는 부족하다 -- 표적 위에서 재야 한다.

    G50 이 두 창의 포함관계를 판정하게 만들자 lab_lean_oddmertens 와
    lab_survivor_selection 이 둘 다 WINDOWS DISJOINT 로 나왔다. 그건
    공개일 뿐 판정이 아니다. 표적 창 위에서 실제로 재 보면 두
    예측자가 서로 다르게 실패한다 -- Modd 는 기울기를 1.73 에서 2.96
    배 과대예측하고 감쇠가 -0.004298 로 평평하며, 사체가중 P 는
    0.6571 에서 0.8676 배로 과소예측하고 -0.230633 으로 μ 의
    -0.167260 보다 2.11 표준오차 빠르게 준다. 일치도도 0.52-0.62 대
    0.74-0.81 로 갈린다. **어느 쪽도 기울기가 재어지는 곳에서 그것을
    나르지 않고, P 가 훨씬 가깝다.**

    그래서 WINDOWS DISJOINT <a> <b> 를 내면 표적 창 위에서 실제로 잰
    CARRIES <a> <b> <최소비> <최대비> 도 있어야 한다. 공개만으로
    통과시키지 않는다.
    """
    if not os.path.isdir(RESULTS):
        return 0
    carries, wins = {}, {}
    for f in sorted(os.listdir(RESULTS)):
        if not f.endswith(".txt"):
            continue
        src = read(os.path.join(RESULTS, f))
        for a, b, lo, hi in CARRIES.findall(src):
            try:
                carries[(a, b)] = (float(lo), float(hi))
            except ValueError:
                pass
        for lab, lo, hi in WINDOW.findall(src):
            wins[lab] = (int(lo), int(hi))
    n = 0
    for f in sorted(os.listdir(RESULTS)):
        if not f.endswith(".txt"):
            continue
        src = read(os.path.join(RESULTS, f))
        for a, b in WINDISJ.findall(src):
            if (a, b) not in carries:
                fails.append(
                    f"G51 results/{f} declares WINDOWS DISJOINT {a} "
                    f"{b} and no file declares CARRIES {a} {b}; "
                    f"disclosing that a mechanism was shown on "
                    f"another window is not measuring it on this one")
                n += 1
                continue
            lo, hi = carries[(a, b)]
            if not (0.0 < lo <= hi):
                fails.append(
                    f"G51 CARRIES {a} {b} declares [{lo:g}, {hi:g}], "
                    f"which is not a ratio range")
                n += 1
    for (a, b) in sorted(carries):
        for lab in (a, b):
            if lab not in wins:
                fails.append(
                    f"G51 CARRIES {a} {b} names {lab}, which declares "
                    f"no WINDOW")
                n += 1
    return n


BOTHPARTS = re.compile(r"move[s]? both parts together|moves both")
ONESIDED = re.compile(r"^ONESIDED (\S+) ([-+]?[\d.eE+-]+) "
                      r"([-+]?[\d.eE+-]+)\s*$", re.M)


def g52_declined_nulls_meet_a_one_sided_randomisation():
    """"둘 다 같이 움직인다"는 한쪽만 깨는 널을 배제하지 않는다.

    lab_predictable_part.py 는 분해 H = beta P + R 의 잔여 몫
    0.6310 … 0.5307 을 발표하면서 널을 사양한다 -- "하나의 잰 합을
    분해하는 것이라 무작위화는 두 부분을 함께 움직인다". 그건 둘 다
    움직이는 무작위화에 대해서만 참이다. P 만 깨고 H 는 건드리지 않는
    것이 둘 있다: (a) 부호를 두고 |P| 를 k 에 걸쳐 순열, (b) |P| 를
    두고 부호를 다시 뽑기. 각각 예측자의 절반씩을 없앤다.

    실제로 돌리니 (a) 는 잰 삭감의 10-15%, (b) 는 0% 를 낸다 -- 오히려
    음수다. 즉 37-47% 의 삭감은 부호 마당만으로도 크기만으로도 나오지
    않고 **k 별 짝맞춤**이 있어야 한다. 사양된 널이 실제로 정보를
    갖고 있었다.

    그래서 널을 "두 부분이 함께 움직인다"는 이유로 사양하는 결과
    파일은 한쪽만 깬 무작위화의 결과를
    ONESIDED <파일이름> <(a) 분율> <(b) 분율> 로 내야 한다. 깨진
    예측자가 온전한 것보다 많이 사면(분율 > 1) 계산이 틀린 것이다.
    """
    if not os.path.isdir(RESULTS):
        return 0
    seen = {}
    for f in sorted(os.listdir(RESULTS)):
        if not f.endswith(".txt"):
            continue
        src = read(os.path.join(RESULTS, f))
        for lab, a, b in ONESIDED.findall(src):
            try:
                seen[lab] = (float(a), float(b))
            except ValueError:
                pass
    n = 0
    for f in sorted(os.listdir(RESULTS)):
        if not f.endswith(".txt"):
            continue
        src = read(os.path.join(RESULTS, f))
        if not BOTHPARTS.search(src):
            continue
        lab = f[:-4]
        if lab not in seen:
            fails.append(
                f"G52 results/{f} declines a null because a "
                f"randomisation would move both parts together and no "
                f"file declares ONESIDED {lab}; the one such decline "
                f"this repository audited hid a null that reaches 15 "
                f"per cent of the effect")
            n += 1
            continue
        a, b = seen[lab]
        if a > 1.05 or b > 1.05:
            fails.append(
                f"G52 ONESIDED {lab} declares fractions "
                f"{a:g} and {b:g}; a randomisation that breaks half "
                f"the predictor cannot buy more of the effect than "
                f"the whole predictor does")
            n += 1
    return n


SHAPESURVIVE = re.compile(r"^SHAPESURVIVE (\S+) (\d+) (\d+) "
                          r"([\d.eE+-]+)\s*$", re.M)
SHAPECURRENT = re.compile(r"^SHAPECURRENT (\S+) (\d+)\s*$", re.M)


def g53_shape_survival_is_redone_on_every_point():
    """모양 판정은 점 개수의 함수다 -- 점이 늘면 다시 해야 한다.

    audit_ladder_shape.py 는 다섯 모양을 원시근사 사다리의 **열한**
    가로대에 맞추고, 최선 r.m.s. 가 자유도 9 에서 23.6% 표준오차를
    가지므로 한 표준오차 안에 세 모양이 살아남아 theta' = 0.56 을
    11.2700, 14.6167, 19.6207 -- 8.35 자릿수 벌어진 자리에 놓는다고
    적었다. 그 뒤 audit_primorial_rung11.py 가 열두 번째 가로대를
    쟀는데 그 판정은 다시 되지 않았다.

    열둘에서 다시 하니 순위는 그대로(직선이 최선)지만 생존이 셋에서
    **둘**로 줄고 폭이 8.3508 에서 2.7845 자릿수로 좁아진다. 가로대
    하나가 산 것이 그만큼이다.

    그래서 SHAPESURVIVE <표적> <점수> <생존수> <폭> 을 내면
    SHAPECURRENT <표적> <점수> 로 어느 것이 현행인지 적어야 하고,
    그 점수는 같은 표적에 대해 선언된 것 중 **최대**여야 한다.
    """
    if not os.path.isdir(RESULTS):
        return 0
    rows, cur = {}, {}
    for f in sorted(os.listdir(RESULTS)):
        if not f.endswith(".txt"):
            continue
        src = read(os.path.join(RESULTS, f))
        for tgt, pts, sv, sp in SHAPESURVIVE.findall(src):
            rows.setdefault(tgt, {})[int(pts)] = (int(sv), float(sp))
        for tgt, pts in SHAPECURRENT.findall(src):
            cur[tgt] = int(pts)
    n = 0
    for tgt in sorted(rows):
        if tgt not in cur:
            fails.append(
                f"G53 SHAPESURVIVE {tgt} is declared at "
                f"{sorted(rows[tgt])} points and no SHAPECURRENT says "
                f"which stands; a shape adjudication is a statement "
                f"about degrees of freedom, not about the shapes")
            n += 1
            continue
        want = max(rows[tgt])
        if cur[tgt] != want:
            fails.append(
                f"G53 SHAPECURRENT {tgt} {cur[tgt]} while "
                f"SHAPESURVIVE {tgt} runs to {want} points; adding a "
                f"rung cut the surviving shapes from 3 to 2 and the "
                f"spread from 8.3508 to 2.7845 decades, so the older "
                f"adjudication may not stand")
            n += 1
    for tgt, pts in sorted(cur.items()):
        if tgt not in rows or pts not in rows[tgt]:
            fails.append(
                f"G53 SHAPECURRENT {tgt} {pts} has no SHAPESURVIVE "
                f"{tgt} {pts} behind it")
            n += 1
    return n


PREDICTOR = re.compile(r"^PREDICTOR (\S+) (\S+) ([\d.eE+-]+) "
                       r"([\d.eE+-]+) ([-+][\d.eE+-]+)\s*$", re.M)
PREDCRIT = re.compile(r"^PREDICTOR CRITERION (\S+) (\S+)\s*$", re.M)
PREDBEST = re.compile(r"^PREDICTOR BEST (\S+) (\S+)\s*$", re.M)


def g54_predictors_are_ranked_on_a_declared_score():
    """두 점수가 반대로 움직이면 이긴 점수를 고를 수 있다.

    기울기를 나르는 초등 예측자를 찾으면서 세 후보가 모였다 -- Modd,
    체 가중 P, 그리고 P 에 log(N-mk) 를 실은 P_log. 마지막은 **기울기
    비에서는 가장 좋고**(0.7018-0.9212 대 P 의 0.6571-0.8676) **부호
    일치도에서는 더 나쁘다**(0.7241-0.7955 대 0.7367-0.8129). 두 점수가
    반대로 움직이므로 기준을 밝히지 않으면 이긴 쪽을 골라 "개선"이라고
    쓸 수 있다.

    그래서 한 표적에 예측자를 둘 이상 실으면
    PREDICTOR <표적> <이름> <일치도> <기울기비> <감쇠기울기> 를 각각
    내고, PREDICTOR CRITERION <표적> <점수> 로 순위 기준을 밝히고,
    PREDICTOR BEST <표적> <이름> 이 그 기준의 최선과 일치해야 한다.
    기준은 agreement(큰 쪽) 또는 leanratio(1 에 가까운 쪽)다.
    """
    if not os.path.isdir(RESULTS):
        return 0
    n = 0
    for f in sorted(os.listdir(RESULTS)):
        if not f.endswith(".txt"):
            continue
        src = read(os.path.join(RESULTS, f))
        rows = {}
        for tgt, nm, ag, rt, sl in PREDICTOR.findall(src):
            try:
                rows.setdefault(tgt, {})[nm] = (float(ag), float(rt))
            except ValueError:
                pass
        if not rows:
            continue
        crit = dict(PREDCRIT.findall(src))
        best = dict(PREDBEST.findall(src))
        for tgt, r in sorted(rows.items()):
            if len(r) < 2:
                continue
            if tgt not in crit or tgt not in best:
                fails.append(
                    f"G54 results/{f} lists {len(r)} predictors for "
                    f"{tgt} and declares no CRITERION or no BEST; the "
                    f"one candidate this repository tested was the "
                    f"best of three on one score and the worst on "
                    f"another")
                n += 1
                continue
            c = crit[tgt]
            if c == "agreement":
                want = max(r, key=lambda nm: r[nm][0])
            elif c == "leanratio":
                want = min(r, key=lambda nm: abs(r[nm][1] - 1.0))
            else:
                fails.append(
                    f"G54 results/{f} declares CRITERION {tgt} {c}, "
                    f"which is not agreement or leanratio")
                n += 1
                continue
            if best[tgt] != want:
                fails.append(
                    f"G54 results/{f} declares BEST {tgt} "
                    f"{best[tgt]} while the criterion {c} makes it "
                    f"{want}")
                n += 1
    return n


LEVEL = re.compile(r"^LEVEL (\S+) (\S+)\s*$", re.M)
UNBOUNDED = re.compile(r"^UNBOUNDED LEVEL (\S+)\s*$", re.M)


def g55_predictors_declare_their_sieve_level():
    """레벨이 `N` 과 함께 자라면 그건 초등 예측자가 아니다.

    기울기를 나르는 예측자를 찾다가 답이 나왔다: 체의 깊이 Q 를
    29 에서 ceil(sqrt(N)) 까지 올리면 부호 일치도가 0.74-0.81 에서
    0.99 로 오르고 기울기 비가 0.9885-0.9942, 감쇠 기울기가 μ 와
    0.05 표준오차 안에 든다. 그런데 Q = sqrt(N) 에서 생존자는 곧
    소수이므로 그건 **유계 모듈러스 대상이 아니다** --
    {#rem:provablehalf} 가 P 를 "모든 조건이 곱셈적이거나 유계
    모듈러스"라고 부른 그 뜻에서 초등이 아니다. 고정 레벨 29 에서는
    일치도가 0.7367-0.8129 이고 `N` 이 커질수록 나빠진다.

    G54 는 점수 기준만 못박으므로, 그 기준의 최선이 레벨이 자라는
    대상이면 "초등 예측자를 찾았다"로 읽힐 수 있다. 그래서
    PREDICTOR 로 실린 이름은 LEVEL <이름> <값> 을 가져야 하고, 값이
    정수가 아니면(레벨이 자라면) BEST 로 지명될 때
    UNBOUNDED LEVEL <이름> 이 함께 있어야 한다.
    """
    if not os.path.isdir(RESULTS):
        return 0
    lev, unb = {}, set()
    for f in sorted(os.listdir(RESULTS)):
        if not f.endswith(".txt"):
            continue
        src = read(os.path.join(RESULTS, f))
        for nm, v in LEVEL.findall(src):
            lev[nm] = v
        unb.update(UNBOUNDED.findall(src))
    n = 0
    for f in sorted(os.listdir(RESULTS)):
        if not f.endswith(".txt"):
            continue
        src = read(os.path.join(RESULTS, f))
        names = set(nm for _t, nm, _a, _r, _s
                    in PREDICTOR.findall(src))
        if not names:
            continue
        for nm in sorted(names):
            if nm not in lev:
                fails.append(
                    f"G55 results/{f} lists the predictor {nm} and no "
                    f"file declares LEVEL {nm}; the one predictor that "
                    f"carried the lean here did so at a level growing "
                    f"like sqrt(N), which is not a bounded modulus")
                n += 1
        for _t, nm in PREDBEST.findall(src):
            if nm not in lev:
                continue
            if not lev[nm].isdigit() and nm not in unb:
                fails.append(
                    f"G55 results/{f} names {nm} best and its LEVEL "
                    f"is {lev[nm]}, which grows with N, without an "
                    f"UNBOUNDED LEVEL {nm}; a predictor whose sieve "
                    f"level grows is not elementary in the sense the "
                    f"programme needs")
                n += 1
    return n


THRESHFROM = re.compile(r"^THRESHOLD FROM (\S+) (\S+)\s*$", re.M)
TSTATA = re.compile(r"^TSTAT slope_(\S+?)_a([\d.]+) ([\d.eE+-]+)\s*$",
                    re.M)


def g56_thresholds_follow_from_the_declared_statistic():
    """문턱은 눈으로 고르는 게 아니라 규칙으로 나와야 한다.

    체 레벨을 Q = N^alpha 로 쓸어 부호 기울기가 언제 미끄러지기를
    멈추는지 물었더니, **일치도**는 어느 레벨에서도 두 표준오차에
    못 미치고(고정 레벨조차 1.10) **기울기 비**는 3.52, 7.09, 4.31,
    1.52, 0.22, 2.39 로 alpha = 0.3 에서 해소되지 않게 된다. 두
    통계량이 다른 답을 주므로 문턱은 어느 것에서 나왔는지 밝히지
    않으면 고를 수 있다. 지난 사이클 {#rem:sievedepth} 가 "둘 다
    나빠진다"고 적은 것도 그 둘을 가르지 않았기 때문이다.

    그래서 LEVEL <이름>_threshold <값> 을 내면
    THRESHOLD FROM <이름> <통계량> 으로 출처를 밝혀야 하고, 그 값은
    그 통계량의 TSTAT slope_<통계량>_a<alpha> 중 **처음으로 2 미만이
    되는 alpha** 와 같아야 한다.
    """
    if not os.path.isdir(RESULTS):
        return 0
    n = 0
    for f in sorted(os.listdir(RESULTS)):
        if not f.endswith(".txt"):
            continue
        src = read(os.path.join(RESULTS, f))
        lv = dict(LEVEL.findall(src))
        thr = dict(THRESHFROM.findall(src))
        fams = {}
        for fam, a, t in TSTATA.findall(src):
            try:
                fams.setdefault(fam, {})[float(a)] = abs(float(t))
            except ValueError:
                pass
        for nm, val in sorted(lv.items()):
            if not nm.endswith("_threshold"):
                continue
            if nm not in thr:
                fails.append(
                    f"G56 results/{f} declares LEVEL {nm} {val} and "
                    f"no THRESHOLD FROM {nm}; the two statistics on "
                    f"this sweep disagree, one resolving at 3.52 and "
                    f"the other at 1.10, so the source decides the "
                    f"answer")
                n += 1
                continue
            fam = thr[nm]
            if fam not in fams:
                fails.append(
                    f"G56 results/{f} says {nm} comes from {fam} and "
                    f"no TSTAT slope_{fam}_a* is declared")
                n += 1
                continue
            unres = sorted(a for a, t in fams[fam].items() if t < 2.0)
            want = "%.1f" % unres[0] if unres else "none"
            if val != want:
                fails.append(
                    f"G56 results/{f} declares LEVEL {nm} {val} while "
                    f"the first unresolved alpha of {fam} is {want}; "
                    f"a threshold is the rule applied, not the level "
                    f"chosen")
                n += 1
    return n


AXIS = re.compile(r"^AXIS (\S+) (\S+) (\S+)\s*$", re.M)
THRDIFFER = re.compile(r"^THRESHOLDS DIFFER (\S+)\s*$", re.M)


def g57_one_axis_carries_every_statistic_s_threshold():
    """한 축 위의 문턱은 통계량마다 다르다 -- 싼 쪽을 옮겨 쓰면 안 된다.

    체 레벨 alpha 축 위에서 두 통계량이 서로 다른 곳에서 꺾인다.
    부호 기울기 비는 alpha = 0.3 에서 미끄러짐이 멈추고
    ({#rem:levelthreshold}), 수요의 잔여 몫은 거기서 0.4866-0.6022 로
    레벨 29 의 0.5307-0.6310 에서 거의 안 내려가며 절반이 되는 건
    alpha = 0.5 에서다. 감당 가능한 레벨이 사는 것은 무한 레벨이 살
    것의 **8-15%** 뿐이다.

    즉 "부호는 alpha = 0.3 에서 붙잡힌다"를 수요로 옮겨 읽으면 경로가
    싸 보인다. 그래서 문턱을 선언하는 축은 그 위에서 잰 **모든**
    통계량에 대해 AXIS <축> <통계량> <문턱> 을 내야 하고, 두 문턱이
    다르면 THRESHOLDS DIFFER <축> 을 붙여야 한다. 같은데 붙여도
    실패한다.
    """
    if not os.path.isdir(RESULTS):
        return 0
    n = 0
    for f in sorted(os.listdir(RESULTS)):
        if not f.endswith(".txt"):
            continue
        src = read(os.path.join(RESULTS, f))
        rows = {}
        for ax, st, th in AXIS.findall(src):
            rows.setdefault(ax, {})[st] = th
        if not rows:
            continue
        diff = set(THRDIFFER.findall(src))
        for ax, r in sorted(rows.items()):
            if len(r) < 2:
                fails.append(
                    f"G57 results/{f} declares one statistic on the "
                    f"axis {ax}; a threshold on an axis is a fact "
                    f"about the statistic, and the one axis this "
                    f"repository swept carries two thresholds that "
                    f"differ")
                n += 1
                continue
            same = len(set(r.values())) == 1
            if same == (ax in diff):
                fails.append(
                    f"G57 results/{f} declares thresholds "
                    f"{sorted(r.items())} on {ax} and the "
                    f"THRESHOLDS DIFFER {ax} marker "
                    f"{'is present anyway' if same else 'is missing'}")
                n += 1
    return n


ACCOUNT = re.compile(r"^ACCOUNT (\S+) ([\d.eE+-]+) ([\d.eE+-]+) "
                     r"([\d.eE+-]+) ([\d.eE+-]+)\s*$", re.M)
ACCUNEXPL = re.compile(r"^ACCOUNT UNEXPLAINED (\S+)\s*$", re.M)


def g58_negligible_remainders_are_accounted_for():
    """0 에 가깝다고 적은 나머지는 무엇으로 이루어졌는지 대야 한다.

    {#rem:leveldemand} 는 alpha = 1/2 에서 잔여 몫 0.2271-0.2525 를
    남기고 그것이 바닥인지 물었다. 답은 바닥이 아니라 **세는 규약의
    값**이었다: 그 레벨에서 생존자가 곧 소수이므로 log(N-mk) 가 log p
    와 같고, 로그 가중 예측자를 쓰면 몫이 0.006403-0.008981 로 떨어진다.
    남는 그 조각도 정체가 있다 -- 체가 소수를 그 배수와 함께 지우므로
    sqrt(N) 이하 소수는 진짜 기여자인데 제거된다. 그 질량을 직접 재면
    0.007108-0.009755 로 잔여와 겹친다.

    0 에 가까운 나머지는 "무시할 만하다"로 넘어가기 쉽고, 그러면 그것이
    규약의 값인지 진짜 바닥인지가 흐려진다. 그래서 상한이 0.02 아래인
    PERN 범위는 ACCOUNT <이름> <잔여 최소> <잔여 최대> <설명 최소>
    <설명 최대> 로 직접 잰 설명을 함께 내야 하고, 두 구간이 겹치지
    않으면 ACCOUNT UNEXPLAINED <이름> 을 적어야 한다.
    """
    if not os.path.isdir(RESULTS):
        return 0
    acc, unexp = {}, set()
    for f in sorted(os.listdir(RESULTS)):
        if not f.endswith(".txt"):
            continue
        src = read(os.path.join(RESULTS, f))
        for lab, a, b, c, d in ACCOUNT.findall(src):
            try:
                acc[lab] = (float(a), float(b), float(c), float(d))
            except ValueError:
                pass
        unexp.update(ACCUNEXPL.findall(src))
    n = 0
    for f in sorted(os.listdir(RESULTS)):
        if not f.endswith(".txt"):
            continue
        src = read(os.path.join(RESULTS, f))
        for lab, _cnt, lo, hi in PERN.findall(src):
            try:
                hiv = float(hi)
            except ValueError:
                continue
            if hiv >= 0.02:
                continue
            if lab not in acc:
                fails.append(
                    f"G58 results/{f} declares PERN {lab} with an "
                    f"upper end of {hiv:g} and no ACCOUNT {lab}; the "
                    f"one near-zero remainder this repository checked "
                    f"was the cost of a counting convention and not a "
                    f"floor")
                n += 1
                continue
            rlo, rhi, clo, chi = acc[lab]
            if not (rlo <= chi and clo <= rhi) and lab not in unexp:
                fails.append(
                    f"G58 ACCOUNT {lab} gives a remainder "
                    f"[{rlo:g}, {rhi:g}] and an account "
                    f"[{clo:g}, {chi:g}] that do not overlap, with no "
                    f"ACCOUNT UNEXPLAINED {lab}")
                n += 1
    return n


SPREADCAP = re.compile(r"a spread of [\d.]+\s+\(cap [\d.]+\)")
CONSTSPREAD = re.compile(r"^CONSTSPREAD (\S+) ([\d.eE+-]+) "
                         r"([\d.eE+-]+)\s*$", re.M)
CONSTDRIFT = re.compile(r"^CONST DRIFTS (\S+)\s*$", re.M)


def g59_capped_spreads_declare_their_sampling_error():
    """추정량 자신의 표집오차보다 좁은 상한을 걸면 반드시 걸린다.

    세 비의 항등식을 확인하면서 그것들을 잇는 상수 c 가 N 에 걸쳐
    상수인지 0.02 상한으로 물었다. 관측된 로그 폭은 0.219062 로 상한을
    열 배 넘었지만, c 는 256 뽑기의 **중앙값**이라 자기 표집오차를
    갖는다. 뽑기를 16 개씩 열여섯 조로 갈라 조별 중앙값의 산포를 재면
    전체 중앙값의 표준오차가 나오고, 그것이 주는 폭이 0.237831 --
    관측된 폭과 같은 크기다. 즉 c 는 흐르지 않고, 상한이 추정량보다
    정밀했을 뿐이다.

    그래서 `a spread of X (cap Y)` 로 폭을 판정하는 결과 파일은
    CONSTSPREAD <이름> <관측 폭> <표집 폭> 을 내야 하고, 관측이
    표집의 두 배를 넘으면 CONST DRIFTS <이름> 을, 안 넘으면 그것을
    붙이지 않아야 한다.
    """
    if not os.path.isdir(RESULTS):
        return 0
    n = 0
    for f in sorted(os.listdir(RESULTS)):
        if not f.endswith(".txt"):
            continue
        src = read(os.path.join(RESULTS, f))
        if not SPREADCAP.search(src):
            continue
        rows = CONSTSPREAD.findall(src)
        if not rows:
            fails.append(
                f"G59 results/{f} judges a spread against a cap and "
                f"declares no CONSTSPREAD; the one such cap this "
                f"repository set was ten times tighter than the "
                f"estimator behind the quantity")
            n += 1
            continue
        drifts = set(CONSTDRIFT.findall(src))
        for lab, obs, samp in rows:
            try:
                o, s = float(obs), float(samp)
            except ValueError:
                continue
            if (o > 2.0 * s) != (lab in drifts):
                fails.append(
                    f"G59 results/{f} declares CONSTSPREAD {lab} "
                    f"{o:g} against a sampling spread {s:g} and the "
                    f"CONST DRIFTS {lab} marker "
                    f"{'is present anyway' if o <= 2.0 * s else 'is missing'}")
                n += 1
    return n


SHAPEGAP = re.compile(r"^SHAPEGAP (\S+) ([\d.eE+-]+) ([\d.eE+-]+)"
                      r"\s*$", re.M)
SHAPETIED = re.compile(r"^SHAPES TIED (\S+)\s*$", re.M)


def g60_shape_gaps_are_read_against_their_own_error():
    """모양 사이의 r.m.s. 차는 r.m.s. 자신의 표준오차로 재야 한다.

    평탄도 F 에 두 모양을 맞추니 멱법칙이 유계 모양보다 r.m.s. 로
    0.000163 앞선다. 그런데 여덟 점 두 모수의 r.m.s. 는 자기
    표준오차가 0.001870 (28.9%) 이고, 그 차는 그것의 **0.09 배**다 --
    자료는 두 모양을 전혀 가르지 못한다. 가른 것은 적합이 아니라
    코시-슈바르츠(F <= 1)였다. 사다리 쪽도 같다: 열두 가로대에서 최선
    둘의 차가 0.000626, 표준오차가 0.000828 이다.

    "무엇이 더 잘 맞는다"는 문장은 그 차가 표준오차보다 클 때만
    뜻이 있다. 그래서 SHAPESURVIVE 를 내는 표적은
    SHAPEGAP <표적> <차> <표준오차> 도 내야 하고, 차가 표준오차 이하면
    SHAPES TIED <표적> 을, 넘으면 그것을 붙이지 않아야 한다.

    판정은 **현행 점수의 것 하나**만 읽는다. 처음엔 SHAPEGAP 을
    파일 순서로 덮어쓰면서 SHAPES TIED 는 파일들의 합집합으로 모았는데,
    그 둘이 어긋난 채로 실제 사건이 왔다: audit_primorial_dense.py 가
    같은 표적 ladder_theta 를 209 점에서 다시 판정해 차 0.000238 대
    표준오차 0.000196 으로 **갈랐는데**, 열두 가로대에서 묶였던 옛
    SHAPES TIED 가 합집합에 남아 현행 판정과 모순됐다. G53 이 이미
    "점이 늘면 다시 하라"고 강제하는 마당에, 그 다시 한 판정을 옛
    표지가 뒤집을 수 있으면 안 된다. 그래서 표적마다 SHAPESURVIVE 의
    점수가 가장 큰 파일을 현행으로 잡고, 그 파일의 SHAPEGAP 과 그
    파일의 SHAPES TIED 만 대조한다. 옛 파일의 표지는 자기 점수에 대한
    역사적 진술로 남되 현행을 구속하지 않는다.
    """
    if not os.path.isdir(RESULTS):
        return 0
    # 표적 -> (점수, 파일, 자기 파일의 gap, 자기 파일의 tied 여부)
    current = {}
    for f in sorted(os.listdir(RESULTS)):
        if not f.endswith(".txt"):
            continue
        src = read(os.path.join(RESULTS, f))
        gaps_f = {}
        for lab, g, e in SHAPEGAP.findall(src):
            try:
                gaps_f[lab] = (float(g), float(e))
            except ValueError:
                pass
        tied_f = set(SHAPETIED.findall(src))
        for t, pts, _s, _sp in SHAPESURVIVE.findall(src):
            pts = int(pts)
            if t in current and current[t][0] >= pts:
                continue
            current[t] = (pts, f, gaps_f.get(t), t in tied_f)
    n = 0
    for t in sorted(current):
        pts, f, ge, istied = current[t]
        if ge is None:
            fails.append(
                f"G60 results/{f} adjudicates SHAPESURVIVE {t} at "
                f"{pts} points, the largest declared, and gives no "
                f"SHAPEGAP {t}; the two shape comparisons this "
                f"repository first ran were separated by 0.09 and 0.76 "
                f"of the r.m.s.'s own standard error, which is to say "
                f"not at all")
            n += 1
            continue
        g, e = ge
        if (g <= e) != istied:
            fails.append(
                f"G60 results/{f} is the current adjudication of {t} "
                f"({pts} points) with SHAPEGAP {g:g} against a "
                f"standard error {e:g}, and its SHAPES TIED {t} marker "
                f"{'is present anyway' if g > e else 'is missing'}")
            n += 1
    return n


GAINSPLIT = re.compile(r"^GAINSPLIT (\S+) ([-+][\d.eE+-]+) "
                       r"([-+][\d.eE+-]+) ([-+][\d.eE+-]+)\s*$", re.M)
GAINEXP = re.compile(r"^TSTAT slope_leanidentity_G ", re.M)


def g61_whole_range_gains_declare_their_split():
    """전체 구간의 상쇄 이득은 부분들의 반대 거동을 가린다.

    확대에 걸친 이득 G 의 지수가 전체에서 +0.153911 인데, 크기 순으로
    가르면 위 10분의 1 에서 +0.077963, 아래 10분의 9 에서 +0.340006
    이다. 즉 작은 항들은 제곱근보다 **잘** 상쇄되고(θ'/2 = 0.28 을
    2.33 표준오차 넘김) 큰 항들은 거의 상쇄되지 않는다 -- 머리에서
    같은 부호를 갖는 비율이 1.0000 에서 0.8274 로, 여덟 N 내내 8할이
    넘는다. 전체 수치 하나만 보면 이 대비가 통째로 사라진다.

    {#rem:nocrossk} 의 T4 는 최상위 십분위의 **질량 몫**(0.3337-0.3587)
    만 보고 "무거운 꼬리가 아니다"라고 읽었는데, 물었어야 할 것은 그
    십분위가 얼마나 상쇄되는가였다.

    그래서 전체 구간의 이득 지수를 발표하면
    GAINSPLIT <이름> <머리 지수> <꼬리 지수> <전체 지수> 도 내야 하고,
    세 값이 머리 <= 전체 <= 꼬리를 만족해야 한다.
    """
    if not os.path.isdir(RESULTS):
        return 0
    rows, seen = {}, False
    for f in sorted(os.listdir(RESULTS)):
        if not f.endswith(".txt"):
            continue
        src = read(os.path.join(RESULTS, f))
        if GAINEXP.search(src):
            seen = True
        for lab, h, t, w in GAINSPLIT.findall(src):
            try:
                rows[lab] = (float(h), float(t), float(w))
            except ValueError:
                pass
    n = 0
    if seen and not rows:
        fails.append(
            "G61 a whole-range cross-k gain exponent is declared and "
            "no file declares a GAINSPLIT; the one such exponent this "
            "repository split ran from +0.078 on the top tenth to "
            "+0.340 on the rest")
        n += 1
    for lab, (h, t, w) in sorted(rows.items()):
        if not (h <= w <= t):
            fails.append(
                f"G61 GAINSPLIT {lab} gives head {h:g}, whole {w:g}, "
                f"tail {t:g}, which do not bracket; the whole range "
                f"must lie between its parts")
            n += 1
    return n


SPLITOVERLAP = re.compile(r"^SPLITOVERLAP (\S+) ([\d.eE+-]+) "
                          r"([\d.eE+-]+)\s*$", re.M)


def g62_mass_splits_declare_their_overlap_with_the_range():
    """질량 순 분할을 구간 제한으로 읽으면 틀린다.

    {#rem:gainsplit} 의 머리는 |a_k| = (log k)|H(N;k)| 의 상위 10분의
    1 이고, |H| 는 안쪽 합이 길수록 크므로 그 머리는 작은 k 쪽일 것
    같다. 실제로 재 보니 가장 작은 10분의 1 과의 겹침이 0.2174 에서
    0.3263 뿐이다 -- log k 인자가 k 와 함께 자라 |H| 의 감소와 맞서기
    때문이고, 머리는 k-순서의 0.22 에서 0.32 지점에 중앙값을 두고
    구간 전체에 퍼져 있다. (가중을 떼면 |H| 상위 십분위와는 0.7794 에서
    0.8863 으로 겹친다.)

    질량 순 분할을 "작은 k 쪽"으로 읽으면 그 부분집합에 대한 산술적
    설명을 찾게 되는데, 그런 설명은 없다. 그래서 GAINSPLIT 을 내는
    이름은 SPLITOVERLAP <이름> <최소 겹침> <최대 겹침> 도 내야 하고,
    겹침은 0 과 1 사이여야 한다.
    """
    if not os.path.isdir(RESULTS):
        return 0
    ov, labs = {}, set()
    for f in sorted(os.listdir(RESULTS)):
        if not f.endswith(".txt"):
            continue
        src = read(os.path.join(RESULTS, f))
        for lab, lo, hi in SPLITOVERLAP.findall(src):
            try:
                ov[lab] = (float(lo), float(hi))
            except ValueError:
                pass
        labs.update(lab for lab, _h, _t, _w in GAINSPLIT.findall(src))
    n = 0
    for lab in sorted(labs):
        if lab not in ov:
            fails.append(
                f"G62 GAINSPLIT {lab} is declared and no "
                f"SPLITOVERLAP {lab}; the one mass split this "
                f"repository checked overlapped the corresponding "
                f"range decile by only 0.22 to 0.33")
            n += 1
            continue
        lo, hi = ov[lab]
        if not (0.0 <= lo <= hi <= 1.0):
            fails.append(
                f"G62 SPLITOVERLAP {lab} declares [{lo:g}, {hi:g}], "
                f"which is not a range of overlaps")
            n += 1
    return n


RESIDSCALE = re.compile(r"^RESIDSCALE (\S+) (\d+) ([\d.eE+-]+) "
                        r"([\d.eE+-]+) ([\d.eE+-]+)\s*$", re.M)
CROSSFLOOR = re.compile(r"^CROSSES FLOOR (\S+)\s*$", re.M)


def g63_residual_spreads_are_reported_at_every_scale():
    """한 규모의 잔차로 "설명한다"를 판정하면 방향을 못 본다.

    {#rem:residuearithmetic} 는 일곱 산술형의 수준 지수를 예산의 로그에
    회귀해 상관 0.97565 를 얻고 예산이 퍼짐을 설명한다고 읽었다 -- 한
    규모에서. 네 규모에서 다시 하니 잔차 폭이 0.0218, 0.0240, 0.0129,
    0.0103 으로 **바닥 0.0133 을 x2 와 x4 사이에서 가로지르고** 상관이
    0.97565 에서 0.99577 로 오른다. 즉 두 번째 변수는 작은 N 에서만
    보이고 커지면 사라진다 -- 한 규모만 보면 그 방향을 알 수 없다.

    그래서 잔차 폭을 `PERN <이름>_residual` 로 발표하면
    RESIDSCALE <이름> <규모 수> <최소> <최대> <바닥> 도 내야 하고,
    바닥이 그 사이에 있으면 CROSSES FLOOR <이름> 을 붙여야 한다.
    """
    if not os.path.isdir(RESULTS):
        return 0
    rows, cross, labs = {}, set(), set()
    for f in sorted(os.listdir(RESULTS)):
        if not f.endswith(".txt"):
            continue
        src = read(os.path.join(RESULTS, f))
        for lab, n_, lo, hi, fl in RESIDSCALE.findall(src):
            try:
                rows[lab] = (int(n_), float(lo), float(hi), float(fl))
            except ValueError:
                pass
        cross.update(CROSSFLOOR.findall(src))
        for lab, _c, _lo, _hi in PERN.findall(src):
            if lab.endswith("_residual"):
                labs.add(lab[:-len("_residual")])
    n = 0
    for lab in sorted(labs):
        if lab not in rows:
            fails.append(
                f"G63 a residual spread {lab}_residual is published "
                f"and no RESIDSCALE {lab}; the one such spread this "
                f"repository followed up the scale crossed its floor "
                f"between the second and third rung")
            n += 1
            continue
        cnt, lo, hi, fl = rows[lab]
        if cnt < 2:
            fails.append(
                f"G63 RESIDSCALE {lab} covers {cnt} scale, which "
                f"cannot show a crossing either way")
            n += 1
            continue
        if (lo < fl < hi) != (lab in cross):
            fails.append(
                f"G63 RESIDSCALE {lab} spans [{lo:g}, {hi:g}] about "
                f"a floor {fl:g} and the CROSSES FLOOR {lab} marker "
                f"{'is present anyway' if not (lo < fl < hi) else 'is missing'}")
            n += 1
    return n


TRUST = re.compile(r"^TRUST (\S+) ([\d.eE+-]+) ([\d.eE+-]+)\s*$",
                   re.M)
OUTSIDETRUST = re.compile(r"^FORECAST OUTSIDE (\S+)\s*$", re.M)


def g64_shape_forecasts_declare_their_trust_range():
    """모양이 갈라지는 자리를 넘어서면 그건 자료가 아니라 모양의 말이다.

    원시근사 사다리의 열두 가로대에서 두 모양이 살아남고 theta' = 0.56
    을 2.78 자릿수 떨어진 곳에 놓는다. 그런데 그 둘이 사다리 자신의
    r.m.s. 0.00370 보다 많이 갈라지기 시작하는 자리는 log10 N = 8.1253
    -- 꼭대기 가로대(7.7889)에서 겨우 **0.34 자릿수** 위다. 두 예보는
    11.0762 와 13.8607 로 그보다 훨씬 멀리 있다. 평탄도 쪽도 같아서
    8.6994 에서 갈라지는데 경계는 28.6782 에 있다.

    "두 모양이 X 자릿수 어긋난다"만 적으면 어디까지가 자료의 말인지가
    빠진다. 그래서 SHAPESURVIVE 를 내는 표적은
    TRUST <표적> <갈라지는 log10 N> <예보 log10 N> 도 내야 하고,
    예보가 그 자리를 넘으면 FORECAST OUTSIDE <표적> 을 붙여야 한다.
    """
    if not os.path.isdir(RESULTS):
        return 0
    tr, out, targets = {}, set(), set()
    for f in sorted(os.listdir(RESULTS)):
        if not f.endswith(".txt"):
            continue
        src = read(os.path.join(RESULTS, f))
        for lab, t, fc in TRUST.findall(src):
            try:
                tr[lab] = (float(t), float(fc))
            except ValueError:
                pass
        out.update(OUTSIDETRUST.findall(src))
        targets.update(t for t, _p, _s, _sp
                       in SHAPESURVIVE.findall(src))
    n = 0
    for t in sorted(targets):
        if t not in tr:
            fails.append(
                f"G64 SHAPESURVIVE {t} is declared and no TRUST {t}; "
                f"the one ladder this repository measured parts only "
                f"0.34 decades above its top rung, far below either "
                f"forecast")
            n += 1
            continue
        trust, fc = tr[t]
        if (fc > trust) != (t in out):
            fails.append(
                f"G64 TRUST {t} puts the shapes parting at {trust:g} "
                f"and the forecast at {fc:g}, and the "
                f"FORECAST OUTSIDE {t} marker "
                f"{'is present anyway' if fc <= trust else 'is missing'}")
            n += 1
    return n


FROZEN = re.compile(r"^FROZEN (\S+) ([-+][\d.eE+-]+) "
                    r"([-+][\d.eE+-]+)\s*$", re.M)
TRENDCONV = re.compile(r"^TREND CONVENTION (\S+)\s*$", re.M)


def g65_frozen_constants_declare_the_other_convention():
    """상수를 쓸기 최댓값에 얼려 두면 추세의 부호가 바뀔 수 있다.

    {#rem:provablehalf} 의 규칙 W3 은 고전적 한계가 예산의 13.98 에서
    19.83 배를 쓰며 "게다가 나빠진다"고 적는다. 그 계산은 함축 상수 A 를
    1 로(사실상 쓸기 최댓값으로) 얼려 둔 것이다. A 를 각 N 자신의
    최댓값으로 두면 A 가 1.2119 에서 0.3487 로 무너지고 기울기가
    +0.125950 에서 **-0.325836** 으로 뒤집힌다 -- 7.03 표준오차로
    나빠지는 게 아니라 좋아진다.

    어느 규약도 틀리지 않았다. 그 remark 는 상한을 원해서 최댓값을
    쓴다고 스스로 적고 있다. 말할 수 없는 것은 "나빠진다"를 규약 없이
    적는 것이다.

    그래서 상수를 얼려 계산한 추세를 발표하면
    FROZEN <이름> <언 기울기> <점별 기울기> 도 내야 하고, 두 부호가
    다르면 TREND CONVENTION <이름> 을 붙여야 한다. 같은데 붙여도
    실패한다.
    """
    if not os.path.isdir(RESULTS):
        return 0
    rows, conv, labs = {}, set(), set()
    for f in sorted(os.listdir(RESULTS)):
        if not f.endswith(".txt"):
            continue
        src = read(os.path.join(RESULTS, f))
        for lab, a, b in FROZEN.findall(src):
            try:
                rows[lab] = (float(a), float(b))
            except ValueError:
                pass
        conv.update(TRENDCONV.findall(src))
        for lab, _c, _lo, _hi in PERN.findall(src):
            if lab.endswith("_A"):
                labs.add(lab[:-2])
    n = 0
    for lab in sorted(labs):
        if lab not in rows:
            fails.append(
                f"G65 a constant range {lab}_A is published and no "
                f"FROZEN {lab}; the one trend this repository "
                f"recomputed with the constant taken per point "
                f"changed sign")
            n += 1
    for lab, (a, b) in sorted(rows.items()):
        if ((a > 0) != (b > 0)) != (lab in conv):
            fails.append(
                f"G65 FROZEN {lab} gives {a:g} frozen and {b:g} per "
                f"point and the TREND CONVENTION {lab} marker "
                f"{'is present anyway' if (a > 0) == (b > 0) else 'is missing'}")
            n += 1
    return n


FORECASTBOTH = re.compile(r"^FORECAST BOTH (\S+) ([\d.eE+-]+) "
                          r"([\d.eE+-]+) ([\d.eE+-]+)\s*$", re.M)
CONVSPLIT = re.compile(r"^FORECAST CONVENTION SPLIT (\S+)\s*$", re.M)


def g66_forecasts_on_frozen_constants_declare_both():
    """언 상수 위의 예보는 상수를 풀었을 때의 값도 함께 대야 한다.

    {#rem:provablehalf} 의 10^5474.8 은 함축 상수 A 를 쓸기 최댓값에
    얼려 푼 것이다. A 는 상수가 아니라 각 N 의 최댓값이고 10.30
    표준오차로 1.2119 에서 0.3487 로 떨어지며, 그 감쇠의 모양은 두
    후보가 r.m.s. 자기 오차의 0.44 배 안에서 묶여 가려지지 않는다.
    A 를 각 모양으로 외삽해 같은 방정식을 풀면 교차가 10^8.96 과
    10^10.29 로 나온다 -- 얼린 값과 5465 자릿수 차이다.

    두 답은 서로 다른 물음의 답이다. A 를 최댓값에 얼리는 것은 지출의
    상한을 원할 때 옳고, 떨어뜨리는 것은 자료가 하는 일의 기술이다.
    그러나 하나만 적으면 읽는 쪽은 그 차이를 볼 수 없다.

    그래서 TREND CONVENTION 을 단 이름이 예보를 먹이면(TRUST 로
    확인된다) FORECAST BOTH <이름> <언 예보> <풀린 최소> <풀린 최대>
    를 내야 하고, 언 예보가 그 구간 밖이면
    FORECAST CONVENTION SPLIT <이름> 을 붙여야 한다.
    """
    if not os.path.isdir(RESULTS):
        return 0
    both, split, conv, trust = {}, set(), set(), set()
    for f in sorted(os.listdir(RESULTS)):
        if not f.endswith(".txt"):
            continue
        src = read(os.path.join(RESULTS, f))
        for lab, a, lo, hi in FORECASTBOTH.findall(src):
            try:
                both[lab] = (float(a), float(lo), float(hi))
            except ValueError:
                pass
        split.update(CONVSPLIT.findall(src))
        conv.update(SHAPETIED.findall(src))
        trust.update(OUTSIDETRUST.findall(src))
    n = 0
    for lab in sorted(conv & trust):
        if lab not in both:
            fails.append(
                f"G66 {lab} forecasts outside its trust range on an "
                f"input whose shape is tied, and declares no "
                f"FORECAST BOTH {lab}; the "
                f"one such forecast this repository unfroze moved by "
                f"5465 decades")
            n += 1
            continue
        a, lo, hi = both[lab]
        inside = lo <= a <= hi
        if inside == (lab in split):
            fails.append(
                f"G66 FORECAST BOTH {lab} puts the frozen answer at "
                f"{a:g} against [{lo:g}, {hi:g}] and the "
                f"FORECAST CONVENTION SPLIT {lab} marker "
                f"{'is present anyway' if inside else 'is missing'}")
            n += 1
    return n


MAXOVER = re.compile(r"max ratio over every")
CROSSAXIS = re.compile(r"^CROSSAXIS (\S+) (\d+) (\d+) ([\d.eE+-]+)"
                       r"\s*$", re.M)
AXISRISE = re.compile(r"^AXIS RISE (\S+)\s*$", re.M)


def g67_maximised_constants_are_checked_off_axis():
    """한 축에서 최대를 취한 상수는 다른 축에서도 봐야 한다.

    {#rem:provablehalf} 는 고전적 모양의 상수 A 를 안쪽 길이 N/k 로
    잘라 가며 최대를 취한다 -- 1.2119, 1.0710, 0.7309, 0.3363. 그건
    고전적 추정이 이미 통제하는 축이다. 통제하지 않는 축은 k 이고,
    "k 에 대해 균일하게"가 바로 그 축의 진술이다.

    안쪽 길이를 옥타브로 고정하고 그 안에서 k 를 (다섯 N 에 걸쳐 인수
    16 만큼) 움직여 보면 축이 평평하지 않다: 짧은 안쪽합에서는 비가
    떨어지고(기울기 -0.555) 긴 쪽에서는 오른다(+0.167 에 3.37 표준오차,
    +0.228 에 2.18). 여섯 옥타브 중 둘에서 해소된 상승이 있고, 그 둘이
    초등 합의 0.0673 을 나른다 -- 드리프트는 실재하지만 질량이 있는
    곳은 아니다.

    그래서 한 축에서 최대를 취한 상수를 발표하는 결과 파일은
    CROSSAXIS <파일이름> <옥타브 수> <상승 수> <상승분의 질량 몫> 으로
    덮여야 하고, 상승이 하나라도 있으면 AXIS RISE <파일이름> 를
    붙여야 한다.
    """
    if not os.path.isdir(RESULTS):
        return 0
    rows, rise = {}, set()
    for f in sorted(os.listdir(RESULTS)):
        if not f.endswith(".txt"):
            continue
        src = read(os.path.join(RESULTS, f))
        for lab, n_, r_, m_ in CROSSAXIS.findall(src):
            try:
                rows[lab] = (int(n_), int(r_), float(m_))
            except ValueError:
                pass
        rise.update(AXISRISE.findall(src))
    n = 0
    for f in sorted(os.listdir(RESULTS)):
        if not f.endswith(".txt"):
            continue
        src = read(os.path.join(RESULTS, f))
        if not MAXOVER.search(src):
            continue
        lab = f[:-4]
        if lab not in rows:
            fails.append(
                f"G67 results/{f} maximises a constant over one axis "
                f"and no file declares CROSSAXIS {lab}; the one such "
                f"constant this repository checked off axis rises "
                f"resolvably at two octaves of six")
            n += 1
            continue
        cnt, ris, mass = rows[lab]
        if (ris > 0) != (lab in rise):
            fails.append(
                f"G67 CROSSAXIS {lab} reports {ris} rising octaves "
                f"of {cnt} and the AXIS RISE {lab} marker "
                f"{'is present anyway' if ris == 0 else 'is missing'}")
            n += 1
        elif not (0.0 <= mass <= 1.0):
            fails.append(
                f"G67 CROSSAXIS {lab} declares a mass share {mass:g}, "
                f"which is not a share")
            n += 1
    return n


SITSIN = re.compile(r"^SITSIN (\S+) ([\d.eE+-]+) ([\d.eE+-]+)\s*$",
                    re.M)


def g68_subset_claims_declare_the_share_they_carry():
    """"큰 항들에 있다"는 그 부분집합이 얼마를 나르는지 대야 한다.

    {#rem:signmass} 는 "상관은 큰 항들에 있다"고 적는다. 세는 부호는
    0.4121-0.4808 로 균형이고 무게를 주면 0.2273-0.3207 이라 간격이
    0.1325-0.1848 인데, 상위 십분위를 빼도 간격이 0.0637-0.1038 남는다 --
    머리가 나르는 것은 전체의 44 에서 52 퍼센트뿐이다. 같은 크기를
    다시 부호매기면 간격이 0.0004-0.0044 로 사라지므로 남은 절반도
    mu 의 것이다.

    "X 에 있다"는 문장은 X 를 빼고도 남는 양을 함께 적을 때만 뜻이
    있다. 그래서 <이름>_whole 과 <이름>_nohead 의 PERN 을 짝지어 RATIO
    를 내는 파일은 SITSIN <이름> <최소 몫> <최대 몫> 으로 그 부분집합이
    나르는 몫을 대야 하고, 그 몫은 1 에서 비를 뺀 값과 맞아야 한다.
    """
    if not os.path.isdir(RESULTS):
        return 0
    n = 0
    for f in sorted(os.listdir(RESULTS)):
        if not f.endswith(".txt"):
            continue
        src = read(os.path.join(RESULTS, f))
        sits = {}
        for lab, lo, hi in SITSIN.findall(src):
            try:
                sits[lab] = (float(lo), float(hi))
            except ValueError:
                pass
        for a, b, lo, hi in RATIO.findall(src):
            if not (a.endswith("_nohead") and b.endswith("_whole")):
                continue
            lab = a[:-len("_nohead")]
            if lab not in sits:
                fails.append(
                    f"G68 results/{f} compares {a} with {b} and "
                    f"declares no SITSIN {lab}; the one subset this "
                    f"repository checked carried only 44 to 52 per "
                    f"cent of the effect it was said to hold")
                n += 1
                continue
            try:
                rlo, rhi = float(lo), float(hi)
            except ValueError:
                continue
            slo, shi = sits[lab]
            if (abs(slo - (1.0 - rhi)) > 1e-3
                    or abs(shi - (1.0 - rlo)) > 1e-3):
                fails.append(
                    f"G68 SITSIN {lab} declares [{slo:g}, {shi:g}] "
                    f"while the RATIO [{rlo:g}, {rhi:g}] makes it "
                    f"[{1.0 - rhi:g}, {1.0 - rlo:g}]")
                n += 1
    return n


REGCORR = re.compile(r"^CORR (\S+)_regressors ([-\d.eE+]+)\s*$",
                     re.M)
NOTSEP = re.compile(r"^COEFF NOT SEPARABLE (\S+)\s*$", re.M)


def g69_two_regressor_fits_declare_their_collinearity():
    """두 회귀변수가 겹치면 계수를 따로 읽으면 안 된다.

    {#rem:ladderderived} 는 유도된 계수 4.7036 을 자유 적합의 5.7691 과
    비교했다. 자유 적합은 r.m.s. 를 0.01419 에서 0.00395 로 떨어뜨리니
    "자료가 더 큰 계수를 원한다"고 읽고 싶어진다. 그런데 두 회귀변수
    loglog N/log N 과 1/log N 은 이 2.78 자릿수 위에서 상관이 0.99883 --
    계수와 상수가 서로를 먹는다. 그 비교는 값이 아니라 축의 짧음을
    잰 것이다.

    그래서 자유 계수를 유도 계수와 대는 파일은 CORR <이름>_regressors
    로 두 회귀변수의 상관을 대야 하고, 그 절댓값이 0.99 이상이면
    COEFF NOT SEPARABLE <이름> 을 붙여 계수를 따로 읽지 말라고 적어야
    한다. 상관이 그 아래인데 그 표지를 붙이는 것도 막는다 -- 표지는
    잰 값을 따라야지 그 반대가 아니다.
    """
    if not os.path.isdir(RESULTS):
        return 0
    n = 0
    for f in sorted(os.listdir(RESULTS)):
        if not f.endswith(".txt"):
            continue
        src = read(os.path.join(RESULTS, f))
        cor = {}
        for lab, v in REGCORR.findall(src):
            try:
                cor[lab] = float(v)
            except ValueError:
                pass
        sep = set(NOTSEP.findall(src))
        if "free coefficient" in src and not cor:
            fails.append(
                f"G69 results/{f} reads a free coefficient against a "
                f"derived one and declares no CORR <label>_regressors; "
                f"the one such fit this repository made had its two "
                f"regressors correlated at 0.99883")
            n += 1
        for lab, v in cor.items():
            if abs(v) >= 0.99 and lab not in sep:
                fails.append(
                    f"G69 results/{f} has CORR {lab}_regressors "
                    f"{v:g} and no COEFF NOT SEPARABLE {lab}")
                n += 1
        for lab in sep:
            if lab in cor and abs(cor[lab]) < 0.99:
                fails.append(
                    f"G69 results/{f} declares COEFF NOT SEPARABLE "
                    f"{lab} while CORR {lab}_regressors is "
                    f"{cor[lab]:g}, below 0.99")
                n += 1
    return n


# ----------------------------------------------------------------- G70
def g70_every_citation_exists(docs):
    """번호 붙은 진술만이 아니라 **모든** evidence 마커가 실재해야 한다.

    G1은 statements() 를 쓴다 -- Theorem·Proposition·Lemma·Corollary·
    Conjecture 만 본다. 그런데 이 작업의 발견은 거의 전부 Remark 에
    있고, cited() 의 docstring 이 그것을 이유로 G10·G12·G14 를 옮겨
    놓았는데 **존재 검사만 statements() 에 남았다.**

    세어 보면 마커 115개 중 G1이 22개를 보고, G10의 백틱이 7개를 더
    잡고, 나머지 **91개는 아무 검사도 보지 않는다.** G12와 G14는 파일이
    없으면 "G1이 보고한다", "G10이 보고한다"며 건너뛰는데 Remark 마커에
    대해서는 둘 다 틀린 가정이다.

    그래서 code/ 가 통째로 사라져도 논문은 그대로 남고 게이트는 조용히
    통과한다. 그건 가상의 사고가 아니라 이 프로그램의 다른 갈래가 지금
    겪고 있는 사고이고, 거기서는 인용된 96개 중 21개만 남아 있다.
    복원은 재작성이지 복구가 아니다 -- 원본이 없으면 논문이 인쇄한
    숫자를 맞히는 것 말고 검증할 방법이 없고, 숫자를 인쇄하지 않은
    주장은 그것조차 안 된다. 그러니 조용히 통과하는 경로를 막는다.
    """
    n = 0
    for path, src in docs:
        for label, ev, ln in cited(src):
            s, r = evidence_paths(ev)
            if not os.path.exists(s):
                fails.append(f"G70 {rel(path)}:{ln} cites code/{ev} for "
                             f"'{label}', which does not exist")
                n += 1
            elif not os.path.exists(r):
                fails.append(f"G70 {rel(path)}:{ln} cites code/{ev} for "
                             f"'{label}' but results/"
                             f"{os.path.basename(r)} is missing")
                n += 1
    return n


# ----------------------------------------------------------------- G71
RESIDLINE = re.compile(r"^\s*residuals about [^\n:]*:\s*"
                       r"([+-][\d.]+(?:\s*,\s*[+-][\d.]+)+)\s*$", re.M)
SIGNRUN = re.compile(r"^SIGNRUN (\S+) (\d+) (\d+)\s*$", re.M)


def g71_residual_signs_are_declared():
    """잔차를 r.m.s.로만 판정하면 부호가 보이지 않는다.

    {#rem:primorialgap} 이 사다리의 열한 가로대 선 주위로 내부 네 점의
    잔차를 냈다: +0.0001, +0.0025, +0.0049, +0.0052. r.m.s. 는 0.0038
    이고 발표된 산포는 0.0037 -- 머리카락 하나 차이라 "산포 안"으로
    읽고 넘어가기 쉽다. 그런데 **넷이 전부 양수이고 단조 증가**다.
    그건 흩어지는 게 아니라 선을 떠나는 것이고, r.m.s. 는 정의상
    그것을 볼 수 없다. 같은 자리에서 실제로 기울기가 +0.006778 에서
    +0.007139 로 움직였다.

    부호가 다 같을 확률은 교환가능 뽑기에서 2^{1-n} 이다 -- 네 점이면
    0.125, 여섯 점이면 0.03. r.m.s. 하나로는 절대 도달하지 못하는
    증거인데 어떤 검사도 읽고 있지 않았다.

    그래서 적합의 잔차 목록을 인쇄하는 결과 파일은 같은 이름으로
    SIGNRUN <이름> <같은 부호 개수> <전체 개수> 를 내야 한다. 세는
    것은 게이트가 아니라 스크립트가 하고, 게이트는 냈는지와 그 수가
    실제 목록과 맞는지만 본다 -- G19 가 손으로 적는 것을 막으므로
    표지만 붙여 통과시킬 수 없다.
    """
    if not os.path.isdir(RESULTS):
        return 0
    n = 0
    for f in sorted(os.listdir(RESULTS)):
        if not f.endswith(".txt"):
            continue
        src = read(os.path.join(RESULTS, f))
        lists = RESIDLINE.findall(src)
        if not lists:
            continue
        declared = {lab: (int(a), int(b)) for lab, a, b
                    in SIGNRUN.findall(src)}
        if not declared:
            fails.append(
                f"G71 results/{f} prints a residual list and declares "
                f"no SIGNRUN line; four residuals that are all one "
                f"sign carry evidence an r.m.s. cannot, and the one "
                f"such fit here moved the slope by 0.000361")
            n += 1
            continue
        vals = [float(v) for v in re.split(r"\s*,\s*", lists[0])]
        pos = sum(1 for v in vals if v > 0)
        run, tot = max(pos, len(vals) - pos), len(vals)
        for lab, (a, b) in declared.items():
            if (a, b) != (run, tot):
                fails.append(
                    f"G71 results/{f} declares SIGNRUN {lab} {a} {b} "
                    f"but its residual list is {run} of {tot} one way")
                n += 1
    return n


# ----------------------------------------------------------------- G72
BETWEENCLAIM = re.compile(r"cross-block|between the blocks|"
                          r"between blocks", re.I)
CROSSSHAREM = re.compile(r"^CROSSSHARE (\S+) (\d+) ([\d.eE+-]+)\s*$",
                         re.M)
RESDEP = re.compile(r"^RESOLUTION DEPENDENT (\S+)\s*$", re.M)


def g72_within_between_splits_declare_their_resolution():
    """'안'과 '사이'로 가른 분해는 분할의 해상도를 함께 내야 한다.

    audit_gain_profile.py 가 이득의 감쇠 0.153911 을 블록 안의 상쇄
    0.098386 과 블록 사이의 상쇄 0.055525 로 갈랐다. 3.67 표준오차로
    해소된 실측이고, 그 자체로는 옳다. 그런데 **그 수는 분할의 해상도에
    달려 있다.** 블록 B 개로 자르고 블록 사이의 상쇄를 금지한 값
    sum_d |sum_d a| / l1 은 B = 1 에서 1/G 이고 B = #k 에서 정확히 1 이라
    N 에 안 움직이므로, "사이"의 몫은 B = 1 에서 0 이고 B = #k 에서 1 이다.
    사이의 몫은 그 사이를 걸어가는 함수이지 한 수가 아니다.

    실제로 재 보니 B = 2, 5, 10, 20, 50 에서 0.2452, 0.3684, 0.3608,
    0.3216, 0.2064 -- 단조도 아니고 1.78 배로 벌어진다. 십분위 하나만
    발표하면 정칙 분해처럼 읽히고, 그 위에서 "상쇄의 3분의 1이 블록
    사이의 것"이라는 문장이 만들어진다. 게이트의 어느 검사도 이의를
    제기하지 않았다 -- G24 는 쓸린 매개변수를 요구하지만 그건 적합
    지수에만 걸리고, 분해의 분할은 매개변수로 보이지 않기 때문이다.

    그래서 결과 파일이 '사이'의 상쇄를 진술하면 같은 이름으로
    CROSSSHARE <이름> <블록수> <몫> 을 **셋 이상**의 해상도에서 내야
    하고, 그 몫이 1.5 배를 넘게 벌어지면 RESOLUTION DEPENDENT <이름> 을
    붙여야 한다. 넘지 않으면 붙이지 않아야 한다 -- SHAPES TIED 와 같은
    규약이고, 붙이는 쪽과 안 붙이는 쪽이 모두 진술이다.
    """
    if not os.path.isdir(RESULTS):
        return 0
    n = 0
    for f in sorted(os.listdir(RESULTS)):
        if not f.endswith(".txt"):
            continue
        src = read(os.path.join(RESULTS, f))
        if not BETWEENCLAIM.search(src):
            continue
        rows = {}
        for lab, b, v in CROSSSHAREM.findall(src):
            try:
                rows.setdefault(lab, {})[int(b)] = float(v)
            except ValueError:
                pass
        dep = set(RESDEP.findall(src))
        if not rows:
            fails.append(
                f"G72 results/{f} attributes part of a decay to "
                f"cancellation between blocks and declares no "
                f"CROSSSHARE line; that share is 0 at one block and 1 "
                f"at one k per block, so without the resolution it is "
                f"not a measurement")
            n += 1
            continue
        for lab, byb in sorted(rows.items()):
            if len(byb) < 3:
                fails.append(
                    f"G72 results/{f} declares CROSSSHARE {lab} at "
                    f"{len(byb)} resolution(s); three are the fewest "
                    f"that can show the dependence, and the one case "
                    f"measured here was not even monotone in it")
                n += 1
                continue
            vals = [v for v in byb.values() if v > 0]
            if len(vals) < 2:
                continue
            spread = max(vals) / min(vals)
            if (spread > 1.5) != (lab in dep):
                fails.append(
                    f"G72 results/{f} declares CROSSSHARE {lab} over a "
                    f"spread of {spread:.4f} and the "
                    f"RESOLUTION DEPENDENT {lab} marker "
                    f"{'is present anyway' if spread <= 1.5 else 'is missing'}")
                n += 1
    return n


# ----------------------------------------------------------------- G73
SIGNAGREE = re.compile(r"at which\s+sign|sign agreement of|"
                       r"agreement of\s+sign")
MARGINAL = re.compile(r"^MARGINAL (\S+) ([\d.eE+-]+)\s*$", re.M)
DEGENERATE = re.compile(r"^DEGENERATE (\S+)\s*$", re.M)
DEGCAP = 0.9


def g73_agreements_declare_the_predictor_s_marginals():
    """일치율은 예측자에 분산이 있을 때만 측정이다.

    audit_oddmertens_range.py 가 `k < N^theta'` 에서 sign H 와
    sign Modd(floor(N/k)) 의 일치를 0.5201 로 재고 "0.70 아래"라고 적었다.
    그 문장은 참이지만 약한 예측자를 기술하는 것처럼 읽힌다. 실제로는
    audit_tail_mertens.py 가 재보니 그 범위에서 예측자가 k 의 0.9829 …
    0.9970 에서 음이고 꼬리에서는 **전부** 음이다 -- 상수다. 상수 예측자는
    어느 집합에서든 그 집합이 그 값을 갖는 비율만큼 일치하므로, 0.5201 은
    예측이 아니라 sign H 의 주변 음수율을 되읽은 값이다. 실제로 꼬리에서는
    일치가 관측 음수 몫과 자릿수까지 같고 순열 기준선이 그것을 따라잡는다.

    lab_lean_mechanism.py 의 NULL 은 이 함정을 알고 있었다 -- "H도 대개
    음이고 M(x)도 대개 음이라 일치가 우연히 높다"고 적고 순열로 보정한다.
    그런데 그 사실이 **수치로 발표되지 않아서** 다음 스크립트가 같은
    함정에 다시 걸렸고, 게이트의 어느 검사도 울지 않았다. G54 는 두
    예측자를 비교할 때의 점수 규약을 강제하지만 예측자 하나의 분산은
    보지 않는다.

    그래서 부호 일치를 선언하는 결과 파일은 같은 이름으로
    MARGINAL <이름> <몫> 을 내야 한다 -- 그 파일이 일치를 보고하는 모든
    집합 위에서 예측자가 다수 부호를 갖는 몫 가운데 **가장 큰 것**
    (예측자에게 가장 불리한 경우)이다. 처음엔 가장 작은 것으로 썼다가
    audit_oddmertens_range 에서 걸렸다: 그 파일은 창 둘을 보고하는데 통제
    창이 0.7928, 문제의 창이 0.9970 이어서, 최솟값을 내면 축퇴한 창이
    표시를 피해 간다. 하나라도 축퇴면 파일이 표시를 달아야 한다. 그 몫이 0.9 를
    넘으면 DEGENERATE <이름> 을 붙이고, 넘지 않으면 붙이지 않아야 한다.
    RESOLUTION DEPENDENT 와 같은 규약이고, 붙이는 쪽과 안 붙이는 쪽이
    모두 진술이다.
    """
    if not os.path.isdir(RESULTS):
        return 0
    n = 0
    for f in sorted(os.listdir(RESULTS)):
        if not f.endswith(".txt"):
            continue
        src = read(os.path.join(RESULTS, f))
        if not SIGNAGREE.search(src):
            continue
        lab = f[:-4]
        rates = {a: b for a, b in MARGINAL.findall(src)}
        deg = set(DEGENERATE.findall(src))
        if lab not in rates:
            fails.append(
                f"G73 results/{f} declares a statistic of the form "
                f"'the fraction at which sign X equals sign Y' and no "
                f"MARGINAL {lab}; the one such predictor this "
                f"repository followed up took one sign on 0.9829 of "
                f"its range, so its agreement was the other side's "
                f"marginal rate read back")
            n += 1
            continue
        try:
            r = float(rates[lab])
        except ValueError:
            fails.append(f"G73 results/{f} has an unparsable MARGINAL "
                         f"value for {lab}")
            n += 1
            continue
        if not (0.0 <= r <= 1.0 + 1e-9):
            fails.append(
                f"G73 results/{f} declares MARGINAL {lab} {r:g}, which "
                f"is not a share; a majority share lies in [0.5, 1]")
            n += 1
            continue
        if (r >= DEGCAP) != (lab in deg):
            fails.append(
                f"G73 results/{f} declares MARGINAL {lab} {r:.4f} and "
                f"the DEGENERATE {lab} marker "
                f"{'is present anyway' if r < DEGCAP else 'is missing'}; "
                f"above {DEGCAP:g} an agreement carries no information "
                f"about the pairing")
            n += 1
    return n


# ----------------------------------------------------------------- G74
CLASSFORM = re.compile(r"every N = \d+\^[a-z](?:\s+\d+\^[a-z])*")
CLASSES = re.compile(r"^COPRIME (\d+)\s*$", re.M)
CLASSESFOR = re.compile(r"^COPRIME FOR (\S+) (\d+)\s*$", re.M)
SPLITMARK = re.compile(r"^FIELD SPLIT (\S+)\s*$", re.M)


def g74_parametric_families_declare_their_coprimality_classes():
    """매개변수 꼴로 생성한 가족은 k-집합 계급 수를 선언해야 한다.

    audit_flatness_fill.py 가 "every N = 2^a 5^b 는 홀근기가 5라
    허용 k-집합과 문턱이 배가열과 똑같이 고정된다"고 적고 70점을
    한 장으로 적합했다. 그 문장은 가족의 양 끝에서 거짓이다 --
    열거가 a = 0 과 b = 0 에서 시작하므로 N = 2^a (홀근기 1, k 가
    2 와만 서로소라 5의 배수 k 가 들어온다)와 N = 5^b (홀수 N, k 가
    5 와만 서로소라 짝수 k 가 들어온다)가 섞인다. 계급 셋이 한
    적합에 들어갔고, 그 열 점의 이득이 나머지 예순의 창 밖에 있어서
    (audit_fill_field.py 의 V2) 천장 초과가 4.98 표준오차에서 0.02 로
    씻겨 보였다.

    G34 는 울지 않았다. 그 검사가 세는 것은 **홀**근기이고 2^a 은
    빈 홀근기를 내므로 파일이 RADICALS 2 를 선언해도 규약에 맞았기
    때문이다. k-집합을 고정하는 것은 홀근기가 아니라 N 의 소인수
    전체다 -- 그리고 그 차이는 정확히 축퇴한 구성원에서만 생기므로
    가족을 손으로 나열하는 동안에는 보이지 않고 꼴로 생성하는 순간
    들어온다.

    그래서 결과 파일이 `every N = <꼴>` 로 가족을 선언하면 그 가족이
    나르는 서로소 계급의 수 COPRIME <m> 을 내야 한다. 자기가 못 내면
    다른 파일이 COPRIME FOR <이름> <m> 으로 대신 낼 수 있다 -- 결함이
    발견된 파일을 다시 돌리지 않고 감사가 판정하는 길이다. 그리고
    m 이 1 을 넘으면 FIELD SPLIT <이름> 이 어딘가 있어야 한다: 계급을
    나눠 다시 적합한 파일이 실재한다는 선언이고, 없으면 그 파일의
    적합은 섞인 장 위의 값이다.
    """
    if not os.path.isdir(RESULTS):
        return 0
    forlab, splits = {}, set()
    for f in sorted(os.listdir(RESULTS)):
        if not f.endswith(".txt"):
            continue
        src = read(os.path.join(RESULTS, f))
        for lab, m in CLASSESFOR.findall(src):
            try:
                forlab[lab] = int(m)
            except ValueError:
                pass
        splits.update(SPLITMARK.findall(src))
    n = 0
    for f in sorted(os.listdir(RESULTS)):
        if not f.endswith(".txt"):
            continue
        src = read(os.path.join(RESULTS, f))
        if not CLASSFORM.search(src):
            continue
        lab = f[:-4]
        own = CLASSES.findall(src)
        if own:
            m = int(own[0])
        elif lab in forlab:
            m = forlab[lab]
        else:
            fails.append(
                f"G74 results/{f} declares a family by a parametric "
                f"form and no COPRIME count, its own or another "
                f"file's COPRIME FOR {lab}; the one family this "
                f"repository generated that way carried three "
                f"coprimality classes while its RADICALS line said "
                f"two, because 2^a contributes the empty odd radical")
            n += 1
            continue
        if m > 1 and lab not in splits:
            fails.append(
                f"G74 results/{f} declares COPRIME {m} and no file "
                f"declares FIELD SPLIT {lab}; above one class the "
                f"fits over that family are fits over a mixture, and "
                f"the mixture washed a 4.98-standard-error excess "
                f"down to 0.02")
            n += 1
    return n


# ----------------------------------------------------------------- G75
TOLJUDGE = re.compile(r"\btol ([\d.]+)")
READSRESULT = re.compile(r"results/[a-z_0-9]+\.txt")
PRINTBOUND = re.compile(r"^PRINTBOUND (\S+) (\d+) ([\d.eE+-]+)\s*$", re.M)
PRINTBOUNDFOR = re.compile(r"^PRINTBOUND FOR (\S+) (\d+) "
                           r"([\d.eE+-]+)\s*$", re.M)
TOLBELOW = re.compile(r"^TOL BELOW PRINT (\S+)\s*$", re.M)
TOLOWN = re.compile(r"^TOL NOT FROM PRINT (\S+)\s*$", re.M)


def g75_tolerances_follow_the_printing_they_judge():
    """인쇄된 표를 판정하는 허용오차는 그 인쇄에서 나와야 한다.

    두 번 물렸다. audit_slope_significance.py 의 M1 이 10^-5 로 통제를
    걸었다가 반증됐고, 재 보니 모든 간극이 인쇄 반올림이 강제하는 자기
    한계 안이었다 -- 그 remark 가 "M1 이 반증하는 것은 허용오차"라고
    적었다. 그리고 audit_shape_power.py 의 P1 이 0.000001 로 걸었다가
    같은 이유로 반증됐다. 표가 소수 넷째 자리까지 찍으면 각 값은 자기를
    만든 값의 0.00005 안이고, r.m.s. 는 그 섭동에 대해 비확대이므로
    그보다 좁은 허용오차는 자료와 무관하게 발화한다.

    두 번 다 같은 실수이고 게이트의 어느 검사도 울지 않았다. G17 은 N 에
    의존하는 문턱을 계산하라고 하지만 인쇄 정밀도에서 나오는 문턱은
    그물 밖이었다.

    그래서 다른 결과 파일을 읽으면서 tol <값> 으로 통제를 판정하는 결과
    파일은 PRINTBOUND <이름> <소수자리> <한계> 를 내야 한다 (다른 파일이
    PRINTBOUND FOR <이름> 으로 대신 낼 수 있다). 그리고 판정에 쓴 가장
    작은 허용오차가 그 한계보다 작으면 TOL BELOW PRINT <이름> 을 붙여야
    한다 -- 붙이는 쪽도 진술이다. 그런 파일은 자기 통제가 인쇄 때문에
    발화한다는 것을 알고 그렇게 적은 것이고, 모르고 지나간 것과는 다르다.
    """
    if not os.path.isdir(RESULTS):
        return 0
    forlab, below, own_only = {}, set(), set()
    for f in sorted(os.listdir(RESULTS)):
        if not f.endswith(".txt"):
            continue
        src = read(os.path.join(RESULTS, f))
        for lab, dec, bd in PRINTBOUNDFOR.findall(src):
            try:
                forlab[lab] = (int(dec), float(bd))
            except ValueError:
                pass
        below.update(TOLBELOW.findall(src))
        own_only.update(TOLOWN.findall(src))
    n = 0
    for f in sorted(os.listdir(RESULTS)):
        if not f.endswith(".txt"):
            continue
        src = read(os.path.join(RESULTS, f))
        tols = TOLJUDGE.findall(src)
        if not tols or not READSRESULT.search(src):
            continue
        lab = f[:-4]
        if lab in own_only:
            # the tolerance judges quantities this file computed, not
            # values it read off a printed table; that is a statement
            # and the file has made it.
            continue
        own = [(int(d), float(b)) for lb, d, b in PRINTBOUND.findall(src)
               if lb == lab]
        if own:
            dec, bound = own[0]
        elif lab in forlab:
            dec, bound = forlab[lab]
        else:
            fails.append(
                f"G75 results/{f} judges a control with a tol and "
                f"reads another result file, and declares neither "
                f"PRINTBOUND {lab} nor TOL NOT FROM PRINT {lab}; "
                f"the two controls this repository set without asking "
                f"what the printed table carries both refuted the "
                f"tolerance rather than the fit")
            n += 1
            continue
        if bound <= 0 or dec <= 0:
            fails.append(
                f"G75 results/{f} declares PRINTBOUND {lab} {dec} "
                f"{bound:g}, which is not a printing bound")
            n += 1
            continue
        try:
            tightest = min(float(t) for t in tols)
        except ValueError:
            continue
        if (tightest < bound) != (lab in below):
            fails.append(
                f"G75 results/{f} judges at tol {tightest:g} against "
                f"a printing bound of {bound:g} and the "
                f"TOL BELOW PRINT {lab} marker "
                f"{'is present anyway' if tightest >= bound else 'is missing'}"
                f"; below the bound a control fires on rounding, and "
                f"that has to be said rather than discovered")
            n += 1
    return n


# ----------------------------------------------------------------- G76
READLINE = re.compile(r"^READ (\S+\.txt) (\S+ \S+) (\S.*?)\s*$", re.M)


def g76_reads_match_their_source():
    """남의 결과에서 읽었다고 선언한 값이 정말 그 줄의 값이어야 한다.

    audit_target_band.py 가 audit_local_floor.txt 에서 창 B 의 기울기를
    re.S 패턴으로 뽑았는데, FIELD 머리글이 "window B" 를 먼저 부르는
    바람에 검색이 거기서 시작해 창 A 의 기울기를 집었다. 등록된 문턱은
    0.0578 인데 보고된 값은 0.4611 이 됐고, **게이트 일흔다섯 검사가
    전부 통과했다** -- G11 은 인쇄된 숫자가 어느 결과 파일엔가 있는지만
    묻지, 그 값이 자기가 읽었다고 말하는 줄에서 나왔는지는 묻지 않는다.

    그래서 읽은 값을 READ <파일> <표지> <값> 으로 선언하게 하고, 그
    파일에 정말 "<표지> <값>" 줄이 있는지 대조한다. 산문에서 긁지 말고
    표지를 통째로 맞추라는 뜻이기도 하다.
    """
    if not os.path.isdir(RESULTS):
        return 0
    n = 0
    for f in sorted(os.listdir(RESULTS)):
        if not f.endswith(".txt"):
            continue
        src = read(os.path.join(RESULTS, f))
        for name, label, val in READLINE.findall(src):
            dep = os.path.join(RESULTS, name)
            if not os.path.exists(dep):
                fails.append(f"G76 results/{f} declares READ from "
                             f"{name}, which does not exist")
                n += 1
                continue
            want = re.compile(r"^%s %s\s*$"
                              % (re.escape(label), re.escape(val)),
                              re.M)
            if not want.search(read(dep)):
                fails.append(
                    f"G76 results/{f} says it read '{label} {val}' "
                    f"from {name}, which prints no such line")
                n += 1
    return n


# ----------------------------------------------------------------- G77
M9TRIG = re.compile(r"REFUTED[^.\n]*(standard error|\|t\|\s*[<>=]"
                    r"|two standard|below two)", re.I)
M9NAMED = re.compile(r"unresolved|too noisy|cannot tell|not resolved"
                     r"|fails to resolve|noisy to tell|is not zero",
                     re.I)
# 검사 이전에 쓰인 것들. 줄어들기만 한다 -- 목록에 있는데 더 이상
# 방아쇠에 안 걸리거나 이미 경우를 명시했으면 그것도 실패다.
M9GRANDFATHER = {
    "audit_budget_gap.py", "audit_deficit_direct.py",
    "audit_lean_extended.py", "audit_lean_floor.py",
    "audit_lean_identity.py", "audit_level_magnitude.py",
    "audit_level_slope_reach.py", "audit_level_threshold.py",
    "audit_logweight_predictor.py", "audit_oddmertens_range.py",
    "audit_predictable_null.py", "audit_primorial_rung16.py",
    "audit_primorial_rung18.py", "audit_provable_forecast.py",
    "audit_provable_share.py", "audit_residue_coin_rank.py",
    "audit_sieve_depth.py", "audit_signed_gain.py",
    "audit_split_value.py", "audit_survivor_range.py",
}


def _refutation_block(src):
    m = re.search(r"REFUTATION RULE(.*?)"
                  r"(?:\n\s*(?:NO NULL|THE NULL|A NULL)|\n\"\"\")",
                  src, re.S)
    return m.group(1) if m else None


def g77_resolution_rules_name_the_unresolved():
    """해소로 판정하는 반증 규칙은 '미해소'를 경우로 지목해야 한다.

    M9 는 README 규칙인데 같은 오류를 세 번 못 막았다 --
    {#rem:thetalaw} 의 U4, {#rem:alphalocal} 의 Z4, 그리고
    {#rem:deficitshape} 의 B5 이고, 마지막은 두 번째를 기록한 바로 그
    세션이 냈다. 안 무는 규칙은 규칙이 아니라 희망이다.

    방아쇠는 좁다: 표준오차나 t 로 판정하는 REFUTED 절. 그런 절은 언제나
    "너무 시끄러워 못 본다"로 깨질 수 있고, 블록이 그걸 한 번도 말하지
    않으면 조건이 갖지 않은 사건을 지목하고 있는 것이다.

    스무 스크립트가 검사보다 먼저 쓰였다. G16 이 체 해시를, G17 이 타이핑된
    문턱을 명단으로 두는 것과 같이 명단에 둔다 -- 재계산 비용 없이 새
    작업에만 물리고, 세 번의 실패가 전부 새 작업에서 났으므로 거기 걸린다.
    """
    n = 0
    seen = set()
    for f in sorted(os.listdir(CODE)):
        if not f.endswith(".py"):
            continue
        blk = _refutation_block(read(os.path.join(CODE, f)))
        if blk is None:
            continue
        trig = bool(M9TRIG.search(blk))
        named = bool(M9NAMED.search(blk))
        if f in M9GRANDFATHER:
            seen.add(f)
            if not trig or named:
                fails.append(
                    f"G77 code/{f} is on the grandfather list and no "
                    f"longer needs to be; the list may only shrink")
                n += 1
            continue
        if trig and not named:
            fails.append(
                f"G77 code/{f} judges a refutation by a standard "
                f"error or a t and never names the unresolved case; "
                f"'too noisy to tell' is always one way that "
                f"condition fails (M9)")
            n += 1
    for f in sorted(M9GRANDFATHER - seen):
        fails.append(f"G77 audit list names code/{f}, which is not in "
                     f"code/; the list is stale")
        n += 1
    return n



# ----------------------------------------------------------------- G78
PINCONST = re.compile(r"^([A-Z][A-Z0-9_]*)\s*=\s*(.+?)\s*(?:#.*)?$", re.M)
PINNUM = re.compile(r"-?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?")
PINSKIP = ("BLOCK", "CHUNK", "WIDTH", "DEC", "DIGITS", "PREC")
PINPATH = ("OUT", "ROOT", "RES", "RESULTS", "CODE", "VERIFY")
PINGRANDFATHER = {
    ("code/audit_beta_optimal.py", "SWEEP"),
    ("code/audit_cR_window.py", "OCT"),
    ("code/audit_flatness_shape.py", "UMAX"),
    ("code/audit_floor_law.py", "HALF"),
    ("code/audit_ladder_shape.py", "UMAX"),
    ("code/audit_ladder_shape12.py", "UMAX"),
    ("code/audit_mask_deepform.py", "WIN"),
    ("code/audit_mask_formreach.py", "WIN"),
    ("code/audit_mask_rivals.py", "WIN"),
    ("code/audit_provable_uniformity.py", "OCT"),
    ("code/audit_residue_arithmetic.py", "SEED"),
    ("code/audit_residue_coin_rank.py", "OCT"),
    ("code/audit_residue_kexponent.py", "OCT"),
    ("code/audit_shape_trust.py", "UMAX"),
    ("code/audit_weightgap_null.py", "JS"),
    ("code/audit_weightgap_pairing.py", "JS"),
    ("code/lab_elementary_provable.py", "OCT"),
    ("code/lab_primorial_ladder.py", "SEED"),
    ("code/lab_residue_cancellation.py", "OCT"),
}


def _pin_consts(src):
    body = src.split('"""', 2)[-1] if src.count('"""') >= 2 else src
    out = {}
    for m in PINCONST.finditer(body):
        name, rhs = m.group(1), m.group(2)
        if name in PINPATH or any(w in name for w in PINSKIP):
            continue
        if rhs.startswith(("os.", 'r"', "r'", '"', "'", "re.")):
            continue
        if "(" in rhs and not rhs.startswith(("(", "[")):
            continue
        toks = PINNUM.findall(rhs)
        if not toks or len(toks) > 24:
            continue
        if re.search(r"[A-Za-z_]", re.sub(r"[eE]", "", rhs)):
            continue
        out[name] = [float(t) for t in toks]
    return out


def _pin_present(v, nums):
    s = ("%r" % v).rstrip("0").rstrip(".")
    d = len(s.split(".")[1]) if "." in s else 0
    tol = 0.5 * 10.0 ** (-d)
    return any(abs(v - x) <= tol for x in nums)


def _pin_result(p):
    if VERIFY in p:
        return os.path.join(os.path.dirname(os.path.dirname(p)),
                            "results",
                            os.path.basename(p)[:-3] + ".txt")
    return os.path.join(RESULTS, os.path.basename(p)[:-3] + ".txt")


def g78_constants_reach_their_result():
    """결과 파일은 홀로 서야 한다 -- G4 가 이미 그렇게 정했다.

    G4 는 STATISTIC: 과 FIELD: 줄의 **존재**만 본다. 마당을 이름 없이
    부르는 FIELD: 는 아무것도 고정하지 않는다. verify/pass2 와 pass3 가
    다른 트리의 스탬프 둘에서 같은 결함을 찾았다 -- 수는 옳은데 인쇄된
    것으로 재구성되지 않는다. 하나는 격자가, 하나는 셀 색인과 소수
    범위가 없었다.

    같은 기준을 여기 안 대는 건 앞뒤가 안 맞아서 audit_field_pinned.py
    로 쟀고, 154 짝 중 22 에 구멍이 있었다. 그중 일곱은 SEED 다 --
    np.random.default_rng(SEED) 를 모는 상수가 자기 결과 파일에 없다.
    씨앗이 안 적힌 널은 그 파일만으로 재현되지 않는다.

    규칙: 스크립트가 모듈 수준에 고정한 수치 상수는 자기 결과 파일
    어딘가에 나타나야 한다. 값이 우연히 맞을 수 있으므로 이 검사는
    한쪽으로만 유효하다 -- 걸린 것은 구멍이고, 통과한 것은 증명이 아니다.

    BLOCK·CHUNK(타일링)와 WIDTH·DEC·DIGITS·PREC(인쇄)는 수를 못 바꾸므로
    이름으로 뺀다. SEED 는 일부러 안 뺀다.

    검사 이전의 것들은 G16·G17·G77 처럼 명단에 둔다. **명단은 줄기만
    한다** -- 새로 넣으려면 그 스크립트를 고쳐 다시 돌리는 쪽이 싸다.
    """
    n = 0
    for p in py_files(CODE) + py_files(VERIFY):
        rel = os.path.relpath(p, ROOT).replace(os.sep, "/")
        r = _pin_result(p)
        if not os.path.exists(r):
            continue
        cons = _pin_consts(read(p))
        if not cons:
            continue
        nums = [float(t) for t in PINNUM.findall(read(r))]
        for name, vals in sorted(cons.items()):
            if (rel, name) in PINGRANDFATHER:
                continue
            gone = [v for v in vals if not _pin_present(v, nums)]
            if gone:
                fails.append(
                    f"G78 {rel} fixes {name} = "
                    f"{', '.join('%g' % v for v in gone)} and no such "
                    f"number is in {os.path.basename(r)}; a result "
                    f"file has to carry the field it was run on")
                n += 1
    return n


# ----------------------------------------------------------------- G79
DEPLOY = os.path.join(os.path.dirname(ROOT), "deploy", "papers")

# 확정된 결함마다 (이름, 다시 나타나면 안 되는 문자열). 두 트리 모두에서
# 사라졌는지 본다. 문자열은 결함의 서명이지 문장 전체가 아니다.
PARITY = [
    ("pass4 F3  the two lists must not share a numerator silently",
     "against $2.1591"),
    # 서명은 고쳐진 형태에 걸리면 안 된다. 처음 쓴 "(\log k)\,\n  \bigl|"
    # 은 mu^2 를 넣은 줄에도 그대로 남아 있어 늘 실패했다 -- 결함의
    # 서명은 "(k,N)=1}} 바로 뒤에 (\log k)" 라는 자리다.
    ("pass4 F4  B(N) carries mu^2(k)",
     "(k,N)=1}} (\\log k)"),
    ("pass4 F11 the coin factor is not a standard deviation",
     "coin standard deviations"),
    ("pass4 F15 the ineligible rounding of the Heath-Brown share",
     "0.833180,"),
    ("pass4 F19c the residual is measured against N",
     "of $\\SS$ only near"),
    ("pass6 F1  one table, one resolution",
     "e & -0.000879"),
    ("pass6 F2  the adjective must match the cell that carries the effect",
     "conservative by about two orders of magnitude"),
    ("pass6 F3  the sampling error is not the same at every depth",
     "the same $0.0013$ throughout"),
    ("pass6 F7  a generic even N does not have that local factor",
     "at a generic even $N$"),
]


def g79_one_finding_both_trees():
    """한 발견은 두 트리 모두에 가야 한다.

    pass6 이 찾은 것: 확정된 결함 열 중 셋은 어느 트리에도 반영되지
    않았고 일곱은 한 트리에만 반영됐다. 커밋이 `paper/` 와 `deploy/`
    를 번갈아 만지는 동안 같은 수정이 한쪽에만 갔고, 어느 검사도 그
    비대칭을 보지 않았다. 배포본이 이미 고친 결함을 원본이 나르면
    다음 원고가 그것을 물려받는다 -- 그중 하나는 이항 골드바흐 등가
    하한을 고전 문헌에 귀속시키고 있었다.

    규칙: 아래 서명 문자열은 `paper/` 에도 `deploy/papers/` 에도
    나타나면 안 된다. 어느 한쪽에만 남아 있으면 반영이 절반만 간
    것이고, 그것이 이 검사가 잡는 것이다.
    """
    n = 0
    trees = [("paper", PAPER)]
    if os.path.isdir(DEPLOY):
        trees.append(("deploy/papers", DEPLOY))
    else:
        notes.append("G79 deploy/papers not present; only paper/ checked")
    for name, stale in PARITY:
        for tag, base in trees:
            for b, _, fs in os.walk(base):
                for f in sorted(fs):
                    if not f.endswith((".md", ".tex")):
                        continue
                    if stale in read(os.path.join(b, f)):
                        fails.append(
                            f"G79 {tag}/{f} still carries the signature of "
                            f"[{name}] -- the fix reached one tree only")
                        n += 1
    return n



# ----------------------------------------------------------------- G80
MEAS_TEX = re.compile(
    r"\\begin\{(measurement|observation)\}(.*?)\\end\{\1\}", re.S)
MEAS_MD = re.compile(
    r"^#### (?:Measurement|Observation)\b[^\n]*\n(.*?)(?=^#{2,4} |\Z)",
    re.S | re.M)
SCRIPT_IN = re.compile(r"([A-Za-z0-9_\\]+\.py)")


def _result_text(base, tree):
    """그 스크립트의 결과 파일 본문. 없으면 None."""
    stem = os.path.splitext(base)[0] + ".txt"
    for cand in (os.path.join(tree, "results", stem),
                 os.path.join(RESULTS, stem)):
        if os.path.exists(cand):
            return read(cand)
    for pas in range(1, 10):
        cand = os.path.join(VERIFY, "pass%d" % pas, "results", stem)
        if os.path.exists(cand):
            return read(cand)
    return None


def g80_cited_script_produces_its_numbers():
    """측정이 인쇄한 수는 그 측정이 이름 부른 스크립트에서 나와야 한다.

    G11 은 인쇄된 소수가 results/ 어딘가에 있는지만 본다. 어딘가에
    있으면서 **가리킨 곳에는 없는** 경우가 남고, pass6 이 그것을 찾았다:
    한 measurement 가 다섯 자리 유효숫자 다섯 개를 인쇄하면서 그 수를
    내지 못하는 스크립트 둘을 가리킨다. 배포 패킷만 받은 독자가 인용된
    스크립트를 돌리면 다른 목록을 얻는다.

    규칙: measurement 나 observation 이 스크립트를 하나 이상 이름으로
    부르면, 그 블록이 인쇄한 소수점 아래 세 자리 이상 리터럴은 전부
    그 스크립트들 중 하나의 결과 파일에 있어야 한다. 스크립트를 아무도
    안 부르면 이 검사는 발동하지 않는다 -- 그때는 G11 이 본다.
    """
    n = 0
    jobs = [(PAPER, MEAS_MD, ROOT, "paper")]
    if os.path.isdir(DEPLOY):
        jobs.append((DEPLOY, MEAS_TEX, os.path.dirname(DEPLOY),
                     "deploy/papers"))
    for base, pat, tree, tag in jobs:
        for b, _, fs in os.walk(base):
            for f in sorted(fs):
                if not f.endswith((".md", ".tex")):
                    continue
                src = read(os.path.join(b, f))
                for m in pat.finditer(src):
                    body = m.group(m.lastindex)
                    names = set()
                    for sm in SCRIPT_IN.finditer(body):
                        names.add(os.path.basename(
                            sm.group(1).replace("\\", "")))
                    if not names:
                        continue
                    texts = [t for t in (_result_text(x, tree)
                                         for x in sorted(names))
                             if t is not None]
                    if not texts:
                        continue
                    blob = "\n".join(texts)
                    # 논문은 결과 파일보다 적은 자리로 인쇄해도 된다.
                    # 그러니 문자열 일치만 보면 옳게 반올림한 자리까지
                    # 걸린다 -- 자기 자릿수로 반올림해 맞는지도 본다.
                    vals = []
                    for tok in re.findall(
                            r"[-+]?\d+\.\d+(?:[eE][-+]?\d+)?", blob):
                        try:
                            vals.append(float(tok))
                        except ValueError:
                            pass
                    missing = []
                    for d in sorted(set(DEC.findall(body))):
                        if d in blob:
                            continue
                        k = len(d.split(".")[1])
                        if any(("%.*f" % (k, v)) == d for v in vals):
                            continue
                        missing.append(d)
                    if missing:
                        fails.append(
                            f"G80 {tag}/{f}: a measurement citing "
                            f"{', '.join(sorted(names))} prints "
                            f"{len(missing)} decimals none of those "
                            f"scripts produce, e.g. "
                            f"{', '.join(missing[:3])}")
                        n += 1
    return n



# ----------------------------------------------------------------- G81
DEPLOY_ROOT = os.path.dirname(DEPLOY)


def _sha(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _repo_evidence_index():
    """저장소가 아는 증거 파일 전부: 이름 -> 해시 집합."""
    idx = {}
    roots = [CODE, RESULTS]
    for pas in sorted(os.listdir(VERIFY)) if os.path.isdir(VERIFY) else []:
        for sub in ("code", "results"):
            d = os.path.join(VERIFY, pas, sub)
            if os.path.isdir(d):
                roots.append(d)
    for r in roots:
        for b, _, fs in os.walk(r):
            if "__pycache__" in b:
                continue
            for f in fs:
                idx.setdefault(f, set()).add(_sha(os.path.join(b, f)))
    return idx


def g81_deployed_evidence_is_accounted_for():
    """배포된 증거는 저장소의 어느 파일이고, 인용된 것은 배포돼 있어야 한다.

    배포 패킷은 저장소에서 잘라 낸 것이므로, 그 안의 모든 코드와 결과는
    저장소 어딘가에 **같은 바이트로** 있어야 한다. 그렇지 않은 파일은
    패킷에서만 자라난 것이고, 저장소가 그것을 재현할 수 없다. 반대
    방향도 본다: 배포본 논문이 이름 부른 스크립트가 패킷 안에 없으면
    패킷만 받은 독자는 그 수를 확인할 길이 없다 -- pass6 이 그 자리를
    하나 찾았다.
    """
    n = 0
    if not os.path.isdir(DEPLOY_ROOT):
        return 0
    idx = _repo_evidence_index()
    for sub in ("code", "results"):
        d = os.path.join(DEPLOY_ROOT, sub)
        if not os.path.isdir(d):
            continue
        for f in sorted(os.listdir(d)):
            p = os.path.join(d, f)
            if not os.path.isfile(p) or f.endswith(".pyc"):
                continue
            if _sha(p) not in idx.get(f, set()):
                fails.append(
                    f"G81 deploy/{sub}/{f} is not in the repository with "
                    f"the same bytes; the packet cannot be regenerated "
                    f"from what the repository keeps")
                n += 1
    named = set()
    for b, _, fs in os.walk(DEPLOY):
        for f in sorted(fs):
            if not f.endswith(".tex"):
                continue
            for m in re.finditer(r"\\texttt\{([A-Za-z0-9_\\]+\.py)\}",
                                 read(os.path.join(b, f))):
                named.add((f, os.path.basename(
                    m.group(1).replace("\\", ""))))
    for tex, base in sorted(named):
        s = os.path.join(DEPLOY_ROOT, "code", base)
        r = os.path.join(DEPLOY_ROOT, "results",
                         os.path.splitext(base)[0] + ".txt")
        if not os.path.exists(s):
            fails.append(f"G81 deploy/papers/{tex} names {base}, which is "
                         f"not in the packet")
            n += 1
        elif not os.path.exists(r) and base != "indep.py":
            fails.append(f"G81 deploy/papers/{tex} names {base} but its "
                         f"result file is not in the packet")
            n += 1
    return n



def main():
    docs = [(p, read(p)) for p in paper_files()]
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
        ("G9 no stray CR", g9_line_endings()),
        ("G10 named scripts exist", g10_named_scripts(docs)),
        ("G11 printed decimals backed", g11_numbers_backed(docs)),
        ("G12 refutations disclosed", g12_refutations_disclosed(docs)),
        ("G13 lab measurements declare a null", g13_lab_nulls()),
        ("G14 evidence names its statement", g14_evidence_relevant(docs)),
        ("G15 cross-references resolve", g15_refs_resolve(docs)),
        ("G16 every sieve variant is audited", g16_sieve_manifest()),
        ("G17 N-dependent thresholds are computed",
         g17_no_typed_thresholds()),
        ("G18 results are not stale", g18_results_not_stale()),
        ("G19 no hand-written measurements in output",
         g19_no_typed_measurements()),
        ("G20 Euler products use a fixed bound",
         g20_euler_bound()),
        ("G21 cross-checked statistics agree",
         g21_cross_check()),
        ("G22 consumers run after what they read",
         g22_consumer_order()),
        ("G23 every source is non-empty and parses",
         g23_sources_compile()),
        ("G24 fitted exponents are swept",
         g24_exponents_swept()),
        ("G25 values from elsewhere are read, not copied",
         g25_no_copied_values()),
        ("G26 declined nulls have an auditor",
         g26_declines_are_audited()),
        ("G27 every pre-registered prediction is adjudicated",
         g27_predictions_adjudicated()),
        ("G28 forecasts carry a bracket",
         g28_forecasts_carry_a_bracket()),
        ("G29 octave fits have no unbounded bin",
         g29_octave_fits_are_bounded()),
        ("G30 octave fits declare their thinnest bin",
         g30_octave_fits_declare_population()),
        ("G31 octave fits declare their correlation",
         g31_octave_fits_declare_correlation()),
        ("G32 reported levels declare their budget",
         g32_levels_declare_their_budget()),
        ("G33 forecasts declare the drift they extrapolate",
         g33_forecasts_declare_drift()),
        ("G34 reported levels declare their arithmetic",
         g34_levels_declare_their_arithmetic()),
        ("G35 slope forecasts declare their scatter",
         g35_slope_forecasts_declare_scatter()),
        ("G36 slope forecasts declare their shape count",
         g36_slope_forecasts_declare_shapes()),
        ("G37 spans declare their noise floor",
         g37_spans_declare_their_noise_floor()),
        ("G38 published slopes declare their significance",
         g38_published_slopes_declare_significance()),
        ("G39 t-ratios declare the range they were measured over",
         g39_tstats_declare_their_range()),
        ("G40 single-point crossings beat their own floor",
         g40_crossings_are_wider_than_their_floor()),
        ("G41 censored crossings are counted and declared",
         g41_censored_crossings_are_declared()),
        ("G42 two per-N ranges require a paired ratio",
         g42_two_ranges_require_a_paired_ratio()),
        ("G43 spans declare how many scales they cover",
         g43_spans_declare_how_many_scales()),
        ("G44 quoted means declare the series behind them",
         g44_quoted_means_declare_their_series()),
        ("G45 cross-radical spans are judged against the sweep",
         g45_cross_radical_spans_are_judged_against_the_declared_sweep()),
        ("G46 counts are judged against an exchangeable null",
         g46_counts_are_judged_against_an_exchangeable_null()),
        ("G47 count references declare the magnitude one",
         g47_count_references_declare_the_magnitude_one()),
        ("G48 reference levels declare their own trend",
         g48_reference_levels_declare_their_own_trend()),
        ("G49 delegated floors cover the range they are used on",
         g49_delegated_floors_cover_the_range_they_are_used_on()),
        ("G50 mechanisms declare the window they were shown on",
         g50_mechanisms_declare_the_window_they_were_shown_on()),
        ("G51 disjoint windows are adjudicated on the target",
         g51_disjoint_windows_are_adjudicated_on_the_target()),
        ("G52 declined nulls meet a one-sided randomisation",
         g52_declined_nulls_meet_a_one_sided_randomisation()),
        ("G53 shape survival is redone on every point",
         g53_shape_survival_is_redone_on_every_point()),
        ("G54 predictors are ranked on a declared score",
         g54_predictors_are_ranked_on_a_declared_score()),
        ("G55 predictors declare their sieve level",
         g55_predictors_declare_their_sieve_level()),
        ("G56 thresholds follow from the declared statistic",
         g56_thresholds_follow_from_the_declared_statistic()),
        ("G57 one axis carries every statistic's threshold",
         g57_one_axis_carries_every_statistic_s_threshold()),
        ("G58 negligible remainders are accounted for",
         g58_negligible_remainders_are_accounted_for()),
        ("G59 capped spreads declare their sampling error",
         g59_capped_spreads_declare_their_sampling_error()),
        ("G60 shape gaps are read against their own error",
         g60_shape_gaps_are_read_against_their_own_error()),
        ("G61 whole-range gains declare their split",
         g61_whole_range_gains_declare_their_split()),
        ("G62 mass splits declare their overlap with the range",
         g62_mass_splits_declare_their_overlap_with_the_range()),
        ("G63 residual spreads are reported at every scale",
         g63_residual_spreads_are_reported_at_every_scale()),
        ("G64 shape forecasts declare their trust range",
         g64_shape_forecasts_declare_their_trust_range()),
        ("G65 frozen constants declare the other convention",
         g65_frozen_constants_declare_the_other_convention()),
        ("G66 forecasts on frozen constants declare both",
         g66_forecasts_on_frozen_constants_declare_both()),
        ("G67 maximised constants are checked off axis",
         g67_maximised_constants_are_checked_off_axis()),
        ("G68 subset claims declare the share they carry",
         g68_subset_claims_declare_the_share_they_carry()),
        ("G69 two-regressor fits declare their collinearity",
         g69_two_regressor_fits_declare_their_collinearity()),
        ("G70 every citation exists, not just the numbered ones",
         g70_every_citation_exists(docs)),
        ("G71 residual lists declare their sign run",
         g71_residual_signs_are_declared()),
        ("G72 within/between splits declare their resolution",
         g72_within_between_splits_declare_their_resolution()),
        ("G73 agreements declare the predictor's marginals",
         g73_agreements_declare_the_predictor_s_marginals()),
        ("G74 parametric families declare their classes",
         g74_parametric_families_declare_their_coprimality_classes()),
        ("G75 tolerances follow the printing they judge",
         g75_tolerances_follow_the_printing_they_judge()),
        ("G76 declared reads match their source",
         g76_reads_match_their_source()),
        ("G77 resolution rules name the unresolved case",
         g77_resolution_rules_name_the_unresolved()),
        ("G78 constants reach their result file",
         g78_constants_reach_their_result()),
        ("G79 one finding, both trees",
         g79_one_finding_both_trees()),
        ("G80 a cited script produces its numbers",
         g80_cited_script_produces_its_numbers()),
        ("G81 deployed evidence is accounted for",
         g81_deployed_evidence_is_accounted_for()),
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
