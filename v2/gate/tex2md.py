# -*- coding: utf-8 -*-
"""
LaTeX -> Markdown. 결정적이고 다시 돌릴 수 있는 변환.

    python gate/tex2md.py paper/wall_v3.tex > paper/wall_v3.md

손으로 옮기지 않는 이유: 이 세션의 오류 상당수가 백슬래시·이스케이프에서
나왔다. 변환은 스크립트가 하고, 결과는 사람이 본다.

출력 규약 (게이트가 이걸 읽는다):

    ### Theorem (w_k = 1) {#thm:A}
    <!-- evidence: analytic -->

수식은 `$...$` / `$$...$$` 그대로 둔다 -- 마크다운에서도 수학은 LaTeX다.
표와 변환이 애매한 덩어리는 삼중 백틱 `latex` 블록으로 **그대로** 남긴다.
조용히 뭉개는 것보다 안 건드린 채 보이는 게 낫다.
"""

import io
import os
import re
import sys

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

KINDS = {
    "theorem": "Theorem", "proposition": "Proposition", "lemma": "Lemma",
    "corollary": "Corollary", "conjecture": "Conjecture",
    "remark": "Remark", "definition": "Definition",
}

ACCENTS = [
    (r'\\"o', "ö"), (r'\\"a', "ä"), (r'\\"u', "ü"), (r'\\"O', "Ö"),
    (r'\\"A', "Ä"), (r'\\"U', "Ü"), (r"\\'e", "é"), (r"\\`e", "è"),
    (r"\\l\{\}", "ł"), (r"\\l\b", "ł"), (r"\\ss\b", "ß"),
    (r"\\c\{c\}", "ç"), (r"\\i\b", "ı"),
]


def protect_math(text, store):
    """수식을 자리표시자로 빼둔다. 그 안에서는 아무것도 치환하지 않는다."""
    def keep(m):
        store.append(m.group(0))
        return f"\x00M{len(store) - 1}\x00"
    text = re.sub(r"\\\[(.*?)\\\]", keep, text, flags=re.S)
    text = re.sub(r"\$\$(.*?)\$\$", keep, text, flags=re.S)
    text = re.sub(r"(?<!\\)\$(?:[^$\\]|\\.)+\$", keep, text)
    text = re.sub(r"\\begin\{(equation|align|gather)\*?\}.*?"
                  r"\\end\{\1\*?\}", keep, text, flags=re.S)
    return text


def restore_math(text, store):
    def put(m):
        raw = store[int(m.group(1))]
        d = re.match(r"\\\[(.*?)\\\]$", raw, re.S)
        if d:
            return "\n$$\n" + d.group(1).strip() + "\n$$\n"
        if raw.startswith("\\begin{"):
            return "\n$$\n" + raw.strip() + "\n$$\n"
        return raw
    return re.sub(r"\x00M(\d+)\x00", put, text)


def keep_verbatim(text, store, env):
    """표처럼 손대면 안 되는 환경을 통째로 보존한다."""
    def keep(m):
        store.append(m.group(0))
        return f"\x00V{len(store) - 1}\x00"
    return re.sub(r"\\begin\{" + env + r"\*?\}.*?\\end\{" + env + r"\*?\}",
                  keep, text, flags=re.S)


def convert(src):
    # 본문만
    m = re.search(r"\\begin\{document\}(.*?)\\end\{document\}", src, re.S)
    body = m.group(1) if m else src
    preamble = src[:m.start()] if m else ""

    verb = []
    for env in ("longtable", "tabular", "table", "center", "verbatim"):
        body = keep_verbatim(body, verb, env)

    math = []
    body = protect_math(body, math)

    # 주석 제거 (evidence 마커는 살린다)
    body = re.sub(r"(?<!\\)%(?!\s*evidence:).*", "", body)

    # 제목 구조
    body = re.sub(r"\\section\*?\{([^}]*)\}\s*\\label\{([^}]*)\}",
                  r"\n## \1 {#\2}\n", body)
    body = re.sub(r"\\section\*?\{([^}]*)\}", r"\n## \1\n", body)
    body = re.sub(r"\\subsection\*?\{([^}]*)\}\s*\\label\{([^}]*)\}",
                  r"\n### \1 {#\2}\n", body)
    body = re.sub(r"\\subsection\*?\{([^}]*)\}", r"\n### \1\n", body)
    body = re.sub(r"\\paragraph\{([^}]*)\}", r"\n**\1.**", body)

    # 진술 환경
    for env, name in KINDS.items():
        def open_stmt(mm, name=name):
            opt = (mm.group(1) or "").strip("[]")
            lab = mm.group(2)
            head = f"{name} ({opt})" if opt else name
            return f"\n#### {head} {{#{lab}}}\n"
        body = re.sub(r"\\begin\{" + env + r"\}(\[[^\]]*\])?\s*"
                      r"\\label\{([^}]*)\}", open_stmt, body)
        body = re.sub(r"\\begin\{" + env + r"\}(\[[^\]]*\])?",
                      lambda mm, n=name: f"\n#### {n}"
                      + (f" ({mm.group(1).strip('[]')})" if mm.group(1)
                         else "") + "\n", body)
        body = re.sub(r"\\end\{" + env + r"\}", "\n", body)

    body = re.sub(r"\\begin\{proof\}(\[[^\]]*\])?",
                  lambda mm: "\n**Proof"
                  + (f" ({mm.group(1).strip('[]')})" if mm.group(1) else "")
                  + ".** ", body)
    body = re.sub(r"\\end\{proof\}", " ∎\n", body)
    body = re.sub(r"\\begin\{abstract\}", "\n## Abstract\n", body)
    body = re.sub(r"\\end\{abstract\}", "\n", body)

    # 목록
    body = re.sub(r"\\begin\{itemize\}", "\n", body)
    body = re.sub(r"\\end\{itemize\}", "\n", body)
    body = re.sub(r"\\begin\{enumerate\}", "\n", body)
    body = re.sub(r"\\end\{enumerate\}", "\n", body)
    body = re.sub(r"\\item\s*\[([^\]]*)\]", r"\n- **\1** ", body)
    body = re.sub(r"\\item\s*", "\n- ", body)

    # 인라인 서식
    body = re.sub(r"\\(?:emph|textit)\{([^{}]*)\}", r"*\1*", body)
    body = re.sub(r"\\textbf\{([^{}]*)\}", r"**\1**", body)
    body = re.sub(r"\\texttt\{([^{}]*)\}", r"`\1`", body)
    body = re.sub(r"\\ref\{([^}]*)\}", r"[\1]", body)
    body = re.sub(r"\\eqref\{([^}]*)\}", r"[\1]", body)
    body = re.sub(r"\\cite\{([^}]*)\}", r"[\1]", body)
    body = re.sub(r"\\label\{([^}]*)\}", r"{#\1}", body)
    body = re.sub(r"\\footnote\{([^{}]*)\}", r" (\1)", body)

    # 서지
    body = re.sub(r"\\begin\{thebibliography\}\{[^}]*\}",
                  "\n## References\n", body)
    body = re.sub(r"\\end\{thebibliography\}", "\n", body)
    body = re.sub(r"\\bibitem\{([^}]*)\}", r"\n- **[\1]** ", body)

    # 남은 것들
    body = re.sub(r"\\(?:maketitle|tableofcontents|newpage|noindent|"
                  r"medskip|smallskip|bigskip|centering|small|hline|"
                  r"endhead)\b", "", body)
    body = re.sub(r"\\(?:title|author|date)\{.*?\}", "", body, flags=re.S)
    for pat, rep in ACCENTS:
        body = re.sub(pat, rep, body)
    body = body.replace("---", "—").replace("--", "–")
    body = body.replace("~", " ")          # LaTeX 비분할 공백
    body = re.sub(r"\\,|\\ |\\;|\\!|\\quad|\\qquad", " ", body)
    body = re.sub(r"``|''", '"', body)
    body = re.sub(r"\\\\(?=\s*\n)", "", body)

    body = restore_math(body, math)

    # 중첩된다: table 안에 center, 그 안에 tabular. 한 번만 치환하면
    # 안쪽 자리표시자가 그대로 남아 널 바이트가 파일에 들어간다.
    def put_verb(mm):
        return "\n```latex\n" + verb[int(mm.group(1))].strip() + "\n```\n"
    for _ in range(8):
        new = re.sub(r"\x00V(\d+)\x00", put_verb, body)
        if new == body:
            break
        body = new
    body = re.sub(r"\n```latex\n(\s*)```latex\n", r"\n```latex\n\1", body)
    body = re.sub(r"\n```\n(\s*)```\n", r"\n```\n", body)

    body = re.sub(r"\n{3,}", "\n\n", body).strip()

    # 매크로는 보존한다 -- 수식이 이것들을 쓴다
    macros = re.findall(r"\\(?:newcommand|renewcommand|DeclareMathOperator)"
                        r"\*?\{[^\n]*", preamble)
    title = re.search(r"\\title\{(.*?)\}\s*\n", src, re.S)
    if title:
        t = re.sub(r"\\\\|\n|\s+", " ", title.group(1)).strip()
        t = t.replace("---", "—").replace("--", "–")
        for pat, rep in ACCENTS:
            t = re.sub(pat, rep, t)
        head = "# " + t
    else:
        head = "# (제목 없음)"

    out = [head, ""]
    if macros:
        out += ["```latex", "% 수식이 쓰는 매크로 — 렌더러/역변환용",
                *macros, "```", ""]
    out += [body, ""]
    return "\n".join(out)


def main():
    if len(sys.argv) < 2:
        print("usage: tex2md.py <file.tex>", file=sys.stderr)
        return 2
    # 줄끝은 읽는 즉시 LF로 눌러 버린다. newline=""로 읽어 CRLF를 안고 가면,
    # 텍스트 모드 stdout이 스크립트가 넣은 \n만 CRLF로 바꾸므로 원래의
    # \r\n이 \r\r\n이 된다. 그러면 파일 안에 홀CR이 남고, splitlines()가
    # 그것도 줄바꿈으로 세서 게이트의 줄 번호가 실제와 어긋난다.
    src = io.open(sys.argv[1], encoding="utf-8", newline="").read()
    src = src.replace("\r\n", "\n").replace("\r", "\n")
    # 그리고 stdout이 아니라 바이너리로 쓴다 -- 리다이렉트가 다시 번역하지
    # 못하게.
    sys.stdout.flush()
    sys.stdout.buffer.write(convert(src).encode("utf-8"))
    return 0


if __name__ == "__main__":
    sys.exit(main())
