# -*- coding: utf-8 -*-
r"""Extract the two sides into comparable inventories.

SOURCE  : v2/paper/{theorem_A,wall_v3}.md   -- markdown, "#### Remark (title) {#rem:key}"
PROJ    : deploy/papers/P{1..5}*.tex        -- LaTeX, \begin{env}[title]\label{key}

Writes machine-readable inventories to B-out/results/ so every later
script reads one file instead of re-parsing prose.

Usage: python extract.py            (writes inventory_source.json etc.)
"""
import io
import json
import os
import re
import sys

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

ROOT = r"z:\업무\goldbach-ln2-real"
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..",
                   "results")
OUT = os.path.abspath(OUT)

SOURCE = ["v2/paper/theorem_A.md", "v2/paper/wall_v3.md"]
PROJ = [
    "deploy/papers/P1-mobius-fixed-class.tex",
    "deploy/papers/P2-no-go-divisor-switch.tex",
    "deploy/papers/P3-wall-second-moment.tex",
    "deploy/papers/P4-coherent-cell-floor.tex",
    "deploy/papers/P5-negative-map.tex",
]

# ---------------------------------------------------------------- source

# "#### Remark (title) {#rem:key}"  /  "#### Proposition (title) {#prop:key}"
MD_HEAD = re.compile(
    r"^(#{2,4})\s+"
    r"(Remark|Proposition|Theorem|Lemma|Corollary|Conjecture|Definition|"
    r"Observation|Measurement|Note)"
    r"\s*(\([^\n]*?\))?\s*(\{#([A-Za-z0-9:_.-]+)\})?\s*$")

MD_ANYHEAD = re.compile(r"^(#{1,6})\s+(.*)$")


def read(rel):
    with io.open(os.path.join(ROOT, rel.replace("/", os.sep)),
                 encoding="utf-8") as f:
        return f.read().split("\n")


def parse_source():
    items = []
    for rel in SOURCE:
        lines = read(rel)
        heads = []
        for i, ln in enumerate(lines):
            m = MD_HEAD.match(ln)
            if m:
                heads.append((i, m))
        for j, (i, m) in enumerate(heads):
            # body runs to the next heading of any level
            end = len(lines)
            for k in range(i + 1, len(lines)):
                if MD_ANYHEAD.match(lines[k]):
                    end = k
                    break
            items.append({
                "file": rel,
                "line": i + 1,
                "kind": m.group(2),
                "title": (m.group(3) or "").strip("()") or "",
                "key": m.group(5) or "",
                "body": "\n".join(lines[i + 1:end]).strip(),
            })
    return items


# ------------------------------------------------------------ projection

TEX_BEGIN = re.compile(
    r"\\begin\{(note|observation|measurement|theorem|proposition|lemma|"
    r"corollary|conjecture|definition|remark)\}")
TEX_LABEL = re.compile(r"\\label\{([^}]*)\}")


def parse_proj():
    items = []
    for rel in PROJ:
        with io.open(os.path.join(ROOT, rel.replace("/", os.sep)),
                     encoding="utf-8") as f:
            lines = f.read().split("\n")
        for i, ln in enumerate(lines):
            m = TEX_BEGIN.search(ln)
            if not m:
                continue
            env = m.group(1)
            # header may wrap: gather until \label or the line after next
            blob = ln
            k = i
            while ("\\label{" not in blob) and (k + 1 < len(lines)) \
                    and (k - i < 3):
                k += 1
                blob += " " + lines[k]
            lab = TEX_LABEL.search(blob)
            # title in [...] possibly spanning lines
            t = re.search(r"\\begin\{%s\}\[(.*?)\]" % env, blob, re.S)
            end = len(lines)
            for q in range(i, len(lines)):
                if ("\\end{%s}" % env) in lines[q]:
                    end = q
                    break
            items.append({
                "file": rel,
                "line": i + 1,
                "env": env,
                "title": re.sub(r"\s+", " ", (t.group(1) if t else "")),
                "key": lab.group(1) if lab else "",
                "body": "\n".join(lines[i:end + 1]),
            })
    return items


def all_labels_proj():
    out = {}
    for rel in PROJ:
        with io.open(os.path.join(ROOT, rel.replace("/", os.sep)),
                     encoding="utf-8") as f:
            txt = f.read()
        for m in TEX_LABEL.finditer(txt):
            out.setdefault(m.group(1), []).append(rel)
    return out


def all_keys_source():
    out = {}
    for rel in SOURCE:
        with io.open(os.path.join(ROOT, rel.replace("/", os.sep)),
                     encoding="utf-8") as f:
            txt = f.read()
        for m in re.finditer(r"\{#([A-Za-z0-9:_.-]+)\}", txt):
            out.setdefault(m.group(1), []).append(rel)
    return out


def main():
    src = parse_source()
    prj = parse_proj()
    j = lambda n, o: io.open(os.path.join(OUT, n), "w",
                             encoding="utf-8", newline="\n").write(
        json.dumps(o, ensure_ascii=False, indent=1))
    j("inventory_source.json", src)
    j("inventory_proj.json", prj)
    j("labels_proj.json", all_labels_proj())
    j("labels_source.json", all_keys_source())

    from collections import Counter
    print("STATISTIC: counts of numbered environments on each side, by kind")
    print("FIELD:     v2/paper/*.md (source) and deploy/papers/*.tex "
          "(projection)")
    print("CONSTANTS: none")
    print("NULL:      n/a -- this is an inventory, not a test")
    print("DENOM:     every heading matching the kind regex in each file")
    print()
    print("SOURCE by kind:")
    for k, v in sorted(Counter(x["kind"] for x in src).items()):
        print("   %-14s %d" % (k, v))
    print("   TOTAL          %d" % len(src))
    print()
    print("SOURCE remarks per file:")
    for rel in SOURCE:
        n = sum(1 for x in src if x["file"] == rel and x["kind"] == "Remark")
        print("   %-24s %d" % (os.path.basename(rel), n))
    print()
    print("PROJECTION by env:")
    for k, v in sorted(Counter(x["env"] for x in prj).items()):
        print("   %-14s %d" % (k, v))
    print("   TOTAL          %d" % len(prj))
    print()
    print("PROJECTION per file:")
    for rel in PROJ:
        n = sum(1 for x in prj if x["file"] == rel)
        print("   %-32s %d" % (os.path.basename(rel), n))
    print()
    print("source anchors {#...}: %d   projection \\label{}: %d"
          % (len(all_keys_source()), len(all_labels_proj())))


if __name__ == "__main__":
    main()
