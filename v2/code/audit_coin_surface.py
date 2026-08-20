# -*- coding: utf-8 -*-
r"""
How much of this paper stands on having escaped a coin?

WHAT IS AT STAKE

rem:cncoindeep found that an iid coin resolved to 512 draws still
excludes mu where a resolved multiplicative ensemble covers it, so the
iid control lem:coin uses is narrower than the right null wherever the
object is multiplicative.  It also fixed the direction of the damage:
a narrower null is conservative for killing -- anything the coin
already covers, a wider null covers too -- so every claim lem:coin
killed stays killed.  What is not safe is the other direction, a
measurement recorded as having escaped the coin.

That bound is only useful with a count.  A repository-wide caveat that
touches three remarks is a footnote; one that touches thirty is a
rewrite.  **Nothing so far has counted.**  This does.

THE CENSUS AND ITS LIMITS

Claims live in paper/, so the paper is what is scanned, split into its
remarks and numbered statements.  A block counts as coin-calibrated if
it mentions a coin or an ensemble, and as claiming an escape if it
also contains one of a fixed list of phrases, set below before the
run:

    no coin draw, none of the, no draw of, 0 of, outside its coin,
    outside the coin, outside its control, cannot make, does not
    reach, does not go where

**This is a keyword census and it is a lower bound.**  A block phrased
some other way is missed; a block that uses a phrase while saying the
opposite is a false hit.  Every hit is printed with the phrase that
caught it and the label it sits in, so the list is a set of pointers
to read rather than a verdict, exactly as the constants census of
G78 was.

BACKS: Remark {#rem:coinsurface} in paper/wall_v3.md.

PRE-REGISTERED PREDICTIONS (written before this script was run)

  H1  THE GATE, a positive and a negative control.  rem:cnclass and
      rem:cncoindeep are found as escapes -- both state one in so many
      words -- and rem:cnshift is not, since it found the field inside
      its control at both octaves.  If the search cannot separate
      those three it cannot separate anything.
  H2  The coin is pervasive: at least 15 blocks of the paper mention
      a coin or an ensemble.
  H3  **The escapes are few: fewer than 10 blocks claim one.**  That
      is the bound on how much of the paper rem:cncoindeep's caveat
      touches.
  H4  And they are concentrated here: at least half of the escapes
      are in the C(N) family, whose labels begin rem:cn.  The surface
      outside this branch is then smaller still.

REFUTATION RULE (fixed before the run)

  H1  REFUTED if either control fails; nothing below is reported,
      because a search that cannot tell an escape from its opposite
      is not measuring what the count is supposed to mean.
  H2  REFUTED below 15.  Then the coin is not pervasive and the
      question this run exists for is smaller than it looked -- which
      would itself be worth knowing and would make H3 easy.
  H3  **REFUTED at 10 or more.**  Then the caveat is not a footnote
      and the paper needs a pass rather than a remark.  This is the
      outcome that would cost the most and it is the one to report
      first if it comes.
  H4  REFUTED below half.  Then the affected claims are spread
      through the paper rather than sitting in one branch, and the
      list matters more than the count.

  WHAT THIS CANNOT DO.  It cannot tell whether an escape is from an
  iid coin or from something else, whether the object in that block
  is multiplicative, or whether the escape still matters after the
  block's own caveats.  Each of those is a reading of the block and
  the blocks are named so they can be read.
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
OUT = os.path.join(ROOT, "results", "audit_coin_surface.txt")

COINWORD = re.compile(r"\bcoin\b|\bensemble\b", re.I)
ESCAPE = (
    "no coin draw",
    "none of the",
    "no draw of",
    "0 of",
    "outside its coin",
    "outside the coin",
    "outside its control",
    "cannot make",
    "does not reach",
    "does not go where",
)
HEADER = re.compile(r"^#{3,5} .*?\{#([A-Za-z:]+[A-Za-z0-9:]*)\}", re.M)
ESCRE = tuple(re.compile(r"(?<![0-9A-Za-z])" + re.escape(p_), re.I)
              for p_ in ESCAPE)


def flatten(t):
    """the paper writes numbers inside $...$; the census must not"""
    return re.sub(r"\s+", " ", t.replace("$", " "))


POSCTL = ("rem:cnclass", "rem:cncoindeep")
NEGCTL = "rem:cnshift"
MINCOIN = 15
MAXESC = 10
FAMILY = "rem:cn"


HEAD = [
    "STATISTIC: per labelled block of paper/, whether it mentions a",
    "           coin or an ensemble and whether it also contains one",
    "           of ten fixed phrases that state an escape from one;",
    "           the counts of each and the list of escapes.",
    "FIELD: every block of every document in v2/paper/ that carries a",
    "       {#label}, split at the labelled headings, at the commit",
    "       this was run on. The phrase list is fixed in the script",
    "       and printed with each hit, so the census is a lower bound",
    "       by construction: a block phrased otherwise is missed and a",
    "       block using a phrase to say the opposite is a false hit.",
    "NULL: none is run and none applies. This is a census of text on",
    "      disk, not a signal against a background.",
    "",
]


def blocks():
    out = []
    for f in sorted(os.listdir(PAPER)):
        if not f.endswith(".md"):
            continue
        src = io.open(os.path.join(PAPER, f), encoding="utf-8").read()
        hits = list(HEADER.finditer(src))
        for i, m in enumerate(hits):
            end = hits[i + 1].start() if i + 1 < len(hits) else len(src)
            out.append((f, m.group(1), src[m.start():end]))
    return out


def main():
    lines = []

    def say(t=""):
        print(t)
        sys.stdout.flush()
        lines.append(t)

    bs = blocks()
    say("  NOTE, disclosed: the first execution failed H1 -- the "
        "positive control")
    say("  rem:cncoindeep was missed because the paper writes its "
        "numbers inside")
    say("  dollar signs and this census read the markup literally. "
        "The fix is to")
    say("  the reader, not the rule: dollars are stripped and "
        "whitespace collapsed")
    say("  before matching, and each phrase is matched at a word "
        "boundary so that")
    say("  \"0 of\" does not fire inside \"10 of\". The phrase list "
        "is exactly the")
    say("  one registered above.")
    say("%d labelled blocks across %d documents"
        % (len(bs), len(set(f for f, _, _ in bs))))
    say("  the phrase list, fixed before the run and unchanged "
        "by the fix below: %s"
        % "; ".join(ESCAPE))

    coined, escapes = [], []
    for f, lab, body in bs:
        if not COINWORD.search(body):
            continue
        coined.append(lab)
        low = flatten(body)
        got = [ESCAPE[i] for i, rx in enumerate(ESCRE)
               if rx.search(low)]
        if got:
            escapes.append((lab, f, got))

    # -------------------------------------------------------------- H1
    say()
    say("H1  can the search tell an escape from its opposite?")
    labs = [l for l, _, _ in escapes]
    pos = [c for c in POSCTL if c in labs]
    neg = NEGCTL in labs
    say("  positive controls found: %s of %s"
        % (", ".join(pos) if pos else "none", ", ".join(POSCTL)))
    say("  negative control %s: %s"
        % (NEGCTL, "FLAGGED, which it should not be" if neg
           else "not flagged, as it should not be"))
    h1 = len(pos) == len(POSCTL) and not neg
    say("  H1 %s   (cap: both positives, neither negative)"
        % ("hold" if h1 else "REFUTED"))
    if not h1:
        io.open(OUT, "w", encoding="utf-8", newline="\n").write(
            "\n".join(HEAD + lines) + "\n")
        raise SystemExit(1)

    # -------------------------------------------------------------- H2
    say()
    say("H2  how pervasive is the coin?")
    say("  %d of %d blocks mention a coin or an ensemble"
        % (len(coined), len(bs)))
    h2 = len(coined) >= MINCOIN
    say("COUNT coinsurface_coined %d" % len(coined))
    say("  H2 %s   (cap: at least %d)"
        % ("hold" if h2 else "REFUTED", MINCOIN))

    # -------------------------------------------------------------- H3
    say()
    say("H3  how many claim an escape?")
    say("  %d of the %d coin blocks do" % (len(escapes), len(coined)))
    say("      label                          document"
        "                  phrase that caught it")
    for lab, f, got in escapes:
        say("      %-30s %-24s %s" % (lab, f, "; ".join(got)))
    say("SCALES 1")
    say("COUNT coinsurface_escapes %d" % len(escapes))
    h3 = len(escapes) < MAXESC
    say("  H3 %s   (cap: fewer than %d)"
        % ("hold" if h3 else "REFUTED", MAXESC))

    # -------------------------------------------------------------- H4
    say()
    say("H4  are they concentrated in this branch?")
    fam = [l for l, _, _ in escapes if l.startswith(FAMILY)]
    share = len(fam) / float(len(escapes)) if escapes else 0.0
    say("  %d of %d begin %s, a share of %.4f"
        % (len(fam), len(escapes), FAMILY, share))
    outside = [l for l, _, _ in escapes if not l.startswith(FAMILY)]
    say("  the ones outside it: %s"
        % (", ".join(outside) if outside else "none"))
    say("SHARE coinsurface_family %.4f" % share)
    h4 = share >= 0.5
    say("  H4 %s   (cap: at least half)"
        % ("hold" if h4 else "REFUTED"))

    say()
    say("=" * 70)
    say("H1 %s  H2 %s  H3 %s  H4 %s"
        % tuple("hold" if v else "REFUTED"
                for v in (h1, h2, h3, h4)))
    say()
    if h3 and h4:
        say("the caveat is a footnote and not a rewrite. The coin is "
            "everywhere in")
        say("this paper and almost everywhere it is killing something, "
            "which a")
        say("narrower null does conservatively. The blocks that stand "
            "on having")
        say("escaped one are few and mostly in the branch that found "
            "the problem;")
        say("the few outside it are named above and are what to read "
            "next.")
    elif not h3:
        say("the caveat is not a footnote. Enough of the paper stands "
            "on escaping")
        say("a coin that the right response is a pass over those "
            "blocks rather")
        say("than a remark, and the list above is that pass's work "
            "list.")
    else:
        say("the escapes are spread through the paper rather than "
            "concentrated")
        say("here, so the list matters more than the count and every "
            "block on it")
        say("has to be read on its own terms.")

    io.open(OUT, "w", encoding="utf-8", newline="\n").write(
        "\n".join(HEAD + lines) + "\n")
    print("\nwrote %s" % os.path.normpath(OUT))
    return 0


if __name__ == "__main__":
    sys.exit(main())
