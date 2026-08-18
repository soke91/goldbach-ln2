# -*- coding: utf-8 -*-
r"""
The other two destinations, and what the data actually bound.

WHAT IS AT STAKE

{#rem:shapepower} measured the power of this repository's shape
discriminator on F and found none: fed synthetic sweeps drawn from
the fitted shape itself, it picks the generating shape in 0.0000 of
trials, and the span that would give power 0.95 puts the top N at
10^15. It closed by naming what it did not do -- {#rem:splitreach}'s
one-sign fraction and {#rem:headaxis}'s axis spreads rest on the same
kind of comparison and were not put through the test. This does that,
and then asks a better question than the one the discriminator asks.

Model selection is not what either open question needs. Nobody wants
to know which of two curves fits; they want to know WHERE the trend
goes -- whether the head's one sign reaches a coin, and whether the
imbalance axis hands over to mass. That is a statement about an
asymptote, and an asymptote has a confidence interval whether or not
a shape contest has a winner. Fitting y = L + b*g(log N) for a bounded
g makes L the destination and gives it a standard error directly.

The catch is that L depends on g, and that dependence is the shape
ambiguity in the one place it can be quantified: the spread of L
across bases, against the error within one. If the spread dominates,
the destination is undetermined and it is the basis choice, not the
noise, that leaves it so.

BACKS: Remark {#rem:destination} in paper/theorem_A.md.

PRE-REGISTERED PREDICTIONS (written before this script was run)

  Q1  The control. The linear slopes refitted here reproduce the
      published ones inside the bound the tables' printing forces,
      computed here rather than assumed -- the lesson
      {#rem:shapepower} paid for.
  Q2  The transfer holds. On both quantities the discriminator's
      power, at their own abscissae and residual r.m.s., is below
      0.95.
  Q3  The destination is undetermined, and by basis and not by
      noise: across the three bounded bases 1/x, 1/sqrt(x) and 1/x^2
      the spread of the fitted asymptote exceeds twice the largest
      within-basis standard error, for each quantity.
  Q4  And no bounded fit sends the head to a coin: on every basis the
      one-sign asymptote's two-sigma interval lies above one half.

REFUTATION RULE (fixed before the run)

  Q1  REFUTED if any refit leaves its own computed bound. Then this
      is not the fit those files published. THIS ONE GATES.
  Q2  REFUTED if either power reaches 0.95. Then that quantity's
      shape question IS answerable at the reach already computed, and
      it should be answered rather than left open.
  Q3  REFUTED if the spread is within twice the largest within-basis
      error. Then the bases agree on the destination, the asymptote
      is determined to that precision, and the open question closes
      with a number instead of a refusal.
  Q4  REFUTED if any basis's interval reaches one half. Then a coin
      is inside what the data allow, {#rem:splitreach}'s "nowhere
      near a coin" is a statement about the measured range only, and
      the head may be dissolving after all.

  Q1 gates. Q2 to Q4 are the measurement and do not gate.

  THE NULL IS THE MEASUREMENT for Q2, as in {#rem:shapepower}: every
  trial draws from a shape known in advance and the number reported is
  how often the test recovers it. For Q3 and Q4 no null applies --
  they are interval estimates of a parameter, and the interval is its
  own statement of what the data exclude. The draws use the fixed SEED
  declared below.
"""

import io
import math
import os
import re
import sys

import numpy as np

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
RES = os.path.join(ROOT, "results")
OUT = os.path.join(RES, "audit_power_reach.txt")

SEED = 20260819
TRIALS = 2000
TARGET = 0.95
HALF = 0.5

SPLITROW = re.compile(r"^  (\d{5,})\s+\d+\s+\d+\s+[\d.]+\s+[\d.]+\s+"
                      r"[\d.]+\s+[\d.]+\s+([\d.]+)\s+[\d.]+\s*$", re.M)
AXISROW = re.compile(r"^  (\d{5,})\s+\d+\s+[\d.]+\s+([+-][\d.]+)\s+"
                     r"([+-][\d.]+)\s+([+-][\d.]+)\s+[\d.]+\s+"
                     r"[\d.]+\s+[+-]?[\d.]+\s*$", re.M)


def decimals(s):
    return len(s.split(".")[1]) if "." in s else 0


def read_inputs():
    src = io.open(os.path.join(RES, "audit_split_reach.txt"),
                  encoding="utf-8").read()
    one, dec1 = [], 0
    for m in SPLITROW.finditer(src):
        one.append((int(m.group(1)), float(m.group(2))))
        dec1 = max(dec1, decimals(m.group(2)))
    p1 = float(re.search(r"least-squares slope in log N = "
                         r"([+-][\d.]+)", src).group(1))
    src2 = io.open(os.path.join(RES, "audit_headaxis_reach.txt"),
                   encoding="utf-8").read()
    spI, spT, dec2 = [], [], 0
    for m in AXISROW.finditer(src2):
        spI.append((int(m.group(1)), float(m.group(2))))
        spT.append((int(m.group(1)), float(m.group(3))))
        dec2 = max(dec2, decimals(m.group(2)))
    p2 = float(re.search(r"least-squares slope in log N = "
                         r"([+-][\d.]+)", src2).group(1))
    p3 = float(re.search(r"the T spread's slope is ([+-][\d.]+)",
                         src2).group(1))
    return one, dec1, p1, spI, spT, dec2, p2, p3


def linfit(x, y):
    a, b = np.polyfit(x, y, 1)
    r = y - (a * x + b)
    se = math.sqrt(float((r ** 2).sum() / (x.size - 2))
                   / float(((x - x.mean()) ** 2).sum()))
    return float(a), float(b), float(np.sqrt((r ** 2).mean())), se


def basisfit(x, y, g):
    """y = L + b g(x); returns L, its standard error, and the r.m.s."""
    u = g(x)
    A = np.column_stack([np.ones_like(u), u])
    c, *_ = np.linalg.lstsq(A, y, rcond=None)
    res = y - A.dot(c)
    n = x.size
    s2 = float((res ** 2).sum()) / (n - 2)
    cov = s2 * np.linalg.inv(A.T.dot(A))
    return (float(c[0]), math.sqrt(float(cov[0, 0])),
            float(np.sqrt((res ** 2).mean())), float(c[1]))


def rms_pair(x, y):
    """the repository's two shapes: linear in log N against L + b/x"""
    a, b, r1, _se = linfit(x, y)
    _L, _sL, r2, _bb = basisfit(x, y, lambda t: 1.0 / t)
    return r1, r2


def verdict(r1, r2, n):
    se = min(r1, r2) / math.sqrt(2.0 * (n - 2))
    return abs(r1 - r2) > se, se


def power(xs, truth, sigma, rng, trials, want_linear):
    hit = 0
    sep = 0
    for _ in range(trials):
        y = truth + rng.normal(0.0, sigma, size=truth.size)
        r1, r2 = rms_pair(xs, y)
        ok, _se = verdict(r1, r2, xs.size)
        if ok:
            sep += 1
            if (r1 < r2) == want_linear:
                hit += 1
    return hit / float(trials), sep / float(trials)


BASES = (("1/x", lambda t: 1.0 / t),
         ("1/sqrt(x)", lambda t: 1.0 / np.sqrt(t)),
         ("1/x^2", lambda t: 1.0 / (t * t)))


def main():
    lines = []

    def say(s=""):
        print(s)
        lines.append(s)

    one, dec1, pub1, spI, spT, dec2, pub2, pub3 = read_inputs()
    say("read %d one-sign fractions from "
        "results/audit_split_reach.txt (slope %+.6f published)"
        % (len(one), pub1))
    say("  and %d axis spreads from "
        "results/audit_headaxis_reach.txt (|I| %+.6f, T %+.6f "
        "published)" % (len(spI), pub2, pub3))
    say("SEED %d" % SEED)
    say("TRIALS %d" % TRIALS)

    sets = [("one-sign fraction", one, dec1, pub1),
            ("|I| spread", spI, dec2, pub2),
            ("T spread", spT, dec2, pub3)]

    rng = np.random.default_rng(SEED)

    # -------------------------------------------------------------- Q1
    say()
    say("Q1  the controls, each against the bound its own printing "
        "forces")
    say("  quantity              slope here   published    bound      "
        "inside")
    q1 = True
    prepared = []
    for nm, rows, dec, pub in sets:
        x = np.log(np.array([r[0] for r in rows], dtype=np.float64))
        y = np.array([r[1] for r in rows])
        a, b, rms, se = linfit(x, y)
        rnd = 0.5 * 10.0 ** (-dec)
        cs = (x - x.mean()) / float(((x - x.mean()) ** 2).sum())
        bound = float(np.abs(cs).sum()) * rnd
        ins = abs(a - pub) <= bound
        if not ins:
            q1 = False
        prepared.append((nm, x, y, a, b, rms, se))
        say("  %-21s %+-12.6f %+-12.6f %-10.6f %s"
            % (nm, a, pub, bound, "yes" if ins else "NO"))
    say("  the printing carries %d and %d decimals, so a single value "
        "is within %.6f and %.6f of what produced it"
        % (dec1, dec2, 0.5 * 10.0 ** (-dec1), 0.5 * 10.0 ** (-dec2)))
    say("  Q1 %s   (cap: each slope's own computed bound)"
        % ("hold" if q1 else "REFUTED"))

    # -------------------------------------------------------------- Q2
    say()
    say("Q2  the power of the discriminator on these two")
    say("  quantity              span     noise      separates  "
        "picks it")
    q2 = True
    for nm, x, y, a, b, rms, se in prepared:
        truth = a * x + b
        pw, sp = power(x, truth, rms, rng, TRIALS, True)
        if pw >= TARGET:
            q2 = False
        say("  %-21s %-8.4f %-10.6f %-10.4f %.4f"
            % (nm, float(x.max() - x.min()), rms, sp, pw))
    say("  drawn from the linear shape itself at the observed "
        "abscissae and residual r.m.s.")
    say("  Q2 %s   (target power %.2f)"
        % ("hold" if q2 else "REFUTED", TARGET))

    # -------------------------------------------------------------- Q3
    say()
    say("Q3  where does each trend actually go, and how well is that "
        "pinned?")
    say("  the destination is the intercept L of y = L + b g(log N), "
        "so it is")
    say("  estimated with a standard error on every basis; the "
        "question is")
    say("  whether the bases agree better than they are each pinned.")
    q3 = True
    dests = {}
    for nm, x, y, a, b, rms, se in prepared:
        say("  %s" % nm)
        say("    basis        L            s.e.       r.m.s.")
        Ls, ses = [], []
        for bn, g in BASES:
            L, sL, r, _bb = basisfit(x, y, g)
            Ls.append(L)
            ses.append(sL)
            say("    %-12s %+-12.6f %-10.6f %.6f" % (bn, L, sL, r))
        spread = max(Ls) - min(Ls)
        worst = max(ses)
        dests[nm] = (Ls, ses)
        if not (spread > 2.0 * worst):
            q3 = False
        say("    spread across bases %.6f against twice the largest "
            "standard error %.6f  %s"
            % (spread, 2.0 * worst,
               "basis dominates" if spread > 2.0 * worst
               else "errors dominate"))
    say("  Q3 %s   (cap: twice the largest within-basis error)"
        % ("hold" if q3 else "REFUTED"))

    # -------------------------------------------------------------- Q4
    say()
    say("Q4  does any bounded fit send the head to a coin?")
    Ls, ses = dests["one-sign fraction"]
    say("  basis        L            two-sigma interval        "
        "clears one half")
    q4 = True
    for (bn, _g), L, sL in zip(BASES, Ls, ses):
        lo, hi = L - 2.0 * sL, L + 2.0 * sL
        cl = lo > HALF
        if not cl:
            q4 = False
        say("    %-12s %+-12.6f [%+.6f, %+.6f]   %s"
            % (bn, L, lo, hi, "yes" if cl else "NO"))
    say("  the cap is one half, which is what a coin gives a "
        "majority share")
    say("  Q4 %s" % ("hold" if q4 else "REFUTED"))

    # ------------------------------------------- not pre-registered
    say()
    say("X1  and the crossing the axis question asks about")
    say("  (written after Q3; not pre-registered). {#rem:headaxis} "
        "found the two")
    say("  spreads moving in opposite directions at near-equal rates "
        "and published")
    say("  no crossing. On each bounded basis the two destinations "
        "are these:")
    say("  basis        |I| goes to   T goes to    gap at the "
        "destination")
    LI, sI = dests["|I| spread"]
    LT, sT = dests["T spread"]
    for (bn, _g), li, lt, s1, s2 in zip(BASES, LI, LT, sI, sT):
        say("    %-12s %+-13.6f %+-12.6f %+.6f  (+- %.6f)"
            % (bn, li, lt, li - lt,
               2.0 * math.sqrt(s1 * s1 + s2 * s2)))
    gaps = [li - lt for li, lt in zip(LI, LT)]
    resolved = [abs(g) > 2.0 * math.sqrt(a * a + b * b)
                for g, a, b in zip(gaps, sI, sT)]
    agree = all(g < 0 for g in gaps) or all(g > 0 for g in gaps)
    say("  a negative gap is a handover completed, a positive one is "
        "not.")
    say("  the bases agree on the sign: %s, and each gap clears its "
        "own two-sigma: %s"
        % ("yes" if agree else "no",
           ", ".join("yes" if r else "no" for r in resolved)))
    say("  the gaps themselves run %+.6f to %+.6f, a factor of %.2f, "
        "so the" % (min(gaps), max(gaps),
                    max(abs(g) for g in gaps)
                    / min(abs(g) for g in gaps)))
    say("  direction of the handover is agreed and its size is not -- "
        "which is Q3's")
    say("  verdict in the units the question was asked in.")
    say("TSTAT slope_audit_power_reach %.2f"
        % abs(prepared[0][3] / prepared[0][6]))
    say("SPREAD slope_audit_power_reach %.4f"
        % float(prepared[0][1].max() - prepared[0][1].min()))
    say("SENSITIVITY powerreach %.6f %.6f"
        % (min(r[5] for r in prepared), max(r[5] for r in prepared)))

    say()
    say("=" * 70)
    say("Q1 %s  Q2 %s  Q3 %s  Q4 %s"
        % tuple("hold" if v else "REFUTED"
                for v in (q1, q2, q3, q4)))

    head = [
        "STATISTIC: for the head's one-sign fraction read from",
        "           results/audit_split_reach.txt and the |I| and T",
        "           decile spreads read from",
        "           results/audit_headaxis_reach.txt, the",
        "           least-squares slope in log N against the",
        "           published one and against the bound each table's",
        "           printing forces; the power of this repository's",
        "           shape discriminator on each, by synthetic draws",
        "           from the fitted linear shape at the observed",
        "           abscissae and residual r.m.s.; and the fitted",
        "           asymptote L of y = L + b g(log N) with its",
        "           standard error on each of the bases 1/x,",
        "           1/sqrt(x) and 1/x^2, with the spread of L across",
        "           bases against the largest within-basis error.",
        "NULL: the null IS the measurement for the power, as in",
        "      {#rem:shapepower}: every trial draws from a shape known",
        "      in advance. For the asymptotes no null applies -- they",
        "      are interval estimates, and the interval states what",
        "      the data exclude. The draws use the fixed SEED below.",
        "FIELD: no arithmetic is computed here. Both series are read",
        "       from files whose own field is the on-field",
        "       N = 2^a 5^b with a, b >= 1 in [2e5, 1.024e8], one",
        "       coprimality class; the one-sign fraction is the",
        "       majority sign share of the top tenth of k by |a_k|,",
        "       and the spreads are top-minus-bottom decile",
        "       differences of the negative share cut on |I| and on T",
        "       as code/audit_head_sign.py defines them.",
        "",
    ]
    io.open(OUT, "w", encoding="utf-8", newline="\n").write(
        "\n".join(head + lines) + "\n")
    print("\nwrote %s" % OUT)
    if not q1:
        raise SystemExit(1)
    return 0


if __name__ == "__main__":
    sys.exit(main())
