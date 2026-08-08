# -*- coding: utf-8 -*-
r"""
What shape is the ladder, and does the shape decide the answer?

WHAT IS AT STAKE

Remark {#rem:primorialrung10} observed the primorial ladder's exponent
cross 1/2 at N = 10^7.4879, and then said where the same line puts the
theta' = 0.56 these papers actually need: 10^11.2680. That second
number is not a measurement. It is a linear fit in log N extrapolated
four decades, and nothing has ever justified the shape.

The shape matters more than the fit. Balancing
sum_{k<K}(log k) c_R sqrt(N/k) against a budget of order N/log N gives
K of order N/log^4 N, whose exponent is 1 - c log log N / log N -- a
curve that rises to one but with a derivative that falls. A saturating
form a + b/log N rises to a and stops. A line rises forever. All three
can look alike over a factor 1024 in N and say entirely different
things about 0.56: one reaches it, one reaches it much later, one
never does.

This costs nothing to test. The eleven exponents are already
published; what is missing is the comparison of shapes on them.

BACKS: Remark {#rem:laddershape} in paper/theorem_A.md.

PRE-REGISTERED PREDICTIONS (written before this script was run)

  H1  The control: the line refitted here reproduces the published
      slope and r.m.s. residual of
      results/audit_primorial_rung10.txt to within 1e-4.
  H2  The line is not the best shape and not distinguishably so: at
      least one alternative fits with a residual scatter no larger
      than the line's.
  H3  And the shapes disagree about theta' = 0.56 by more than one
      decade in log10 N -- or one of them never reaches it at all.
      That is the one that matters: it would withdraw the 10^11.2680
      as a statement about anything.
  H4  But they agree where the data are: every shape that fits within
      the line's scatter puts the crossing of 1/2 inside the bracket
      results/audit_primorial_rung10.txt publishes.

REFUTATION RULE (fixed before the run)

  H1  REFUTED at 1e-4.
  H2  REFUTED if every alternative scatters more than the line, which
      would mean the linear shape is not merely convenient but
      preferred by the data.
  H3  REFUTED if all shapes agree on 0.56 to within one decade, which
      would make the extrapolation robust and the 10^11.2680 worth
      quoting.
  H4  REFUTED if any adequately fitting shape puts 1/2 outside the
      published bracket.

  All four gate.

  NO NULL IS RUN and none applies. Several deterministic shapes are
  least-squares fitted to eleven measured numbers and compared by
  residual; there is no background to detect against. The coin arms
  for these exponents were run in lab_primorial_ladder.py and
  lab_primorial_share.py, which established that both the rise and its
  scatter are facts about magnitudes and not about mu.
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
OUT = os.path.join(ROOT, "results", "audit_ladder_shape.txt")

HALF = 0.5
THETA = 0.56
UMAX = 400.0            # log N ceiling for the crossing search


def read_rungs():
    """the eleven exponents, the line and its bracket -- read"""
    p = os.path.join(ROOT, "results", "audit_primorial_rung10.txt")
    src = io.open(p, encoding="utf-8").read()
    i = src.index("N            log10 N   exponent   fitted     "
                  "residual")
    ex = {}
    for ln in src[i:].splitlines()[1:]:
        f = ln.split()
        if len(f) < 5 or not f[0].isdigit():
            break
        ex[int(f[0])] = float(f[2])
    m = re.search(r"slope %s([\d.]+), r\.m\.s\. residual ([\d.]+)"
                  % re.escape("+"), src)
    sl, rms = float(m.group(1)), float(m.group(2))
    b = re.search(r"BRACKET log10_N_primorial_reaches_half_v3 "
                  r"([\d.]+) ([\d.]+) ([\d.]+)", src)
    return ex, sl, rms, (float(b.group(1)), float(b.group(2)),
                         float(b.group(3)))


def main():
    lines = []

    def say(s=""):
        print(s)
        lines.append(s)

    ex, pubsl, pubrms, (bpt, blo, bhi) = read_rungs()
    ns = sorted(ex)
    say("read %d rungs, the published slope %+.6f, r.m.s. %.4f and"
        % (len(ns), pubsl, pubrms))
    say("  the bracket 10^%.4f [10^%.4f, 10^%.4f] from"
        % (bpt, blo, bhi))
    say("  results/audit_primorial_rung10.txt")

    u = np.array([math.log(n) for n in ns])
    y = np.array([ex[n] for n in ns])
    lu = np.log(u)

    # each shape: a name, the design matrix columns, and a function
    # giving the fitted exponent at a given log N
    def design(name, u_):
        if name == "line":
            return np.vstack([np.ones_like(u_), u_]).T
        if name == "saturating":
            return np.vstack([np.ones_like(u_), 1.0 / u_]).T
        if name == "loglog":
            return np.vstack([np.ones_like(u_), np.log(u_)]).T
        if name == "heuristic2":
            return np.vstack([np.ones_like(u_),
                              np.log(u_) / u_]).T
        raise KeyError(name)

    SHAPES = ["line", "saturating", "loglog", "heuristic2"]
    LABEL = {"line": "a + b log N",
             "saturating": "a + b / log N",
             "loglog": "a + b log log N",
             "heuristic2": "a + b log log N / log N"}

    fits = {}
    for nm in SHAPES:
        X = design(nm, u)
        c, *_ = np.linalg.lstsq(X, y, rcond=None)
        r = y - X @ c
        fits[nm] = (c, float(np.sqrt((r ** 2).mean())), r)

    # the one-parameter heuristic shape, 1 - c loglog/log
    z = lu / u
    c1 = float(((1.0 - y) * z).sum() / (z * z).sum())
    r1 = y - (1.0 - c1 * z)
    fits["heuristic1"] = (np.array([c1]),
                          float(np.sqrt((r1 ** 2).mean())), r1)
    SHAPES.append("heuristic1")
    LABEL["heuristic1"] = "1 - c log log N / log N"

    def value(nm, uu):
        c = fits[nm][0]
        if nm == "heuristic1":
            return 1.0 - c[0] * math.log(uu) / uu
        X = design(nm, np.array([uu]))
        return float(X @ c)

    def crossing(nm, target):
        lo, hi = u[-1], UMAX
        if value(nm, hi) < target:
            return None
        lo2 = 2.0
        for _ in range(300):
            mid = 0.5 * (lo2 + hi)
            if value(nm, mid) < target:
                lo2 = mid
            else:
                hi = mid
        return 0.5 * (lo2 + hi) / math.log(10.0)

    # ------------------------------------------------------------- H1
    say()
    c, rms, _ = fits["line"]
    d1 = abs(c[1] - pubsl)
    d2 = abs(rms - pubrms)
    h1 = d1 < 1e-4 and d2 < 1e-4
    say("H1  the control: the line refitted here")
    say("  slope %+.6f against the published %+.6f, diff %.6f"
        % (c[1], pubsl, d1))
    say("  r.m.s. %.4f against the published %.4f, diff %.6f"
        % (rms, pubrms, d2))
    say("  H1 %s" % ("hold" if h1 else "REFUTED"))

    # ------------------------------------------------------------- H2
    say()
    say("H2  how well each shape fits")
    say("  shape                        r.m.s. residual   vs the line")
    base = fits["line"][1]
    h2 = False
    for nm in SHAPES:
        rr = fits[nm][1]
        if nm != "line" and rr <= base:
            h2 = True
        say("  %-28s %-17.5f %s"
            % (LABEL[nm], rr,
               "-" if nm == "line" else "%+.1f%%"
               % (100.0 * (rr / base - 1.0))))
    say("  H2 %s" % ("hold" if h2 else "REFUTED"))

    # ---------------------------------------------------------- H3/H4
    say()
    say("H3/H4  where each shape puts the two levels")
    say("  shape                        log10 N at 0.50   at 0.56")
    ok56, ok50 = [], []
    for nm in SHAPES:
        a = crossing(nm, HALF)
        b = crossing(nm, THETA)
        adequate = fits[nm][1] <= base * 1.0000001 or nm == "line"
        if adequate:
            if a is not None:
                ok50.append(a)
            ok56.append(b)
        say("  %-28s %-17s %s"
            % (LABEL[nm],
               "%.4f" % a if a is not None else "never",
               "%.4f" % b if b is not None else "never"))
    fin = [v for v in ok56 if v is not None]
    never = any(v is None for v in ok56)
    spread = (max(fin) - min(fin)) if len(fin) > 1 else 0.0
    h3 = never or spread > 1.0
    say("  among the shapes that fit at least as well as the line:")
    say("    0.56 reached at %s%s"
        % (", ".join("%.4f" % v for v in fin) or "none",
           "; and NEVER by at least one" if never else ""))
    say("    spread %.4f decades" % spread)
    say("  H3 %s   (floor 1 decade, or one shape never reaching it)"
        % ("hold" if h3 else "REFUTED"))
    h4 = all(blo <= v <= bhi for v in ok50) if ok50 else False
    say("  1/2 reached at %s, published bracket [%.4f, %.4f]"
        % (", ".join("%.4f" % v for v in ok50), blo, bhi))
    say("  H4 %s" % ("hold" if h4 else "REFUTED"))

    say()
    say("  the arithmetic and the budget are those of the ladder,")
    say("  declared by the file this reads; nothing new is computed")
    say("  here beyond least squares on eleven published numbers.")
    say("RADICALS 1")

    say()
    say("  DIAGNOSTIC (post hoc). The residuals, shape by shape, so")
    say("  that 'fits as well' can be seen and not taken on trust:")
    say("  N            " + "  ".join("%-9s" % nm for nm in SHAPES))
    for i, n in enumerate(ns):
        say("  %-12d %s"
            % (n, "  ".join("%+9.4f" % fits[nm][2][i]
                            for nm in SHAPES)))
    say()
    say("  H2 and H3 asked for a shape that fits AT LEAST AS WELL as")
    say("  the line, and none does; the line is the best of the five.")
    dof = len(ns) - 2
    say("  That criterion is too strict and the numbers say why. With")
    say("  %d points and 2 parameters the r.m.s. residual is itself"
        % len(ns))
    say("  estimated from %d degrees of freedom, so its own standard"
        % dof)
    say("  error is r.m.s./sqrt(2 df):")
    se = base / math.sqrt(2.0 * dof)
    say("    best r.m.s. %.5f, its standard error %.5f, i.e. %.1f%%"
        % (base, se, 100.0 * se / base))
    say("  A shape is not excluded unless it sits more than that above")
    say("  the best. Which do:")
    say("  shape                        r.m.s.    (r-best)/s.e.  "
        "excluded?")
    keep = []
    for nm in SHAPES:
        rr = fits[nm][1]
        t = (rr - base) / se
        exc = t > 1.0
        if not exc:
            keep.append(nm)
        say("  %-28s %-9.5f %-14.2f %s"
            % (LABEL[nm], rr, t, "yes" if exc else "no"))
    say("  %d of %d shapes survive: %s"
        % (len(keep), len(SHAPES),
           ", ".join(LABEL[n] for n in keep)))
    say("  and where they put the two levels:")
    say("  shape                        0.50      0.56")
    k50, k56 = [], []
    for nm in keep:
        a = crossing(nm, HALF)
        b = crossing(nm, THETA)
        if a is not None:
            k50.append(a)
        k56.append(b)
        say("  %-28s %-9s %s"
            % (LABEL[nm],
               "%.4f" % a if a is not None else "never",
               "%.4f" % b if b is not None else "never"))
    fk = [v for v in k56 if v is not None]
    if len(fk) > 1:
        say("  0.56 spans %.4f decades among the surviving shapes;"
            % (max(fk) - min(fk)))
    if k50:
        say("  0.50 spans %.4f decades among them, against a published"
            % (max(k50) - min(k50)))
        say("  bracket %.4f decades wide" % (bhi - blo))

    say()
    say("  and what each shape does in the limit:")
    for nm in SHAPES:
        big = value(nm, 1e6)
        say("  %-28s exponent at log N = 1e6: %s"
            % (LABEL[nm],
               "%.4f" % big if abs(big) < 100 else "%.3e" % big))

    say()
    say("=" * 70)
    ok = h1 and h2 and h3 and h4
    say("the shapes agree where the data are and disagree entirely "
        "about theta'" if ok else "REFUTED")

    head = [
        "STATISTIC: least-squares fits of five shapes to the eleven",
        "           published exponents of the primorial ladder --",
        "           a + b log N, a + b/log N, a + b log log N,",
        "           a + b log log N / log N, and the one-parameter",
        "           1 - c log log N / log N -- their r.m.s. residuals,",
        "           and the log10 N at which each reaches 0.50 and 0.56.",
        "NULL: none is run and none applies. Deterministic shapes are",
        "      fitted to eleven measured numbers and compared by",
        "      residual; there is no background to detect against. The",
        "      coin arms for these exponents were run in",
        "      lab_primorial_ladder.py and lab_primorial_share.py.",
        "FIELD: the eleven exponents of N = 30030 * 2^j, j = 0..10, read",
        "       from results/audit_primorial_rung10.txt, together with",
        "       the line and bracket that file publishes. No new",
        "       arithmetic is done here; the k-range, sieve weight and",
        "       budget are those declared there.",
        "",
    ]
    io.open(OUT, "w", encoding="utf-8", newline="\n").write(
        "\n".join(head + lines) + "\n")
    print("\nwrote %s" % OUT)
    if not ok:
        raise SystemExit(1)
    return 0


if __name__ == "__main__":
    sys.exit(main())
