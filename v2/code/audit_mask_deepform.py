# -*- coding: utf-8 -*-
r"""
The same form question, asked where the mask has not sunk.

WHAT IS AT STAKE

{#rem:maskformreach} tried to measure what range would separate
N^-a from (log N)^-b for the mask and could not, because at the
depths it used the weights collapse: the effective sample froze at
5.01 while the band count went to 1793.  It named the reason -- the
weight is (dm/se)^2 and falls like N to twice the difference of the
amplitude and error exponents -- and it named what it had not tried:
**a design that measures where the mask has not yet sunk below its
own error.**

That design is in the same table.  At depths 3, 4 and 5 every band
clears |z| = 2, running from 2.04 to 12.19, and the exponent
difference that kills the weight is small: {#rem:maskdmse} measured
0.0317 at depth 5 against 0.4659 and 0.2309 at the shallower ones.
So at depth 5 the weight falls like N^-0.063 -- almost flat -- and
extending the range should add bands that carry information rather
than bands that carry none.

If the forms separate there at a reachable range, the mask's form is
measurable and only the shallow depths are hopeless.  If they do not
separate even where the signal is strongest and the weights do not
collapse, then it is the two forms that are close over any range this
design can offer, and OPEN.md item 3's limitation is about the
question rather than about the reach.

Nothing about the mask is measured.  Every input is read from v1's
table.

BACKS: Remark {#rem:maskdeepform} in paper/wall_v3.md.

PRE-REGISTERED PREDICTIONS (written before this script was run)

  I1  The control.  The power-law exponents refitted here reproduce
      depths 3, 4 and 5 to the four decimals v1 prints.
  I2  The weights do not collapse: at depth 5 the effective sample
      size grows by more than a factor of four as the range is
      extended, against the 4.99 to 5.01 {#rem:maskformreach}
      recorded at depth 0.
  I3  v1's statement holds at these depths too: at the observed
      range the power law wins in fewer than 95 per cent of draws at
      all three.
  I4  **And here extending the range works**: at depth 5 the power
      law reaches 95 per cent at some extension tried.
  I5  But not cheaply: the factor in N at which it does is above
      10^6.

REFUTATION RULE (fixed before the run)

  I1  REFUTED outside the printed decimals.  THIS ONE GATES.
  I2  REFUTED if the effective sample grows by four or less.  Then
      this design collapses like the last one and I3 to I5 are
      uninterpretable in the same way -- which is the outcome that
      would say the collapse is not about the depth but about the
      weighting itself.
  I3  REFUTED if the forms separate already at the observed range at
      any of the three.  Then v1's limitation was a statement about
      the shallow depths only and the mask's form is measured, not
      assumed, where the signal is strong.
  I4  REFUTED if depth 5 never reaches 95 per cent.  With I2 holding
      that cannot be blamed on the weights, so it would say the two
      forms are genuinely close over this design at any range --
      **and "does not reach 95 per cent" must be read against the
      printed table**: a fraction still rising at the last extension
      means the search stopped early, a fraction flat means it did
      not.
  I5  REFUTED below 10^6.  Then the mask's form is settleable by a
      computation of a size this program has already run on other
      axes, and item 3 becomes a task.

  I1 gates.  I2 to I5 are the measurement and do not gate.

  THE NULL IS RUN, and it is the whole of I2 to I5: data are drawn
  from the fitted power law with Gaussian errors of the modelled se,
  both forms are fitted by the same weighted least squares, and the
  fraction of draws in which the power law has the smaller weighted
  r.m.s. is what "separated" counts.
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
V1 = os.path.abspath(os.path.join(ROOT, "..", "v1", "results", "wall"))
OUT = os.path.join(RES, "audit_mask_deepform.txt")

SEED = 20260823
DRAWS = 2000
DEC = 4
DEEP = (3, 4, 5)
WIN = 0.95
TARGET = 6.0
EXT = (1.0, 2.0, 4.0, 8.0, 16.0, 32.0, 64.0, 128.0)

ROW = (r"^  depth (\d)\s+\([^)]*\)\s+(\d+) points\s*\n[^\n]*N\s+n_c"
       r"[^\n]*\n((?:^\s+[\d.eE+-]+\s+\d+\s+[-+]?[\d.]+\s+[\d.]+\s+"
       r"[-+]?[\d.]+\s*\n)+)")


def read_bands():
    src = io.open(os.path.join(V1, "lab_mask_amplitude_law.txt"),
                  encoding="utf-8").read()
    return {int(d): np.array([[float(t) for t in ln.split()]
                              for ln in body.strip().splitlines()])
            for d, _n, body in re.findall(ROW, src, re.M)}


def read_published():
    src = io.open(os.path.join(V1, "lab_mask_exponent_se.txt"),
                  encoding="utf-8").read()
    return {int(m.group(1)): float(m.group(3))
            for m in re.finditer(
                r"^\s+(\d)\s+(\d+)\s+([\d.]+)\s+([\d.]+)\s+"
                r"([\d.]+)\s+([\d.]+)\s*$", src, re.M)}


def wfit(u, y, w):
    A = np.column_stack([np.ones_like(u), u])
    W = np.diag(w)
    c = np.linalg.inv(A.T.dot(W).dot(A)).dot(A.T).dot(W).dot(y)
    r = y - A.dot(c)
    return -float(c[1]), float(math.sqrt((w * r * r).sum() / w.sum()))


def main():
    lines = []

    def say(t=""):
        print(t)
        sys.stdout.flush()
        lines.append(t)

    bands = read_bands()
    pub = read_published()
    rnd = 0.5 * 10.0 ** (-DEC)
    say("read the depth tables from "
        "v1/results/wall/lab_mask_amplitude_law.txt; nothing about "
        "the mask")
    say("  is measured here")
    say("SEED %d" % SEED)
    say("DRAWS %d" % DRAWS)
    say("PRINTBOUND audit_mask_deepform %d %.8f" % (DEC, rnd))

    # -------------------------------------------------------------- I1
    say()
    say("I1  the control on the depths where every band clears")
    i1 = True
    base = {}
    for d in DEEP:
        R = bands[d]
        keep = R[:, 2] != 0.0
        x = np.log(R[keep, 0])
        dm, se = np.abs(R[keep, 2]), R[keep, 3]
        w = (dm / se) ** 2
        a, _r = wfit(x, np.log(dm), w)
        sa, _r2 = wfit(x, np.log(se), np.ones_like(se))
        base[d] = (x, dm, se, a, sa)
        ok = abs(a - pub[d]) <= rnd
        i1 = i1 and ok
        z = np.abs(R[keep, 4])
        say("  depth %d  a = %.4f against v1's %.4f; se exponent "
            "%.4f; |z| runs %.2f to %.2f  %s"
            % (d, a, pub[d], sa, z.min(), z.max(),
               "yes" if ok else "NO"))
        say("  the weight falls like N to minus %.4f" % (2 * (a - sa)))
    say("  I1 %s   (cap: the four decimals v1 prints)"
        % ("hold" if i1 else "REFUTED"))
    if not i1:
        io.open(OUT, "w", encoding="utf-8", newline="\n").write(
            "\n".join(lines) + "\n")
        raise SystemExit(1)

    # --------------------------------------------------------- I2..I5
    say()
    say("I2/I3/I4/I5  extending the range where the signal is strong")
    say("  depth  extension  factor in N   bands   power wins  "
        "effective n")
    rng = np.random.default_rng(SEED)
    need, neff0, neffL, first = {}, {}, {}, {}
    for d in DEEP:
        x0, dm0, se0, a, sa = base[d]
        lo, hi = float(x0.min()), float(x0.max())
        dens = (x0.size - 1) / (hi - lo)
        c0 = math.log(se0[0]) + sa * lo
        amp0 = math.log(dm0[0]) + a * lo
        for e in EXT:
            top = lo + e * (hi - lo)
            n = int(round(dens * (top - lo))) + 1
            x = np.linspace(lo, top, n)
            se = np.exp(c0 - sa * x)
            mu = np.exp(amp0 - a * x)
            w = (mu / se) ** 2
            neff = float((w.sum() ** 2) / (w * w).sum())
            wins = 0
            for _ in range(DRAWS):
                y = mu + rng.normal(0.0, se)
                ok = y > 0
                if ok.sum() < 4:
                    continue
                _a1, r1 = wfit(x[ok], np.log(y[ok]), w[ok])
                _b1, r2 = wfit(np.log(x[ok]), np.log(y[ok]), w[ok])
                if r1 < r2:
                    wins += 1
            frac = wins / float(DRAWS)
            if e == EXT[0]:
                neff0[d], first[d] = neff, frac
            neffL[d] = neff
            say("  %-6d %-10.0f %-13.3e %-7d %-11.4f %.2f"
                % (d, e, math.exp(top - lo), n, frac, neff))
            if frac >= WIN and d not in need:
                need[d] = math.exp(top - lo)
        say("NULL mask_deepform_depth%d %.4f" % (d, frac))

    say()
    i2 = neffL[5] / neff0[5] > 4.0
    say("  at depth 5 the effective sample goes %.2f to %.2f, a "
        "factor of %.1f" % (neff0[5], neffL[5], neffL[5] / neff0[5]))
    say("  I2 %s   (cap: a factor of four)"
        % ("hold" if i2 else "REFUTED"))
    i3 = all(first[d] < WIN for d in DEEP)
    say("  at the observed range the power law wins %s"
        % ", ".join("%.4f" % first[d] for d in DEEP))
    say("  I3 %s   (cap: none reaching 95 per cent)"
        % ("hold" if i3 else "REFUTED"))
    i4 = 5 in need
    if i4:
        say("FORMREACH mask_deep_depth5 %.2f" % math.log10(need[5]))
        say("  depth 5 separates at a factor of 10^%.2f in N"
            % math.log10(need[5]))
    else:
        say("FORMREACH mask_deep_depth5 0.00")
        say("  depth 5 does not separate at any extension tried")
    say("  I4 %s   (cap: 95 per cent at some extension)"
        % ("hold" if i4 else "REFUTED"))
    i5 = i4 and math.log10(need[5]) > TARGET
    say("  I5 %s   (cap: a factor above 10^%.0f)"
        % ("hold" if i5 else "REFUTED", TARGET))

    say()
    say("what this settles")
    if i2 and i4:
        say("  where the mask has not sunk the weights do not "
            "collapse and the forms do")
        say("  separate, so the limitation is the shallow depths' "
            "and not the question's")
    elif i2 and not i4:
        say("  the weights do not collapse and the forms still do "
            "not separate, so the two")
        say("  are close over this design at any range and item 3's "
            "limitation is the")
        say("  question's, not the reach's")
    else:
        say("  the weights collapse here too, so the collapse is the "
            "weighting's and not")
        say("  the depth's, and nothing below I2 can be read")

    say()
    say("=" * 70)
    say("I1 %s  I2 %s  I3 %s  I4 %s  I5 %s"
        % tuple("hold" if v else "REFUTED"
                for v in (i1, i2, i3, i4, i5)))

    head = [
        "STATISTIC: the fraction of simulated fields in which N^-a",
        "           beats (log N)^-b on weighted r.m.s. at the three",
        "           depths where every band clears |z| = 2, as the",
        "           log-range is extended with the band density and",
        "           the error law held fixed; the effective sample",
        "           size at each extension; and the factor in N at",
        "           which the fraction first reaches 95 per cent.",
        "NULL: RUN, and it is the whole measurement. Data are drawn",
        "      from the fitted power law with Gaussian errors of the",
        "      modelled se, both forms are fitted by the same",
        "      weighted least squares, and 'separated' counts the",
        "      draws in which the true form has the smaller weighted",
        "      r.m.s. The fixed SEED governs every draw.",
        "FIELD: the per-depth band tables of",
        "       v1/results/wall/lab_mask_amplitude_law.txt at depths",
        "       3, 4 and 5, where |z| runs from 2.04 to 12.19 and no",
        "       band fails to clear -- the region",
        "       {#rem:maskformreach} named as untried. Fifteen",
        "       half-octave bands over a factor 160 in N, with the",
        "       exact errors of Lemma lem:cellmom. Nothing about the",
        "       mask is measured; the exponents are refitted only as",
        "       I1's control against",
        "       v1/results/wall/lab_mask_exponent_se.txt.",
        "",
    ]
    io.open(OUT, "w", encoding="utf-8", newline="\n").write(
        "\n".join(head + lines) + "\n")
    print("\nwrote %s" % OUT)
    if not i1:
        raise SystemExit(1)
    return 0


if __name__ == "__main__":
    sys.exit(main())
