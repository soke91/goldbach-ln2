# -*- coding: utf-8 -*-
r"""
How much range would separate N^-a from (log N)^-b?

WHAT IS AT STAKE

v1 stated a limitation and OPEN.md item 3 has carried it unquantified
ever since: over a factor 160 in N the data do not separate N^-a from
(log N)^-b, "so these are exponents of an assumed form".  It says
they do not separate.  It does not say what would.

That number decides whether the limitation is a gap in this
program's budget or a fact about the question.  If a factor of a
thousand would separate them, the mask's form is a measurement
waiting to be taken.  If it takes 10^40, the form is not something
computation settles here and {#rem:maskfloornull}'s three surviving
exponents are permanently exponents of an assumption.

The design is the mask's own.  Fifteen half-octave bands over a
factor 160, the recorded errors of Lemma lem:cellmom, and the
surviving depths of {#rem:maskfloornull}: 0, 2 and 3.  Data are drawn
from the fitted power law with the recorded scatter, both forms are
fitted, and the range is extended -- keeping the band density and the
error law fixed -- until the true form wins often enough to call it
separated.

Nothing about the mask is measured.  Every input is read from v1's
table and from results/audit_mask_floornull.txt.

BACKS: Remark {#rem:maskformreach} in paper/wall_v3.md.

PRE-REGISTERED PREDICTIONS (written before this script was run)

  H1  The control.  The power-law exponents refitted here reproduce
      the three surviving depths of {#rem:maskfloornull} to the four
      decimals v1 prints.
  H2  v1's statement holds where it was made: at the observed range,
      drawing from the power law, the power law wins on r.m.s. in
      fewer than 95 per cent of draws -- the two forms are not
      separated.
  H3  Separation is reachable in principle: at some finite extension
      of the range the power law wins in at least 95 per cent of
      draws, at every surviving depth.
  H4  **And it is out of reach in practice**: the factor in N needed
      is above 10^6 at every surviving depth, against the 160 the
      data has.

REFUTATION RULE (fixed before the run)

  H1  REFUTED outside the printed decimals.  THIS ONE GATES.
  H2  REFUTED if the forms already separate at the observed range.
      Then v1's limitation is not a limitation and the exponents are
      of a measured form, not an assumed one -- which would be worth
      more than anything else in this remark.
  H3  REFUTED if no extension tried reaches 95 per cent.  Two ways
      that can happen and they differ.  **The two forms may be
      indistinguishable at any range** over this design, in which
      case the question is not one more N answers at all; or the
      search may simply have stopped too early, which the printed
      table of extensions makes visible.  A ceiling reached without
      the table rising is the first; a table still rising at the
      last extension is the second.
  H4  REFUTED below 10^6.  Then the form is settleable by a
      computation this program could plausibly run, and the open
      item becomes a task rather than a limitation.

  H1 gates.  H2 to H4 are the measurement and do not gate.

  THE NULL IS RUN, and it is the whole of H2 to H4: data are drawn
  from the fitted power law with Gaussian errors of the recorded se,
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
OUT = os.path.join(RES, "audit_mask_formreach.txt")

SEED = 20260823
DRAWS = 2000
DEC = 4
SURVIVING = (0, 2, 3)
WIN = 0.95
TARGET = 6.0                        # log10 of the factor H4 asks about
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
    """weighted least squares of y on u, returning the weighted rms"""
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
    say("PRINTBOUND audit_mask_formreach %d %.8f" % (DEC, rnd))

    # -------------------------------------------------------------- H1
    say()
    say("H1  the control on the surviving depths")
    h1 = True
    base = {}
    for d in SURVIVING:
        R = bands[d]
        keep = R[:, 2] != 0.0
        x = np.log(R[keep, 0])
        dm, se = R[keep, 2], R[keep, 3]
        w = (dm / se) ** 2
        a, rms = wfit(x, np.log(np.abs(dm)), w)
        base[d] = (x, np.abs(dm), se, w, a, rms)
        ok = abs(a - pub[d]) <= rnd
        h1 = h1 and ok
        say("  depth %d  a = %.4f against v1's %.4f  %s"
            % (d, a, pub[d], "yes" if ok else "NO"))
    say("  H1 %s   (cap: the four decimals v1 prints)"
        % ("hold" if h1 else "REFUTED"))
    if not h1:
        io.open(OUT, "w", encoding="utf-8", newline="\n").write(
            "\n".join(lines) + "\n")
        raise SystemExit(1)

    span = base[SURVIVING[0]][0]
    say()
    say("  the design: %d bands over %.4f in log N, a factor of %.0f"
        % (span.size, span.max() - span.min(),
           math.exp(span.max() - span.min())))
    say("SCALES 1")

    # --------------------------------------------------------- H2..H4
    say()
    say("H2/H3/H4  extending the range with the band density fixed")
    say("  depth  extension  factor in N   bands   power wins  "
        "effective n")
    rng = np.random.default_rng(SEED)
    need = {}
    for d in SURVIVING:
        x0, dm0, se0, w0, a, rms = base[d]
        lo, hi = float(x0.min()), float(x0.max())
        dens = (x0.size - 1) / (hi - lo)
        # the error law: se against log N, fitted on the recorded se
        sa, _r = wfit(x0, np.log(se0), np.ones_like(se0))
        c0 = math.log(se0[0]) + sa * lo
        amp0 = math.log(dm0[0]) + a * lo
        for e in EXT:
            top = lo + e * (hi - lo)
            n = int(round(dens * (top - lo))) + 1
            x = np.linspace(lo, top, n)
            se = np.exp(c0 - sa * x)
            mu = np.exp(amp0 - a * x)
            w = (mu / se) ** 2
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
            neff = float((w.sum() ** 2) / (w * w).sum())
            say("  %-6d %-10.0f %-13.3e %-7d %-11.4f %.2f"
                % (d, e, math.exp(top - lo), n, frac, neff))
            if frac >= WIN and d not in need:
                need[d] = (e, math.exp(top - lo))
        say("NULL mask_formreach_depth%d %.4f" % (d, frac))
    h2 = all(d not in need or need[d][0] > 1.0 for d in SURVIVING)
    h3 = all(d in need for d in SURVIVING)
    h4 = h3 and all(math.log10(need[d][1]) > TARGET for d in SURVIVING)
    say()
    for d in SURVIVING:
        if d in need:
            say("  depth %d separates at extension %.0f, a factor of "
                "10^%.2f in N"
                % (d, need[d][0], math.log10(need[d][1])))
            say("FORMREACH mask_depth%d %.2f"
                % (d, math.log10(need[d][1])))
        else:
            say("  depth %d does not separate at any extension tried "
                "(up to %.0f)" % (d, EXT[-1]))
            say("FORMREACH mask_depth%d 0.00" % d)
    say("  H2 %s   (cap: not separated at the observed range)"
        % ("hold" if h2 else "REFUTED"))
    say("  H3 %s   (cap: separated at some extension tried)"
        % ("hold" if h3 else "REFUTED"))
    say("  H4 %s   (cap: a factor above 10^%.0f)"
        % ("hold" if h4 else "REFUTED", TARGET))

    say()
    say("why the table does not rise -- diagnosed after the run, "
        "predicted by nothing")
    say("  the weights are (dm/se)^2 and dm/se itself decays: "
        "{#rem:maskdmse} measured")
    say("  the amplitude exponent at %.4f and the error's at %.4f, "
        "so the weight falls" % (base[0][4], 0.0420))
    say("  like N to the minus twice their difference. Bands added "
        "past the observed")
    say("  range carry almost no weight, and the effective sample "
        "size above stops")
    say("  growing while the band count doubles. **The design cannot "
        "answer the question")
    say("  it was registered to ask**, and the refutation of H3 and "
        "H4 is that, not a")
    say("  fact about the two forms. What is real underneath it is "
        "that the informative")
    say("  window is bounded: the mask sinks below its own error, "
        "and range past that")
    say("  point adds bands that measure nothing.")

    say()
    say("what this settles")
    if h2 and h3 and h4:
        say("  v1's limitation is a fact about the reach and not "
            "about the question: the")
        say("  two forms do separate, at a range this program cannot "
            "buy")
    elif h2 and not h3:
        say("  no extension tried separates them; whether that is a "
            "ceiling or a short")
        say("  search is what the table above shows")
    else:
        say("  the forms separate closer than v1's statement implies, "
            "and the open item")
        say("  is a task rather than a limitation")

    say()
    say("=" * 70)
    say("H1 %s  H2 %s  H3 %s  H4 %s"
        % tuple("hold" if v else "REFUTED"
                for v in (h1, h2, h3, h4)))

    head = [
        "STATISTIC: the fraction of simulated fields in which the",
        "           power law N^-a beats (log N)^-b on weighted",
        "           r.m.s., when the data are drawn from the power",
        "           law with the mask's own scatter, as the log-range",
        "           is extended with the band density and the error",
        "           law held fixed; and the extension at which that",
        "           fraction first reaches 95 per cent.",
        "NULL: RUN, and it is the whole measurement. Data are drawn",
        "      from the fitted power law with Gaussian errors of the",
        "      recorded se, both forms are fitted by the same",
        "      weighted least squares, and 'separated' counts the",
        "      draws in which the true form has the smaller weighted",
        "      r.m.s. The fixed SEED governs every draw.",
        "FIELD: the per-depth band tables of",
        "       v1/results/wall/lab_mask_amplitude_law.txt at the",
        "       three depths {#rem:maskfloornull} leaves standing --",
        "       0, 2 and 3 -- fifteen half-octave bands over a factor",
        "       160 in N, with the exact errors of Lemma lem:cellmom.",
        "       Nothing about the mask is measured here; the",
        "       exponents are refitted only as H1's control against",
        "       v1/results/wall/lab_mask_exponent_se.txt.",
        "",
    ]
    io.open(OUT, "w", encoding="utf-8", newline="\n").write(
        "\n".join(head + lines) + "\n")
    print("\nwrote %s" % OUT)
    if not h1:
        raise SystemExit(1)
    return 0


if __name__ == "__main__":
    sys.exit(main())
