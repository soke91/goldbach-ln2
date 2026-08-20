# -*- coding: utf-8 -*-
r"""
Is 10^7.72 the price of the form, or the price against one rival?

WHAT IS AT STAKE

{#rem:maskdeepform} priced the mask's form question: at depth 5,
where the weights do not collapse, N^-a separates from (log N)^-b at
a factor of 10^7.72 in N.  It closed with a caution that is really a
question -- **the separation measured is between two named forms and
says nothing about a third.**

That matters for how the number should be read.  If every reasonable
two-parameter rival separates at about the same range, 10^7.72 is the
price of establishing the form.  If the rivals differ by orders of
magnitude, it is the price against (log N)^-b in particular, and the
honest figure for "the form" is whatever the hardest rival costs.

Three rivals are put against the power law, each with two free
parameters so the comparison of weighted r.m.s. is fair, and each a
shape this literature actually writes:

    power     log y = A - a x            (the fitted form, x = log N)
    log-power log y = A - b log x        ({#rem:maskdeepform}'s rival)
    stretched log y = A - c sqrt(x)      (a stretched exponential in N)
    root-log  log y = A - d x / log x    (the shape a sieve bound takes)

Nothing about the mask is measured.  Every input is read from v1's
table, and the design -- band density, error law, extensions -- is
{#rem:maskdeepform}'s unchanged.

BACKS: Remark {#rem:maskrivals} in paper/wall_v3.md.

PRE-REGISTERED PREDICTIONS (written before this script was run)

  J1  The control.  Depth 5's exponent reproduces v1's to four
      decimals, and the reach against (log N)^-b reproduces
      {#rem:maskdeepform}'s 10^7.72 to the two decimals it prints.
  J2  At the observed range no rival is separated: the power law
      wins in fewer than 95 per cent of draws against all three.
  J3  **The reaches differ by at least a factor of ten in N**, so
      10^7.72 is a price against one rival and not the price of the
      form.
  J4  And the hardest rival costs more than (log N)^-b did: its
      reach exceeds 10^7.72.

REFUTATION RULE (fixed before the run)

  J1  REFUTED outside the printed decimals.  THIS ONE GATES.
  J2  REFUTED if any rival is already separated at the observed
      range.  Then the mask's form is measured against that rival
      today and only the others are open.
  J3  REFUTED if all three reaches fall within a factor of ten.
      That is the outcome that would let the number be quoted as the
      price of the form rather than of a matchup -- the stronger
      reading, and not the one predicted.
  J4  REFUTED if every rival separates at or below 10^7.72, in which
      case {#rem:maskdeepform}'s figure is already the worst case.
      **A rival that separates at no extension tried is not evidence
      for J4 and must be read from the table**: a fraction still
      rising at the last extension is a short search, a flat one is
      a rival this design cannot separate at any range, and the
      second would make "the price of the form" undefined rather
      than large.

  J1 gates.  J2 to J4 are the measurement and do not gate.

  THE NULL IS RUN, and it is the whole measurement: data are drawn
  from the fitted power law with Gaussian errors of the modelled se,
  every form is fitted by the same weighted least squares, and the
  fraction of draws in which the power law has the smallest weighted
  r.m.s. against a given rival is what "separated" counts.
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
OUT = os.path.join(RES, "audit_mask_rivals.txt")

SEED = 20260823
DRAWS = 2000
DEC = 4
DEPTH = 5
WIN = 0.95
EXT = (1.0, 2.0, 4.0, 8.0, 16.0, 32.0, 64.0, 128.0)

RIVALS = (("log-power", lambda x: np.log(x)),
          ("stretched", lambda x: np.sqrt(x)),
          ("root-log", lambda x: x / np.log(x)))

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


def read_reach():
    src = io.open(os.path.join(RES, "audit_mask_deepform.txt"),
                  encoding="utf-8").read()
    m = re.search(r"^FORMREACH mask_deep_depth5 ([\d.]+)\s*$",
                  src, re.M)
    if not m:
        raise SystemExit("no FORMREACH mask_deep_depth5 marker")
    return float(m.group(1))


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
    prev = read_reach()
    rnd = 0.5 * 10.0 ** (-DEC)
    say("read the depth-%d table from "
        "v1/results/wall/lab_mask_amplitude_law.txt and the reach"
        % DEPTH)
    say("READ audit_mask_deepform.txt FORMREACH mask_deep_depth5 "
        "%.2f" % prev)
    say("SEED %d" % SEED)
    say("DRAWS %d" % DRAWS)
    say("PRINTBOUND audit_mask_rivals %d %.8f" % (DEC, rnd))

    R = bands[DEPTH]
    keep = R[:, 2] != 0.0
    x0 = np.log(R[keep, 0])
    dm0, se0 = np.abs(R[keep, 2]), R[keep, 3]
    w0 = (dm0 / se0) ** 2
    a, _r = wfit(x0, np.log(dm0), w0)
    sa, _r2 = wfit(x0, np.log(se0), np.ones_like(se0))
    say()
    say("  depth %d: a = %.4f against v1's %.4f, se exponent %.4f"
        % (DEPTH, a, pub[DEPTH], sa))
    say("SCALES 1")

    lo, hi = float(x0.min()), float(x0.max())
    dens = (x0.size - 1) / (hi - lo)
    c0 = math.log(se0[0]) + sa * lo
    amp0 = math.log(dm0[0]) + a * lo
    rng = np.random.default_rng(SEED)

    say()
    say("  extension  factor in N   %s"
        % "  ".join("%-11s" % nm for nm, _f in RIVALS))
    need = {}
    firsts = {}
    for e in EXT:
        top = lo + e * (hi - lo)
        n = int(round(dens * (top - lo))) + 1
        x = np.linspace(lo, top, n)
        se = np.exp(c0 - sa * x)
        mu = np.exp(amp0 - a * x)
        w = (mu / se) ** 2
        wins = {nm: 0 for nm, _f in RIVALS}
        for _ in range(DRAWS):
            y = mu + rng.normal(0.0, se)
            ok = y > 0
            if ok.sum() < 4:
                continue
            ly = np.log(y[ok])
            _p, r0 = wfit(x[ok], ly, w[ok])
            for nm, f in RIVALS:
                _q, r1 = wfit(f(x[ok]), ly, w[ok])
                if r0 < r1:
                    wins[nm] += 1
        fr = {nm: wins[nm] / float(DRAWS) for nm, _f in RIVALS}
        if e == EXT[0]:
            firsts = dict(fr)
        say("  %-10.0f %-13.3e %s"
            % (e, math.exp(top - lo),
               "  ".join("%-11.4f" % fr[nm] for nm, _f in RIVALS)))
        for nm, _f in RIVALS:
            if fr[nm] >= WIN and nm not in need:
                need[nm] = math.exp(top - lo)

    # -------------------------------------------------------------- J1
    say()
    say("J1  the control")
    j1a = abs(a - pub[DEPTH]) <= rnd
    lp = math.log10(need["log-power"]) if "log-power" in need else 0.0
    j1b = abs(lp - prev) <= 0.005
    j1 = j1a and j1b
    say("  exponent %.4f against v1's %.4f; the (log N)^-b reach "
        "10^%.2f against 10^%.2f" % (a, pub[DEPTH], lp, prev))
    say("  J1 %s   (cap: the printed decimals)"
        % ("hold" if j1 else "REFUTED"))
    if not j1:
        io.open(OUT, "w", encoding="utf-8", newline="\n").write(
            "\n".join(lines) + "\n")
        raise SystemExit(1)

    # ---------------------------------------------------------- J2..J4
    say()
    say("J2  is anything separated at the observed range?")
    j2 = all(firsts[nm] < WIN for nm, _f in RIVALS)
    say("  %s" % ", ".join("%s %.4f" % (nm, firsts[nm])
                           for nm, _f in RIVALS))
    say("  J2 %s   (cap: none reaching 95 per cent)"
        % ("hold" if j2 else "REFUTED"))

    say()
    say("J3/J4  how far apart are the rivals?")
    got = []
    for nm, _f in RIVALS:
        if nm in need:
            v = math.log10(need[nm])
            got.append(v)
            say("RIVALREACH mask_%s %.2f" % (nm.replace("-", ""), v))
            say("  %-11s separates at 10^%.2f" % (nm, v))
        else:
            say("RIVALREACH mask_%s 0.00" % nm.replace("-", ""))
            say("  %-11s does not separate at any extension tried"
                % nm)
    unresolved = [nm for nm, _f in RIVALS if nm not in need]
    if got:
        spread = max(got) - min(got)
        say("SPREAD mask_rivals %.2f" % spread)
        say("  the separated rivals span %.2f in log10 of the factor"
            % spread)
    else:
        spread = 0.0
    j3 = spread >= 1.0 or bool(unresolved)
    j4 = (bool(unresolved)
          or (got and max(got) > prev + 0.005))
    say("  J3 %s   (cap: a factor of ten between reaches)"
        % ("hold" if j3 else "REFUTED"))
    say("  J4 %s   (cap: above 10^%.2f)"
        % ("hold" if j4 else "REFUTED", prev))
    if unresolved:
        say("  %d rival(s) unseparated at any extension: %s -- read "
            "the table above"
            % (len(unresolved), ", ".join(unresolved)))
        say("  a column still rising at the last extension is a "
            "short search; a flat one")
        say("  is a rival this design cannot separate at any range")

    say()
    say("what this settles")
    if j3:
        say("  the reaches are rival-specific, so 10^%.2f is the "
            "price against (log N)^-b" % prev)
        say("  and not the price of establishing the form; the "
            "honest figure is the hardest")
        say("  rival's, and it is %s"
            % ("undefined while a rival stays unseparated"
               if unresolved else "10^%.2f" % max(got)))
    else:
        say("  every rival costs about the same, so the number may "
            "be read as the price of")
        say("  the form rather than of a matchup")

    say()
    say("=" * 70)
    say("J1 %s  J2 %s  J3 %s  J4 %s"
        % tuple("hold" if v else "REFUTED"
                for v in (j1, j2, j3, j4)))

    head = [
        "STATISTIC: at depth 5, the fraction of simulated fields in",
        "           which N^-a beats each of three two-parameter",
        "           rivals on weighted r.m.s. as the log-range is",
        "           extended, and the factor in N at which each",
        "           fraction first reaches 95 per cent.",
        "NULL: RUN, and it is the whole measurement. Data are drawn",
        "      from the fitted power law with Gaussian errors of the",
        "      modelled se, every form is fitted by the same",
        "      weighted least squares on two parameters, and",
        "      'separated' counts the draws in which the true form",
        "      has the smaller weighted r.m.s. The fixed SEED",
        "      governs every draw.",
        "FIELD: the depth-5 band table of",
        "       v1/results/wall/lab_mask_amplitude_law.txt -- the",
        "       cell where {#rem:maskdeepform} found the weights do",
        "       not collapse -- fourteen half-octave bands, with the",
        "       exact errors of Lemma lem:cellmom. The rivals are",
        "       log y = A - b log x, A - c sqrt(x) and A - d x/log x",
        "       with x = log N, each two free parameters like the",
        "       power law itself. Nothing about the mask is",
        "       measured; the exponent and the (log N)^-b reach are",
        "       refitted only as J1's control against",
        "       v1/results/wall/lab_mask_exponent_se.txt and",
        "       results/audit_mask_deepform.txt.",
        "",
    ]
    io.open(OUT, "w", encoding="utf-8", newline="\n").write(
        "\n".join(head + lines) + "\n")
    print("\nwrote %s" % OUT)
    if not j1:
        raise SystemExit(1)
    return 0


if __name__ == "__main__":
    sys.exit(main())
