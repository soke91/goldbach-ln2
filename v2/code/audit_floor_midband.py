# -*- coding: utf-8 -*-
r"""
The floor law, asked for a window it was not fitted on.

WHAT IS AT STAKE

{#rem:floorlaw} fitted log(floor) on log(N) over the four windows
measured so far and got a slope of -0.321639, roughly N^(-1/3).  It
also said, on the line beneath, why the printed t = -9.60 is not nine
standard errors: the four windows sit in two clusters 2.6 decades
apart, so the slope is in substance a two-point determination.  A
two-point line through two clusters is not a law; it is a line.

The way to tell them apart is a third place.  10^8.5 sits almost
exactly between the clusters -- between rung 13 at 10^8.3909 and rung
14 at 10^8.6920 -- and at N near 3*10^8 a window costs minutes rather
than hours.  So the law is asked for a number before the number is
measured.

Fitted on the four, it puts the floor at the geometric mean N of the
window below at **0.000767**.  The four-window fit's own r.m.s.
residual is 0.1334 in log, which is also about the scatter between
two floors measured at the same N (0.001883 against 0.001662 at one
cluster, 0.000308 against 0.000264 at the other).  So a new window
landing within twice that -- a factor of 1.31 -- is the law
predicting; landing outside it is two clusters and a line drawn
through them.

Nothing about the ladder changes either way.  What changes is whether
a window can be sized before the compute is spent, which is the thing
{#rem:targetband} showed decides between resolving a crossing and
sampling a band from inside it.

The sieve, the packing and the statistic are imported verbatim from
code/audit_primorial_rung18.py.  The four floors and the window
abscissae are read the way {#rem:floorlaw} reads them: whole marker
lines, declared with READ so G76 checks them, and window definitions
imported from the scripts that made them.

BACKS: Remark {#rem:floormid} in paper/theorem_A.md.

PRE-REGISTERED PREDICTIONS (written before this script was run)

  Q1  The ladder control.  Rung 13 returns the K*_R and exponent
      results/audit_ladder_cap.txt prints at the uniform cap, 27343
      and 0.528766, exactly and to the six decimals printed.
  Q2  The radical control.  The threshold S(N)(1-A(N))/N is the same
      at every N measured, to the six decimals printed.
  Q3  **The law predicts it.**  The floor measured in this window is
      within a factor 1.31 of 0.000767 -- that is, |log(measured /
      predicted)| is below 0.2668, twice the four-window fit's own
      r.m.s. residual.
  Q4  And it lies between the clusters: the measured floor is below
      0.001662 and above 0.000308.
  Q5  Refitting on five windows keeps the slope negative and
      resolved at |t| > 2.

REFUTATION RULE (fixed before the run)

  Q1  REFUTED by a K*_R that differs or an exponent differing in the
      six decimals printed.  THIS ONE GATES.
  Q2  REFUTED by a threshold differing in the decimals printed; the
      points would not be on one ladder.  THIS ONE GATES.
  Q3  REFUTED outside a factor 1.31.  **That is the outcome worth
      having and it is not unlikely**: it would say the four floors
      are four numbers, that the two clusters happen to lie on a
      line, and that no window anywhere can be sized in advance --
      every one would have to be measured and then judged, which is
      what {#rem:primorialgap} did and paid for.
  Q4  REFUTED if the floor falls outside the clusters' own values.
      Then it is not even monotone in N and the fit was fitting
      noise.
  Q5  REFUTED if the slope turns positive or |t| drops to 2.  With
      the middle filled the t is no longer inflated by clustering,
      so this is the first honest reading of it -- and a slope that
      stops resolving once the clustering is broken would mean the
      earlier t was the clustering and nothing else.

  Q1 and Q2 gate.  Q3 to Q5 are the measurement and do not gate.

  NO NULL IS RUN and none applies.  K*_R is a deterministic integer
  once N is fixed; there is no sampling noise and no background to
  detect against.  What the floor measures is arithmetic fluctuation
  across N, which is why it is estimated from neighbouring N.  The
  coin arms for this statistic were run in lab_primorial_ladder.py
  and lab_primorial_share.py.
"""

import importlib.util
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
CODE = os.path.join(ROOT, "code")
RES = os.path.join(ROOT, "results")
OUT = os.path.join(RES, "audit_floor_midband.txt")

CLIM = 4_000_000                    # the fixed Euler bound (G20)
CAP = 1_000_000                     # the uniform cap of {#rem:laddercap}
DEC = 6                             # the decimals the ladder prints


def module(name):
    p = os.path.join(CODE, name + ".py")
    spec = importlib.util.spec_from_file_location(name, p)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


R11 = module("audit_primorial_rung11")
R18 = module("audit_primorial_rung18")
LF = module("audit_local_floor")
FI = module("audit_rung18_fill")
TB = module("audit_target_band")
primes_upto = R11.primes_upto
BASE = R11.BASE
CONTROL = BASE * (1 << 13)          # 246005760, the rung 13 point
MID = tuple(BASE * m for m in (10080, 10206, 10395, 10560, 10725,
                               10890, 11025))
TOP = max(MID + (CONTROL,))


def smooth_ok(n):
    return set(q for q in R11.factor_set(n) if q > 2) == {3, 5, 7, 11,
                                                          13}


def marker(fname, label):
    src = io.open(os.path.join(RES, fname), encoding="utf-8").read()
    m = re.search(r"^%s ([-+]?[\d.]+)\s*$" % re.escape(label), src,
                  re.M)
    if not m:
        raise SystemExit("no line '%s ...' in %s" % (label, fname))
    return float(m.group(1))


def gmean(ns):
    return math.exp(sum(math.log(float(n)) for n in ns) / len(ns))


def local_floor(ns, es):
    xs = np.log(np.array(ns, dtype=np.float64))
    ys = np.array(es)
    sl, ic = np.polyfit(xs, ys, 1)
    r = ys - (sl * xs + ic)
    return float(sl), float(np.sqrt((r ** 2).mean()))


def linefit(x, y):
    sl, ic = np.polyfit(x, y, 1)
    n = x.size
    r = y - (sl * x + ic)
    sse = float((r ** 2).sum()) / (n - 2)
    sxx = float(((x - x.mean()) ** 2).sum())
    se = math.sqrt(sse / sxx)
    return sl, ic, se, float(np.sqrt((r ** 2).mean()))


def main():
    lines = []

    def say(t=""):
        print(t)
        sys.stdout.flush()
        lines.append(t)

    floA = marker("audit_local_floor.txt", "FLOOR local_window_A")
    floB = marker("audit_local_floor.txt", "FLOOR local_window_B")
    floF = marker("audit_rung18_fill.txt", "FLOOR rung18_fill_local")
    floT = marker("audit_target_band.txt", "FLOOR target_band_local")
    scat4 = marker("audit_floor_law.txt", "SCATTER slope_audit_floor_law")

    say("the four floors the law was fitted on, each read from a "
        "whole marker line:")
    say("READ audit_local_floor.txt FLOOR local_window_A %.6f" % floA)
    say("READ audit_local_floor.txt FLOOR local_window_B %.6f" % floB)
    say("READ audit_target_band.txt FLOOR target_band_local %.6f"
        % floT)
    say("READ audit_rung18_fill.txt FLOOR rung18_fill_local %.6f"
        % floF)
    say("READ audit_floor_law.txt SCATTER slope_audit_floor_law %.4f"
        % scat4)

    old = [(gmean(LF.WIN_A), floA), (gmean(LF.WIN_B), floB),
           (gmean(TB.BAND), floT),
           (gmean(FI.FILL + (FI.RUNG18,)), floF)]
    x4 = np.array([math.log(n) for n, _ in old])
    y4 = np.array([math.log(f) for _, f in old])
    sl4, ic4, se4, rms4 = linefit(x4, y4)
    gN = gmean(MID)
    pred = math.exp(sl4 * math.log(gN) + ic4)
    tol = 2.0 * scat4
    say()
    say("the four-window law is log(floor) = %+.6f log N %+.6f, "
        "r.m.s. %.4f" % (sl4, ic4, rms4))
    say("  at this window's geometric mean N = %d (log10 %.4f) it "
        "predicts %.6f" % (int(gN), math.log10(gN), pred))
    say("PREDICTED floor_midband %.6f" % pred)
    say("  the registered tolerance is twice the four-window r.m.s., "
        "%.4f in log," % tol)
    say("  a factor of %.2f" % math.exp(tol))

    for N in MID:
        assert smooth_ok(N), "the odd radical moved at %d" % N
    say()
    say("the seven N in the window, all of odd radical 3*5*7*11*13")
    say("   N            log10 N")
    for N in MID:
        say("   %-12d %.4f" % (N, math.log10(N)))
    say("  spanning %.4f decades"
        % (math.log10(MID[-1]) - math.log10(MID[0])))
    say("SCALES 1")

    qs = [int(q) for q in primes_upto(R11.QSIEVE) if q > 2]
    say()
    say("sieving to %d on half-indices, once for every N below it"
        % TOP)
    kind, mu = R18.kind_and_mu_odd(TOP)
    vmask = R18.mask5_odd(TOP, qs)
    pv, plg = R18.power_table(TOP)
    say("BYTES resident_arrays %d"
        % (kind.nbytes + mu.nbytes + vmask.nbytes))
    artin, twin = 1.0, 2.0
    assert CLIM == R11.CLIM
    for p in primes_upto(CLIM):
        p = int(p)
        artin *= 1.0 - 1.0 / (p * (p - 1.0))
        if p > 2:
            twin *= 1.0 - 1.0 / (p - 1.0) ** 2
    say("  the Euler products at the fixed bound %d: Artin %.9f, "
        "twin %.9f" % (CLIM, artin, twin))

    got = {}
    say()
    say("  N            log10 N   K*_R     exponent")
    for N in (CONTROL,) + MID:
        out = R18.measure_kind(N, kind, mu, vmask, pv, plg, qs, artin,
                               twin, CAP)
        if out is None:
            say("  %-12d no crossing below k = %d" % (N, CAP))
            continue
        kstar, e, bpn, beta, nk = out
        got[N] = (kstar, e, bpn)
        say("  %-12d %-9.4f %-8d %.6f"
            % (N, math.log10(N), kstar, e))
        say("BUDGET kstar_R_S1AN_N%d %.6f" % (N, bpn))
    say("RADICALS 1")

    # -------------------------------------------------------------- Q1
    rnd = 0.5 * 10.0 ** (-DEC)
    say()
    say("Q1  the ladder control at rung 13")
    say("PRINTBOUND audit_floor_midband %d %.8f" % (DEC, rnd))
    src = io.open(os.path.join(RES, "audit_ladder_cap.txt"),
                  encoding="utf-8").read()
    pk = pe = None
    for ln in src.splitlines():
        m = re.match(r"^  (\d+)\s+(\d+)\s+(.*)$", ln)
        if m and int(m.group(2)) == CONTROL:
            f = [t for t in m.group(3).split() if t != "cap-invariant"]
            pk, pe = int(f[6]), float(f[7])
    k_, e_, _ = got[CONTROL]
    q1 = pk is not None and k_ == pk and abs(e_ - pe) <= rnd
    say("  K*_R here %d against %s, exponent %.6f against %s, %s"
        % (k_, pk, e_, pe, "equal" if q1 else "DIFFERENT"))
    say("  Q1 %s   (cap: exact on K*_R, the printing bound)"
        % ("hold" if q1 else "REFUTED"))

    # -------------------------------------------------------------- Q2
    say()
    say("Q2  the radical control: one threshold everywhere")
    thrs = [got[N][2] for N in (CONTROL,) + MID if N in got]
    q2 = max(thrs) - min(thrs) <= rnd
    say("  the threshold reads %.6f at all %d N, spread %.8f"
        % (thrs[0], len(thrs), max(thrs) - min(thrs)))
    say("  Q2 %s   (cap: the printing bound)"
        % ("hold" if q2 else "REFUTED"))
    if not (q1 and q2):
        io.open(OUT, "w", encoding="utf-8", newline="\n").write(
            "\n".join(lines) + "\n")
        raise SystemExit(1)

    # -------------------------------------------------------------- Q3
    say()
    say("Q3  does the law predict this window's floor?")
    slope, flo = local_floor(MID, [got[N][1] for N in MID])
    d = math.log(flo / pred)
    q3 = abs(d) < tol
    say("  the window's slope is %+.6f per log unit and its r.m.s. "
        "residual %.6f" % (slope, flo))
    say("FLOOR midband_local %.6f" % flo)
    say("  measured %.6f against predicted %.6f: log ratio %+.4f, "
        "a factor of %.2f" % (flo, pred, d, math.exp(abs(d))))
    say("  against the registered tolerance %.4f in log" % tol)
    say("  Q3 %s   (cap: twice the four-window r.m.s.)"
        % ("hold" if q3 else "REFUTED"))

    # -------------------------------------------------------------- Q4
    say()
    say("Q4  does it lie between the clusters?")
    q4 = floF < flo < floB
    say("  %.6f (rung 18 window) < %.6f (here) < %.6f (the gap "
        "window)?  %s" % (floF, flo, floB, "yes" if q4 else "NO"))
    say("  Q4 %s   (cap: the two clusters' own floors)"
        % ("hold" if q4 else "REFUTED"))

    # -------------------------------------------------------------- Q5
    say()
    say("Q5  does the slope survive breaking the clustering?")
    x5 = np.append(x4, math.log(gN))
    y5 = np.append(y4, math.log(flo))
    sl5, ic5, se5, rms5 = linefit(x5, y5)
    t5 = sl5 / se5
    q5 = sl5 < 0.0 and abs(t5) > 2.0
    say("  on five windows the slope is %+.6f +- %.6f, t = %.2f, "
        "r.m.s. %.4f" % (sl5, se5, t5, rms5))
    say("  the four-window values were %+.6f and r.m.s. %.4f"
        % (sl4, rms4))
    say("TSTAT floor_law5_slope %.2f" % t5)
    if abs(t5) < 2.0:
        say("UNRESOLVED SIGN floor_law5_slope")
    say("SPREAD floor_law5_slope %.4f" % (x5.max() - x5.min()))
    say("SCATTER slope_audit_floor_midband %.4f" % rms5)
    say("  Q5 %s   (cap: negative and |t| > 2)"
        % ("hold" if q5 else "REFUTED"))

    say()
    say("=" * 70)
    say("Q1 %s  Q2 %s  Q3 %s  Q4 %s  Q5 %s"
        % tuple("hold" if v else "REFUTED"
                for v in (q1, q2, q3, q4, q5)))

    head = [
        "STATISTIC: the local floor -- the r.m.s. residual of seven",
        "           same-radical N about a line in log N -- in a",
        "           window at 10^8.5, between the two clusters the",
        "           floor law of {#rem:floorlaw} was fitted on;",
        "           against what that law, fitted without this",
        "           window, predicts there; and the slope refitted",
        "           on five windows with the clustering broken.",
        "           Rung 13 is recomputed as a control.",
        "NULL: none is run and none applies. K*_R is a deterministic",
        "      integer once N is fixed; there is no sampling noise",
        "      and no background to detect against. What the floor",
        "      measures is arithmetic fluctuation across N, which is",
        "      why it is estimated from neighbouring N. The coin",
        "      arms for this statistic were run in",
        "      lab_primorial_ladder.py and lab_primorial_share.py.",
        "FIELD: N = 30030*m for m = 10080, 10206, 10395, 10560,",
        "       10725, 10890, 11025, every m composed only of 2, 3,",
        "       5, 7, 11 and 13 so the odd radical is 3*5*7*11*13 and",
        "       the threshold is the ladder's; and N = 30030*2^13 as",
        "       the control. k squarefree and coprime to N with",
        "       2 <= k < 1000000, beta fitted on the same range; m",
        "       odd, squarefree and coprime to k, m < N/k; the sieve",
        "       weight over the odd primes below 30; the Euler",
        "       products at the fixed bound 4000000. One odd",
        "       radical, as the RADICALS line declares. The sieve,",
        "       the packing and the statistic are imported from",
        "       code/audit_primorial_rung18.py; the four earlier",
        "       floors and the fit's scatter are read from",
        "       results/audit_local_floor.txt,",
        "       results/audit_target_band.txt,",
        "       results/audit_rung18_fill.txt and",
        "       results/audit_floor_law.txt, and the window",
        "       abscissae are imported from the scripts that defined",
        "       them; rung 13 comes from",
        "       results/audit_ladder_cap.txt.",
        "",
    ]
    io.open(OUT, "w", encoding="utf-8", newline="\n").write(
        "\n".join(head + lines) + "\n")
    print("\nwrote %s" % OUT)
    if not (q1 and q2):
        raise SystemExit(1)
    return 0


if __name__ == "__main__":
    sys.exit(main())
