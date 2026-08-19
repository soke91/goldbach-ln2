# -*- coding: utf-8 -*-
r"""
The floor at 1/2, measured where the question lives.

WHAT IS AT STAKE

[rem:rung18fill] measured the arithmetic fluctuation between
neighbouring N of the ladder's odd radical, near rung 18, and found
0.000308 -- twelve times below the ladder's scatter 0.0037. The
ladder's scatter is the r.m.s. residual about a curve fitted across
ten decades, so at local scale it is almost all shape misfit rather
than arithmetic noise, and every single-point verdict this branch has
made against it has been judged against a floor too big for the
question.

Two such verdicts are at 1/2, and both retreated.

[rem:primorialrung10] read the barrier crossed at rung 10 by 0.0023
and had to withdraw the reading to a forecast, because 0.0023 is
below 0.0037. [rem:primorialgap] filled the void below it and found
the crossing inside the interval, but of its four interior points
only one, at +0.0040, cleared the floor; the point that brackets the
crossing from below is short by 0.0002 and the one above clears by
0.0030, so it had to call its own interval "a point-estimate bracket
and not a resolved one".

Neither retreat is wrong under the rule it was judged by. Both may be
unnecessary. This measures the local floor at both places, the same
way [rem:rung18fill] did: same-radical N placed in a narrow window,
a line through them, the r.m.s. residual as the floor.

**The floor is not assumed to be small here.** Arithmetic fluctuation
need not behave at 3*10^7 as it does at 8*10^9, and if it is larger
here the retreats stand and are confirmed rather than lifted. That is
why it is measured at both windows rather than carried over.

Everything is at the uniform cap 10^6 of [rem:laddercap], which is
the ladder the current rungs live on; the published gap numbers were
computed at the earlier cap, so the points are recomputed rather than
compared.

The sieve, the packing and the statistic are imported verbatim from
code/audit_primorial_rung18.py.

BACKS: Remark {#rem:localfloor} in paper/theorem_A.md.

PRE-REGISTERED PREDICTIONS (written before this script was run)

  C1  The ladder control. Rungs 9 and 10 return the K*_R and
      exponents results/audit_ladder_cap.txt prints at the uniform
      cap, 3551 with 0.494008 and 5779 with 0.502394, exactly and to
      the six decimals printed.
  C2  The radical control. The threshold S(N)(1-A(N))/N is the same
      at every N measured, to the six decimals printed.
  L1  Both local floors are below the ladder's 0.0037.
  L2  Rung 10's crossing of 1/2 is resolved locally: its clearance
      exceeds the local floor of the window around it. This is the
      retreat of {#rem:primorialrung10}.
  L3  The upper end of the gap bracket is resolved locally: at
      N = 20180160 the clearance over 1/2 exceeds the local floor of
      the window around it.
  L4  And so is the lower end: at N = 18018000 the shortfall below
      1/2 exceeds that same floor. **This is the one that is likely
      to fail** -- the published shortfall there is 0.0002, and a
      floor would have to be smaller than that to resolve it.

REFUTATION RULE (fixed before the run)

  C1  REFUTED by a K*_R that differs or an exponent differing in the
      six decimals printed. THIS ONE GATES.
  C2  REFUTED by a threshold differing in the decimals printed; the
      points would not be on one ladder. THIS ONE GATES.
  L1  REFUTED if either local floor reaches 0.0037. Then arithmetic
      fluctuation at 3*10^7 really is what the ladder scatter says,
      the retreats at 1/2 were necessary, and {#rem:rung18fill}'s
      finding does not generalise down the ladder -- which would be
      worth more than confirming it, because it would mean the floor
      is scale-dependent and every verdict needs its own.
  L2  REFUTED if rung 10's clearance is at or below its window's
      floor. {#rem:primorialrung10}'s retreat then stands as made.
  L3  REFUTED if the clearance at 20180160 is at or below the floor.
  L4  REFUTED if the shortfall at 18018000 is at or below the floor.
      Then the crossing of 1/2 is bracketed from below only at some
      smaller N, and the resolved bracket is wider than the
      point-estimate one -- which is the honest form of what
      {#rem:primorialgap} already suspected.

  C1 and C2 gate. L1 to L4 are the measurement and do not gate.

  NO NULL IS RUN and none applies. K*_R is a deterministic integer
  once N is fixed; there is no sampling noise and no background. What
  the floor measures is arithmetic fluctuation across N, which is
  why it is estimated from neighbouring N. The coin arms for this
  statistic were run in lab_primorial_ladder.py and
  lab_primorial_share.py.
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
OUT = os.path.join(RES, "audit_local_floor.txt")

CLIM = 4_000_000                    # the fixed Euler bound (G20)
TARGET = 0.5
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
primes_upto = R11.primes_upto
BASE = R11.BASE

RUNG9 = BASE * (1 << 9)             # 15375360
RUNG10 = BASE * (1 << 10)           # 30750720
GAP_LOW = BASE * 600                # 18018000, the published last short
GAP_HIGH = BASE * 672               # 20180160, the published first over

# window A: around rung 10, the crossing {#rem:primorialrung10} read
WIN_A = tuple(BASE * m for m in (975, 990, 1008, 1024, 1040, 1056))
# window B: around the gap bracket of {#rem:primorialgap}
WIN_B = tuple(BASE * m for m in (600, 616, 630, 640, 650, 660, 672,
                                 686))
TOP = max(WIN_A + WIN_B)


def smooth_ok(n):
    return set(q for q in R11.factor_set(n) if q > 2) == {3, 5, 7, 11,
                                                          13}


def read_uniform():
    """rungs 9 and 10 as the uniform ladder prints them"""
    src = io.open(os.path.join(RES, "audit_ladder_cap.txt"),
                  encoding="utf-8").read()
    got = {}
    for ln in src.splitlines():
        m = re.match(r"^  (\d+)\s+(\d+)\s+(.*)$", ln)
        if not m:
            continue
        f = [t for t in m.group(3).split() if t != "cap-invariant"]
        if len(f) < 8 or f[6] == "none":
            continue
        got[int(m.group(2))] = (int(f[6]), float(f[7]))
    s18 = io.open(os.path.join(RES, "audit_primorial_rung18.txt"),
                  encoding="utf-8").read()
    scat = float(re.search(r"^FLOOR primorial_rung18 ([\d.]+)\s*$",
                           s18, re.M).group(1))
    sf = io.open(os.path.join(RES, "audit_rung18_fill.txt"),
                 encoding="utf-8").read()
    far = float(re.search(r"^FLOOR rung18_fill_local ([\d.]+)\s*$",
                          sf, re.M).group(1))
    return got, scat, far


def local_floor(ns, es):
    """the r.m.s. residual of a window about a line in log N"""
    xs = np.log(np.array(ns, dtype=np.float64))
    ys = np.array(es)
    sl, ic = np.polyfit(xs, ys, 1)
    r = ys - (sl * xs + ic)
    return float(sl), float(np.sqrt((r ** 2).mean()))


def main():
    lines = []

    def say(t=""):
        print(t)
        sys.stdout.flush()
        lines.append(t)

    pub, scat, far = read_uniform()
    rnd = 0.5 * 10.0 ** (-DEC)
    say("read rungs 9 and 10 from results/audit_ladder_cap.txt, the "
        "ladder's floor")
    say("  %.4f from results/audit_primorial_rung18.txt, and the "
        "floor %.6f" % (scat, far))
    say("  measured near rung 18 from results/audit_rung18_fill.txt")
    say("PRINTBOUND audit_local_floor %d %.8f" % (DEC, rnd))

    for N in WIN_A + WIN_B + (RUNG9,):
        assert smooth_ok(N), "the odd radical moved at %d" % N

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

    todo = sorted(set(WIN_A + WIN_B + (RUNG9, RUNG10)))
    got = {}
    say()
    say("  N              log10 N   K*_R     exponent   over 1/2")
    for N in todo:
        out = R18.measure_kind(N, kind, mu, vmask, pv, plg, qs, artin,
                               twin, CAP)
        if out is None:
            say("  %-14d %-9.4f no crossing below k = %d" % (N, 0, CAP))
            continue
        kstar, e, bpn, beta, nk = out
        got[N] = (kstar, e, bpn)
        say("  %-14d %-9.4f %-8d %-10.6f %+.6f"
            % (N, math.log10(N), kstar, e, e - TARGET))
        say("BUDGET kstar_R_S1AN_N%d %.6f" % (N, bpn))
    say("RADICALS 1")

    # -------------------------------------------------------------- C1
    say()
    say("C1  the ladder control at rungs 9 and 10")
    c1 = True
    for j, N in ((9, RUNG9), (10, RUNG10)):
        pk, pe = pub[N]
        k_, e_, _ = got[N]
        ok = k_ == pk and abs(e_ - pe) <= rnd
        c1 = c1 and ok
        say("  rung %-3d K*_R here %d against %d, exponent %.6f "
            "against %.6f, %s"
            % (j, k_, pk, e_, pe, "equal" if ok else "DIFFERENT"))
    say("  C1 %s   (cap: exact on K*_R, the printing bound)"
        % ("hold" if c1 else "REFUTED"))

    # -------------------------------------------------------------- C2
    say()
    say("C2  the radical control: one threshold everywhere")
    thrs = [got[N][2] for N in todo if N in got]
    c2 = max(thrs) - min(thrs) <= rnd
    say("  the threshold reads %.6f at all %d N, spread %.8f"
        % (thrs[0], len(thrs), max(thrs) - min(thrs)))
    say("  C2 %s   (cap: the printing bound)"
        % ("hold" if c2 else "REFUTED"))
    if not (c1 and c2):
        io.open(OUT, "w", encoding="utf-8", newline="\n").write(
            "\n".join(lines) + "\n")
        raise SystemExit(1)

    # -------------------------------------------------------------- L1
    say()
    say("L1  the local floor in each window")
    slA, floA = local_floor(WIN_A, [got[N][1] for N in WIN_A])
    slB, floB = local_floor(WIN_B, [got[N][1] for N in WIN_B])
    say("  window A, %d points over %.4f decades around rung 10:"
        % (len(WIN_A), math.log10(max(WIN_A)) - math.log10(min(WIN_A))))
    say("    slope %+.6f per log unit, r.m.s. residual %.6f"
        % (slA, floA))
    say("FLOOR local_window_A %.6f" % floA)
    say("  window B, %d points over %.4f decades around the gap "
        "bracket:"
        % (len(WIN_B), math.log10(max(WIN_B)) - math.log10(min(WIN_B))))
    say("    slope %+.6f per log unit, r.m.s. residual %.6f"
        % (slB, floB))
    say("FLOOR local_window_B %.6f" % floB)
    l1 = floA < scat and floB < scat
    say("  against the ladder's floor %.4f, and the %.6f measured "
        "near rung 18" % (scat, far))
    say("  L1 %s   (cap: the ladder's floor)"
        % ("hold" if l1 else "REFUTED"))

    # -------------------------------------------------------------- L2
    say()
    say("L2  is rung 10's crossing of 1/2 resolved locally?")
    c10 = got[RUNG10][1] - TARGET
    l2 = c10 > floA
    say("  rung 10 clears 1/2 by %.6f against the local floor %.6f, "
        "a ratio of %.2f" % (c10, floA, c10 / floA))
    say("MARGIN audit_local_floor_rung10 %.6f %.6f" % (c10, floA))
    if c10 <= floA:
        say("INSIDE FLOOR audit_local_floor_rung10")
    say("  L2 %s   (cap: window A's floor)"
        % ("hold" if l2 else "REFUTED"))

    # --------------------------------------------------------- L3, L4
    say()
    say("L3/L4  are the two ends of the gap bracket resolved?")
    hi = got[GAP_HIGH][1] - TARGET
    lo = TARGET - got[GAP_LOW][1]
    l3 = hi > floB
    l4 = lo > floB
    say("  N = %d clears 1/2 by %+.6f, ratio %.2f"
        % (GAP_HIGH, hi, hi / floB))
    say("  N = %d falls short by %+.6f, ratio %.2f"
        % (GAP_LOW, lo, lo / floB))
    say("MARGIN audit_local_floor_gaphigh %.6f %.6f" % (hi, floB))
    if hi <= floB:
        say("INSIDE FLOOR audit_local_floor_gaphigh")
    say("MARGIN audit_local_floor_gaplow %.6f %.6f" % (lo, floB))
    if lo <= floB:
        say("INSIDE FLOOR audit_local_floor_gaplow")
    say("  L3 %s   (cap: window B's floor)"
        % ("hold" if l3 else "REFUTED"))
    say("  L4 %s   (cap: window B's floor)"
        % ("hold" if l4 else "REFUTED"))

    # ---------------------------------------------- the resolved pair
    say()
    say("the resolved bracket for 1/2, from window B")
    below = [N for N in WIN_B
             if TARGET - got[N][1] > floB]
    above = [N for N in WIN_B if got[N][1] - TARGET > floB]
    if below and above:
        b, a = max(below), min(above)
        say("  highest N resolved below: %d at 10^%.4f, exponent "
            "%.6f" % (b, math.log10(b), got[b][1]))
        say("  lowest N resolved above:  %d at 10^%.4f, exponent "
            "%.6f" % (a, math.log10(a), got[a][1]))
        say("BRACKET half_resolved %.4f %.4f %.4f"
            % ((math.log10(b) + math.log10(a)) / 2.0,
               math.log10(b), math.log10(a)))
        say("DRIFT half_resolved %.4f"
            % abs(math.log10(a) - math.log10(b)))
        say("SCATTER slope_audit_local_floor %.6f" % floB)
        say("SHAPES 1")
    else:
        say("  no pair of window-B points is resolved on both sides; "
            "the bracket stays")
        say("  a point-estimate one")
    say("SCALES 2")

    say()
    say("=" * 70)
    say("C1 %s  C2 %s  L1 %s  L2 %s  L3 %s  L4 %s"
        % tuple("hold" if v else "REFUTED"
                for v in (c1, c2, l1, l2, l3, l4)))

    head = [
        "STATISTIC: the truncation K*_R and its exponent",
        "           log K*_R / log N at same-radical N in two narrow",
        "           windows, one around rung 10 and one around the",
        "           interval {#rem:primorialgap} bracketed; the",
        "           r.m.s. residual of each window about a line in",
        "           log N, as a locally measured floor; and the",
        "           clearance or shortfall of three published",
        "           verdicts at 1/2 against those floors. Rungs 9",
        "           and 10 are recomputed as a control.",
        "NULL: none is run and none applies. K*_R is a deterministic",
        "      integer once N is fixed; there is no sampling noise",
        "      and no background to detect against. What the floor",
        "      measures is arithmetic fluctuation across N, which is",
        "      why it is estimated from neighbouring N rather than",
        "      from a fit across the whole ladder. The coin arms for",
        "      this statistic were run in lab_primorial_ladder.py",
        "      and lab_primorial_share.py.",
        "FIELD: N = 30030*m, every m composed only of 2, 3, 5, 7, 11",
        "       and 13 so the odd radical is 3*5*7*11*13 and the",
        "       threshold is the ladder's: m = 975, 990, 1008, 1024,",
        "       1040, 1056 in window A and m = 600, 616, 630, 640,",
        "       650, 660, 672, 686 in window B, with m = 512 and",
        "       1024 as controls. k squarefree and coprime to N with",
        "       2 <= k < 1000000, beta fitted on the same range; m",
        "       odd, squarefree and coprime to k, m < N/k; the sieve",
        "       weight over the odd primes below 30; the Euler",
        "       products at the fixed bound 4000000. One odd",
        "       radical, as the RADICALS line declares. The sieve,",
        "       the packing and the statistic are imported from",
        "       code/audit_primorial_rung18.py; rungs 9 and 10 come",
        "       from results/audit_ladder_cap.txt, the ladder floor",
        "       from results/audit_primorial_rung18.txt and the",
        "       rung-18 local floor from",
        "       results/audit_rung18_fill.txt.",
        "",
    ]
    io.open(OUT, "w", encoding="utf-8", newline="\n").write(
        "\n".join(head + lines) + "\n")
    print("\nwrote %s" % OUT)
    if not (c1 and c2):
        raise SystemExit(1)
    return 0


if __name__ == "__main__":
    sys.exit(main())
