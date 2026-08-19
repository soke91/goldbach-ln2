# -*- coding: utf-8 -*-
r"""
Is rung 18's crossing of 0.56 an arithmetic fluctuation?

WHAT IS AT STAKE

[rem:rung18] measured 0.562768 at rung 18 and declared the level
crossed -- and on the same line declared `INSIDE FLOOR`, because the
clearance 0.002768 is below the ladder's scatter 0.0037.  So the
0.56 bracket (10^9.5951, 10^9.8961] has a resolved lower end and an
unresolved upper one, and G37 and G40 say that a single point
clearing a barrier by less than the scatter has not cleared it.

What that floor measures is worth being exact about, because it
decides what would resolve the upper end.  K*_R is a deterministic
integer once N is fixed; there is no sampling noise in it at all.
The ladder's scatter is the spread of its exponents about a fitted
curve -- **arithmetic fluctuation across N, not measurement error.**
So "the clearance is inside the floor" says: a neighbouring N of the
same odd radical might give an exponent 0.0037 lower, and that would
be below 0.56.

That is a statement about neighbouring N, and it is tested by
measuring neighbouring N.  [rem:primorialgap] did exactly this for
1/2: it filled the void between two rungs with N of the same odd
radical and closed the interval by the sign of measurements rather
than by a fit.

Four N are placed here within 0.021 decades of rung 18 on either
side, all of the form 30030*m with m composed only of 2, 3, 5, 7, 11
and 13 so that the odd radical, and hence the threshold
S(N)(1-A(N)), is the same as on the ladder.  Across a window that
narrow the trend contributes about 0.001, a quarter of the floor, so
what the five points (these four and rung 18) spread by is very
nearly the arithmetic fluctuation itself.

**That replaces the ladder-wide floor with a locally measured one**,
which is the right yardstick for "could a nearby N be lower", and it
is the only thing that can resolve the upper end without another
doubling.

The statistic, the sieve and the packing are imported verbatim from
code/audit_primorial_rung18.py, whose C1, C2 and C3 compared them
against the production route elementwise.

BACKS: Remark {#rem:rung18fill} in paper/theorem_A.md.

PRE-REGISTERED PREDICTIONS (written before this script was run)

  C1  The harness control.  Run through this script, rung 16 returns
      the K*_R and exponent that results/audit_ladder_cap.txt prints,
      126079 and 0.548808, exactly and to the six decimals printed.
  C2  The radical control.  The threshold S(N)(1-A(N))/N is the same
      at all five N to the six decimals printed, which is what
      "the same odd radical" has to mean if these points are on the
      ladder at all.
  F1  Every one of the four new N has exponent above 0.56.
  F2  The local scatter is smaller than the ladder's floor.  The
      r.m.s. residual of the five exponents about a line in log N
      fitted to them is below 0.0037.
  F3  **The crossing is resolved locally.**  The lowest of the five
      exponents clears 0.56 by more than that local scatter.

REFUTATION RULE (fixed before the run)

  C1  REFUTED by a K*_R that differs or an exponent differing in the
      six decimals printed.  THIS ONE GATES.
  C2  REFUTED by a threshold differing in the decimals printed.
      Then the five N are not on one ladder and F1 to F3 compare
      quantities that are not comparable.  THIS ONE GATES.
  F1  REFUTED by a single new N at or below 0.56.  **That is the
      outcome that would settle it the other way**: it would show
      directly that a neighbouring N can fall below the level, that
      rung 18's clearance is the fluctuation `INSIDE FLOOR` warned
      of, and that the upper end of the bracket is not merely
      unresolved but wrong.
  F2  REFUTED if the local scatter reaches 0.0037.  Then the ladder
      floor was the right yardstick after all, nothing has been
      gained by measuring locally, and the upper end stays where
      {#rem:rung18} left it.
  F3  REFUTED if the lowest of the five clears 0.56 by the local
      scatter or less.  Then the crossing is not resolved even
      against a local floor, and what remains is a doubling -- rung
      19 -- not more points.

  C1 and C2 gate.  F1 to F3 are the measurement and do not gate.

  NO NULL IS RUN and none applies.  A deterministic integer is
  located against a computed threshold at five N; there is no
  background to detect against.  The coin arms for this statistic
  were run in lab_primorial_ladder.py and lab_primorial_share.py.
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
OUT = os.path.join(RES, "audit_rung18_fill.txt")

CLIM = 4_000_000                    # the fixed Euler bound (G20)
TARGET = 0.56
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
CONTROL = BASE * (1 << 16)          # 1968046080, the rung 16 point
RUNG18 = BASE * (1 << 18)           # 7872184320
# m composed only of 2,3,5,7,11,13, so the odd radical is 3*5*7*11*13
FILL = (BASE * 250047,              # 7508911410, -0.0205 decades
        BASE * 256000,              # 7687680000, -0.0103
        BASE * 268800,              # 8072064000, +0.0109
        BASE * 275000)              # 8258250000, +0.0208
TOP = max(FILL)


def smooth_ok(n):
    """n's odd radical is 3*5*7*11*13 and nothing else"""
    return set(q for q in R11.factor_set(n) if q > 2) == {3, 5, 7, 11, 13}


def read_published():
    """rung 16 and rung 18 as the uniform ladder printed them"""
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
    m = re.search(r"^  N = " + str(RUNG18) +
                  r"\s+thr ([\d.]+)\s+#k \d+\s+beta [\d.]+\s+"
                  r"K\*_R (\d+)\s+exp ([\d.]+)\s*$", s18, re.M)
    got[RUNG18] = (int(m.group(2)), float(m.group(3)))
    thr18 = float(m.group(1))
    scat = float(re.search(r"^FLOOR primorial_rung18 ([\d.]+)\s*$",
                           s18, re.M).group(1))
    return got, thr18, scat


def main():
    lines = []

    def say(t=""):
        print(t)
        sys.stdout.flush()
        lines.append(t)

    pub, thr18, scat = read_published()
    say("read rung 16 and rung 18 from results/audit_ladder_cap.txt "
        "and")
    say("  results/audit_primorial_rung18.txt, with the ladder's "
        "floor %.4f" % scat)
    say("PRINTBOUND audit_rung18_fill %d %.8f"
        % (DEC, 0.5 * 10.0 ** (-DEC)))

    say()
    say("the five N in the window, all of odd radical 3*5*7*11*13")
    say("   N              log10 N   decades from rung 18")
    for N in sorted(FILL + (RUNG18,)):
        assert smooth_ok(N), "the odd radical moved at %d" % N
        say("   %-14d %.4f    %+.4f"
            % (N, math.log10(N),
               math.log10(N) - math.log10(RUNG18)))
    say("SCALES 1")
    say("  the window spans %.4f decades; the ladder's slope over it "
        "contributes"
        % (math.log10(max(FILL)) - math.log10(min(FILL))))
    say("  far less than the floor, so what these spread by is very "
        "nearly the")
    say("  arithmetic fluctuation itself")

    qs = [int(q) for q in primes_upto(R11.QSIEVE) if q > 2]
    say()
    say("sieving to %d on half-indices, once for every N below it"
        % TOP)
    kind, mu = R18.kind_and_mu_odd(TOP)
    vmask = R18.mask5_odd(TOP, qs)
    pv, plg = R18.power_table(TOP)
    say("BYTES resident_arrays %d"
        % (kind.nbytes + mu.nbytes + vmask.nbytes))
    say("  three arrays of %.2f GB, and a power table of %d entries"
        % ((kind.nbytes + mu.nbytes + vmask.nbytes) / 2.0 ** 30,
           pv.size))
    artin, twin = 1.0, 2.0
    assert CLIM == R11.CLIM
    for p in primes_upto(CLIM):
        p = int(p)
        artin *= 1.0 - 1.0 / (p * (p - 1.0))
        if p > 2:
            twin *= 1.0 - 1.0 / (p - 1.0) ** 2
    say("  the Euler products at the fixed bound %d: Artin %.9f, "
        "twin %.9f" % (CLIM, artin, twin))

    say()
    got = {}
    for N in (CONTROL,) + FILL:
        out = R18.measure_kind(N, kind, mu, vmask, pv, plg, qs, artin,
                               twin, CAP)
        if out is None:
            say("  N = %-14d no crossing below k = %d" % (N, CAP))
            continue
        kstar, e, bpn, beta, nk = out
        got[N] = (kstar, e, bpn)
        say("  N = %-14d thr %.6f  #k %-7d beta %.6f  K*_R %-8d "
            "exp %.6f" % (N, bpn, nk, beta, kstar, e))
        say("BUDGET kstar_R_S1AN_N%d %.6f" % (N, bpn))
    say("RADICALS 1")

    # -------------------------------------------------------------- C1
    say()
    say("C1  the harness control at rung 16")
    rnd = 0.5 * 10.0 ** (-DEC)
    pk, pe = pub[CONTROL]
    k_, e_, _ = got[CONTROL]
    c1 = k_ == pk and abs(e_ - pe) <= rnd
    say("  K*_R here %d against the published %d, exponent %.6f "
        "against %.6f, %s"
        % (k_, pk, e_, pe, "equal" if c1 else "DIFFERENT"))
    say("  C1 %s   (cap: exact on K*_R, the printing bound on the "
        "exponent)" % ("hold" if c1 else "REFUTED"))

    # -------------------------------------------------------------- C2
    say()
    say("C2  the radical control: one threshold at all five N")
    thrs = [thr18] + [got[N][2] for N in FILL if N in got]
    c2 = max(thrs) - min(thrs) <= rnd
    say("  thresholds %s, and rung 18's %.6f"
        % (", ".join("%.6f" % got[N][2] for N in FILL if N in got),
           thr18))
    say("  spread %.8f against the printing bound %.8f"
        % (max(thrs) - min(thrs), rnd))
    say("  C2 %s   (cap: the printing bound)"
        % ("hold" if c2 else "REFUTED"))
    if not (c1 and c2):
        io.open(OUT, "w", encoding="utf-8", newline="\n").write(
            "\n".join(lines) + "\n")
        raise SystemExit(1)

    # -------------------------------------------------------------- F1
    say()
    say("F1  does every new N clear 0.56?")
    es = [got[N][1] for N in FILL if N in got]
    below = [N for N in FILL if N in got and got[N][1] <= TARGET]
    f1 = len(es) == len(FILL) and not below
    for N in FILL:
        if N in got:
            say("  N = %-14d exp %.6f, %s 0.56 by %.6f"
                % (N, got[N][1],
                   "over" if got[N][1] > TARGET else "UNDER",
                   abs(got[N][1] - TARGET)))
    say("  F1 %s   (cap: the level 0.56)"
        % ("hold" if f1 else "REFUTED"))

    # -------------------------------------------------------------- F2
    say()
    say("F2  what is the local scatter?")
    allN = sorted(FILL + (RUNG18,))
    xs = np.array([math.log(N) for N in allN])
    ys = np.array([got[N][1] if N in got else pub[N][1]
                   for N in allN])
    sl, ic = np.polyfit(xs, ys, 1)
    resid = ys - (sl * xs + ic)
    local = float(np.sqrt((resid ** 2).mean()))
    f2 = local < scat
    say("  five exponents %s"
        % ", ".join("%.6f" % v for v in ys))
    say("  a line through them has slope %+.6f per log unit; the "
        "r.m.s. residual is %.6f" % (sl, local))
    say("FLOOR rung18_fill_local %.6f" % local)
    say("  against the ladder's floor %.4f" % scat)
    say("  F2 %s   (cap: the ladder's floor)"
        % ("hold" if f2 else "REFUTED"))

    # -------------------------------------------------------------- F3
    say()
    say("F3  is the crossing resolved against the local floor?")
    worst = float(ys.min())
    clear = worst - TARGET
    f3 = clear > local
    say("  the lowest of the five is %.6f, clearing 0.56 by %.6f"
        % (worst, clear))
    say("MARGIN audit_rung18_fill %.6f %.6f" % (clear, local))
    if clear <= local:
        say("INSIDE FLOOR audit_rung18_fill")
    say("  against the local floor %.6f, a ratio of %.2f"
        % (local, clear / local if local > 0 else float("inf")))
    say("  F3 %s   (cap: the local floor)"
        % ("hold" if f3 else "REFUTED"))

    say()
    say("=" * 70)
    say("C1 %s  C2 %s  F1 %s  F2 %s  F3 %s"
        % tuple("hold" if v else "REFUTED"
                for v in (c1, c2, f1, f2, f3)))

    head = [
        "STATISTIC: the truncation K*_R and its exponent",
        "           log K*_R / log N at four N of the ladder's own",
        "           odd radical placed within 0.021 decades of rung",
        "           18 on either side, together with rung 18 itself;",
        "           the r.m.s. residual of the five about a line in",
        "           log N, as a locally measured floor; and the",
        "           clearance of the lowest of them over 0.56",
        "           against that floor.  Rung 16 is recomputed as a",
        "           control on the harness.",
        "NULL: none is run and none applies. K*_R is a deterministic",
        "      integer once N is fixed; there is no sampling noise",
        "      and no background to detect against. What the floor",
        "      measures is arithmetic fluctuation across N, which is",
        "      why it is estimated here from neighbouring N rather",
        "      than from a fit across the whole ladder. The coin",
        "      arms for this statistic were run in",
        "      lab_primorial_ladder.py and lab_primorial_share.py.",
        "FIELD: N = 30030*m for m = 250047, 256000, 268800, 275000",
        "       and 262144, every m composed only of 2, 3, 5, 7, 11",
        "       and 13 so that the odd radical is 3*5*7*11*13 and",
        "       the threshold S(N)(1-A(N)) is the ladder's; and",
        "       N = 30030*2^16 as the control. k squarefree and",
        "       coprime to N with 2 <= k < 1000000, beta fitted on",
        "       the same range; m odd, squarefree and coprime to k,",
        "       m < N/k; the sieve weight over the odd primes below",
        "       30; the Euler products at the fixed bound 4000000.",
        "       One odd radical, as the RADICALS line declares. The",
        "       sieve, the packing and the statistic are imported",
        "       from code/audit_primorial_rung18.py; rung 16 comes",
        "       from results/audit_ladder_cap.txt and rung 18 with",
        "       the floor from results/audit_primorial_rung18.txt.",
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
