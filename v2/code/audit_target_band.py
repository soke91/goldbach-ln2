# -*- coding: utf-8 -*-
r"""
Is 0.56 crossed once, or does the ladder sit on it as it does on 1/2?

WHAT IS AT STAKE

[rem:localfloor] put eight same-radical N across the interval
[rem:primorialgap] had bracketed at 1/2 and found the sign of
exponent - 1/2 changing five times, with six of the eight inside the
local floor.  It read that as an oscillation, and corrected
[rem:primorialgap]'s "the crossing is inside the interval" to "the
ladder sits on 1/2 within the fluctuation".

The arithmetic behind that is worth stating, because it says where to
look next rather than leaving the correction as a warning.  A point
is indeterminate when |exponent - level| is under the local floor, so
the indeterminate band is about 2*floor/slope wide.  At 1/2 the
window's slope was 0.024973 per log unit and its floor 0.001662,
putting the band at 0.0578 decades -- and the window spanned 0.0582.
**The window was the band.** Every point in it had to be
indeterminate; the oscillation was not a discovery about the ladder
so much as the shape of sampling a band from inside.

The same arithmetic at 0.56 gives a narrower band.  Rung 17 at
10^9.5951 reads 0.555925 and the lowest point of {#rem:rung18fill} at
10^9.8756 reads 0.563206, a chord of 0.025957 per decade; against the
0.000308 floor measured near rung 18 that is 0.0237 decades, and
against the flatter local slope there, 0.0465.  The band is somewhere
between, and the chord puts the level itself at 10^9.7521.

So seven same-radical N are placed across 0.0641 decades centred
there, spaced about 0.0107 -- fine enough to sample a band of that
width three to five times.  The question is whether the sign of
exponent - 0.56 turns over once, as a crossing, or repeatedly, as at
1/2.

This is the interior [rem:rung18fill] left unmeasured, and the only
thing that decides whether "0.56 is crossed in
(10^9.5951, 10^9.8756]" is a statement about a crossing or about a
region.

The sieve, the packing and the statistic are imported verbatim from
code/audit_primorial_rung18.py.

BACKS: Remark {#rem:targetband} in paper/theorem_A.md.

PRE-REGISTERED PREDICTIONS (written before this script was run)

  C1  The ladder control.  Rung 16 returns the K*_R and exponent
      results/audit_ladder_cap.txt prints, 126079 and 0.548808,
      exactly and to the six decimals printed.
  C2  The radical control.  The threshold S(N)(1-A(N))/N is the same
      at every N measured, to the six decimals printed.
  W1  The floor keeps shrinking.  The local floor here is below the
      0.001662 measured at 1/2 by {#rem:localfloor}.
  W2  The window straddles the level: at least one of the seven is
      below 0.56 and at least one is above.
  W3  **The sign turns over exactly once.**  Ordered by N, the sign
      of exponent - 0.56 changes once across the seven points.
  W4  And the band here is narrower than at 1/2: 2*floor/slope,
      from this window's own floor and slope, is below the 0.0578
      decades that arithmetic gives at 1/2.

REFUTATION RULE (fixed before the run)

  C1  REFUTED by a K*_R that differs or an exponent differing in the
      six decimals printed.  THIS ONE GATES.
  C2  REFUTED by a threshold differing in the decimals printed; the
      points would not be on one ladder.  THIS ONE GATES.
  W1  REFUTED if the floor reaches 0.001662.  Then the fluctuation
      does not shrink with N the way {#rem:localfloor} read it, and
      every floor on this ladder has to be measured rather than
      expected.
  W2  REFUTED if all seven fall on one side.  The level is then
      outside this window, which would mean the chord between two
      resolved measurements misplaces it by more than 0.03 decades
      -- itself a fact about how far from a line the ladder runs
      across the gap.
  W3  REFUTED if the sign changes more than once.  **That is the
      outcome that matters**: it would put 0.56 in the same position
      as 1/2, where the ladder sits on the level rather than passing
      through it, and "0.56 is crossed in
      (10^9.5951, 10^9.8756]" would have to be restated as a region
      the ladder is indistinguishable from the level on, exactly as
      {#rem:localfloor} had to restate {#rem:primorialgap}.  A sign
      that changes zero times is W2's failure, not this one's.
  W4  REFUTED if the band is 0.0578 decades or wider.  Then nothing
      about 0.56 is better determined than 1/2 was, and W3 holding
      would be luck rather than resolution.

  C1 and C2 gate.  W1 to W4 are the measurement and do not gate.

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
OUT = os.path.join(RES, "audit_target_band.txt")

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
BAND = tuple(BASE * m for m in (174960, 179200, 183708, 188160,
                                193050, 198000, 202800))
TOP = max(BAND)


def smooth_ok(n):
    return set(q for q in R11.factor_set(n) if q > 2) == {3, 5, 7, 11,
                                                          13}


def read_context():
    """rung 16, the floor at 1/2, and the two resolved endpoints"""
    src = io.open(os.path.join(RES, "audit_ladder_cap.txt"),
                  encoding="utf-8").read()
    pub = {}
    for ln in src.splitlines():
        m = re.match(r"^  (\d+)\s+(\d+)\s+(.*)$", ln)
        if not m:
            continue
        f = [t for t in m.group(3).split() if t != "cap-invariant"]
        if len(f) < 8 or f[6] == "none":
            continue
        pub[int(m.group(2))] = (int(f[6]), float(f[7]))
    sl = io.open(os.path.join(RES, "audit_local_floor.txt"),
                 encoding="utf-8").read()
    half = float(re.search(r"^FLOOR local_window_B ([\d.]+)\s*$",
                           sl, re.M).group(1))
    halfslope = float(re.search(r"^SLOPE local_window_B "
                                r"([-+][\d.]+)\s*$",
                                sl, re.M).group(1))
    sf = io.open(os.path.join(RES, "audit_rung18_fill.txt"),
                 encoding="utf-8").read()
    far = float(re.search(r"^FLOOR rung18_fill_local ([\d.]+)\s*$",
                          sf, re.M).group(1))
    s17 = io.open(os.path.join(RES, "audit_primorial_rung17.txt"),
                  encoding="utf-8").read()
    lo = float(re.search(r"^  N = " + str(BASE * (1 << 17)) +
                         r"\s+thr [\d.]+\s+#k \d+\s+beta [\d.]+\s+"
                         r"K\*_R \d+\s+exp ([\d.]+)\s*$",
                         s17, re.M).group(1))
    return pub, half, halfslope, far, lo


def local_floor(ns, es):
    xs = np.log(np.array(ns, dtype=np.float64))
    ys = np.array(es)
    sl, ic = np.polyfit(xs, ys, 1)
    r = ys - (sl * xs + ic)
    return float(sl), float(np.sqrt((r ** 2).mean()))


def turns(vals):
    """how many times a sequence of signs changes"""
    s = [1 if v > 0 else -1 for v in vals]
    return sum(1 for i in range(1, len(s)) if s[i] != s[i - 1])


def main():
    lines = []

    def say(t=""):
        print(t)
        sys.stdout.flush()
        lines.append(t)

    pub, halffloor, halfslope, farfloor, e17 = read_context()
    rnd = 0.5 * 10.0 ** (-DEC)
    say("read rung 16 from results/audit_ladder_cap.txt, the floor "
        "%.6f at 1/2 from" % halffloor)
    say("  results/audit_local_floor.txt, the floor %.6f near rung 18 "
        "from" % farfloor)
    say("  results/audit_rung18_fill.txt, and rung 17's exponent "
        "%.6f from" % e17)
    say("  results/audit_primorial_rung17.txt")
    say("PRINTBOUND audit_target_band %d %.8f" % (DEC, rnd))
    say("READ audit_local_floor.txt FLOOR local_window_B %.6f"
        % halffloor)
    say("READ audit_local_floor.txt SLOPE local_window_B %+.6f"
        % halfslope)
    say("READ audit_rung18_fill.txt FLOOR rung18_fill_local %.6f"
        % farfloor)

    for N in BAND:
        assert smooth_ok(N), "the odd radical moved at %d" % N

    say()
    say("the seven N across the band, all of odd radical 3*5*7*11*13")
    say("   N              log10 N")
    for N in BAND:
        say("   %-14d %.4f" % (N, math.log10(N)))
    say("  spanning %.4f decades, spaced about %.4f"
        % (math.log10(BAND[-1]) - math.log10(BAND[0]),
           (math.log10(BAND[-1]) - math.log10(BAND[0]))
           / (len(BAND) - 1)))
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
    say("  N              log10 N   K*_R     exponent   over 0.56")
    for N in (CONTROL,) + BAND:
        out = R18.measure_kind(N, kind, mu, vmask, pv, plg, qs, artin,
                               twin, CAP)
        if out is None:
            say("  %-14d no crossing below k = %d" % (N, CAP))
            continue
        kstar, e, bpn, beta, nk = out
        got[N] = (kstar, e, bpn)
        say("  %-14d %-9.4f %-8d %-10.6f %+.6f"
            % (N, math.log10(N), kstar, e, e - TARGET))
        say("BUDGET kstar_R_S1AN_N%d %.6f" % (N, bpn))
    say("RADICALS 1")

    # -------------------------------------------------------------- C1
    say()
    say("C1  the ladder control at rung 16")
    pk, pe = pub[CONTROL]
    k_, e_, _ = got[CONTROL]
    c1 = k_ == pk and abs(e_ - pe) <= rnd
    say("  K*_R here %d against %d, exponent %.6f against %.6f, %s"
        % (k_, pk, e_, pe, "equal" if c1 else "DIFFERENT"))
    say("  C1 %s   (cap: exact on K*_R, the printing bound)"
        % ("hold" if c1 else "REFUTED"))

    # -------------------------------------------------------------- C2
    say()
    say("C2  the radical control: one threshold everywhere")
    thrs = [got[N][2] for N in (CONTROL,) + BAND if N in got]
    c2 = max(thrs) - min(thrs) <= rnd
    say("  the threshold reads %.6f at all %d N, spread %.8f"
        % (thrs[0], len(thrs), max(thrs) - min(thrs)))
    say("  C2 %s   (cap: the printing bound)"
        % ("hold" if c2 else "REFUTED"))
    if not (c1 and c2):
        io.open(OUT, "w", encoding="utf-8", newline="\n").write(
            "\n".join(lines) + "\n")
        raise SystemExit(1)

    es = [got[N][1] for N in BAND]
    over = [e - TARGET for e in es]

    # -------------------------------------------------------------- W1
    say()
    say("W1  the local floor in this window")
    slope, flo = local_floor(BAND, es)
    perdec = slope * math.log(10.0)
    w1 = flo < halffloor
    say("  slope %+.6f per log unit (%+.6f per decade), r.m.s. "
        "residual %.6f" % (slope, perdec, flo))
    say("FLOOR target_band_local %.6f" % flo)
    say("  against %.6f at 1/2 and %.6f near rung 18"
        % (halffloor, farfloor))
    say("  W1 %s   (cap: the floor at 1/2)"
        % ("hold" if w1 else "REFUTED"))

    # -------------------------------------------------------------- W2
    say()
    say("W2  does the window straddle 0.56?")
    nbelow = sum(1 for v in over if v <= 0.0)
    nabove = len(over) - nbelow
    w2 = nbelow > 0 and nabove > 0
    say("  %d of %d below the level, %d above"
        % (nbelow, len(over), nabove))
    say("  W2 %s   (cap: one on each side)"
        % ("hold" if w2 else "REFUTED"))

    # -------------------------------------------------------------- W3
    say()
    say("W3  how many times does the sign turn over?")
    t = turns(over)
    w3 = t == 1
    say("  signs in N order: %s"
        % " ".join("+" if v > 0 else "-" for v in over))
    say("SIGNRUN target_band_turns %d" % t)
    say("  the sign changes %d time%s" % (t, "" if t == 1 else "s"))
    inside = sum(1 for v in over if abs(v) <= flo)
    say("  %d of %d points are inside the local floor, so their side "
        "is not determined" % (inside, len(over)))
    say("  W3 %s   (cap: exactly one turn)"
        % ("hold" if w3 else "REFUTED"))

    # -------------------------------------------------------------- W4
    say()
    say("W4  how wide is the indeterminate band here?")
    band = 2.0 * flo / abs(perdec) if perdec else float("inf")
    HALFBAND = 2.0 * halffloor / abs(halfslope * math.log(10.0))
    w4 = band < HALFBAND
    say("  2*floor/slope is %.4f decades, against %.4f at 1/2, "
        "whose slope" % (band, HALFBAND))
    say("  %+.6f per log unit and floor %.6f are read from "
        "results/audit_local_floor.txt" % (halfslope, halffloor))
    say("  the window spans %.4f decades, so it %s the band"
        % (math.log10(BAND[-1]) - math.log10(BAND[0]),
           "covers" if (math.log10(BAND[-1]) - math.log10(BAND[0]))
           > band else "sits inside"))
    say("SPREAD target_band %.4f"
        % (math.log10(BAND[-1]) - math.log10(BAND[0])))
    say("SCATTER slope_audit_target_band %.6f" % flo)
    say("  W4 %s   (cap: the band at 1/2)"
        % ("hold" if w4 else "REFUTED"))

    # ---------------------------------------------- the resolved pair
    say()
    say("the resolved pair for 0.56, from this window")
    below = [N for N, v in zip(BAND, over) if -v > flo]
    above = [N for N, v in zip(BAND, over) if v > flo]
    if below and above:
        b, a = max(below), min(above)
        say("  highest N resolved below: %d at 10^%.4f, exponent "
            "%.6f" % (b, math.log10(b), got[b][1]))
        say("  lowest N resolved above:  %d at 10^%.4f, exponent "
            "%.6f" % (a, math.log10(a), got[a][1]))
        say("BRACKET target_056_resolved %.4f %.4f %.4f"
            % ((math.log10(b) + math.log10(a)) / 2.0,
               math.log10(b), math.log10(a)))
        say("DRIFT target_056_resolved %.4f"
            % abs(math.log10(a) - math.log10(b)))
        say("SHAPES 1")
    else:
        say("  no pair in this window is resolved on both sides")
    say("SCALES 2")

    say()
    say("=" * 70)
    say("C1 %s  C2 %s  W1 %s  W2 %s  W3 %s  W4 %s"
        % tuple("hold" if v else "REFUTED"
                for v in (c1, c2, w1, w2, w3, w4)))

    head = [
        "STATISTIC: the truncation K*_R and its exponent",
        "           log K*_R / log N at seven same-radical N spanning",
        "           0.0641 decades about the place the chord between",
        "           two resolved measurements puts 0.56; the number",
        "           of times the sign of exponent - 0.56 turns over",
        "           across them; the r.m.s. residual of the seven",
        "           about a line in log N as a local floor; and the",
        "           width 2*floor/slope of the band in which a point",
        "           cannot be told from the level. Rung 16 is",
        "           recomputed as a control.",
        "NULL: none is run and none applies. K*_R is a deterministic",
        "      integer once N is fixed; there is no sampling noise",
        "      and no background to detect against. What the floor",
        "      measures is arithmetic fluctuation across N, which is",
        "      why it is estimated from neighbouring N rather than",
        "      from a fit across the whole ladder. The coin arms for",
        "      this statistic were run in lab_primorial_ladder.py",
        "      and lab_primorial_share.py.",
        "FIELD: N = 30030*m for m = 174960, 179200, 183708, 188160,",
        "       193050, 198000, 202800, every m composed only of 2,",
        "       3, 5, 7, 11 and 13 so the odd radical is 3*5*7*11*13",
        "       and the threshold is the ladder's; and N = 30030*2^16",
        "       as the control. k squarefree and coprime to N with",
        "       2 <= k < 1000000, beta fitted on the same range; m",
        "       odd, squarefree and coprime to k, m < N/k; the sieve",
        "       weight over the odd primes below 30; the Euler",
        "       products at the fixed bound 4000000. One odd",
        "       radical, as the RADICALS line declares. The sieve,",
        "       the packing and the statistic are imported from",
        "       code/audit_primorial_rung18.py; rung 16 comes from",
        "       results/audit_ladder_cap.txt, the floor at 1/2 from",
        "       results/audit_local_floor.txt, the floor near rung 18",
        "       from results/audit_rung18_fill.txt and rung 17's",
        "       exponent from results/audit_primorial_rung17.txt.",
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
