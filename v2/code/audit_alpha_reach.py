# -*- coding: utf-8 -*-
r"""
The sign axis, asked four octaves further out.

WHAT IS AT STAKE

This is the axis the proof is stuck on.  {#rem:leanidentity} left one
requirement: |sum a| must reach the order of its own l2 norm, a gap of
+0.134019 in exponent.  Every constructive route to it is closed --
{#rem:levelmagnitude} in principle, {#rem:filter} because the inverse
filter grows like 1.916413^m, {#rem:meanonly} because the denominator
is a correlation and not a mean.  What is left is one computational
question: **does the gap close as N grows?**

{#rem:alphalocal} answered it "no, nothing is closing", and the
answer is honest but weak, and weak in a way that matters.  The
deficit's slope against log N is +0.014678 with a standard error of
0.017249, t = 0.85.  A drift two standard errors inside that -- 0.027
per log unit -- would close a gap of 0.134019 in about five log units,
two and a bit decades.  **So the measurement as it stands cannot
exclude the budget route; it only fails to see it.**  Its sign is the
wrong way for the route besides: the deficit is drifting up, not
down, at a size the errors cannot resolve.

The field stops at N = 1.024e8, nine octaves.  That bound is an
artefact of the sieve it was measured with: {#rem:alphalocal} uses
the float64 route of code/audit_gain_split.py, which holds Lambda at
eight bytes per index and an int32 cofactor beside it.  The packing
built for the level axis in {#rem:rung18} stores one byte per index
and no prime at all, and it has carried a sieve to 8.26e9.  **The
sign axis has never been given it.**

So the field is extended to 8e9 -- fifteen octaves instead of nine,
the same reach the level axis now has -- and the two slopes are
refitted.  The half-index trick of
{#rem:rung17} does not apply here: m runs over every integer coprime
to k, not the odd ones, so N - mk takes both parities and the whole
index is addressed.  What does apply is the kind byte, which is what
lifts the memory from eight bytes an index to two.

BACKS: Remark {#rem:alphareach} in paper/wall_v3.md.

PRE-REGISTERED PREDICTIONS (written before this script was run)

  R1  The octave control.  Recomputed through this route, the top
      published octave returns the alpha_local and deficit
      results/audit_alpha_local.txt prints for octave 8, +0.664374
      and +0.108604, to the six decimals printed.
  R2  The field control.  The extended field has one coprimality
      class, as the published COPRIME line says of the shorter one.
  R3  The alpha drift stays unresolved: over the extended set of
      octaves the slope of alpha_local against mid log N has
      |t| < 2.
  R4  **The deficit grows, resolved.**  Its slope over the extended
      set is positive with |t| > 2.  The budget route would then be
      closed by measurement rather than by failure to see.
  R5  And the extension buys the resolution it was run for: the
      deficit slope's standard error falls below 0.010, from the
      0.017249 nine octaves gave.  The lever arm grows by about 1.7
      and the octave count by two thirds, so 0.010 is a real test
      rather than an arithmetic certainty.

REFUTATION RULE (fixed before the run)

  R1  REFUTED if either differs in the six decimals printed.  The
      sums here are split into blocks where the published route took
      them whole, so this is the check that the association does not
      move the exponent.  THIS ONE GATES.
  R2  REFUTED by more than one class; the points would not be one
      field and the octaves would not be comparable.  THIS ONE
      GATES.
  R3  REFUTED if |t| reaches 2.  A resolved alpha drift would be a
      result in its own right and would have to be read together
      with e(l2)'s, since only their difference is the requirement.
  R4  REFUTED two ways, and they mean opposite things, so both are
      named here as M9 requires.  **(a) The slope resolves
      negative**: the deficit is closing, the budget route is alive,
      and the branch would have a number for how far N must be
      pushed -- this is the outcome that would matter most and it is
      not the one predicted.  **(b) The slope stays unresolved**:
      nothing is learned about the route except that four more
      octaves were not enough, and the honest report is that the
      exclusion is still weak.  "Too noisy to tell" is not "the
      deficit is flat".
  R5  REFUTED at or above 0.010.  Then the extension did not buy
      what it cost, and reaching further is the only remaining move
      on this axis.

  R1 and R2 gate.  R3 to R5 are the measurement and do not gate.

  NO NULL IS RUN and none applies.  |sum a| and its l2 norm are
  deterministic once N is fixed; there is no sampling noise and no
  background to detect against.  The coin arms for the sign axis were
  run in lab_primorial_share.py and are what established that the
  quantity is about mu and not about magnitudes.
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
OUT = os.path.join(RES, "audit_alpha_reach.txt")

BLOCK = 1 << 24
DEC = 6
NEWTOP = 8_000_000_000              # the extended bound


def module(name):
    p = os.path.join(CODE, name + ".py")
    spec = importlib.util.spec_from_file_location(name, p)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


SPL = module("audit_gain_split")
AL = module("audit_alpha_local")
THETA = SPL.THETA
LO, HI = AL.LO, AL.HI


def family(lo, hi):
    return sorted(set(
        2 ** a * 5 ** b
        for a in range(1, 64) for b in range(1, 64)
        if lo <= 2 ** a * 5 ** b <= hi))


def power_table(n):
    """the prime powers p^j, j >= 2, below n, and their log p"""
    vs, lg = [], []
    for p in SPL.primes_upto(int(math.isqrt(n))):
        p = int(p)
        lp = math.log(float(p))
        q = p * p
        while q <= n:
            vs.append(q)
            lg.append(lp)
            if q > n // p:
                break
            q *= p
    o = np.argsort(np.array(vs, dtype=np.int64))
    return (np.array(vs, dtype=np.int64)[o],
            np.array(lg, dtype=np.float64)[o])


def kind_and_mu(n, block=BLOCK):
    """the kind byte and Moebius on whole indices

    kind[v] is 1 when v is prime, 2 when it is a higher power of a
    prime, 0 otherwise, so Lambda(v) is log v, log p from the power
    table, or zero.  Two bytes an index against the eight the
    published route holds for Lambda alone, and no prime is stored,
    so nothing here is capped at 2^32.
    """
    root = int(math.isqrt(n))
    pr = [int(p) for p in SPL.primes_upto(root)]
    kind = np.zeros(n + 1, dtype=np.uint8)
    mu = np.empty(n + 1, dtype=np.int8)
    for lo in range(0, n + 1, block):
        hi = min(lo + block, n + 1)
        w = hi - lo
        vals = np.arange(lo, hi, dtype=np.int64)
        rem = vals.copy()
        m = np.ones(w, dtype=np.int8)
        for p in pr:
            s = (-lo) % p
            if s < w:
                m[s::p] = -m[s::p]
                rem[s::p] //= p
            q = p * p
            s = (-lo) % q
            if s < w:
                m[s::q] = 0
            pk = p * p
            while pk <= n:
                s = (-lo) % pk
                if s < w:
                    rem[s::pk] //= p
                if pk > n // p:
                    break
                pk *= p
        big = rem > 1
        m[big] = -m[big]
        mu[lo:hi] = m
        idx = np.flatnonzero(big & (rem == vals))
        if idx.size:
            kind[lo + idx] = 1
        del vals, rem, m, big, idx
    mu[0] = 0
    for p in pr:
        kind[p] = 1
    pv, _ = power_table(n)
    kind[pv] = 2
    return kind, mu


def weighted(N, kind, mu, pv, plg, block=BLOCK):
    """(log k) H(N;k) over squarefree k < N^theta coprime to N

    The set of k and of m is the published route's; only the sum over
    m is split into blocks, and Lambda is read off the kind byte with
    the logarithm taken at the point of use of the same float64
    integer.  R1 is what says neither moves the exponent.
    """
    PN = SPL.factor_set(N)
    K = int(N ** THETA)
    ks, Hs = [], []
    for k in range(2, K):
        if mu[k] == 0:
            continue
        if any(k % q == 0 for q in PN):
            continue
        M = (N - 1) // k
        if M < 1:
            continue
        drop = sorted(SPL.factor_set(k))
        h = 0.0
        for lo in range(1, M + 1, block):
            hi = min(lo + block, M + 1)
            ms = np.arange(lo, hi, dtype=np.int64)
            for q in drop:
                ms = ms[ms % q != 0]
            if ms.size == 0:
                continue
            vals = N - ms * k
            g = mu[ms].astype(np.float64)
            kd = kind[vals]
            pri = kd == 1
            if pri.any():
                h += float((np.log(vals[pri].astype(np.float64))
                            * g[pri]).sum())
            pwr = kd == 2
            if pwr.any():
                w = vals[pwr]
                h += float((plg[np.searchsorted(pv, w)] * g[pwr]).sum())
        ks.append(k)
        Hs.append(h)
    ks = np.array(ks, dtype=np.int64)
    return ks, np.log(ks.astype(np.float64)) * np.array(Hs)


def fit(x, y):
    a, b = np.polyfit(x, y, 1)
    r = y - (a * x + b)
    n = x.size
    se = math.sqrt(float((r ** 2).sum() / (n - 2))
                   / float(((x - x.mean()) ** 2).sum()))
    return float(a), float(se)


def octave_marker(fname, j):
    src = io.open(os.path.join(RES, fname), encoding="utf-8").read()
    m = re.search(r"^OCTAVE alphalocal_%d ([\d.]+) ([-+][\d.]+) "
                  r"([\d.]+) ([-+][\d.]+) ([\d.]+)\s*$" % j, src, re.M)
    if not m:
        raise SystemExit("no OCTAVE alphalocal_%d line" % j)
    return tuple(float(g) for g in m.groups())


def main():
    lines = []

    def say(t=""):
        print(t)
        sys.stdout.flush()
        lines.append(t)

    pub = [octave_marker("audit_alpha_local.txt", j) for j in range(9)]
    say("the nine published octaves, each read from a whole marker "
        "line:")
    for j, r in enumerate(pub):
        say("READ audit_alpha_local.txt OCTAVE alphalocal_%d %.4f "
            "%+.6f %.6f %+.6f %.6f" % ((j,) + r))
    rnd = 0.5 * 10.0 ** (-DEC)
    say("PRINTBOUND audit_alpha_reach %d %.8f" % (DEC, rnd))

    NS_old = family(LO, HI)
    NS_new = [N for N in family(LO, NEWTOP) if N > HI]
    # the published binning, continued: octave j is
    # [LO*2^j, LO*2^(j+1)), the last one closed at its top
    ctrl = [N for N in NS_old if LO * 2 ** 8 <= N <= HI]
    if len(ctrl) < 3:
        raise SystemExit("the control octave is empty; the binning "
                         "does not match code/audit_alpha_local.py")
    say()
    say("the published field is %d points in [%d, %d]; the extension "
        "adds %d" % (len(NS_old), LO, HI, len(NS_new)))
    say("  to %d, %.1f octaves against %.1f"
        % (NEWTOP, math.log2(max(NS_new) / float(LO)),
           math.log2(HI / float(LO))))
    say("  the control octave holds %d of the published points"
        % len(ctrl))
    classes = sorted(set(tuple(sorted(SPL.factor_set(N)))
                         for N in NS_old + NS_new))
    say("RADICALS %d" % len(classes))
    say("COPRIME %d" % len(classes))
    r2 = len(classes) == 1
    say("  R2 %s   (cap: one class)" % ("hold" if r2 else "REFUTED"))

    TOP = max(NS_new)
    say()
    say("sieving to %d on whole indices with the kind byte" % TOP)
    kind, mu = kind_and_mu(TOP)
    pv, plg = power_table(TOP)
    say("BYTES resident_arrays %d" % (kind.nbytes + mu.nbytes))
    say("BYTES published_route %d" % (8 * (TOP + 1) + (TOP + 1)
                                      + 4 * (TOP + 1)))
    say("  two arrays of %.2f GB, against %.2f GB for the float64 "
        "route of" % ((kind.nbytes + mu.nbytes) / 2.0 ** 30,
                      (13 * (TOP + 1)) / 2.0 ** 30))
    say("  code/audit_gain_split.py at this N; the power table holds "
        "%d entries" % pv.size)

    def measure(NS, tag):
        S, L = [], []
        for N in NS:
            _ks, a = weighted(N, kind, mu, pv, plg)
            S.append(abs(float(a.sum())))
            L.append(float(np.sqrt((a ** 2).sum())))
            say("POINT alphareach_%d %.10e %.10e"
                % (N, S[-1], L[-1]))
            say("  %-8s N = %-12d |sum a| = %.6e   l2 = %.6e"
                % (tag, N, S[-1], L[-1]))
        return np.array(S), np.array(L)

    # -------------------------------------------------------------- R1
    say()
    say("R1  the control: the top published octave through this route")
    Sc, Lc = measure(ctrl, "ctrl")
    xc = np.log(np.array(ctrl, dtype=np.float64))
    ac, sac = fit(xc, np.log(Sc))
    lc, slc = fit(xc, np.log(Lc))
    dc = ac - lc
    r1 = (abs(ac - pub[8][1]) <= rnd and abs(dc - pub[8][3]) <= rnd)
    say("  alpha_local here %+.6f against the published %+.6f"
        % (ac, pub[8][1]))
    say("  deficit here %+.6f against the published %+.6f"
        % (dc, pub[8][3]))
    say("  R1 %s   (cap: the printing bound)"
        % ("hold" if r1 else "REFUTED"))
    if not (r1 and r2):
        io.open(OUT, "w", encoding="utf-8", newline="\n").write(
            "\n".join(lines) + "\n")
        raise SystemExit(1)

    # ------------------------------------------------ the new octaves
    say()
    say("the extension")
    Sn, Ln = measure(NS_new, "new")
    xn = np.log(np.array(NS_new, dtype=np.float64))
    newocts = []
    j = 9
    while LO * 2 ** j < NEWTOP:
        lo, hi = LO * 2 ** j, LO * 2 ** (j + 1)
        sel = np.array([(lo <= N < hi) or (N == min(hi, NEWTOP))
                        for N in NS_new])
        j += 1
        if sel.sum() < 3:
            continue
        xo = xn[sel]
        a1, s1 = fit(xo, np.log(Sn[sel]))
        a2, s2 = fit(xo, np.log(Ln[sel]))
        dd = a1 - a2
        sd = math.sqrt(s1 * s1 + s2 * s2)
        newocts.append((float(xo.mean()), a1, s1, dd, sd,
                        int(sel.sum())))
        say("OCTAVE alphareach_%d %.4f %+.6f %.6f %+.6f %.6f"
            % (j - 1, xo.mean(), a1, s1, dd, sd))
        say("  octave %-3d %2d points  mid log N %.4f  alpha %+.6f "
            "(%.6f)  deficit %+.6f (%.6f)"
            % (j - 1, sel.sum(), xo.mean(), a1, s1, dd, sd))
    say("POP alphareach_octave %d"
        % (min(o[5] for o in newocts) if newocts else 0))

    # --------------------------------------------------------- R3, R4
    npts = [0] * len(pub) + [o[5] for o in newocts]
    mids = np.array([r[0] for r in pub] + [o[0] for o in newocts])
    alph = np.array([r[1] for r in pub] + [o[1] for o in newocts])
    defi = np.array([r[3] for r in pub] + [o[3] for o in newocts])
    m9 = np.array([r[0] for r in pub])
    sa9, sea9 = fit(m9, np.array([r[1] for r in pub]))
    sd9, sed9 = fit(m9, np.array([r[3] for r in pub]))
    sa, sea = fit(mids, alph)
    sd_, sed = fit(mids, defi)
    ta, td = sa / sea, sd_ / sed
    r3 = abs(ta) < 2.0
    r4 = sd_ > 0.0 and abs(td) > 2.0
    say()
    say("R3/R4  the two slopes over %d octaves" % mids.size)
    say("  alpha   slope %+.6f +- %.6f, t = %.2f" % (sa, sea, ta))
    say("TSTAT alphareach_alpha_slope %.2f" % ta)
    if abs(ta) < 2.0:
        say("UNRESOLVED SIGN alphareach_alpha_slope")
    say("  deficit slope %+.6f +- %.6f, t = %.2f" % (sd_, sed, td))
    say("TSTAT alphareach_deficit_slope %.2f" % td)
    if abs(td) < 2.0:
        say("UNRESOLVED SIGN alphareach_deficit_slope")
    say("SPREAD alphareach_alpha_slope %.4f"
        % (mids.max() - mids.min()))
    say("SPREAD alphareach_deficit_slope %.4f"
        % (mids.max() - mids.min()))
    say("SCATTER slope_audit_alpha_reach %.6f"
        % float(np.sqrt(((defi - (sd_ * mids + (defi.mean()
                                                - sd_ * mids.mean())))
                         ** 2).mean())))
    say("  the nine published octaves alone give %+.6f +- %.6f "
        "for alpha" % (sa9, sea9))
    say("  and %+.6f +- %.6f for the deficit, refitted here from "
        "their marker lines" % (sd9, sed9))
    swa, swd = [], []
    for i in range(mids.size):
        keep = np.ones(mids.size, dtype=bool)
        keep[i] = False
        swa.append(fit(mids[keep], alph[keep])[0])
        swd.append(fit(mids[keep], defi[keep])[0])
    swa, swd = np.array(swa), np.array(swd)
    say("  dropping one octave at a time moves the alpha slope over "
        "%+.6f to %+.6f" % (swa.min(), swa.max()))
    say("SWEPT alphareach_alpha_slope leave-one-octave %.6f"
        % (swa.max() - swa.min()))
    say("  and the deficit slope over %+.6f to %+.6f"
        % (swd.min(), swd.max()))
    say("SWEPT alphareach_deficit_slope leave-one-octave %.6f"
        % (swd.max() - swd.min()))
    worst = int(np.argmax(np.abs(swa - sa)))
    say("  the octave that moves it most is the one at mid log N "
        "%.4f, which holds" % mids[worst])
    say("  %d points; the field was cut at %d inside an octave, so "
        "the last one is partial"
        % (npts[worst], NEWTOP))
    say("  R3 %s   (cap: |t| = 2)" % ("hold" if r3 else "REFUTED"))
    say("  R4 %s   (cap: positive and |t| > 2)"
        % ("hold" if r4 else "REFUTED"))

    # -------------------------------------------------------------- R5
    say()
    say("R5  did the extension buy resolution?")
    r5 = sed < 0.010
    say("  the deficit slope's standard error is %.6f, against "
        "%.6f on nine octaves" % (sed, sed9))
    say("  R5 %s   (cap: 0.010)" % ("hold" if r5 else "REFUTED"))

    say()
    say("what would close the gap")
    gap = 0.134019
    if sd_ < 0:
        say("  at %+.6f per log unit the deficit %.6f would close in "
            "%.1f log units" % (sd_, gap, gap / abs(sd_)))
        say("  which is %.1f decades past the top of this field"
            % (gap / abs(sd_) / math.log(10.0)))
    else:
        say("  the slope is not negative, so no closure follows from "
            "it at any N")
    say("  the exclusion this field can support: a drift smaller "
        "than %.6f" % (2.0 * sed))
    say("  is invisible here, and such a drift would close the gap "
        "in %.1f decades"
        % (gap / (2.0 * sed) / math.log(10.0)))

    say()
    say("=" * 70)
    say("R1 %s  R2 %s  R3 %s  R4 %s  R5 %s"
        % tuple("hold" if v else "REFUTED"
                for v in (r1, r2, r3, r4, r5)))

    head = [
        "STATISTIC: the exponents of |sum a| and of its l2 norm",
        "           against N, fitted inside each octave of the field",
        "           N = 2^a 5^b extended from 1.024e8 to 8e9; their",
        "           difference, the deficit that must reach zero; and",
        "           the slopes of both against mid log N over the",
        "           octaves that result, against the nine",
        "           results/audit_alpha_local.txt reports. The top",
        "           published octave is recomputed as a control.",
        "NULL: none is run and none applies. |sum a| and its l2 norm",
        "      are deterministic once N is fixed; there is no",
        "      sampling noise and no background to detect against.",
        "      The coin arms for the sign axis were run in",
        "      lab_primorial_share.py, which established that the",
        "      quantity is about mu and not about magnitudes.",
        "FIELD: N = 2^a 5^b with both a >= 1 and b >= 1 in",
        "       [2e5, 8e9], one coprimality class as COPRIME says",
        "       -- the class {2,5}, k coprime to 10 and N even; k",
        "       squarefree with 2 <= k < N^0.56; m over 1 <= m < N/k",
        "       with (m,k) = 1. Lambda and mu come from a kind-byte",
        "       sieve on whole indices, two bytes per index, and the",
        "       per-k sums are split into blocks; the k-set, the",
        "       m-set and theta' are code/audit_gain_split.py's, and",
        "       R1 is the check that neither the packing nor the",
        "       splitting moves the exponent. The nine published",
        "       octaves are read from marker lines in",
        "       results/audit_alpha_local.txt and declared with READ.",
        "",
    ]
    io.open(OUT, "w", encoding="utf-8", newline="\n").write(
        "\n".join(head + lines) + "\n")
    print("\nwrote %s" % OUT)
    if not (r1 and r2):
        raise SystemExit(1)
    return 0


if __name__ == "__main__":
    sys.exit(main())
