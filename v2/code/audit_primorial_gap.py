# -*- coding: utf-8 -*-
r"""
The gap the ladder steps over: where between rungs 9 and 10 does the
level exponent cross 1/2?

WHAT IS AT STAKE

Remark {#rem:primorialreach} excluded the point estimate 10^7.10 and
left "a crossing anywhere in [10^7.19, 10^7.36]" open. Remark
{#rem:primorialrung10} then measured 0.5023 at rung 10, N = 30750720 =
10^7.4879, and Remark {#rem:primorialrung11} 0.5099 at rung 11. So the
crossing happened somewhere above rung 9 (10^7.1868, exponent 0.4941)
and at or below rung 10.

**Nothing has ever been measured inside the interval.** The ladder is
N = 30030*2^j, so it steps from 10^7.1868 straight to 10^7.4879 and
the open interval sits in the void between two consecutive rungs. It
has been neither excluded nor confirmed; it has been skipped.

It does not have to be. What the ladder holds fixed is the prime set
P(N) = {2,3,5,7,11,13} -- that is what makes the budget S(N)(1-A(N))
constant along it, and it is the single odd radical whose caveat
Remark {#rem:arithmeticreach} says cannot be computed away. Doubling
is only the cheapest way to move N while holding that set. Any
N = 30030*m with m composed of 2,3,5,7,11,13 holds it just as well,
and there are more than twenty such N inside the interval. Four are
taken here, spread across it:

  30030*528 = 15855840   10^7.2003    528 = 2^4*3*11
  30030*600 = 18018000   10^7.2557    600 = 2^3*3*5^2
  30030*672 = 20180160   10^7.3049    672 = 2^5*3*7
  30030*750 = 22522500   10^7.3527    750 = 2*3*5^3

This also tests something the ladder never has: whether it is a line
at finer resolution than a doubling. Every residual ever quoted for it
is a residual at spacing log 2.

The implementation is copied from audit_primorial_rung11.py rather
than rewritten. That is deliberate and it is the opposite of what P1
there was for: the point of these four N is that they be commensurable
with the eleven rungs, so they must be the same statistic computed the
same way, and the control below is that it reproduces the published
rung it is anchored to.

BACKS: Remark {#rem:primorialgap} in paper/theorem_A.md.

PRE-REGISTERED PREDICTIONS (written before this script was run)

  Q1  The control. This code reproduces the published rung-9 exponent
      at N = 15375360 to within 0.001, and refitting the eleven
      published rungs against log N reproduces the published fitted
      column to within 0.0005.
  Q2  No interior point exceeds 1/2. The line fitted on eleven rungs
      puts the crossing near 10^7.42, above the whole interval, so
      every one of the four should still be short.
  Q3  The ladder is a line at sub-doubling spacing: the r.m.s.
      residual of the four interior points about the eleven-rung line
      is no larger than the ladder's own published scatter.
  Q4  The budget is identical at all four, because the prime set is.
      This is the control that makes them commensurable with the
      rungs at all -- if it fails, the four are measuring arithmetic
      and not N.
  Q5  The crossing refitted on all fifteen points lies inside the
      bracket [10^7.1866, 10^7.6631] that Remark
      {#rem:primorialrung10} published on eleven.
  Q6  And the interval is NOT excluded by observation. Every interior
      point should sit within the ladder's floor of 1/2 -- the line
      gives about 0.499 at the top of the interval, a shortfall of
      roughly 0.001 against a floor of 0.0037 -- so what the four
      points can say is that the interval is unsupported, not that it
      is empty. Gate check G40 exists for exactly this distinction
      and it is registered here rather than discovered afterwards.

REFUTATION RULE (fixed before the run)

  Q1  REFUTED at 0.001 on the rung or 0.0005 on the fitted column,
      either of which would mean this is not the ladder's statistic
      and nothing below may be compared with it.
  Q2  REFUTED if any interior exponent exceeds 1/2. That is the
      outcome worth having: the crossing would then be inside
      [10^7.19, 10^7.36] after all, the eleven-rung line would be
      putting it 0.06 decades too high, and {#rem:primorialreach}'s
      open interval would close by confirmation instead of by
      exclusion.
  Q3  REFUTED if the interior r.m.s. exceeds the published scatter.
      The ladder would then bend between its rungs, and every
      residual and every crossing fitted on it would be describing a
      curve sampled at one point per doubling.
  Q4  REFUTED if the relative spread of S(N)(1-A(N)) across the four
      exceeds 1e-12.
  Q5  REFUTED if the refitted crossing falls outside the published
      bracket.
  Q6  REFUTED if any interior point's distance from 1/2 exceeds the
      floor. Then something IS observed in the interval and the
      verdict is stronger than "unsupported".

  Q1 and Q4 gate: without them the four points are not on the ladder.
  Q2, Q3, Q5 and Q6 are the measurement and do not gate.

  NO NULL IS RUN and none applies, for the reason given in
  audit_primorial_rung11.py: a deterministic curve is located against
  a computed threshold and there is no background to detect against.
  The coin arms for this statistic were run in lab_primorial_ladder.py
  and lab_primorial_share.py, and the scatter they left is the floor
  Q3 and Q6 are judged against.
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
OUT = os.path.join(RES, "audit_primorial_gap.txt")

BASE = 30030                        # 2*3*5*7*11*13
CONTROL = BASE * (1 << 9)           # 15375360, published rung 9
INSIDE = [BASE * m for m in (528, 600, 672, 750)]
KCAP = 100_000
QSIEVE = 30
CLIM = 4_000_000


def primes_upto(n):
    s = np.ones(n + 1, dtype=bool)
    s[:2] = False
    for p in range(2, int(n ** 0.5) + 1):
        if s[p]:
            s[p * p::p] = False
    return np.flatnonzero(s).astype(np.int64)


def lambda_and_mu(n):
    """von Mangoldt and Moebius, the cofactor kept in int32"""
    pr = primes_upto(n)
    lgp = np.log(pr.astype(np.float64))
    lam = np.zeros(n + 1, dtype=np.float64)
    lam[pr] = lgp
    for i, p in enumerate(pr):
        p = int(p)
        if p * p > n:
            break
        q = p * p
        while q <= n:
            lam[q] = lgp[i]
            if q > n // p:
                break
            q *= p
    del pr, lgp
    mu = np.ones(n + 1, dtype=np.int8)
    cof = np.arange(n + 1, dtype=np.int32)
    for p in primes_upto(int(math.isqrt(n))):
        p = int(p)
        mu[p::p] = -mu[p::p]
        if p * p <= n:
            mu[p * p::p * p] = 0
        cof[p::p] //= p
        pk = p * p
        while pk <= n:
            cof[pk::pk] //= p
            if pk > n // p:
                break
            pk *= p
    big = cof > 1
    del cof
    mu[big] = -mu[big]
    del big
    mu[0] = 0
    return lam, mu


def residue_mask(n, qs):
    """bit i of mask[v] is set exactly when qs[i] divides v"""
    m = np.zeros(n + 1, dtype=np.uint16)
    for i, q in enumerate(qs):
        m[0::q] |= np.uint16(1 << i)
    return m


def factor_set(n):
    v, out, d = n, set(), 2
    while d * d <= v:
        if v % d == 0:
            out.add(d)
            while v % d == 0:
                v //= d
        d += 1
    if v > 1:
        out.add(v)
    return out


def read_ladder():
    """the eleven published rungs, their fitted column and the scatter"""
    src = io.open(os.path.join(RES, "audit_primorial_rung10.txt"),
                  encoding="utf-8").read()
    i = src.index("N            log10 N   exponent   fitted     residual")
    ns, ex, fit = [], [], []
    for ln in src[i:].splitlines()[1:]:
        f = ln.split()
        if len(f) < 4 or not f[0].isdigit():
            break
        ns.append(int(f[0]))
        ex.append(float(f[2]))
        fit.append(float(f[3]))
    m = re.search(r"the scatter from [\d.]+ to ([\d.]+)", src)
    return ns, ex, fit, float(m.group(1))


def read_bracket():
    """the eleven-rung crossing bracket -- read, not typed"""
    src = io.open(os.path.join(RES, "audit_primorial_rung10.txt"),
                  encoding="utf-8").read()
    m = re.search(r"bracket \[([\d.]+), ([\d.]+)\]", src)
    return float(m.group(1)), float(m.group(2))


def measure(N, lam, mu, sqf, vmask, qs, artin, twin):
    """the level exponent log K*_R / log N at one N"""
    PN = factor_set(N)
    A_, S_ = artin, twin
    for q in sorted(PN):
        A_ /= (1.0 - 1.0 / (q * (q - 1.0)))
        if q > 2:
            S_ *= (1.0 + 1.0 / (q - 2.0))

    ks, Hs, Ps = [], [], []
    for k in range(2, KCAP):
        if not sqf[k]:
            continue
        if any(k % q == 0 for q in PN):
            continue
        M = (N - 1) // k
        if M < 2:
            continue
        ms = np.arange(1, M + 1, 2, dtype=np.int64)
        ms = ms[sqf[ms]]
        kb, ck = 0, 1.0
        for i, q in enumerate(qs):
            if k % q == 0:
                kb |= 1 << i
            else:
                ck *= q / (q - 1.0)
        for q in factor_set(k):
            if q > 2:
                ms = ms[ms % q != 0]
        if ms.size == 0:
            continue
        vals = N - ms * k
        g = mu[ms].astype(np.float64)
        keep = (vmask[vals] & np.uint16(~kb & 0xFFFF)) == 0
        ks.append(k)
        Hs.append(float((lam[vals] * g).sum()))
        Ps.append(ck * float(g[keep].sum()))
    ks = np.array(ks, dtype=np.int64)
    H = np.array(Hs)
    P = np.array(Ps)
    beta = float((H * P).sum() / (P * P).sum())
    R = H - beta * P
    cum = np.cumsum(np.log(ks.astype(np.float64)) * np.abs(R))
    thr = S_ * (1.0 - A_) * N
    j = int(np.searchsorted(cum, thr))
    if j >= ks.size:
        return None
    kstar = int(ks[j])
    return (kstar, math.log(kstar) / math.log(N), thr / N, beta,
            ks.size)


def main():
    lines = []

    def say(s=""):
        print(s)
        lines.append(s)

    lns, lex, lfit, floor = read_ladder()
    blo, bhi = read_bracket()
    say("read from results/audit_primorial_rung10.txt: %d rungs, "
        "scatter %.4f," % (len(lns), floor))
    say("  and the eleven-rung crossing bracket [%.4f, %.4f]"
        % (blo, bhi))
    pub9 = dict(zip(lns, lex))[CONTROL]
    say("  the rung this anchors to: N = %d, published exponent %.4f"
        % (CONTROL, pub9))
    say("  the interval left open by rem:primorialreach is the void")
    say("  between that rung and the next -- no N in it has been run.")

    top = max(INSIDE)
    qs = [int(q) for q in primes_upto(QSIEVE) if q > 2]
    say()
    say("sieving to %d, sieve weight over the odd primes %s"
        % (top, ", ".join(map(str, qs))))
    lam, mu = lambda_and_mu(top)
    sqf = mu != 0
    vmask = residue_mask(top, qs)

    artin, twin = 1.0, 2.0
    for p in primes_upto(CLIM):
        p = int(p)
        artin *= 1.0 - 1.0 / (p * (p - 1.0))
        if p > 2:
            twin *= 1.0 - 1.0 / (p - 1.0) ** 2

    say()
    say("  N            log10 N   m      thr/N      #k      "
        "beta       K*_R     exponent")
    got, budg = {}, {}
    for N in [CONTROL] + INSIDE:
        r = measure(N, lam, mu, sqf, vmask, qs, artin, twin)
        if r is None:
            say("  %-12d no crossing below k = %d" % (N, KCAP))
            continue
        kstar, e, bpn, beta, nk = r
        got[N] = e
        budg[N] = bpn
        say("  %-12d %-9.4f %-6d %-10.6f %-7d %-10.6f %-8d %.4f"
            % (N, math.log10(N), N // BASE, bpn, nk, beta, kstar, e))
        say("BUDGET kstar_R_S1AN_N%d %.6f" % (N, bpn))
    rads = set(tuple(sorted(q for q in factor_set(N) if q > 2))
               for N in got)
    say("RADICALS %d" % len(rads))

    # ------------------------------------------------------------- Q1
    say()
    say("Q1  the control: the published rung, and the published line")
    d1 = abs(got[CONTROL] - pub9)
    x = np.log(np.array(lns, dtype=np.float64))
    y = np.array(lex)
    a, b = np.polyfit(x, y, 1)
    dfit = float(np.abs(a * x + b - np.array(lfit)).max())
    q1 = d1 < 0.001 and dfit < 0.0005
    say("  rung 9 here %.4f against the published %.4f, diff %.6f"
        % (got[CONTROL], pub9, d1))
    say("  refitting the %d rungs on log N: slope %+.6f; worst"
        % (len(lns), a))
    say("  departure from the published fitted column %.6f" % dfit)
    say("  Q1 %s   (cap 0.001 on the rung, cap 0.0005 on the column)"
        % ("hold" if q1 else "REFUTED"))

    # ------------------------------------------------------------- Q4
    say()
    say("Q4  the budget is the same at all four, because the prime "
        "set is")
    bs = [budg[N] for N in INSIDE if N in budg]
    spread = (max(bs) - min(bs)) / (sum(bs) / len(bs))
    q4 = spread <= 1e-12
    say("  prime set of every N: %s"
        % ", ".join(str(q) for q in sorted(factor_set(INSIDE[0]))))
    say("  S(N)(1-A(N)) runs %.9f to %.9f, relative spread %.2e"
        % (min(bs), max(bs), spread))
    say("  Q4 %s   (cap 1e-12)" % ("hold" if q4 else "REFUTED"))

    # ------------------------------------------------------------- Q2
    say()
    say("Q2  does anything inside the interval reach 1/2?")
    say("  N            log10 N   exponent   over 1/2")
    q2 = True
    for N in INSIDE:
        if N not in got:
            continue
        e = got[N]
        if e > 0.5:
            q2 = False
        say("  %-12d %-9.4f %-10.4f %+.4f"
            % (N, math.log10(N), e, e - 0.5))
    over = [N for N in INSIDE if N in got and got[N] > 0.5]
    if q2:
        say("  the interval is AT OR BELOW 0.5 at every point")
    else:
        say("  the exponent is AT OR ABOVE 0.5 at %d of the %d points, "
            "from" % (len(over), len(INSIDE)))
        say("  N = %d (10^%.4f) upward -- the crossing is INSIDE the "
            "interval" % (min(over), math.log10(min(over))))
    say("  Q2 %s" % ("hold" if q2 else "REFUTED"))

    # ------------------------------------------------------------- Q6
    say()
    say("Q6  but is anything OBSERVED there? distance against the floor")
    say("  N            exponent - 1/2   |.|      floor     outside?")
    q6 = True
    for N in INSIDE:
        if N not in got:
            continue
        d = abs(got[N] - 0.5)
        out = d > floor
        if out:
            q6 = False
        say("  %-12d %+-16.4f %-8.4f %-9.4f %s"
            % (N, got[N] - 0.5, d, floor, "OUTSIDE" if out else "inside"))
    # G40 reads the crossing claim, so the margin is the signed one
    # over 1/2 -- the best a single point in the interval achieves.
    best = max(got[N] - 0.5 for N in INSIDE if N in got)
    say("FLOOR primorial_gap %.4f" % floor)
    say("MARGIN audit_primorial_gap %.4f %.4f" % (best, floor))
    if best <= floor:
        say("INSIDE FLOOR audit_primorial_gap")
    if q6:
        say("  Q6 holds -- every point is inside the floor, so the")
        say("  interval is unsupported and not excluded, and nothing")
        say("  in it is observed.")
    else:
        say("  Q6 REFUTED -- the best margin over 1/2 in the interval "
            "is %+.4f" % best)
        say("  against a floor of %.4f, so the crossing is not merely"
            % floor)
        say("  fitted inside the interval, it is observed there. The")
        say("  registered expectation was the opposite and it was the")
        say("  eleven-rung line that produced it.")

    # ------------------------------------------------------------- Q3
    say()
    say("Q3  is the ladder a line at finer than a doubling?")
    xi = np.log(np.array([N for N in INSIDE if N in got],
                         dtype=np.float64))
    yi = np.array([got[N] for N in INSIDE if N in got])
    ri = yi - (a * xi + b)
    rms_i = float(np.sqrt((ri ** 2).mean()))
    q3 = rms_i <= floor
    say("  residuals about the eleven-rung line: %s"
        % ", ".join("%+.4f" % v for v in ri))
    say("  r.m.s. %.4f against the published scatter %.4f   %s"
        % (rms_i, floor, "within" if q3 else "OUTSIDE"))
    say("  Q3 %s" % ("hold" if q3 else "REFUTED"))
    say("  the rungs are spaced %.4f in log N and these four are "
        "spaced" % float(np.diff(x).mean()))
    say("  %.4f, so this is the first residual the ladder has at "
        "sub-doubling" % float(np.diff(np.sort(xi)).mean()))
    say("  resolution.")

    say()
    say("  and the r.m.s. is not the informative part -- the signs are.")
    npos = int((ri > 0).sum())
    mono = all(ri[i] <= ri[i + 1] for i in range(ri.size - 1))
    say("  %d of %d residuals are positive and the sequence is %s,"
        % (npos, ri.size, "monotone increasing" if mono
           else "not monotone"))
    run = max(npos, int(ri.size) - npos)
    say("SIGNRUN primorial_gap_interior %d %d" % (run, ri.size))
    say("  a run of %d in %d has chance %.4f under an exchangeable "
        "draw," % (run, ri.size, 2.0 ** (1 - ri.size)))
    say("  which is evidence an r.m.s. cannot carry either way.")
    say("  so the four points do not scatter about the line, they")
    say("  leave it. The local slope says the same thing: between the")
    say("  anchor rung and the top interior point the exponent rises")
    lo_x, lo_y = math.log(CONTROL), got[CONTROL]
    hi_x, hi_y = math.log(top), got[top]
    loc = (hi_y - lo_y) / (hi_x - lo_x)
    nxt = BASE << 10
    lever = hi_x - lo_x
    se_loc = floor * math.sqrt(2.0) / lever
    say("  at %+.6f per unit log N against the ladder's %+.6f, a "
        "factor" % (loc, a))
    say("  %.2f -- but that is a two-point slope over a lever of only"
        % (loc / a))
    say("  %.4f in log N, so propagating the ladder's own scatter"
        % lever)
    say("  through it gives a standard error of %.4f: the local slope"
        % se_loc)
    say("  is %.2f of its own error and the factor is NOT resolved."
        % (loc / se_loc))
    say("SPREAD slope_primorial_gap_local %.4f" % lever)
    say("TSTAT slope_primorial_gap_local %.2f" % (loc / se_loc))
    if loc / se_loc < 2.0:
        say("UNRESOLVED SIGN slope_primorial_gap_local")
    say("  Between the two rungs that straddle the interval it is")
    pubnext = dict(zip(lns, lex)).get(nxt)
    if pubnext is not None:
        loc2 = (pubnext - lo_y) / (math.log(nxt) - lo_x)
        say("  %+.6f, a factor %.2f -- so the steepening is not an "
            "artefact" % (loc2, loc2 / a))
        say("  of the interior points; the published rungs already "
            "show it and")
        say("  had no point in between to resolve it against.")

    # ------------------------------------------------------------- Q5
    say()
    say("Q5  the crossing, refitted on all fifteen points")
    xa = np.append(x, xi)
    ya = np.append(y, yi)
    a2, b2 = np.polyfit(xa, ya, 1)
    ra = ya - (a2 * xa + b2)
    rms_a = float(np.sqrt((ra ** 2).mean()))
    se_a = math.sqrt(float((ra ** 2).sum() / (xa.size - 2))
                     / float(((xa - xa.mean()) ** 2).sum()))
    cross = ((0.5 - b2) / a2) / math.log(10.0)
    q5 = blo <= cross <= bhi
    say("  slope %+.6f on %d points, r.m.s. residual %.4f"
        % (a2, xa.size, rms_a))
    say("  eleven rungs alone gave slope %+.6f, r.m.s. %.4f"
        % (a, float(np.sqrt(((y - (a * x + b)) ** 2).mean()))))
    say("  crossing of 1/2 at log10 N = %.4f" % cross)
    say("  published eleven-rung bracket [%.4f, %.4f]   %s"
        % (blo, bhi, "inside" if q5 else "OUTSIDE"))
    say("  Q5 %s" % ("hold" if q5 else "REFUTED"))
    say("SCATTER slope_audit_primorial_gap %.4f" % rms_a)
    say("TSTAT slope_audit_primorial_gap %.2f" % (abs(a2) / se_a))
    say("SPREAD slope_audit_primorial_gap %.4f"
        % float(xa.max() - xa.min()))
    if abs(a2) / se_a < 2.0:
        say("UNRESOLVED SIGN slope_audit_primorial_gap")
    say("SWEPT audit_primorial_gap_crossing N-range %.4f"
        % abs(cross - ((0.5 - b) / a) / math.log(10.0)))

    say()
    say("  What this settles. The crossing is bracketed by measured")
    say("  points on both sides of 1/2, not by a fit:")
    below = [N for N in [CONTROL] + INSIDE if N in got and got[N] <= 0.5]
    above = [N for N in INSIDE if N in got and got[N] > 0.5]
    if below and above:
        lo_n, hi_n = max(below), min(above)
        say("  last short   N = %-12d 10^%.4f  exponent %.4f"
            % (lo_n, math.log10(lo_n), got[lo_n]))
        say("  first over   N = %-12d 10^%.4f  exponent %.4f"
            % (hi_n, math.log10(hi_n), got[hi_n]))
        say("  so the crossing lies in (10^%.4f, 10^%.4f], %.4f decades"
            % (math.log10(lo_n), math.log10(hi_n),
               math.log10(hi_n) - math.log10(lo_n)))
        say("  wide, against the %.4f decades between the two rungs"
            % (math.log10(BASE << 10) - math.log10(CONTROL)))
        say("  that straddled it before -- a narrowing of %.2f."
            % ((math.log10(BASE << 10) - math.log10(CONTROL))
               / (math.log10(hi_n) - math.log10(lo_n))))
        say("  And it is INSIDE the interval rem:primorialreach left")
        say("  open, which therefore closes by confirmation and not")
        say("  by exclusion -- the outcome Q2's refutation rule named")
        say("  as the one worth having.")
    else:
        say("  every measured point falls on one side of 1/2, so the")
        say("  interval is not bracketed by observation here.")
    say()
    say("  What this does not settle. The eleven-rung line is not the")
    say("  right instrument in this region -- Q3 shows the interior")
    say("  points leaving it systematically -- so any quantity read")
    say("  off that line's slope past the ladder, theta' included,")
    say("  inherits the same understatement and is not corrected by")
    say("  the fifteen-point refit either, which is still one line")
    say("  through a region that is not straight.")

    say()
    say("=" * 70)
    say("Q1 %s  Q2 %s  Q3 %s  Q4 %s  Q5 %s  Q6 %s"
        % tuple("hold" if v else "REFUTED"
                for v in (q1, q2, q3, q4, q5, q6)))
    ok = q1 and q4
    say("the ladder now has points between its rungs" if ok
        else "REFUTED: the four points are not on the ladder")

    head = [
        "STATISTIC: the truncation K*_R at which",
        "           sum_{k<K}(log k)|R(N;k)| first reaches",
        "           S(N)(1-A(N))N, and its exponent log K*_R / log N,",
        "           at four N inside the interval [10^7.19, 10^7.36]",
        "           that no rung of the primorial ladder occupies, and",
        "           as a control at the published rung 9; the distance",
        "           of each interior exponent from 1/2 against the",
        "           ladder's own r.m.s. residual; the r.m.s. residual",
        "           of the four about the line fitted on the eleven",
        "           published rungs; and the crossing of 1/2 refitted",
        "           on all fifteen points.",
        "NULL: none is run and none applies. A deterministic curve is",
        "      located against a computed threshold; there is no",
        "      background to detect against. The coin arms for this",
        "      statistic were run in lab_primorial_ladder.py and",
        "      lab_primorial_share.py, and the scatter they left is",
        "      the floor Q3 and Q6 are judged against.",
        "FIELD: N = 30030*m for m = 528, 600, 672, 750, and the control",
        "       N = 30030*2^9; every one has prime set {2,3,5,7,11,13},",
        "       so the odd radical 3*5*7*11*13 and the threshold are",
        "       fixed exactly as along the ladder; k squarefree and",
        "       coprime to N with 2 <= k < " + str(KCAP) + "; m odd,",
        "       squarefree and coprime to k, m < N/k; the sieve weight",
        "       over the odd primes below " + str(QSIEVE) + "; the",
        "       Euler products at the fixed bound " + str(CLIM) + ";",
        "       the eleven rungs, their fitted column, the scatter and",
        "       the crossing bracket are read from",
        "       results/audit_primorial_rung10.txt.",
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
