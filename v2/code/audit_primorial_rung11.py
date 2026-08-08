# -*- coding: utf-8 -*-
r"""
One more rung: is the crossing of 1/2 observed, or inside its noise?

WHAT IS AT STAKE

Remark {#rem:primorialrung10} reports the primorial ladder's level
exponent reaching 0.5023 at N = 30030*2^10 = 30750720 and reads that
as the square-root barrier crossed at the hard arithmetic -- observed,
not extrapolated. It is one rung, and the ladder's own scatter about
its fitted line is 0.0037 r.m.s. **A clearance of 0.0023 is less than
that scatter.** By the standard this project adopted in gate checks
G37 and G38 -- a claim judged against a threshold declares the noise
floor of the same statistic -- the rung-10 crossing has never been
put against its own floor, and it does not survive that comparison on
its own.

What settles it is the next rung, N = 30030*2^11 = 61501440. The line
fitted on eleven rungs predicts 0.5010 + 0.006780 = 0.5078 there,
which is 0.0078 above 1/2 -- twice the scatter. If the exponent lands
there the crossing is two independent points and a margin outside the
floor; if it falls back to 1/2 the rung-10 reading was a fluctuation
and has to be withdrawn.

The implementation is the bitmask one of audit_level_slope_reach.py,
independent of audit_primorial_rung10.py's, so P1 doubles as a
cross-check of the published rung.

BACKS: Remark {#rem:primorialrung11} in paper/theorem_A.md.

PRE-REGISTERED PREDICTIONS (written before this script was run)

  P1  The control: this independent implementation reproduces the
      published rung-10 exponent to within 0.001.
  P2  The barrier stays crossed: at N = 61501440 the exponent is
      above 1/2.
  P3  And the crossing is now outside its own noise: the margin over
      1/2 exceeds the ladder's published r.m.s. residual, which
      rung 10's 0.0023 did not.
  P4  The ladder is still a line: rung 11 lies within one published
      r.m.s. residual of the line fitted on the eleven published
      rungs.

REFUTATION RULE (fixed before the run)

  P1  REFUTED at 0.001, which would mean this is not the same
      statistic and nothing below may be compared with the ladder.
  P2  REFUTED if the exponent is at or below 1/2. That is the one
      that matters: {#rem:primorialrung10}'s G3 would have to be
      withdrawn and the crossing would go back to being a forecast.
  P3  REFUTED if the margin is at most the published r.m.s. residual.
      The crossing would then be real as a fitted trend but not
      observable at any single rung, and the paper would have to say
      so.
  P4  REFUTED if the residual exceeds the published r.m.s., which
      would say the ladder bends where it matters and the fitted
      crossing of theta' means nothing.

  All four gate.

  NO NULL IS RUN and none applies. A deterministic curve is located
  against a computed threshold; there is no background to detect
  against. The coin arms for this statistic were run in
  lab_primorial_ladder.py and lab_primorial_share.py, which
  established that both the rise and its scatter are facts about
  magnitudes and not about mu, and the scatter they left is exactly
  the floor P3 is judged against.
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
OUT = os.path.join(RES, "audit_primorial_rung11.txt")

BASE = 30030                        # 2*3*5*7*11*13
CONTROL = BASE * (1 << 10)          # 30750720, the top published rung
NEW = BASE * (1 << 11)              # 61501440
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
    """the eleven published rungs and the ladder's scatter"""
    src = io.open(os.path.join(RES, "audit_primorial_rung10.txt"),
                  encoding="utf-8").read()
    i = src.index("N            log10 N   exponent   fitted     residual")
    ns, ex = [], []
    for ln in src[i:].splitlines()[1:]:
        f = ln.split()
        if len(f) < 3 or not f[0].isdigit():
            break
        ns.append(int(f[0]))
        ex.append(float(f[2]))
    m = re.search(r"the scatter from [\d.]+ to ([\d.]+)", src)
    return ns, ex, float(m.group(1))


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

    lns, lex, floor = read_ladder()
    say("read %d published rungs and the ladder scatter %.4f from "
        "results/audit_primorial_rung10.txt" % (len(lns), floor))
    pub10 = dict(zip(lns, lex))[CONTROL]

    qs = [int(q) for q in primes_upto(QSIEVE) if q > 2]
    say("sieving to %d, sieve weight over the odd primes %s"
        % (NEW, ", ".join(map(str, qs))))
    lam, mu = lambda_and_mu(NEW)
    sqf = mu != 0
    vmask = residue_mask(NEW, qs)

    artin, twin = 1.0, 2.0
    for p in primes_upto(CLIM):
        p = int(p)
        artin *= 1.0 - 1.0 / (p * (p - 1.0))
        if p > 2:
            twin *= 1.0 - 1.0 / (p - 1.0) ** 2

    got = {}
    for N in (CONTROL, NEW):
        r = measure(N, lam, mu, sqf, vmask, qs, artin, twin)
        if r is None:
            say("  N = %-12d no crossing below k = %d" % (N, KCAP))
            continue
        kstar, e, bpn, beta, nk = r
        got[N] = e
        say("  N = %-12d thr %.6f  #k %-7d beta %.6f  K*_R %-8d "
            "exp %.4f" % (N, bpn, nk, beta, kstar, e))
        say("BUDGET kstar_R_S1AN_N%d %.6f" % (N, bpn))
    rads = set(tuple(sorted(q for q in factor_set(N) if q > 2))
               for N in got)
    say("RADICALS %d" % len(rads))

    # ------------------------------------------------------------- P1
    say()
    say("P1  the control at N = %d" % CONTROL)
    d1 = abs(got[CONTROL] - pub10)
    p1 = d1 < 0.001
    say("  %.4f here against the published %.4f, diff %.6f"
        % (got[CONTROL], pub10, d1))
    say("  P1 %s   (cap 0.001)" % ("hold" if p1 else "REFUTED"))

    # ---------------------------------------------------------- P2/P3
    say()
    say("P2/P3  the new rung against 1/2 and against the floor")
    e11 = got[NEW]
    marg = e11 - 0.5
    p2 = e11 > 0.5
    p3 = marg > floor
    say("  N            log10 N   exponent   margin over 1/2  floor")
    say("  %-12d %-9.4f %-10.4f %-16.4f %.4f"
        % (NEW, math.log10(NEW), e11, marg, floor))
    say("  and the rung this extends, for comparison:")
    say("  %-12d %-9.4f %-10.4f %-16.4f %.4f"
        % (CONTROL, math.log10(CONTROL), got[CONTROL],
           got[CONTROL] - 0.5, floor))
    say("FLOOR primorial_rung11 %.4f" % floor)
    for stem, mg in (("audit_primorial_rung10", got[CONTROL] - 0.5),
                     ("audit_primorial_rung11", marg)):
        say("MARGIN %s %.4f %.4f" % (stem, mg, floor))
        if mg <= floor:
            say("INSIDE FLOOR %s" % stem)
    say("  P2 the exponent is above 1/2 (%.4f)   %s"
        % (e11, "hold" if p2 else "REFUTED"))
    say("  P3 the margin %.4f exceeds the floor %.4f   %s"
        % (marg, floor, "hold" if p3 else "REFUTED"))

    # ------------------------------------------------------------- P4
    say()
    say("P4  is the ladder still a line?")
    x = np.log(np.array(lns, dtype=np.float64))
    y = np.array(lex)
    x0 = math.log(NEW)
    a, b = np.polyfit(x, y, 1)
    pred = a * x0 + b
    resid = e11 - pred
    p4 = abs(resid) <= floor
    say("  line on the %d published rungs, against log N: slope %+.6f"
        % (len(lns), a))
    say("  at N = %d it predicts %.4f, measured %.4f, off by %+.4f"
        % (NEW, pred, e11, resid))
    say("  published r.m.s. residual %.4f   %s"
        % (floor, "within" if p4 else "OUTSIDE"))
    say("  P4 %s" % ("hold" if p4 else "REFUTED"))

    say()
    say("  DIAGNOSTIC on P4 (post hoc). The rule compared an")
    say("  out-of-sample point with an in-sample r.m.s., and those")
    say("  are not the same width. A least-squares line predicts a")
    say("  new abscissa with standard error")
    say("  s*sqrt(1 + 1/n + (x0-xbar)^2/sum (x-xbar)^2), s^2 being")
    say("  the residual variance on n-2 degrees of freedom, and x0")
    say("  here is outside the fitted range, which inflates it")
    say("  further. That width is computable:")
    r = y - (a * x + b)
    n = x.size
    s2 = float((r ** 2).sum() / (n - 2))
    sxx = float(((x - x.mean()) ** 2).sum())
    sepred = math.sqrt(s2 * (1.0 + 1.0 / n + (x0 - x.mean()) ** 2 / sxx))
    say("  in-sample r.m.s.        %.4f" % floor)
    say("  prediction s.e. at x0   %.4f" % sepred)
    say("  measured departure      %.4f   %s"
        % (abs(resid),
           "inside" if abs(resid) <= sepred else "OUTSIDE"))
    say("  So P4 fails as registered and the ladder is not thereby")
    say("  shown to bend: the departure is %.2f prediction standard"
        % (abs(resid) / sepred))
    say("  errors. What P4 should have been compared against is this")
    say("  width, and the twelve-rung refit below is the check that")
    say("  matters -- if the ladder bent, the scatter would rise.")

    say()
    say("  the twelve rungs refitted, which is what any later")
    say("  extrapolation must use:")
    x12 = np.append(x, x0)
    y12 = np.append(y, e11)
    a12, b12 = np.polyfit(x12, y12, 1)
    r12 = y12 - (a12 * x12 + b12)
    rms12 = float(np.sqrt((r12 ** 2).mean()))
    se12 = math.sqrt(float((r12 ** 2).sum() / (x12.size - 2))
                     / float(((x12 - x12.mean()) ** 2).sum()))
    say("  slope %+.6f, r.m.s. residual %.4f, standard error %.6f, "
        "t = %.2f" % (a12, rms12, se12, abs(a12) / se12))
    say("SCATTER slope_audit_primorial_rung11 %.4f" % rms12)
    say("TSTAT slope_audit_primorial_rung11 %.2f" % (abs(a12) / se12))
    say("SPREAD slope_audit_primorial_rung11 %.4f"
        % float(x12.max() - x12.min()))
    if abs(a12) / se12 < 2.0:
        say("UNRESOLVED SIGN slope_audit_primorial_rung11")
    say("  no crossing is forecast from this: the barrier is now")
    say("  behind the ladder rather than ahead of it, and where the")
    say("  line meets theta' is the shape question Remark")
    say("  {#rem:laddershape} left open and this rung does not touch.")

    say()
    say("=" * 70)
    ok = p1 and p2 and p3 and p4
    say("the crossing is observed at two rungs and outside its floor"
        if ok else "REFUTED")

    head = [
        "STATISTIC: the truncation K*_R at which",
        "           sum_{k<K}(log k)|R(N;k)| first reaches",
        "           S(N)(1-A(N))N, and its exponent log K*_R / log N,",
        "           at N = 30030*2^11 = 61501440 and, as a control, at",
        "           N = 30030*2^10; the margin of the new exponent",
        "           over 1/2 against the ladder's own r.m.s. residual;",
        "           and the residual of the new rung from the line",
        "           fitted on the eleven published rungs.",
        "NULL: none is run and none applies. A deterministic curve is",
        "      located against a computed threshold; there is no",
        "      background to detect against. The coin arms for this",
        "      statistic were run in lab_primorial_ladder.py and",
        "      lab_primorial_share.py, and the scatter they left is",
        "      the floor the margin is judged against here.",
        "FIELD: N = 30750720 and 61501440, the odd radical 3*5*7*11*13",
        "       fixed so the threshold is constant; k squarefree and",
        "       coprime to N with 2 <= k < " + str(KCAP) + "; m odd,",
        "       squarefree and coprime to k, m < N/k; the sieve weight",
        "       over the odd primes below " + str(QSIEVE) + "; the",
        "       Euler products at the fixed bound " + str(CLIM) + ".",
        "       One odd radical, as the RADICALS line declares. The",
        "       eleven published rungs and the scatter are read from",
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
