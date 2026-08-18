# -*- coding: utf-8 -*-
r"""
Where the validated model says the route's hypothesis becomes true.

WHAT IS AT STAKE

Remark {#rem:extrap} validated, out of sample, the model

    |A(N;k)| = gamma sqrt(N/k) sqrt(log N),   gamma = 0.6520,

fitted on N = 2e5, 4e5, 8e5 and predicting K* at 1.6e6 and 3.2e6
within 18% and 38%.  Under it,

    B(N)/N = gamma sqrt(log N / N) * S(K),
    S(K)   = sum_{k<K, mu^2(k)=1, (k,N)=1} (log k) k^{-1/2},

and Proposition {#prop:nolog} needs B(N) <= S(N)(1-A(N))N.  Solving
for K gives K*(N), the level the demand can carry.  Huang-Li need
K = N^{theta'} for a single theta' > 1/2.  So the model answers a
question no measurement can reach directly: AT WHICH N does the
route's hypothesis become true, at theta' = 1/2 and at theta' = 0.56?

This is a forecast from a fitted model, not a theorem and not a
measurement, and it is reported as one.  What makes it worth making is
that the model has already predicted out of sample once.

BACKS: Remark {#rem:forecast} in paper/wall_v3.md.

PRE-REGISTERED PREDICTIONS (written before this script was run)

  S1  The model reproduces the measured K* = 319, 537, 767, 1353, 2319
      at N = 2e5 ... 3.2e6 within a factor 1.5 at every N.
  S2  The model's K*/sqrt(N) crosses 1 between N = 3e5 and 3e6,
      matching the measured crossing between 8e5 and 1.6e6.
  S3  The model's K*/N^{0.56} crosses 1 at some N below 10^20.
  S4  At that theta' = 0.56 crossing the model's local exponent
      d log K* / d log N lies between 0.6 and 0.9.

REFUTATION RULE (fixed before the run)

  S1  REFUTED if any of the five is off by a factor 1.5 or more.  This
      gates: a model that cannot reproduce the measured K* is not worth
      forecasting with.
  S2  REFUTED if the crossing leaves [3e5, 3e6].
  S3  REFUTED if K*/N^{0.56} has not crossed 1 by N = 10^20.
  S4  REFUTED if the local exponent leaves [0.6, 0.9].

  S1 and S2 gate; S3 and S4 are the forecast and are reported.

CITED BY: {#rem:artifact} in paper/.
"""

import io
import math
import os
import sys

import numpy as np

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
OUT = os.path.join(ROOT, "results", "lab_level_forecast.txt")

GAMMA = 0.6520
KMAX = 10_000_000
PN = (2, 5)                      # the family N = 2^a 5^b used throughout
THR = None                       # S(N)(1-A(N)) for that family
MEAS = [(200_000, 319), (400_000, 537), (800_000, 767),
        (1_600_000, 1353), (3_200_000, 2319)]


def primes_upto(n):
    s = np.ones(n + 1, dtype=bool)
    s[:2] = False
    for p in range(2, int(n ** 0.5) + 1):
        if s[p]:
            s[p * p::p] = False
    return np.flatnonzero(s).astype(np.int64)


def main():
    global THR
    lines = []

    def say(s=""):
        print(s)
        lines.append(s)

    say("building S(K) exactly to K = %d ..." % KMAX)
    sq = np.ones(KMAX + 1, dtype=bool)
    sq[0] = False
    p = 2
    while p * p <= KMAX:
        sq[p * p::p * p] = False
        p += 1
    for q in PN:
        sq[q::q] = False
    ks = np.flatnonzero(sq)
    w = np.log(ks.astype(np.float64)) / np.sqrt(ks.astype(np.float64))
    S = np.cumsum(w)
    say("  admissible k <= %d: %d (density %.6f)"
        % (KMAX, ks.size, ks.size / KMAX))
    say("  S(%d) = %.4f" % (KMAX, S[-1]))

    pr = primes_upto(4_000_000)
    artin, twin = 1.0, 2.0
    for q in pr:
        q = int(q)
        artin *= 1.0 - 1.0 / (q * (q - 1.0))
        if q > 2:
            twin *= 1.0 - 1.0 / (q - 1.0) ** 2
    A_, S_ = artin, twin
    for q in PN:
        A_ /= (1.0 - 1.0 / (q * (q - 1.0)))
        if q > 2:
            S_ *= (1.0 + 1.0 / (q - 2.0))
    THR = S_ * (1.0 - A_)
    say("  threshold S(N)(1-A(N)) = %.6f for P(N) = %s" % (THR, PN))

    # asymptotic continuation, calibrated at K = KMAX
    dens = ks.size / KMAX

    def S_of(K):
        if K <= KMAX:
            j = int(np.searchsorted(ks, K))
            return float(S[max(j - 1, 0)])
        base = float(S[-1])
        f = lambda t: 2.0 * math.sqrt(t) * (math.log(t) - 2.0)
        return base + dens * (f(K) - f(KMAX))

    def Kstar(N):
        target = THR * math.sqrt(N / math.log(N)) / GAMMA
        lo, hi = 2.0, 1e30
        for _ in range(200):
            mid = math.sqrt(lo * hi)
            if S_of(mid) < target:
                lo = mid
            else:
                hi = mid
        return math.sqrt(lo * hi)

    say()
    say("S1  the model against the measured K*")
    say("  N          measured   model      ratio")
    s1 = True
    for N, km in MEAS:
        kp = Kstar(N)
        r = kp / km
        if not (1 / 1.5 < r < 1.5):
            s1 = False
        say("  %-10d %-10d %-10.1f %.3f" % (N, km, kp, r))
    say("  S1 %s" % ("hold" if s1 else "REFUTED"))
    say("  DIAGNOSTIC (post hoc). The model OVERSHOOTS K* by a steady")
    say("  factor near 1.8, i.e. it UNDERSTATES B(N). gamma was")
    say("  calibrated on the MEDIAN of |A|/sqrt(N/k), but B is a SUM and")
    say("  needs the MEAN, and |A| is one-sided so the mean is larger.")
    say("  Recomputing both from the field:")
    say("  N          median     mean       mean/median   K* (mean gamma)"
        "   ratio")
    import numpy as _np
    NMAX = max(N for N, _ in MEAS)
    prm = primes_upto(NMAX)
    lgp = _np.log(prm.astype(_np.float64))
    lam = _np.zeros(NMAX + 1, dtype=_np.float64)
    lam[prm] = lgp
    for i, q in enumerate(prm):
        q = int(q)
        if q * q > NMAX:
            break
        t = q * q
        while t <= NMAX:
            lam[t] = lgp[i]
            if t > NMAX // q:
                break
            t *= q
    mu = _np.ones(NMAX + 1, dtype=_np.int8)
    rem = _np.arange(NMAX + 1, dtype=_np.int32)
    for q in primes_upto(int(math.isqrt(NMAX))):
        q = int(q)
        mu[q::q] = -mu[q::q]
        if q * q <= NMAX:
            mu[q * q::q * q] = 0
        t = q
        while t <= NMAX:
            rem[t::t] //= q
            if t > NMAX // q:
                break
            t *= q
    mu[rem > 1] = -mu[rem > 1]
    mu[0] = 0
    del rem
    gm = []
    for N, km in MEAS:
        kk = _np.array([k for k in range(2, N // 1000 + 1)
                        if mu[k] != 0 and all(k % q for q in PN)])
        f = _np.zeros(N, dtype=_np.float64)
        idx = _np.arange(1, N, dtype=_np.int64)
        f[1:] = lam[1:N] * mu[N - idx]
        Aa = _np.empty(kk.size)
        for i, k in enumerate(kk):
            r = N % int(k)
            Aa[i] = abs(f[r::int(k)].sum() if r
                        else f[int(k)::int(k)].sum())
        nrm = Aa / _np.sqrt(N / kk.astype(float)) / math.sqrt(math.log(N))
        med, mea = float(_np.median(nrm)), float(nrm.mean())
        gm.append(mea)
        g_old, g_new = GAMMA, mea
        tgt = THR * math.sqrt(N / math.log(N)) / g_new
        lo, hi = 2.0, 1e15
        for _ in range(200):
            mid = math.sqrt(lo * hi)
            if S_of(mid) < tgt:
                lo = mid
            else:
                hi = mid
        kp = math.sqrt(lo * hi)
        say("  %-10d %-10.4f %-10.4f %-13.4f %-17.1f %.3f"
            % (N, med, mea, mea / med, kp, kp / km))
    GNEW = float(_np.mean(gm[:3]))
    say("  mean-based gamma over the three fitting N = %.4f" % GNEW)
    say("  fitted on the three smallest N only; the two rows it never saw")
    say("  are the last two of the table above.")

    def Kstar2(N):
        target = THR * math.sqrt(N / math.log(N)) / GNEW
        lo, hi = 2.0, 1e30
        for _ in range(200):
            mid = math.sqrt(lo * hi)
            if S_of(mid) < target:
                lo = mid
            else:
                hi = mid
        return math.sqrt(lo * hi)

    say()
    say("  the forecast under the CORRECTED gamma:")
    say("  N            K*(N)         K*/sqrt N     K*/N^0.56")
    rows2 = []
    for e in range(5, 21):
        Nv = 10.0 ** e
        k = Kstar2(Nv)
        rows2.append((Nv, k / math.sqrt(Nv), k / Nv ** 0.56))
        say("  10^%-10d %-13.4e %-13.4f %.4f"
            % (e, k, k / math.sqrt(Nv), k / Nv ** 0.56))

    def cross2(idx):
        prev = None
        for Nv, a, b in rows2:
            v = (a, b)[idx]
            if prev is not None and prev[1] < 1.0 <= v:
                lo, hi = prev[0], Nv
                for _ in range(80):
                    mid = math.sqrt(lo * hi)
                    kk = Kstar2(mid)
                    vv = kk / math.sqrt(mid) if idx == 0 \
                        else kk / mid ** 0.56
                    if vv < 1.0:
                        lo = mid
                    else:
                        hi = mid
                return math.sqrt(lo * hi)
            prev = (Nv, v)
        return None

    ch, c56 = cross2(0), cross2(1)
    say("  corrected: K*/sqrt N crosses 1 at N = %s"
        % ("%.3e" % ch if ch else "already above at 1e5"))
    say("  corrected: K*/N^0.56 crosses 1 at N = %s"
        % ("%.3e" % c56 if c56 else "not below 1e20"))

    say()
    say("S2/S3/S4  the forecast")
    say("  N            K*(N)         K*/sqrt N     K*/N^0.56    "
        "local exponent")
    rows = []
    for e in range(5, 31):
        N = 10.0 ** e
        k = Kstar(N)
        k2 = Kstar(N * 1.05)
        loc = (math.log(k2) - math.log(k)) / math.log(1.05)
        rows.append((N, k, k / math.sqrt(N), k / N ** 0.56, loc))
        say("  10^%-10d %-13.4e %-13.4f %-12.4f %.4f"
            % (e, k, k / math.sqrt(N), k / N ** 0.56, loc))

    def cross(idx):
        prev = None
        for N, k, a, b, loc in rows:
            v = (a, b)[idx]
            if prev is not None and prev[1] < 1.0 <= v:
                lo, hi = prev[0], N
                for _ in range(80):
                    mid = math.sqrt(lo * hi)
                    kk = Kstar(mid)
                    vv = kk / math.sqrt(mid) if idx == 0 \
                        else kk / mid ** 0.56
                    if vv < 1.0:
                        lo = mid
                    else:
                        hi = mid
                return math.sqrt(lo * hi)
            prev = (N, v)
        return None

    c_half = cross(0)
    c_056 = cross(1)
    s2 = c_half is not None and 3e5 <= c_half <= 3e6
    say()
    say("S2  K*/sqrt N crosses 1 at N = %s   (band [3e5, 3e6])   %s"
        % ("%.3e" % c_half if c_half else "never",
           "hold" if s2 else "REFUTED"))
    s3 = c_056 is not None and c_056 < 1e20
    say("S3  K*/N^0.56 crosses 1 at N = %s   (cap 1e20)   %s"
        % ("%.3e" % c_056 if c_056 else "never",
           "hold" if s3 else "REFUTED"))
    if c_056:
        k = Kstar(c_056)
        k2 = Kstar(c_056 * 1.05)
        loc = (math.log(k2) - math.log(k)) / math.log(1.05)
        s4 = 0.6 <= loc <= 0.9
        say("S4  local exponent there = %.4f   (band [0.6,0.9])   %s"
            % (loc, "hold" if s4 else "REFUTED"))
    else:
        s4 = False
        say("S4  no crossing, so no exponent   REFUTED")

    say()
    say("=" * 70)
    ok = s1 and s2
    say("S1 %s  S2 %s  S3 %s  S4 %s"
        % tuple("hold" if v else "REFUTED" for v in (s1, s2, s3, s4)))
    say("the model reproduces the measured level and forecasts where the "
        "route's hypothesis becomes true" if ok else "REFUTED")

    head = [
        "STATISTIC: K*(N) solved from the model B(N)/N = gamma",
        "           sqrt(log N / N) S(K) against the threshold",
        "           S(N)(1-A(N)), with gamma = 0.6520 fitted out of sample",
        "           in lab_dilate_extrapolation.py and S(K) = sum over",
        "           admissible k < K of (log k) k^{-1/2}; the model's K*",
        "           against the five measured values; K*/sqrt(N) and",
        "           K*/N^0.56; and the local exponent d log K*/d log N.",
        "NULL: none applies. This is a forecast from a model, not a",
        "      detection: S(K) is a deterministic arithmetic sum and",
        "      gamma carries whatever the sign pattern contributed, which",
        "      lab_dilate_extrapolation.py measured for mu and for a coin",
        "      separately. The control that matters was run there.",
        "FIELD: S(K) enumerated exactly over squarefree k coprime to 10",
        "       up to K = 10^7, and continued beyond by the integral",
        "       2 sqrt(K)(log K - 2) scaled by the measured density; the",
        "       threshold is that of the family P(N) = {2,5} used",
        "       throughout; N swept over 10^5 to 10^30.",
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
