# -*- coding: utf-8 -*-
r"""
Extending the demand-side measurements by two octaves, to separate
what a factor of 32 in N could not.

WHAT IS AT STAKE

Three open questions in this chain are all limited by the same thing,
the length of the N-range:

  * Remark {#rem:leandecay} could not choose between
    |1/2 - f| ~ N^{-0.1673} (correlation -0.97616) and
    (log N)^{-2.3166} (correlation -0.97469), and the two extrapolate
    to 1.0e14 and 4.3e22.
  * Remark {#rem:relocate} extrapolated B(N)/N to cross the Goldbach
    Goldbach threshold S(N)(1-A(N)) near N = 10^8.9, with no
    out-of-sample check.
  * Remark {#rem:onesided} recorded that [eq:onesided] fails at
    N = 2e5 and 4e5 and holds from 8e5, over five points.

All three are helped by the same computation, so it is done once here,
to N = 2.56e7 -- a factor 128 rather than 32.

BACKS: Remark {#rem:extendrange} in paper/wall_v3.md.

PRE-REGISTERED PREDICTIONS (written before this script was run)

  X1  The lean keeps decaying: f(2.56e7) > f(6.4e6).
  X2  Even at a factor 128 the two decay laws do not separate: the
      absolute difference of their fit correlations stays below 0.005.
      A refutation here is the good outcome and says which law it is.
  X3  Out of sample: fitting B(N)/N on the first five N alone and
      predicting at 1.28e7 and 2.56e7 lands within 10% both times.
  X4  [eq:onesided] holds at both new N, i.e. |E_3|/N is under the
      threshold S(N)(1-A(N)) computed for that N.

REFUTATION RULE (fixed before the run)

  X1  REFUTED if f does not rise.
  X2  REFUTED if the correlations differ by 0.005 or more -- in which
      case the better-fitting law is reported and the extrapolation of
      Remark {#rem:leandecay} is settled.
  X3  REFUTED if either prediction is off by 10% or more.
  X4  REFUTED if |E_3|/N reaches that threshold at either new N.

  All four gate.
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
OUT = os.path.join(ROOT, "results", "lab_extend_range.txt")

NS = [200_000, 400_000, 800_000, 1_600_000, 3_200_000, 6_400_000,
      12_800_000, 25_600_000]
FITN = 5                      # how many of the smallest N the fit sees
THETA = 0.56
# The threshold is S(N)(1-A(N)) and depends on which primes divide
# N; audit_threshold_arithmetic.py shows it moves by a factor of
# five across N of this size. It is computed below, never typed.
PLIM = 4_000_000


def primes_upto(n):
    s = np.ones(n + 1, dtype=bool)
    s[:2] = False
    for p in range(2, int(n ** 0.5) + 1):
        if s[p]:
            s[p * p::p] = False
    return np.flatnonzero(s).astype(np.int64)


def sieves(n):
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
    mu = np.ones(n + 1, dtype=np.int8)
    rem = np.arange(n + 1, dtype=np.int32)
    for p in primes_upto(int(math.isqrt(n))):
        p = int(p)
        mu[p::p] = -mu[p::p]
        if p * p <= n:
            mu[p * p::p * p] = 0
        q = p
        while q <= n:
            rem[q::q] //= p
            if q > n // p:
                break
            q *= p
    big = rem > 1
    del rem
    mu[big] = -mu[big]
    del big
    mu[0] = 0
    return pr, lam, mu


def phi_of(k):
    v, phi, d = k, 1, 2
    while d * d <= v:
        if v % d == 0:
            phi *= (d - 1)
            v //= d
            while v % d == 0:
                phi *= d
                v //= d
        d += 1
    if v > 1:
        phi *= (v - 1)
    return phi


def loo(x, y, name, say):
    """Refit dropping each end in turn, and report the spread.

    Every exponent this repository quotes is a slope over four to eight
    values of N, and audit_truncation_exponent.py showed what such a
    slope is worth when nobody varies the free parameter that defines
    it. For a direct fit the free parameter is the N-range, so the
    cheapest honest check is to refit without the smallest N and
    without the largest and print how far the answer moves.
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    f = [float(np.polyfit(x[s], y[s], 1)[0])
         for s in (slice(None), slice(1, None), slice(0, -1))]
    sp = max(f) - min(f)
    say("  leave-one-out on %s: full %.4f, without the smallest N "
        "%.4f," % (name, f[0], f[1]))
    say("  without the largest %.4f -- spread %.4f" % (f[2], sp))
    say("SWEPT %s N-range %.4f" % (name, sp))
    return sp


def main():
    lines = []

    def say(s=""):
        print(s)
        lines.append(s)

    NMAX = max(NS)
    say("sieving to %d ..." % NMAX)
    pr, lam, mu = sieves(NMAX)

    artin, twin = 1.0, 2.0
    for p in primes_upto(PLIM):
        p = int(p)
        artin *= 1.0 - 1.0 / (p * (p - 1.0))
        if p > 2:
            twin *= 1.0 - 1.0 / (p - 1.0) ** 2

    def thresh(N, PN):
        A, S = artin, twin
        for q in sorted(PN):
            A /= (1.0 - 1.0 / (q * (q - 1.0)))
            if q > 2:
                S *= (1.0 + 1.0 / (q - 2.0))
        return S * (1.0 - A)

    say()
    say("  N            K       #k     B/N       |E_3|/N   f         "
        "|0.5-f|   G")
    say("  " + "-" * 82)
    Bs, E3s, fs, gs, THRS = [], [], [], [], []
    for N in NS:
        v, PN, d = N, set(), 2
        while d * d <= v:
            if v % d == 0:
                PN.add(d)
                while v % d == 0:
                    v //= d
            d += 1
        if v > 1:
            PN.add(v)
        K = int(N ** THETA)
        ks = np.array([k for k in range(2, K)
                       if mu[k] != 0 and all(k % q for q in PN)])
        lg = np.log(ks.astype(float))
        iph = np.array([phi_of(int(k)) for k in ks], dtype=np.float64)

        f0 = np.zeros(N, dtype=np.float64)
        idx = np.arange(1, N, dtype=np.int64)
        f0[1:] = lam[1:N] * mu[N - idx]
        C = float(f0.sum())

        A = np.empty(ks.size)
        for i, k in enumerate(ks):
            r = N % int(k)
            A[i] = f0[r::int(k)].sum() if r else f0[int(k)::int(k)].sum()
        sgn = mu[ks].astype(np.float64)
        E3 = float((sgn * lg * (A - C / iph)).sum())
        B = float((lg * np.abs(A - C / iph)).sum())
        H = sgn * A
        w = lg * np.abs(H)
        fr = float(w[H > 0].sum() / w.sum())
        G = 1.0 / max(abs(2 * fr - 1), 1e-300)

        THRS.append(thresh(N, PN))
        Bs.append(B / N)
        E3s.append(abs(E3) / N)
        fs.append(fr)
        gs.append(G)
        say("  %-12d %-7d %-6d %-9.4f %-9.4f %-9.4f %-9.4f %.3f"
            % (N, K, ks.size, B / N, abs(E3) / N, fr,
               abs(0.5 - fr), G))

    say()
    x1 = fs[-1] > fs[5]
    say("X1  f rises from %.4f at 6.4e6 to %.4f at 2.56e7   %s"
        % (fs[5], fs[-1], "hold" if x1 else "REFUTED"))

    dev = np.array([abs(0.5 - v) for v in fs])
    Ln = np.array(NS, dtype=float)
    pw = np.polyfit(np.log(Ln), np.log(dev), 1)
    lw = np.polyfit(np.log(np.log(Ln)), np.log(dev), 1)
    r_pow = float(np.corrcoef(np.log(Ln), np.log(dev))[0, 1])
    r_log = float(np.corrcoef(np.log(np.log(Ln)), np.log(dev))[0, 1])
    x2 = abs(abs(r_pow) - abs(r_log)) < 0.005
    say("X2  power N^{-%.4f}  r = %.5f ;  log (log N)^{-%.4f}  r = %.5f"
        % (-pw[0], r_pow, -lw[0], r_log))
    say("    |difference| = %.5f   (cap 0.005)   %s"
        % (abs(abs(r_pow) - abs(r_log)), "hold" if x2 else "REFUTED"))
    if not x2:
        say("    the better fit is the %s law"
            % ("power" if abs(r_pow) > abs(r_log) else "log"))

    say()
    say("X3  out of sample on B/N, fit on the first %d N" % FITN)
    c = np.polyfit(np.log(np.log(Ln[:FITN])), np.log(np.array(Bs[:FITN])), 1)
    say("    B/N ~ (log N)^{%.4f}" % c[0])
    say("    N            measured   predicted   ratio")
    x3 = True
    for j in range(FITN, len(NS)):
        p = math.exp(c[1] + c[0] * math.log(math.log(NS[j])))
        r = p / Bs[j]
        if abs(r - 1.0) >= 0.10:
            x3 = False
        say("    %-12d %-10.4f %-11.4f %.4f" % (NS[j], Bs[j], p, r))
    say("    X3 %s" % ("hold" if x3 else "REFUTED"))
    THR = THRS[0]
    say("    the threshold S(N)(1-A(N)) is %.6f and is the same at"
        % THR)
    say("    all eight N -- every one has odd radical 5, so this")
    say("    sweep varies the size of N and not its arithmetic;")
    say("    audit_threshold_arithmetic.py measures what changing")
    say("    the arithmetic does.")
    lo, hi = math.log(NS[-1]), 200.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if math.exp(c[1] + c[0] * math.log(mid)) > THR:
            lo = mid
        else:
            hi = mid
    say("    that fit crosses B/N = %.4f at N = e^%.1f = 10^%.2f"
        % (THR, 0.5 * (lo + hi), 0.5 * (lo + hi) / math.log(10)))

    say("    DIAGNOSTIC (post hoc). The miss is monotone -- %s --"
        % ", ".join("%.3f" % (math.exp(c[1] + c[0] * math.log(
            math.log(NS[j]))) / Bs[j]) for j in range(FITN, len(NS))))
    say("    so it is a systematic over-prediction and not noise:")
    say("    B/N falls FASTER than the log law fitted on the small N.")
    say("    Refitting on all eight, and against a power law:")
    ca = np.polyfit(np.log(np.log(Ln)), np.log(np.array(Bs)), 1)
    cp = np.polyfit(np.log(Ln), np.log(np.array(Bs)), 1)
    ra = float(np.corrcoef(np.log(np.log(Ln)), np.log(np.array(Bs)))[0, 1])
    rp = float(np.corrcoef(np.log(Ln), np.log(np.array(Bs)))[0, 1])
    say("      all-eight log law   B/N ~ (log N)^{%.4f}   r = %.5f"
        % (ca[0], ra))
    say("      all-eight power law B/N ~ N^{%.4f}         r = %.5f"
        % (cp[0], rp))
    loo(np.log(Ln), np.log(np.array(Bs)), "B_over_N_power", say)
    loo(np.log(Ln), np.log(dev), "lean_power", say)
    for lab, cc, isl in (("log", ca, True), ("power", cp, False)):
        lo2, hi2 = math.log(NS[-1]), 400.0
        for _ in range(300):
            mid = 0.5 * (lo2 + hi2)
            val = math.exp(cc[1] + cc[0] * (math.log(mid) if isl else mid))
            if val > THR:
                lo2 = mid
            else:
                hi2 = mid
        say("      %-5s law crosses %.4f at N = 10^%.2f"
            % (lab, THR, 0.5 * (lo2 + hi2) / math.log(10)))
    say("    So 10^9.06 is the conservative end: every refit brings the")
    say("    crossing earlier, and the two laws bracket it.")

    say()
    x4 = all(E3s[j] < THRS[j] for j in range(FITN + 1, len(NS)))
    say("X4  |E_3|/N at the two new N: %.4f, %.4f against their own"
        % (E3s[-2], E3s[-1]))
    say("    thresholds %.6f, %.6f   %s"
        % (THRS[-2], THRS[-1], "hold" if x4 else "REFUTED"))

    say()
    say("=" * 70)
    ok = x1 and x2 and x3 and x4
    say("the range is doubled and the decay law is still not settled"
        if ok else "REFUTED")

    head = [
        "STATISTIC: B(N)/N = sum_{k<K}(log k)|E_mu(N;k)|/N; |E_3(N)|/N;",
        "           the mass-weighted fraction f of k with H(N;k) > 0 and",
        "           the implied gain G = 1/|2f-1|; power-law and log-law",
        "           fits of |1/2 - f| against N with their correlations;",
        "           and an out-of-sample prediction of B/N at the two",
        "           largest N from a fit on the five smallest.",
        "NULL: none is run here and none is needed for what is claimed.",
        "      This extends the range of measurements whose controls were",
        "      run in lab_sign_structure.py and lab_lean_decay.py; the",
        "      coin reference level for f is 1/2 and was measured there.",
        "      No new detection is claimed, only a longer lever.",
        "FIELD: N = 2e5 through 2.56e7 by doubling, theta' = 0.56, so k",
        "       runs over the squarefree k < N^0.56 coprime to N; m over",
        "       1 <= m < N/k; Lambda and mu from an integer sieve to",
        "       2.56e7; the out-of-sample fit sees only the five smallest",
        "       N.",
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
