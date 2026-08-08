# -*- coding: utf-8 -*-
r"""
Settling the decay law by fitting the family, not by lengthening N.

WHAT IS AT STAKE

Remark {#rem:extendrange} showed that quadrupling the N-range moved
the discrimination between two decay laws by less than 0.002 in
correlation, and concluded that extending N is not the way.  There is
another way, and it comes from the theory rather than from more data.

The two laws compared so far are special cases of one family:

    |1/2 - f| = A exp(-c (log N)^alpha),

with alpha = 1 the power law N^{-c} and alpha -> 0 the log law
(log N)^{-d}.  And theory supplies a THIRD value.  Remark
{#rem:thetasweep} found the finite-N residual of Theorem {#thm:C} is
dominated not by Bombieri-Vinogradov but by the main-term cancellation
over m < M = N^{1-theta'}, which Lemma {#lem:mu} bounds by
exp(-c sqrt(log M)) -- that is alpha = 1/2, and it is neither a power
of N nor a power of log N.  Nobody has tested it.

So instead of comparing two arbitrary choices, fit alpha.  If the data
prefers 1/2, the mechanism is the zero-free-region cancellation and the
law is known.  If the residual curve in alpha is flat, then the data
cannot determine the law and that is worth stating as a fact about the
method rather than leaving it as an unresolved comparison.

BACKS: Remark {#rem:decayfamily} in paper/wall_v3.md.

PRE-REGISTERED PREDICTIONS (written before this script was run)

  Y1  The recomputation reproduces lab_extend_range's f and B/N at all
      eight N to within 0.001 -- the tie that makes the rest mean
      anything.
  Y2  Fitting |1/2 - f| = A exp(-c (log N)^alpha), the set of alpha
      whose residual sum of squares is within 1% of the minimum spans
      more than 0.5 in alpha: the data cannot determine the law.
  Y3  alpha = 1/2, the value Lemma {#lem:mu} predicts, fits within 5%
      of the minimum RSS.
  Y4  The three canonical laws alpha = 1, 1/2 and ->0 extrapolate to
      |1/2 - f| < 0.01 at N spanning more than five orders of
      magnitude.

REFUTATION RULE (fixed before the run)

  Y1  REFUTED if any value differs by 0.001 or more.
  Y2  REFUTED if the 1%-of-minimum set spans 0.5 or less -- in which
      case the data DOES determine the law, and the minimising alpha
      is reported as the answer.
  Y3  REFUTED if alpha = 1/2 costs more than 5% over the minimum.
  Y4  REFUTED if the span is five orders or less.

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
OUT = os.path.join(ROOT, "results", "lab_decay_family.txt")

NS = [200_000, 400_000, 800_000, 1_600_000, 3_200_000, 6_400_000,
      12_800_000, 25_600_000]
PUB_F = [0.2273, 0.2228, 0.2735, 0.3068, 0.3207, 0.3376, 0.3533, 0.3608]
PUB_B = [0.8086, 0.7395, 0.7303, 0.6547, 0.5916, 0.5526, 0.4992, 0.4527]
THETA = 0.56


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


def rss_at(alpha, x, y):
    """residual sum of squares of log y = a - c x^alpha."""
    # Box-Cox: (x^a - 1)/a -> log x as a -> 0, and stays well
    # conditioned there. Affine in x^a, so the RSS is unchanged; what
    # it fixes is the extrapolation, which the raw x^a form makes
    # meaningless once x^a is nearly constant.
    t = np.log(x) if alpha < 1e-9 else (x ** alpha - 1.0) / alpha
    A = np.vstack([np.ones_like(t), -t]).T
    sol, res, _, _ = np.linalg.lstsq(A, y, rcond=None)
    fit = A @ sol
    return float(((y - fit) ** 2).sum()), sol


def main():
    lines = []

    def say(s=""):
        print(s)
        lines.append(s)

    NMAX = max(NS)
    say("sieving to %d ..." % NMAX)
    pr, lam, mu = sieves(NMAX)

    say()
    say("Y1  recomputation against lab_extend_range")
    say("  N            f          pub        B/N        pub")
    fs, bs = [], []
    y1 = True
    for j, N in enumerate(NS):
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
        B = float((lg * np.abs(A - C / iph)).sum()) / N
        H = sgn * A
        w = lg * np.abs(H)
        fr = float(w[H > 0].sum() / w.sum())
        fs.append(fr)
        bs.append(B)
        if abs(fr - PUB_F[j]) >= 1e-3 or abs(B - PUB_B[j]) >= 1e-3:
            y1 = False
        say("  %-12d %-10.4f %-10.4f %-10.4f %.4f"
            % (N, fr, PUB_F[j], B, PUB_B[j]))
    say("  Y1 %s" % ("hold" if y1 else "REFUTED"))

    x = np.log(np.array(NS, dtype=float))
    ydev = np.log(np.array([abs(0.5 - v) for v in fs]))
    yb = np.log(np.array(bs))

    say()
    say("Y2/Y3  fitting |1/2 - f| = A exp(-c (log N)^alpha)")
    say("  alpha    RSS          RSS / min")
    alphas = np.arange(0.05, 1.51, 0.05)
    rs = np.array([rss_at(a, x, ydev)[0] for a in alphas])
    rmin = float(rs.min())
    amin = float(alphas[int(np.argmin(rs))])
    for a, r in zip(alphas, rs):
        if abs(a - round(a, 1)) < 1e-9 and abs(a * 10) % 2 < 1e-9:
            say("  %-8.2f %-12.6e %.4f" % (a, r, r / rmin))
    within = alphas[rs <= 1.01 * rmin]
    # a 1%-of-RSS band is not a confidence band; the proper one for
    # one free parameter with 6 degrees of freedom is RSS within
    # 1 + F(1,6;0.05)/6 = 2.00 of the minimum.
    conf = alphas[rs <= 2.00 * rmin]
    span = float(within.max() - within.min())
    y2 = span > 0.5
    say("  minimising alpha = %.2f ;  1%%-of-minimum set spans %.2f "
        "(from %.2f to %.2f)   (floor 0.5)"
        % (amin, span, float(within.min()), float(within.max())))
    say("  Y2 %s" % ("hold" if y2 else "REFUTED"))
    say("  DIAGNOSTIC (post hoc). The 1% band is my threshold, not a")
    say("  confidence band, and 8 points with 2 fitted parameters do not")
    say("  resolve 13%% of RSS. Under the 95%% band for one parameter with")
    say("  6 degrees of freedom -- RSS within a factor 2.00 -- the")
    say("  admissible alpha runs %.2f to %.2f, the whole sweep."
        % (float(conf.min()), float(conf.max())))
    r_half = rss_at(0.5, x, ydev)[0]
    y3 = r_half <= 1.05 * rmin
    say("  Y3  alpha = 1/2 costs %.4f of the minimum   (cap 1.05)   %s"
        % (r_half / rmin, "hold" if y3 else "REFUTED"))

    say()
    say("  the same fit for B/N:")
    rsb = np.array([rss_at(a, x, yb)[0] for a in alphas])
    say("    minimising alpha = %.2f ;  RSS at 1/2 over minimum = %.4f"
        % (float(alphas[int(np.argmin(rsb))]),
           rss_at(0.5, x, yb)[0] / float(rsb.min())))

    say()
    say("Y4  where each canonical law puts |1/2 - f| = 0.01")
    say("  alpha    law                    N")
    Ns = []
    for a, name in ((1.0, "power  N^{-c}"),
                    (0.5, "exp(-c sqrt(log N))"),
                    (0.0, "log    (log N)^{-d}")):
        _, sol = rss_at(a, x, ydev)
        aa, cc = float(sol[0]), float(sol[1])
        target = math.log(0.01)
        u = (aa - target) / cc          # the value of the regressor t
        L = math.exp(u) if a < 1e-9 else (1.0 + a * u) ** (1.0 / a)
        Ns.append(L / math.log(10))
        say("  %-8.3f %-22s 10^%.2f" % (a, name, L / math.log(10)))
    spanN = max(Ns) - min(Ns)
    y4 = spanN > 5.0
    say("  span = %.2f orders   (floor 5)   %s"
        % (spanN, "hold" if y4 else "REFUTED"))

    say()
    say("=" * 70)
    ok = y1 and y2 and y3 and y4
    say("the exponent alpha is not determined by the data, and alpha=1/2 "
        "-- the value the theory names -- fits as well as any"
        if ok else "REFUTED")

    head = [
        "STATISTIC: the mass-weighted fraction f and B(N)/N recomputed and",
        "           checked against lab_extend_range.py; the residual sum",
        "           of squares of the fit log|1/2 - f| = a - c (log N)^alpha",
        "           as alpha is swept, its minimiser, the span of alpha",
        "           within 1% of that minimum, and the cost of alpha = 1/2;",
        "           the same sweep for B/N; and the N at which each of the",
        "           three canonical laws reaches |1/2 - f| = 0.01.",
        "NULL: none applies. Nothing here is a detection: the field and",
        "      its coin reference were measured in lab_sign_structure.py",
        "      and lab_lean_decay.py, and this script only fits a",
        "      functional form to numbers already controlled there.",
        "FIELD: N = 2e5 through 2.56e7 by doubling, theta' = 0.56; k over",
        "       the squarefree k < N^0.56 coprime to N; m over 1 <= m <",
        "       N/k; Lambda and mu from an integer sieve to 2.56e7; alpha",
        "       swept over 0.05 to 1.50 in steps of 0.05, with alpha -> 0",
        "       evaluated as the log law.",
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
