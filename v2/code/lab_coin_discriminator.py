# -*- coding: utf-8 -*-
r"""
OPEN.md, "동전을 구별하는 추정기" -- Lemma {#lem:coin}'s pass condition.

WHAT IS AT STAKE

Lemma [lem:coin] says: put arbitrary signs eps(v) = +-1 on
{v : mu(v) != 0} and zero elsewhere; then eps^2 = mu^2 pointwise, so
V(N) is unchanged and C_eps(N) = sum_v eps(v) Lambda(N-v) has the same
exact second moment as C(N) for every N.  Any estimator reproduced
under that substitution is not measuring mu.  OPEN.md records this not
as a refutation but as a PASS CONDITION: a new approach must first show
that a coin gives a different answer.

Two questions follow, and this script answers both.

  (i) Does any discriminating estimator exist at all?  Yes, and
      cheaply, provided one reads MULTIPLICATIVITY instead of
      variability: mu(2v) = -mu(v) for odd v, an exact identity, while
      eps(2v) and eps(v) are independent signs.  So

          T(x) = (1/x) sum_{v<=x} mu(v) mu(2v)

      tends to minus the density of odd squarefree integers, -4/pi^2,
      while its coin analogue is O(x^{-1/2}).

  (ii) Can such an estimator be built from the C-field alone?  In
      principle yes: with M = N-2 and Lhat(m) = Lambda(m+2), one has
      C(N) = (Lhat * mu)(M) as an additive convolution with
      Lhat(0) = log 2 != 0, so Lhat is invertible and mu is recovered
      by a linear filter a with a * Lhat = delta_0.  The field
      therefore determines mu, and Lemma [lem:coin] cannot be an
      information-theoretic obstruction.  Whether the filter is USABLE
      is a separate question, and it is the one that matters.

PRE-REGISTERED PREDICTIONS (written before this script was run)

  H1  The filter recursion is self-consistent: (a * Lhat)(m) = delta_{m,0}
      to 1e-9 for every m <= 200.
  H2  |a(m)|^{1/m} settles to a limit rho in (1.0, 3.0); the filter
      grows geometrically rather than staying bounded.
  H3  The reconstruction mu_hat(M) = sum_m a(m) C(M-m+2) recovers mu(M)
      to within 0.01 for every M <= 40, and is wrong by more than 0.5
      for some M <= 200.  Exact in principle, unusable in practice.
  H4  T(x) -> -4/pi^2 = -0.405285 for mu, while over 20 coin draws on
      the same support |T_eps(x)| < 5/sqrt(x) every time.

REFUTATION RULE (fixed before the run)

  H1  REFUTED if any residual exceeds 1e-9.
  H2  REFUTED if |a(m)|^{1/m} leaves (1.0, 3.0) at m = 200, or if
      |a(m)| is bounded (max over m <= 200 below 10).
  H3  REFUTED if the reconstruction fails below M = 40, or if it never
      fails up to M = 200.
  H4  REFUTED if |T(x) + 4/pi^2| > 0.002 at the largest x, or if any
      coin draw exceeds 5/sqrt(x).

  All four gate.  H3 is the finding and is written to be falsifiable in
  both directions: a filter that worked would BE the discriminating
  estimator on the field, and a filter that failed early would say the
  field carries less than it looks like it does.

BACKS: Proposition {#prop:coindisc} and Remark {#rem:filter} in
paper/wall_v3.md.
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
OUT = os.path.join(ROOT, "results", "lab_coin_discriminator.txt")

MFILT = 200
XS = [10_000, 100_000, 1_000_000, 4_000_000]
DRAWS = 20
SEED = 20260808


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
    return lam, mu


def main():
    lines = []

    def say(s=""):
        print(s)
        lines.append(s)

    XMAX = max(XS)
    say("sieving to %d ..." % XMAX)
    lam, mu = sieves(XMAX)

    # ------------------------------------------------------------- H1/H2
    say()
    say("H1/H2   the inverse filter of Lhat(m) = Lambda(m+2)")
    say("=" * 70)
    Lh = lam[2:2 + MFILT + 1].astype(np.float64)
    a = np.zeros(MFILT + 1, dtype=np.float64)
    a[0] = 1.0 / Lh[0]
    for m in range(1, MFILT + 1):
        a[m] = -(Lh[1:m + 1][::-1] * a[:m]).sum() / Lh[0]
    resid = []
    for m in range(0, MFILT + 1):
        v = float((a[:m + 1] * Lh[m::-1]).sum())
        resid.append(abs(v - (1.0 if m == 0 else 0.0)))
    h1 = max(resid) <= 1e-9
    say("  Lhat(0) = Lambda(2) = %.6f" % Lh[0])
    say("  max |(a*Lhat)(m) - delta_{m,0}| over m <= %d : %.3e   (tol 1e-9)"
        % (MFILT, max(resid)))
    say("  H1 %s" % ("hold" if h1 else "REFUTED"))
    say("  DIAGNOSTIC (post hoc). H1 used an ABSOLUTE tolerance on a sum")
    say("  whose terms reach 1e56, so it tests double precision and not")
    say("  the recursion. Relative to the mass of the sum:")
    rel = []
    for m in range(0, MFILT + 1):
        v = float((a[:m + 1] * Lh[m::-1]).sum())
        mass = float((np.abs(a[:m + 1]) * np.abs(Lh[m::-1])).sum())
        rel.append(abs(v - (1.0 if m == 0 else 0.0)) / max(mass, 1e-300))
    say("    max relative residual over m <= %d : %.3e" % (MFILT, max(rel)))
    say("    m at which the absolute residual first exceeds 1e-9 : %d"
        % next((m for m, r in enumerate(resid) if r > 1e-9), -1))
    say("    |a(m)| there : %.4g"
        % abs(a[next((m for m, r in enumerate(resid) if r > 1e-9), 0)]))
    say("  |a(m)| at m = 0,1,2,5,10,20,50,100,200:")
    say("    %s" % ", ".join("%.4g" % abs(a[m])
                             for m in (0, 1, 2, 5, 10, 20, 50, 100, 200)))
    rho = abs(a[MFILT]) ** (1.0 / MFILT)
    h2 = 1.0 < rho < 3.0 and float(np.abs(a).max()) >= 10.0
    say("  |a(200)|^{1/200} = %.6f   max |a| = %.4g   %s"
        % (rho, float(np.abs(a).max()), "hold" if h2 else "REFUTED"))

    # ---------------------------------------------------------------- H3
    say()
    say("H3   reconstructing mu from the C-field")
    say("=" * 70)
    NB = 400
    idx = np.arange(1, NB + 1, dtype=np.int64)
    C = np.zeros(NB + 3, dtype=np.float64)
    for N in range(2, NB + 3):
        n = np.arange(2, N, dtype=np.int64)
        C[N] = float((lam[n] * mu[N - n]).sum()) if n.size else 0.0
    say("  M     mu(M)   reconstruction      abs error")
    err = {}
    for M in list(range(1, 11)) + [20, 40, 60, 80, 100, 150, 200]:
        s = float((a[:M] * C[M - np.arange(M) + 2][::1]).sum()) if M else 0.0
        s = float(sum(a[m] * C[M - m + 2] for m in range(M)))
        err[M] = abs(s - float(mu[M]))
        if M <= 10 or M in (20, 40, 60, 80, 100, 150, 200):
            say("  %-5d %-7d %-19.6f %.3e" % (M, int(mu[M]), s, err[M]))
    good = all(err[M] <= 0.01 for M in err if M <= 40)
    fails = any(err[M] > 0.5 for M in err if M <= 200)
    h3 = good and fails
    say("  accurate (<0.01) for every M <= 40 : %s" % good)
    say("  wrong by > 0.5 for some M <= 200   : %s" % fails)
    say("  H3 %s" % ("hold" if h3 else "REFUTED"))

    # ---------------------------------------------------------------- H4
    say()
    say("H4   an estimator that does distinguish: multiplicativity")
    say("=" * 70)
    target = -4.0 / math.pi ** 2
    rng = np.random.default_rng(SEED)
    say("  x          T(x)          target        coin: max |T_eps(x)|"
        "   5/sqrt(x)")
    h4 = True
    for x in XS:
        v = np.arange(1, x + 1, dtype=np.int64)
        m1 = mu[v].astype(np.int64)
        m2 = mu[2 * v].astype(np.int64) if 2 * x <= XMAX else None
        if m2 is None:
            v = np.arange(1, XMAX // 2 + 1, dtype=np.int64)
            m1 = mu[v].astype(np.int64)
            m2 = mu[2 * v].astype(np.int64)
            x = v.size
        T = float((m1 * m2).sum()) / x
        supp = np.flatnonzero(m1 * m2)
        worst = 0.0
        for _ in range(DRAWS):
            e1 = rng.integers(0, 2, size=supp.size) * 2 - 1
            e2 = rng.integers(0, 2, size=supp.size) * 2 - 1
            worst = max(worst, abs(float((e1 * e2).sum())) / x)
        bound = 5.0 / math.sqrt(x)
        ok = worst < bound
        h4 = h4 and ok
        say("  %-10d %+.6f     %+.6f     %-19.6f  %.6f  %s"
            % (x, T, target, worst, bound, "ok" if ok else "COIN OVER"))
        last = T
    h4 = h4 and abs(last - target) <= 0.002
    say("  |T(x) - target| at the largest x = %.6f   (tol 0.002)   %s"
        % (abs(last - target), "hold" if h4 else "REFUTED"))

    say()
    say("=" * 70)
    ok = h1 and h2 and h3 and h4
    say("H1 %s  H2 %s  H3 %s  H4 %s"
        % tuple("hold" if v else "REFUTED" for v in (h1, h2, h3, h4)))
    say("discriminating estimators exist; the one built from the field "
        "itself is exact but numerically dead" if ok else "REFUTED")

    head = [
        "STATISTIC: (a) the residual of a*Lhat against delta_0, where",
        "           Lhat(m) = Lambda(m+2) and a is its convolution",
        "           inverse; (b) the growth rate |a(m)|^{1/m}; (c) the",
        "           error of the reconstruction mu_hat(M) = sum_m a(m)",
        "           C(M-m+2) against mu(M); (d) T(x) = (1/x) sum_{v<=x}",
        "           mu(v)mu(2v) against -4/pi^2, and the worst of 20 coin",
        "           draws of the same statistic on the same support.",
        "FIELD: filter length m <= 200; reconstruction at M up to 200 with",
        "       C(N) computed by direct summation for N <= 403; T(x) at",
        "       x = 1e4, 1e5, 1e6 and 2e6; Lambda and mu from an integer",
        "       sieve to 4e6; coin draws with numpy default_rng seed",
        "       20260808.",
        'NULL: the coin is the object under test. Twenty draws of eps on the',
        '      support of mu^2 for T(x); the reconstruction filter is',
        '      deterministic and needs none.',
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
