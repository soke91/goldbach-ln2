# -*- coding: utf-8 -*-
"""
Re-verification of the E1 normalisation of §3.1 and of the Gaussian
stamps of Conjecture 11 (`conj:L`) of v1/paper/wall_v1.tex.

THE STATEMENTS UNDER TEST.

  (a) §3.1: "The normalisation is M_k^2, i.e. the trivial bound, not
      M_k, the square-root scale ... at the square-root normalisation
      the demand would be sum|D|^2 / sum M_k << 8.8e-6, which our own
      measurements refute directly (0.305, 0.310, 0.319 at N = 1e8)."

  (b) §3.1: "At the correct normalisation the target is cleared by
      square-root cancellation with margin (N/K)(log N)^{-2A-2} -> inf,
      though that margin is invisible at accessible N."

  (c) Table 1: E1 band ratios 0.966 / 0.950 / 0.922 at N = 1e9 over
      K ~ 1e3 / 3e3 / 1e4, "no growth"; and the Conjecture-11 stamps
      "kurtosis 2.99-3.03", "free class 0.97-1.02".

with
      D(k) = sum_{sqrt N < m <= N/k} mu(m) mu(N - m k),
      M_k  = the length of that m-range.

METHOD HERE. Written from the statement. D(k) is computed by direct
enumeration from a sieved mu; nothing is taken from v1's code. The
random-sign second moment of D(k) is computed EXACTLY as
      E_k = sum_m mu^2(m) mu^2(N - m k),
which is the correct denominator for a band ratio: it is the coin
control of Lemma 18 applied to this field, and it is NOT M_k, because
the sum only sees m with both mu's nonvanishing.

PRE-REGISTRATION (written before the run).

  (1) RULE for (a). sum|D|^2 / sum M_k must land near 0.31 and must be
      of the same order as the joint-squarefree density
      prod_q (1 - 2/q^2) = 0.32264, since square-root cancellation on
      a support of that density gives exactly that. A number far from
      it would mean either the field is not square-root or the range
      convention differs from v1's.

  (2) RULE for (b), the arithmetic. Report the margin
      (N/K)(log N)^{-2A-2} at the tested (N, K) for A = 1 and A = 2.
      The paper says it diverges but is "invisible at accessible N";
      if it is in fact BELOW 1 at the sizes where the supporting
      measurements were taken, that should be stated as a number, not
      as a word.

  (3) RULE for (c). Band ratio = sum|D|^2 / sum E_k. Under
      Conjecture 11 this is 1. Report it, the excess kurtosis of
      D(k)/sqrt(E_k), and the half-normal ratio
      mean|D|/sqrt(E_k) against sqrt(2/pi) = 0.7979.

  (4) PREDICTION, recorded so it cannot be reported as a surprise.
      I predict the band ratio sits slightly BELOW 1 and drifts down
      with K, reproducing the shape of Table 1's 0.966/0.950/0.922,
      and that the kurtosis sits at 3 within its own error bar
      sqrt(24/n_k). I expect no finding here; this is a corroboration
      run, and the informative outcome would be a band ratio that
      moves with K far more than Table 1 admits.
"""
import sys
import math

import numpy as np

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass


def sieve_mu(X):
    """mu on [0, X] in int8, memory-lean."""
    mu = np.ones(X + 1, dtype=np.int8)
    comp = np.zeros(X + 1, dtype=bool)
    lim = int(X ** 0.5) + 1
    for p in range(2, lim):
        if not comp[p]:
            comp[p * p::p] = True
    primes = np.nonzero(~comp)[0]
    primes = primes[primes >= 2]
    for p in primes:
        p = int(p)
        mu[p::p] = -mu[p::p]
        pp = p * p
        if pp <= X:
            mu[pp::pp] = 0
    return mu


def band(N, K, mu, nk=300):
    """D(k), M_k, E_k for k in [K, 2K), at most nk of them, (k,N)=1."""
    rt = int(math.isqrt(N))
    ks = [k for k in range(K, 2 * K) if math.gcd(k, N) == 1]
    if len(ks) > nk:
        step = len(ks) / nk
        ks = [ks[int(i * step)] for i in range(nk)]
    D = np.empty(len(ks))
    M = np.empty(len(ks))
    E = np.empty(len(ks))
    for i, k in enumerate(ks):
        m = np.arange(rt + 1, N // k + 1, dtype=np.int64)
        a = mu[m]
        b = mu[N - m * k]
        D[i] = float(np.dot(a.astype(np.float64), b.astype(np.float64)))
        M[i] = len(m)
        E[i] = float(np.count_nonzero((a != 0) & (b != 0)))
    return np.array(ks), D, M, E


def main():
    print("Re-verification of the E1 normalisation and the Conjecture 11")
    print("Gaussian stamps.")
    print()

    DENS2 = 1.0
    for q in [p for p in range(2, 100000)
              if all(p % r for r in range(2, int(p ** 0.5) + 1))]:
        DENS2 *= 1.0 - 2.0 / (q * q)
    print(f"joint-squarefree density prod(1-2/q^2) = {DENS2:.5f}")
    print()

    for N in (10_000_000, 40_000_000):
        N -= N % 2
        mu = sieve_mu(N)
        lgN = math.log(N)
        print(f"N = {N}   log N = {lgN:.3f}   N^(1/3) = {N**(1/3):.0f}")
        hdr = (f"    {'K':>7} {'#k':>5} {'sum|D|^2/sum M_k':>17} "
               f"{'sum|D|^2/sum M_k^2':>19} {'band ratio':>11} "
               f"{'exc kurt':>9} {'+-':>6} {'E|D|/sqrt(E_k)':>15}")
        print(hdr)
        print("    " + "-" * (len(hdr) - 4))
        for K in (100, 300, 1000):
            if 2 * K > N ** (1 / 3) * 1.05:
                pass
            ks, D, M, E = band(N, K, mu)
            r_sqrt = float((D ** 2).sum() / M.sum())
            r_triv = float((D ** 2).sum() / (M ** 2).sum())
            ratio = float((D ** 2).sum() / E.sum())
            g = D / np.sqrt(E)
            kurt = float(((g - g.mean()) ** 4).mean()
                         / ((g - g.mean()) ** 2).mean() ** 2) - 3.0
            se = math.sqrt(24.0 / len(g))
            hn = float(np.abs(g).mean())
            print(f"    {K:>7} {len(ks):>5} {r_sqrt:>17.4f} "
                  f"{r_triv:>19.3e} {ratio:>11.4f} {kurt:>9.3f} "
                  f"{se:>6.3f} {hn:>15.4f}")
        print(f"    (half-normal target sqrt(2/pi) = "
              f"{math.sqrt(2/math.pi):.4f})")
        print()
        print(f"    (2) the margin (N/K)(log N)^(-2A-2) at this N:")
        for K in (100, 300, 1000):
            for A in (1, 2):
                mg = (N / K) * lgN ** (-2 * A - 2)
                print(f"        K = {K:>5}, A = {A}: {mg:.4g}"
                      f"{'   -- BELOW 1' if mg < 1 else ''}")
        print()
        del mu

    print("(1) v1 quotes 0.305, 0.310, 0.319 for sum|D|^2/sum M_k at")
    print("    N = 1e8. Compare the third column, and compare both")
    print("    against the joint-squarefree density above.")
    print("(3) v1 quotes band ratios 0.966/0.950/0.922 at N = 1e9 over")
    print("    K ~ 1e3/3e3/1e4, and kurtosis 2.99-3.03 (excess -0.01 to")
    print("    +0.03). Compare the last three columns.")
    print("DONE")


if __name__ == "__main__":
    main()
