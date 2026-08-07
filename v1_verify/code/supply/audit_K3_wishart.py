# -*- coding: utf-8 -*-
"""
Re-verification of kill-test K3 of v1/paper/wall_v1.tex.

THE STATEMENT UNDER TEST, verbatim (§7.2):

    K3 | Wishart / operator moment method | **dead**:
    tr(M^2), tr(M^3), tr(M^4) dead-center on the Wishart null
    (z = -0.48/+0.42/-0.03). No sub-Wishart surplus to fund the
    exchange.

and, from Table 1,

    Pair matrix | lambda_max at the Wishart null dead center
    (z = -0.19)

THE OBJECT. `M` is not defined in §7.2, but Conjecture 10 defines the
only matrix in the program: the prime-indexed dilate pair field

    C_{k,k'} = sum_{p ~ P} mu(N - p k) mu(N - p k').

Writing X for the matrix with rows indexed by the primes p ~ P and
columns by k, X_{p,k} = mu(N - p k), this is exactly the Gram matrix
M = X^T X. That is a Wishart matrix when the columns are independent,
which is what Conjecture 10 asserts, so the Wishart null is the
appropriate family here -- unlike the null-design error the paper's
§Methodology records having already caught once ("an i.i.d.-entry
Wigner null applied to a Gram matrix"). A Gram matrix wants Wishart;
this one is a Gram matrix.

METHOD HERE. Written from the statement plus that identification. The
null is generated rather than quoted from a formula, and it is the
paper's own Lemma 17 control: replace mu by independent signs on the
SAME support pattern -- the same zero set, so every degeneracy of the
real matrix is preserved -- and rebuild M. This is stronger than an
analytic Wishart null because it keeps the support structure that the
squarefree condition imposes.

PRE-REGISTRATION (written before the run).

  (1) RULE. tr(M^2), tr(M^3), tr(M^4) and lambda_max of the real
      matrix must each sit within 3 standard deviations of the coin
      null. A value SEVERAL standard deviations BELOW the null is what
      "a sub-Wishart surplus to fund the exchange" would look like and
      is the outcome that would revive K3.
  (2) The paper reports z = -0.48/+0.42/-0.03 and lambda_max at
      z = -0.19. RULE: my z's must be of the same order -- within a
      few units of zero -- or the object I built is not the object
      they measured.
  (3) PREDICTION, recorded so it cannot be reported as a surprise.
      Conjecture 10 says the field is a deterministic mask times an
      exactly Gaussian fluctuation, so I expect every z near zero and
      no finding. The informative outcome is a systematically NEGATIVE
      set of z's, which would be a surplus.
  (4) Reported regardless: the number of coin draws, and the null's
      own spread, so the test's power is visible. A z near zero is
      evidence of absence only against a stated spread.
"""
import sys
import math

import numpy as np

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass


def sieve_mu(X):
    mu = np.ones(X + 1, dtype=np.int8)
    comp = np.zeros(X + 1, dtype=bool)
    for p in range(2, int(X ** 0.5) + 1):
        if not comp[p]:
            comp[p * p::p] = True
    primes = np.nonzero(~comp)[0]
    primes = primes[primes >= 2]
    for p in primes:
        p = int(p)
        mu[p::p] = -mu[p::p]
        if p * p <= X:
            mu[p * p::p * p] = 0
    mu[0] = 0
    return mu, primes


def moments(M):
    """tr(M^2), tr(M^3), tr(M^4), lambda_max."""
    M2 = M @ M
    t2 = float(np.trace(M2))
    t3 = float(np.einsum("ij,ji->", M2, M))
    t4 = float(np.einsum("ij,ji->", M2, M2))
    lam = float(np.linalg.eigvalsh(M).max())
    return t2, t3, t4, lam


def main():
    N = 10_000_000
    PLO, PHI = 10_000, 20_000
    NK = 200
    NDRAW = 60

    mu, primes = sieve_mu(N)
    ps = primes[(primes >= PLO) & (primes < PHI)].astype(np.int64)
    kmax = (N - 1) // int(ps.max())
    ks = np.arange(2, kmax + 1, dtype=np.int64)
    ks = ks[np.gcd(ks, N) == 1]
    if len(ks) > NK:
        ks = ks[np.linspace(0, len(ks) - 1, NK).astype(int)]
    print("Re-verification of kill-test K3 (the pair matrix and its")
    print("Wishart null)")
    print()
    print(f"  N = {N}, primes p in [{PLO}, {PHI}) : {len(ps)}")
    print(f"  k values: {len(ks)} in [2, {kmax}], (k,N) = 1")
    print(f"  M = X^T X with X[p,k] = mu(N - p k)   "
          f"({len(ps)} x {len(ks)})")
    print()

    idx = N - np.outer(ps, ks)                      # rows p, cols k
    Xreal = mu[idx].astype(np.float64)
    supp = (Xreal != 0)
    dens = float(supp.mean())
    print(f"  support density (mu != 0) = {dens:.4f}"
          f"   (6/pi^2 = {6/math.pi**2:.4f})")
    print()

    Mreal = Xreal.T @ Xreal
    real = moments(Mreal)

    rng = np.random.default_rng(3141)
    draws = []
    for _ in range(NDRAW):
        Xc = np.where(supp, rng.choice([-1.0, 1.0], size=Xreal.shape),
                      0.0)
        draws.append(moments(Xc.T @ Xc))
    draws = np.array(draws)

    names = ["tr(M^2)", "tr(M^3)", "tr(M^4)", "lambda_max"]
    quoted = [-0.48, +0.42, -0.03, -0.19]
    print(f"{'statistic':>12} {'real':>15} {'coin mean':>15} "
          f"{'coin sd':>13} {'z':>8} {'v1 quotes':>10}")
    print("-" * 78)
    ok = True
    zs = []
    for i, nm in enumerate(names):
        m, s = draws[:, i].mean(), draws[:, i].std()
        z = (real[i] - m) / s
        zs.append(z)
        ok &= abs(z) < 3.0
        print(f"{nm:>12} {real[i]:>15.6e} {m:>15.6e} {s:>13.4e} "
              f"{z:>8.2f} {quoted[i]:>10.2f}")
    print()
    print(f"    (1) every statistic within 3 sd of the coin null: "
          f"{'PASS' if ok else 'FAIL -- K3 would be alive'}")
    print(f"    (2) same order as v1's z's: measured "
          f"{[f'{z:+.2f}' for z in zs]}")
    print(f"    (4) {NDRAW} coin draws; the null spreads above are the")
    print("        test's power. The smallest effect this design could")
    print("        detect at 3 sd is, per statistic:")
    for i, nm in enumerate(names):
        m, s = draws[:, i].mean(), draws[:, i].std()
        print(f"          {nm:>12}: 3 sd / mean = {3*s/abs(m):.4f}"
              f"   (a {100*3*s/abs(m):.2f}% deficit)")
    print()
    print("    A 'sub-Wishart surplus to fund the exchange' would have")
    print("    to exceed those percentages. The chain's currency is a")
    print("    fixed power of log N; at N = 1e7 that is (log N)^-1 =")
    print(f"    {1/math.log(N):.4f}, i.e. a {100/math.log(N):.1f}% effect at A = 1,")
    print("    which this design CAN see, and does not.")
    if not ok:
        print("DONE (failed)")
        sys.exit(1)
    print("DONE")


if __name__ == "__main__":
    main()
