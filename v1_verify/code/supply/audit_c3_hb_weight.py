# -*- coding: utf-8 -*-
"""
Re-verification of §`sec:c3` item (2) of v1/paper/wall_v1.tex --- the
Heath-Brown weight table --- and, first, the recovery of the
definition the paper does not give.

THE STATEMENT UNDER TEST, verbatim:

    "That bookkeeping needs the Mobius-side a to satisfy
     a <= y^{O(1)} = M^{o(1)}. Measuring the absolute Heath-Brown
     weight W(a) = sum_j binom(J,j) A_j(a) D_j(M/a) across a-sizes,
     the fraction with a > M^{0.05} is 0.939, 0.949, 0.960, 0.947 at
     J = 3, 4, 6, 8, rising to 0.961 and 0.969 at M = 10^6. ...
     the identity with cut z needs z^J >= x, and its j-th term has
     a <= z^j, which at j = J is x for every admissible (z,J) ---
     while the weight concentrates in exactly those high-j terms
     (0.824 at j=3 when J=3; 0.847 at j in {6,7,8} when J=8)."

WHAT IS MISSING. `A_j` and `D_j` are never defined, so the table cannot
be reproduced from the paper. This script reconstructs them from the
identity they must come from and then checks whether the numbers
follow. If they do, the definition is recovered and can be written
down; if they do not, the reconstruction is wrong or the paper's is,
and the disagreement is the finding.

THE RECONSTRUCTION. Heath-Brown's identity of level J with cut
z = x^{1/J} expands a Mobius- or Lambda-type coefficient at n <= x as

    sum_{j=1}^{J} (-1)^{j-1} binom(J,j)
        sum_{m_1...m_j n_1...n_j = n,  m_i <= z} mu(m_1)...mu(m_j) f(n_j)

so the j-th term splits every n into a "Mobius side"
a = m_1...m_j with every m_i <= z, and a complementary side
b = n_1...n_j. Taking absolute values throughout --- the paper says
ABSOLUTE weight --- the natural readings are

    A_j(a) = #{(m_1,...,m_j) : prod m_i = a, m_i <= z, mu(m_i) != 0}
    D_j(y) = sum_{b <= y} tau_j(b)

so that W(a) = sum_j binom(J,j) A_j(a) D_j(M/a) is exactly the number
of terms of the identity whose Mobius side equals a, and
sum_a W(a) is its total term count. That is the only reading under
which "the fraction of the identity's weight with a > M^{0.05}" is a
fraction of anything.

PRE-REGISTRATION (written before the run).

  (1) RULE. Under the reconstruction above, the fraction of
      sum_a W(a) carried by a > M^{0.05} must reproduce
      0.939 / 0.949 / 0.960 / 0.947 at J = 3,4,6,8 and M = 10^5, and
      0.961 / 0.969 at M = 10^6, to the printed three decimals give or
      take rounding. Anything further off means the definition is not
      the one used, and the finding is that the paper does not
      determine its own table.
  (2) SECOND RULE. The per-j concentration must reproduce: 0.824 at
      j = 3 when J = 3, and 0.847 summed over j in {6,7,8} when J = 8.
  (3) PREDICTION, recorded so it cannot be reported as a surprise.
      binom(J,j) A_j D_j grows steeply in j -- A_j counts ordered
      factorisations and D_j sums tau_j -- so the weight will sit at
      high j, and high j means large a. I therefore expect the
      fractions to come out near 0.95 whatever the exact convention,
      and the informative outcome is not the headline number but
      whether the THIRD decimal is reproducible at all. If it is not,
      the table is not a measurement a reader can check.
  (4) VARIANTS. Two conventions are defensible and both are reported:
      z = ceil(M^{1/J}) and z = floor(M^{1/J}). If the table moves in
      the third decimal between them, that alone settles (3).
"""
import sys
import math

import numpy as np

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass


def smallest_prime_factor(M):
    spf = np.zeros(M + 1, dtype=np.int32)
    for i in range(2, int(M ** 0.5) + 1):
        if spf[i] == 0:
            sl = spf[i * i::i]
            sl[sl == 0] = i
    idx = np.arange(M + 1, dtype=np.int32)
    z = spf == 0
    spf[z] = idx[z]
    spf[0] = spf[1] = 0
    return spf


def squarefree(M):
    mu2 = np.ones(M + 1, dtype=bool)
    mu2[0] = False
    for p in range(2, int(M ** 0.5) + 1):
        mu2[p * p::p * p] = False
    return mu2


def tau_j_table(M, J, spf):
    """tau_j(n) for j = 1..J, as a list of float64 arrays."""
    # factor once: n = spf^e * rest
    e = np.zeros(M + 1, dtype=np.int32)
    rest = np.arange(M + 1, dtype=np.int64)
    for n in range(2, M + 1):
        p = int(spf[n])
        m, k = n, 0
        while m % p == 0:
            m //= p
            k += 1
        e[n] = k
        rest[n] = m
    out = []
    for j in range(1, J + 1):
        t = np.zeros(M + 1, dtype=np.float64)
        t[1] = 1.0
        for n in range(2, M + 1):
            t[n] = math.comb(int(e[n]) + j - 1, j - 1) * t[int(rest[n])]
        out.append(t)
    return out


def A_tables(M, z, J, mu2):
    """A_j(a) for j = 1..J: ordered factorisations of a into j
    squarefree factors each <= z, truncated at a <= M."""
    A1 = np.zeros(M + 1, dtype=np.float64)
    ds = [d for d in range(1, z + 1) if mu2[d]]
    for d in ds:
        A1[d] = 1.0
    out = [A1]
    cur = A1
    for _ in range(2, J + 1):
        nxt = np.zeros(M + 1, dtype=np.float64)
        for d in ds:
            # nxt[k*d] += cur[k] for k = 1 .. M//d
            nxt[d::d] += cur[1: M // d + 1]
        out.append(nxt)
        cur = nxt
    return out


def run(M, J, zmode):
    spf = smallest_prime_factor(M)
    mu2 = squarefree(M)
    zf = M ** (1.0 / J)
    z = int(math.ceil(zf)) if zmode == "ceil" else max(1, int(zf))
    taus = tau_j_table(M, J, spf)
    A = A_tables(M, z, J, mu2)
    a = np.arange(M + 1)
    W = np.zeros(M + 1, dtype=np.float64)
    perj = []
    for j in range(1, J + 1):
        Dc = np.cumsum(taus[j - 1])
        # D_j(M/a) = Dc[floor(M/a)]
        q = np.zeros(M + 1, dtype=np.int64)
        q[1:] = M // a[1:]
        contrib = math.comb(J, j) * A[j - 1] * Dc[q]
        contrib[0] = 0.0
        W += contrib
        perj.append(float(contrib.sum()))
    tot = float(W.sum())
    cut = M ** 0.05
    frac = float(W[a > cut].sum()) / tot
    return z, frac, [p / tot for p in perj], tot


def main():
    print("Recovering the definition behind the Heath-Brown weight table")
    print()
    print("  A_j(a) = #{(m_1..m_j): prod m_i = a, m_i <= z, mu(m_i)!=0}")
    print("  D_j(y) = sum_{b<=y} tau_j(b)")
    print("  W(a)   = sum_j binom(J,j) A_j(a) D_j(M/a)")
    print("  z      = M^(1/J)")
    print()
    quoted = {(100000, 3): 0.939, (100000, 4): 0.949,
              (100000, 6): 0.960, (100000, 8): 0.947,
              (1000000, 6): 0.961, (1000000, 8): 0.969}
    hdr = (f"{'M':>9} {'J':>3} {'zmode':>6} {'z':>4} "
           f"{'frac a > M^0.05':>16} {'v1 quotes':>10} {'diff':>8}")
    print(hdr)
    print("-" * len(hdr))
    store = {}
    for M in (100_000, 1_000_000):
        for J in (3, 4, 6, 8):
            for zmode in ("ceil", "floor"):
                z, frac, perj, tot = run(M, J, zmode)
                store[(M, J, zmode)] = (frac, perj)
                q = quoted.get((M, J))
                qs = f"{q:.3f}" if q is not None else "--"
                ds = f"{frac - q:+.3f}" if q is not None else "--"
                print(f"{M:>9} {J:>3} {zmode:>6} {z:>4} {frac:>16.4f} "
                      f"{qs:>10} {ds:>8}")
    print()
    print("(2) per-j concentration, z = ceil(M^(1/J)), M = 1e5")
    for J in (3, 8):
        _, perj = store[(100_000, J, "ceil")]
        s = "  ".join(f"j={j+1}: {v:.3f}" for j, v in enumerate(perj))
        print(f"    J = {J}:  {s}")
        if J == 3:
            print(f"      v1 quotes 0.824 at j=3; here {perj[2]:.3f}")
        else:
            hi = sum(perj[5:8])
            print(f"      v1 quotes 0.847 at j in {{6,7,8}}; here "
                  f"{hi:.3f}")
    print()
    print("(4) how much does the convention move the table?")
    for M in (100_000, 1_000_000):
        for J in (3, 4, 6, 8):
            c = store[(M, J, 'ceil')][0]
            f = store[(M, J, 'floor')][0]
            print(f"    M = {M:>8}, J = {J}: ceil {c:.4f} vs floor "
                  f"{f:.4f}   (moves {abs(c-f):.4f})")
    print("DONE")


if __name__ == "__main__":
    main()
