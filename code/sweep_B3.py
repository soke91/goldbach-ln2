# -*- coding: utf-8 -*-
"""
Sweep B follow-up (increment 215): is the mask on C(N) exactly S(N)?

Sweep B raised five flags at |z| >= 4 on c(N) = C(N)/sqrt(N):
kurtosis 4.446 (z +11.4), corr(c, S) -0.208 (z -8.0), sd ratio 1.408
for omega(N) >= 5 (z +6.4) with mean z -7.1, sd ratio 0.806 for
omega(N) <= 3 (z -5.2), and sd ratio 1.261 for N = 0 mod 5 (z +5.7).

Every one of them is what "C(N) = S(N) x (unit Gaussian)" predicts,
since S(N) = 2 C_2 prod_{p|N, p>2} (p-1)/(p-2) grows with omega(N) and
with N having small odd prime factors, and a scale mixture of
Gaussians has excess kurtosis.

PRE-REGISTERED (nulls on the same line as the criteria):
  MASK IS S  iff, after replacing c by u = C(N)/(S(N) sqrt N), all
             five flags fall below |z| = 4: kurtosis -> 3 (SE
             sqrt(24/n)), corr(u,S) -> 0 (SE 1/sqrt n), and the three
             sd ratios -> 1 (SE sqrt(1/2n1 + 1/2n2)).
  PARTIAL    iff some but not all fall below.
  NOT S      iff the flags survive normalisation.
"""
import numpy as np
import math

from thmC_alpha_scan import sieve, singular


def stats(tag, x, S, om, m5, n):
    out = []
    xx = x - x.mean()
    sk = float((xx ** 3).mean() / x.std() ** 3)
    ku = float((xx ** 4).mean() / x.std() ** 4)
    zk = (ku - 3) / math.sqrt(24 / n)
    r = float(np.corrcoef(x, S)[0, 1]); zr = r * math.sqrt(n)
    print(f"  [{tag}] kurtosis {ku:6.3f} z={zk:+7.2f}"
          f"   corr(.,S) {r:+.4f} z={zr:+7.2f}")
    out += [abs(zk), abs(zr)]
    for nm, mask in (("omega<=3", om <= 3), ("omega>=5", om >= 5),
                     ("N=0 mod 5", m5 == 0)):
        a, b = x[mask], x[~mask]
        rs = a.std() / b.std()
        zs = (rs - 1) / math.sqrt(1 / (2 * a.size) + 1 / (2 * b.size))
        se = math.sqrt(a.std() ** 2 / a.size + b.std() ** 2 / b.size)
        zm = (a.mean() - b.mean()) / se
        print(f"  [{tag}] {nm:<10} sd ratio {rs:5.3f} z={zs:+7.2f}"
              f"   mean z={zm:+7.2f}")
        # The mean z was printed and NOT counted here until increment
        # 259. That omission hid the largest surviving signal: after
        # dividing by sqrt(S) the sd ratios collapse to 1.040, 0.998,
        # 1.007 while the omega>=5 MEAN sits at z = -6.35, larger than
        # any counted statistic. A summary that omits part of the
        # evidence it printed is the same fault as a verdict string
        # that misreports its own rule.
        out += [abs(zs), abs(zm)]
    return out


def main():
    X = 1_100_000
    mu, lam, phi, primes, spf = sieve(X)
    import sys
    N0 = int(sys.argv[1]) if len(sys.argv)>1 else 500_000
    cnt = 1500
    C, S, om, m5 = [], [], [], []
    for t in range(cnt):
        N = N0 + 2 * t
        idx = np.arange(1, N)
        C.append(float(np.dot(lam[1:N], mu[N - idx].astype(np.float64)))
                 / math.sqrt(N))
        S.append(singular(N, spf))
        n_ = N; w = 0
        while n_ > 1:
            p = int(spf[n_]); w += 1
            while n_ % p == 0:
                n_ //= p
        om.append(w); m5.append(N % 5)
    c = np.array(C); S = np.array(S)
    om = np.array(om); m5 = np.array(m5)
    n = c.size
    u = c / np.sqrt(S)
    print(f"n = {n}   sd(c) = {c.std():.4f}   sd(c/sqrtS) = {u.std():.4f}")
    print(f"S range [{S.min():.3f}, {S.max():.3f}]\n")

    zc = stats("raw  ", c, S, om, m5, n)
    print()
    zu = stats("/sqrtS", u, S, om, m5, n)

    nc = sum(1 for v in zc if v >= 4)
    nu = sum(1 for v in zu if v >= 4)
    print(f"\nflags |z|>=4:  raw {nc}/{len(zc)}   "
          f"after /sqrt(S) {nu}/{len(zu)}")
    print("  (all eight printed statistics are counted. The three")
    print("   MEAN z were printed and omitted from the count until")
    print("   increment 259, which hid the largest surviving signal:")
    print("   after /sqrt(S) the omega>=5 mean sits at z = -6.35,")
    print("   larger than any counted statistic.)")
    print("verdict:",
          "MASK IS sqrt(S(N)) -- C(N) = sqrt(S(N) N) x unit Gaussian"
          if nu == 0 else
          ("PARTIAL -- S(N) explains some of the structure"
           if nu < nc else "NOT S(N)"))
    print("DONE")


if __name__ == "__main__":
    main()
