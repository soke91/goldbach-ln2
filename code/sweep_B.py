# -*- coding: utf-8 -*-
"""
Sweep B (increment 214): seven hypotheses on the wall's scalar
C(N) = Sum_{n<N} Lambda(n) mu(N-n). Nulls stated with the criteria.

c(N) = C(N)/sqrt(N); under Conjecture L the c(N) over even N are
iid mean-zero.

 B1  lag 1..4 autocorrelation of c   null: 0, SE 1/sqrt(n)
 B2  skewness of c                   null: 0, SE sqrt(6/n)
 B3  kurtosis of c                   null: 3, SE sqrt(24/n)
 B4  sign balance of C(N)            null: 1/2, SE 1/(2 sqrt n)
 B5  corr(c, S(N))                   null: 0, SE 1/sqrt(n)
 B6  omega(N) bins                   null: flat
 B7  N mod 8 and N mod 5 bins        null: flat

ALIVE iff |z| >= 4 anywhere (family of ~20 statistics).

Note on circularity, stated because it limits what this can show: the
null used here is Conjecture L's own (iid mean-zero). A pass therefore
confirms nothing independently -- it only fails to detect a deviation.
"""
import numpy as np
import math

from thmC_alpha_scan import sieve, singular


def split_report(name, d, mask):
    a, b = d[mask], d[~mask]
    if a.size < 30 or b.size < 30:
        print(f"  {name:<24} (too few)"); return 0
    se = math.sqrt(a.std() ** 2 / a.size + b.std() ** 2 / b.size)
    zm = (a.mean() - b.mean()) / se if se > 0 else 0.0
    rs = a.std() / b.std() if b.std() > 0 else float('nan')
    zr = (rs - 1) / math.sqrt(1 / (2 * a.size) + 1 / (2 * b.size))
    f = int(abs(zm) >= 4 or abs(zr) >= 4)
    print(f"  {name:<24} n={a.size:>5}/{b.size:<5} mean z={zm:+6.2f} "
          f" sd ratio={rs:5.3f} (z={zr:+6.2f}) {'<<<' if f else ''}")
    return f


def main():
    X = 1_100_000
    mu, lam, phi, primes, spf = sieve(X)
    N0, cnt = 500_000, 1500
    cs, Ss, oms, mods8, mods5 = [], [], [], [], []
    Ns = []
    for t in range(cnt):
        N = N0 + 2 * t
        idx = np.arange(1, N)
        C = float(np.dot(lam[1:N], mu[N - idx].astype(np.float64)))
        cs.append(C / math.sqrt(N))
        Ss.append(singular(N, spf))
        n = N; w = 0
        while n > 1:
            p = int(spf[n]); w += 1
            while n % p == 0:
                n //= p
        oms.append(w); mods8.append(N % 8); mods5.append(N % 5)
        Ns.append(N)
    c = np.array(cs); S = np.array(Ss)
    om = np.array(oms); m8 = np.array(mods8); m5 = np.array(mods5)
    n = c.size
    print(f"N in [{Ns[0]}, {Ns[-1]}] even, n = {n}")
    print(f"mean c = {c.mean():+.4f}  sd c = {c.std():.4f}\n")

    flags = 0
    cc = c - c.mean(); v = float(np.dot(cc, cc))
    for l in range(1, 5):
        r = float(np.dot(cc[:-l], cc[l:]) / v); z = r * math.sqrt(n)
        print(f"  B1 lag {l}   r={r:+.4f}  z={z:+6.2f}"
              f" {'<<<' if abs(z) >= 4 else ''}")
        flags += int(abs(z) >= 4)

    sk = float((cc ** 3).mean() / c.std() ** 3)
    ku = float((cc ** 4).mean() / c.std() ** 4)
    zs = sk / math.sqrt(6 / n); zk = (ku - 3) / math.sqrt(24 / n)
    print(f"  B2 skewness {sk:+.4f} (null 0)  z={zs:+6.2f}"
          f" {'<<<' if abs(zs) >= 4 else ''}")
    print(f"  B3 kurtosis {ku:.4f} (null 3)  z={zk:+6.2f}"
          f" {'<<<' if abs(zk) >= 4 else ''}")
    flags += int(abs(zs) >= 4) + int(abs(zk) >= 4)

    pos = float(np.count_nonzero(c > 0)) / n
    zp = (pos - 0.5) / (1 / (2 * math.sqrt(n)))
    print(f"  B4 P(C>0)  {pos:.4f} (null 0.5)  z={zp:+6.2f}"
          f" {'<<<' if abs(zp) >= 4 else ''}")
    flags += int(abs(zp) >= 4)

    r = float(np.corrcoef(c, S)[0, 1]); z = r * math.sqrt(n)
    print(f"  B5 corr(c,S) {r:+.4f} (null 0)  z={z:+6.2f}"
          f" {'<<<' if abs(z) >= 4 else ''}")
    flags += int(abs(z) >= 4)

    print("\n  B6-B7 splits:")
    flags += split_report("B6 omega(N)<=3", c, om <= 3)
    flags += split_report("B6 omega(N)>=5", c, om >= 5)
    flags += split_report("B7 N=2 mod 8", c, m8 == 2)
    flags += split_report("B7 N=0 mod 8", c, m8 == 0)
    flags += split_report("B7 N=0 mod 5", c, m5 == 0)

    print(f"\nSWEEP B: {flags} flags at |z| >= 4")
    print("DONE")


if __name__ == "__main__":
    main()
