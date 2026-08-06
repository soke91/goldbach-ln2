# -*- coding: utf-8 -*-
"""
Sweep A (increment 214): ten hypotheses on the local structure of the
dilate field D(k). Every criterion carries its null on the same line.

Normalised field d(k) = D(k)/sqrt(supp(k)); under Conjecture L the
d(k) are iid mean-zero with sd ~ 1 (half-normal |d| mean 0.798).

 A1  gcd(k,N) > 1 vs = 1        null: mean-z ~ N(0,1), sd ratio ~ 1
 A2  omega(k) bins              null: flat
 A3  v_2(k) bins                null: flat
 A4  smooth k (P+(k) <= sqrt k) null: flat
 A5  k mod 3 and k mod 4 bins   null: flat
 A6  skewness of d              null: 0, SE sqrt(6/n)
 A7  kurtosis of d              null: 3, SE sqrt(24/n)
 A8  runs test on sign(d)       null: n/2 runs, SE sqrt(n)/2
 A9  lag 1..8 autocorrelation   null: 0, SE 1/sqrt(n)
 A10 long lag (n/4, n/2)        null: 0, SE 1/sqrt(n)

ALIVE for any test iff |z| >= 4 (family of ~30 statistics, so 4 sigma
is the Bonferroni-aware bar); DEAD otherwise.
"""
import numpy as np
import math

from e1_forge_r4 import mobius_upto


def split_report(name, d, mask):
    a, b = d[mask], d[~mask]
    if a.size < 30 or b.size < 30:
        print(f"  {name:<26} (too few)"); return 0
    se = math.sqrt(a.std() ** 2 / a.size + b.std() ** 2 / b.size)
    zm = (a.mean() - b.mean()) / se if se > 0 else 0.0
    rs = a.std() / b.std() if b.std() > 0 else float('nan')
    se_r = math.sqrt(1 / (2 * a.size) + 1 / (2 * b.size))
    zr = (rs - 1) / se_r
    flag = int(abs(zm) >= 4 or abs(zr) >= 4)
    print(f"  {name:<26} n={a.size:>5}/{b.size:<5} "
          f"mean z={zm:+6.2f}  sd ratio={rs:5.3f} (z={zr:+6.2f}) "
          f"{'<<<' if flag else ''}")
    return flag


def main():
    N = 9_999_998
    mu = mobius_upto(N)
    SQ = int(N ** 0.5)
    ks = np.arange(300, SQ, dtype=np.int64)     # k < sqrt(N) or empty
    D = np.zeros(len(ks)); S = np.zeros(len(ks))
    for i, k in enumerate(ks):
        k = int(k)
        ms = np.arange(SQ + 1, N // k + 1, dtype=np.int64)
        a = mu[ms].astype(np.int64); b = mu[N - k * ms].astype(np.int64)
        D[i] = float((a * b).sum()); S[i] = float(np.count_nonzero(a * b))
    keep = S > 0
    ks, D, S = ks[keep], D[keep], S[keep]
    d = D / np.sqrt(S)
    n = d.size
    print(f"N = {N}, k in [{ks[0]}, {ks[-1]}], n = {n}")
    print(f"mean d = {d.mean():+.4f}  sd d = {d.std():.4f}  "
          f"mean|d| = {np.abs(d).mean():.4f} (null 0.798)\n")

    # smallest prime factor table for k
    lim = int(ks[-1]) + 1
    spf = np.zeros(lim, dtype=np.int32)
    for i in range(2, int(lim ** 0.5) + 1):
        if spf[i] == 0:
            sl = spf[i * i::i]; sl[sl == 0] = i
    for i in range(2, lim):
        if spf[i] == 0:
            spf[i] = i

    def facts(k):
        f = []
        while k > 1:
            p = int(spf[k]); f.append(p)
            while k % p == 0:
                k //= p
        return f

    om = np.array([len(facts(int(k))) for k in ks])
    pmax = np.array([max(facts(int(k))) for k in ks])
    v2 = np.array([(int(k) & -int(k)).bit_length() - 1 for k in ks])
    g = np.array([math.gcd(int(k), N) for k in ks])

    flags = 0
    print("A1-A5 splits:")
    flags += split_report("A1 gcd(k,N)>1", d, g > 1)
    flags += split_report("A2 omega(k)<=2", d, om <= 2)
    flags += split_report("A2 omega(k)>=4", d, om >= 4)
    flags += split_report("A3 v2(k)=0", d, v2 == 0)
    flags += split_report("A3 v2(k)>=2", d, v2 >= 2)
    flags += split_report("A4 smooth P+<=sqrt k", d,
                          pmax <= np.sqrt(ks.astype(float)))
    flags += split_report("A5 k=0 mod 3", d, ks % 3 == 0)
    flags += split_report("A5 k=1 mod 4", d, ks % 4 == 1)

    print("\nA6-A10 shape and memory:")
    sk = float(((d - d.mean()) ** 3).mean() / d.std() ** 3)
    ku = float(((d - d.mean()) ** 4).mean() / d.std() ** 4)
    zs = sk / math.sqrt(6 / n); zk = (ku - 3) / math.sqrt(24 / n)
    print(f"  A6 skewness  {sk:+.4f} (null 0)   z={zs:+6.2f}"
          f" {'<<<' if abs(zs) >= 4 else ''}")
    print(f"  A7 kurtosis  {ku:.4f} (null 3)   z={zk:+6.2f}"
          f" {'<<<' if abs(zk) >= 4 else ''}")
    flags += int(abs(zs) >= 4) + int(abs(zk) >= 4)

    s = np.sign(d); s = s[s != 0]
    runs = 1 + int(np.count_nonzero(np.diff(s) != 0))
    m_ = s.size
    zr = (runs - m_ / 2) / (math.sqrt(m_) / 2)
    print(f"  A8 runs      {runs} (null {m_/2:.0f})   z={zr:+6.2f}"
          f" {'<<<' if abs(zr) >= 4 else ''}")
    flags += int(abs(zr) >= 4)

    dd = d - d.mean(); v = float(np.dot(dd, dd))
    for l in list(range(1, 9)) + [n // 4, n // 2]:
        r = float(np.dot(dd[:-l], dd[l:]) / v)
        z = r * math.sqrt(n)
        tag = "A9 " if l <= 8 else "A10"
        print(f"  {tag} lag {l:<5} r={r:+.4f}  z={z:+6.2f}"
              f" {'<<<' if abs(z) >= 4 else ''}")
        flags += int(abs(z) >= 4)

    print(f"\nSWEEP A: {flags} flags at |z| >= 4")
    print("DONE")


if __name__ == "__main__":
    main()
