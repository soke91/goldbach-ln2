# -*- coding: utf-8 -*-
"""
The location mask, derived (increment 247): R_A(N) = prod_{q not| N}
(q-3)/(q-1), and the sharp prediction that 3 not| N kills the mask.

Increment 246 established the identity

    C(N) = ( Sum_j T_j ) * R_A(N),
    R_A(N) = Sum_j (-1)^j T_j / Sum_j T_j,

with T_j the Lambda-weighted count of shifted primes v = N - p having
omega(v) = j, and Sum_j T_j = Sum_p log p mu^2(N-p) the trivial bound.
So R_A(N) IS the mask. It also corrected the mechanism: for q not
dividing N, q | v forces N - v = N mod q, which is nonzero, so a small
prime factor of v HELPS N - v be prime.

THE DERIVATION THAT FOLLOWS. Condition on N - v being prime. Then
v must avoid the class v = N mod q, leaving q - 1 classes, of which
exactly one is v = 0 (and 0 = N mod q is impossible when q does not
divide N). Hence

    P(q | v) = 1/(q-1)     for q not dividing N,
    P(q | v) = 0           for q dividing N.

On the support of mu, mu = (-1)^omega, so treating the divisibility
indicators as independent,

    R_A(N) = E[(-1)^omega(v)] = prod_q (1 - 2 P(q|v))
           = prod_{q <= z, q not| N} (q-3)/(q-1).

TWO PREDICTIONS, BOTH SHARP, BOTH MADE BEFORE THE MEASUREMENT.
 (1) The factor at q = 3 is (3-3)/(3-1) = 0 EXACTLY. So if 3 does not
     divide N the mask vanishes to this order. This is not a fitted
     shape; it is a zero.
 (2) Restricted to 3 | N, the per-prime factor for q dividing N is the
     reciprocal of the omitted term:
       f(q) = (q-1)/(q-3):  2.000 at q=5, 1.500 at 7, 1.250 at 11,
       1.200 at 13, 1.143 at 17, 1.125 at 19, 1.100 at 23.

CRITERIA, with nulls on the same line.
 * For (1): mean |R_A| over cells with 3 | N against cells with 3 not|
   N. NULL is a ratio of 1 (no difference). CONFIRMED iff the ratio
   exceeds 5.
 * For (2): fitted f(q) from the cell means, restricted to 3 | N,
   against (q-1)/(q-3). CONFIRMED iff within 20 percent.
 * A permuted-label control is run on the same statistic so the fitted
   spread has a floor. It is printed whether or not it is flattering.

TRUNCATION. The product over q <= z diverges to zero as z grows, like
(log z)^{-2}, since sum 2/(q-1) diverges. So the ABSOLUTE size of R_A
is not predicted here, only its dependence on rad(N); one overall
constant is fitted and that is the only free quantity.
"""
import numpy as np
import math
import time


def sieve(X):
    spf = np.zeros(X + 1, dtype=np.int32)
    for i in range(2, int(X ** 0.5) + 1):
        if spf[i] == 0:
            sl = spf[i * i::i]; sl[sl == 0] = i
    for i in range(2, X + 1):
        if spf[i] == 0:
            spf[i] = i
    mu = np.zeros(X + 1, dtype=np.int8); mu[1] = 1
    for i in range(2, X + 1):
        p = int(spf[i]); j = i // p
        mu[i] = 0 if j % p == 0 else -mu[j]
    primes = np.nonzero(spf[2:] == np.arange(2, X + 1))[0] + 2
    lam = np.zeros(X + 1)
    for p in primes:
        q = int(p); lp = math.log(int(p))
        while q <= X:
            lam[q] = lp; q *= int(p)
    return mu, lam, primes


def main():
    X = 4_000_000
    lo = 100_000
    QS = [3, 5, 7, 11, 13, 17, 19, 23]
    t0 = time.time()
    mu, lam, primes = sieve(X)

    n_fft = 1
    while n_fft < 2 * (X + 1):
        n_fft *= 2
    F = np.zeros(n_fft); G = np.zeros(n_fft)
    F[: X + 1] = mu; G[: X + 1] = lam
    C = np.fft.irfft(np.fft.rfft(F) * np.fft.rfft(G), n_fft)[: X + 1]
    F[: X + 1] = np.abs(mu)
    T = np.fft.irfft(np.fft.rfft(F) * np.fft.rfft(G), n_fft)[: X + 1]
    del F, G
    print(f"convolutions t={time.time()-t0:.0f}s", flush=True)

    Ns = np.arange(lo, X + 1, 2)
    RA = C[Ns] / T[Ns]                     # the mask, by the identity
    cell = np.zeros(len(Ns), dtype=np.int64)
    for i, q in enumerate(QS):
        cell |= ((Ns % q == 0).astype(np.int64) << i)
    ncell = 1 << len(QS)

    div3 = (Ns % 3 == 0)
    m3, m3n = float(RA[div3].mean()), float(RA[~div3].mean())
    s3 = float(RA[div3].std() / math.sqrt(div3.sum()))
    s3n = float(RA[~div3].std() / math.sqrt((~div3).sum()))
    print(f"\n(1) the zero at q = 3 -- does 3 not| N kill the mask?")
    print(f"  3 | N      mean R_A = {m3:+.6f} +- {s3:.6f}"
          f"   n = {int(div3.sum())}")
    print(f"  3 not| N   mean R_A = {m3n:+.6f} +- {s3n:.6f}"
          f"   n = {int((~div3).sum())}")
    print(f"  ratio |3|N| / |3 not|N| = {abs(m3/m3n):.2f}"
          f"    (null 1.00, criterion > 5)")
    print(f"  predicted: the q = 3 factor is (3-3)/(3-1) = 0 exactly")

    print(f"\n(2) per-prime factors, restricted to 3 | N")
    print(f"{'q':>4} {'mean R_A, q|N':>15} {'mean R_A, q not|N':>18} "
        f"{'f(q) meas':>10} {'(q-1)/(q-3)':>12} {'perm f':>8}")
    rng = np.random.default_rng(20260806)
    perm = rng.permutation(len(Ns))
    RAp = RA[perm]
    for q in QS[1:]:
        base = div3
        a = base & (Ns % q == 0)
        b = base & (Ns % q != 0)
        if a.sum() < 500:
            continue
        ma, mb = float(RA[a].mean()), float(RA[b].mean())
        pa, pb = float(RAp[a].mean()), float(RAp[b].mean())
        f = ma / mb if mb else float('nan')
        fp = pa / pb if pb else float('nan')
        print(f"{q:>4} {ma:>15.6f} {mb:>18.6f} {f:>10.3f} "
              f"{(q-1)/(q-3):>12.3f} {fp:>8.3f}")

    print(f"\n(3) the whole formula against the cell means (3 | N only)")
    pred = np.ones(ncell)
    for c in range(ncell):
        if not (c & 1):                    # 3 does not divide N
            pred[c] = 0.0
            continue
        for i, q in enumerate(QS):
            if not (c >> i & 1) and q > 3:
                pred[c] *= (q - 3.0) / (q - 1.0)
    rows = []
    for c in range(ncell):
        idx = np.nonzero(cell == c)[0]
        if len(idx) < 200 or not (c & 1):
            continue
        rows.append((float(RA[idx].mean()), pred[c], len(idx), c))
    a = np.array([r[0] for r in rows]); b = np.array([r[1] for r in rows])
    w = np.array([r[2] for r in rows], dtype=np.float64)
    k0 = float(np.dot(w * a, b) / np.dot(w * b, b))
    r2 = 1 - float((w * (a - k0 * b) ** 2).sum()) / \
        float((w * (a - np.average(a, weights=w)) ** 2).sum())
    print(f"  one fitted constant k0 = {k0:+.6f},  weighted R^2 "
          f"= {r2:+.4f}   ({len(rows)} cells)")
    rows.sort()
    print(f"{'primes | N':>26} {'count':>7} {'measured':>10} "
          f"{'predicted':>10} {'ratio':>7}")
    for (mc, pc, cnt, c) in rows[:10]:
        lab = "*".join(str(q) for i, q in enumerate(QS) if c >> i & 1)
        print(f"{lab:>26} {cnt:>7} {mc:>10.5f} {k0*pc:>10.5f} "
              f"{(k0*pc)/mc if mc else float('nan'):>7.3f}")
    print("DONE")


if __name__ == "__main__":
    main()
