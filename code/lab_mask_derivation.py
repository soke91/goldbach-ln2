# -*- coding: utf-8 -*-
"""
Deriving the location mask (increment 244): is the wall's deterministic
term the Mobius bias over numbers coprime to rad(N)?

Increment 243 closed the guessing phase: additive, Euler-multiplicative
and singular-series-power models all fail to reproduce the enumerated
mask. What was missing is a derivation, and there is one available that
has not been tried.

THE MECHANISM, STATED EXACTLY. If n is prime and q | N, then q | N - n
would force q | n, impossible for n > q. Hence, up to the omega(N)
terms with n = q,

    C(N) = Sum_{v < N, (v, rad N) = 1} mu(v) Lambda(N - v)
           + O(omega(N) log N),

so the Mobius variable is confined to the integers coprime to rad(N).
That is not a neutral restriction. Sieving out the small primes leaves
the ROUGH numbers, among which the primes -- where mu = -1, with no
sign variation at all -- are over-represented relative to their share
of all integers. The Mobius sum over such a set carries a bias, and the
bias grows as more small primes are removed.

THE HYPOTHESIS. Averaging Lambda(N - v) over the surviving v gives a
factor that is itself computable, not a free constant. The numerator is

    Sum_{v<N, (v,rad N)=1} Lambda(N-v) = Sum_{n<N} Lambda(n)
        - Sum_{n : q | N-n for some q | N} Lambda(n),

and for q | N the condition q | N - n forces q | n, i.e. n = q^k, so the
subtracted part is O(omega(N) log^2 N) and the numerator is psi(N) ~ N.
The denominator is the count of such v, namely N prod_{q|N}(1 - 1/q).
Hence

    kappa(N) = prod_{q | N} q/(q-1),      C(N) ~= kappa_0 kappa(N) M_{rad N}(N),
    M_P(x) := Sum_{v <= x, (v,P) = 1} mu(v),

with kappa_0 a single absolute constant and no free shape at all.

This is a prediction for EVERY N, not for a cell mean, and M_P(N) is a
concrete computable number. If it holds, the mask is not a fitted table
but the partial sums of the Mobius function over rough numbers, and the
enumeration of increment 240 was measuring exactly that.

WHAT IS TESTED, WITH NULLS ON THE SAME LINE.
 (A) corr(C(N), M(N)) over all even N in the range. NULL: the same
     correlation after permuting M within each cell, which destroys the
     N-by-N link while preserving every cell mean -- so a nonzero
     permuted correlation would mean the cell structure alone produces
     it. SE of the null is 1/sqrt(n).
 (B) the fitted kappa per cell. The hypothesis says kappa is a property
     of the Lambda average, hence roughly constant across cells;
     PREDICTION made before running: within a factor 1.5.
 (C) how much of the CELL-MEAN structure M reproduces: the cell means
     of kappa*M against the cell means of C, on the cells increment 240
     enumerated. This is the direct comparison with the fitted mask.
 (D) the deepest cells one by one, since those are where every previous
     model broke.

TRUNCATION, STATED. rad(N) is replaced by its part supported on
q <= 17. Larger prime factors of N impose a coprimality of density
1 - 1/q, close to no restriction, so the truncation should cost little
-- but it is a truncation and (D) will show if the deepest cells suffer
from it.
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
    QS = [3, 5, 7, 11, 13, 17]          # 2 is always in rad(N) here
    t0 = time.time()
    mu, lam, primes = sieve(X)

    n_fft = 1
    while n_fft < 2 * (X + 1):
        n_fft *= 2
    F = np.zeros(n_fft); F[: X + 1] = mu
    G = np.zeros(n_fft); G[: X + 1] = lam
    C = np.fft.irfft(np.fft.rfft(F) * np.fft.rfft(G), n_fft)[: X + 1]
    del F, G
    print(f"convolution t={time.time()-t0:.0f}s", flush=True)

    Ns = np.arange(lo, X + 1, 2)
    cell = np.zeros(len(Ns), dtype=np.int64)
    for i, q in enumerate(QS):
        cell |= ((Ns % q == 0).astype(np.int64) << i)

    v = np.arange(X + 1)
    muf = mu.astype(np.float64)
    Mval = np.zeros(len(Ns))
    ncell = 1 << len(QS)
    for c in range(ncell):
        sel = cell == c
        if not sel.any():
            continue
        keep = (v % 2 == 1)                       # 2 | N always
        for i, q in enumerate(QS):
            if c >> i & 1:
                keep &= (v % q != 0)
        Ms = np.cumsum(muf * keep)
        Mval[sel] = Ms[Ns[sel] - 1]
    print(f"M_P computed for {ncell} cells  t={time.time()-t0:.0f}s",
          flush=True)

    Cv = C[Ns]
    n = len(Cv)
    r = float(np.corrcoef(Cv, Mval)[0, 1])

    rng = np.random.default_rng(20260806)
    Mp = Mval.copy()
    for c in range(ncell):
        idx = np.nonzero(cell == c)[0]
        if len(idx) > 1:
            Mp[idx] = Mval[idx][rng.permutation(len(idx))]
    r0 = float(np.corrcoef(Cv, Mp)[0, 1])

    print(f"\n(A) does M_P predict C, N by N?")
    print(f"  corr(C, M)            = {r:+.4f}")
    print(f"  corr(C, M) permuted   = {r0:+.4f}   "
          f"(null 0, SE {1/math.sqrt(n):.5f})")
    print(f"  n = {n}")

    # kappa(N) = prod_{q|N} q/(q-1), derived above, not fitted
    kapN = np.full(len(Ns), 2.0)              # q = 2 always divides N
    for i, q in enumerate(QS):
        kapN = np.where(((cell >> i) & 1).astype(bool),
                        kapN * q / (q - 1.0), kapN)
    W = kapN * Mval
    Wp = kapN * Mp
    r_w = float(np.corrcoef(Cv, W)[0, 1])
    r_w0 = float(np.corrcoef(Cv, Wp)[0, 1])
    print(f"\n(A2) with the derived kappa(N) = prod q/(q-1) included")
    print(f"  corr(C, kappa*M)          = {r_w:+.4f}")
    print(f"  corr(C, kappa*M) permuted = {r_w0:+.4f}")

    kap = float(np.dot(Cv, W) / np.dot(W, W))
    print(f"\n(B) fitted kappa_0 overall = {kap:.4f}")
    print(f"{'primes | N':>24} {'count':>8} {'kappa_c':>9} "
          f"{'corr_c':>8}")
    ks = []
    for c in np.argsort(-np.bincount(cell, minlength=ncell)):
        idx = np.nonzero(cell == c)[0]
        if len(idx) < 500:
            continue
        cc, mm = Cv[idx], Mval[idx]
        if np.dot(mm, mm) == 0:
            continue
        ww = W[idx]
        if np.dot(ww, ww) == 0:
            continue
        k = float(np.dot(cc, ww) / np.dot(ww, ww))
        rc = float(np.corrcoef(cc, ww)[0, 1])
        lab = "*".join(str(q) for i, q in enumerate(QS) if c >> i & 1)
        ks.append(k)
        if len(ks) <= 12:
            print(f"{lab or '(none)':>24} {len(idx):>8} {k:>9.4f} "
                  f"{rc:>8.4f}")
    if ks:
        print(f"  kappa range over cells: [{min(ks):.4f}, "
              f"{max(ks):.4f}]  ratio {max(ks)/min(ks):.2f}"
              f"   (predicted within 1.5)")

    print(f"\n(C) does it reproduce the CELL MEANS the mask enumerates?")
    print(f"{'primes | N':>24} {'count':>8} {'mean C':>12} "
          f"{'k0*kap*mean M':>14} {'ratio':>7}")
    rows = []
    for c in range(ncell):
        idx = np.nonzero(cell == c)[0]
        if len(idx) < 100:
            continue
        mc = float(Cv[idx].mean()); mm = kap * float(W[idx].mean())
        rows.append((mc, mm, len(idx), c))
    rows.sort()
    for (mc, mm, cnt, c) in rows[:12]:
        lab = "*".join(str(q) for i, q in enumerate(QS) if c >> i & 1)
        print(f"{lab or '(none)':>24} {cnt:>8} {mc:>12.1f} "
              f"{mm:>13.1f} {mm/mc if mc else float('nan'):>7.3f}")
    a = np.array([r[0] for r in rows]); b = np.array([r[1] for r in rows])
    w = np.array([r[2] for r in rows], dtype=np.float64)
    r2 = 1 - float((w * (a - b) ** 2).sum()) / \
        float((w * (a - np.average(a, weights=w)) ** 2).sum())
    print(f"  weighted R^2 of kappa*M against the cell means: {r2:+.4f}")
    print(f"  (the fitted models of increment 243 gave +0.219 additive,")
    print(f"   -714 multiplicative, +0.910 for a power of Sig with an")
    print(f"   unnatural exponent 9.7 -- and this has NO free shape,")
    print(f"   only the single constant kappa)")
    print("DONE")


if __name__ == "__main__":
    main()
