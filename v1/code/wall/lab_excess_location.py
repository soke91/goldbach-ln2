# -*- coding: utf-8 -*-
"""
Where does the deep-N excess live? (increment 268)

Increment 267 established, against an exact sign-randomised null, that
at N divisible by 510510 the per-prime sums satisfy

    Sum_p log p |D_p|  =  1.584 x (random-sign level),

while at all even N the same ratio is 0.959. So at deep N each |D_p| is
inflated well above what random signs would give -- yet the signed
total Sum_p log p D_p = C(N) - Lambda(N-1) stays at sqrt(N log N)
scale. The systematic components must therefore cancel across p almost
exactly. This asks where they sit.

THE DISCRIMINATION. Increment 263 measured which p carry a REMOVABLE
deterministic part: subtracting a split-estimated per-p mean reduced
the demand at p <= 8192 (reductions 0.42 to 0.68) and increased it
above that. So the mask, as far as it can be estimated, lives at small
and middle p. If the 58 percent excess sits there too, the excess IS
the mask. If it sits at large p -- where the mass is, and where
increment 263's subtraction only added noise -- then the excess is a
second phenomenon and the mask is not the whole story.

PREDICTION, fixed before the run: the excess sits at LARGE p. Reason:
for p > N/2 the sum D_p is a single term mu(p)Lambda(N-p), and at deep
N those terms are not sign-random at all -- v = p is prime, so
mu(p) = -1 always. A one-term sum with a forced sign is the extreme
case of "worse than random", and it costs nothing to check whether that
is what the ratio is measuring.

NULLS. The sign-randomised null is computed per band from the same
data, so each band's ratio carries its own reference; the spread over
draws is printed. Bands with fewer than 3 primes are printed with that
fact and not dropped.
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
    lam = np.zeros(X + 1, dtype=np.float64)
    for p in primes:
        q = int(p); lg = math.log(int(p))
        while q <= X:
            lam[q] = lg; q *= int(p)
    del spf
    return mu, lam, primes


def bandsums(N, mu, lam, ps, lp, edges, rng, ndraw):
    v = np.arange(1, N, dtype=np.int64)
    base = lam[N - v] / np.log(np.maximum(v, 2).astype(np.float64))
    mv = mu[1:N].astype(np.float64)
    t = mv * base; t[0] = 0.0
    a = np.abs(mv) * base; a[0] = 0.0
    sup = a > 0
    nb = len(edges) - 1
    S = np.zeros(nb); M = np.zeros(nb)
    D = np.zeros((ndraw, nb))
    bidx = np.searchsorted(edges, ps, side='right') - 1
    for i in range(len(ps)):
        p = int(ps[i]); b = bidx[i]
        if b < 0 or b >= nb:
            continue
        sl = slice(p - 1, None, p)
        S[b] += lp[i] * abs(float(t[sl].sum()))
        M[b] += lp[i] * float(a[sl].sum())
    ts = np.empty_like(t)
    for d in range(ndraw):
        ts[:] = 0.0
        ts[sup] = a[sup] * rng.choice(np.array([-1.0, 1.0]),
                                      size=int(sup.sum()))
        for i in range(len(ps)):
            p = int(ps[i]); b = bidx[i]
            if b < 0 or b >= nb:
                continue
            D[d, b] += lp[i] * abs(float(ts[p - 1::p].sum()))
    return S, M, D


def main():
    X = 1_600_000
    NDRAW = 2
    t0 = time.time()
    mu, lam, primes = sieve(X)
    ps = primes[primes < X]
    lp = np.log(ps.astype(np.float64))
    rng = np.random.default_rng(20260806)
    edges = [2 ** k for k in range(1, 22)]
    edges = [e for e in edges if e <= X] + [X + 1]
    print(f"sieve, {len(ps)} primes  t={time.time()-t0:.0f}s", flush=True)

    lo = X // 2
    groups = [("deep  k*30030", [k * 30030
                                 for k in range(lo // 30030 + 1,
                                                X // 30030 + 1)][:6]),
              ("all even N", list(range(lo + 2, X + 1,
                                        2 * ((X - lo) // 12)))[:6])]
    out = {}
    for name, Ns in groups:
        S = M = None
        D = None
        for N in Ns:
            s, m, d = bandsums(int(N), mu, lam, ps, lp, edges,
                               rng, NDRAW)
            S = s if S is None else S + s
            M = m if M is None else M + m
            D = d if D is None else D + d
        out[name] = (S, M, D, len(Ns))
        print(f"  {name}: {len(Ns)} values  t={time.time()-t0:.0f}s",
              flush=True)

    print(f"\n{'p band':>18} {'#p':>7} {'mass share':>11} "
          f"{'deep S/S_sign':>14} {'even S/S_sign':>14}")
    Sd, Md, Dd, _ = out["deep  k*30030"]
    Se, Me, De, _ = out["all even N"]
    for b in range(len(edges) - 1):
        npb = int(((ps >= edges[b]) & (ps < edges[b + 1])).sum())
        if npb == 0:
            continue
        rd = Sd[b] / Dd[:, b].mean() if Dd[:, b].mean() > 0 else float('nan')
        re_ = Se[b] / De[:, b].mean() if De[:, b].mean() > 0 else float('nan')
        tag = f"{npb}" + ("*" if npb < 3 else "")
        print(f"{edges[b]:>8}-{edges[b+1]:>9} {tag:>7} "
              f"{Md[b]/Md.sum():>11.4f} {rd:>14.4f} {re_:>14.4f}")
    print("    * marks a band with fewer than 3 primes")
    print(f"\n  totals: deep {Sd.sum()/Dd.sum(axis=1).mean():.4f}, "
          f"all even {Se.sum()/De.sum(axis=1).mean():.4f}")
    print("  increment 263 found the REMOVABLE part at p <= 8192; if")
    print("  the excess is instead at large p, it is a second thing")
    print("DONE")


if __name__ == "__main__":
    main()
