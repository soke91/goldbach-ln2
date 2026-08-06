# -*- coding: utf-8 -*-
"""
Does removing the mask reduce transform P's demand? (increment 263)

Proposition P.4 asks for S_abs(N) = Sum_p log p |D_p(N)| = o(N), and
increments 236 and 256 measured S_abs/N at 0.30 and falling, with the
limit undecidable. The location mask gives a new handle on it: if

    D_p(N) = m_p(N) + f_p(N)

with m_p a deterministic part computable from rad(N), then a strictly
weaker sufficient condition than P.4 is

    | Sum_p log p m_p(N) |  +  Sum_p log p |f_p(N)|  =  o(N),

because the first term is the mask itself -- of order sqrt(N log N),
hence already o(N) and free -- while the second may be materially
smaller than Sum_p log p |D_p|. The triangle inequality is applied to
the fluctuation only, which is where it costs least.

WHAT IS MEASURED. Two groups of N of matched size, one deep in the
mask (N = k*30030, so 2,3,5,7,11,13 all divide N) and one shallow
(N = 2q, q prime). Within each group the per-p mean

    m_p = mean over the group of D_p(N)

is estimated from the group itself, and then for every N in the group

    S_raw(N) = Sum_p log p |D_p(N)|        (P.4's demand)
    S_dem(N) = Sum_p log p |D_p(N) - m_p|  (demand after de-masking)

are compared. Also reported by dyadic p-band, since sessions 6 and 10
put all of transform P's difficulty at large p and the question is
whether the mask lives there too.

NULLS AND CRITERIA, on the same line.
 * Estimating m_p from the same N that are then de-masked REMOVES
   variance by construction: with n values per group, subtracting the
   sample mean shrinks Sum|D - mean| by a factor of order 1 - c/n even
   if m_p is truly zero. The honest control is therefore a SPLIT
   estimate -- m_p from the odd-indexed half of the group, applied to
   the even-indexed half -- and a PERMUTATION control in which the
   group labels are shuffled so that any true m_p is destroyed while
   the estimation noise is preserved. Both are computed and printed.
 * CRITERION: de-masking helps iff the split-estimate reduction on the
   deep group exceeds the permuted reduction by a clear margin. Any
   reduction not exceeding the permuted control is estimation
   shrinkage and not the mask.
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


def dvec(N, mu, lam, ps):
    """D_p(N) for every p < N, as in Proposition P.1."""
    v = np.arange(1, N)
    t = (mu[1:N].astype(np.float64) * lam[N - v]
         / np.log(np.maximum(v, 2).astype(np.float64)))
    t[0] = 0.0
    out = np.zeros(len(ps))
    for i, p in enumerate(ps):
        p = int(p)
        out[i] = t[p - 1::p].sum()
    return out


def main():
    X = 1_200_000
    t0 = time.time()
    mu, lam, primes = sieve(X)
    ps = primes[primes < X]
    lp = np.log(ps.astype(np.float64))
    print(f"sieve, {len(ps)} primes  t={time.time()-t0:.0f}s",
          flush=True)

    CORE = 30030
    deep = [k * CORE for k in range(1, X // CORE + 1)]
    cand = primes[np.searchsorted(primes, X // 3):]
    shal = [2 * int(q) for q in cand[:: max(1, len(cand) // 60)]
            if 2 * int(q) <= X][:len(deep)]
    groups = [("deep  k*30030", deep), ("shallow  2q", shal)]

    store = {}
    for name, Ns in groups:
        M = np.zeros((len(Ns), len(ps)))
        for i, N in enumerate(Ns):
            M[i] = dvec(int(N), mu, lam, ps)
        store[name] = (np.array(Ns), M)
        print(f"  {name}: {len(Ns)} values  t={time.time()-t0:.0f}s",
              flush=True)

    rng = np.random.default_rng(20260806)
    print(f"\n{'group':>16} {'n':>4} {'S_raw/N':>9} {'in-sample':>10} "
          f"{'split':>9} {'permuted':>10}")
    for name, Ns in groups:
        Ns_, M = store[name]
        n = len(Ns_)
        raw = (np.abs(M) @ lp) / Ns_
        # in-sample: mean from all, applied to all (optimistic)
        ins = (np.abs(M - M.mean(0)) @ lp) / Ns_
        # split: mean from odd rows, applied to even rows
        a, b = M[1::2], M[0::2]
        Na = Ns_[0::2]
        spl = (np.abs(b - a.mean(0)) @ lp) / Na
        # permuted: rows shuffled across the two groups' union is not
        # available here, so shuffle SIGNS of each row, which destroys
        # any common m_p while preserving each row's magnitude profile
        sgn = rng.choice([-1.0, 1.0], size=(len(a), 1))
        prm = (np.abs(b - (a * sgn).mean(0)) @ lp) / Na
        print(f"{name:>16} {n:>4} {raw.mean():>9.4f} {ins.mean():>10.4f} "
              f"{spl.mean():>9.4f} {prm.mean():>10.4f}")

    print("\n  in-sample is optimistic by construction; split is the")
    print("  honest number; permuted is the floor that estimation")
    print("  shrinkage alone produces")

    # P.4 is a statement about EVERY N, not about an average over N.
    # The deep group's raw demand is far above the all-N figure of
    # about 0.30 recorded at increments 236 and 256, so the question
    # is whether it declines with N there too.
    print(f"\nP.4's demand within each group, by N-tercile")
    print(f"{'group':>16} {'N range':>22} {'n':>4} {'S_raw/N':>9}")
    for name, _ in groups:
        Ns_, M = store[name]
        raw = (np.abs(M) @ lp) / Ns_
        order = np.argsort(Ns_)
        k = max(1, len(order) // 3)
        for j, sl in enumerate((order[:k], order[k:2*k], order[2*k:])):
            if len(sl) == 0:
                continue
            print(f"{name if j == 0 else '':>16} "
                  f"{int(Ns_[sl].min()):>10}-{int(Ns_[sl].max()):>11} "
                  f"{len(sl):>4} {raw[sl].mean():>9.4f}")
    print("    all-N figures for comparison: 0.3295 at 1e5, 0.2905 at")
    print("    8e5, 0.2592 at 6.4e6 (increment 236)")

    print(f"\nwhere the removed part lives, by dyadic p-band (deep group)")
    Ns_, M = store["deep  k*30030"]
    a, b = M[1::2], M[0::2]
    mhat = a.mean(0)
    print(f"{'p range':>16} {'#p':>7} {'raw share':>10} "
          f"{'demasked share':>15} {'reduction':>10}")
    bnd = 2
    while bnd < X:
        hi = min(2 * bnd, X)
        sel = (ps >= bnd) & (ps < hi)
        if sel.sum():
            r = float((np.abs(b[:, sel]) @ lp[sel]).mean())
            d = float((np.abs(b[:, sel] - mhat[sel]) @ lp[sel]).mean())
            tot = float((np.abs(b) @ lp).mean())
            print(f"{bnd:>7}-{hi:>8} {int(sel.sum()):>7} "
                  f"{r/tot:>10.4f} {d/tot:>15.4f} "
                  f"{1 - d/max(r, 1e-9):>10.4f}")
        bnd *= 2
    print("DONE")


if __name__ == "__main__":
    main()
