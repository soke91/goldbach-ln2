# -*- coding: utf-8 -*-
"""
Why the cell means do not concentrate (increment 304).

THE OPEN QUESTION, #113. A coin's between-cell variance is
B/T ~ 0.043 and does NOT fall with sample size: 0.054 at n = 5e4,
0.040 at n = 1.6e6. Independent samples would give (k-1)/n, i.e. 0.0045
and 0.00014 -- a factor 32 fall across the range, and 290x below what
is seen at the top. Increment 303 left the mechanism unidentified
rather than guessing.

THE MECHANISM, derived rather than guessed. For a coin,
Z(N) = C(N)/sqrt(V(N)) with C = eps * Lambda, and

    Cov(Z(N), Z(N'))
      = Sum_v Lambda(N-v) Lambda(N'-v) mu^2(v)
        / sqrt(V(N) V(N')).

Put h = N'-N. The numerator is Sum_w Lambda(w)Lambda(w+h) mu^2(N-w) --
a PRIME-PAIR count, of order S_2(h) N. The denominator is of order
A(N) N log N, because V carries the extra log from Lambda^2. So

    Cov(Z(N), Z(N'))  ~  c S_2(h) / log N,

which is O(1/log N) for EVERY pair, at every distance. A mean over any
subset of N therefore has variance ~rho_bar rather than ~1/n:

    Var(cell mean) = (1/n_c^2)[n_c Var(Z) + n_c(n_c-1) rho_bar Var(Z)]
                   -> rho_bar Var(Z)   as n_c grows.

The cell means cannot concentrate. That predicts BOTH observations at
once: the size (1/log N ~ 0.06-0.08 over this range) and the flatness
in n.

PRE-REGISTRATION (fixed before the run).

  Three quantities per band, for R = 12 coins:
    rho_bar  computed exactly from
             [(Sum Z)^2 - Sum Z^2] / (n(n-1) Var(Z)),
             which needs no model and no pairing loop;
    B/T      the between-cell variance share of #112;
    1/log N  the predicted shape.

  DECISION RULES.
    (a) MECHANISM: B/T agrees with rho_bar to within 20% in every band.
        If the cell floor is the pairwise correlation, these are the
        same number by the algebra above; if not, the derivation is
        wrong.
    (b) SHAPE: rho_bar * log N is constant across bands to within 15%.
        That is what "Cov ~ c/log N" asserts and it is what the two
        endpoint values 0.054 and 0.040 already hint at, since
        0.054/0.040 = 1.35 against log(1.6e7)/log(1.4e5) = 1.39.

  WHAT WOULD REFUTE. B/T far from rho_bar (the floor is not the
  pairwise correlation), or rho_bar * log N drifting (the 1/log N shape
  is wrong). Both are possible outcomes.

WHY IT MATTERS BEYOND #113. If this is right, EVERY statistic in this
program that averages Z over a subset of N -- every cell mean, every
de-masking, every band mean -- carries an irreducible floor of order
1/log N that does not shrink with more data. That is the same floor
that produced the bias of #106 and the contamination of #112, and it
would explain all three from one cause.
"""
import math
import time

import numpy as np

QS = [3, 5, 7, 11, 13, 17, 19, 23]


def sieve(X):
    spf = np.zeros(X + 1, dtype=np.int32)
    for i in range(2, int(X ** 0.5) + 1):
        if spf[i] == 0:
            sl = spf[i * i::i]
            sl[sl == 0] = i
    for i in range(2, X + 1):
        if spf[i] == 0:
            spf[i] = i
    mu = np.zeros(X + 1, dtype=np.int8)
    mu[1] = 1
    for i in range(2, X + 1):
        p = int(spf[i])
        j = i // p
        mu[i] = 0 if j % p == 0 else -mu[j]
    primes = np.nonzero(spf[2:] == np.arange(2, X + 1))[0] + 2
    lam = np.zeros(X + 1, dtype=np.float64)
    for p in primes:
        q = int(p)
        lg = math.log(int(p))
        while q <= X:
            lam[q] = lg
            q *= int(p)
    return mu, lam


def conv(X, a, b, nfft):
    A = np.zeros(nfft); A[: X + 1] = a
    B = np.zeros(nfft); B[: X + 1] = b
    return np.fft.irfft(np.fft.rfft(A) * np.fft.rfft(B), nfft)[: X + 1]


def stats(sig, lam, V, Ns, key, sels, X, nfft):
    C = conv(X, sig, lam, nfft)[Ns]
    Z = C / np.sqrt(V[Ns])
    out = []
    for sel in sels:
        z = Z[sel]
        n = len(z)
        s1 = float(z.sum())
        s2 = float((z * z).sum())
        var = s2 / n - (s1 / n) ** 2
        # exact mean pairwise correlation, no pairing loop needed
        rho = (s1 * s1 - s2) / (n * (n - 1)) / var
        k = key[sel]
        uniq, inv = np.unique(k, return_inverse=True)
        cnt = np.bincount(inv, minlength=len(uniq)).astype(np.float64)
        tot = np.bincount(inv, weights=z, minlength=len(uniq))
        cm = tot / cnt
        gm = s1 / n
        B = float((cnt * (cm - gm) ** 2).sum()) / n
        out.append((rho, B / var))
    return out


def main():
    X = 16_000_000
    lo = 100_000
    t0 = time.time()
    mu, lam = sieve(X)
    nfft = 1
    while nfft < 2 * (X + 1):
        nfft *= 2
    supp = (mu != 0)
    V = conv(X, supp.astype(np.float64), lam ** 2, nfft)
    Ns = np.arange(lo, X + 1, 2)
    key = np.zeros(len(Ns), dtype=np.int32)
    for i, q in enumerate(QS):
        key |= ((Ns % q) == 0).astype(np.int32) << i
    print(f"sieve + V  t={time.time()-t0:.0f}s", flush=True)

    sels, lab = [], []
    b = lo
    while b < X:
        hi = min(2 * b, X)
        sel = (Ns >= b) & (Ns < hi)
        if int(sel.sum()) > 1000:
            sels.append(sel)
            lab.append((b, hi, int(sel.sum()),
                        math.log(math.sqrt(b * hi))))
        b = hi

    R = 12
    rng = np.random.default_rng(304)
    idx = np.nonzero(supp)[0]
    rhos = np.empty((R, len(sels)))
    bts = np.empty((R, len(sels)))
    for r in range(R):
        eps = np.zeros(X + 1, dtype=np.float64)
        eps[idx] = rng.integers(0, 2, size=len(idx)) * 2.0 - 1.0
        st = stats(eps, lam, V, Ns, key, sels, X, nfft)
        rhos[r] = [x[0] for x in st]
        bts[r] = [x[1] for x in st]
        if (r + 1) % 4 == 0:
            print(f"  coin {r+1}/{R}  t={time.time()-t0:.0f}s", flush=True)

    print(f"\ncoin: mean pairwise correlation against the cell floor")
    print(f"{'band':>21} {'n':>9} {'rho_bar':>9} {'B/T':>9} "
          f"{'B/T / rho':>10} {'1/logN':>8} {'rho*logN':>9}")
    rb = rhos.mean(axis=0)
    bt = bts.mean(axis=0)
    okA = True
    prod = []
    for i, (b, hi, n, Lm) in enumerate(lab):
        ratio = bt[i] / rb[i]
        okA &= abs(ratio - 1.0) <= 0.20
        prod.append(rb[i] * Lm)
        print(f"{b:>9}-{hi:>11} {n:>9} {rb[i]:>9.5f} {bt[i]:>9.5f} "
              f"{ratio:>10.3f} {1.0/Lm:>8.5f} {rb[i]*Lm:>9.4f}")
    prod = np.array(prod)
    spread = (prod.max() - prod.min()) / prod.mean()
    okB = spread <= 0.15

    print(f"\n    (a) B/T within 20% of rho_bar in every band: "
          f"{'PASS' if okA else 'FAIL'}")
    print(f"    (b) rho_bar * log N constant to 15%: "
          f"{'PASS' if okB else 'FAIL'}"
          f"   (spread {spread:.1%}, mean {prod.mean():.4f})")
    if okA and okB:
        v = ("the cell floor IS the pairwise correlation, and that "
             "correlation is c/log N with c = "
             f"{prod.mean():.3f}")
    elif okA:
        v = ("the floor is the pairwise correlation, but its shape is "
             "not 1/log N")
    elif okB:
        v = ("the 1/log N shape holds but the floor is not the pairwise "
             "correlation")
    else:
        v = "neither clause holds; the derivation is wrong"
    print(f"    {v}")
    print(f"\n    Consequence if (a) holds: every statistic in this")
    print(f"    program that averages Z over a subset of N carries an")
    print(f"    irreducible floor of this size, whatever the sample.")
    print(f"    That is one cause for the bias of #106, the")
    print(f"    contamination of #112, and #113 itself.")
    print("DONE")


if __name__ == "__main__":
    main()
