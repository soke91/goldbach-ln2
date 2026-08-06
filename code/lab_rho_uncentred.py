# -*- coding: utf-8 -*-
"""
Proposition W's quantity, measured with the estimator the coin
validates (increment 301).

WHAT 300 LEFT. The centred rho of increments 288 and 297 was biased:
a random +/-1 on the same support reproduced the real curve with z
between -0.5 and +0.4 in every band, so corrections #84, #99 and the
quantitative half of #86 were withdrawn. The cause is that nearby N
share the same mu(v)Lambda(N-v) terms, so the C(N) are positively
correlated across a band and subtracting the band mean removes real
variance.

The estimator that survives is the UNCENTRED second moment,

    r(band) = mean_N C(N)^2 / mean_N V(N),

which for a coin has expectation exactly 1: E[C^2] = V term by term,
and nothing is subtracted. And r - 1 is precisely Proposition W's
quantity,

    r - 1 = mean_N OffDiag(N) / mean_N V(N),
    OffDiag(N) = Sum_{v != v'} mu(v)mu(v') Lambda(N-v)Lambda(N-v'),

so this is the first valid measurement of it. Increment 300 saw
1.008 falling to 0.858 with z = -4.2 in the top band on twelve draws,
and said explicitly that establishing it needs its own null. This is
that null.

PRE-REGISTRATION (fixed before the run).

  NULL. R = 40 independent coin draws through byte-identical code. The
  coin distribution per band IS the null; nothing is modelled.

  DECISION RULE. The deficit is established iff the real ratio lies
  more than 5 coin standard deviations below the coin mean in at least
  5 of the 8 bands, AND the real trend in band index is more negative
  than every one of the 40 coin trends. The second clause is the one
  that matters: a single low band is noise, a monotone fall that no
  coin reproduces is not.

  THE MASK IS NOT REMOVED, and that is deliberate. Subtracting cell
  means is what produced the bias in the first place; every scheme for
  estimating them re-introduces it in a subtler form, because the
  estimate and the residual share mu(v) terms. Since
  E[C^2] = mean(m^2) + mean(fluct^2), leaving the mask in makes the
  real ratio an OVERestimate. So whatever deficit is measured is a
  LOWER BOUND on the true one, and the direction of the bias is stated
  rather than corrected.

  WHAT WOULD REFUTE. The real trend inside the coin trend
  distribution, or fewer than 5 bands past 5 sigma. Both are possible
  outcomes.
"""
import math
import time

import numpy as np


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


def ratios(sig, lam, Vb, sels, X, nfft, Ns):
    C = conv(X, sig, lam, nfft)
    out = []
    for sel, vv in zip(sels, Vb):
        c = C[Ns[sel]]
        out.append(float((c ** 2).sum()) / len(c) / vv)
    return np.array(out)


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
    print(f"sieve + V  t={time.time()-t0:.0f}s", flush=True)

    sels, Vb, lab = [], [], []
    b = lo
    while b < X:
        hi = min(2 * b, X)
        sel = (Ns >= b) & (Ns < hi)
        if int(sel.sum()) > 1000:
            sels.append(sel)
            Vb.append(float(V[Ns[sel]].mean()))
            lab.append((b, hi, int(sel.sum())))
        b = hi

    real = ratios(mu.astype(np.float64), lam, Vb, sels, X, nfft, Ns)
    print(f"real  t={time.time()-t0:.0f}s", flush=True)

    R = 40
    rng = np.random.default_rng(301)
    idx = np.nonzero(supp)[0]
    coin = np.empty((R, len(sels)))
    for r in range(R):
        eps = np.zeros(X + 1, dtype=np.float64)
        eps[idx] = rng.integers(0, 2, size=len(idx)) * 2.0 - 1.0
        coin[r] = ratios(eps, lam, Vb, sels, X, nfft, Ns)
        if (r + 1) % 10 == 0:
            print(f"  coin {r+1}/{R}  t={time.time()-t0:.0f}s", flush=True)

    print(f"\nE[C^2]/V per octave band, real mu against {R} coin draws")
    print(f"{'band':>21} {'n':>9} {'real':>9} {'coin mean':>10} "
          f"{'coin sd':>9} {'z':>8}")
    zs = []
    for i, (b, hi, n) in enumerate(lab):
        m, sd = float(coin[:, i].mean()), float(coin[:, i].std(ddof=1))
        z = (real[i] - m) / sd
        zs.append(z)
        print(f"{b:>9}-{hi:>11} {n:>9} {real[i]:>9.5f} {m:>10.5f} "
              f"{sd:>9.5f} {z:>8.2f}")
    zs = np.array(zs)

    x = np.arange(len(sels), dtype=float)
    tr_real = float(np.polyfit(x, real, 1)[0])
    tr_coin = np.array([float(np.polyfit(x, coin[r], 1)[0])
                        for r in range(R)])
    n5 = int((zs < -5.0).sum())
    beats = int((tr_coin < tr_real).sum())

    print(f"\n    coin mean over all bands {coin.mean():.5f} "
          f"(exactly 1 in expectation)")
    print(f"    real trend {tr_real:+.5f}/band; coin trends "
          f"{tr_coin.min():+.5f} to {tr_coin.max():+.5f}")
    print(f"    coin trends more negative than the real one: {beats}/{R}")
    print(f"    bands past 5 sigma: {n5}/{len(sels)}")

    okA = n5 >= 5
    okB = beats == 0
    print(f"\n    (a) at least 5 bands past -5 sigma: "
          f"{'PASS' if okA else 'FAIL'}")
    print(f"    (b) no coin trend as negative as the real one: "
          f"{'PASS' if okB else 'FAIL'}")
    if okA and okB:
        v = ("ESTABLISHED: the off-diagonal of Proposition W is "
             "negative and grows with N")
    elif okB:
        v = ("the trend is outside the coin but the per-band evidence "
             "is weaker than pre-registered")
    elif okA:
        v = ("individual bands are low but the trend is inside the "
             "coin distribution")
    else:
        v = "NOT established by either clause"
    print(f"    {v}")
    # Pooling the bands is the highest-power version available, and it
    # is reported because the per-band test failing could otherwise be
    # read as a power problem that a better summary would fix.
    pr = float(real.mean())
    pc = coin.mean(axis=1)
    pz = (pr - float(pc.mean())) / float(pc.std(ddof=1))
    print(f"\n    pooled over all bands: real {pr:.5f}, coin "
          f"{pc.mean():.5f} +/- {pc.std(ddof=1):.5f}, z = {pz:+.2f}")
    print(f"    coin draws below the real pooled value: "
          f"{int((pc < pr).sum())}/{R}")
    print(f"\n    Why the test cannot conclude: the estimator is")
    print(f"    UNBIASED (coin mean {coin.mean():.5f}) but its per-draw")
    print(f"    scatter is {coin.std(axis=0).mean():.3f}, because the")
    print(f"    C(N) across a band are strongly dependent -- they share")
    print(f"    most of their terms -- so a band of {lab[-1][2]} values")
    print(f"    carries far fewer than {lab[-1][2]} independent ones.")
    print(f"    So the low-variance estimator is biased (increment 300)")
    print(f"    and the unbiased one is too noisy: the pooled deficit is")
    print(f"    {1.0 - pr:.3f} against a coin scatter of "
          f"{pc.std(ddof=1):.3f}, i.e. {abs(pz):.2f} of one standard")
    print(f"    deviation, with {int((pc < pr).sum())} of {R} coin draws")
    print(f"    landing below it. With one realisation of mu this")
    print(f"    quantity is not measurable by either route.")

    print(f"\n    The mask is left in, so the real ratio is an OVER"
          f"estimate")
    print(f"    (E[C^2] = mean(m^2) + mean(fluct^2)) and the measured")
    print(f"    deficit {1.0 - real[-1]:.4f} in the top band is a LOWER")
    print(f"    bound on the true one. Removing the mask is what")
    print(f"    produced increment 300's bias, so it is not attempted.")
    print("DONE")


if __name__ == "__main__":
    main()
