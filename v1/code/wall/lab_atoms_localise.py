# -*- coding: utf-8 -*-
"""
Localising the factor-5 contradiction of #169 (increment 328)

THE CONTRADICTION. For C(N) over even N in [2e5, 8e6]:

    atomic energy SHARE      real 0.0649   coin 0.3346   ratio 5.155
    per-frequency |sum|^2/Q  real 1.1174   coin 1.0156   ratio 0.91

The spectrum of C is exactly |mu-hat(f)|^2 |Lambda-hat(f)|^2 and the
coin's is |eps-hat(f)|^2 |Lambda-hat(f)|^2, so a per-frequency ratio of
1 with a share ratio of 5 is impossible unless the two disagree
somewhere else. Increment 327 said so and stopped rather than guess a
third mechanism.

A SHARE IS A RATIO, AND I ONLY EVER LOOKED AT THE RATIO. The share is
numerator over denominator, and the two have never been reported
separately. If the numerators agree and the denominators differ by 5,
the effect is in the TOTAL energy over the window and has nothing to do
with the atoms; if the numerators differ, it is in the atoms after all
and the per-frequency measurement is measuring something else.

That is one run and it decides which.

PRE-REGISTRATION (fixed before the run).

  (L1) SPLIT THE SHARE. Report, absolutely and not as a ratio:
       numerator = sum over atomic bins of |C-hat|^2,
       denominator = sum over all non-zero bins,
       for the real and for the coin. RULE: none -- this is the
       measurement that was missing.

  (L2) WHICH ONE CARRIES THE FACTOR 5. Compute num_coin/num_real and
       den_coin/den_real. Exactly one of them should be near 5.155 and
       the other near 1. RULE: the two ratios multiply to the share
       ratio to within 1%, which is arithmetic and confirms the split
       is being read correctly.

  (L3) WHERE THE DENOMINATOR LIVES. Split the full band into
       low (bins 1..1000), atomic, and the rest, and report each side's
       share of its own total. C(N) grows like sqrt(N log N) across a
       window spanning a factor 40 in N, so a large low-frequency mass
       is expected in BOTH and its ratio is the thing to look at.

  WHAT THIS CANNOT DO. It localises; it does not explain. If the
  denominators differ by 5 the next question is why, and that is a
  separate run.
"""
import math
import sys
import time

import numpy as np

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

MOD = 30030
QP = [3, 5, 7, 11, 13]


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


def divisors(ps):
    out = [1]
    for p in ps:
        out += [d * p for d in out]
    return sorted(out)


def main():
    X = 8_000_000
    lo = 200_000
    t0 = time.time()
    mu, lam = sieve(X)
    nf = 1
    while nf < 2 * (X + 1):
        nf *= 2
    F_lam = np.fft.rfft(np.pad(lam, (0, nf - X - 1)))
    C = np.fft.irfft(np.fft.rfft(np.pad(mu.astype(np.float64),
                                        (0, nf - X - 1))) * F_lam,
                     nf)[: X + 1]
    rng = np.random.default_rng(328)
    idx = np.nonzero(mu != 0)[0]
    eps = np.zeros(nf)
    eps[idx] = rng.integers(0, 2, size=len(idx)) * 2.0 - 1.0
    Cc = np.fft.irfft(np.fft.rfft(eps) * F_lam, nf)[: X + 1]

    Ns = np.arange(lo, X + 1, 2)
    n = (len(Ns) // MOD) * MOD
    Ns = Ns[:n]
    atom = np.zeros(n // 2 + 1, dtype=bool)
    for q in divisors(QP):
        if q == 1:
            continue
        for j in range(1, q):
            if math.gcd(j, q) == 1:
                b = (j * n) // q
                if 0 < b < len(atom):
                    atom[b] = True
    print(f"n = {n}, {int(atom.sum())} atomic bins  "
          f"t={time.time()-t0:.0f}s", flush=True)

    def parts(y):
        y = y - y.mean()
        P = np.abs(np.fft.rfft(y)) ** 2
        num = float(P[atom].sum())
        den = float(P[1:].sum())
        lowm = float(P[1:1001].sum())
        rest = den - num - lowm + float(P[1:1001][atom[1:1001]].sum())
        return num, den, lowm, rest

    nr, dr, lr, rr = parts(C[Ns])
    nc, dc, lc, rc = parts(Cc[Ns])

    print(f"\n(L1) the share, split into its parts (absolute)")
    print(f"{'':<10} {'numerator':>14} {'denominator':>14} "
          f"{'share':>10}")
    print(f"{'real':<10} {nr:>14.5e} {dr:>14.5e} {nr/dr:>10.5f}")
    print(f"{'coin':<10} {nc:>14.5e} {dc:>14.5e} {nc/dc:>10.5f}")

    rn, rd = nc / nr, dc / dr
    sh = (nc / dc) / (nr / dr)
    print(f"\n(L2) which side carries the factor")
    print(f"    numerator   coin/real = {rn:>8.3f}")
    print(f"    denominator coin/real = {rd:>8.3f}")
    print(f"    product of the two    = {rn/rd:>8.3f}   "
          f"share ratio = {sh:.3f}")
    okL2 = abs((rn / rd) / sh - 1.0) < 0.01
    print(f"    (L2) the split multiplies back to the share ratio: "
          f"{'PASS' if okL2 else 'FAIL'}")

    print(f"\n(L3) where the energy sits, each side as a fraction of "
          f"its own total")
    print(f"{'':<10} {'low (1..1000)':>15} {'atoms':>10} {'rest':>10}")
    print(f"{'real':<10} {lr/dr:>15.5f} {nr/dr:>10.5f} {rr/dr:>10.5f}")
    print(f"{'coin':<10} {lc/dc:>15.5f} {nc/dc:>10.5f} {rc/dc:>10.5f}")
    print(f"    low-band coin/real = {lc/lr:>8.3f}   "
          f"rest coin/real = {rc/rr:>8.3f}")

    al = alias_test(mu, eps, QP, divisors)
    okL4 = al["odd v only"] > 2.0

    if okL2 and abs(rn - 1) < 0.5 and rd > 3:
        v = (f"the ATOMS agree ({rn:.2f}x) and the DENOMINATORS differ "
             f"({rd:.2f}x). The factor 5 is not in the atomic bins at "
             f"all -- it is in the total energy over the window, and "
             f"the per-frequency measurement of #167/#169 was right. "
             f"What has to be explained is why the real C carries "
             f"{1/rd:.2f}x the coin's total energy here, which is a "
             f"different question from the one #166 asked")
    elif okL2 and rn > 3 and okL4:
        v = (f"the ATOMS differ by {rn:.2f}x while the denominators "
             f"agree to {rd:.2f}x, and the mechanism is the EVEN-N "
             f"SUBSAMPLING. Each bin aliases f and f+1/2; Lambda lives "
             f"on odd numbers so its two aliases nearly cancel, and "
             f"what survives is 2*Lambda-hat(f)*sum_{{v odd}} mu(v) "
             f"e(vf). Restricted to odd v the coin exceeds mu by "
             f"{al['odd v only']:.2f}x, against {al['all v']:.2f}x over "
             f"all v. #166 is NOT 'mu suppresses the periodic "
             f"covariance' -- it is mu's ODD-restricted sum being "
             f"smaller than a coin's, seen through a subsampling I "
             f"chose")
    elif okL2 and rn > 3:
        v = (f"the ATOMS differ by {rn:.2f}x while the denominators "
             f"agree, and the even-N aliasing does not account for it "
             f"either: restricted to odd v the ratio is "
             f"{al['odd v only']:.2f}x. Third mechanism refuted; the "
             f"contradiction stands")
    elif okL2:
        v = (f"both sides move -- numerator {rn:.2f}x, denominator "
             f"{rd:.2f}x -- so the factor 5 is not localised to either "
             f"and the contradiction stands undiagnosed")
    else:
        v = ("the split does not multiply back to the share ratio, so "
             "the parts are not being read correctly and nothing here "
             "is usable")
    print(f"\n    {v}")
    print("DONE")




def alias_test(mu, eps, QP, divisors):
    """(L4, increment 328b) The periodogram runs over EVEN N, so each
    bin ALIASES two frequencies: f = b/(2n) and f + 1/2. Lambda is
    supported on primes, all odd but one, so Lambda-hat(f + 1/2) is
    very nearly -Lambda-hat(f), and the two aliases combine into

        Lambda-hat(f) * [ mu-hat(f) - mu-hat(f + 1/2) ]
            = 2 Lambda-hat(f) * sum_{v ODD} mu(v) e(v f).

    For a coin the two aliases are independent draws and no such
    coherent combination occurs. So the quantity the periodogram
    actually sees at an atom is the ODD-RESTRICTED sum, not mu-hat
    itself -- which is why #167 and #169 measured 1 where the
    periodogram measures 4.4."""
    print("\n(L4) the odd-restricted sums, which is what an even-N "
          "periodogram sees")
    sup = (mu != 0)
    vv = np.nonzero(sup)[0]
    mv = mu[sup].astype(np.float64)
    ev = eps[vv]
    odd = (vv % 2) == 1
    Qo = float(odd.sum())
    print(f"{'evaluated at':<22} {'mu':>12} {'coin':>12} {'coin/mu':>10}")
    out = {}
    for tag, use_odd in (("all v", False), ("odd v only", True)):
        sel = odd if use_odd else np.ones(len(vv), dtype=bool)
        vs, ms, es = vv[sel], mv[sel], ev[sel]
        Qn = float(len(vs))
        rm, re_ = [], []
        for q in divisors(QP):
            if q == 1:
                continue
            M = 2 * q
            r = (vs % M).astype(np.int64)
            am = np.bincount(r, weights=ms, minlength=M)
            ae = np.bincount(r, weights=es, minlength=M)
            for jj in range(1, q):
                if math.gcd(jj, q) != 1:
                    continue
                w = np.exp(2j * np.pi * jj * np.arange(M) / M)
                rm.append(abs(complex(np.dot(am, w))) ** 2 / Qn)
                re_.append(abs(complex(np.dot(ae, w))) ** 2 / Qn)
        rm, re_ = np.array(rm), np.array(re_)
        out[tag] = float(re_.mean() / rm.mean())
        print(f"{tag:<22} {rm.mean():>12.4f} {re_.mean():>12.4f} "
              f"{re_.mean()/rm.mean():>10.2f}")
    return out


if __name__ == "__main__":
    main()
