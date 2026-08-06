# -*- coding: utf-8 -*-
"""
Re-audit of sweep_B's five flags against the location mask
(increment 253).

sweep_B fired five flags at |z| >= 4 on C(N) and MEASUREMENTS section 12
reads all five as one thing: "All five are what a S-dependent SCALE
predicts... So the mask is sqrt(S)". Dividing by sqrt(S) then brings the
sd ratios to 0.997/1.145/1.083 -> 1.047/0.966/0.961 and the kurtosis to
2.983, which is a genuine fix for those.

But one of the five is not a scale statistic. B5 is

    corr( C(N), S(N) ) = -0.2075,  z = -8.04,

and a scale mask CANNOT produce a correlation with S. Dividing by
sqrt(S) rescales the fluctuation; it does not move a mean. The recorded
table shows corr(., S) going -0.194, -0.138, -0.070 across three N and
that declining trend was read as convergence to zero. Increments
239-252 establish that C(N) has a LOCATION mask indexed by rad(N), which
would produce exactly a negative corr(C, S) that sqrt(S)-division leaves
behind.

WHAT IS TESTED, WITH NULLS ON THE SAME LINE.
 (A) corr(C, S) and corr(C/sqrt(S), S) over every even N in
     [1e5, 4e6], n = 1.95e6, SE = 1/sqrt(n) = 7.2e-4 -- against
     sweep_B's n = 1500, SE 2.6e-2. If the sqrt(S)-corrected
     correlation is still far from zero, the scale-only reading is
     refuted.
 (B) The same by dyadic N band, to see whether it declines with N as
     the recorded table suggested or holds.
 (C) The decisive separation: divide by the EXACT per-N second moment
     sqrt(V(N)) = sqrt(Sum_v mu^2(v) Lambda(N-v)^2) instead of
     sqrt(S). That removes every scale effect exactly, with no model
     and no fitted power, so any surviving correlation with S is
     location and nothing else. NULL: 0, SE 1/sqrt(n).
 (D) And the converse check the scale reading deserves: after dividing
     by sqrt(V), is the sd still S-dependent? If not, the scale part of
     sweep_B's reading was right and only its extension to B5 was
     wrong.
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


def conv(a, b, n_fft, X):
    A = np.zeros(n_fft); A[: X + 1] = a
    B = np.zeros(n_fft); B[: X + 1] = b
    return np.fft.irfft(np.fft.rfft(A) * np.fft.rfft(B), n_fft)[: X + 1]


def main():
    X = 4_000_000
    lo = 100_000
    t0 = time.time()
    mu, lam, primes = sieve(X)
    n_fft = 1
    while n_fft < 2 * (X + 1):
        n_fft *= 2
    C = conv(mu.astype(np.float64), lam, n_fft, X)
    V = conv(np.abs(mu).astype(np.float64), lam ** 2, n_fft, X)
    print(f"convolutions t={time.time()-t0:.0f}s", flush=True)

    C2 = 0.6601618158468696
    S = np.full(X + 1, 2 * C2)
    for p in primes:
        p = int(p)
        if p > 2:
            S[p::p] *= (p - 1) / (p - 2)

    Ns = np.arange(lo, X + 1, 2)
    Cv, Sv, Vv = C[Ns], S[Ns], V[Ns]
    n = len(Cv)
    se = 1 / math.sqrt(n)

    def r(a, b):
        return float(np.corrcoef(a, b)[0, 1])

    print(f"\n(A) the B5 statistic at full power   n = {n}, "
          f"SE = {se:.5f}")
    print(f"{'statistic':>38} {'corr with S':>12} {'z':>10}")
    for tag, x in (("C            [sweep_B's B5]", Cv),
                   ("C / sqrt(S)  [sweep_B's fix]", Cv / np.sqrt(Sv)),
                   ("C / S", Cv / Sv),
                   ("C / sqrt(V)  [exact scale removal]",
                    Cv / np.sqrt(Vv))):
        rr = r(x, Sv)
        print(f"{tag:>38} {rr:>12.4f} {rr/se:>10.1f}")
    print(f"  sweep_B ran at n = 1500, SE 2.6e-2, and read")
    print(f"  corr = -0.2075 at z = -8.04")

    print(f"\n(B) does it decline with N, as the recorded table read?")
    print(f"{'band':>22} {'n':>9} {'corr(C,S)':>11} "
          f"{'corr(C/sqrtS,S)':>17} {'corr(C/sqrtV,S)':>17}")
    b = lo
    while b < X:
        hi = min(2 * b, X)
        sel = (Ns >= b) & (Ns < hi)
        if sel.sum() > 5000:
            c_, s_, v_ = Cv[sel], Sv[sel], Vv[sel]
            print(f"{b:>10}-{hi:>11} {int(sel.sum()):>9} "
                  f"{r(c_, s_):>11.4f} "
                  f"{r(c_/np.sqrt(s_), s_):>17.4f} "
                  f"{r(c_/np.sqrt(v_), s_):>17.4f}")
        b *= 2
    print("    a scale mask cannot make any of these nonzero; only a")
    print("    location mask can")

    print(f"\n(C) the scale half of sweep_B's reading, checked")
    print(f"{'S band':>16} {'n':>9} {'sd C':>11} {'sd C/sqrtS':>12} "
          f"{'sd C/sqrtV':>12} {'mean C/sqrtV':>13}")
    edges = [1.2, 1.5, 2.0, 2.5, 3.0, 4.0, 6.0]
    for i in range(len(edges) - 1):
        sel = (Sv >= edges[i]) & (Sv < edges[i + 1])
        if sel.sum() < 2000:
            continue
        c_, s_, v_ = Cv[sel], Sv[sel], Vv[sel]
        print(f"{edges[i]:>6.2f}-{edges[i+1]:<9.2f} {int(sel.sum()):>9} "
              f"{c_.std():>11.1f} {(c_/np.sqrt(s_)).std():>12.1f} "
              f"{(c_/np.sqrt(v_)).std():>12.4f} "
              f"{(c_/np.sqrt(v_)).mean():>13.4f}")
    print("    if sd C/sqrt(V) is flat across S bands the scale is")
    print("    fully removed; the mean column is then pure location")
    print("DONE")


if __name__ == "__main__":
    main()
