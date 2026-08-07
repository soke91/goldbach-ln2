# -*- coding: utf-8 -*-
"""
Every ratio this program quotes, measured both ways (increment 332)

WHY. Hazard 9 (inc. 330): two summaries of one object are not
comparable until each one's weight is stated. Increment 331 swept it at
rho and found three definitions differing by 10% at small N, with one
recorded figure needing restatement. rho is not the only ratio this
program quotes over a band, and the others have never been checked.

For a per-N pair (num, den) there are two natural band summaries:

    ratio of means   mean(num) / mean(den)
    mean of ratios   mean(num / den)

They agree only when den is constant or num/den is uncorrelated with
den. Every quantity below is quoted somewhere without saying which.

WHAT IS SWEPT, and where each is used:

    A(N) = V/W          Proposition V's local factor; increment 309
                        quotes 0.8106 with "spread 0.00%"
    rho  = C^2/V        Proposition W, #121, #175 -- included as the
                        control, since its answer is already known
    Q/N                 the squarefree density, 6/pi^2; #133 was a
                        normalisation fault on exactly this
    C^2/(W*Q)           the Cauchy-Schwarz deficit of increment 309,
                        inverted so it is a ratio of per-N quantities
    |C|/sqrt(V)         the half-normal arm of verify_all

PRE-REGISTRATION (fixed before the run).

  (A1) WHICH DISAGREE. For each ratio, per octave band, report both
       summaries and their relative gap. RULE: flag any whose gap
       exceeds 2% in any band. A flag is not an error -- it means the
       quantity has two values and every quotation of it must say
       which.

  (A2) THE CONTROL BEHAVES. rho must reproduce increment 331's
       numbers: gap about 1.5% at the top band and larger at the
       bottom. RULE: rho's gap at the smallest band exceeds its gap at
       the largest. If the control does not reproduce, the harness is
       wrong and no other row reads.

  (A3) A(N) IS THE ONE THAT MATTERS. Proposition V asserts
       V(N) = W(N) A(N) POINTWISE, so the pointwise summary -- the
       mean of ratios -- is the one that tests it, and increment 309's
       0.8106 came from the ratio of means. RULE: none; the gap is
       reported, and if it is small the two agree and 309 stands as
       quoted.
"""
import math
import sys
import time

import numpy as np

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass


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


def main():
    X = 16_000_000
    lo = 100_000
    t0 = time.time()
    mu, lam = sieve(X)
    nf = 1
    while nf < 2 * (X + 1):
        nf *= 2
    V = np.fft.irfft(np.fft.rfft(np.pad((mu != 0).astype(np.float64),
                                        (0, nf - X - 1)))
                     * np.fft.rfft(np.pad(lam ** 2, (0, nf - X - 1))),
                     nf)[: X + 1]
    C = np.fft.irfft(np.fft.rfft(np.pad(mu.astype(np.float64),
                                        (0, nf - X - 1)))
                     * np.fft.rfft(np.pad(lam, (0, nf - X - 1))),
                     nf)[: X + 1]
    W = np.cumsum(lam ** 2)
    Q = np.cumsum((mu != 0).astype(np.float64))
    Ns = np.arange(lo, X + 1, 2)
    print(f"sieve  t={time.time()-t0:.0f}s", flush=True)

    bands = []
    b = lo
    while b < X:
        hi = min(2 * b, X)
        sel = (Ns >= b) & (Ns < hi)
        if int(sel.sum()) > 1000:
            bands.append((b, hi, Ns[sel]))
        b = hi

    def pair(name, num_f, den_f):
        gaps = []
        print(f"\n  {name}")
        print(f"{'band':>21} {'ratio of means':>16} "
              f"{'mean of ratios':>16} {'gap':>8}")
        for b0, hi, Nb in bands:
            nu, de = num_f(Nb), den_f(Nb)
            rm = float(nu.mean() / de.mean())
            mr = float((nu / de).mean())
            g = abs(mr / rm - 1.0)
            gaps.append(g)
            print(f"{b0:>9}-{hi:>11} {rm:>16.6f} {mr:>16.6f} "
                  f"{g:>8.2%}")
        return np.array(gaps)

    G = {}
    G["A(N) = V/W"] = pair("A(N) = V / W   (Proposition V's local factor)",
                           lambda Nb: V[Nb], lambda Nb: W[Nb])
    G["rho = C^2/V"] = pair("rho = C^2 / V   (control, #175)",
                            lambda Nb: C[Nb] ** 2, lambda Nb: V[Nb])
    G["Q/N"] = pair("Q(N) / N   (squarefree density, #133)",
                    lambda Nb: Q[Nb],
                    lambda Nb: Nb.astype(np.float64))
    G["C^2/(W Q)"] = pair("C^2 / (W Q)   (the Cauchy-Schwarz deficit, "
                          "inverted)",
                          lambda Nb: C[Nb] ** 2,
                          lambda Nb: W[Nb] * Q[Nb])
    G["|C|/sqrt V"] = pair("|C| / sqrt(V)   (the half-normal arm)",
                           lambda Nb: np.abs(C[Nb]),
                           lambda Nb: np.sqrt(V[Nb]))

    print(f"\n(A1) which ratios have two values")
    print(f"{'quantity':<16} {'max gap':>9} {'top band':>10} "
          f"{'bottom band':>12}  flag")
    flagged = []
    for k, g in G.items():
        fl = g.max() > 0.02
        if fl:
            flagged.append(k)
        print(f"{k:<16} {g.max():>9.2%} {g[-1]:>10.2%} "
              f"{g[0]:>12.2%}  {'FLAG' if fl else ''}")

    rg = G["rho = C^2/V"]
    okA2 = rg[0] > rg[-1]
    print(f"\n    (A2) the control reproduces (rho's gap larger at the "
          f"bottom): {'PASS' if okA2 else 'FAIL'}  "
          f"({rg[0]:.2%} against {rg[-1]:.2%})")
    ag = G["A(N) = V/W"]
    print(f"    (A3) A(N)'s gap is {ag.max():.2%} at worst and "
          f"{ag[-1]:.2%} at the top band")

    if okA2 and not flagged:
        v = ("no ratio in this sweep has two values at the 2% level, "
             "so rho was the exception and the rest of the program's "
             "quoted ratios are unambiguous")
    elif okA2:
        v = (f"{len(flagged)} of {len(G)} quoted ratios have two values "
             f"at the 2% level: {', '.join(flagged)}. Each is a "
             f"quantity whose every quotation must say which summary "
             f"it is, and A(N)'s gap of {ag.max():.2%} is the one that "
             f"touches Proposition V")
    else:
        v = ("the control does not reproduce increment 331's ordering, "
             "so this harness is not measuring what that one measured "
             "and no row reads")
    print(f"\n    {v}")
    print("DONE")


if __name__ == "__main__":
    main()
