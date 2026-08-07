# -*- coding: utf-8 -*-
"""
sec:coin's major-arc factors 8.40 (q=3) and 15.16 (q=5): reachable?
(v1_verify2, Phase 1, blind.)

audit_propE_majorarcs.py measured, at N = 2^20, a ratio of the coin's
scale to |S_mu(j/q)| of 2.40 at q=3 and 2.76 at q=5, against the paper's
8.40 and 15.16. The paper states neither the N at which those were
measured nor which "a coin's" scale it means (rms, or mean modulus,
which differ by sqrt(pi/2)).

This script scans N so the reader can see whether 8.40 and 15.16 are
reachable anywhere in the range this program works in, and under which
convention.

PRE-REGISTRATION.  Decision rule: report, for N = 2^14 .. 2^24 and
q = 2,3,5,7, the ratio under both coin conventions.
  REACHABLE   : some (N, convention) in range returns 8.40 at q=3 and
                15.16 at q=5 together.
  NOT REACHABLE : none does; then the two figures cannot be placed and
                the finding is that the statistic is undefined as
                printed.

Prediction: NOT REACHABLE together. |S_mu(j/q)| for fixed small q is a
Mertens-type sum of size roughly N^{1/2+o(1)}, and the coin's scale is
exactly (6N/pi^2)^{1/2}, so the ratio should drift only slowly and stay
of order a few. I predict it does not reach 8.4 by 2^24, and that q=5
does not exceed q=3 by the factor 1.8 the paper's pair implies.
"""

import sys

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

import numpy as np


def main():
    X = 1 << 24
    print("audit_majorarc_scan   (v1_verify2 Phase 1, blind)")
    print(f"sieving mu to {X:,} ...")
    spf = np.zeros(X + 1, dtype=np.int32)
    spf[2::2] = 2
    i = 3
    while i * i <= X:
        if spf[i] == 0:
            sl = spf[i * i:: 2 * i]
            sl[sl == 0] = i
        i += 2
    rest = np.nonzero(spf == 0)[0]
    rest = rest[rest >= 2]
    spf[rest] = rest
    primes = np.nonzero(np.arange(X + 1, dtype=np.int64) == spf)[0]
    primes = primes[primes >= 2]
    del spf, rest
    mu = np.ones(X + 1, dtype=np.int8)
    mu[0] = 0
    for p in primes:
        p = int(p)
        mu[p::p] *= -1
        if p * p <= X:
            mu[p * p:: p * p] = 0
    m = mu[1:].astype(np.float64)
    n = np.arange(1, X + 1)
    sf = np.cumsum(m != 0)

    print("=" * 78)
    print("ratio = (coin scale) / |S_mu(j/q)|, rms over j coprime to q")
    print("  coin rms  = sqrt(#squarefree <= N)")
    print("  coin mean = sqrt(2/pi) * coin rms")
    print()
    Ns = [1 << e for e in range(14, 25)]
    for q in (3, 5):
        print(f"  q = {q}   [paper: "
              f"{'8.40' if q == 3 else '15.16'}]")
        print(f"    {'N':>12}{'|S_mu| rms':>13}{'coin rms':>11}"
              f"{'ratio(rms)':>12}{'ratio(mean)':>13}")
        parts = []
        for j in range(1, q):
            if np.gcd(j, q) != 1:
                continue
            ph = np.exp(2j * np.pi * (j / q) * n)
            parts.append(np.cumsum(m * ph))
        for N in Ns:
            vals = [abs(p[N - 1]) for p in parts]
            v = float(np.sqrt(np.mean(np.array(vals) ** 2)))
            rms = float(np.sqrt(sf[N - 1]))
            print(f"    {N:>12,}{v:>13.1f}{rms:>11.1f}"
                  f"{rms / v:>12.2f}{np.sqrt(2 / np.pi) * rms / v:>13.2f}")
        del parts
        print()
    print("  For reference, the same at q=2 and q=7 at N=2^24:")
    for q in (2, 7):
        vals = []
        for j in range(1, q):
            if np.gcd(j, q) != 1:
                continue
            ph = np.exp(2j * np.pi * (j / q) * n)
            vals.append(abs(complex(np.dot(m, ph.real), np.dot(m, ph.imag))))
        v = float(np.sqrt(np.mean(np.array(vals) ** 2)))
        rms = float(np.sqrt(sf[X - 1]))
        print(f"    q={q:>3}   |S_mu| rms = {v:>10.1f}   "
              f"ratio(rms) = {rms / v:>7.2f}   "
              f"ratio(mean) = {np.sqrt(2 / np.pi) * rms / v:>7.2f}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
