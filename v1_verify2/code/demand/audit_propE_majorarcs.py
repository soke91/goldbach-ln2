# -*- coding: utf-8 -*-
"""
prop:E's circle-method margins, and sec:coin's major-arc factors.
(v1_verify2, Phase 1, blind.)

STATEMENTS UNDER TEST, verbatim:

  prop:E  "The measured margin N/(sup|S_mu| ||S_Lambda||_1) is
           0.168, 0.175, 0.158, 0.152 at N = 2^14, ..., 2^20 --- below 1
           and decaying."
          "Cauchy-Schwarz gives ~(6/pi^2)^{1/2} N (log N)^{1/2}"

  sec:coin "On the rationals j/q with small q, the Mobius exponential
            sums are markedly smaller than a coin's --- by a factor 8.40
            at q=3 and 15.16 at q=5."

PRE-REGISTRATION.

  Decision rule.
   (a) compute sup_alpha |S_mu(alpha)| and ||S_Lambda||_1 exactly on a
       fine grid (the sup over a grid is a lower bound for the true sup,
       so the computed margin is an UPPER bound; state that), and report
       the margin at every N = 2^14..2^20. REPRODUCED if four of the
       seven land on the quoted values.
   (b) form the ratio |S_mu(j/q)| against a coin's expected size on the
       same support, for q = 3 and q = 5, and report it. REPRODUCED if
       8.40 and 15.16 come out.
       The coin's scale must be stated before the ratio is quoted: for
       random signs on the squarefree support, E|S_eps(j/q)|^2
       = #{v <= N : mu(v) != 0} = (6/pi^2)N, so the coin's rms is
       sqrt(6N/pi^2). Both the rms and the half-normal mean
       sqrt(2/pi)*rms are reported, since the paper does not say which
       "a coin's" means -- that is exactly the failure mode its own rule
       4 names.

  Predictions written before running.
   (a) The margins REPRODUCE in size (order 0.15-0.18) but I predict the
       sequence is not monotone, since the paper's own four quoted values
       already rise before falling.
   (b) I predict the ratio depends on which coin scale is used, by a
       factor sqrt(pi/2) = 1.2533, and that the paper's two numbers pick
       one of the two without saying which.
"""

import sys

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

import numpy as np


def sieve_mu_lambda(X):
    spf = np.zeros(X + 1, dtype=np.int64)
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
    primes = np.nonzero(np.arange(X + 1) == spf)[0]
    primes = primes[primes >= 2]
    mu = np.ones(X + 1, dtype=np.int8)
    mu[0] = 0
    for p in primes:
        p = int(p)
        mu[p::p] *= -1
        if p * p <= X:
            mu[p * p:: p * p] = 0
    lam = np.zeros(X + 1, dtype=np.float64)
    lam[primes] = np.log(primes.astype(np.float64))
    for p in primes[primes * primes <= X]:
        q = int(p) * int(p)
        lp = float(np.log(p))
        while q <= X:
            lam[q] = lp
            q *= int(p)
    return mu, lam


def main():
    print("audit_propE_majorarcs   (v1_verify2 Phase 1, blind)")
    print("=" * 74)
    print()
    print("--- (a) prop:E, the measured margin ------------------------------")
    print("    margin := N / (sup_alpha|S_mu| * ||S_Lambda||_1)")
    print("    sup and L1 are computed on a grid of 8N points; a grid sup")
    print("    UNDER-estimates the true sup, so the margin printed is an")
    print("    upper bound on the true margin.")
    print()
    print(f"    {'N':>10}{'sup|S_mu|':>14}{'sup/sqrt(N)':>13}"
          f"{'||S_L||_1':>13}{'margin':>10}{'paper':>9}")
    quoted = {14: 0.168, 16: 0.175, 18: 0.158, 20: 0.152}
    for e in range(14, 21):
        X = 1 << e
        mu, lam = sieve_mu_lambda(X)
        G = 8 * X
        smu = np.abs(np.fft.rfft(mu[1:X + 1].astype(np.float64), G))
        sup = float(smu.max())
        del smu
        sl = np.abs(np.fft.fft(lam[1:X + 1], G))
        l1 = float(sl.mean())
        del sl
        margin = X / (sup * l1)
        q = quoted.get(e, None)
        print(f"    {X:>10,}{sup:>14.1f}{sup / np.sqrt(X):>13.4f}"
              f"{l1:>13.2f}{margin:>10.4f}"
              f"{('' if q is None else f'{q:.3f}'):>9}")
    print()
    print("    the paper names seven exponents (2^14,...,2^20) and gives")
    print("    four values, so the abscissa of each is not recoverable.")
    print()

    print("--- (b) sec:coin, the major arcs ---------------------------------")
    X = 1 << 20
    mu, lam = sieve_mu_lambda(X)
    m = mu[1:X + 1].astype(np.float64)
    Q = float((m != 0).sum())
    rms = np.sqrt(Q)                       # coin rms on the same support
    halfnormal = np.sqrt(2 / np.pi) * rms  # coin mean |.|
    n = np.arange(1, X + 1)
    print(f"    N = {X:,},  squarefree support = {Q:,.0f},  "
          f"(6/pi^2)N = {6 / np.pi ** 2 * X:,.0f}")
    print(f"    coin rms        = {rms:.1f}")
    print(f"    coin mean |.|   = {halfnormal:.1f}   "
          f"(= sqrt(2/pi) x rms)")
    print()
    print(f"    {'q':>4}{'|S_mu(j/q)| rms over j':>24}{'ratio vs coin rms':>20}"
          f"{'ratio vs coin mean':>21}{'paper':>9}")
    paper = {3: 8.40, 5: 15.16}
    for q in (2, 3, 4, 5, 6, 7, 11, 143):
        vals = []
        for j in range(1, q):
            if np.gcd(j, q) != 1:
                continue
            s = np.exp(2j * np.pi * (j / q) * n)
            vals.append(abs(float(np.dot(m, s.real)) +
                            1j * float(np.dot(m, s.imag))))
        v = np.sqrt(np.mean(np.array(vals) ** 2))
        p = paper.get(q, None)
        print(f"    {q:>4}{v:>24.2f}{rms / v:>20.2f}{halfnormal / v:>21.2f}"
              f"{('' if p is None else f'{p:.2f}'):>9}")
    print()
    print("    the two coin conventions differ by sqrt(pi/2) = "
          f"{np.sqrt(np.pi / 2):.4f}; the paper states neither.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
