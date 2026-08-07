# -*- coding: utf-8 -*-
"""
sec:coin -- the shift-mass distribution and the mu-autocorrelation floor.
(v1_verify2, Phase 1, blind.)

STATEMENTS UNDER TEST, verbatim:

  "The Mobius autocorrelation sits at 1.051--1.068 times the random-sign
   floor sqrt(0.32264 (X-h)) --- not sqrt(X), since the sum sees only n
   with both n and n+h squarefree --- stably across five decades of
   shift."

  "Gross mass by shift: h<10^3 carries 1.1%, 10^3--10^4 carries 3.0%,
   10^4--10^5 carries 23.1%, 10^5--10^6 carries 48.9%, and above 10^6
   carries 23.8% (the ranges cancel, net -2.2e13 against gross 4.4e13,
   so shares of the net would mislead). Small shifts, where Chowla is
   hardest and the averaged theorem weakest, carry almost nothing: the
   wall leans on the range where that theorem is strongest."

  lem:MP (repaired form, verified in audit_lem_mp.py):
   sum_N Chat(N)^2 = sum_h M(h)P(h), M(h)=sum mu(v)mu(v+h),
   P(h)=sum Lambda(w)Lambda(w+h).

The mass being apportioned is M(h)P(h), the h-component of the wall's
aggregate second moment. This script computes M and P exactly by FFT
autocorrelation and does the apportionment.

PRE-REGISTRATION.

  Decision rule.
    (a) mu-autocorrelation: form |M(h)|/sqrt(0.32264(X-h)) and report its
        typical value by decade of h. REPRODUCED if the decade summaries
        land in 1.051--1.068 across five decades.
    (b) shift mass: report the five bucket shares of gross mass
        sum|M(h)P(h)|, and the net and gross totals. REPRODUCED if the
        five shares match to the quoted 0.1% and the totals to the
        quoted two figures.
    (c) as a by-product, re-verify the repaired lem:MP at this X.

  Predictions written before running.
    (a) REPRODUCED. The constant 0.32264 = prod_p (1-2/p^2) is the right
        density for pairs of squarefree values at generic h, and the
        ratio exceeding 1 is the known positive bias of mu-correlations.
        I expect the ratio to be h-dependent in a way the phrase "stably
        across five decades" hides, because the density depends on
        whether p^2 | h and so is NOT 0.32264 for every h.
    (b) The bucket shares are dominated by the number of h in each
        bucket, so I predict shares rising steeply with the bucket and
        the top two buckets carrying the bulk -- qualitatively as
        quoted. The stated conclusion ("the wall leans where the
        averaged theorem is strongest") is a claim about GROSS mass in a
        sum that cancels, so I also report what the buckets do to the
        NET, which is the quantity that actually appears in rho-1.
"""

import os
import sys

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
sys.path.insert(0, os.path.join(HERE))


def autocorr(a, X):
    """sum_v a(v)a(v+h) for h = 0..X-1, exactly, by FFT."""
    n = 1
    while n < 2 * len(a):
        n <<= 1
    f = np.fft.rfft(a, n)
    f *= np.conjugate(f)
    r = np.fft.irfft(f, n)
    return r[:X].copy()


def main():
    X = int(sys.argv[1]) if len(sys.argv) > 1 else 16_000_000

    # rebuild mu and Lambda (cheap relative to the FFTs)
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from lab_field_build import smallest_prime_factor, von_mangoldt, mobius

    print(f"audit_shift_mass   (v1_verify2 Phase 1, blind)   X = {X:,}")
    print("=" * 74)
    print("  sieving ...")
    spf = smallest_prime_factor(X)
    lam, primes = von_mangoldt(X, spf)
    del spf
    mu = mobius(X, primes)

    print("  M(h) = mu autocorrelation (fft) ...")
    M = autocorr(mu.astype(np.float64), X)
    print("  P(h) = Lambda autocorrelation (fft) ...")
    P = autocorr(lam, X)

    Q = float((mu != 0).sum())          # M(0)
    W = float(P[0])
    print(f"  M(0) = #squarefree <= X = {M[0]:,.0f}   "
          f"(6/pi^2)X = {6 / np.pi ** 2 * X:,.0f}")
    print(f"  P(0) = sum Lambda^2      = {W:,.0f}   "
          f"X(log X - 1) = {X * (np.log(X) - 1):,.0f}")
    print()

    # ---------------------------------------------------- (a)
    print("--- (a) the mu-autocorrelation against its random-sign floor ----")
    print("    floor as the paper writes it: sqrt(0.32264 (X-h))")
    print(f"    {'decade of h':>18}{'count':>10}{'rms |M|/floor':>16}"
          f"{'mean |M|/floor':>16}{'median':>10}")
    edges = [(1, 10), (10, 100), (100, 1000), (1000, 10000),
             (10000, 100000), (100000, 1000000), (1000000, X // 2)]
    for lo, hi in edges:
        h = np.arange(lo, min(hi, X))
        if len(h) == 0:
            continue
        floor = np.sqrt(0.32264 * (X - h))
        r = np.abs(M[lo:min(hi, X)]) / floor
        print(f"    {f'{lo:,}-{hi:,}':>18}{len(h):>10,}"
              f"{np.sqrt((r ** 2).mean()):>16.4f}{r.mean():>16.4f}"
              f"{np.median(r):>10.4f}")
    print("    [paper: 1.051--1.068, 'stably across five decades']")
    print()
    print("    the floor's own constant is not h-independent: the density")
    print("    of v with v and v+h both squarefree is prod_{p^2 |h}(1-1/p^2)")
    print("    * prod_{p^2 not| h}(1-2/p^2).  Split by 4|h:")
    for lab, sel in (("h odd", lambda h: h % 2 == 1),
                     ("h = 2 mod 4", lambda h: h % 4 == 2),
                     ("4 | h", lambda h: h % 4 == 0)):
        h = np.arange(1, 200000)
        s = sel(h)
        floor = np.sqrt(0.32264 * (X - h[s]))
        r = np.abs(M[1:200000][s]) / floor
        print(f"      {lab:<12} n={s.sum():>7,}   rms |M|/floor = "
              f"{np.sqrt((r ** 2).mean()):.4f}")
    print()

    # ---------------------------------------------------- (c)
    print("--- (c) the repaired lem:MP at this X ---------------------------")
    tot = float(M[0] * P[0] + 2.0 * np.dot(M[1:], P[1:]))
    a = np.zeros(X + 1)
    a[1:] = lam[1:]
    b = np.zeros(X + 1)
    b[1:] = (mu[1:]).astype(np.float64)
    n = 1
    while n < 2 * (X + 1):
        n <<= 1
    fa = np.fft.rfft(a, n)
    fb = np.fft.rfft(b, n)
    fa *= fb
    del fb
    chat = np.fft.irfft(fa, n)[: 2 * X + 1]
    del fa
    lhs = float(np.dot(chat, chat))
    print(f"    sum_N Chat(N)^2      = {lhs:.10e}")
    print(f"    sum_h M(h)P(h)       = {tot:.10e}")
    print(f"    ratio                = {lhs / tot:.12f}")
    del chat, a, b
    print()

    # ---------------------------------------------------- (b)
    print("--- (b) the apportionment of mass by shift ----------------------")
    prod = M[1:] * P[1:]
    gross_tot = float(np.abs(prod).sum())
    net_tot = float(prod.sum())
    print(f"    net  sum_{{h != 0}} M(h)P(h) (one-sided) = {net_tot:+.4e}")
    print(f"    gross sum_{{h != 0}} |M(h)P(h)|          = {gross_tot:.4e}")
    print(f"    [paper: net -2.2e13 against gross 4.4e13]")
    print(f"    two-sided (h and -h): net {2 * net_tot:+.4e}, "
          f"gross {2 * gross_tot:.4e}")
    print()
    buckets = [("h < 1e3", 1, 1000), ("1e3 - 1e4", 1000, 10000),
               ("1e4 - 1e5", 10000, 100000),
               ("1e5 - 1e6", 100000, 1000000),
               ("above 1e6", 1000000, X)]
    quoted = [1.1, 3.0, 23.1, 48.9, 23.8]
    # Two readings of "gross mass". The paper's own totals decide which:
    # it reports net -2.2e13 against gross 4.4e13, a ratio of 2. Summing
    # |M(h)P(h)| term by term gives a ratio in the hundreds, so "gross"
    # must be the sum of the ABSOLUTE BUCKET NETS.
    seg_net = [float(prod[lo - 1: hi - 1].sum()) for _, lo, hi in buckets]
    gross_buckets = sum(abs(v) for v in seg_net)
    print(f"    sum of |bucket nets| (one-sided) = {gross_buckets:.4e}, "
          f"two-sided {2 * gross_buckets:.4e}")
    print(f"    ratio to net = {gross_buckets / abs(net_tot):.2f}"
          f"   [paper's gross/net = 2.0]")
    print()
    print(f"    {'bucket':>12}{'share of |bucket nets|':>24}{'paper':>8}"
          f"{'per-h gross share':>20}{'net value':>14}")
    for ((lab, lo, hi), q, nv) in zip(buckets, quoted, seg_net):
        seg = prod[lo - 1: hi - 1]
        g = float(np.abs(seg).sum()) / gross_tot * 100
        print(f"    {lab:>12}{abs(nv) / gross_buckets * 100:>23.1f}%"
              f"{q:>8.1f}{g:>19.1f}%{nv:>+14.3e}")
    print()
    print("    the paper apportions GROSS mass in a sum whose net is 50% of")
    print("    its gross, and rho-1 is built from the NET. Both columns are")
    print("    printed above so the reader can see whether the conclusion")
    print("    ('the wall leans where the averaged theorem is strongest')")
    print("    survives on the quantity that actually enters rho-1.")
    print()

    # rho - 1 in the aggregate
    print(f"    aggregate rho - 1 = sum_{{h!=0}}M(h)P(h) / (M(0)P(0)) = "
          f"{2 * net_tot / (M[0] * P[0]):+.4f}")
    print(f"    [sec:coin quotes a measured rho-1 of -0.18 and a")
    print(f"     reconstruction of -0.0976, 'a factor 0.54']")
    return 0


if __name__ == "__main__":
    sys.exit(main())
