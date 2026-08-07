# -*- coding: utf-8 -*-
"""
How much is left in the uncertainty bound? (increment 321)

WHY. Asked whether quantum uncertainty and entanglement could help
with the wall, the honest answer is that both already appear in this
program's results, in their exact mathematical forms, and both are on
the wall's side:

  UNCERTAINTY is Parseval. Heisenberg's inequality is a statement about
  a Fourier pair, and the Fourier pair here gives
      sup_a |S_mu(a)| >= ||S_mu||_2 = sqrt(6N/pi^2),
  which is Proposition E's route (ii) -- the inequality that CLOSES the
  pointwise route rather than opening it. THEOREM_A.md states the
  consequence in its own words: "no improvement in Mobius
  exponential-sum technology can push sup_a|S_mu| below
  (6/pi^2)^{1/2} N^{1/2}".

  ENTANGLEMENT has a precise counterpart. C(N) = sum_v mu(v)Lambda(N-v)
  is a correlation of two arithmetic functions at a FIXED SUM, and if
  they were independent -- unentangled -- the cancellation ratio
  rho = Var C / V would be exactly 1. Proposition W identifies the
  deficit: rho - 1 = (1/V) sum_{h != 0} c(h) S(h), a prime-pair-weighted
  Chowla correlation. So rho - 1 IS this program's measure of the
  non-factorisation, and it is measured: 0.810 +/- 0.018 at N ~ 1e8.

  But increment 320 showed rho and A(N) CANCEL out of the comparison
  that decides whether a method reaches the target. A complete theory
  of the entanglement therefore moves nothing.

THE ONE THING NOT YET MEASURED. Route (ii) needs sup_a|S_mu| small, and
Parseval floors it at ||S_mu||_2. How much slack is there between the
actual sup and its floor? If the sup already sits at the floor, the
route is closed tightly and no technology can move it. If there is
slack, the question is how much saturating it would buy -- and the
target needs a GROWING log power, so a bounded gain buys nothing
whatever its size.

PRE-REGISTRATION (fixed before the run).

  (X1) THE PARSEVAL FLOOR IS EXACT, a self-test: ||S_mu||_2/sqrt(N)
       must reproduce sqrt(6/pi^2) = 0.7797 to three decimals at every
       size. It is an identity, not a measurement.

  (X2) THE SLACK IS BOUNDED. sup_a|S_mu| / ||S_mu||_2 measured across
       N = 2^14 .. 2^22. RULE: the ratio does not grow with N -- fitted
       slope of log(ratio) against log N within 0.05 of zero. If it
       grew, saturating the floor would buy a growing factor and the
       route would deserve another look.

  (X3) WHAT SATURATION WOULD BUY, in the units of increment 320: the
       best conceivable pointwise bound is
       ||S_mu||_2 * ||S_Lambda||_1, and it is compared against the
       target N(log N)^{-A}. Reported, not thresholded.

  WHAT WOULD REFUTE. (X2) failing -- a growing ratio -- would mean the
  gap between the sup and its floor is itself a resource, and the
  pointwise route is not closed as tightly as Proposition E says.
"""
import math
import sys
import time

import numpy as np

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

C2PI = 6.0 / math.pi ** 2


def sieve(X):
    mu = np.ones(X + 1, dtype=np.int8)
    pm = np.ones(X + 1, dtype=bool)
    pm[:2] = False
    for p in range(2, int(X ** 0.5) + 1):
        if pm[p]:
            pm[p * p::p] = False
            mu[p::p] *= -1
            mu[p * p::p * p] = 0
    val = np.arange(X + 1, dtype=np.int64)
    for p in range(2, int(X ** 0.5) + 1):
        if pm[p]:
            val[p::p] //= p
            pp = p * p
            while pp <= X:
                val[pp::pp] //= p
                pp *= p
    mu[val > 1] *= -1
    lam = np.zeros(X + 1, dtype=np.float64)
    for p in range(2, X + 1):
        if pm[p]:
            q = p
            lg = math.log(p)
            while q <= X:
                lam[q] = lg
                q *= p
    return mu, lam


def main():
    t0 = time.time()
    XM = 1 << 22
    mu, lam = sieve(XM)
    print(f"sieve  t={time.time()-t0:.0f}s", flush=True)

    print(f"\n(X1)(X2) the sup against its Parseval floor")
    print(f"{'N':>10} {'||S_mu||2/sqrtN':>16} {'sup|S_mu|/sqrtN':>16} "
          f"{'slack':>8} {'||S_L||1/sqrtN':>15}")
    okX1 = True
    rows = []
    for e in range(14, 23):
        N = 1 << e
        g = 4 * N
        sm = np.fft.fft(mu[:N].astype(np.float64), g)
        sl = np.fft.fft(lam[:N], g)
        l2 = float(np.sqrt((np.abs(sm) ** 2).mean()))
        sup = float(np.abs(sm).max())
        l1lam = float(np.abs(sl).mean())
        r2 = l2 / math.sqrt(N)
        okX1 &= abs(r2 - math.sqrt(C2PI)) < 5e-3
        rows.append((N, r2, sup / math.sqrt(N), sup / l2,
                     l1lam / math.sqrt(N)))
        print(f"{N:>10} {r2:>16.4f} {sup/math.sqrt(N):>16.4f} "
              f"{sup/l2:>8.3f} {l1lam/math.sqrt(N):>15.4f}",
              flush=True)

    Nn = np.array([r[0] for r in rows], dtype=float)
    sl_ = np.array([r[3] for r in rows])
    slope = float(np.polyfit(np.log(Nn), np.log(sl_), 1)[0])
    okX2 = abs(slope) <= 0.05
    print(f"\n    (X1) ||S_mu||2/sqrt(N) = sqrt(6/pi^2) = "
          f"{math.sqrt(C2PI):.4f} at every size: "
          f"{'PASS' if okX1 else 'FAIL'}")
    print(f"    (X2) the slack sup/||.||_2 does not grow with N: "
          f"{'PASS' if okX2 else 'FAIL'}  "
          f"(fitted slope {slope:+.4f}, ratio "
          f"{sl_.min():.2f} to {sl_.max():.2f})")

    print(f"\n(X3) what saturating the uncertainty bound would buy")
    N, r2, supr, slack, l1 = rows[-1]
    best = r2 * math.sqrt(N) * l1 * math.sqrt(N)
    actual = supr * math.sqrt(N) * l1 * math.sqrt(N)
    lg = math.log(N)
    print(f"    at N = {N}:")
    print(f"      actual   sup|S_mu| * ||S_L||_1  = {actual:>12.4e}"
          f"   = {actual/N:>7.3f} N")
    print(f"      best possible ||S_mu||_2 * ||S_L||_1 "
          f"= {best:>12.4e}   = {best/N:>7.3f} N")
    print(f"      gain from full saturation        "
          f"{actual/best:>7.3f}x  -- a BOUNDED constant")
    for A in (1, 2):
        need = lg ** A
        print(f"      target needs a saving of (log N)^{A} = "
              f"{need:>9.2f} over trivial;")
        print(f"        the best pointwise product is still "
              f"{best/N:.3f}x trivial, i.e. it SAVES nothing")
    if okX1 and okX2:
        v = ("the uncertainty bound is already the binding constraint "
             f"and its slack is a bounded {sl_.mean():.1f}x that does "
             f"not grow. Saturating it exactly would improve the "
             f"pointwise product by {actual/best:.1f}x and leave it "
             f"still ABOVE trivial, while the target needs a growing "
             f"log power. Uncertainty is on the wall's side, and "
             f"entanglement -- rho - 1, Proposition W -- cancels out of "
             f"the comparison entirely (inc. 320). Neither is a "
             f"resource here")
    elif okX1:
        v = ("the slack grows with N, so the gap between the sup and "
             "its Parseval floor is itself a resource and Proposition "
             "E's route (ii) is not closed as tightly as recorded")
    else:
        v = ("the Parseval floor does not reproduce; nothing here reads")
    print(f"\n    {v}")
    print("DONE")


if __name__ == "__main__":
    main()
