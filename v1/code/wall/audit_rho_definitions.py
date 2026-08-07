# -*- coding: utf-8 -*-
"""
One name, three definitions: rho (increment 331)

WHY. Hazard 9 was named at increment 330 -- two summaries of one object
are not comparable until each one's weight is stated -- and this
program's own rule is that naming a hazard does not prevent it; only a
check does (#285, #296). So it gets swept, and the first place to look
is the quantity the most claims rest on.

rho = Var C / V appears in Proposition W, in #121's rho = 0.810 at
N ~ 1e8, in increment 309's re-derivation of Proposition E, and in
increment 318's per-cell share. Reading the code, it is computed three
different ways:

    A   mean(C^2) / mean(V)          reaudit_prop_E.py
    B   mean(C^2 / V) = var(C/sqrt V)  lab_mask_share_percell.py
    C   (pi/2) * mean(|C|/sqrt V)^2    verify_all.py's half-normal arm,
                                        which is where 0.810 came from

A is a ratio of means; B is a mean of ratios; C is a mean of square
roots squared. They coincide only if V is constant across the band and
C/sqrt(V) is exactly half-normal, and neither holds: V varies by a
factor 2 within an octave, and C carries the location mask.

Nobody has checked that they agree, and every recorded rho is quoted
without saying which one it is.

PRE-REGISTRATION (fixed before the run).

  (R1) DO THEY AGREE? Compute all three per octave band. RULE: all
       three within 2% of each other in every band. If they are, the
       name is unambiguous and this is a clean bill. If not, the
       program has been quoting different quantities under one name.

  (R2) BY HOW MUCH, AND IN WHICH DIRECTION. Report the spread band by
       band and the systematic ordering, so that any recorded figure
       can be corrected rather than merely doubted.

  (R3) WHICH ONE IS THE RIGHT ONE FOR PROPOSITION W. Proposition W is
       rho - 1 = (1/V) sum_{h != 0} c(h) S(h), which is a statement
       about the RATIO OF TOTALS: sum_N C^2 over sum_N V. That is
       definition A. RULE: none -- this is a reading of the
       proposition, and it is stated so it can be disputed.

  WHAT WOULD REFUTE. (R1) passing would mean the three names one
  quantity and nothing needs correcting.
"""
import math
import sys
import time

import numpy as np

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

# The half-normal relation is E|Z| = sqrt(2 rho / pi), so
# rho = (pi/2) (E|Z|)^2. The first draft multiplied by 2/pi
# instead -- the reciprocal -- and produced a 172% 'disagreement'
# that was entirely the constant, off by pi^2/4 = 2.467.
HN = math.pi / 2.0


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
    Ns = np.arange(lo, X + 1, 2)
    print(f"sieve  t={time.time()-t0:.0f}s", flush=True)

    print(f"\n(R1)(R2) three definitions of rho, per octave")
    print(f"{'band':>21} {'A ratio-of-means':>17} "
          f"{'B mean-of-ratios':>17} {'C half-normal':>14} "
          f"{'spread':>8}")
    okR1 = True
    rows = []
    b = lo
    while b < X:
        hi = min(2 * b, X)
        sel = (Ns >= b) & (Ns < hi)
        if int(sel.sum()) > 1000:
            Nb = Ns[sel]
            c, v = C[Nb], V[Nb]
            A = float((c * c).mean() / v.mean())
            B = float((c * c / v).mean())
            Cc = float(HN * (np.abs(c) / np.sqrt(v)).mean() ** 2)
            sp = max(A, B, Cc) / min(A, B, Cc) - 1.0
            okR1 &= sp <= 0.02
            rows.append((b, hi, A, B, Cc, sp))
            print(f"{b:>9}-{hi:>11} {A:>17.5f} {B:>17.5f} "
                  f"{Cc:>14.5f} {sp:>8.2%}")
        b = hi

    sp = np.array([r[5] for r in rows])
    As = np.array([r[2] for r in rows])
    Bs = np.array([r[3] for r in rows])
    Cs = np.array([r[4] for r in rows])
    print(f"\n    (R1) all three within 2% in every band: "
          f"{'PASS' if okR1 else 'FAIL'}  "
          f"(max spread {sp.max():.2%})")
    print(f"    (R2) ordering, averaged: A {As.mean():.4f}, "
          f"B {Bs.mean():.4f}, C {Cs.mean():.4f}")
    print(f"         B/A = {float((Bs/As).mean()):.4f}   "
          f"C/A = {float((Cs/As).mean()):.4f}   "
          f"C/B = {float((Cs/Bs).mean()):.4f}")
    print(f"\n    (R3) Proposition W states rho - 1 = (1/V) sum_h c(h)S(h),")
    print(f"         which is a ratio of TOTALS -- definition A. The")
    print(f"         0.810 of #121 came from the half-normal arm, i.e.")
    print(f"         definition C, and increment 318's per-cell share")
    print(f"         used definition B.")

    if okR1:
        v = ("the three definitions agree to 2% everywhere, so rho "
             "names one quantity and no recorded figure needs "
             "correcting")
    else:
        v = (f"the three differ by up to {sp.max():.1%}, so rho has "
             f"named up to three quantities. C sits "
             f"{abs(float((Cs/As).mean())-1):.1%} from A "
             f"systematically, which is the half-normal assumption "
             f"failing on a field that carries a mask. Every recorded "
             f"rho must say which definition it is: Proposition W "
             f"wants A, #121's 0.810 is C, #318's share is B")
    print(f"\n    {v}")
    print("DONE")


if __name__ == "__main__":
    main()
