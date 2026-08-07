# -*- coding: utf-8 -*-
"""
prop:Dpp -- the two pieces of CP_2, and the mean-zero tuning.
(v1_verify2, Phase 1, blind.)

STATEMENTS UNDER TEST, verbatim:

  "Measured at N=10^6, 4*10^6, 1.6*10^7: the two pieces of CP_2 are the
   same order, with ratio 0.771, 0.790, 0.810 drifting towards 1, and
   CP_2/(N log N) = 2.886, 2.949, 2.997. The closure rests on
   nonnegativity, not on one piece dominating. The one analytically
   canonical tuning, f(x)=x^2-2 gamma x, which makes b mean-zero by
   killing the pole of zeta W at s=1, moves the complete part by about
   five percent."

DERIVATION USED (from the statements, not from v1's code).  The complete
part is sum_{u<N} Lambda(N-u) mu^2(u) b_u with b = mu * w. For the
monomial w_k = log^2 k, b = Lambda_2 = Lambda*log + Lambda*Lambda. The
factor mu^2(u) restricts to squarefree u, and Lambda_2 is supported on
omega(u) <= 2, so exactly two pieces survive:

   r = 1 : u = p prime,  Lambda_2(p) = (log p)^2
   r = 2 : u = p q, p != q,  Lambda_2(pq) = 2 log p log q

so  CP_2 = sum_{p<N} (log p)^2 Lambda(N-p)
          + 2 sum_{pq<N, p<q} log p log q Lambda(N-pq).

For f(x) = x^2 - 2 gamma x one has b = Lambda_2 - 2 gamma Lambda, so the
tuned complete part is CP_2 - 2 gamma CP_1 with
CP_1 = sum_{p<N} log p Lambda(N-p), the binary Goldbach sum.

PRE-REGISTRATION.

  Decision rule. Compute both pieces at the three stated N.
    REPRODUCED iff the piece ratio gives 0.771/0.790/0.810 and
    CP_2/(N log N) gives 2.886/2.949/2.997.
    Report which orientation of the ratio (r=1 over r=2, or the reverse)
    the paper's numbers correspond to, since it does not say.
  Also report the tuned shift as a fraction, against the quoted "about
  five percent".

  Prediction written before running.  REPRODUCED. This is an
  unambiguous finite sum with no threshold and no null, the kind of
  quantity that has reproduced everywhere in this pass so far. I predict
  the "drifting towards 1" claim is an extrapolation the data cannot
  carry: three points moving by ~0.02 per 4x in N reach 1 only near
  N ~ 1e13, and I predict the increments do not accelerate.
"""

import sys

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

import numpy as np

EULER_GAMMA = 0.5772156649015329


def sieve(X):
    s = np.ones(X + 1, dtype=bool)
    s[:2] = False
    for i in range(2, int(X ** 0.5) + 1):
        if s[i]:
            s[i * i:: i] = False
    primes = np.nonzero(s)[0]
    lam = np.zeros(X + 1, dtype=np.float64)
    lam[primes] = np.log(primes.astype(np.float64))
    for p in primes[primes * primes <= X]:
        q = int(p) * int(p)
        lp = float(np.log(p))
        while q <= X:
            lam[q] = lp
            q *= int(p)
    return lam, primes


def main():
    print("audit_propDpp_cp2   (v1_verify2 Phase 1, blind)")
    print("=" * 74)
    Xmax = 16_000_000
    print(f"  sieving to {Xmax:,} ...")
    lam, primes = sieve(Xmax)
    logp = np.log(primes.astype(np.float64))

    print()
    print(f"  {'N':>12}{'piece r=1':>16}{'piece r=2':>16}{'r1/r2':>9}"
          f"{'r2/r1':>9}{'CP_2':>16}{'CP_2/(N logN)':>15}")
    rows = []
    for N in (1_000_000, 4_000_000, 16_000_000):
        sel = primes < N
        p = primes[sel]
        lp = logp[sel]
        # r = 1 : u = p
        piece1 = float(np.dot(lp * lp, lam[N - p]))
        # r = 2 : u = p*q, p < q
        piece2 = 0.0
        root = int(N ** 0.5)
        for i, pp in enumerate(primes[primes < root]):
            pp = int(pp)
            hi = N // pp
            q = primes[(primes > pp) & (primes < hi)]
            if len(q) == 0:
                continue
            u = pp * q
            piece2 += 2.0 * float(np.log(pp)) * float(
                np.dot(np.log(q.astype(np.float64)), lam[N - u]))
        cp2 = piece1 + piece2
        cp1 = float(np.dot(lp, lam[N - p]))
        rows.append((N, piece1, piece2, cp2, cp1))
        print(f"  {N:>12,}{piece1:>16.4e}{piece2:>16.4e}"
              f"{piece1 / piece2:>9.4f}{piece2 / piece1:>9.4f}"
              f"{cp2:>16.4e}{cp2 / (N * np.log(N)):>15.4f}")
    print(f"  [paper: ratio 0.771, 0.790, 0.810; "
          f"CP_2/(N log N) = 2.886, 2.949, 2.997]")

    print()
    print("  the mean-zero tuning f(x) = x^2 - 2 gamma x:")
    print(f"  {'N':>12}{'CP_2':>16}{'CP_1':>16}{'tuned':>16}{'shift':>10}")
    for N, p1, p2, cp2, cp1 in rows:
        tuned = cp2 - 2 * EULER_GAMMA * cp1
        print(f"  {N:>12,}{cp2:>16.4e}{cp1:>16.4e}{tuned:>16.4e}"
              f"{(tuned - cp2) / cp2 * 100:>9.2f}%")
    print(f"  [paper: 'moves the complete part by about five percent']")

    print()
    print("  'drifting towards 1': the increments, and where they arrive")
    r = [x[1] / x[2] for x in rows]
    print(f"    ratios      : {r[0]:.4f}, {r[1]:.4f}, {r[2]:.4f}")
    print(f"    increments  : {r[1] - r[0]:+.4f}, {r[2] - r[1]:+.4f} "
          f"per 4x in N")
    if r[2] > r[1]:
        need = (1.0 - r[2]) / max(r[2] - r[1], 1e-12)
        print(f"    at the last observed rate, reaching 1 needs "
              f"{need:.1f} further 4x steps,")
        print(f"    i.e. N ~ 1e{np.log10(1.6e7 * 4 ** need):.0f}")
    print(f"    all measurements in this paper are at N <= 1.6e7, "
          f"with two arms at 1e8.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
