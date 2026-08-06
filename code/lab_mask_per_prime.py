# -*- coding: utf-8 -*-
"""
The location mask's per-prime factor (increment 243): can the table
be turned into a formula?

Increment 242 showed the mask compounds MULTIPLICATIVELY across the
primes dividing N -- the 9699690 cell sits 1.26 to 1.35 times deeper
than the 510510 cell, for the single added prime 19 -- and that is why
an additive one-coefficient-per-prime model missed the primorials by a
factor 18 while full modular enumeration worked. Enumeration captures
the structure without knowing its form. This asks for the form.

WHAT IS MEASURED. For a fixed core M and a prime q not dividing M,
compare the deterministic term over the family N = k*M split by whether
q divides k:

    f(q; M) := mean( C(N)/sqrt(N) | q | k ) / mean( C(N)/sqrt(N) | q not| k )

Everything else is matched: same core, same N range, same scale. If the
mask is multiplicative, f(q; M) does not depend on M, and testing that
IS the test of multiplicativity. Two cores are used, M = 30030
(2*3*5*7*11*13) and M = 210 (2*3*5*7).

THE PREDICTION, made before running. The mechanism is that n prime
forces v = N - n coprime to rad(N), so v is drawn from the integers
with no prime factor in rad(N). For a Mobius sum restricted that way,

    Sum_{(v,q)=1} mu(v) v^{-s} = zeta(s)^{-1} / (1 - q^{-s}),

so removing q from the support multiplies the sum by 1/(1 - 1/q) =
q/(q-1) at s = 1. PREDICTED f(q) = q/(q-1): 1.063 at q = 17, 1.056 at
19, 1.045 at 23, 1.036 at 29. The competing candidate, from the
singular series rather than the Mobius sum, is (q-1)/(q-2): 1.067,
1.059, 1.048, 1.037 -- numerically almost the same, so this experiment
can confirm the SCALE of the effect but not distinguish those two.

  CONFIRMED  iff f(q) sits within its error bar of q/(q-1), and is
             the same across the two cores.
  REFUTED    iff f(q) is far from q/(q-1) -- in which case the ratio
             itself is the new datum, and the Euler-factor reading of
             the mechanism is wrong.

NULLS. Each mean carries the standard error sd/sqrt(n) of its own
subgroup, and the ratio's error is propagated from both. A control
prime is included: q chosen LARGER than any that should matter, where
the predicted factor is ~1 and any measured departure is noise, which
calibrates how much apparent structure this design manufactures.
"""
import numpy as np
import math
import sys
import time

from lab_location_mask_large import sieve_mu_and_primepowers, CV


def group_stat(Ns, mu, primes, lpr):
    vals = []
    for N in Ns:
        C, V = CV(int(N), mu, primes, lpr)
        vals.append(C / math.sqrt(N))
    a = np.array(vals)
    return float(a.mean()), float(a.std(ddof=1) / math.sqrt(len(a))), \
        len(a)


def main():
    X = int(sys.argv[1]) if len(sys.argv) > 1 else 30_000_000
    t0 = time.time()
    mu, primes = sieve_mu_and_primepowers(X)
    lpr = np.log(primes.astype(np.float64))
    print(f"sieve to {X}: t={time.time()-t0:.0f}s", flush=True)

    QS = [17, 19, 23, 29, 31, 37, 101]       # 101 is the control prime
    CORES = [("M = 30030 (2*3*5*7*11*13)", 30030),
             ("M = 210   (2*3*5*7)", 210)]

    for label, M in CORES:
        kmax = X // M
        ks = np.arange(1, kmax + 1)
        if len(ks) > 1100:                    # keep the runtime finite
            ks = ks[:: len(ks) // 1100][:1100]
        print(f"\n=== {label}   k up to {kmax}, using {len(ks)} of them")
        print(f"{'q':>5} {'n(q|k)':>7} {'mean|q':>9} {'SE':>7} "
              f"{'n(q!|k)':>8} {'mean|not':>9} {'SE':>7} "
              f"{'f(q)':>7} {'SE':>6} {'q/(q-1)':>8} {'z vs pred':>10}")
        for q in QS:
            A = [int(k) * M for k in ks if int(k) % q == 0]
            B = [int(k) * M for k in ks if int(k) % q != 0]
            if len(A) < 12:
                continue
            if len(B) > 700:
                B = B[:: len(B) // 700][:700]
            ma, sa, na = group_stat(A, mu, primes, lpr)
            mb, sb, nb = group_stat(B, mu, primes, lpr)
            f = ma / mb
            sf = abs(f) * math.sqrt((sa / ma) ** 2 + (sb / mb) ** 2)
            pred = q / (q - 1.0)
            print(f"{q:>5} {na:>7} {ma:>9.3f} {sa:>7.3f} "
                  f"{nb:>8} {mb:>9.3f} {sb:>7.3f} "
                  f"{f:>7.3f} {sf:>6.3f} {pred:>8.3f} "
                  f"{(f-pred)/sf:>10.2f}")
            print(f"      t={time.time()-t0:.0f}s", flush=True)

    print("\nreading: q = 101 is the control -- its predicted factor is")
    print("1.010, so whatever f(101) reads is the design's own noise.")
    print("If f(17..37) sit near q/(q-1) in BOTH cores, the mask is")
    print("multiplicative with that Euler factor and the table becomes")
    print("a formula. If they sit well above it, the Euler-factor")
    print("reading of the mechanism is wrong and the measured ratios")
    print("are the datum to explain.")
    print("DONE")


if __name__ == "__main__":
    main()
