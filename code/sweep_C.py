# -*- coding: utf-8 -*-
"""
Sweep C (increment 214): five variants of the wall's scalar, to see
whether any relative of C(N) is smaller than C(N) itself. Nulls stated
with the criteria.

 C1  pure Mobius pair   P(N) = Sum_n mu(n) mu(N-n)
 C2  primes only        Q(N) = Sum_{p<N} log p * mu(N-p)
 C3  mask-only          Z(N) = Sum_n Lambda(n) (mu^2(N-n) - 6/pi^2)
 C4  Mertens calibration M(N) = Sum_{n<N} mu(n)      [known object]
 C5  twisted            T(N) = Sum_n Lambda(n) mu(N-n) * (-1)^n

Each is measured over four octave groups and its exponent fitted:
|X(N)| ~ N^alpha. The reference is C(N) itself at alpha = 0.503.

  ALIVE iff some variant has alpha <= 0.40, i.e. genuinely smaller
        than square-root -- that variant would be a softer target.
  DEAD  iff every alpha sits in [0.45, 0.60], the square-root band.

C4 is a calibration line, not a hypothesis: Mertens is known to be
square-root-ish, so an alpha far from 0.5 there would indict the
measurement rather than the mathematics.
"""
import numpy as np
import math

from thmC_alpha_scan import sieve


def main():
    X = 1_100_000
    mu, lam, phi, primes, spf = sieve(X)
    isprime = np.zeros(X + 1, dtype=bool)
    isprime[primes] = True

    groups = [(140_000, 60), (280_000, 60), (560_000, 60),
              (1_050_000, 60)]
    names = ["C1 mu*mu", "C2 primes", "C3 mask-only",
             "C4 Mertens", "C5 twisted"]
    acc = {k: [] for k in names}
    Ns = []
    for N0, cnt in groups:
        vals = {k: [] for k in names}
        for t in range(cnt):
            N = N0 + 2 * t
            idx = np.arange(1, N)
            mrev = mu[N - idx].astype(np.float64)
            lv = lam[1:N]
            vals["C1 mu*mu"].append(
                abs(float(np.dot(mu[1:N].astype(np.float64), mrev))))
            pv = np.where(isprime[1:N], lv, 0.0)
            vals["C2 primes"].append(abs(float(np.dot(pv, mrev))))
            vals["C3 mask-only"].append(
                abs(float(np.dot(lv, (mrev ** 2) - 6 / math.pi ** 2))))
            vals["C4 Mertens"].append(
                abs(float(mu[1:N].astype(np.float64).sum())))
            sgn = np.where(idx % 2 == 0, 1.0, -1.0)
            vals["C5 twisted"].append(abs(float(np.dot(lv * sgn, mrev))))
        for k in names:
            acc[k].append(float(np.mean(vals[k])))
        Ns.append(float(N0))

    print(f"{'variant':<14} " + " ".join(f"{int(x):>10}" for x in Ns)
          + f" {'alpha':>7}")
    flags = 0
    for k in names:
        a = np.array(acc[k])
        al = np.polyfit(np.log(Ns), np.log(np.maximum(a, 1e-9)), 1)[0]
        f = int(al <= 0.40 and not k.startswith("C4"))
        flags += f
        print(f"{k:<14} " + " ".join(f"{v:>10.1f}" for v in a)
              + f" {al:>7.3f} {'<<<' if f else ''}")
    print(f"\nreference: C(N) itself has alpha = 0.503")
    print(f"SWEEP C: {flags} variants below alpha 0.40")
    print("DONE")


if __name__ == "__main__":
    main()
