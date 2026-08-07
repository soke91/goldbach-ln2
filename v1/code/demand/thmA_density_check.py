# -*- coding: utf-8 -*-
"""
Theorem A load-bearing check (increment 191): exact rational
verification of the main-term density identity

    Sum_{g | m} mu(g) / ( phi(m/g) * g * phi(g) )  =  1/m      (m squarefree)

This single line carries Theorem A: density exponent exactly 1 gives
PNT-strength cancellation of Sum_m mu(m) lambda(m)/m; any other
exponent would leave only (log M)^{-c} (Selberg-Delange) and
Theorem A would be false.
"""
from fractions import Fraction
from sympy import mobius, totient, divisors, factorint


def squarefree(n):
    return all(e == 1 for e in factorint(n).values())


def main(limit=400):
    bad = []
    tested = 0
    for m in range(1, limit):
        if m > 1 and not squarefree(m):
            continue
        tested += 1
        s = Fraction(0)
        for g in divisors(m):
            mu = int(mobius(g))
            if mu == 0:
                continue
            s += Fraction(mu, 1) / (int(totient(m // g)) * g * int(totient(g)))
        if s != Fraction(1, m):
            bad.append((m, s, Fraction(1, m)))
    print(f"squarefree m tested: {tested} (m < {limit})")
    print(f"mismatches: {len(bad)}")
    if bad:
        print(bad[:5])
    else:
        print("identity holds exactly (rational arithmetic)")


if __name__ == "__main__":
    main()
