# -*- coding: utf-8 -*-
"""
Transform Lab, session 5 (increment 224): which pairings admit a
simultaneous dilation, and which do not.

Tao's reduction and Helfgott-Radziwill's expansion both rest on one
mechanism: the correlation is invariant under dilating BOTH factors at
once. For lambda,

    lambda(pn) lambda(pn + p) = lambda(p)^2 lambda(n) lambda(n+1)
                              = lambda(n) lambda(n+1),

so the correlation at scale pn is the correlation at scale n. That is
the self-similarity the whole machine runs on.

Our wall is C(N) = Sum_v mu(v) Lambda(N-v): a pairing of TWO DIFFERENT
functions, only one of which is completely multiplicative. This script
checks three pairings for a simultaneous dilation, exactly:

 (1) lambda-lambda, additive shift    : the mechanism, as a control
 (2) mu-Lambda, additive constraint N : our wall
 (3) mu-mu, additive constraint N     : the wall with Lambda replaced,
     where BOTH factors are completely multiplicative -- and where the
     dilation turns out to relate different N rather than acting
     within one N

NULLS AND CRITERIA. Each test is an exact identity check, so the null
is zero mismatches; a pairing either admits the dilation or does not.
For (3) the identity checked is
    mu(pv) mu(pN' - pv) = mu(v) mu(N' - v)
for p dividing neither v nor N' - v -- BOTH side conditions are
needed, since each factor separately requires one of them. An earlier
run enforced only the first and failed on 5.2% of draws, exactly the
density of p | (N' - v); a test fault, recorded. The identity
relates the correlation at N = pN' to the one at N'.
"""
import numpy as np
import math


def mobius_upto(X):
    mu = np.ones(X + 1, dtype=np.int8)
    pm = np.ones(X + 1, dtype=bool); pm[:2] = False
    for p in range(2, int(X ** 0.5) + 1):
        if pm[p]:
            pm[p * p::p] = False
            mu[p::p] *= -1
            mu[p * p::p * p] = 0
    val = np.arange(X + 1, dtype=np.int64)
    for p in range(2, int(X ** 0.5) + 1):
        if pm[p]:
            val[p::p] //= p
    mu[val > 1] *= -1
    mu[0] = 0
    return mu


def lam_upto(X):
    sieve = np.ones(X + 1, dtype=bool); sieve[:2] = False
    for i in range(2, int(X ** 0.5) + 1):
        if sieve[i]:
            sieve[i * i::i] = False
    lam = np.zeros(X + 1)
    for p in np.nonzero(sieve)[0]:
        q = int(p); lp = math.log(int(p))
        while q <= X:
            lam[q] = lp; q *= int(p)
    return lam


def liouville_upto(X):
    om = np.zeros(X + 1, dtype=np.int32)
    for p in range(2, X + 1):
        if om[p] == 0:                      # p prime
            q = p
            while q <= X:
                om[q::q] += 1
                q *= p
                if q > X:
                    break
    return np.where(om % 2 == 0, 1, -1)


def main():
    X = 300_000
    mu = mobius_upto(X)
    lam = lam_upto(X)
    lio = liouville_upto(X)
    rng = np.random.default_rng(20260910)
    PS = [2, 3, 5, 7, 11, 13]

    print("(1) control -- lambda(pn)lambda(pn+p) = lambda(n)lambda(n+1)")
    bad = t = 0
    for _ in range(200000):
        p = int(rng.choice(PS)); n = int(rng.integers(1, 20000))
        if p * n + p > X:
            continue
        t += 1
        if int(lio[p * n]) * int(lio[p * n + p]) != \
           int(lio[n]) * int(lio[n + 1]):
            bad += 1
    print(f"    tested {t}, mismatches {bad}   "
          f"{'HOLDS -- the mechanism' if bad == 0 else 'FAILS'}")

    print("\n(2) our wall -- is mu(pv)Lambda(N-pv) tied to "
          "mu(v)Lambda(N'-v) for any N'?")
    print("    Lambda is not completely multiplicative: Lambda(pm) is")
    print("    log p when m is a power of p and 0 otherwise, so no")
    print("    factor of the form Lambda(p) x Lambda(m) exists. The")
    print("    check below reports how often Lambda(p*m) equals")
    print("    anything proportional to Lambda(m).")
    nz = prop = 0
    for _ in range(200000):
        p = int(rng.choice(PS)); m = int(rng.integers(2, 20000))
        if p * m > X:
            continue
        if lam[m] != 0.0:
            nz += 1
            if lam[p * m] != 0.0:
                prop += 1
    print(f"    Lambda(m) nonzero in {nz} draws; Lambda(pm) also "
          f"nonzero in {prop} ({prop/max(nz,1):.4f})")
    print("    => the second factor does not survive dilation: NO "
          "simultaneous dilation")

    print("\n(3) mu-mu -- mu(pv)mu(pN'-pv) = mu(v)mu(N'-v), p not| v")
    bad = t = 0
    for _ in range(200000):
        p = int(rng.choice(PS))
        Np = int(rng.integers(10, 4000)) * 2
        v = int(rng.integers(1, Np))
        # BOTH side conditions: p not| v AND p not| (N'-v)
        if p * Np > X or v % p == 0 or (Np - v) % p == 0:
            continue
        t += 1
        if int(mu[p * v]) * int(mu[p * (Np - v)]) != \
           int(mu[v]) * int(mu[Np - v]):
            bad += 1
    print(f"    tested {t}, mismatches {bad}   "
          f"{'HOLDS' if bad == 0 else 'FAILS'}")
    print("    but note what it relates: the correlation at N = p N'")
    print("    restricted to p | v, against the correlation at N'.")
    print("    The dilation moves N, so it acts ACROSS the N-family,")
    print("    not within a single N -- which is the N-average again.")

    print("\nConclusion: the Tao/HR mechanism needs both factors to")
    print("transform under the SAME dilation. Our wall pairs mu with")
    print("Lambda and only mu does. Replacing Lambda by mu restores the")
    print("dilation but moves N, landing back on the N-average.")
    print("DONE")


if __name__ == "__main__":
    main()
