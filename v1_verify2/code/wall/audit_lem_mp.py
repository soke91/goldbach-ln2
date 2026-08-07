# -*- coding: utf-8 -*-
"""
Is Lemma MP (lem:MP) an identity?  (v1_verify2, Phase 1, blind)

THE STATEMENT UNDER TEST, verbatim from v1/paper/wall_v1.tex:

    Lemma (the second-moment identity).  For any X,
        sum_{N<=X} C(N)^2  =  sum_{|h|<X} M(h) P(h),
        M(h) = sum_v mu(v) mu(v+h),   P(h) = sum_w Lambda(w) Lambda(w+h),
    the inner sums being over the ranges that keep both arguments below X.

    "The identity is exact and unconditional."

WHY IT IS SUSPECT.  Expanding C(N)^2 as a double sum over (n,v),(n',v')
with n+v = n'+v' = N and summing over N <= X leaves three free indices
(n, v, h) subject to n + v = N <= X.  The Lambda-pair carries (n, n-h)
and the mu-pair carries (v, v+h), but the two are still tied together by
the SIMPLEX constraint n + v <= X.  Factoring into M(h) * P(h) replaces
that simplex by the BOX {n < X} x {v < X}.  The proof in the paper does
the collection step and then says "collecting the Lambda-pairs and the
mu-pairs separately gives the two factors" -- that separation is the
step that is not available.

PRE-REGISTRATION (fixed before this ran).

  Decision rule.  Compute both sides exactly, in integer/exact float
  arithmetic, by brute force, at X = 200, 400, 800, 1600, 3200.
    PASS  : |RHS/LHS - 1| < 1e-9 at every X.  Lemma MP is an identity.
    FAIL  : the ratio is bounded away from 1 at every X.  Then report
            the ratio and whether it is stable in X (a stable ratio
            means a fixed missing factor; a drifting one means the
            defect is not a constant).

  Prediction written before running.  FAIL, with RHS/LHS in the
  neighbourhood of 2, because the simplex n+v <= X is (to first order)
  half of the box {n<X} x {v<X} and the summands do not conspire.
  I expect the ratio to drift slowly rather than sit exactly at 2,
  because Lambda and mu are not flat across the range.

  What would refute the finding.  RHS/LHS -> 1 as X grows, or a reading
  of "the ranges that keep both arguments below X" under which the two
  sides agree.  The second is tested explicitly: three candidate
  readings of the range convention are all evaluated.
"""

import sys

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

import numpy as np


def mobius(n):
    """mu(0..n) as int8."""
    mu = np.ones(n + 1, dtype=np.int8)
    primes = []
    is_c = np.zeros(n + 1, dtype=bool)
    for i in range(2, n + 1):
        if not is_c[i]:
            primes.append(i)
            mu[i] = -1
        for p in primes:
            if i * p > n:
                break
            is_c[i * p] = True
            if i % p == 0:
                mu[i * p] = 0
                break
            mu[i * p] = -mu[i]
    mu[0] = 0
    return mu


def vonmangoldt(n):
    """Lambda(0..n) as float64."""
    lam = np.zeros(n + 1, dtype=np.float64)
    sieve = np.ones(n + 1, dtype=bool)
    sieve[:2] = False
    for i in range(2, int(n ** 0.5) + 1):
        if sieve[i]:
            sieve[i * i:: i] = False
    for p in np.nonzero(sieve)[0]:
        q = int(p)
        while q <= n:
            lam[q] = np.log(p)
            q *= int(p)
    return lam


def lhs_exact(X, lam, mu):
    """sum_{N<=X} C(N)^2 with C(N) = sum_{n<N} Lambda(n) mu(N-n)."""
    tot = 0.0
    per_N = np.zeros(X + 1)
    for N in range(2, X + 1):
        n = np.arange(1, N)
        c = float(np.dot(lam[1:N], mu[N - 1:0:-1].astype(np.float64)))
        per_N[N] = c
        tot += c * c
    return tot, per_N


def rhs_box(X, lam, mu, strict):
    """sum_{|h|<X} M(h) P(h), both arguments constrained by `strict`.

    strict = 'lt'  : arguments < X   (the paper's words, read literally)
    strict = 'le'  : arguments <= X  (the other reading of "below X")
    """
    top = X - 1 if strict == "lt" else X
    lamf = lam[: top + 1].astype(np.float64)
    muf = mu[: top + 1].astype(np.float64)
    tot = 0.0
    for h in range(-(X - 1), X):
        if h >= 0:
            a = lamf[1: top + 1 - h]
            b = lamf[1 + h: top + 1]
            c = muf[1: top + 1 - h]
            d = muf[1 + h: top + 1]
        else:
            g = -h
            a = lamf[1 + g: top + 1]
            b = lamf[1: top + 1 - g]
            c = muf[1 + g: top + 1]
            d = muf[1: top + 1 - g]
        P = float(np.dot(a, b))
        M = float(np.dot(c, d))
        tot += M * P
    return tot


def rhs_simplex(X, lam, mu):
    """The identity that IS true: keep the coupling n + v <= X.

    sum_{N<=X} C(N)^2 = sum_h sum_{n+v<=X, n-h>=1, v+h>=1}
                            Lambda(n)Lambda(n-h) mu(v)mu(v+h)
    Evaluated directly, as a check on the derivation itself.
    """
    lamf = lam.astype(np.float64)
    muf = mu.astype(np.float64)
    tot = 0.0
    for h in range(-(X - 1), X):
        s = 0.0
        for n in range(max(1, 1 + h), X):
            w = lamf[n] * lamf[n - h]
            if w == 0.0:
                continue
            vmax = X - n
            vlo = max(1, 1 - h)
            if vlo > vmax:
                continue
            v = np.arange(vlo, vmax + 1)
            s += w * float(np.dot(muf[v], muf[v + h]))
        tot += s
    return tot


def repair(X, lam, mu):
    """The identity that IS available with a factorized right-hand side.

    Let Chat(N) = sum_{n+v=N, n<=X, v<=X} Lambda(n) mu(v)  --  the
    TRUNCATED convolution.  Then, summing over ALL N (which runs to 2X,
    not to X),

        sum_{N<=2X} Chat(N)^2 = sum_{|h|<X} M(h) P(h)

    with M, P over arguments <= X.  The box is now the honest domain on
    both sides.  What cannot be done is to keep sum_{N<=X} C(N)^2 on the
    left.
    """
    lamf = lam[: X + 1].astype(np.float64)
    muf = mu[: X + 1].astype(np.float64)
    a = np.zeros(X + 1)
    a[1:] = lamf[1:]
    b = np.zeros(X + 1)
    b[1:] = muf[1:]
    chat = np.convolve(a, b)
    return float(np.dot(chat, chat))


def main():
    print("audit_lem_mp -- is Lemma MP an identity?  (v1_verify2 Phase 1)")
    print("pre-registered rule: PASS iff |RHS/LHS - 1| < 1e-9 at every X")
    print("prediction on record: FAIL, ratio near 2, slowly drifting")
    print()

    Xs = [200, 400, 800, 1600, 3200]
    Nmax = max(Xs) + 2
    mu = mobius(Nmax)
    lam = vonmangoldt(Nmax)

    print(f"{'X':>6} {'LHS':>18} {'RHS(<X)':>18} {'RHS(<=X)':>18} "
          f"{'RHS/LHS(<X)':>13} {'RHS/LHS(<=X)':>13}")
    ratios = []
    for X in Xs:
        L, _ = lhs_exact(X, lam, mu)
        Rlt = rhs_box(X, lam, mu, "lt")
        Rle = rhs_box(X, lam, mu, "le")
        ratios.append((X, L, Rlt, Rle, Rlt / L, Rle / L))
        print(f"{X:>6} {L:>18.6e} {Rlt:>18.6e} {Rle:>18.6e} "
              f"{Rlt / L:>13.6f} {Rle / L:>13.6f}")

    print()
    print("cross-check of the derivation: the simplex-coupled form")
    print("(what the expansion actually gives) against the LHS")
    for X in [200, 400, 800]:
        L, _ = lhs_exact(X, lam, mu)
        S = rhs_simplex(X, lam, mu)
        print(f"  X={X:>5}  LHS={L:.10e}  simplex={S:.10e}  "
              f"ratio={S / L:.12f}")

    print()
    print("the repair: truncated convolution, summed over ALL N (to 2X)")
    print("            against the same box-range RHS")
    for X in [200, 400, 800, 1600]:
        R = rhs_box(X, lam, mu, "le")
        L2 = repair(X, lam, mu)
        print(f"  X={X:>5}  sum_N Chat(N)^2={L2:.10e}  RHS={R:.10e}  "
              f"ratio={L2 / R:.12f}")

    print()
    ok = all(abs(r[4] - 1.0) < 1e-9 for r in ratios)
    okle = all(abs(r[5] - 1.0) < 1e-9 for r in ratios)
    print(f"VERDICT (paper's literal range, args < X) : "
          f"{'PASS' if ok else 'FAIL'}")
    print(f"VERDICT (alternative range, args <= X)    : "
          f"{'PASS' if okle else 'FAIL'}")
    if not (ok or okle):
        rs = [r[4] for r in ratios]
        print(f"  ratio RHS/LHS over X={Xs}: "
              + ", ".join(f"{v:.4f}" for v in rs))
        print(f"  drift across the range: {rs[-1] / rs[0]:.4f}")
        print("  => Lemma MP is NOT an identity under either reading.")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
