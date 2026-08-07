# -*- coding: utf-8 -*-
"""
Is the mask the whole of the atomic structure at its own moduli?
(increment 326)

WHY. Two objects in this program were derived independently and never
compared:

  the LOCATION MASK, from M.1 -- a theorem about divisibility: n prime
    and q | N force q not| (N-n). Its estimator is the mean of
    Z = C/sqrt(V) over the 32 cells of the divisibility pattern of
    {3,5,7,11,13}.
  the ATOMIC SPECTRUM, from the covariance -- increments 324 and 325:
    Z's spectral measure is atomic on the rationals with
    Hardy-Littlewood weights.

Both live at the same frequencies. The cell means are functions of
N mod 30030, so they sit inside the span of the atoms with q | 30030.
The question nobody has asked is whether they are the SAME thing:

    is the atomic structure at q | 30030 exhausted by the 32-dimensional
    divisibility-cell subspace, or is there more there?

If exhausted, the mask IS the atomic component and one object has been
found twice. If not, there is periodic structure at those moduli that
M.1 does not reach -- and the coin control says whether it belongs to
mu or to Lambda.

PRE-REGISTRATION (fixed before the run).

  (S1) CONTAINMENT, a structural self-test. The cell-mean projection
       P_cell Z is a function of N mod 30030, so almost all of its
       energy must lie in the atomic bins with q | 30030. RULE: at
       least 95% does. If not, the two spans are not nested and
       nothing below compares like with like.

  (S2) IS THERE MORE? Report ||P_atom Z||^2 / ||P_cell Z||^2 for the
       real Z. A ratio near 1 means the mask exhausts the atomic
       structure at its own moduli; a ratio well above 1 means it does
       not.

  (S3) WHOSE IS THE EXCESS? The same ratio for a coin, which has no
       mask but has the same covariance. RULE: if the real's ratio and
       the coin's agree within 20%, the excess is Lambda's and the
       mask is simply the part of it that M.1 explains. If the real's
       is materially larger, there is atomic structure that needs mu --
       which would be the first such finding in this program.

  WHAT WOULD REFUTE. (S1) failing invalidates the comparison. (S3)
  coming out with a real excess would be a genuine finding and would
  need its own replication before anything is claimed.
"""
import math
import sys
import time

import numpy as np

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

MOD = 30030
QP = [3, 5, 7, 11, 13]


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


def divisors(ps):
    out = [1]
    for p in ps:
        out += [d * p for d in out]
    return sorted(out)


def main():
    X = 8_000_000
    lo = 200_000
    t0 = time.time()
    mu, lam = sieve(X)
    nf = 1
    while nf < 2 * (X + 1):
        nf *= 2
    F_lam = np.fft.rfft(np.pad(lam, (0, nf - X - 1)))
    V = np.fft.irfft(np.fft.rfft(np.pad((mu != 0).astype(np.float64),
                                        (0, nf - X - 1)))
                     * np.fft.rfft(np.pad(lam ** 2, (0, nf - X - 1))),
                     nf)[: X + 1]
    C = np.fft.irfft(np.fft.rfft(np.pad(mu.astype(np.float64),
                                        (0, nf - X - 1))) * F_lam,
                     nf)[: X + 1]
    Ns = np.arange(lo, X + 1, 2)
    n = (len(Ns) // MOD) * MOD
    Ns = Ns[:n]
    key = np.zeros(n, dtype=np.int32)
    for i, q in enumerate(QP):
        key |= ((Ns % q) == 0).astype(np.int32) << i
    uq, inv = np.unique(key, return_inverse=True)
    cnt = np.bincount(inv).astype(np.float64)
    print(f"n = {n}, {len(uq)} cells   t={time.time()-t0:.0f}s",
          flush=True)

    atom = np.zeros(n // 2 + 1, dtype=bool)
    for q in divisors(QP):
        if q == 1:
            continue
        for j in range(1, q):
            if math.gcd(j, q) == 1:
                b = (j * n) // q
                if 0 < b < len(atom):
                    atom[b] = True
    print(f"    {int(atom.sum())} atomic bins at q | {MOD}", flush=True)

    def analyse(Z, tag):
        Z = Z - Z.mean()
        tot = float((Z * Z).sum())
        cm = np.bincount(inv, weights=Z) / cnt
        Pcell = cm[inv]
        F = np.fft.rfft(Z)
        Fa = np.where(atom, F, 0.0)
        Patom = np.fft.irfft(Fa, n)
        # containment: how much of P_cell lies in the atomic bins
        Fc = np.fft.rfft(Pcell)
        num = float((np.abs(np.where(atom, Fc, 0.0)) ** 2).sum())
        den = float((np.abs(Fc[1:]) ** 2).sum())
        contain = num / den if den > 0 else float("nan")
        ec = float((Pcell * Pcell).sum()) / tot
        ea = float((Patom * Patom).sum()) / tot
        print(f"    {tag:<6} cell {ec:>9.5f}   atom {ea:>9.5f}   "
              f"atom/cell {ea/ec:>8.3f}   containment {contain:>7.4f}")
        return ec, ea, contain

    print(f"\n(S1)(S2)(S3) energy shares of Z, and the two projections")
    print(f"    {'':<6} {'cell':>14} {'atom':>14} {'ratio':>18} "
          f"{'containment':>19}")
    ecr, ear, conr = analyse(C[Ns] / np.sqrt(V[Ns]), "real")
    rng = np.random.default_rng(326)
    idx = np.nonzero(mu != 0)[0]
    eps = np.zeros(nf)
    eps[idx] = rng.integers(0, 2, size=len(idx)) * 2.0 - 1.0
    Cc = np.fft.irfft(np.fft.rfft(eps) * F_lam, nf)[: X + 1]
    ecc, eac, conc = analyse(Cc[Ns] / np.sqrt(V[Ns]), "coin")

    okS1 = conr > 0.95 and conc > 0.95
    rr, rc = ear / ecr, eac / ecc
    okS3 = abs(rr / rc - 1.0) <= 0.20
    print(f"\n    (S1) the cell projection lies inside the atomic span "
          f"(>95%): {'PASS' if okS1 else 'FAIL'}  "
          f"(real {conr:.4f}, coin {conc:.4f})")
    print(f"    (S2) atom/cell energy: real {rr:.3f}, coin {rc:.3f}")
    print(f"    (S3) the two ratios agree within 20%: "
          f"{'PASS' if okS3 else 'FAIL'}  "
          f"(real/coin = {rr/rc:.3f})")
    # (S4) THE MECHANISM, checked rather than guessed. The coin
    # carries MORE atomic energy than the real, which is the opposite
    # of what a mask would do. The candidate: at a rational j/q the
    # prime side has a major arc while mu's exponential sum is small
    # (Davenport), and a random sign sequence has no such cancellation
    # -- |sum eps(v) e(v f)|^2 sits at the support count for every f.
    # If that is it, the ratio of the two exponential sums at the
    # atomic frequencies should reproduce the energy gap.
    print("\n(S4) the two sign sequences at the atomic frequencies")
    sup = (mu != 0)
    vv = np.nonzero(sup)[0]
    mv = mu[sup].astype(np.float64)
    ev = eps[vv]
    Qn = float(len(vv))
    # e(v j/q) depends only on v mod q, so one bincount per modulus
    # replaces one length-Q exponential per frequency. The first draft
    # did the latter -- 255 frequencies x 4.9e6 terms -- and had to be
    # killed.
    rm, re_ = [], []
    for q in divisors(QP):
        if q == 1:
            continue
        r = (vv % q).astype(np.int64)
        am = np.bincount(r, weights=mv, minlength=q)
        ae = np.bincount(r, weights=ev, minlength=q)
        for jj in range(1, q):
            if math.gcd(jj, q) != 1:
                continue
            w = np.exp(2j * np.pi * jj * np.arange(q) / q)
            rm.append(abs(complex(np.dot(am, w))) ** 2 / Qn)
            re_.append(abs(complex(np.dot(ae, w))) ** 2 / Qn)
    rm, re_ = np.array(rm), np.array(re_)
    gap = float(re_.mean() / rm.mean())
    print(f"    mean |sum mu(v) e(v j/q)|^2 / Q  = {rm.mean():>12.4f}")
    print(f"    mean |sum eps(v) e(v j/q)|^2 / Q = {re_.mean():>12.4f}")
    print(f"    coin / real                      = {gap:>12.2f}x")
    print(f"    energy gap to explain (coin/real atom-to-cell) "
          f"= {(eac/ecc)/(ear/ecr):>7.2f}x")
    okS4 = gap > 3.0
    print(f"    (S4) mu's exponential sum is small at rationals where "
          f"a coin's is not: {'PASS' if okS4 else 'FAIL'}  "
          f"({gap:.1f}x)")

    if okS1 and okS3:
        v = (f"the mask does NOT exhaust the atomic structure at its own "
             f"moduli -- there is {rr:.1f}x more energy in the atoms "
             f"than the 32 divisibility cells capture -- but the coin "
             f"shows the same factor, {rc:.1f}x, so the excess is "
             f"Lambda's through the shift and not mu's. The mask is the "
             f"part of that structure M.1 explains, and the rest is the "
             f"same covariance seen at finer frequencies")
    elif okS1 and okS4:
        v = (f"the coin carries {rc/rr:.1f}x the real's atom-to-cell "
             f"ratio, and the mechanism is measured, not guessed: at "
             f"the atomic frequencies mu's exponential sum is "
             f"{gap:.0f}x smaller than a random sign sequence's. That "
             f"is Davenport's bound appearing in the spectrum. ⚠ It "
             f"also means THE COIN IS A POOR NULL FOR RATIONAL-"
             f"FREQUENCY STRUCTURE: a random sign sequence has no "
             f"major-arc cancellation while mu does, so the coin "
             f"OVERSTATES atomic energy several-fold")
    elif okS1:
        v = (f"the coin carries {rc/rr:.2f}x the real's atom-to-cell "
             f"ratio and the exponential-sum mechanism does not "
             f"account for it; the difference needs its own "
             f"explanation before anything is claimed")
    else:
        v = ("the cell projection does not lie inside the atomic span, "
             "so the two objects are not nested and the comparison "
             "does not read")
    print(f"\n    {v}")
    print("DONE")


if __name__ == "__main__":
    main()
