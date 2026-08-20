# -*- coding: utf-8 -*-
"""C-05.  Which object does P2's exponent -0.2658 belong to?

c02_flatsum_k1.txt prints "fitted log-weighted exponent : -0.2696" and P2
{#meas:flatsum} prints -0.2658 for what it calls |T_log|/N.  The deploy
session asked whether that is a fitting-convention difference or a defect.

It is neither: the two exponents belong to two DIFFERENT quantities that
lab_weight_gap.txt prints in two different places.

    sum_k (log k) H(N;k)      -- the ratio's denominator, Z4 table j=1.00
    T_log = E_3 = sum_k mu(k)(log k)[A(N;k) - C(N)/phi(k)]
                              -- the |E_3|/N column

They differ by C(N) * B_log(K), which does NOT vanish (B_log -> -S(N)),
unlike the w=1 branch where B_1 -> 0.  This script fits both, from the
target's own printed columns and from an independent recomputation, and
identifies which is which.

It also settles the second question: whether |T_1|/N = 0.01425 is the
k>=1 flat sum.  It is not -- it is T_1, which is invariant under the k=1
question, and it merely sits close to the k>=1 flat sum because B_1(K)
is small.
"""
import io
import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from indep import Lambda_bf, mu_bf, phi_bf  # noqa: E402

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                   "..", "results", "c05_exponent_convention.txt")
NS = [200_000, 400_000, 800_000, 1_600_000, 3_200_000, 6_400_000,
      12_800_000, 25_600_000]
THETA = 0.56
# the two columns lab_weight_gap.txt prints, transcribed from that file
TGT_E3 = [0.4377, 0.3837, 0.3172, 0.2608, 0.2073, 0.1855, 0.1459, 0.1245]
TGT_LOGH = [0.43948, 0.39668, 0.32530, 0.25645, 0.21066, 0.18167,
            0.14628, 0.12565]
PUB_EXP = -0.2658
lines = []

HEAD = [
    "STATISTIC: the fitted decay exponent of two quantities that P2",
    "           {#meas:flatsum} and lab_weight_gap.txt both print --",
    "           |sum_k (log k) H(N;k)|/N and |T_log|/N = |E_3|/N -- fitted",
    "           the same way (polyfit of log y on log N), from the target's",
    "           own printed columns and from an independent recomputation;",
    "           their difference C(N) B_log(K)/N; and the quantities that",
    "           the sentence 'at the top |T_1|/N = 0.01425 against",
    "           |T_log|/N = 0.1245' places next to the ratio list.",
    "FIELD: N = 2e5 through 2.56e7 by doubling; theta' = THETA = 0.56, so k",
    "       runs over the squarefree k < floor(N^0.56) coprime to N, both",
    "       from k = 1 and from k = 2; Lambda and mu from an independent",
    "       integer sieve to 2.56e7 (mu by an omega-counter plus a",
    "       squarefull mask); phi by trial-division factorisation; sums in",
    "       float64; exponents by np.polyfit of log y on log N, 8 points.",
    "CONSTANTS: THETA = 0.56; NS = 2e5 ... 2.56e7 by doubling; sieve bound",
    "           NMAX = 25600000; PUB_EXP = -0.2658 and the two transcribed",
    "           columns TGT_E3, TGT_LOGH are quoted from",
    "           the audited results/lab_weight_gap.txt, not recomputed from it.",
    "NULL: none applies. Every quantity is a deterministic sum over one",
    "      integer sieve and the question is which of two labelled objects a",
    "      published exponent belongs to. No random number is drawn; there",
    "      is no seed.",
    "DENOM: both fitted quantities are divided by N before the log-log fit.",
    "",
]


def say(s=""):
    print(s)
    lines.append(s)


def my_sieves(n):
    s = np.ones(n + 1, dtype=bool)
    s[:2] = False
    for p in range(2, int(n ** 0.5) + 1):
        if s[p]:
            s[p * p::p] = False
    pr = np.flatnonzero(s).astype(np.int64)
    omega = np.zeros(n + 1, dtype=np.int8)
    sqf = np.ones(n + 1, dtype=bool)
    for p in pr:
        p = int(p)
        omega[p::p] += 1
        if p * p <= n:
            sqf[p * p::p * p] = False
    mu = np.where(sqf, np.where(omega & 1, -1, 1), 0).astype(np.int8)
    mu[0] = 0
    lam = np.zeros(n + 1, dtype=np.float64)
    lg = np.log(pr.astype(np.float64))
    lam[pr] = lg
    for i, p in enumerate(pr):
        p = int(p)
        if p * p > n:
            break
        q = p * p
        while q <= n:
            lam[q] = lg[i]
            if q > n // p:
                break
            q *= p
    return pr, lam, mu


def fit(ys):
    x = np.log(np.array(NS, dtype=float))
    return float(np.polyfit(x, np.log(np.array(ys, dtype=float)), 1)[0])


say("C-05  which object owns the exponent -0.2658")
say("=" * 78)
say()
say("PART 1.  fitted from lab_weight_gap.txt's own printed columns")
say("-" * 78)
say("  |E_3|/N column          -> exponent %.4f   (P2 prints %.4f)"
    % (fit(TGT_E3), PUB_EXP))
say("  |sum (log k) H|/N (j=1) -> exponent %.4f" % fit(TGT_LOGH))
say("  c02_flatsum_k1.txt printed -0.2696 for the SECOND of these.")
say()

pr, lam, mu = my_sieves(max(NS))
bad = [n for n in range(10_001)
       if int(mu[n]) != mu_bf(n) or abs(float(lam[n]) - Lambda_bf(n)) > 1e-12]
say("PART 2.  independent recomputation (sieve self-check: %d mismatches)"
    % len(bad))
say("-" * 78)
say("  N            |E_3|/N   |S(log k)H|/N  C*B_log/N   |T_1|/N   "
    "flat k>=1  flat k>=2")
rows = []
for N in NS:
    v, PN, d = N, set(), 2
    while d * d <= v:
        if v % d == 0:
            PN.add(d)
            while v % d == 0:
                v //= d
        d += 1
    if v > 1:
        PN.add(v)
    K = int(N ** THETA)
    ks = np.array([k for k in range(1, K)
                   if mu[k] != 0 and all(k % q for q in PN)])
    f0 = np.zeros(N, dtype=np.float64)
    idx = np.arange(1, N, dtype=np.int64)
    f0[1:] = lam[1:N] * mu[N - idx]
    C = float(f0.sum())
    A = np.empty(ks.size)
    for i, k in enumerate(ks):
        k = int(k)
        r = N % k
        A[i] = f0[r::k].sum() if r else f0[k::k].sum()
    sgn = mu[ks].astype(np.float64)
    H = sgn * A
    lg = np.log(ks.astype(float))
    iph = np.array([phi_bf(int(k)) for k in ks], dtype=np.float64)

    logH = float((lg * H).sum())
    Blog = float((sgn * lg / iph).sum())
    E3 = float((sgn * lg * (A - C / iph)).sum())
    B1 = float((sgn / iph).sum())
    T1 = float((sgn * (A - C / iph)).sum())
    flat1 = float(H.sum())
    flat2 = float(H[ks > 1].sum())
    rows.append((N, abs(E3) / N, abs(logH) / N, C * Blog / N, abs(T1) / N,
                 abs(flat1) / N, abs(flat2) / N, C / N, B1, Blog))
    say("  %-12d %-9.5f %-14.5f %+-11.5f %-9.5f %-10.6f %.6f"
        % rows[-1][:7])

say()
say("  exponent of |E_3|/N            : %.4f" % fit([r[1] for r in rows]))
say("  exponent of |sum (log k) H|/N  : %.4f" % fit([r[2] for r in rows]))
say()
say("  the two differ by C(N) B_log(K), and B_log -> -S(N) does NOT vanish:")
for r in rows:
    say("    N = %-10d  B_log(K) = %+.5f   C(N)/N = %+.6f   "
        "C*B_log/N = %+.6f" % (r[0], r[9], r[7], r[3]))

say()
say("PART 3.  what |T_1|/N = 0.01425 is, and what it is not")
say("-" * 78)
say("  T_1 = sum_k mu(k)[A(N;k) - C(N)/phi(k)] and its k=1 term is")
say("  A(N;1) - C(N)/phi(1) = C(N) - C(N) = 0 exactly, so T_1 does not")
say("  depend on the index-set question at all.  It is CLOSE to the k>=1")
say("  flat sum only because T_1 = sum_{k>=1} H - C(N) B_1(K) and B_1 is")
say("  small -- that is Huang-Li's Lemma 1, and it is the whole reason the")
say("  w=1 branch is the proved one.")
say()
say("  N            |T_1|/N    flat k>=1   gap = |C*B_1|/N   B_1(K)")
for r in rows:
    say("  %-12d %-10.6f %-11.6f %-17.2e %+.6f"
        % (r[0], r[4], r[5], abs(r[7] * r[8]), r[8]))
say()
say("  So at N = 2.56e7: |T_1|/N = %.6f and the k>=1 flat sum is %.6f."
    % (rows[-1][4], rows[-1][5]))
say("  They agree to four decimals (both round to 0.01425) but they are")
say("  not the same number, and the k>=2 flat sum %.6f is a third one."
    % rows[-1][6])

io.open(OUT, "w", encoding="utf-8", newline="\n").write(
    "\n".join(HEAD + lines) + "\n")
print("\nwrote", os.path.abspath(OUT))
