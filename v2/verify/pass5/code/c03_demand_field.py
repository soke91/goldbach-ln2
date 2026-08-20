# -*- coding: utf-8 -*-
"""C-03.  Measurement {#meas:relocate} on an N-field that is not all 2^a 5^b.

P1's Measurement "what the weakening relocates" reports

    B(N)/N = 0.8086 ... 0.5916  against a threshold of 0.3745 N,
    so B(N)/(S(N)(1-A(N))N) falls 2.159 ... 1.580,

and concludes the object Proposition {#prop:nolog} constrains "is already
within a factor 1.6 of the threshold".  Every N in that field --
2e5, 4e5, 8e5, 1.6e6, 3.2e6 -- is 2^a * 5^b, so S(N)(1-A(N)) is the SAME
number 0.3745 at all five points: the sweep varies size only, never the
arithmetic that sets the threshold.  P1's own Measurement {#meas:margin}
puts the threshold's minimum at 0.060890 (primorials), a factor 6 below.

This script recomputes B(N)/N and B(N)/threshold on the same range with
primorial-type N included.
"""
import io
import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from indep import Lambda_bf, mu_bf, phi_bf  # noqa: E402

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                   "..", "results", "c03_demand_field.txt")
THETA = 0.56
PLIM = 4_000_000
CLIM = 4_000_000
# the target's field, then primorial-type N inside the same range
NS_TARGET = [200_000, 400_000, 800_000, 1_600_000, 3_200_000]
NS_EXTRA = [210_210, 510_510, 570_570, 1_021_020, 2_042_040, 2_072_070,
            223_092, 3_063_060]
lines = []

HEAD = [
    "STATISTIC: B(N) = sum_{k<K,(k,N)=1} (log k) |E_mu(N;k)| of P1 [eq:Bsum],",
    "           divided by N and by the one-sided threshold",
    "           S(N)(1-A(N))N of {#prop:onesided}/{#prop:nolog}; alongside",
    "           |E_3(N)|/N on the same k-range, so that |E_3|/B -- what the",
    "           signs mu(k) buy -- can be read off. Reported for the target's",
    "           own five N and for primorial-type N inside the same range.",
    "FIELD: theta' = THETA = 0.56, so K = floor(N^0.56) and k runs over the",
    "       squarefree k in [2,K) coprime to N (k=1 contributes log 1 = 0 to",
    "       both B and E_3, so the index set is not at issue here -- see",
    "       c02_flatsum_k1 for the sum where it is); N = the target's",
    "       2e5, 4e5, 8e5, 1.6e6, 3.2e6, every one of them 2^a 5^b, plus",
    "       eight N of other prime support in the same range; Lambda and mu",
    "       from an independent integer sieve to PLIM (mu by an omega-counter",
    "       plus a squarefull mask); phi by trial-division factorisation;",
    "       the Euler products for 2C_2 and Artin truncated at CLIM.",
    "CONSTANTS: THETA = 0.56; PLIM = 4000000 (sieve bound);",
    "           CLIM = 4000000 (Euler-product truncation for 2C_2 and",
    "           Artin -- audit_constants.py puts the tail there at ~2e-8,",
    "           four orders below the sixth printed decimal of S(1-A));",
    "           NS_TARGET = 2e5, 4e5, 8e5, 1.6e6, 3.2e6;",
    "           NS_EXTRA = 210210, 510510, 570570, 1021020, 2042040,",
    "           2072070, 223092, 3063060; brute-force cross-check bound 10^4.",
    "NULL: none applies. B(N) and its threshold are deterministic arithmetic",
    "      and the comparison between two N-fields is a comparison of two",
    "      deterministic quantities; there is no sign pattern to permute and",
    "      no random number is drawn. There is no seed.",
    "DENOM: B/N and thresh/N are divided by N; B/thresh is divided by",
    "       S(N)(1-A(N))N, i.e. by thresh itself.",
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


say("C-03  B(N) against its threshold, with the arithmetic of N varied")
say("=" * 86)
pr, lam, mu = my_sieves(PLIM)
bad = [n for n in range(10_001)
       if int(mu[n]) != mu_bf(n) or abs(float(lam[n]) - Lambda_bf(n)) > 1e-12]
say("  self-check of this script's sieve against trial division: %d "
    "mismatches" % len(bad))

artin, twin = 1.0, 2.0
for p in pr:
    p = int(p)
    if p > CLIM:
        break
    artin *= 1.0 - 1.0 / (p * (p - 1.0))
    if p > 2:
        twin *= 1.0 - 1.0 / (p - 1.0) ** 2
say("  Artin %.9f   2C_2 %.9f  (Euler products over p < %d)"
    % (artin, twin, CLIM))
say()
say("  N          P(N)                  S(N)     1-A(N)  thresh/N  B/N     "
    "  B/thresh  |E_3|/N")
say("  " + "-" * 86)

rows = []
for N in NS_TARGET + NS_EXTRA:
    v, PN, d = N, [], 2
    while d * d <= v:
        if v % d == 0:
            PN.append(d)
            while v % d == 0:
                v //= d
        d += 1
    if v > 1:
        PN.append(v)
    A, S = artin, twin
    for q in PN:
        A /= (1.0 - 1.0 / (q * (q - 1.0)))
        if q > 2:
            S *= (1.0 + 1.0 / (q - 2.0))
    K = int(N ** THETA)
    f = np.zeros(N, dtype=np.float64)
    idx = np.arange(1, N, dtype=np.int64)
    f[1:] = lam[1:N] * mu[N - idx]
    C = float(f.sum())
    B = 0.0
    E3 = 0.0
    for k in range(2, K):
        if mu[k] == 0 or any(k % q == 0 for q in PN):
            continue
        ph = phi_bf(k)
        r = N % k
        inner = float(f[r::k].sum()) if r else float(f[k::k].sum())
        e = inner - C / ph
        B += math.log(k) * abs(e)
        E3 += int(mu[k]) * math.log(k) * e
    thr = S * (1.0 - A) * N
    rows.append((N, PN, S, 1.0 - A, thr / N, B / N, B / thr, abs(E3) / N))
    say("  %-10d %-21s %-8.4f %-7.4f %-9.4f %-8.4f %-9.2f %.4f"
        % (N, ",".join(str(q) for q in PN), S, 1.0 - A, thr / N,
           B / N, B / thr, abs(E3) / N))

say()
say("  the target's five N all carry P(N) = {2,5}, so thresh/N is the")
say("  same 0.3745 at all five; B/thresh reproduces 2.159 ... 1.580.")
tgt = [r for r in rows if r[0] in NS_TARGET]
ext = [r for r in rows if r[0] in NS_EXTRA]
say("  B/thresh on the target field : %s"
    % ", ".join("%.3f" % r[6] for r in tgt))
say("  B/thresh on primorial-type N : %s"
    % ", ".join("%.2f" % r[6] for r in ext))
say("  worst on the target field %.2f;  worst here %.2f  (factor %.1f)"
    % (max(r[6] for r in tgt), max(r[6] for r in ext),
       max(r[6] for r in ext) / max(r[6] for r in tgt)))
say()
say("  B/N itself barely moves with the arithmetic of N -- it is the")
say("  THRESHOLD that collapses, so the ratio the Measurement quotes is a")
say("  property of the field it was measured on:")
say("    B/N   target %.4f-%.4f   primorial-type %.4f-%.4f"
    % (min(r[5] for r in tgt), max(r[5] for r in tgt),
       min(r[5] for r in ext), max(r[5] for r in ext)))
say("    thr/N target %.4f          primorial-type %.4f-%.4f"
    % (tgt[0][4], min(r[4] for r in ext), max(r[4] for r in ext)))

io.open(OUT, "w", encoding="utf-8", newline="\n").write(
    "\n".join(HEAD + lines) + "\n")
print("\nwrote", os.path.abspath(OUT))
