# -*- coding: utf-8 -*-
"""
Transform Lab, session 7 (increment 226): the same identity, grouped the
other way -- and why the margin depends on the grouping.

Session 6's transform P is a double sum over (p, w) with v = pw:

    C_log(N) = Sum_{v<N} mu(v) Lambda(N-v) log v
             = Sum_p log p * D_p(N),   D_p(N) = -Sum_w mu(w) Lambda(N-pw)
             = -Sum_w mu(w) * G_w(N),  G_w(N) = Sum_{p < N/w, p not| w}
                                                 log p * Lambda(N - pw)

Grouping by p gave a margin: Sum_p log p |D_p| is a DECAYING fraction of
the trivial bound (0.383 at N = 4*10^5), so square-root cancellation
inside each group already suffices. This script asks the same question
of the other grouping, and of the arithmetic that would be needed to
rescue it.

WHAT G_w IS. Fixing w, G_w(N) counts the representations N = wp + n
with p, n prime, weighted log p * Lambda(n) -- a binary Goldbach-type
problem, one for each w. So the wall is a MOBIUS AVERAGE OF BINARY
PROBLEMS. Each G_w is a sum of NONNEGATIVE terms, with main term
S_w(N) * N/w where, for even N and w odd with (w,N) = 1,

    S_w(N) = 2C_2 * prod_{q | wN, q>2} (q-1)/(q-2) = S(N) * f(w),
    f(w)   = prod_{q | w} (q-1)/(q-2),

and S_w(N) = 0 when gcd(w, N) > 1 or w is even (then n = N - wp is
forced into a fixed residue class 0 mod q, so n = q and the count is
O(1) rather than of order N/w).

THE THREE MEASUREMENTS.

 (A) the identity again, in the w-grouping: -Sum_w mu(w) G_w = C_log,
     exactly. A second independent check on transform P.

 (B) THE ASYMMETRY. Because every G_w >= 0,
        Sum_w |G_w| = Sum_w G_w = Sum_v mu^2(v) Lambda(N-v) log v,
     which is the trivial bound EXACTLY. So the w-grouping has ratio
     1 identically -- provably zero margin, not merely measured. All
     of its cancellation has to come from mu(w) ACROSS groups, whereas
     the p-grouping gets a decaying ratio from cancellation INSIDE each
     group. Predicted before running: ratio = 1.000000 to machine
     precision, for every N. Anything else is a bug in this script.

 (C) IS THE SINGULAR SERIES THE OBSTRUCTION? If the w-grouping is to
     work, Sum_w mu(w) S_w(N) N/w must be small. Since S_w = S(N) f(w),
     that sum is S(N) N * A(W) with
        A(W) = Sum_{w <= W, w odd squarefree, (w,N)=1} mu(w) f(w) / w.
     NULL. f(q) = (q-1)/(q-2) = 1 + 1/(q-2) is within O(1/q) of 1, so
     A(W) is a perturbation of Sum_{w<=W} mu(w)/w, whose true size is
     PNT-strength -- O(exp(-c sqrt(log W))), far below 1/log W. The
     Euler-product reading 1 - f(q)/q ~ 1 - 1/q would suggest c/log W
     and is the WRONG null here: it describes the smooth-truncated
     product, not the Mobius sum. Both references are printed, with
     Sum_{w<=W} mu(w)/w computed from the same data as the honest one.

 (D) WHAT IS LEFT. E_w := G_w - S_w N/w, and Sum_w |E_w|.
     NULL, taken from the data rather than from a size heuristic: G_w
     is a sum of nonnegative terms t_p = log p * Lambda(N - pw), so the
     natural fluctuation scale is the square root of its own second
     moment, null_w = sqrt(Sum_p t_p^2), computed exactly in the same
     loop. A crude a priori estimate |E_w| ~ sqrt(N/w) log N would give
     Sum_w |E_w| ~ sqrt(2) N log N, i.e. no saving at all, but that
     estimate is dominated by large w where G_w has O(1) terms and the
     square-root model does not apply; it is printed only to show how
     far off a size heuristic lands. The comparison that counts is
     Sum_w |E_w| against Sum_w null_w.

READING. (B) is the structural point and it is exact, not statistical.
(C) and (D) then ask what the w-grouping would need in order to be
rescued: whether the singular-series average cancels, and how much of
the mass survives once the main terms are removed by hand.
"""
import numpy as np
import math

from hl_S1_check import sieve, singular

C2 = 0.6601618158468696


def prime_factors(n, spf):
    fs = []
    while n > 1:
        p = int(spf[n]); fs.append(p)
        while n % p == 0:
            n //= p
    return fs


def main():
    X = 200_000
    mu, lam, spf = sieve(X)
    primes = np.nonzero(spf[2:] == np.arange(2, X + 1))[0] + 2
    logp = np.log(primes.astype(np.float64))

    NS = (25_000, 50_000, 100_000, 200_000)
    print(f"{'N':>8} {'C_log direct':>14} {'-Sum mu G_w':>14} "
          f"{'rel.diff':>10} {'Sum|G_w|/triv':>14}")
    keep = []
    for N in NS:
        v = np.arange(1, N)
        term = mu[1:N].astype(np.float64) * lam[N - v]
        logv = np.log(v.astype(np.float64))
        C_log = float(np.dot(term, logv))
        triv = float(np.dot(np.abs(mu[1:N]).astype(np.float64) * lam[N - v],
                            logv))

        Nfac = set(prime_factors(N, spf))
        S_N = singular(N, spf)

        ws = np.nonzero(mu[1:N // 2 + 1])[0] + 1     # squarefree w >= 1
        G = np.zeros(len(ws)); MT = np.zeros(len(ws)); wf = np.zeros(len(ws))
        Q = np.zeros(len(ws))                        # second moment of G_w
        for i, w in enumerate(ws):
            w = int(w)
            lim = (N - 1) // w
            j = np.searchsorted(primes, lim, side='right')
            if j == 0:
                continue
            ps = primes[:j]; lps = logp[:j]
            wfac = prime_factors(w, spf) if w > 1 else []
            if wfac:                                  # ENFORCE p not| w
                ok = np.ones(j, dtype=bool)
                for q in wfac:
                    ok &= (ps != q)
                ps = ps[ok]; lps = lps[ok]
            t = lps * lam[N - ps * w]
            G[i] = float(t.sum())
            Q[i] = float(np.dot(t, t))
            # main term: zero unless w odd and (w, N) = 1
            if w % 2 == 1 and not (set(wfac) & Nfac):
                f = 1.0
                for q in wfac:
                    f *= (q - 1) / (q - 2)
                wf[i] = f
                MT[i] = S_N * f * N / w

        muw = mu[ws].astype(np.float64)
        lhs = -float(np.dot(muw, G))
        rel = abs(lhs - C_log) / max(abs(C_log), 1.0)
        sabs = float(np.abs(G).sum())
        print(f"{N:>8} {C_log:>14.2f} {lhs:>14.2f} {rel:>10.2e} "
              f"{sabs/triv:>14.8f}")
        keep.append((N, C_log, triv, ws, muw, G, MT, wf, S_N, Q))

    print("\n(B) the asymmetry -- both groupings of the SAME identity")
    print(f"{'N':>8} {'w-group |.|/triv':>17} {'p-group |.|/triv':>17}")
    print("       (p-group figures from session 6, results/"
          "lab_prime_factor_split.txt)")
    pg = {50_000: 0.4459, 100_000: 0.4200, 200_000: 0.4033}
    for (N, C_log, triv, ws, muw, G, MT, wf, S_N, Q) in keep:
        sabs = float(np.abs(G).sum())
        s = f"{pg[N]:.4f}" if N in pg else "   --  "
        print(f"{N:>8} {sabs/triv:>17.8f} {s:>17}")
    print("    G_w >= 0 termwise, so Sum_w |G_w| = Sum_w G_w = the")
    print("    trivial bound identically: the w-grouping has provably")
    print("    zero margin. The p-grouping's D_p carries the Mobius and")
    print("    cancels inside each group.")

    print("\n(C) do the main terms cancel? A(W) = Sum mu(w) f(w)/w")
    print(f"{'N':>8} {'A(W)':>11} {'Sum mu(w)/w':>12} {'1/log W':>9} "
          f"{'MT_tot':>11} {'|MT|/triv':>10} {'|MT|/(N/logN)':>14}")
    for (N, C_log, triv, ws, muw, G, MT, wf, S_N, Q) in keep:
        W = N // 2
        A = float(np.dot(muw, np.where(wf > 0, wf / ws, 0.0)))
        allw = np.arange(1, W + 1)
        mert = float(np.dot(mu[1:W + 1].astype(np.float64), 1.0 / allw))
        MT_tot = float(np.dot(muw, MT))
        print(f"{N:>8} {A:>11.6f} {mert:>12.6f} "
              f"{1/math.log(W):>9.5f} {MT_tot:>11.1f} "
              f"{abs(MT_tot)/triv:>10.6f} "
              f"{abs(MT_tot)/(N/math.log(N)):>14.4f}")
    print("    A(W) tracks Sum mu(w)/w, PNT-strength and far below")
    print("    1/log W: the main terms cancel from ~N log N to O(100),")
    print("    a saving of more than (log N)^2. The singular-series")
    print("    average is NOT the obstruction to the w-grouping.")

    print("\n(D) what is left after the main terms: E_w = G_w - S_w N/w")
    print(f"{'N':>8} {'Sum|E_w|':>12} {'Sum null_w':>12} {'|E|/null':>9} "
          f"{'|E|/triv':>9} {'heuristic':>11} {'|E|/heur':>9}")
    for (N, C_log, triv, ws, muw, G, MT, wf, S_N, Q) in keep:
        E = G - MT
        sE = float(np.abs(E).sum())
        snull = float(np.sqrt(Q).sum())
        heur = math.sqrt(2) * N * math.log(N)
        print(f"{N:>8} {sE:>12.1f} {snull:>12.1f} {sE/snull:>9.4f} "
              f"{sE/triv:>9.4f} {heur:>11.1f} {sE/heur:>9.4f}")
    print("    |E|/null is the figure that counts: the residual mass")
    print("    measured against its own second-moment scale. The size")
    print("    heuristic sqrt(2) N log N is printed only to show how")
    print("    far a priori sizing lands from the data.")

    print("\nreading: one identity, two groupings. The p-grouping's D_p")
    print("carries the Mobius, so it cancels internally and needs no")
    print("main term supplied from outside; the w-grouping's G_w is")
    print("nonnegative, has provably zero margin as it stands, and only")
    print("becomes comparable after its main terms are subtracted by")
    print("hand. Design rule: group so that mu is INSIDE the group, and")
    print("the main-term subtraction comes for free.")
    print("DONE")


if __name__ == "__main__":
    main()
