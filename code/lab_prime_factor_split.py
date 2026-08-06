# -*- coding: utf-8 -*-
"""
Transform Lab, session 6 (increment 225): the prime-factor split, and
whether it leaves any margin.

THE CANDIDATE. mu vanishes off the squarefree numbers, so on the
support of mu the identity log v = Sum_{p | v} log p is EXACT (no
higher prime powers to account for). Feeding it into the wall,

    C_log(N) := Sum_{v<N} mu(v) Lambda(N-v) log v
              = Sum_p log p * D_p(N),
    D_p(N)   := Sum_{v<N, p | v} mu(v) Lambda(N-v).

An exact identity with no error term. Equivalently, writing v = pw
(p not| w automatic on the support of mu), D_p(N) = -Sum_w mu(w)
Lambda(N - pw): the Mobius variable is dilated and the additive form
becomes N - pw, i.e. Lambda restricted to the progression N mod p.

WHY THIS CANDIDATE. The three-clause test from sessions 1-5 asks
whether a transform (1) moves a finite pencil vertex to infinity,
(2) removes a divisibility restriction we inserted ourselves, or
(3) requires the two factors to dilate together. This one does none:
it is not a divisor switch, the restriction p | v is intrinsic to v
rather than inserted, and only mu is dilated -- Lambda is redirected to
a progression instead of being asked to transform.

WHAT THIS SCRIPT MEASURES. Being lossless is necessary, not
sufficient; the question is whether the split leaves MARGIN. Write

    M_p(N) := Sum_{v<N, p | v} mu^2(v) Lambda(N-v)      (trivial mass)

so that |D_p| <= M_p termwise and Sum_p log p M_p = Sum_v mu^2(v)
Lambda(N-v) log v, the trivial bound for C_log. The wall needs
C_log(N) = o(N log N). Three things get measured, in order:

 (A) the identity, exactly, as a check on the transform itself;
 (B) losslessness: Sum_p log p M_p against the trivial bound;
 (C) THE DECISIVE ONE -- the absolute aggregate
        S_abs(N) = Sum_p log p |D_p(N)|,
     which uses NO cancellation across p at all. If S_abs = o(N log N)
     then per-p cancellation alone suffices and the split has margin;
     if S_abs is a fixed fraction of the trivial bound, it has none.

NULLS, TAKEN FROM THE DATA AND NOT FROM A SIZE HEURISTIC. D_p is a
signed sum of the terms t_v = mu(v) Lambda(N-v) over v divisible by p,
so under a random-sign model its scale is the square root of its own
second moment, which is computable in the same loop:

    V_p    := Sum_{v<N, p|v} mu^2(v) Lambda(N-v)^2,
    null_p := sqrt(V_p),          rho_p := |D_p| / M_p.

    S_null(N) := Sum_p log p * sqrt(V_p).

A size heuristic would instead write null_p = sqrt(M_p log N), on the
grounds that D_p has ~ M_p/log N terms of size ~ log N; that form is
also printed, as a check on how far such reasoning lands, but it is
NOT the null the verdict is read against (hazard 4, CLOSURE_REAUDIT).

Asymptotically V_p ~ c N log N / p, so S_null ~ sqrt(c N log N) *
Sum_{p<N} log p / sqrt(p) ~ 2 sqrt(c) N sqrt(log N) and the PREDICTION
under square-root-per-p is

    S_null / (N log N)  ~  const / sqrt(log N),

decaying, with sqrt(log N) to spare over what the wall needs. Both
S_abs and S_null are printed beside the trivial bound so the measured
ratio is read against its own null rather than a bare threshold.

LIMIT OF THE TEST, STATED UP FRONT. Over N in [5*10^4, 4*10^5] the
factor sqrt(log N) changes by only sqrt(12.9/10.8) = 1.093, so a
predicted 9% decay is all that is on offer. This measurement can
therefore distinguish "decaying like the null" from "flat or growing",
and nothing finer. The fitted exponent is reported with that caveat
attached, not as a confirmation of -1/2.
"""
import numpy as np
import math


def sieve(X):
    spf = np.zeros(X + 1, dtype=np.int32)
    for i in range(2, int(X ** 0.5) + 1):
        if spf[i] == 0:
            sl = spf[i * i::i]; sl[sl == 0] = i
    for i in range(2, X + 1):
        if spf[i] == 0:
            spf[i] = i
    mu = np.zeros(X + 1, dtype=np.int8); mu[1] = 1
    for i in range(2, X + 1):
        p = int(spf[i]); j = i // p
        mu[i] = 0 if j % p == 0 else -mu[j]
    primes = np.nonzero(spf[2:] == np.arange(2, X + 1))[0] + 2
    lam = np.zeros(X + 1)
    for p in primes:
        q = int(p); lp = math.log(int(p))
        while q <= X:
            lam[q] = lp; q *= int(p)
    return mu, lam, spf, primes


def split(N, mu, lam, primes):
    """D_p and M_p for every p < N, plus the two global sums."""
    v = np.arange(1, N)
    muv = mu[1:N].astype(np.float64)
    lamr = lam[N - v]
    term = muv * lamr
    absterm = np.abs(muv) * lamr
    logv = np.log(v.astype(np.float64))
    C = float(term.sum())
    C_log = float(np.dot(term, logv))
    M_log = float(np.dot(absterm, logv))

    ps = primes[primes < N]
    D = np.empty(len(ps)); M = np.empty(len(ps)); V = np.empty(len(ps))
    for i, p in enumerate(ps):
        idx = np.arange(int(p), N, int(p))
        m = mu[idx].astype(np.float64); L = lam[N - idx]
        t = m * L
        D[i] = t.sum()
        M[i] = np.dot(np.abs(m), L)
        V[i] = np.dot(t, t)                 # second moment, exact
    return C, C_log, M_log, ps, D, M, V


def main():
    X = 400_000
    mu, lam, spf, primes = sieve(X)
    NS = (50_000, 100_000, 200_000, 400_000)

    print("(A) the identity  C_log(N) = Sum_p log p D_p(N)")
    print(f"{'N':>8} {'C(N)':>12} {'C_log(N)':>14} {'Sum lp*D_p':>14} "
          f"{'rel.diff':>10}")
    rows = []
    for N in NS:
        C, C_log, M_log, ps, D, M, V = split(N, mu, lam, primes)
        lp = np.log(ps.astype(np.float64))
        lhs = float(np.dot(lp, D))
        rel = abs(lhs - C_log) / max(abs(C_log), 1.0)
        print(f"{N:>8} {C:>12.2f} {C_log:>14.2f} {lhs:>14.2f} "
              f"{rel:>10.2e}")
        rows.append((N, C, C_log, M_log, ps, D, M, V, lp))
    print("    exact identity: no error term, and mu kills the")
    print("    non-squarefree v where log v = Sum_{p|v} log p would fail")

    print("\n(B) losslessness -- Sum_p log p M_p vs the trivial bound")
    print(f"{'N':>8} {'Sum lp*M_p':>14} {'M_log':>14} {'rel.diff':>10} "
          f"{'/(N log N)':>11}")
    for (N, C, C_log, M_log, ps, D, M, V, lp) in rows:
        s = float(np.dot(lp, M))
        rel = abs(s - M_log) / M_log
        print(f"{N:>8} {s:>14.2f} {M_log:>14.2f} {rel:>10.2e} "
              f"{s/(N*math.log(N)):>11.4f}")
    print("    the split neither gains nor loses at the trivial scale:")
    print("    it redistributes the same mass over p")

    print("\n(C) the decisive one -- S_abs uses NO cancellation over p")
    print(f"{'N':>8} {'S_abs':>13} {'S_null':>13} {'triv':>13} "
          f"{'S_abs/triv':>11} {'S_null/triv':>12} {'S_abs/S_null':>13} "
          f"{'heur/S_null':>12}")
    meas = []
    for (N, C, C_log, M_log, ps, D, M, V, lp) in rows:
        logN = math.log(N)
        S_abs = float(np.dot(lp, np.abs(D)))
        S_null = float(np.dot(lp, np.sqrt(V)))
        S_heur = float(np.dot(lp, np.sqrt(M * logN)))
        triv = float(np.dot(lp, M))
        meas.append((N, S_abs, S_null, triv))
        print(f"{N:>8} {S_abs:>13.1f} {S_null:>13.1f} {triv:>13.1f} "
              f"{S_abs/triv:>11.4f} {S_null/triv:>12.4f} "
              f"{S_abs/S_null:>13.4f} {S_heur/S_null:>12.4f}")

    # fitted exponent in log N, reported with its caveat
    xs = np.array([math.log(math.log(N)) for (N, a, b, t) in meas])
    ya = np.array([math.log(a / t) for (N, a, b, t) in meas])
    yn = np.array([math.log(b / t) for (N, a, b, t) in meas])
    ea = float(np.polyfit(xs, ya, 1)[0])
    en = float(np.polyfit(xs, yn, 1)[0])
    print(f"\n    fitted: S_abs/triv ~ (log N)^{ea:+.3f}, "
          f"S_null/triv ~ (log N)^{en:+.3f}   (null exponent -0.5)")
    print("    caveat: over this N range sqrt(log N) moves by 9%, so")
    print("    these exponents separate 'decaying' from 'flat', no more")

    print("\n(D) where the mass sits -- dyadic profile at N = 400000")
    (N, C, C_log, M_log, ps, D, M, V, lp) = rows[-1]
    logN = math.log(N)
    triv = float(np.dot(lp, M))
    print(f"{'p range':>16} {'#p':>7} {'mass frac':>10} {'abs frac':>9} "
          f"{'mean rho':>9} {'null rho':>9} {'ratio':>7}")
    b = 2
    while b < N:
        sel = (ps >= b) & (ps < min(2 * b, N))
        if sel.sum():
            m = M[sel]; d = D[sel]; l = lp[sel]
            wm = float(np.dot(l, m))
            wa = float(np.dot(l, np.abs(d)))
            rho = float(np.dot(l * m, np.abs(d) / np.maximum(m, 1e-9))
                        ) / max(wm, 1e-9)
            v = V[sel]
            rn = float(np.dot(l * m, np.sqrt(v) / np.maximum(m, 1e-9))
                       ) / max(wm, 1e-9)
            print(f"{b:>7}-{min(2*b, N):>8} {int(sel.sum()):>7} "
                  f"{wm/triv:>10.4f} {wa/triv:>9.4f} {rho:>9.4f} "
                  f"{rn:>9.4f} {rho/max(rn,1e-9):>7.3f}")
        b *= 2
    print("    mass frac = share of the trivial bound carried by the")
    print("    range; rho = |D_p|/M_p averaged with weight log p * M_p")

    print("\n(E) how much of C_log survives -- the wall itself")
    print(f"{'N':>8} {'C(N)':>12} {'C/sqrt(N)':>11} {'C_log':>12} "
          f"{'C_log/(N logN)':>15} {'S_abs/(N logN)':>15}")
    for (N, C, C_log, M_log, ps, D, M, V, lp), (n2, S_abs, S_null, triv) \
            in zip(rows, meas):
        print(f"{N:>8} {C:>12.2f} {C/math.sqrt(N):>11.4f} "
              f"{C_log:>12.2f} {C_log/(N*math.log(N)):>15.6f} "
              f"{S_abs/(N*math.log(N)):>15.4f}")

    print("\nreading: (C) is the verdict. If S_abs/triv tracks the null")
    print("exponent, per-p square-root cancellation alone would give")
    print("C_log = o(N log N) with sqrt(log N) to spare, and the split")
    print("is the first transform in this campaign with positive margin.")
    print("If S_abs/triv is flat, the split fragments the wall without")
    print("creating room, and it joins the closed classes.")
    print("DONE")


if __name__ == "__main__":
    main()
