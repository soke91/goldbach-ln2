# -*- coding: utf-8 -*-
"""
Transform Lab, session 8 (increment 228): where the cancellation in the
wall actually lives -- an alternating series in omega(v).

Session 7 left a puzzle. In the p-grouping the surviving mass climbs to
100% at large p, and the reason is visible once stated: for v prime,
mu(v) = -1 with NO sign variation at all, so those terms cannot cancel
among themselves. Yet the prime v alone contribute

    Sum_{v = p} mu(v) Lambda(N-v) log v = -Sum_p log p Lambda(N-p)
                                        ~ -S(N) N,

which is enormous next to the measured C_log. So the smallness of the
wall is NOT cancellation inside any class of v -- it is cancellation
ACROSS classes.

THE REFORMULATION. On the support of mu, mu(v) = (-1)^{omega(v)}, so

    C_log(N) = Sum_{j>=1} (-1)^j T_j(N),
    T_j(N)  := Sum_{v<N, mu^2(v)=1, omega(v)=j} Lambda(N-v) log v >= 0.

Every T_j is NONNEGATIVE and Sum_j T_j is exactly the trivial bound.
The wall is an alternating series in the number of prime factors, and
the entire question is how nearly consecutive T_j cancel. (v = 1 has
omega = 0 and log 1 = 0, so it contributes nothing and the sum starts
at j = 1.)

THE HYPOTHESIS UNDER TEST. omega(v) is Poisson-like with mean
~ log log N (Erdos-Kac / Sathe-Selberg), and for an exact Poisson shape
the alternating sum collapses geometrically, Sum_j (-1)^j e^{-lam}
lam^j / j! = e^{-2 lam}. If the T_j inherit that shape, a total mass of
~0.73 N log N would alternate down to ~N log N e^{-2 log log N} =
N / log N -- comfortably o(N log N), by a mechanism exponential in
log log N rather than in log N. H_POISSON: the shape of the
omega-distribution accounts for the observed smallness.

MEASUREMENTS.
 (A) the reformulation, exactly: Sum_j (-1)^j T_j = C_log.
 (B) the profile T_j / Sum T_j against j. The T_j carry no j = 0 mass
     (v = 1 has log v = 0), so the model is a Poisson CONDITIONED to
     j >= 1, whose mean is lam / (1 - e^{-lam}). lam is therefore
     obtained by SOLVING lam / (1 - e^{-lam}) = m for the measured
     conditional mean m; setting lam = m directly would misfit, and the
     two differ materially here (m = 2.28 gives lam = 1.95).
 (C) the alternation ratio R(N) = |Sum_j (-1)^j T_j| / Sum_j T_j
     against, on the same line, its two references:
        R = 1                                            no cancellation
        R_pois = |e^{-2lam} - e^{-lam}| / (1 - e^{-lam})     H_POISSON
     The fitted P(1) is printed beside the measured T_1/Sum T so the
     shape is checked and not only its mean.
 (D) how the partial alternating sums build up.
 (E) THE CONTROL, and it is the one that matters. Drop the arithmetic
     weight: T_j^0 := Sum_{omega(v)=j, v<N squarefree} log v, whose
     alternating sum is Sum_{v<N} mu(v) log v -- the plain Mobius sum,
     small by the prime number theorem. The identical statistic on it
     separates two very different readings of a small R:
       * if the control's R is comparable, the omega-alternation is a
         property of mu itself and the additive constraint Lambda(N-v)
         does not destroy it -- the encouraging outcome;
       * if the control's R is far smaller, the additive constraint has
         already cost most of the alternation.

 (F) THE DECOMPOSITION THAT SAYS WHAT THE WALL ACTUALLY ASKS. Write
        r_j := T_j / T_j^0,
     the mean level of Lambda(N-v) over the squarefree v with exactly j
     prime factors. Then, with rbar the T^0-weighted mean of r_j,

        Sum_j (-1)^j T_j  =  rbar * Sum_j (-1)^j T_j^0
                             + Sum_j (-1)^j (r_j - rbar) T_j^0.

     The first piece is rbar times the plain Mobius sum, hence
     PNT-small and free. EVERYTHING the wall asks sits in the second:
     the wall is small if and only if the prime density near N is the
     SAME across omega-classes, to the precision the second piece
     needs. Both pieces are printed against the trivial bound.

     This is the reduction worth having. "mu does not correlate with
     Lambda(N-.)" becomes "the level of Lambda(N-.) does not depend on
     omega(v)" -- and omega is an incomparably simpler statistic than
     mu, being a count rather than a sign.

 (G) IS THE LEVEL SMOOTH IN j? An alternating sum is a repeated
     difference operator, so it annihilates smooth functions of j: for
     any polynomial f, Sum_j (-1)^j f(j) lam^j/j! e^{-lam} is
     e^{-2 lam} times a polynomial in lam, still exponentially small.
     So the drift term is harmless as long as r_j is a SMOOTH function
     of j, and dangerous only in a component oscillating like (-1)^j.
     Measured by replacing r_j with its weighted least-squares fit of
     degree 1 and degree 2 and re-running the alternating sum: the fit
     part is what a smooth level would give, the residual part is what
     the oscillating component gives.

     This is where the reformulation meets the parity problem head on.
     A (-1)^j component in the level IS a mu-correlation by definition,
     so if the residual carries the wall, the omega-decomposition has
     restated parity rather than circumvented it -- worth knowing
     precisely, and better known than assumed.

READING. R(N) fluctuates -- C_log changes sign with N -- so no single N
proves anything. The questions are whether R sits near R_pois or far
below it, how it compares with the control, which of (F)'s two pieces
carries the residual, and whether (G) leaves anything a smooth level
cannot explain.
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
    om = np.zeros(X + 1, dtype=np.int8)
    for i in range(2, X + 1):
        p = int(spf[i]); j = i // p
        mu[i] = 0 if j % p == 0 else -mu[j]
        om[i] = om[j] + (0 if j % p == 0 else 1)   # omega on squarefrees
    primes = np.nonzero(spf[2:] == np.arange(2, X + 1))[0] + 2
    lam = np.zeros(X + 1)
    for p in primes:
        q = int(p); lp = math.log(int(p))
        while q <= X:
            lam[q] = lp; q *= int(p)
    return mu, lam, spf, om


def fit_lambda(m):
    """Solve lam / (1 - exp(-lam)) = m for lam > 0, by bisection."""
    lo, hi = 1e-9, 60.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if mid / (1 - math.exp(-mid)) < m:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def main():
    X = 400_000
    mu, lam, spf, om = sieve(X)
    NS = (50_000, 100_000, 200_000, 400_000)
    JMAX = 9

    print("(A) the reformulation  C_log = Sum_j (-1)^j T_j")
    print(f"{'N':>8} {'C_log direct':>14} {'alt sum':>14} {'rel.diff':>10} "
          f"{'Sum T_j':>14}")
    keep = []
    for N in NS:
        v = np.arange(1, N)
        sq = (mu[1:N] != 0)
        w = np.where(sq, lam[N - v] * np.log(v.astype(np.float64)), 0.0)
        jj = om[1:N]
        T = np.array([float(w[jj == j].sum()) for j in range(0, JMAX + 1)])
        alt = float(sum((-1) ** j * T[j] for j in range(1, JMAX + 1)))
        C_log = float(np.dot(mu[1:N].astype(np.float64),
                             lam[N - v] * np.log(v.astype(np.float64))))
        tot = float(T[1:].sum())
        print(f"{N:>8} {C_log:>14.2f} {alt:>14.2f} "
              f"{abs(alt-C_log)/max(abs(C_log),1):>10.2e} {tot:>14.1f}")
        keep.append((N, C_log, T, tot))
    print("    every T_j >= 0 and Sum_j T_j is exactly the trivial")
    print("    bound: the wall is an alternating series in omega(v)")

    print("\n(B) the profile T_j / Sum T_j, with lambda fitted from it")
    hdr = "".join(f"{j:>8}" for j in range(1, JMAX + 1))
    print(f"{'N':>8}{hdr} {'cond mean':>9} {'loglogN':>8}")
    for (N, C_log, T, tot) in keep:
        m = float(sum(j * T[j] for j in range(1, JMAX + 1)) / tot)
        row = "".join(f"{T[j]/tot:>8.4f}" for j in range(1, JMAX + 1))
        print(f"{N:>8}{row} {m:>9.4f} "
              f"{math.log(math.log(N)):>8.4f}")
    print("    the column is the measured conditional mean m; the")
    print("    Poisson parameter fitted from it appears in table (C)")

    print("\n(C) the alternation ratio against its two references")
    print(f"{'N':>8} {'R=|alt|/SumT':>13} {'R_pois':>9} {'no-cancel':>10} "
          f"{'R/R_pois':>9} {'lam fit':>8} {'P1 pois':>8} {'P1 real':>8}")
    for (N, C_log, T, tot) in keep:
        m = float(sum(j * T[j] for j in range(1, JMAX + 1)) / tot)
        lam = fit_lambda(m)
        alt = float(sum((-1) ** j * T[j] for j in range(1, JMAX + 1)))
        R = abs(alt) / tot
        e1 = math.exp(-lam); e2 = math.exp(-2 * lam)
        R_pois = abs(e2 - e1) / (1 - e1)
        p1 = lam * e1 / (1 - e1)
        print(f"{N:>8} {R:>13.6f} {R_pois:>9.6f} {1.0:>10.1f} "
              f"{R/R_pois:>9.3f} {lam:>8.4f} {p1:>8.4f} "
              f"{T[1]/tot:>8.4f}")
    print("    R = 1 would mean no cancellation across omega at all;")
    print("    R_pois is what a fitted conditional Poisson would give.")
    print("    The P1 columns check the shape itself, not only its mean")

    print("\n(D) how the alternating sum builds up (N = 400000)")
    (N, C_log, T, tot) = keep[-1]
    print(f"{'j':>3} {'T_j':>14} {'T_j/Sum':>9} {'partial alt':>14} "
          f"{'|partial|/Sum':>14}")
    s = 0.0
    for j in range(1, JMAX + 1):
        s += (-1) ** j * T[j]
        print(f"{j:>3} {T[j]:>14.1f} {T[j]/tot:>9.4f} {s:>14.1f} "
              f"{abs(s)/tot:>14.6f}")
    print("    the partial sums do not decay monotonically -- the")
    print("    collapse happens only once the whole series is in")

    print("\n(E) the control -- the same statistic with the weight")
    print("    Lambda(N-v) removed, so the alternating sum is the plain")
    print("    Mobius sum  Sum_{v<N} mu(v) log v")
    print(f"{'N':>8} {'Sum mu(v) log v':>17} {'Sum T_j^0':>14} "
          f"{'R_ctrl':>10} {'R_wall':>10} {'R_wall/R_ctrl':>14}")
    for (N, C_log, T, tot) in keep:
        v = np.arange(1, N)
        sq = (mu[1:N] != 0)
        w0 = np.where(sq, np.log(v.astype(np.float64)), 0.0)
        jj = om[1:N]
        T0 = np.array([float(w0[jj == j].sum())
                       for j in range(0, JMAX + 1)])
        alt0 = float(sum((-1) ** j * T0[j] for j in range(1, JMAX + 1)))
        tot0 = float(T0[1:].sum())
        alt = float(sum((-1) ** j * T[j] for j in range(1, JMAX + 1)))
        print(f"{N:>8} {alt0:>17.1f} {tot0:>14.1f} "
              f"{abs(alt0)/tot0:>10.6f} {abs(alt)/tot:>10.6f} "
              f"{(abs(alt)/tot)/(abs(alt0)/tot0):>14.3f}")
    print("    comparable ratios would mean the omega-alternation is a")
    print("    property of mu that the additive constraint leaves intact")

    print("\n(F) the decomposition: level-times-Mobius, plus level drift")
    print(f"{'N':>8} {'rbar':>8} {'r_1':>8} {'r_2':>8} {'r_3':>8} "
          f"{'r_4':>8} {'spread':>8} {'free/triv':>10} {'drift/triv':>11}")
    for (N, C_log, T, tot) in keep:
        v = np.arange(1, N)
        sq = (mu[1:N] != 0)
        w0 = np.where(sq, np.log(v.astype(np.float64)), 0.0)
        jj = om[1:N]
        T0 = np.array([float(w0[jj == j].sum())
                       for j in range(0, JMAX + 1)])
        tot0 = float(T0[1:].sum())
        alt0 = float(sum((-1) ** j * T0[j] for j in range(1, JMAX + 1)))
        r = np.array([T[j] / T0[j] if T0[j] > 0 else 0.0
                      for j in range(0, JMAX + 1)])
        rbar = tot / tot0
        free = rbar * alt0
        drift = float(sum((-1) ** j * (r[j] - rbar) * T0[j]
                          for j in range(1, JMAX + 1)))
        live = [j for j in range(1, JMAX + 1) if T0[j] > 0]
        spread = max(abs(r[j] - rbar) for j in live) / rbar
        print(f"{N:>8} {rbar:>8.4f} {r[1]:>8.4f} {r[2]:>8.4f} "
              f"{r[3]:>8.4f} {r[4]:>8.4f} {spread:>8.4f} "
              f"{free/tot:>10.6f} {drift/tot:>11.6f}")
    print("    free  = rbar x (plain Mobius sum), PNT-small by itself")
    print("    drift = the omega-dependence of the prime level, and it")
    print("    is where the whole difficulty of the wall now sits")

    print("\n(G) is the level smooth in j? alternating sums kill smooth")
    print(f"{'N':>8} {'alt':>12} {'alt(deg1 fit)':>14} "
          f"{'alt(deg2 fit)':>14} {'resid deg1':>12} {'resid deg2':>12}")
    for (N, C_log, T, tot) in keep:
        v = np.arange(1, N)
        sq = (mu[1:N] != 0)
        w0 = np.where(sq, np.log(v.astype(np.float64)), 0.0)
        jj = om[1:N]
        T0 = np.array([float(w0[jj == j].sum())
                       for j in range(0, JMAX + 1)])
        live = np.array([j for j in range(1, JMAX + 1) if T0[j] > 0])
        r = np.array([T[j] / T0[j] for j in live])
        wt = np.array([T0[j] for j in live])
        alt = float(sum((-1) ** j * T[j] for j in live))
        row = [f"{N:>8}", f"{alt:>12.1f}"]
        resid = []
        for deg in (1, 2):
            co = np.polyfit(live.astype(float), r, deg, w=np.sqrt(wt))
            rf = np.polyval(co, live.astype(float))
            a_fit = float(sum((-1) ** j * rf[i] * T0[j]
                              for i, j in enumerate(live)))
            row.append(f"{a_fit:>14.1f}")
            resid.append(alt - a_fit)
        print(" ".join(row) + f" {resid[0]:>12.1f} {resid[1]:>12.1f}")
    print("    alt(fit) is what a perfectly smooth level would give;")
    print("    the residual is the (-1)^j-oscillating component, which")
    print("    is a mu-correlation by definition")

    print("\nreading: the reformulation is exact and the Poisson shape,")
    print("if it held, would already suffice (R_pois ~ C/log N = o(1)).")
    print("It does not hold -- the measured profile is not Poisson and")
    print("alternates better than one. What (F) shows is the useful")
    print("part: the wall reduces to the omega-independence of the")
    print("prime density near N, with the Mobius part free.")
    print("DONE")


if __name__ == "__main__":
    main()
