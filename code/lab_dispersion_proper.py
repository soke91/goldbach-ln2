# -*- coding: utf-8 -*-
"""
Transform Lab, session 12 (increment 232): the dispersion method proper
-- subtract the main term, bound only the variance.

Session 11 showed that plain Cauchy-Schwarz on the range-grouped form
lands at exactly the demand (Sum_R CS_R = 1.45 N, i.e. 1.83 x the
trivial bound), while the DIAGONAL of the expansion sits at 0.57 N and
falls. The damage is the positive off-diagonal main terms, which are
explicit. So the specification was: do not bound |Sum|^2, compute it.

THE ARGUMENT THIS SCRIPT TESTS. Write

    C_R = -Sum_w mu(w) H_{R,w},
    H_{R,w} = Sum_{p in R, p not| w} log p Lambda(N-pw)/log(pw),

and let h_w be the MAIN TERM of H_{R,w}. For fixed w, as p runs over
R = [P, 2P) the sum Sum_{p in R} log p Lambda(N-pw) counts the
representations N = wp + n in primes with p in R, whose main term is
S_w(N) x |R| with S_w(N) the singular series of that binary equation --
the same S_w that appeared in session 7. Pulling out the smooth weight,

    h_w := kappa_R * S_w(N) / log(P w),

with kappa_R fitted by weighted least squares from the same data (one
parameter; the shape S_w/log(Pw) is derived, only the constant is
fitted). Then split rather than bound:

    C_R = -Sum_w mu(w) h_w  -  Sum_w mu(w) (H_{R,w} - h_w),

  * the FIRST piece is a Mobius sum against an explicit multiplicative
    shape. Session 7 measured its analogue A(W) = Sum mu(w) f(w)/w at
    PNT strength, ~1e-3 against 1/log W ~ 0.09, so it is expected to be
    free. Measured here directly, not assumed.
  * the SECOND piece is where Cauchy-Schwarz is now applied, and it
    costs sqrt(n_sqf * Var_R) with Var_R = Sum_w mu^2(w)(H-h)^2. This is
    the variance, and the whole question is how much smaller it is than
    the raw second moment Sum_w mu^2(w) H^2 that session 11 used.

CRITERION, with the references on the same line. The wall needs
Sum_R |C_R| = o(N).
    NO MARGIN iff Sum_R disp_R / N >= 1, as plain Cauchy-Schwarz gave
              (1.429, 1.450, 1.446).
    MARGIN    iff Sum_R disp_R / N is below 1 and falling, where
              disp_R := |first piece| + sqrt(n_sqf * Var_R).
Also reported: the variance-explained fraction 1 - Var_R/Sum_w H^2, so
that a small disp_R can be traced to the main term actually fitting
rather than to an accounting slip.

HONEST LIMIT, stated before the numbers. Fitting kappa_R from the data
is legitimate for a main term whose SHAPE is derived, but it is one
free parameter per range, and with ~17 ranges that is 17 parameters
fitted against 17 reported numbers. The variance-explained column is
what guards against this: if the derived shape is right, it should
explain most of the second moment, and the fit should be doing very
little. If the shape is wrong, kappa absorbs nothing and Var stays near
the raw second moment.
"""
import numpy as np
import math

from lab_prime_factor_split import sieve

C2 = 0.6601618158468696


def prime_factors(n, spf):
    fs = []
    while n > 1:
        p = int(spf[n]); fs.append(p)
        while n % p == 0:
            n //= p
    return fs


def singular_w(N, W, spf, Nfac):
    """S_w(N) for the equation N = w p + n, w <= W. Zero unless w is
    odd and coprime to N (else n is forced into a fixed class)."""
    S = np.zeros(W + 1)
    C2N = 2 * C2
    for q in Nfac:
        if q > 2:
            C2N *= (q - 1) / (q - 2)
    for w in range(1, W + 1):
        if w % 2 == 0:
            continue
        fs = prime_factors(w, spf) if w > 1 else []
        if any(q in Nfac for q in fs):
            continue
        v = C2N
        for q in fs:
            v *= (q - 1) / (q - 2)
        S[w] = v
    return S


def main():
    X = 200_000
    mu, lam, spf, primes = sieve(X)
    NS = (50_000, 100_000, 200_000)

    agg = []
    for N in NS:
        Nfac = set(prime_factors(N, spf))
        v = np.arange(1, N)
        triv = float((np.abs(mu[1:N].astype(np.float64)) * lam[N - v]).sum())
        ps = primes[primes < N]
        rows = []
        b = 2
        while b < N:
            hi = min(2 * b, N)
            sel = ps[(ps >= b) & (ps < hi)]
            if len(sel) == 0:
                b *= 2
                continue
            W = (N - 1) // b
            H = np.zeros(W + 1)
            for p in sel:
                p = int(p)
                wmax = min(W, (N - 1) // p)
                if wmax < 1:
                    continue
                ww = np.arange(1, wmax + 1)
                t = (math.log(p) * lam[N - p * ww]
                     / np.log((p * ww).astype(np.float64)))
                H[1:wmax + 1] += np.where(ww % p != 0, t, 0.0)

            Sw = singular_w(N, W, spf, Nfac)
            ww = np.arange(1, W + 1)
            g = Sw[1:] / np.log(np.maximum(b * ww, 2).astype(np.float64))
            muw = mu[1:W + 1].astype(np.float64)
            sq = (muw != 0).astype(np.float64)
            Hh = H[1:]
            den = float(np.dot(sq * g, g))
            kappa = float(np.dot(sq * Hh, g)) / den if den > 0 else 0.0
            h = kappa * g
            piece1 = -float(np.dot(muw, h))
            resid = Hh - h
            Var = float(np.dot(sq * resid, resid))
            raw = float(np.dot(sq * Hh, Hh))
            n_sqf = float(sq.sum())
            C_R = -float(np.dot(muw, Hh))
            disp = abs(piece1) + math.sqrt(n_sqf * Var)
            CS = math.sqrt(n_sqf * raw)
            rows.append((b, hi, len(sel), C_R, piece1, Var, raw, n_sqf,
                         disp, CS, kappa))
            b *= 2
        agg.append((N, triv, rows))

    print("(A) per range at N = 200000 -- main term fitted, then split")
    N, triv, rows = agg[-1]
    print(f"{'p range':>16} {'#w':>8} {'|C_R|':>10} {'|piece1|':>10} "
          f"{'var expl':>9} {'sqrt(nV)':>11} {'disp_R':>11} "
          f"{'disp/CS':>8}")
    for (b, hi, npz, C_R, p1, Var, raw, n, disp, CS, kap) in rows:
        flag = "  <- overfit" if n < 40 else ""
        print(f"{b:>7}-{hi:>8} {int(n):>8} {abs(C_R):>10.1f} "
              f"{abs(p1):>10.1f} {1-Var/raw:>9.4f} "
              f"{math.sqrt(n*Var):>11.1f} {disp:>11.1f} "
              f"{disp/CS:>8.4f}{flag}")
    print("    #w is the number of squarefree w the one-parameter fit")
    print("    sees. Where it is small the fit interpolates and the")
    print("    reported var-expl of 1.0000 means overfitting, not a")
    # verdict-ok: structural: overfitting at small #w, independent of the numbers
    print("    good model -- those rows prove nothing.")

    print("\n(B) the verdict -- does the dispersion method leave margin?")
    print(f"{'N':>8} {'Sum|C_R|':>10} {'Sum disp':>11} {'Sum CS':>11} "
          f"{'disp/N':>8} {'CS/N':>7} {'disp/CS':>8} {'Sum|p1|/N':>10}")
    for (N, triv, rows) in agg:
        c = sum(abs(r[3]) for r in rows)
        d = sum(r[8] for r in rows)
        s = sum(r[9] for r in rows)
        p1 = sum(abs(r[4]) for r in rows)
        print(f"{N:>8} {c:>10.1f} {d:>11.1f} {s:>11.1f} {d/N:>8.3f} "
              f"{s/N:>7.3f} {d/s:>8.4f} {p1/N:>10.4f}")
    print("    the wall needs Sum_R |C_R| = o(N). disp/N >= 1 means the")
    print("    dispersion method has spent the budget too")

    print("\n(C) how much of the second moment the derived shape explains")
    print(f"{'N':>8} {'mean var expl':>14} {'min':>8} {'max':>8} "
          f"{'Var/raw total':>14}")
    for (N, triv, rows) in agg:
        ve = [1 - r[5] / r[6] for r in rows]
        tv = sum(r[5] for r in rows) / sum(r[6] for r in rows)
        print(f"{N:>8} {float(np.mean(ve)):>14.4f} {min(ve):>8.4f} "
              f"{max(ve):>8.4f} {tv:>14.4f}")
    print("    if this is near zero the fitted main term is absorbing")
    print("    nothing and disp_R is not really a dispersion bound")
    print("DONE")


if __name__ == "__main__":
    main()
