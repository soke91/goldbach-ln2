# -*- coding: utf-8 -*-
"""
Which decay regime does rho actually sit in, and does Theorem D get
stronger or weaker when the input is strengthened? (increment 278)

BACKGROUND. Theorem D's no-go rests on one analytic input, Huang-Li's
Lemma 1 (Goldston-Yildirim):

    rho_n(x) = Sum_{j < x, (j,n)=1} mu(j)/phi(j)  <<  exp(-c sqrt(log x))

and the no-go's strength is the ratio ||b||_1 / |B_w|, which is exactly
1 / max_d |rho|. So the quality of the bound on rho IS the quality of
the theorem, and one question decides whether Theorem D is a real
obstruction or an artifact of a weak unconditional estimate:

    is exp(-c sqrt(log x)) the truth, or only what the classical
    zero-free region can prove?

It is only what the classical zero-free region can prove. The Dirichlet
series has 1/zeta to the FIRST power (THEOREM_A.md, ingredient 3), so
under RH the same argument gives rho_n(x) << x^{-1/2+eps}, a power
saving. Substituting that into Theorem D's inequality replaces
exp(c sqrt(log N)) by a POWER of N. If that is right, Theorem D is not
fragile: strengthening the hypothesis strengthens the no-go. That is
the opposite of how an artifact behaves, and it is worth knowing,
because every other route in this campaign has closed on parity and one
would like to know whether this one closed for a real reason.

PRE-REGISTRATION (written before the run; hazard 2 is that a threshold
gets chosen after seeing the data, hazard 4 is that a null gets set by
a size heuristic instead of from data).

  Two competing shapes for the observed decay of |rho_n(x)| in x:
      H_exp    log|rho| = a - c*sqrt(log x)      (classical)
      H_pow    log|rho| = a - beta*log x         (RH-type power saving)
  These are not nested and are distinguished by which regression is
  straighter. DECISION RULE, fixed now:
      - fit both by least squares on the same points,
      - report R^2 for each and the residual RMS,
      - declare a regime only if one model's residual RMS is below half
        the other's; otherwise declare INDETERMINATE.
  A fitted beta is reported with its spread across the n tested, and no
  claim is made that a measured beta at x <= 2*10^7 identifies the true
  exponent -- rho is a tail of a conditionally convergent series and
  the asymptotic regime may lie far beyond reach. This is stated in
  advance so that a clean-looking fit is not later read as more than it
  is (hazard 5: a trend toward a null is not evidence of the null).

WHAT THIS CANNOT DO. It cannot prove the RH-conditional bound; that is
a substitution into a published argument, not a measurement. It can
only say which regime the accessible range is consistent with, and
whether the measured decay is at least as fast as the unconditional
bound -- which is the part that matters for Theorem D, since a FASTER
decay makes the no-go stronger.
"""
import numpy as np
import math
import time


def sieve_mu_phi(X):
    """mu and phi on [0, X] by smallest-prime-factor recursion."""
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
    phi = np.zeros(X + 1, dtype=np.int64)
    phi[1] = 1
    for i in range(2, X + 1):
        p = int(spf[i])
        j = i // p
        mu[i] = 0 if j % p == 0 else -mu[j]
        phi[i] = phi[j] * (p if j % p == 0 else p - 1)
    return mu, phi, spf


def radset(n, spf):
    s = []
    while n > 1:
        p = int(spf[n])
        s.append(p)
        while n % p == 0:
            n //= p
    return s


def fit(xs, ys):
    """least squares ys = a + b*xs; returns a, b, R^2, residual RMS."""
    xs = np.asarray(xs, dtype=np.float64)
    ys = np.asarray(ys, dtype=np.float64)
    b, a = np.polyfit(xs, ys, 1)
    pred = a + b * xs
    resid = ys - pred
    ss_res = float(np.dot(resid, resid))
    ss_tot = float(np.dot(ys - ys.mean(), ys - ys.mean()))
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else float('nan')
    return a, b, r2, math.sqrt(ss_res / len(xs))


def main():
    X = 20_000_000
    t0 = time.time()
    mu, phi, spf = sieve_mu_phi(X)
    print(f"sieve  t={time.time()-t0:.0f}s", flush=True)

    # rho_n needs (j, n) = 1. N is even in the Goldbach setting and the
    # modulus in Theorem D is dN, so 2 always divides it; n = 2 is the
    # shallowest relevant case and 30030 a deep one.
    NS = [2, 6, 30, 2310, 30030]

    # sample points, geometric, fixed before the run
    XS = [10 ** k for k in range(3, 8)]
    XS += [3 * 10 ** k for k in range(3, 7)]
    XS = sorted(x for x in XS if x <= X)

    term = np.zeros(X + 1, dtype=np.float64)
    nz = np.nonzero(mu[1:] != 0)[0] + 1
    term[nz] = mu[nz].astype(np.float64) / phi[nz].astype(np.float64)

    print(f"\n{'n':>7} {'x':>11} {'rho_n(x)':>14} {'|rho|':>11}")
    curves = {}
    for n in NS:
        R = radset(n, spf)
        keep = np.ones(X + 1, dtype=bool)
        for q in R:
            keep[q::q] = False
        t = np.where(keep, term, 0.0)
        cs = np.cumsum(t)
        pts = []
        for x in XS:
            r = float(cs[x - 1])          # j < x
            pts.append((x, r))
            print(f"{n:>7} {x:>11} {r:>+14.6e} {abs(r):>11.4e}")
        curves[n] = pts
        print()

    print("=" * 68)
    print("regime discrimination, decision rule fixed before the run:")
    print("  declare a regime only if one model's residual RMS is below")
    print("  half the other's; otherwise INDETERMINATE")
    print("=" * 68)
    print(f"\n{'n':>7} {'RMS exp':>10} {'RMS pow':>10} {'R2 exp':>8} "
          f"{'R2 pow':>8} {'beta':>7}  verdict")
    betas = []
    verdicts = []
    for n in NS:
        pts = [(x, r) for x, r in curves[n] if r != 0.0]
        if len(pts) < 4:
            print(f"{n:>7}  too few nonzero points")
            continue
        lx = [math.log(x) for x, _ in pts]
        ly = [math.log(abs(r)) for _, r in pts]
        sx = [math.sqrt(v) for v in lx]
        _, ce, r2e, rmse = fit(sx, ly)     # H_exp
        _, be, r2p, rmsp = fit(lx, ly)     # H_pow
        if rmse < 0.5 * rmsp:
            v = "H_exp"
        elif rmsp < 0.5 * rmse:
            v = "H_pow"
        else:
            v = "INDETERMINATE"
        verdicts.append(v)
        betas.append(-be)
        print(f"{n:>7} {rmse:>10.4f} {rmsp:>10.4f} {r2e:>8.4f} "
              f"{r2p:>8.4f} {-be:>7.4f}  {v}")

    if betas:
        print(f"\nfitted power-law exponent beta across n: "
              f"mean {np.mean(betas):.4f}, "
              f"spread [{min(betas):.4f}, {max(betas):.4f}]")
        print("  (RH-type substitution into Theorem D would want "
              "beta ~ 1/2;")
        print("   a measured beta at x <= 2e7 does NOT identify the true")
        print("   exponent -- rho is a tail of a conditionally "
              "convergent series)")

    print("\n" + "=" * 68)
    print("what this means for Theorem D")
    print("=" * 68)
    print("Theorem D's loss factor is 1/max_d |rho_{dN}(K/d)|, with")
    print("K = N^theta' and b supported on d <= N^{theta'-1/2-delta},")
    print("so K/d >= N^{1/2+delta}. Substituting each regime:")
    for lbl, f in (("unconditional  exp(-c sqrt(log y))",
                    lambda L: math.exp(0.5 * math.sqrt(0.5 * L))),
                   ("RH-type        y^{-1/2}",
                    lambda L: math.exp(0.5 * 0.5 * L))):
        row = "  ".join(f"N=1e{k}: {f(k*math.log(10)):.3e}"
                        for k in (10, 50, 200))
        print(f"  {lbl:>34}   {row}")
    print("\n  The RH-type column is a POWER of N and the unconditional")
    print("  one is not. Strengthening the input makes the no-go")
    print("  stronger, which is the opposite of how an artifact of a")
    print("  weak estimate behaves.")
    print("DONE")


if __name__ == "__main__":
    main()
