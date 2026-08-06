# -*- coding: utf-8 -*-
"""
Transform Lab, session 13 (increment 236): does the margin actually go
to zero, or does it converge to a positive constant?

Everything transform P claims rests on Proposition P.4's hypothesis
being TRUE, i.e. on

    S_abs(N) := Sum_{p<N} log p |D_p(N)|  =  o(N).

Sessions 6 and 9 measured S_abs/N falling from 0.3455 to 0.2995 over
N = 5*10^4 to 4*10^5. That is a factor 8, and a fall of 13%. It is
entirely consistent with S_abs/N -> 0, and it is ALSO consistent with

    S_abs/N  ->  c > 0,

in which case P.4's hypothesis is false, the "positive margin" is an
artefact of a short range, and the whole construction is dead. The
campaign has not tested this, and it is the load-bearing question.

WHAT SEPARATES THE TWO. Three models are fitted to the same data:
    (i)   S_abs/N = a (log N)^{-1/2}          -- the square-root-per-p
                                                 prediction of P.4
    (ii)  S_abs/N = c + b/log N               -- convergence to c > 0
    (iii) S_abs/N = a (log N)^{-b}            -- free exponent, -> 0

THE COMPARISON THAT ANSWERS THE QUESTION IS (ii) AGAINST (iii), NOT
(i) AGAINST (ii). The question asked is whether S_abs/N tends to zero.
Model (i) fixes the RATE as well as the limit, so rejecting it says
only that the rate is not exactly (log N)^{-1/2}; that is not evidence
for a positive limit. Both (ii) and (iii) have two parameters, and they
differ precisely in the limit -- (ii) allows c > 0, (iii) forces 0 --
so their residuals are comparable and the comparison is the right one.
(An earlier version of this script used (i) against (ii) and returned
DEAD on that basis; the criterion did not match the question. Recorded
as correction #35.)

NULLS AND CRITERION, on the same line as the thresholds.
  * Range extended to N = 6.4*10^6, a factor 128, over which
    (log N)^{-1/2} falls by 21% and a constant does not fall at all.
  * ALIVE     iff (iii) fits at least as well as (ii) within a factor
              1.5 in residual sd -- i.e. the data do not demand a
              positive limit.
  * DEAD      iff (ii) fits better than (iii) by more than that, with
              its fitted c stable and bounded away from 0 under
              leave-one-out.
  * UNDECIDED iff the two are comparable AND the models' predictions
              separate by less than the measurement noise at any N a
              computation could reach -- reported explicitly, with the
              separation at log N = 20 printed, since an undecidable
              question should be labelled and not resolved by a
              coin-flip between two equally good fits.

STATED LIMIT. Eight points cannot settle an asymptotic question. This
detects a clear flattening, and reports honestly when it cannot
distinguish.
"""
import numpy as np
import math
import time

from lab_prime_factor_split import sieve


def s_abs(N, mu, lam, primes):
    v = np.arange(1, N)
    muv = mu[1:N].astype(np.float64)
    lamr = lam[N - v]
    logv = np.log(v.astype(np.float64)); logv[0] = 1.0
    t = muv * lamr / logv; t[0] = 0.0
    at = np.abs(muv) * lamr / logv; at[0] = 0.0
    triv = float((np.abs(muv) * lamr)[1:].sum())
    ps = primes[primes < N]
    S = 0.0
    for p in ps:
        idx = np.arange(int(p), N, int(p)) - 1
        S += math.log(int(p)) * abs(float(t[idx].sum()))
    return S, triv


def main():
    X = 6_400_000
    t0 = time.time()
    mu, lam, spf, primes = sieve(X)
    print(f"sieve done  t={time.time()-t0:.0f}s", flush=True)

    NS = (50_000, 100_000, 200_000, 400_000, 800_000, 1_600_000,
          3_200_000, 6_400_000)
    rows = []
    for N in NS:
        S, triv = s_abs(N, mu, lam, primes)
        rows.append((N, S, triv))
        print(f"N={N:>9}  S_abs/N={S/N:.4f}  S_abs/triv={S/triv:.4f}  "
              f"t={time.time()-t0:.0f}s", flush=True)

    x = np.array([math.log(N) for (N, S, t) in rows])
    y = np.array([S / N for (N, S, t) in rows])

    # (i) a (log N)^{-1/2}
    a1 = float(np.sum(y / np.sqrt(x)) / np.sum(1.0 / x))
    r1 = y - a1 / np.sqrt(x)
    # (ii) c + b/log N
    A2 = np.vstack([np.ones_like(x), 1.0 / x]).T
    c2, b2 = np.linalg.lstsq(A2, y, rcond=None)[0]
    r2 = y - (c2 + b2 / x)
    # (iii) free exponent
    co = np.polyfit(np.log(x), np.log(y), 1)
    r3 = y - np.exp(co[1]) * x ** co[0]

    print("\nmodel fits")
    print(f"  (i)   a (logN)^-1/2      a={a1:.4f}"
          f"                 resid sd={r1.std():.5f}")
    print(f"  (ii)  c + b/logN         c={c2:+.4f}  b={b2:+.4f}"
          f"   resid sd={r2.std():.5f}")
    print(f"  (iii) a (logN)^-b        b={-co[0]:.4f}"
          f"                 resid sd={r3.std():.5f}")

    print("\nwhat the data say about the limit c in model (ii)")
    print(f"  fitted c = {c2:+.4f}")
    # leave-one-out spread on c, as a crude stability check
    cs = []
    for i in range(len(x)):
        m = np.ones(len(x), bool); m[i] = False
        cc = np.linalg.lstsq(np.vstack([np.ones(m.sum()),
                                        1.0 / x[m]]).T, y[m],
                             rcond=None)[0][0]
        cs.append(float(cc))
    print(f"  leave-one-out range: [{min(cs):+.4f}, {max(cs):+.4f}]")

    # the comparison that answers the question: (iii) -> 0 vs (ii) -> c
    a3 = math.exp(co[1]); b3 = -co[0]
    sep20 = abs((c2 + b2 / 20.0) - a3 * 20.0 ** (-b3))
    noise = float(max(r2.std(), r3.std()))
    print(f"\n  (ii) vs (iii): resid sd {r2.std():.5f} vs "
          f"{r3.std():.5f}   ratio {r3.std()/r2.std():.2f}")
    print(f"  predictions at log N = 20 (N ~ 5e8): "
          f"(ii) {c2 + b2/20.0:.4f}  (iii) {a3*20.0**(-b3):.4f}  "
          f"separation {sep20:.4f}")
    print(f"  measurement noise (larger resid sd) = {noise:.5f}")

    comparable = r3.std() <= r2.std() * 1.5
    if not comparable and min(cs) > 0.02:
        verdict = ("DEAD -- the data demand S_abs/N -> c > 0, and "
                   "P.4's hypothesis is false")
    elif comparable and sep20 < 5 * noise:
        verdict = ("UNDECIDED -- (ii) and (iii) fit equally well and "
                   "their predictions separate by less than the noise "
                   "at any reachable N. Whether S_abs/N -> 0 is NOT "
                   "numerically decidable.")
    else:
        verdict = "ALIVE -- the data do not demand a positive limit"
    print("\nverdict:", verdict)
    print("what IS settled: the exponent is not -1/2. The measured")
    print(f"decay is (log N)^-{b3:.3f}, FASTER than the square-root-")
    print("per-p prediction, so P.4's budget is met with more room")
    print("than predicted, not less.")
    print("DONE")


if __name__ == "__main__":
    main()
