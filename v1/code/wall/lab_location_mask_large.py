# -*- coding: utf-8 -*-
"""
Does the location mask survive at large N? (increment 241)

Increment 240 established that the mean of Z(N) = C(N)/sqrt(V(N)) is a
deterministic function of which small primes divide N, and recorded the
limitation that matters: the mask's DEEPEST cells -- N divisible by many
small primes -- are exactly its LEAST POPULATED ones. At X = 4*10^6
there are only seven multiples of 510510 = 2*3*5*7*11*13*17 and none at
all of 9699690 = 2*3*5*7*11*13*17*19. A mean estimated from seven points
is not a measurement of a mask; it is a rumour about one.

Sweeping the whole range by FFT cannot fix this -- memory grows with X
and the deep cells stay rare. But nothing forces a sweep. C(N) and V(N)
for a SINGLE N cost one pass over the primes below N,

    C(N) = Sum_{p^k < N} log p * mu(N - p^k),
    V(N) = Sum_{p^k < N} (log p)^2 * mu^2(N - p^k),

so a few hundred targeted N can be taken far beyond where a sweep
stops. That is what this does.

DESIGN.
  * DEEP:      N = k * 510510, every k with N <= X.
  * DEEPER:    N = k * 9699690 (adds 19), every k with N <= X.
  * DEEPEST:   N = k * 223092870 (adds 23) if any fit.
  * CONTROL:   even N drawn at random from the top decade, giving the
    local mean and sd of Z at this scale, against which the deep cells
    are read. Without it there is no null, only a number.
  * A SHALLOW comparison group: N = 2*q with q prime, the emptiest
    possible cell, to check that the mask's other end behaves too.

NULLS AND CRITERION, on the same line.
  * The control sample supplies mean(Z) and sd(Z) at this N. NULL for
    a deep cell is the control mean, and the deviation is quoted in
    control sd units.
  * PREDICTION, made before running, from increment 240 at N ~ 10^6:
    the 510510 cell sits near -8.5 to -9 in units where the control has
    sd 1. CONFIRMED iff the deep cells at X = 3*10^7 still sit within
    about 2 of that, i.e. the mask does not wash out with N.
  * REFUTED iff the deep-cell mean drifts toward the control mean as N
    grows, which would make increment 240 a small-N artefact.
  * Reported alongside: C(N)/sqrt(N) and C(N)/N, so it stays visible
    that the effect is a sqrt(N)-scale main term and not a threat to
    C(N) = o(N).
"""
import numpy as np
import math
import sys
import time


def sieve_mu_and_primepowers(X):
    """Vectorised Mobius up to X, plus the list of prime powers <= X
    with their log p."""
    mu = np.ones(X + 1, dtype=np.int8)
    isp = np.ones(X + 1, dtype=bool); isp[:2] = False
    r = int(X ** 0.5)
    for p in range(2, r + 1):
        if isp[p]:
            isp[p * p::p] = False
            mu[p::p] = -mu[p::p]
            mu[p * p::p * p] = 0
    val = np.arange(X + 1, dtype=np.int32)
    for p in range(2, r + 1):
        if isp[p] or (p <= r and val[p] == p):
            pass
    # recompute primality properly (isp above was destroyed for p<=r)
    isp = np.ones(X + 1, dtype=bool); isp[:2] = False
    for p in range(2, r + 1):
        if isp[p]:
            isp[p * p::p] = False
    for p in np.nonzero(isp[: r + 1])[0]:
        val[int(p)::int(p)] //= int(p)
    mu[val > 1] = -mu[val > 1]
    mu[0] = 0
    primes = np.nonzero(isp)[0].astype(np.int64)
    return mu, primes


def CV(N, mu, primes, lpr):
    """C(N) and V(N) by one pass over the prime powers below N."""
    j = int(np.searchsorted(primes, N))
    p = primes[:j]; lg = lpr[:j]
    m = mu[N - p].astype(np.float64)
    C = float(np.dot(lg, m))
    V = float(np.dot(lg * lg, np.abs(m)))
    # prime powers p^k, k >= 2
    for q in primes[: int(np.searchsorted(primes, int(N ** 0.5)) + 1)]:
        q = int(q); e = q * q
        while e < N:
            mm = float(mu[N - e])
            C += math.log(q) * mm
            V += math.log(q) ** 2 * abs(mm)
            e *= q
    return C, V


def main():
    X = int(sys.argv[1]) if len(sys.argv) > 1 else 30_000_000
    t0 = time.time()
    mu, primes = sieve_mu_and_primepowers(X)
    lpr = np.log(primes.astype(np.float64))
    print(f"sieve to {X}: {len(primes)} primes, t={time.time()-t0:.0f}s",
          flush=True)

    P6 = 510510            # 2*3*5*7*11*13*17
    P7 = 9699690           # * 19
    P8 = 223092870         # * 23
    rng = np.random.default_rng(20260806)

    groups = {}
    def thin(vals, cap):
        if len(vals) <= cap:
            return vals
        step = len(vals) // cap
        return vals[::step][:cap]
    groups["deep  510510"] = thin(
        [k * P6 for k in range(1, X // P6 + 1)], 60)
    groups["deeper 9699690"] = thin(
        [k * P7 for k in range(1, X // P7 + 1)], 40)
    if X >= P8:
        groups["deepest 223092870"] = thin(
            [k * P8 for k in range(1, X // P8 + 1)], 20)
    lo = X // 2
    nctl = 300 if X <= 40_000_000 else 150
    ctrl = np.unique(rng.integers(lo // 2, X // 2,
                                  size=2 * nctl) * 2)
    groups["control (random even)"] = [int(v) for v in ctrl[:nctl]]
    sh = []
    cand = primes[np.searchsorted(primes, X // 4):]
    for q in cand[:: max(1, len(cand) // 120)][:60]:
        if 2 * int(q) <= X:
            sh.append(2 * int(q))
    groups["shallow  2*q"] = sh

    res = {}
    for name, Ns in groups.items():
        zs = []
        for N in Ns:
            if N < 1000 or N > X:
                continue
            C, V = CV(int(N), mu, primes, lpr)
            zs.append((int(N), C, V, C / math.sqrt(V)))
        res[name] = zs
        print(f"  {name}: {len(zs)} values  t={time.time()-t0:.0f}s",
              flush=True)

    cz = np.array([r[3] for r in res["control (random even)"]])
    cm, cs = float(cz.mean()), float(cz.std())
    print(f"\ncontrol at N ~ {X}: mean {cm:+.4f}  sd {cs:.4f}  "
          f"n = {len(cz)}")
    print("  (this is the null every other row is read against)")

    print(f"\n{'group':>24} {'n':>5} {'mean Z':>9} {'sd':>7} "
          f"{'(mean-ctrl)/sd_ctrl':>20}")
    for name, zs in res.items():
        z = np.array([r[3] for r in zs])
        if len(z) == 0:
            continue
        print(f"{name:>24} {len(z):>5} {z.mean():>9.3f} "
              f"{z.std():>7.3f} {(z.mean()-cm)/cs:>20.2f}")

    print(f"\nthe deep cells one by one -- does it hold at every N?")
    print(f"{'N':>12} {'Z':>9} {'C':>14} {'C/sqrt(N)':>11} {'C/N':>10}")
    for name in ("deep  510510", "deeper 9699690",
                 "deepest 223092870"):
        if name not in res:
            continue
        print(f"  -- {name}")
        for (N, C, V, z) in res[name][:14]:
            print(f"{N:>12} {z:>9.3f} {C:>14.1f} "
                  f"{C/math.sqrt(N):>11.3f} {C/N:>10.5f}")

    # Z = C/sqrt(V) puts a sqrt(log N) in the denominator, because
    # V ~ kappa S N log N. If the deterministic term is really of size
    # sqrt(N), then Z must decay like 1/sqrt(log N) FOR THAT REASON
    # ALONE, and reading that decay as the mask weakening would be the
    # normalisation hazard a fourth time. Test the scale-free form.
    print(f"\n(E) which invariant is actually constant?")
    print(f"{'group':>24} {'n':>4} {'mean C/sqrt(N)':>15} "
          f"{'1st half':>10} {'2nd half':>10} {'mean Z':>9} "
          f"{'Z*sqrt(logN)':>13}")
    for name, zs in res.items():
        if len(zs) < 3:
            continue
        cn = np.array([r[1] / math.sqrt(r[0]) for r in zs])
        zz = np.array([r[3] for r in zs])
        zl = np.array([r[3] * math.sqrt(math.log(r[0])) for r in zs])
        h = len(cn) // 2
        print(f"{name:>24} {len(cn):>4} {cn.mean():>15.3f} "
              f"{cn[:h].mean():>10.3f} {cn[h:].mean():>10.3f} "
              f"{zz.mean():>9.3f} {zl.mean():>13.3f}")
    print("    if C/sqrt(N) is flat across the halves while Z is not,")
    print("    the invariant is C/sqrt(N) and Z's decay is the")
    print("    normaliser's sqrt(log N), not the mask weakening")

    dz = np.array([r[3] for r in res["deep  510510"]])
    dev = (dz.mean() - cm) / cs
    print(f"\nverdict:")
    if dev < -6:
        print(f"  CONFIRMED -- at N up to {X} the 510510 cell still")
        print(f"  sits {dev:.1f} control sd below the null. The mask")
        print(f"  does not wash out with N.")
    elif dev < -2:
        print(f"  WEAKENED -- {dev:.1f} sd, against about -8.5 at")
        print(f"  N ~ 10^6. The mask survives but shrinks; the rate")
        print(f"  matters and needs a third scale.")
    else:
        print(f"  REFUTED -- {dev:.1f} sd. Increment 240 was a small-N")
        print(f"  artefact.")
    print("  note C/N in the table: the effect is a sqrt(N)-scale main")
    print("  term, so it is o(N) and does not threaten the wall.")
    print("DONE")


if __name__ == "__main__":
    main()
