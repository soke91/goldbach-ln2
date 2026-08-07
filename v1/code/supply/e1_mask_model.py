# -*- coding: utf-8 -*-
"""
Increment 150: closed-form mask model, blind test.

The factorization law (inc. 144-149) says the only structure in the
dispersion field is a local mask. Here the mask is COMPUTED A PRIORI,
per pair, by finite modular enumeration, and tested blind against the
observed zero-density and the observed |C|.

Model: term t_p = mu(N-pk) mu(N-pk') is zero iff q^2 | (N-pk) or
q^2 | (N-pk') for some prime q. For q <= QCUT, the density rho_q of
that event over p uniform in the units mod q^2 is enumerated exactly.
For q > QCUT, both arguments are ~ random: tail factor
prod_{q>QCUT} (1 - 2/q^2) (independence approximation), computed
numerically. Predicted survivor density: s_pred = prod_q (1 - rho_q)
* tail. Blind comparisons on fresh pairs:
  (1) s_pred vs observed support fraction (scatter, corr, mean ratio),
  (2) r1_raw = |C|/sqrt(n_p) vs 0.798 * sqrt(s_pred).
"""
import numpy as np, time
from math import gcd

QCUT = 50

def mobius_upto(X):
    mu = np.ones(X + 1, dtype=np.int8)
    pm = np.ones(X + 1, dtype=bool); pm[:2] = False
    for p in range(2, int(X ** 0.5) + 1):
        if pm[p]:
            pm[p*p::p] = False
            mu[p::p] *= -1
            mu[p*p::p*p] = 0
    val = np.arange(X + 1, dtype=np.int64)
    for p in range(2, int(X ** 0.5) + 1):
        if pm[p]:
            val[p::p] //= p
    mu[val > 1] *= -1
    return mu

def primes_in(lo, hi):
    sieve = np.ones(hi - lo + 1, dtype=bool)
    if lo <= 1:
        sieve[:max(0, 2 - lo)] = False
    for p in range(2, int(hi ** 0.5) + 1):
        start = max(p * p, ((lo + p - 1) // p) * p)
        sieve[start - lo::p] = False
    return np.nonzero(sieve)[0] + lo

def small_primes(n):
    ps = []
    s = np.ones(n + 1, dtype=bool); s[:2] = False
    for i in range(2, n + 1):
        if s[i]:
            ps.append(i)
            s[i*i::i] = False
    return ps

SP = small_primes(QCUT)

def rho_q(q, N, k, kp):
    """Exact density over p in units mod q^2 of q^2 | (N-pk) or
    q^2 | (N-pk')."""
    q2 = q * q
    cnt = 0
    tot = 0
    for p in range(q2):
        if p % q == 0:
            continue
        tot += 1
        if (N - p * k) % q2 == 0 or (N - p * kp) % q2 == 0:
            cnt += 1
    return cnt / tot

def tail_factor():
    # prod over primes q in (QCUT, 100000] of (1 - 2/q^2)
    ps = small_primes(100000)
    f = 1.0
    for q in ps:
        if q > QCUT:
            f *= (1 - 2 / (q * q))
    return f

def main():
    rng = np.random.default_rng(20260811)
    N = 199_999_998
    t0 = time.time()
    mu = mobius_upto(N)
    print(f"mu ready {time.time()-t0:.0f}s", flush=True)
    tf = tail_factor()
    print(f"tail factor (q>{QCUT}) = {tf:.5f}", flush=True)
    K0, K1 = 2000, 4000
    P0 = N // (2 * K1); P1 = 2 * P0
    ps = primes_in(P0, P1)

    target = 4000
    s_pred = np.zeros(target); s_obs = np.zeros(target)
    r1_raw = np.zeros(target)
    done = 0
    while done < target:
        k = int(rng.integers(K0, K1)); kp = int(rng.integers(K0, K1))
        if k == kp:
            continue
        pmax = (N - 2) // max(k, kp)
        pp = ps[ps <= pmax]
        if len(pp) < 200:
            continue
        t = (mu[N - pp * k].astype(np.int64) *
             mu[N - pp * kp].astype(np.int64))
        pred = tf
        for q in SP:
            pred *= (1 - rho_q(q, N, k, kp))
        s_pred[done] = pred
        s_obs[done] = np.count_nonzero(t) / len(pp)
        r1_raw[done] = abs(t.sum()) / np.sqrt(len(pp))
        done += 1
        if done % 500 == 0:
            print(f"{done}/{target}  t={time.time()-t0:.0f}s", flush=True)

    corr = np.corrcoef(s_pred, s_obs)[0, 1]
    ratio = np.mean(s_obs / np.maximum(s_pred, 1e-12))
    maxerr = np.max(np.abs(s_obs - s_pred))
    print(f"(1) mask blind: corr(s_pred, s_obs) = {corr:.4f}  "
          f"mean ratio = {ratio:.4f}  max|err| = {maxerr:.4f}",
          flush=True)
    # bucket check
    for lo, hi in ((0, .1), (.1, .3), (.3, .5), (.5, .7), (.7, 1.01)):
        m = (s_pred >= lo) & (s_pred < hi)
        if m.sum() < 30:
            continue
        print(f"    s_pred [{lo:.1f},{hi:.1f}): n={m.sum():5d}  "
              f"pred={s_pred[m].mean():.3f}  obs={s_obs[m].mean():.3f}",
              flush=True)
    pred_r1 = 0.798 * np.sqrt(s_pred)
    m = s_pred > 0.01
    rr = np.mean(r1_raw[m]) / np.mean(pred_r1[m])
    print(f"(2) amplitude blind: mean r1_raw / mean 0.798*sqrt(s_pred) "
          f"= {rr:.4f}  (target 1.0)", flush=True)
    print("DONE", flush=True)

if __name__ == "__main__":
    main()
