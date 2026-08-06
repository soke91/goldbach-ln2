# -*- coding: utf-8 -*-
"""
Forge round-2 kill-test R1 (increment 178): is the zero spectrum
VISIBLE in the dilate field?

The explicit formula writes Mobius partial sums as sums of x^{rho}
oscillations. If that structure survives the pairing with mu(N-mk),
the field D(k) = Sum_m mu(m) mu(N-mk) has a component aligned with
the zero-templates
    T_gamma(k) = Sum_m m^{i gamma} mu(N-mk),
and the explicit formula becomes an external handle on that
component (a linearization channel for part of the energy).

Kill-test: least-squares projection of {D(k)} onto
span{Re T_gamma, Im T_gamma : first 30 zeros} over 300 k, vs the
same-dimension projection onto random-frequency templates
(6 independent draws of 30 frequencies uniform in [10, 105]).
Pre-registered: ALIVE iff R2_zeros >= 2 x mean(R2_random); DEAD
otherwise (the zero oscillations are invisible to the paired field
-- consistent with L's featurelessness, and closing the direct
explicit-formula channel).
"""
import numpy as np, time

ZEROS = [14.134725, 21.022040, 25.010858, 30.424876, 32.935062,
         37.586178, 40.918719, 43.327073, 48.005151, 49.773832,
         52.970321, 56.446248, 59.347044, 60.831779, 65.112544,
         67.079811, 69.546402, 72.067158, 75.704691, 77.144840,
         79.337375, 82.910381, 84.735493, 87.425275, 88.809111,
         92.491899, 94.651344, 95.870634, 98.831194, 101.317851]

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

def build(mu, N, ks, freqs, t0):
    """Return y (D values) and X (templates, len(ks) x 2*len(freqs))."""
    SQ = int(N ** 0.5)
    y = np.zeros(len(ks))
    X = np.zeros((len(ks), 2 * len(freqs)))
    for i, k in enumerate(ks):
        k = int(k)
        ms = np.arange(SQ + 1, N // k + 1, dtype=np.int64)
        w = mu[N - k * ms].astype(np.float64)
        y[i] = float((mu[ms] * mu[N - k * ms]).astype(np.int64).sum())
        lm = np.log(ms.astype(np.float64))
        for j, g in enumerate(freqs):
            ph = g * lm
            X[i, 2*j] = float(np.dot(np.cos(ph), w))
            X[i, 2*j+1] = float(np.dot(np.sin(ph), w))
        if i % 50 == 49:
            print(f"  k {i+1}/{len(ks)}  t={time.time()-t0:.0f}s",
                  flush=True)
    return y, X

def r2_of(y, X):
    coef, _, _, _ = np.linalg.lstsq(X, y, rcond=None)
    resid = y - X @ coef
    return 1 - float(np.sum(resid**2)) / float(np.sum((y - y.mean())**2))

def main():
    rng = np.random.default_rng(20260902)
    N = 199_999_998
    t0 = time.time()
    mu = mobius_upto(N)
    print(f"mu ready {time.time()-t0:.0f}s", flush=True)
    ks = np.arange(2000, 2300)  # 300 k

    print("building zero-templates...", flush=True)
    y, Xz = build(mu, N, ks, ZEROS, t0)
    R2z = r2_of(y, Xz)
    print(f"R2_zeros = {R2z:.4f}  (60 regressors / 300 points; "
          f"chance ~ 60/300 = 0.20)", flush=True)

    r2rs = []
    for d in range(6):
        freqs = rng.uniform(10, 105, size=len(ZEROS))
        print(f"random draw {d+1}/6...", flush=True)
        _, Xr = build(mu, N, ks, freqs, t0)
        r2rs.append(r2_of(y, Xr))
        print(f"  R2_random[{d}] = {r2rs[-1]:.4f}", flush=True)
    m = float(np.mean(r2rs)); s = float(np.std(r2rs))
    print(f"R2_random = {m:.4f} +- {s:.4f}", flush=True)
    print("=== R1 KILL-TEST ===", flush=True)
    print(f"ratio = {R2z/max(m,1e-9):.2f}  (pre-registered ALIVE "
          f"threshold: >= 2)", flush=True)
    print("verdict:",
          "ALIVE (zero spectrum visible)" if R2z >= 2*m else
          "DEAD (zero oscillations invisible to the paired field; "
          "the direct explicit-formula channel closes)", flush=True)
    print("DONE", flush=True)

if __name__ == "__main__":
    main()
