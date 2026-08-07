# -*- coding: utf-8 -*-
"""
C-III task T-S0b (increment 187) -- corrected sampling: measure Lemma S's object BEFORE any
theory -- does mu against Bessel/square-root phases actually exhibit
square-root cancellation, uniformly?

Object: S(u', c) = Sum_{w ~ W} mu(w) e(2 sqrt(u'(N-w)) / c).
Sample (u', c) over ranges with meaningful oscillation counts
(osc ~ W sqrt(u')/(c sqrt(N)) spanning ~3 to ~3000), 240 pairs.
Report r = |S|/sqrt(support) against half-normal 0.798; flag any
family-level failure (mean r >> 0.8 or heavy tail of large r).
Pre-registered: Lemma-S object HEALTHY iff family mean r in
[0.6, 1.0] and max r < 4; else DEAD ON ARRIVAL and the C-III tree
closes at T-S0.
"""
import numpy as np, time

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

def main():
    rng = np.random.default_rng(20260910)
    N = 199_999_998
    t0 = time.time()
    mu = mobius_upto(N)
    print(f"mu ready {time.time()-t0:.0f}s", flush=True)
    W0, W1 = 50_000_000, 100_000_000
    ws = np.arange(W0, W1, dtype=np.int64)
    muw = mu[ws].astype(np.float64)
    sup = int(np.count_nonzero(muw))
    root = np.sqrt((N - ws).astype(np.float64))
    print(f"W-range {W0}..{W1}, support {sup}", flush=True)

    rs = []
    worst = 0.0
    for i in range(240):
        # T-S0b fix: sample the oscillation count log-uniformly in
        # [3, 3000] and back-solve u' -- guarantees valid pairs
        # (the first run's ranges left only 9/240 usable; increment
        # 187 correction)
        c = float(np.exp(rng.uniform(np.log(2_000), np.log(300_000))))
        osc = float(np.exp(rng.uniform(np.log(3.0), np.log(3000.0))))
        up = (osc * c * np.sqrt(N) / (W1 - W0)) ** 2
        if up < 1:
            continue
        ph = (2.0 * np.sqrt(up) / c) * root
        S = np.sum(muw * np.exp(2j * np.pi * ph))
        r = abs(S) / np.sqrt(sup)
        rs.append((r, osc, up, c))
        worst = max(worst, r)
        if (i + 1) % 40 == 0:
            print(f"{i+1}/240  t={time.time()-t0:.0f}s", flush=True)

    rr = np.array([x[0] for x in rs])
    print(f"pairs used: {len(rr)}", flush=True)
    print(f"mean r = {rr.mean():.3f}  (half-normal 0.798)", flush=True)
    print(f"max r  = {rr.max():.3f}", flush=True)
    # oscillation-band breakdown
    oscs = np.array([x[1] for x in rs])
    for lo, hi in ((3, 30), (30, 300), (300, 1e9)):
        m = (oscs >= lo) & (oscs < hi)
        if m.sum() >= 20:
            print(f"  osc [{lo:.0f},{hi:.0f}): n={m.sum():4d}  "
                  f"mean r={rr[m].mean():.3f}  max={rr[m].max():.3f}",
                  flush=True)
    healthy = (0.6 <= rr.mean() <= 1.0) and (rr.max() < 4)
    print("=== T-S0 ===", flush=True)
    print("verdict:",
          "HEALTHY (Lemma S's object shows uniform square-root "
          "cancellation -- theory tasks T-S1..T-S3 proceed)"
          if healthy else
          "DEAD ON ARRIVAL (nature refuses; the C-III tree closes "
          "at T-S0)", flush=True)
    print("DONE", flush=True)

if __name__ == "__main__":
    main()
