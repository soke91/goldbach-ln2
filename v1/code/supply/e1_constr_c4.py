# -*- coding: utf-8 -*-
"""
Construction experiment C4 (increment 184): manufactured modularity --
the last measurable candidate.

Does the generating function Phi(z) = Sum_m t_m e(mz) of the pair
field satisfy ANY approximate Fricke-type transformation law
    |Phi(-1/(Q^2 z))| ~ |Q z|^kappa |Phi(z)|
for some level Q and weight kappa? A discovered approximate law would
be the seed of a C-IV representation; its absence closes the last
candidate with a numeric probe (C-III, the Motohashi-type class,
requires an automorphic realization of 1/zeta -- a known-open
problem with no finite test).

Method: for k0 = 2001 (t_m on m in (sqrt N, N/k0], M ~ 8.6e4 terms),
evaluate Phi at 40 z-points near |z| ~ 1/Q (where both z and the
involuted point have workable decay), per Q in {1..6}. Defect(Q) =
median over the grid of | log(|Phi(wz)| / (|Qz|^kappa |Phi(z)|)) |
minimized over kappa. Compare against 6 mask-null draws.
Pre-registered: C-IV ALIVE iff defect_real <= 0.5 x defect_null at
>= 2 levels Q; else CLOSED.
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

def phi_vals(t, ms, zs):
    out = np.empty(len(zs), dtype=np.complex128)
    mf = ms.astype(np.float64)
    for i, z in enumerate(zs):
        out[i] = np.sum(t * np.exp(2j * np.pi * mf * z))
    return out

def defect(t, ms, Q, rng):
    # grid near |z| ~ 1/Q with y-scale set for effective decay
    M = float(ms[-1])
    y0 = 3.0 / M          # keeps ~M/3 effective terms
    zs, ws = [], []
    tries = 0
    while len(zs) < 40 and tries < 4000:
        tries += 1
        x = rng.uniform(0.2, 1.0) / Q
        y = y0 * rng.uniform(1.0, 4.0)
        z = complex(x, y)
        w = -1.0 / (Q * Q * z)
        if w.imag >= y0:      # involuted point must also converge
            zs.append(z); ws.append(w)
    if len(zs) < 15:
        return None
    P1 = np.abs(phi_vals(t, ms, zs))
    P2 = np.abs(phi_vals(t, ms, ws))
    az = np.abs(np.array(zs)) * Q
    good = (P1 > 1e-9) & (P2 > 1e-9)
    if good.sum() < 10:
        return None
    L = np.log(P2[good] / P1[good])
    X = np.log(az[good])
    # best kappa: minimize median |L - kappa X|; scan
    kaps = np.linspace(-6, 6, 241)
    meds = [np.median(np.abs(L - kp * X)) for kp in kaps]
    return float(np.min(meds))

def main():
    rng = np.random.default_rng(20260908)
    N = 199_999_998
    t0 = time.time()
    mu = mobius_upto(N)
    print(f"mu ready {time.time()-t0:.0f}s", flush=True)
    k0 = 2001
    SQ = int(N ** 0.5)
    ms = np.arange(SQ + 1, N // k0 + 1, dtype=np.int64)
    t = (mu[ms].astype(np.int64) * mu[N - k0 * ms]).astype(np.float64)
    sup = t != 0
    print(f"M = {len(ms)}", flush=True)

    alive_levels = 0
    for Q in range(1, 7):
        dr = defect(t, ms, Q, rng)
        dns = []
        for _ in range(6):
            s = rng.choice([-1.0, 1.0], size=t.shape)
            tn = np.where(sup, s, 0.0)
            d = defect(tn, ms, Q, rng)
            if d is not None:
                dns.append(d)
        if dr is None or not dns:
            print(f"Q={Q}: insufficient grid", flush=True)
            continue
        mn = float(np.mean(dns))
        ok = dr <= 0.5 * mn
        alive_levels += int(ok)
        print(f"Q={Q}: defect_real={dr:.3f}  defect_null={mn:.3f} "
              f"(x{dr/max(mn,1e-9):.2f})"
              f"{'  <-- LAW CANDIDATE' if ok else ''}  "
              f"t={time.time()-t0:.0f}s", flush=True)

    print(f"=== C4: {alive_levels} law-candidate levels ===", flush=True)
    print("verdict:",
          "C-IV ALIVE (approximate transformation law found)"
          if alive_levels >= 2 else
          "C-IV CLOSED (no approximate Fricke-type law at any tested "
          "level; the generating function is as modularity-free as "
          "its mask-null)", flush=True)
    print("DONE", flush=True)

if __name__ == "__main__":
    main()
