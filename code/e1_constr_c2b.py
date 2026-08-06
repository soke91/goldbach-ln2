# -*- coding: utf-8 -*-
"""
Construction C2b (increment 183): adversarial verification of the C2
hits BEFORE any claim.

C2 flagged 4/300 k with inverse-domain special-frequency excess
(max z = +10.97). Known failure modes from this program's history:
imprecise null SE (8 draws), wrong-null category errors (the z=9
Wishart precedent), degenerate special frequencies, heavy-tailed
statistics. Verification battery on the 4 hit k's + 4 random
control k's:
  (1) 64-draw nulls (accurate SE) + fresh recomputation;
  (2) per-line breakdown: which frequency carries the excess
      (a = 1 / N / Nbar), and is N % k degenerate (near +-1, tiny,
      or sharing factors with k)?
  (3) replication at a second N (fresh arithmetic; the special
      frequencies move with N -- a real inverse-domain structure
      tied to N should move with them, an artifact of k should not);
  (4) permutation null (shuffle the pairing p <-> value instead of
      randomizing signs): kills any sign structure AND any
      value-support correlation while keeping both marginals.
Pre-registered: the hit SURVIVES iff z_64 >= 4 at the same k AND the
second-N replication shows the family-level effect (>= 1 of the 4
transported tests with z >= 3); else the C2 verdict downgrades.
"""
import numpy as np, time
from math import gcd

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

def lines(vals, invs, k, N):
    h = np.zeros(k)
    np.add.at(h, invs, vals)
    g = np.fft.rfft(h)
    E = np.abs(g) ** 2
    out = {}
    Nk = N % k
    Nb = pow(Nk, -1, k) if gcd(Nk, k) == 1 else None
    for name, a in (("a=1", 1), ("a=N", Nk),
                    ("a=Nbar", Nb if Nb is not None else 0)):
        aa = min(a % k, (-a) % k)
        out[name] = float(E[aa]) if aa < len(E) else 0.0
    out["total_special"] = sum(out.values())
    return out

def analyze(mu, N, k, ps, rng, ndraw, tag):
    pmax = (N - 2) // k
    pp = ps[(ps <= pmax) & (ps % k != 0)]
    vals = mu[N - pp * k].astype(np.float64)
    invs = np.array([pow(int(p), -1, k) for p in pp], dtype=np.int64)
    real = lines(vals, invs, k, N)
    sup = vals != 0
    nullsE = {key: [] for key in real}
    for _ in range(ndraw):
        s = rng.choice([-1.0, 1.0], size=vals.shape)
        vn = np.where(sup, s, 0.0)
        nl = lines(vn, invs, k, N)
        for key in nl:
            nullsE[key].append(nl[key])
    # permutation null: shuffle values across support positions
    permE = []
    idx = np.nonzero(sup)[0]
    for _ in range(ndraw):
        vperm = np.zeros_like(vals)
        vperm[idx] = vals[idx][rng.permutation(len(idx))]
        permE.append(lines(vperm, invs, k, N)["total_special"])
    rep = [f"[{tag}] k={k}  N%k={N%k}  gcd(N,k)={gcd(N, k)}"]
    for key in real:
        m = float(np.mean(nullsE[key])); s = float(np.std(nullsE[key]))
        z = (real[key] - m) / max(s, 1e-9)
        rep.append(f"    {key:14s} real={real[key]:.3e}  z={z:+.2f}")
    mp = float(np.mean(permE)); sp = float(np.std(permE))
    zp = (real["total_special"] - mp) / max(sp, 1e-9)
    rep.append(f"    perm-null      z={zp:+.2f}")
    print("\n".join(rep), flush=True)
    m = float(np.mean(nullsE["total_special"]))
    s = float(np.std(nullsE["total_special"]))
    return (real["total_special"] - m) / max(s, 1e-9), zp

def main():
    rng = np.random.default_rng(20260907)
    t0 = time.time()
    N1 = 199_999_998
    mu = mobius_upto(N1)
    print(f"mu(N1) ready {time.time()-t0:.0f}s", flush=True)
    P0 = N1 // (2 * 4200); P1 = 2 * P0
    ps = primes_in(P0, P1)

    HITS = [3859, 3275, 2827, 2993]
    CTRL = [2005, 2469, 3141, 3771]
    z64 = {}
    print("=== (1)(2)(4): 64-draw + per-line + permutation, N1 ===",
          flush=True)
    for k in HITS:
        z, zp = analyze(mu, N1, k, ps, rng, 64, "HIT")
        z64[k] = (z, zp)
    for k in CTRL:
        analyze(mu, N1, k, ps, rng, 64, "ctrl")
    del mu

    print("=== (3): replication at N2 ===", flush=True)
    N2 = 200_000_002  # fresh even N
    mu2 = mobius_upto(N2)
    print(f"mu(N2) ready {time.time()-t0:.0f}s", flush=True)
    ps2 = primes_in(N2 // (2*4200), 2*(N2 // (2*4200)))
    rep_hits = 0
    for k in HITS:
        z, zp = analyze(mu2, N2, k, ps2, rng, 64, "N2")
        if z >= 3:
            rep_hits += 1
    print("=== C2b SUMMARY ===", flush=True)
    for k in HITS:
        print(f"  k={k}: z64(N1)={z64[k][0]:+.2f}  perm={z64[k][1]:+.2f}",
              flush=True)
    print(f"  N2 replication hits (z>=3): {rep_hits}/4", flush=True)
    survive = any(z64[k][0] >= 4 for k in HITS) and rep_hits >= 1
    print("verdict:",
          "SURVIVES verification (genuine inverse-domain structure)"
          if survive else
          "DOWNGRADED (fails accurate-null or replication -- treat C2 "
          "as a statistical artifact pending deeper study)", flush=True)
    print("DONE", flush=True)

if __name__ == "__main__":
    main()
