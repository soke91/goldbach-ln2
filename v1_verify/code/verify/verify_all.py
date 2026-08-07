# -*- coding: utf-8 -*-
"""
One-shot reproduction stamp for v1_verify.

Every finding in `paper/ADVERSARIAL_FINDINGS.md` rests on a number, and
this re-derives each of those numbers at a size small enough to run in
a few minutes, judges it against a pre-registered interval, and exits
nonzero if any row fails. The full-size runs are the individual
`audit_*.py` scripts; this is the gate.

Each interval below was fixed from the corresponding full run's own
spread, not chosen to pass. Where a row's quantity is exact arithmetic
(an identity), the interval is machine precision.

    python v1_verify/code/verify/verify_all.py
"""
import sys
import math

import numpy as np

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass


def sieve(X):
    mu = np.ones(X + 1, dtype=np.int64)
    is_p = np.zeros(X + 1, dtype=bool)
    rem = np.arange(X + 1, dtype=np.int64)
    mu[0] = 0
    for p in range(2, X + 1):
        if rem[p] == p:
            is_p[p] = True
            mu[p::p] *= -1
            rem[p::p] //= p
            if p * p <= X:
                mu[p * p::p * p] = 0
    lam = np.zeros(X + 1, dtype=np.float64)
    for p in range(2, X + 1):
        if is_p[p]:
            lg = math.log(p)
            q = p
            while q <= X:
                lam[q] = lg
                q *= p
    return mu.astype(np.float64), lam


def conv(a, b, X):
    nf = 1
    while nf < 2 * (X + 2):
        nf *= 2
    return np.fft.irfft(np.fft.rfft(np.pad(a, (0, nf - len(a))))
                        * np.fft.rfft(np.pad(b, (0, nf - len(b)))),
                        nf)[: X + 1]


def autocorr_box(a, X):
    nf = 1
    while nf < 2 * (X + 2):
        nf *= 2
    F = np.fft.rfft(np.pad(a[: X + 1], (0, nf - X - 1)))
    return np.fft.irfft(F * np.conj(F), nf)[:X]


ROWS = []


def row(name, value, lo, hi, note=""):
    ROWS.append((name, value, lo, hi, note))


def main():
    X = 1_000_000
    mu, lam = sieve(X)
    mu2 = (mu != 0).astype(np.float64)

    # ---- 1. Lemma 13 as stated is false; the repaired form is exact --
    for XX in (2000, 8000):
        m, l = mu[: XX + 1], lam[: XX + 1]
        C = np.zeros(XX + 1)
        for n in range(2, XX + 1):
            if l[n]:
                C[n + 1:XX + 1] += l[n] * m[1:XX - n + 1]
        lhs_paper = float((C ** 2).sum())
        Ct = np.zeros(2 * XX + 1)
        for n in range(2, XX + 1):
            if l[n]:
                Ct[n + 1:n + XX + 1] += l[n] * m[1:XX + 1]
        lhs_rep = float((Ct ** 2).sum())
        M = autocorr_box(m, XX)
        P = autocorr_box(l, XX)
        rhs = float(M[0] * P[0] + 2.0 * np.dot(M[1:], P[1:]))
        row(f"Lemma 13 as stated, LHS/RHS at X={XX}",
            lhs_paper / rhs, 0.55, 0.70, "must not be 1")
        row(f"Lemma 13 repaired, |LHS/RHS - 1| at X={XX}",
            abs(lhs_rep / rhs - 1.0), 0.0, 1e-12, "exact")

    # ---- 2. Proposition 11: V = A(N) W ------------------------------
    V = conv(mu2, lam ** 2, X)
    W = float((lam[:X] ** 2).sum())
    Vx = float(V[X])
    A = 1.0
    for q in range(2, 100000):
        if all(q % r for r in range(2, int(q ** 0.5) + 1)):
            if X % q:
                A *= 1.0 - 1.0 / (q * (q - 1.0))
    row("Proposition 11: W(X)/V(X) against 1/A(N)",
        (W / Vx) / (1.0 / A), 0.995, 1.005, f"1/A = {1/A:.5f}")

    # ---- 3. the Gamma of Proposition 15 -----------------------------
    for N in (100_000, 1_000_000):
        lg = lam[:N]
        csum = float(lg.sum()) ** 2 - float((lg * lg).sum())
        v = np.arange(1, N)
        VN = float(((mu[1:N] ** 2) * (lam[N - v] ** 2)).sum())
        G = csum / VN
        row(f"Gamma(N) log N / N at N={N}", G * math.log(N) / N,
            1.30, 1.50, "-> 1/A = 1.270")
    row("Gamma growth over a factor 10 in N",
        (ROWS[-1][1] * 1_000_000 / math.log(1_000_000))
        / (ROWS[-2][1] * 100_000 / math.log(100_000)),
        7.0, 11.0, "must diverge, not be O(1)")

    # ---- 4. Proposition 20's three terms cancel ---------------------
    band = np.arange(X // 2 + (X // 2) % 2, X + 1, 2)
    band = band[V[band] > 0]
    n = len(band)
    depth = np.zeros(X + 1, dtype=np.int8)
    for p in (3, 5, 7, 11, 13):
        depth[::p] += 1
    nf = 1
    while nf < 2 * (X + 2):
        nf *= 2
    FL = np.fft.rfft(np.pad(lam, (0, nf - X - 1)))

    def u_of(sel):
        g = np.zeros(nf)
        g[sel] = 1.0 / np.sqrt(V[sel])
        return np.fft.irfft(np.fft.rfft(g) * np.conj(FL), nf)[: X + 1]

    ua = u_of(band)
    Q = lambda a, b: float((mu2[: X + 1] * a * b).sum())
    Qaa = Q(ua, ua)
    cell = band[depth[band] == 1]
    nc = len(cell)
    uc = u_of(cell)
    t1 = Q(uc, uc) / nc ** 2
    var = t1 - 2.0 * Q(uc, ua) / (nc * n) + Qaa / n ** 2
    row("Proposition 20: Var / (Q_cc/n_c^2) at depth 1",
        var / t1, 0.010, 0.030, "the three terms cancel")

    # ---- 5. the coin control on rho ---------------------------------
    C = conv(mu, lam, X)
    Ns = np.arange(200_002, X + 1, 2)
    key = np.zeros(len(Ns), dtype=np.int32)
    for i, q in enumerate((3, 5, 7, 11, 13)):
        key |= ((Ns % q) == 0).astype(np.int32) << i
    uq, inv = np.unique(key, return_inverse=True)
    cnt = np.bincount(inv, minlength=len(uq)).astype(np.float64)

    def rho_of(Cf):
        cv = Cf[Ns]
        tot = np.bincount(inv, weights=cv, minlength=len(uq))
        cc = cv - (tot / cnt)[inv]
        return float((cc ** 2).mean() / V[Ns].mean())

    rng = np.random.default_rng(5)
    idx = np.nonzero(mu2)[0]
    coins = []
    for _ in range(12):
        eps = np.zeros(X + 1)
        eps[idx] = rng.choice([-1.0, 1.0], size=len(idx))
        coins.append(rho_of(conv(eps, lam, X)))
    coins = np.array(coins)
    z = (rho_of(C) - coins.mean()) / coins.std()
    row("Lemma 17: |z| of the real rho against a coin", abs(z),
        0.0, 3.0, "a coin reproduces it")

    # ---- 6. the zeta null: white vs autocorrelation-preserving ------
    cv = C[Ns]
    tot = np.bincount(inv, weights=cv, minlength=len(uq))
    Z = (cv - (tot / cnt)[inv]) / np.sqrt(V[Ns])
    b = 200_000
    while b < X:
        hi = min(2 * b, X)
        s = (Ns >= b) & (Ns < hi)
        if int(s.sum()) > 1000:
            Z[s] = (Z[s] - Z[s].mean()) / Z[s].std(ddof=1)
        b = hi
    Z = Z - Z.mean()
    L = np.log(Ns.astype(np.float64))
    G10 = [14.134725, 21.022040, 25.010858, 30.424876, 32.935062,
           37.586178, 40.918719, 43.327073, 48.005151, 49.773832]
    cols = []
    for t in G10:
        cols.append(np.cos(t * L))
        cols.append(np.sin(t * L))
    B = np.stack(cols, axis=1)
    B -= B.mean(axis=0)
    Gi = np.linalg.inv(B.T @ B)

    def r2(y):
        yc = y - y.mean()
        bty = B.T @ yc
        return float((Gi @ bty) @ bty / (yc @ yc))

    nn = len(Z)
    perm = np.mean([r2(Z[rng.permutation(nn)]) for _ in range(30)])
    rot = np.mean([r2(np.roll(Z, int(rng.integers(1, nn))))
                   for _ in range(30)])
    row("conj:wall item 4: rotation null / permutation null",
        rot / perm, 20.0, 400.0, "the white null is wrong by ~10^2")

    # ---- 7. R4's exact identity -------------------------------------
    bad = 0
    for N in (5000, 20000):
        m = mu[: N + 1]
        tot = 0
        for k in range(1, N):
            j = np.arange(1, (N - 1) // k + 1, dtype=np.int64)
            if not len(j):
                break
            tot += int(np.dot(m[j].astype(np.int64),
                              m[N - j * k].astype(np.int64)))
        bad += (tot != int(m[N - 1]))
    row("R4: the full-range switch identity", float(bad), 0.0, 0.0,
        "= mu(N-1) exactly")

    # ---- report -----------------------------------------------------
    print("v1_verify reproduction stamp")
    print(f"{'row':<52} {'value':>12} {'interval':>20} {'':>6}")
    print("-" * 94)
    nfail = 0
    for name, val, lo, hi, note in ROWS:
        ok = lo <= val <= hi
        nfail += not ok
        print(f"{name:<52} {val:>12.6f} "
              f"[{lo:>8.4f}, {hi:>8.4f}] {'ok' if ok else 'FAIL':>6}"
              f"   {note}")
    print()
    print(f"{len(ROWS)} rows, {nfail} failed")
    if nfail:
        print("DONE (failed)")
        sys.exit(1)
    print("DONE")


if __name__ == "__main__":
    main()
