# -*- coding: utf-8 -*-
"""
Is the location mask multiplicative, and with what factor?
(increment 243)

Increment 242 showed the mask compounds across primes rather than
adding. Increment 243's first attempt measured the per-prime factor by
splitting targeted families, and it was too noisy to settle anything:
the two cores disagreed by 2-3 standard errors (f(29) read 1.11 with
one core and 1.63 with the other), and -- recorded as a fault -- the
pre-registered control prime q = 101 was SILENTLY SKIPPED by an
n >= 12 guard, so the design's own noise was never calibrated. A
control that does not run and does not say so is worse than no control.

This uses all the data instead of a slice of it. The FFT sweep gives
Z(N) for every even N <= 4 x 10^6, band-standardised; the 256 cells of
the finite modular enumeration over q in {3,...,23} give cell means
m_c with counts n_c; and the two candidate laws are fitted to those
cell means directly:

    MULTIPLICATIVE   -m_c  =  A * prod_{q in c} f(q)
    ADDITIVE         m_c   =  a + sum_{q in c} b_q

The multiplicative fit is linear in logs, log(-m_c) = log A + sum
log f(q), weighted by n_c. Cells with too few members, or with m_c not
clearly negative, are excluded and the exclusion is reported.

PREDICTIONS, stated before fitting.
  * From the mechanism -- v = N - n is forced coprime to rad(N), and
    a Mobius sum restricted away from q gains an Euler factor -- the
    prediction is f(q) = q/(q-1): 1.500, 1.250, 1.167, 1.100, 1.083,
    1.063, 1.056, 1.045 for q = 3, 5, 7, 11, 13, 17, 19, 23.
  * The competing shape from the singular series is (q-1)/(q-2):
    2.000, 1.333, 1.200, 1.111, 1.091, 1.067, 1.059, 1.048. These
    differ sharply at q = 3 (1.50 against 2.00) and hardly at all
    above q = 11, so q = 3 is where the two can be told apart.
  * CONFIRMED iff the fitted f(q) track one of these within their
    errors AND the multiplicative fit beats the additive one on the
    same cells.

NULL FOR THE FIT. R^2 on the cell means, against an intercept-only
model. A control is built in and cannot be silently skipped: the same
fit is run on a PERMUTED assignment of cells to N, where by
construction f(q) = 1 for every q, and the fitted spread there is the
noise floor for every number in the real fit.
"""
import numpy as np
import math
import time


def sieve(X):
    spf = np.zeros(X + 1, dtype=np.int32)
    for i in range(2, int(X ** 0.5) + 1):
        if spf[i] == 0:
            sl = spf[i * i::i]; sl[sl == 0] = i
    for i in range(2, X + 1):
        if spf[i] == 0:
            spf[i] = i
    mu = np.zeros(X + 1, dtype=np.int8); mu[1] = 1
    for i in range(2, X + 1):
        p = int(spf[i]); j = i // p
        mu[i] = 0 if j % p == 0 else -mu[j]
    primes = np.nonzero(spf[2:] == np.arange(2, X + 1))[0] + 2
    lam = np.zeros(X + 1)
    for p in primes:
        q = int(p); lp = math.log(int(p))
        while q <= X:
            lam[q] = lp; q *= int(p)
    return mu, lam, primes


def cell_means(Zb, cell, ncell, minn):
    m = np.zeros(ncell); n = np.zeros(ncell, dtype=np.int64)
    s = np.zeros(ncell)
    for c in range(ncell):
        sel = cell == c
        n[c] = int(sel.sum())
        if n[c] >= 2:
            v = Zb[sel]
            m[c] = float(v.mean())
            s[c] = float(v.std(ddof=1) / math.sqrt(n[c]))
    return m, s, n


def fit_both(m, s, n, QS, minn, tag):
    # EVALUATE on every cell with enough members, with NO sign filter.
    # Filtering cells by mean < 0 selects negative fluctuations and
    # manufactures R^2 even on permuted labels -- the first version of
    # this script did exactly that and its "noise floor" of 0.79 was
    # the selection, not the noise.
    ncell = len(m)
    idx = np.nonzero(n >= minn)[0]
    if len(idx) < len(QS) + 2:
        print(f"  {tag}: only {len(idx)} usable cells, cannot fit")
        return None
    A = np.zeros((len(idx), 1 + len(QS)))
    A[:, 0] = 1.0
    for r, c in enumerate(idx):
        for i in range(len(QS)):
            A[r, 1 + i] = 1.0 if (c >> i) & 1 else 0.0
    w = np.sqrt(n[idx].astype(np.float64))
    y = m[idx]

    ba, *_ = np.linalg.lstsq(A * w[:, None], y * w, rcond=None)
    pred_add = A @ ba

    # multiplicative: fitted in logs on the negative cells only (logs
    # require it), then EVALUATED on all of them alongside the additive
    neg = y < -0.02
    if neg.sum() >= len(QS) + 2:
        bm, *_ = np.linalg.lstsq((A[neg] * w[neg, None]),
                                 np.log(-y[neg]) * w[neg], rcond=None)
        pred_mul = -np.exp(A @ bm)
    else:
        bm = np.zeros(1 + len(QS)); pred_mul = np.zeros(len(idx))

    def wr2(p):
        num = float((w * (y - p) ** 2).sum())
        den = float((w * (y - np.average(y, weights=w)) ** 2).sum())
        return 1 - num / den

    r2m, r2a = wr2(pred_mul), wr2(pred_add)
    print(f"  {tag}: {len(idx)} cells (n >= {minn}), no sign filter; "
          f"{int(neg.sum())} negative")
    print(f"    weighted R^2 on cell means:  multiplicative {r2m:+.4f}"
          f"   additive {r2a:+.4f}")
    return bm, ba, idx, r2m, r2a


def main():
    X = 4_000_000
    lo = 100_000
    t0 = time.time()
    mu, lam, primes = sieve(X)
    n_fft = 1
    while n_fft < 2 * (X + 1):
        n_fft *= 2
    F = np.zeros(n_fft); F[: X + 1] = mu
    G = np.zeros(n_fft); G[: X + 1] = lam
    C = np.fft.irfft(np.fft.rfft(F) * np.fft.rfft(G), n_fft)[: X + 1]
    F[: X + 1] = np.abs(mu); G[: X + 1] = lam ** 2
    V = np.fft.irfft(np.fft.rfft(F) * np.fft.rfft(G), n_fft)[: X + 1]
    del F, G
    print(f"convolutions t={time.time()-t0:.0f}s", flush=True)

    Ns = np.arange(lo, X + 1, 2)
    Z = C[Ns] / np.sqrt(V[Ns])
    Zb = np.empty_like(Z)
    b = lo
    while b < X:
        hi = min(2 * b, X)
        sel = (Ns >= b) & (Ns < hi)
        if sel.sum() > 1000:
            v = Z[sel]; Zb[sel] = (v - v.mean()) / v.std()
        else:
            Zb[sel] = 0.0
        b *= 2

    QS = [3, 5, 7, 11, 13, 17, 19, 23]
    ncell = 1 << len(QS)
    cell = np.zeros(len(Ns), dtype=np.int64)
    for i, q in enumerate(QS):
        cell |= ((Ns % q == 0).astype(np.int64) << i)

    m, s, n = cell_means(Zb, cell, ncell, 100)
    print(f"\n(A) fitting the two laws to the cell means")
    out = fit_both(m, s, n, QS, 100, "real")

    # the control that cannot be skipped: permute the cell labels, so
    # every f(q) is 1 by construction, and refit identically
    rng = np.random.default_rng(20260806)
    perm = rng.permutation(len(Ns))
    m0, s0, n0 = cell_means(Zb[perm], cell, ncell, 100)
    print(f"\n(B) the same fit on permuted labels -- the noise floor")
    out0 = fit_both(m0, s0, n0, QS, 100, "permuted")

    if out is None:
        print("DONE"); return
    bm, ba, idx, r2m, r2a = out
    print(f"\n(C) the fitted per-prime factors")
    print(f"{'q':>4} {'f(q) fitted':>12} {'q/(q-1)':>9} "
          f"{'(q-1)/(q-2)':>12} {'permuted f':>11}")
    for i, q in enumerate(QS):
        fq = math.exp(bm[1 + i])
        f0 = math.exp(out0[0][1 + i]) if out0 is not None else float('nan')
        print(f"{q:>4} {fq:>12.4f} {q/(q-1.0):>9.4f} "
              f"{(q-1.0)/(q-2.0):>12.4f} {f0:>11.4f}")
    print(f"  A = {math.exp(bm[0]):.4f}   (the base depth at rad(N) = 2)")

    print(f"\n(D) how well the multiplicative law predicts the deep cells")
    print(f"{'primes | N':>28} {'count':>7} {'measured':>9} "
          f"{'mult pred':>10} {'add pred':>9}")
    order = np.argsort(m)
    shown = 0
    for c in order:
        if n[c] < 100:
            continue
        lab = "*".join(str(q) for i, q in enumerate(QS) if c >> i & 1)
        row = np.zeros(1 + len(QS)); row[0] = 1.0
        for i in range(len(QS)):
            row[1 + i] = 1.0 if (c >> i) & 1 else 0.0
        pm = -math.exp(float(row @ bm))
        pa = float(row @ ba)
        print(f"{lab or '(none)':>28} {n[c]:>7} {m[c]:>9.3f} "
              f"{pm:>10.3f} {pa:>9.3f}")
        shown += 1
        if shown >= 12:
            break

    # Neither closed form works, so ask the model-free question: is the
    # cell mean a function of the cell's singular-series factor alone?
    # Sig(c) = prod_{q in c} (q-1)/(q-2) is what S(N) contributes from
    # the primes the enumeration sees.
    print(f"\n(E) is the cell mean a power of the singular-series factor?")
    sig = np.ones(ncell)
    for c in range(ncell):
        for i, q in enumerate(QS):
            if c >> i & 1:
                sig[c] *= (q - 1.0) / (q - 2.0)
    sel = (n >= 100) & (m < -0.02)
    x = np.log(sig[sel]); y = np.log(-m[sel])
    w = np.sqrt(n[sel].astype(np.float64))
    Amat = np.vstack([np.ones_like(x), x]).T
    co, *_ = np.linalg.lstsq(Amat * w[:, None], y * w, rcond=None)
    pred = Amat @ co
    r2 = 1 - float((w * (y - pred) ** 2).sum()) / \
        float((w * (y - np.average(y, weights=w)) ** 2).sum())
    print(f"  fit  log(-mean) = {co[0]:+.4f} + {co[1]:.4f} log(Sig)")
    print(f"  i.e. mean = -{math.exp(co[0]):.4f} * Sig^{co[1]:.3f}"
          f"    weighted R^2 = {r2:.4f}   ({int(sel.sum())} cells)")
    # the same on permuted labels, as the floor
    sel0 = (n0 >= 100) & (m0 < -0.02)
    if sel0.sum() > 4:
        x0 = np.log(sig[sel0]); y0 = np.log(-m0[sel0])
        w0 = np.sqrt(n0[sel0].astype(np.float64))
        A0 = np.vstack([np.ones_like(x0), x0]).T
        c0, *_ = np.linalg.lstsq(A0 * w0[:, None], y0 * w0, rcond=None)
        p0 = A0 @ c0
        r20 = 1 - float((w0 * (y0 - p0) ** 2).sum()) / \
            float((w0 * (y0 - np.average(y0, weights=w0)) ** 2).sum())
        print(f"  permuted floor: exponent {c0[1]:+.3f}, R^2 {r20:.4f}"
              f"   ({int(sel0.sum())} cells)")
    print(f"\n  {'primes | N':>26} {'Sig':>7} {'measured':>9} "
          f"{'power-law':>10}")
    shown = 0
    for c in np.argsort(m):
        if not sel[c]:
            continue
        lab = "*".join(str(q) for i, q in enumerate(QS) if c >> i & 1)
        pv = -math.exp(co[0]) * sig[c] ** co[1]
        print(f"  {lab:>26} {sig[c]:>7.3f} {m[c]:>9.3f} {pv:>10.3f}")
        shown += 1
        if shown >= 10:
            break
    print("DONE")


if __name__ == "__main__":
    main()
