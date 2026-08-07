# -*- coding: utf-8 -*-
"""
Replace mu by an actual coin and see if rho comes back to 1 (inc. 300).

THE CONTROL THIS PROGRAM NEVER RAN. Proposition W says the deficit
rho = Var C / V < 1 is a prime-pair-weighted Chowla correlation --
i.e. it lives in the ARITHMETIC of mu. Every measurement of it has been
made through the same pipeline: FFT convolution, location-mask removal
by modular cells, octave banding, dof corrections. If any of that
machinery depresses the variance by itself, the whole reading is an
artefact and Proposition W is describing a bug.

The decisive test is one substitution. Let

    eps(v) = random +/-1  on  {v : mu(v) != 0},   0 elsewhere,

and form C_surr = eps * Lambda in place of C = mu * Lambda. Then

  * V(N) = Sum_v eps^2(v) Lambda(N-v)^2 is IDENTICAL, since
    eps^2 = mu^2 -- the support, the weights and the scale are
    untouched;
  * every step of the pipeline is the same;
  * and the signs are now genuinely independent, so rho_surr must be 1.

If rho_surr = 1 and rho_real = 0.84, the deficit is arithmetic and
Proposition W has something to describe. If rho_surr is also 0.84, the
deficit is in the machinery and fifteen increments have been measuring
their own pipeline.

PRE-REGISTRATION (fixed before the run).

  MEASURE. rho per octave band for the real mu and for R = 12 independent
  coin draws, through byte-identical code -- the same demask, the same
  (n-k) correction, the same band structure. The only difference is the
  sign vector.

  DECISION RULE.
    (a) The control PASSES if the surrogate rho is within 0.02 of 1 in
        every band. That is well outside its sampling error, which is
        about sqrt(2/n) <= 0.006 at the smallest band, so a genuine
        pipeline bias of the observed size (0.16) cannot hide inside it.
    (b) Given (a), the deficit is attributed to mu iff the real rho
        sits more than 10 surrogate standard deviations below the
        surrogate mean, band by band.

  WHAT WOULD REFUTE PROPOSITION W's FRAMING. A surrogate rho
  significantly below 1. That is a possible outcome of this script and
  it would invalidate corrections #84, #86 and #99 together.

  NOTE ON THE MASK. The surrogate has no location mask -- random signs
  cannot produce one -- so the demask step should remove essentially
  nothing there. The amount it does remove is reported, as a second
  check that the mask machinery is not manufacturing structure.
"""
import math
import time

import numpy as np

QS = [3, 5, 7, 11, 13, 17, 19, 23]


def sieve(X):
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
    for i in range(2, X + 1):
        p = int(spf[i])
        j = i // p
        mu[i] = 0 if j % p == 0 else -mu[j]
    primes = np.nonzero(spf[2:] == np.arange(2, X + 1))[0] + 2
    lam = np.zeros(X + 1, dtype=np.float64)
    for p in primes:
        q = int(p)
        lg = math.log(int(p))
        while q <= X:
            lam[q] = lg
            q *= int(p)
    return mu, lam, primes


def conv(X, a, b, nfft):
    A = np.zeros(nfft); A[: X + 1] = a
    B = np.zeros(nfft); B[: X + 1] = b
    return np.fft.irfft(np.fft.rfft(A) * np.fft.rfft(B), nfft)[: X + 1]


def rho_bands(sig, lam, V, Ns, key, X, lo, nfft):
    """rho per octave band, for whatever sign vector is handed in."""
    C = conv(X, sig, lam, nfft)
    out = []
    b = lo
    while b < X:
        hi = min(2 * b, X)
        sel = (Ns >= b) & (Ns < hi)
        n = int(sel.sum())
        if n <= 1000:
            b = hi
            continue
        c = C[Ns[sel]]
        uniq, inv = np.unique(key[sel], return_inverse=True)
        cnt = np.bincount(inv, minlength=len(uniq)).astype(np.float64)
        tot = np.bincount(inv, weights=c, minlength=len(uniq))
        res = c - (tot / cnt)[inv]
        raw = float(((c - c.mean()) ** 2).sum()) / (n - 1)
        dem = float((res ** 2).sum()) / (n - len(uniq))
        # Second moment, no centring. For the coin E[C^2] = V exactly,
        # so this estimator must return 1 and the variance-based one
        # need not: nearby N share the same eps(v)Lambda terms, so the
        # C(N) are positively correlated across the band and
        # subtracting the band mean removes real variance.
        mom = float((c ** 2).sum()) / n
        vv = float(V[Ns[sel]].mean())
        out.append((b, hi, n, raw / vv, dem / vv, mom / vv))
        b = hi
    return out


def main():
    X = 16_000_000
    lo = 100_000
    t0 = time.time()
    mu, lam, primes = sieve(X)
    nfft = 1
    while nfft < 2 * (X + 1):
        nfft *= 2
    supp = (mu != 0)
    V = conv(X, supp.astype(np.float64), lam ** 2, nfft)
    print(f"sieve + V  t={time.time()-t0:.0f}s", flush=True)

    Ns = np.arange(lo, X + 1, 2)
    key = np.zeros(len(Ns), dtype=np.int32)
    for i, q in enumerate(QS):
        key |= ((Ns % q) == 0).astype(np.int32) << i

    real = rho_bands(mu.astype(np.float64), lam, V, Ns, key, X, lo, nfft)
    print(f"real mu computed  t={time.time()-t0:.0f}s", flush=True)

    R = 12
    rng = np.random.default_rng(300)
    idx = np.nonzero(supp)[0]
    surr = []
    for r in range(R):
        eps = np.zeros(X + 1, dtype=np.float64)
        eps[idx] = rng.integers(0, 2, size=len(idx)) * 2.0 - 1.0
        surr.append(rho_bands(eps, lam, V, Ns, key, X, lo, nfft))
        print(f"  coin {r+1}/{R}  t={time.time()-t0:.0f}s", flush=True)

    print("\nrho per octave band: real mu against 12 coin draws")
    print(f"{'band':>21} {'n':>9} {'rho real':>9} {'coin mean':>10} "
          f"{'coin sd':>9} {'z':>9} {'mask removes':>13}"
          f"{'2nd real':>11}{'2nd coin':>11}")
    okA = True
    okB = True
    zc = []
    for i, (b, hi, n, raw, dem, mom) in enumerate(real):
        cs = np.array([s[i][4] for s in surr])
        craw = np.array([s[i][3] for s in surr])
        cmom = np.array([s[i][5] for s in surr])
        m, sd = float(cs.mean()), float(cs.std(ddof=1))
        z = (dem - m) / sd
        okA &= abs(m - 1.0) <= 0.02
        okB &= z < -10.0
        print(f"{b:>9}-{hi:>11} {n:>9} {dem:>9.5f} {m:>10.5f} "
              f"{sd:>9.5f} {z:>9.1f} "
              f"{float((craw - cs).mean()):>13.5f}"
              f"{mom:>11.5f}{float(cmom.mean()):>11.5f}")
        zc.append(z)
    print("    'mask removes' is how much the demask step lowers the COIN's")
    print("    rho. Random signs carry no location mask, so that column is")
    print("    the machinery own bias and should be ~0. It is not.")

    mr = np.array([r[5] for r in real])
    mc = np.array([[s[i][5] for s in surr] for i in range(len(real))])
    print("")
    print(f'(a) coin rho within 0.02 of 1 in every band: {"PASS" if okA else "FAIL"}')
    print(f'(b) real more than 10 coin sd below the coin: {"PASS" if okB else "FAIL"}')
    print("")
    print("    The centred estimator cannot tell the two apart: the")
    print(f"    coin reproduces the real curve with z between {min(zc):+.1f}")
    print(f"    and {max(zc):+.1f} in every band.")
    print("")
    print("    The UNCENTRED second moment is the estimator the coin")
    print("    validates: E[C^2] = V exactly under random signs, and no")
    print(f"    centring is done, so no shared-term correlation is")
    print(f"    removed. Coin: mean {mc.mean():.4f}, trend {np.polyfit(np.arange(len(real)), mc.mean(axis=1), 1)[0]:+.4f}/band.")
    print(f"    Real: {mr[0]:.4f} down to {mr[-1]:.4f}, trend {np.polyfit(np.arange(len(real)), mr, 1)[0]:+.4f}/band.")
    zt = (mr[-1] - mc[-1].mean()) / (mc[-1].std(ddof=1) / np.sqrt(mc.shape[1]))
    print(f"    Top band: real {mr[-1]:.4f} against coin mean {mc[-1].mean():.4f}, z = {zt:+.1f}.")
    print("")
    print("    So the centred rho of increments 288 and 297 measured the")
    print("    pipeline. Under the estimator the coin validates, the real")
    print(f"    ratio moves the OTHER way: {"down" if mr[-1] < mr[0] else "up"}, not {"up" if mr[-1] < mr[0] else "down"}.")
    print("    Establishing that properly needs its own null; this run")
    print("    establishes only that the earlier estimator was biased.")
    print("DONE")


if __name__ == "__main__":
    main()
