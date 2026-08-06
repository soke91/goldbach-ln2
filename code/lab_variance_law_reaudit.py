# -*- coding: utf-8 -*-
"""
The variance law, re-fitted with the mask removed, and asked whether
its exponent is identified at all (increment 280).

TWO FAULTS, FOUND BY LOOKING AT DATES.

(1) ORDER OF DISCOVERY. The law `Var C = kappa*S(N)*N*(log N)^alpha`
was measured at increments 236-238 (commit 014d491). The LOCATION MASK
-- the deterministic term m(N) indexed by which small primes divide N
-- was found at increment 240 (commit fa05ea9), two increments LATER.
The fit computes `G = C/sqrt(S*N)` and then `g.std()`, which removes
one band-wide mean and nothing else. So m(N)'s variation across N sits
INSIDE the measured variance, and the law has never been re-fitted
since the mask was known. If m(N) scales differently from
sqrt(S*N*log N), its contamination TILTS alpha.

(2) THE VERDICT OVERSTATES THE FIT. The recorded output declares
"CONFIRMED -- Var C(N) = 0.465 * S(N) * N * log N" while its own fitted
exponent is alpha = 0.8953. The pre-registered criterion was "fitted
alpha within 0.15 of 1", which 0.895 satisfies, so the criterion was
not violated -- but the sentence that got copied into CONJECTURE_L.md
and MEASUREMENTS.md asserts alpha = 1 exactly. Whether the data can
tell 0.895 from 1 was never asked.

PRE-REGISTRATION (fixed before the run).

  A. Mask contamination.
     H0: removing the mask leaves the exponent alone, |d alpha| < 0.02.
     H1: |d alpha| >= 0.02, i.e. the recorded exponent is contaminated.
     The mask is removed by the same finite modular enumeration that
     defined it: cells keyed by the subset of {3,5,7,11,13,17,19,23}
     dividing N, per-cell means subtracted WITHIN each band, with an
     exact (n - k) degrees-of-freedom correction so that removing k
     cell means cannot by itself shrink the variance.

  B. Is alpha identified?
     The noise floor is not a guess: for a variance estimated from n
     samples the sampling CV is sqrt(2/(n-1)), which each band supplies
     from its own count. alpha is declared IDENTIFIED only if the
     spread of alpha implied by that floor excludes the rival value.
     Reported as a standard error from weighted least squares, not as
     an eyeball.

  C. Is the LOG identified, or would a pure power do?
     Fit `sd^2 ~ (log N)^alpha` and `sd^2 ~ N^eps` on the same points.
     Both have two parameters, so their residuals are comparable. If
     N^eps fits as well, then over this range the data does not
     distinguish a log factor from a slightly larger power of N, and
     the law's FORM is consistent with the data rather than shown by
     it. Decision rule fixed now: the log is declared identified only
     if the power fit's residual RMS exceeds the log fit's by 2x.

WHY THIS MATTERS. `C(N) = m(N) + sqrt(kappa*S(N)*N*log N)*G(N)` is
Conjecture L's central formula and the campaign quotes kappa = 0.465
and the bare `log N` as though both were measured. This asks how much
of that is measurement.
"""
import numpy as np
import math
import sys
import time

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


def wls(x, y, w):
    """weighted least squares y = a + b x; returns b and its s.e."""
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    w = np.asarray(w, float)
    W = w.sum()
    mx = (w * x).sum() / W
    my = (w * y).sum() / W
    sxx = (w * (x - mx) ** 2).sum()
    sxy = (w * (x - mx) * (y - my)).sum()
    b = sxy / sxx
    a = my - b * mx
    resid = y - (a + b * x)
    dof = max(len(x) - 2, 1)
    s2 = (w * resid ** 2).sum() / dof
    se = math.sqrt(s2 / sxx)
    rms = math.sqrt((resid ** 2).mean())
    return a, b, se, rms


def main():
    X = int(sys.argv[1]) if len(sys.argv) > 1 else 4_000_000
    t0 = time.time()
    mu, lam, primes = sieve(X)

    n_fft = 1
    while n_fft < 2 * (X + 1):
        n_fft *= 2
    A = np.zeros(n_fft); A[: X + 1] = mu
    B = np.zeros(n_fft); B[: X + 1] = lam
    C = np.fft.irfft(np.fft.rfft(A) * np.fft.rfft(B), n_fft)[: X + 1]
    del A, B
    print(f"convolution  t={time.time()-t0:.0f}s", flush=True)

    C2 = 0.6601618158468696
    S = np.full(X + 1, 2 * C2)
    for p in primes:
        p = int(p)
        if p > 2:
            S[p::p] *= (p - 1) / (p - 2)

    lo = 100_000
    Ns = np.arange(lo, X + 1, 2)
    G = C[Ns] / np.sqrt(S[Ns] * Ns)

    # cell key: which of QS divide N -- the mask's own enumeration
    key = np.zeros(len(Ns), dtype=np.int32)
    for i, q in enumerate(QS):
        key |= ((Ns % q) == 0).astype(np.int32) << i

    print(f"\n{'band':>21} {'count':>8} {'cells':>6} {'sd^2 raw':>9} "
          f"{'sd^2 demask':>12} {'mask share':>11} {'logNmid':>8}")
    rows = []
    b = lo
    while b < X:
        hi = min(2 * b, X)
        sel = (Ns >= b) & (Ns < hi)
        n = int(sel.sum())
        if n <= 1000:
            b = hi
            continue
        g = G[sel]
        k = key[sel]
        v_raw = float(((g - g.mean()) ** 2).sum()) / (n - 1)

        # per-cell means within the band, then an exact dof correction
        uniq, inv = np.unique(k, return_inverse=True)
        cnt = np.bincount(inv, minlength=len(uniq))
        tot = np.bincount(inv, weights=g, minlength=len(uniq))
        cm = tot / cnt
        resid = g - cm[inv]
        kcells = len(uniq)
        v_dem = float((resid ** 2).sum()) / (n - kcells)

        mid = math.sqrt(b * hi)
        L = math.log(mid)
        share = 1.0 - v_dem / v_raw
        rows.append((n, kcells, v_raw, v_dem, L))
        print(f"{b:>9}-{hi:>11} {n:>8} {kcells:>6} {v_raw:>9.4f} "
              f"{v_dem:>12.4f} {share:>10.2%} {L:>8.3f}")
        b = hi

    ns = np.array([r[0] for r in rows], float)
    vr = np.array([r[2] for r in rows])
    vd = np.array([r[3] for r in rows])
    Ls = np.array([r[4] for r in rows])
    # sampling CV of a variance from n samples: sqrt(2/(n-1)).
    cvs = np.sqrt(2.0 / (ns - 1.0))
    w = 1.0 / cvs ** 2          # weights in log space

    print(f"\nsampling noise floor per band (CV of a variance):")
    print("  " + "  ".join(f"{c:.3%}" for c in cvs))

    print("\n(A) mask contamination: does removing it move alpha?")
    out = {}
    for lab, v in (("raw   ", vr), ("demask", vd)):
        _, al, se, rms = wls(np.log(Ls), np.log(v), w)
        out[lab.strip()] = (al, se, rms)
        print(f"  {lab}  alpha = {al:.4f} +/- {se:.4f}   "
              f"resid RMS(log) = {rms:.5f}")
    da = abs(out["raw"][0] - out["demask"][0])
    print(f"  |d alpha| = {da:.4f}   "
          f"{'H1: contaminated' if da >= 0.02 else 'H0: alpha unmoved'}")
    print(f"  mask share of the variance: "
          f"{1 - vd.mean()/vr.mean():.2%} on average")

    print("\n(B) is alpha identified? test the recorded claim alpha = 1")
    for lab in ("raw", "demask"):
        al, se, _ = out[lab]
        z = (al - 1.0) / se
        print(f"  {lab:>6}  alpha = {al:.4f} +/- {se:.4f}   "
              f"z vs alpha=1 : {z:+.2f}   "
              f"{'alpha=1 EXCLUDED' if abs(z) > 3 else 'alpha=1 consistent'}")

    print("\n(C) is the LOG identified, or would a pure power of N do?")
    print("    decision rule fixed before the run: the log is identified")
    print("    only if the power fit's residual RMS exceeds the log")
    print("    fit's by a factor of 2")
    lognmid = Ls           # log N
    for lab, v in (("raw   ", vr), ("demask", vd)):
        _, al, _, rms_log = wls(np.log(Ls), np.log(v), w)
        _, ep, _, rms_pow = wls(lognmid, np.log(v), w)
        ratio = rms_pow / rms_log if rms_log > 0 else float('inf')
        verdict = ("log identified" if ratio > 2.0
                   else "NOT identified -- a power of N fits as well")
        print(f"  {lab}  (log N)^{al:.3f}: RMS {rms_log:.5f}   "
              f"N^{ep:.5f}: RMS {rms_pow:.5f}   "
              f"ratio {ratio:.2f}   {verdict}")

    # (C2) range stability. A coefficient that is real should not walk
    # as the fitting window grows. This is hazard 5 turned on the
    # estimate itself: look for a trend in the ESTIMATE before quoting
    # the estimate.
    print("\n(C2) range stability: alpha from the first j bands only")
    print(f"{'j':>3} {'logN max':>9} {'alpha raw':>20} "
          f"{'alpha demask':>22}")
    for j in range(3, len(Ls) + 1):
        _, ar, ser, _ = wls(np.log(Ls[:j]), np.log(vr[:j]), w[:j])
        _, ad, sed, _ = wls(np.log(Ls[:j]), np.log(vd[:j]), w[:j])
        print(f"{j:>3} {Ls[j-1]:>9.3f} {ar:>12.4f} +/- {ser:<6.4f} "
              f"{ad:>13.4f} +/- {sed:<6.4f}")
    print("    A value that walks with j is not a measured exponent,")
    print("    whatever its standard error says.")

    # (C3) The mask's OWN scaling law, which nobody has measured. The
    # quantity removed above is Var(m)/(S*N) in Z units; if it decays
    # like N^-g then m(N) ~ sqrt(S) * N^{(1-g)/2}, and the mask is
    # o(sqrt(S*N)) exactly when g > 0. That decides whether the
    # deterministic term threatens C(N) = o(N) or is lower order.
    print("\n(C3) the mask's own scaling: Var_mask in Z units vs N")
    vm = vr - vd
    logN_band = Ls
    _, gm, seg, rmsg = wls(logN_band, np.log(vm), w)
    print(f"{'band logN':>10} {'Var_mask(Z)':>13}")
    for L, v in zip(logN_band, vm):
        print(f"{L:>10.3f} {v:>13.5f}")
    print(f"    fit Var_mask ~ N^g :  g = {gm:+.4f} +/- {seg:.4f}"
          f"   RMS {rmsg:.5f}")
    print(f"    => m(N) ~ sqrt(S(N)) * N^{(1 + gm) / 2:.4f}")
    print(f"    g < 0 means the mask DECAYS relative to the")
    print(f"    fluctuation, i.e. it is lower order and does not")
    print(f"    threaten C(N) = o(N).")
    print("    stability of g across the window:")
    for j in range(3, len(Ls) + 1):
        _, gj, sej, _ = wls(logN_band[:j], np.log(vm[:j]), w[:j])
        print(f"      j={j}  logN<={Ls[j-1]:6.3f}   "
              f"g = {gj:+.4f} +/- {sej:.4f}")

    print("\n(D) how much lever arm is there, really?")
    print(f"    N spans a factor {X/lo:.0f}, but log N spans only")
    print(f"    {Ls[-1]/Ls[0]:.4f} -- and it is log N that carries the")
    print(f"    exponent. log(log N) range = {math.log(Ls[-1]/Ls[0]):.4f}.")
    need = out["demask"][1] * math.sqrt(2.0)
    print(f"    To halve the s.e. on alpha ({out['demask'][1]:.4f}) the")
    print(f"    log(log N) range must DOUBLE, i.e. log N must reach")
    print(f"    {Ls[0] * (Ls[-1]/Ls[0])**2:.1f}, i.e. N ~ 1e"
          f"{Ls[0]*(Ls[-1]/Ls[0])**2/math.log(10):.0f}. "
          f"(se would go to ~{need/2:.4f})")
    print("    This is why the exponent is hard to pin: the observable")
    print("    grows like log log N, and no reachable N moves it much.")
    print("DONE")


if __name__ == "__main__":
    main()
