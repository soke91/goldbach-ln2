# -*- coding: utf-8 -*-
"""
lem:cellmom, prop:coh, and sec:floor's decay table.
(v1_verify2, Phase 1, blind.)

STATEMENTS UNDER TEST, verbatim from v1/paper/wall_v1.tex:

  sec:floor  "Throughout, cells are indexed by depth d, the number of
              3,5,7,11,13 dividing N, and Z(N)=C(N)/sqrt(V(N))."

  lem:cellmom  with u_c(v) = sum_{N in c} Lambda(N-v)/sqrt(V(N)) and
              Q_cd = sum_v mu^2(v) u_c(v) u_d(v),
              Var(m_c - mbar) = Q_cc/n_c^2 - 2 Q_ca/(n_c n) + Q_aa/n^2
              under independent signs. "Every term is computable exactly
              by one convolution; no simulation is involved."

  conj:wall item 2  "Against the exactly computed floor of
              Lemma cellmom it clears Bonferroni in every octave to
              1.6e7, at max_c |z_c| = 8.4 in the top band."

  prop:coh   Q_cc/n_c^2 ~ (6/pi^2)N / (A(N) N log N) ~ 1/log N, and
             "fitting se ~ N^{-b} gives b = 0.1105, 0.0465, 0.0385,
              0.0379, 0.0378, 0.0379 at depths 5,4,3,2,1,0"

  sec:floor  the decay table a_d = 0.1434, 0.2152, 0.2713, 0.3686,
             0.0437, 0.6289 at depths 5,4,3,2,1,0 with s.e. 0.0155,
             0.0065, 0.0040, 0.0052, 0.0556, 0.0121

  sec:floor  "Rarity. Its per-cell share reaches 0.94 in the deepest
             cell at N ~ 1.4e7, against a pooled share of 0.018."

PRE-REGISTRATION (fixed before this ran).

  Decision rule.
   (1) Verify lem:cellmom ITSELF by Monte Carlo: draw independent signs
       on the squarefree support, form Z_eps, and compare the empirical
       variance of m_c - mbar against the closed form.
       CONFIRMED iff the two agree to within the Monte Carlo error at
       every depth.
   (2) Compute the exact floor and z_c per octave; report max_c|z_c| in
       the top band and whether Bonferroni is cleared in every octave.
   (3) Fit se ~ N^{-b} per depth across octaves; compare the six b.
   (4) Fit |delta_d| ~ N^{-a_d} per depth, weighted by the exact errors,
       with each exponent's s.e. from the fit covariance; compare the
       six a_d and their s.e.
   (5) Report the per-cell share of the mask in the deepest cell and
       pooled.

  Predictions written before running.
   (1) CONFIRMED. The proof is two lines of bilinear expansion and I
       checked it by hand; E[eps(v)eps(v')] = mu^2(v) delta_{vv'} gives
       Cov(n_c m_c, n_d m_d) = Q_cd exactly.
   (2) I predict the floor is right and the mask is real -- the depth
       means measured earlier (-0.357 at d=3, -1.412 at d=4, -3.729 at
       d=5 in the top octave) are far too large to be sampling noise
       under any bar. I predict max|z_c| comes out much LARGER than 8.4,
       because with 6 depth cells the deepest cell's mean is enormous;
       8.4 looks like a figure for a different cell definition.
   (3) and (4): no prediction on the individual values. I do predict
       that depth 1's exponent is the least stable, since the paper
       itself reports it as unmeasurable.

  What would refute (2): max|z_c| = 8.4 on the stated cell definition.
"""

import os
import sys

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
CACHE = os.path.join(ROOT, "v1_verify2_log", "cache")
sys.path.insert(0, HERE)


def next_pow2(x):
    n = 1
    while n < x:
        n <<= 1
    return n


def correlate_to(w, fl, B, nfft):
    """u(v) = sum_N w[N] L[N-v], for v = 0..B, by circular correlation.

    fl is the precomputed conjugated rfft of L at length nfft, so the
    Lambda transform is paid for once per band and not once per cell.
    """
    fw = np.fft.rfft(w, nfft)
    fw *= fl
    r = np.fft.irfft(fw, nfft)
    del fw
    return r[: B + 1].copy()


def main():
    X = int(sys.argv[1]) if len(sys.argv) > 1 else 16_000_000
    z = np.load(os.path.join(CACHE, f"field_{X}.npz"))
    goodm = (z["V"] > 0) & (z["W"] > 0)
    Ni = z["N"][goodm]
    Vv = z["V"][goodm]
    Cc = z["C"][goodm]
    depth = z["depth"][goodm]
    Z = Cc / np.sqrt(Vv)

    print("audit_cell_floor   (v1_verify2 Phase 1, blind)")
    print("cells = the six depths, as sec:floor defines them")
    print("=" * 78)

    from lab_field_build import smallest_prime_factor, von_mangoldt, mobius
    print("  sieving ...")
    spf = smallest_prime_factor(X)
    lam, primes = von_mangoldt(X, spf)
    del spf
    mu = mobius(X, primes)
    mu2 = (mu != 0).astype(np.float64)
    del mu

    # map even N -> index
    pos = np.zeros(X + 1, dtype=np.int64)
    pos[Ni] = np.arange(len(Ni))

    def band_analysis(lo, hi, ntrial=0):
        """Exact floor for the band (lo, hi], cells = depths."""
        sel = (Ni > lo) & (Ni <= hi)
        Nb = Ni[sel]
        Zb = Z[sel]
        db = depth[sel]
        Vb = Vv[sel]
        n = len(Nb)
        B = int(hi)
        nfft = next_pow2(2 * B + 4)
        L = np.zeros(B + 1)
        L[1:] = lam[1: B + 1]
        fl = np.conjugate(np.fft.rfft(L, nfft))
        del L

        # whole band
        wa = np.zeros(B + 1)
        wa[Nb] = 1.0 / np.sqrt(Vb)
        ua = correlate_to(wa, fl, B, nfft)
        m2 = mu2[: B + 1]
        Qaa = float(np.dot(m2 * ua, ua))

        out = []
        for d in range(6):
            m = db == d
            nc = int(m.sum())
            if nc < 3:
                out.append(None)
                continue
            wc = np.zeros(B + 1)
            wc[Nb[m]] = 1.0 / np.sqrt(Vb[m])
            uc = correlate_to(wc, fl, B, nfft)
            Qcc = float(np.dot(m2 * uc, uc))
            Qca = float(np.dot(m2 * uc, ua))
            var = Qcc / nc ** 2 - 2 * Qca / (nc * n) + Qaa / n ** 2
            delta = float(Zb[m].mean() - Zb.mean())
            out.append((d, nc, delta, np.sqrt(max(var, 0.0)), Qcc / nc ** 2))
            del wc, uc
        del wa, ua, fl

        mc_res = None
        if ntrial:
            # ---- (1) Monte Carlo check of lem:cellmom -----------------
            from lab_field_build import fftconv_prefix
            rng = np.random.default_rng(99)
            supp = mu2 > 0
            acc = {d: [] for d in range(6)}
            for _ in range(ntrial):
                eps = np.zeros(X + 1)
                idx = np.nonzero(supp)[0]
                eps[idx] = rng.integers(0, 2, size=len(idx)) * 2.0 - 1.0
                Ce = fftconv_prefix(lam, eps, X)
                Ze = Ce[Nb] / np.sqrt(Vb)
                for d in range(6):
                    m = db == d
                    if m.sum() >= 3:
                        acc[d].append(Ze[m].mean() - Ze.mean())
                del eps, Ce, Ze
            mc_res = {d: (np.std(v, ddof=1) if len(v) > 2 else np.nan)
                      for d, v in acc.items()}
        return n, out, mc_res

    # ------------------------------------------------ (1)
    print()
    print("--- (1) lem:cellmom against Monte Carlo, band (1e5, 2e5] --------")
    n, out, mc = band_analysis(1e5, 2e5, ntrial=60)
    print(f"    band size n = {n:,};  60 independent-sign draws")
    print(f"    {'depth':>6}{'n_c':>8}{'closed form se':>17}"
          f"{'Monte Carlo sd':>17}{'ratio':>9}")
    for r in out:
        if r is None:
            continue
        d, nc, delta, se, _ = r
        m = mc[d]
        print(f"    {d:>6}{nc:>8,}{se:>17.6f}{m:>17.6f}{m / se:>9.3f}")
    print("    (ratio ~ 1 +- 1/sqrt(2*60) = 1 +- 0.091 confirms the lemma)")

    # ------------------------------------------------ (2)-(5)
    print()
    print("--- (2) the exact floor, and z_c, by octave ---------------------")
    bands = []
    b = X
    while b > 1e5:
        bands.append((b / 2.0, b))
        b /= 2.0
    bands = bands[::-1]
    store = {}
    for lo, hi in bands:
        n, out, _ = band_analysis(lo, hi)
        store[hi] = (n, out)
        zs = [abs(r[2] / r[3]) for r in out if r is not None and r[3] > 0]
        bonf = 0.05 / max(len([r for r in out if r is not None]), 1)
        # two-sided normal critical value for Bonferroni
        from math import erfc, sqrt
        crit = 0.0
        lo_, hi_ = 0.0, 20.0
        for _ in range(200):
            mid = 0.5 * (lo_ + hi_)
            if erfc(mid / sqrt(2.0)) > bonf:
                lo_ = mid
            else:
                hi_ = mid
        crit = 0.5 * (lo_ + hi_)
        print(f"    N in ({lo:>10,.0f}, {hi:>10,.0f}]  n={n:>9,}   "
              f"max|z_c| = {max(zs):>8.1f}   Bonferroni crit = {crit:.2f}   "
              f"{'clears' if max(zs) > crit else 'DOES NOT CLEAR'}")
    print(f"    [paper: 'clears Bonferroni in every octave to 1.6e7, at")
    print(f"     max_c|z_c| = 8.4 in the top band']")

    print()
    print("--- the top band in detail --------------------------------------")
    n, out = store[X]
    print(f"    {'depth':>6}{'n_c':>10}{'delta_c':>12}{'se (exact)':>13}"
          f"{'z_c':>10}{'Qcc/n_c^2':>12}")
    for r in out:
        if r is None:
            continue
        d, nc, delta, se, q = r
        print(f"    {d:>6}{nc:>10,}{delta:>+12.5f}{se:>13.6f}"
              f"{delta / se:>+10.1f}{q:>12.3e}")
    print()
    print(f"    prop:coh predicts Qcc/n_c^2 ~ (6/pi^2)/(A log N):")
    A_even = 0.3739558136 / 0.5
    pred = (6 / np.pi ** 2) / (A_even * np.log(X))
    print(f"      (6/pi^2)/(A log N) at N={X:,} = {pred:.3e}")

    # ------------------------------------------------ (3),(4)
    print()
    print("--- (3) se ~ N^-b, and (4) |delta_d| ~ N^-a_d -------------------")
    print(f"    {'depth':>6}{'b':>10}{'paper b':>10}{'a_d':>10}"
          f"{'s.e.':>9}{'paper a_d':>11}{'paper s.e.':>11}")
    pb = {5: 0.1105, 4: 0.0465, 3: 0.0385, 2: 0.0379, 1: 0.0378, 0: 0.0379}
    pa = {5: (0.1434, 0.0155), 4: (0.2152, 0.0065), 3: (0.2713, 0.0040),
          2: (0.3686, 0.0052), 1: (0.0437, 0.0556), 0: (0.6289, 0.0121)}
    for d in (5, 4, 3, 2, 1, 0):
        xs, ses, dels, errs = [], [], [], []
        for lo, hi in bands:
            n, out = store[hi]
            r = out[d]
            if r is None or r[3] <= 0:
                continue
            xs.append(np.sqrt(lo * hi))
            ses.append(r[3])
            dels.append(abs(r[2]))
            errs.append(r[3])
        if len(xs) < 3:
            print(f"    {d:>6}  too few bands")
            continue
        lx = np.log(np.array(xs))
        bfit = -np.polyfit(lx, np.log(np.array(ses)), 1)[0]
        # weighted fit of log|delta| on log N, weights from exact errors
        y = np.log(np.array(dels))
        sy = np.array(errs) / np.array(dels)
        w = 1.0 / sy ** 2
        Sw, Sx = w.sum(), (w * lx).sum()
        Sxx, Sy, Sxy = (w * lx * lx).sum(), (w * y).sum(), (w * lx * y).sum()
        den = Sw * Sxx - Sx * Sx
        slope = (Sw * Sxy - Sx * Sy) / den
        slope_se = np.sqrt(Sw / den)
        print(f"    {d:>6}{bfit:>10.4f}{pb[d]:>10.4f}{-slope:>10.4f}"
              f"{slope_se:>9.4f}{pa[d][0]:>11.4f}{pa[d][1]:>11.4f}")
    print()
    print(f"    prop:coh's predicted apparent exponent 1/(2<log N>) over")
    print(f"    1e5..1.4e7 = {1 / (2 * 0.5 * (np.log(1e5) + np.log(1.4e7))):.4f}")

    # ------------------------------------------------ (5)
    print()
    print("--- (5) rarity: the mask's per-cell share ------------------------")
    n, out = store[X]
    tot_var = float(Z[(Ni > X / 2)].var(ddof=1))
    for r in out:
        if r is None:
            continue
        d, nc, delta, se, q = r
        m = (Ni > X / 2) & (depth == d)
        share = delta ** 2 / float(Z[m].var(ddof=1) + delta ** 2)
        print(f"    depth {d}: n_c={nc:>9,}  delta={delta:+.4f}  "
              f"share of that cell's mean square = {share:.4f}")
    pooled = sum((r[2] ** 2) * r[1] for r in out if r is not None) / \
        sum(r[1] for r in out if r is not None) / tot_var
    print(f"    pooled share = {pooled:.4f}   "
          f"[paper: 0.94 in the deepest cell, 0.018 pooled]")
    return 0


if __name__ == "__main__":
    sys.exit(main())
