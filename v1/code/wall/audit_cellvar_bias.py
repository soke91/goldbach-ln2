# -*- coding: utf-8 -*-
"""
Where the +0.020 came from (increment 319)

#152 recorded a residual in the closed-form cell variance: over 16 coin
draws, `mean(d^2 - sd_c^2) = +0.0204 +/- 0.0029`, `z = +7.11`, source
not identified. That formula --

    Var(m_c - gm) = Q_cc/n_c^2 - 2 Q_ca/(n_c n) + Q_aa/n^2

-- is what increments 305, 314 and 318 lean on, so an unexplained bias
in it is not something to leave standing.

THE ALGEBRA IS EXACT. Sum_{N in c} Z = sum_v eps(v) u_c(v) with
u_c(v) = sum_{N in c} Lambda(N-v)/sqrt(V(N)), and since eps is
independent +/-1 on {mu != 0}, E[(sum_c Z)(sum_a Z)] = Q_ca exactly.
There is no approximation to be wrong.

WHICH LEAVES THE STANDARD ERROR. The 4096 numbers averaged in #152 are
16 draws x 8 bands x 32 cells, and the cells within one draw are a
PARTITION OF THE SAME REALISATION of eps. They are not independent, and
dividing by sqrt(4096) treats them as if they were. If the effective
sample is the 16 draws, the standard error is larger by up to
sqrt(4096/16) = 16 and z falls from 7.11 to about 0.4.

That is the same species as #41 and as the E1 seam estimator: an SE
written as std/sqrt(n) on samples that share a source.

PRE-REGISTRATION (fixed before the run).

  (U1) THE CLUSTERED TEST. Average the debiased square over all cells
       and bands WITHIN each draw, then test those R per-draw numbers.
       RULE: |z_clustered| < 3. If it passes, #152's z = 7.11 was an
       artefact of the standard error and the closed form is vindicated;
       if it fails, the bias is real and the formula is wrong.

  (U2) THE EFFECTIVE SAMPLE SIZE, measured rather than argued:
       n_eff = (sd_naive / sd_clustered)^2 * R. RULE: report it. If it
       comes out near 4096 the cells really were independent and (U1)
       cannot excuse anything; if near R, they were one sample.

  (U3) SIGN CHECK ON A SECOND STATISTIC. #152 noted that increment
       305's z-based self-test gave sd(z) = 0.9936 < 1, which implies
       Var(d) BELOW sd_c^2 -- the opposite sign to the +0.020. Both are
       computed here from the same draws. RULE: if the clustered test
       passes, the two must be consistent with zero simultaneously.

  WHAT WOULD REFUTE. (U1) failing with (U2) near 4096 means the
  closed form has a real positive bias and #152 stands as written.
"""
import math
import sys
import time

import numpy as np

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

QS = [3, 5, 7, 11, 13]
R = 24


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
    return mu, lam


def main():
    X = 16_000_000
    lo = 100_000
    t0 = time.time()
    mu, lam = sieve(X)
    nfft = 1
    while nfft < 2 * (X + 1):
        nfft *= 2
    supp = (mu != 0)
    suppf = supp.astype(np.float64)
    F_supp = np.fft.rfft(np.pad(suppf, (0, nfft - X - 1)))
    F_lam = np.fft.rfft(np.pad(lam, (0, nfft - X - 1)))
    Fl_c = np.conj(F_lam)
    V = np.fft.irfft(F_supp * np.fft.rfft(
        np.pad(lam ** 2, (0, nfft - X - 1))), nfft)[: X + 1]
    Ns = np.arange(lo, X + 1, 2)
    key = np.zeros(len(Ns), dtype=np.int32)
    for i, q in enumerate(QS):
        key |= ((Ns % q) == 0).astype(np.int32) << i
    invV = 1.0 / np.sqrt(V[Ns])
    muw = suppf[: X + 1]
    print(f"sieve + V  t={time.time()-t0:.0f}s", flush=True)

    sels = []
    b = lo
    while b < X:
        hi = min(2 * b, X)
        sel = (Ns >= b) & (Ns < hi)
        if int(sel.sum()) > 1000:
            sels.append(sel)
        b = hi

    SD = []
    for i, sel in enumerate(sels):
        n = int(sel.sum())
        kb = key[sel]
        u_all = np.zeros(X + 1)
        us, ns = [], []
        for c in np.unique(kb):
            m = sel.copy()
            m[sel] = (kb == c)
            w = np.zeros(nfft)
            w[Ns] = np.where(m, invV, 0.0)
            u = np.fft.irfft(Fl_c * np.fft.rfft(w), nfft)[: X + 1]
            us.append(u)
            ns.append(int(m.sum()))
            u_all += u
        mu_all = muw * u_all
        Qaa = float(np.dot(mu_all, u_all))
        sd = np.zeros(len(us))
        for j, u in enumerate(us):
            nc = ns[j]
            var = (float(np.dot(muw * u, u)) / nc ** 2
                   - 2 * float(np.dot(mu_all, u)) / (nc * n)
                   + Qaa / n ** 2)
            sd[j] = math.sqrt(max(var, 0.0))
        SD.append(sd)
        print(f"  band {i+1}/{len(sels)}  t={time.time()-t0:.0f}s",
              flush=True)

    idx = np.nonzero(supp)[0]
    rng = np.random.default_rng(319)
    per_cell, per_draw, zs_all = [], [], []
    for r in range(R):
        eps = np.zeros(nfft)
        eps[idx] = rng.integers(0, 2, size=len(idx)) * 2.0 - 1.0
        C = np.fft.irfft(np.fft.rfft(eps) * F_lam, nfft)[: X + 1]
        Z = C[Ns] * invV
        acc = []
        for i, sel in enumerate(sels):
            z = Z[sel]
            k = key[sel]
            _uq, inv = np.unique(k, return_inverse=True)
            cnt = np.bincount(inv).astype(np.float64)
            tot = np.bincount(inv, weights=z)
            d = tot / cnt - float(z.mean())
            db = d ** 2 - SD[i] ** 2
            per_cell.append(db)
            acc.append(db)
            zs_all.append(d / SD[i])
        per_draw.append(float(np.concatenate(acc).mean()))

    flat = np.concatenate(per_cell)
    pd = np.array(per_draw)
    m_naive = float(flat.mean())
    se_naive = float(flat.std(ddof=1)) / math.sqrt(len(flat))
    m_clu = float(pd.mean())
    se_clu = float(pd.std(ddof=1)) / math.sqrt(R)
    z_naive = m_naive / se_naive
    z_clu = m_clu / se_clu
    n_eff = (se_naive / se_clu) ** 2 * R if se_clu > 0 else float("nan")

    zz = np.concatenate(zs_all)
    zmean, zsd = float(zz.mean()), float(zz.std(ddof=1))

    print(f"\n(U1) the same residual, two standard errors")
    print(f"    naive, {len(flat)} cell-draws treated as independent:")
    print(f"        {m_naive:+.5f} +/- {se_naive:.5f}   z = {z_naive:+.2f}")
    print(f"    clustered, {R} draws (cells within a draw partition one")
    print(f"    realisation of eps and are not independent):")
    print(f"        {m_clu:+.5f} +/- {se_clu:.5f}   z = {z_clu:+.2f}")
    okU1 = abs(z_clu) < 3.0
    print(f"    (U1) |z_clustered| < 3: {'PASS' if okU1 else 'FAIL'}")

    print(f"\n(U2) effective sample size = "
          f"(se_naive/se_clustered)^2 * R = {n_eff:.1f}")
    print(f"     against {len(flat)} nominal and {R} draws -- "
          f"{'one sample per draw' if n_eff < 4 * R else 'partly independent'}")

    print(f"\n(U3) the z-based statistic on the same draws: "
          f"mean {zmean:+.4f}, sd {zsd:.4f}")
    # 첫 판은 여기서 "sd below 1"을 단언했다 — 그건 증분 305의
    # 숫자였고 이 실행의 것이 아니다. 부호는 데이터에서 계산한다.
    direc = "BELOW" if zsd < 1.0 else "ABOVE"
    print(f"     sd {'<' if zsd < 1.0 else '>'} 1 means Var(d) is "
          f"{direc} sd_c^2 by {abs(zsd**2 - 1):.1%};")
    print(f"     with n_eff = {n_eff:.1f} the sampling error on a "
          f"dispersion is large, so this")
    print(f"     is not a 3-sigma statement either way.")
    okU3 = abs(zsd - 1.0) < 0.05
    print(f"     (U3) sd(z) within 0.05 of 1: "
          f"{'PASS' if okU3 else 'FAIL'}")

    if okU1 and n_eff < 4 * R:
        v = (f"#152's z = +7.11 was an artefact of the standard error, "
             f"not a bias in the closed form. The 4096 cell-draws are "
             f"{n_eff:.0f} independent numbers, not 4096: the cells "
             f"within a draw partition one realisation of eps. "
             f"Clustered, the residual is {m_clu:+.4f} +/- {se_clu:.4f}, "
             f"z = {z_clu:+.2f}. The formula stands and every share "
             f"computed from it at 318 stands with it")
    elif okU1:
        v = ("the clustered test passes but the cell-draws are not far "
             "from independent, so the standard error is not the whole "
             "story and the residual needs a second explanation")
    else:
        v = ("the residual survives clustering: the closed form has a "
             "real positive bias and #152 stands as written")
    print(f"\n    {v}")
    print("DONE")


if __name__ == "__main__":
    main()
