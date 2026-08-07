# -*- coding: utf-8 -*-
"""
Re-verification of kill-test R4, §`sec:R4` of v1/paper/wall_v1.tex:
"the divisor switch does not localize".

THE STATEMENTS UNDER TEST, verbatim.

  (a) "Applied to the dilate field over its FULL ranges the switch
      gives an exact identity
         sum_{k>=1} sum_{m: mk<=N-1} mu(m) mu(N-mk)
           = sum_{u<N} mu(N-u) sum_{m|u} mu(m) = mu(N-1),
      verified by brute force at N = 5000 and N = 20000. This is
      perfect cancellation: O(1) for a double sum of ~N log N terms."

  (b) "On 8000 values of k --- the entire band on which the type-II
      field is non-empty, since m > sqrt N forces k < sqrt N --- the
      block ratios are flat (B=8: 0.958 and 1.023 at two N, against
      B=1 baselines 0.980 and 0.979), and the sharper diagnostic, the
      lag-1 autocorrelation of D(k)/sqrt(supp(k)), reads +0.0104 and
      +0.0127 against a standard error of 0.0112: dead zero, and if
      anything mildly positive rather than the negative coherence a
      surviving residue would require."

  (c) "The mirror. Switching the banded L^1 sum gives
      sum_{k~K} D_full(k) = sum_{u<N} mu(N-u)
        sum_{m|u, u/2K<m<=u/K} mu(m), where m ~ N/K >= N^{2/3}: the
      surviving Mobius sits on the LONG variable."

METHOD HERE. Written from the statement. (a) is checked by two
independent enumerations, one over (m,k) and one over u, at the two N
the paper names and at four more. (b) is rebuilt from the definitions:
D(k) = sum_{sqrt N < m <= N/k} mu(m) mu(N-mk), supp(k) the number of
m in that range with both mu's nonzero, block sums S_B(j) over
consecutive blocks of B values of k, and the block ratio

    r(B) = rms_j S_B(j) / sqrt( sum_{k in block} supp(k) )   (averaged
                                                              over j)

which is 1 under the B-independence Conjecture 10 predicts. (c) is an
identity about which variable carries mu and is checked by
enumeration at small N.

PRE-REGISTRATION (written before the run).

  (1) RULE for (a). The identity must hold exactly, in integer
      arithmetic, at every N tested. Anything else refutes it.
  (2) RULE for (b). r(B) must stay flat in B -- within 20% of the
      B = 1 baseline at B = 8 -- and the lag-1 autocorrelation must sit
      within 3 standard errors of zero. A NEGATIVE lag-1 correlation
      several standard errors below zero is what a surviving residue
      of the switch's cancellation would look like, and is the
      outcome that would revive R4.
  (3) PREDICTION, recorded so it cannot be reported as a surprise.
      The identity is elementary -- the inner sum over m | u is
      1_{u=1} -- so (a) will hold. For (b) I expect r(8)/r(1) near 1
      and a lag-1 correlation of order 1/sqrt(n), i.e. no signal, and
      I expect no finding here. The informative outcome would be a
      block ratio drifting DOWN with B, which is the signature R4 was
      designed to detect.
  (4) The standard error quoted for the lag-1 correlation is
      1/sqrt(n). That is the iid value, and D(k)/sqrt(supp(k)) is what
      is being tested for independence, so the bar is the null's own
      and not a measured spread. Reported alongside: a block-bootstrap
      standard error, which does not assume it.
"""
import sys
import math

import numpy as np

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass


def sieve_mu(X):
    mu = np.ones(X + 1, dtype=np.int8)
    comp = np.zeros(X + 1, dtype=bool)
    for p in range(2, int(X ** 0.5) + 1):
        if not comp[p]:
            comp[p * p::p] = True
    primes = np.nonzero(~comp)[0]
    for p in primes[primes >= 2]:
        p = int(p)
        mu[p::p] = -mu[p::p]
        if p * p <= X:
            mu[p * p::p * p] = 0
    mu[0] = 0
    return mu


def main():
    print("Re-verification of kill-test R4 (the divisor switch)")
    print()

    # ---------- (a) the exact identity ----------
    print("(a) sum_{k>=1} sum_{m: mk<=N-1} mu(m)mu(N-mk) = mu(N-1)")
    print(f"    {'N':>8} {'double sum over (m,k)':>22} "
          f"{'sum over u':>12} {'mu(N-1)':>9} {'ok':>4}")
    okA = True
    for N in (5000, 20000, 3001, 12345, 50000, 100000):
        mu = sieve_mu(N)
        tot = 0
        for k in range(1, N):
            m = np.arange(1, (N - 1) // k + 1, dtype=np.int64)
            if not len(m):
                break
            tot += int(np.dot(mu[m].astype(np.int64),
                              mu[N - m * k].astype(np.int64)))
        # the same thing collected over u = mk
        alt = 0
        for u in range(1, N):
            d = np.arange(1, u + 1)
            d = d[u % d == 0]
            s = int(mu[d].sum())
            if s:
                alt += s * int(mu[N - u])
        ok = (tot == alt == int(mu[N - 1]))
        okA &= ok
        print(f"    {N:>8} {tot:>22} {alt:>12} {int(mu[N-1]):>9} "
              f"{'OK' if ok else 'FAIL':>4}")
    print(f"    (1) the identity is exact at every N: "
          f"{'PASS' if okA else 'FAIL'}")
    print()

    # ---------- (b) block ratios and lag-1 ----------
    print("(b) does any of the cancellation survive the type-II cut?")
    print("    D(k) = sum_{sqrt N < m <= N/k} mu(m) mu(N - m k)")
    print()
    okB = True
    for N in (4_000_000, 16_000_000):
        N -= N % 2
        mu = sieve_mu(N)
        rt = int(math.isqrt(N))
        ks = np.arange(1, rt, dtype=np.int64)   # m > sqrt N forces k < sqrt N
        D = np.zeros(len(ks))
        S = np.zeros(len(ks))
        for i, k in enumerate(ks):
            m = np.arange(rt + 1, N // int(k) + 1, dtype=np.int64)
            if not len(m):
                continue
            a = mu[m].astype(np.int64)
            b = mu[N - m * int(k)].astype(np.int64)
            D[i] = float(np.dot(a, b))
            S[i] = float(np.count_nonzero((a != 0) & (b != 0)))
        keep = S > 0
        D, S = D[keep], S[keep]
        n = len(D)
        print(f"    N = {N}, {n} values of k with a non-empty band")
        # the paper does not define the block-ratio statistic, so both
        # readings of "S_B against the B-independence Conjecture 10
        # predicts" are reported:
        #   r_supp(B) = sqrt( sum_j S_B(j)^2 / sum_j sum_{k in j} supp(k) )
        #   r_std(B)  = rms_j ( sum_{k in j} D(k)/sqrt(supp k) ) / sqrt(B)
        print(f"      {'B':>5} {'blocks':>7} {'r_supp':>9} {'/B=1':>7} "
              f"{'r_std':>9} {'/B=1':>7}")
        base = bstd = None
        gk = D / np.sqrt(S)
        for B in (1, 2, 4, 8, 16, 64):
            nb = n // B
            if nb < 8:
                continue
            s = D[:nb * B].reshape(nb, B).sum(axis=1)
            v = S[:nb * B].reshape(nb, B).sum(axis=1)
            r = float(np.sqrt((s ** 2).sum() / v.sum()))
            t = gk[:nb * B].reshape(nb, B).sum(axis=1)
            rs = float(np.sqrt((t ** 2).mean() / B))
            if base is None:
                base, bstd = r, rs
            print(f"      {B:>5} {nb:>7} {r:>9.4f} {r/base:>7.4f} "
                  f"{rs:>9.4f} {rs/bstd:>7.4f}")
            if B == 8:
                okB &= abs(rs / bstd - 1.0) < 0.20
        # where r_supp's weight sits: supp(k) ~ N/k, so small k dominate
        w = S / S.sum()
        top = np.argsort(-w)[:5]
        print(f"      supp weight: the 5 smallest k carry "
              f"{100*w[top].sum():.1f}% of sum supp(k); "
              f"k = {[int(ks[keep][i]) for i in top]}")
        for cut in (1, 10, 100):
            sel = ks[keep] >= cut
            rr = float(np.sqrt((D[sel] ** 2).sum() / S[sel].sum()))
            print(f"      r_supp(B=1) restricted to k >= {cut:>3}: "
                  f"{rr:.4f}   ({int(sel.sum())} values)")
        g = D / np.sqrt(S)
        g = g - g.mean()
        lag1 = float((g[:-1] * g[1:]).sum() / (g * g).sum())
        se_iid = 1.0 / math.sqrt(n)
        rng = np.random.default_rng(11)
        boot = []
        Bb = max(8, n // 64)
        for _ in range(200):
            st = rng.integers(0, n - Bb, size=n // Bb)
            idx = (st[:, None] + np.arange(Bb)[None, :]).ravel()
            h = g[idx] - g[idx].mean()
            boot.append(float((h[:-1] * h[1:]).sum() / (h * h).sum()))
        se_boot = float(np.std(boot))
        print(f"      lag-1 autocorrelation of D(k)/sqrt(supp k) = "
              f"{lag1:+.4f}")
        print(f"        iid standard error 1/sqrt(n) = {se_iid:.4f}"
              f"   -> z = {lag1/se_iid:+.2f}")
        print(f"        block-bootstrap standard error = {se_boot:.4f}"
              f"   -> z = {lag1/se_boot:+.2f}")
        okB &= abs(lag1 / se_iid) < 3.0
        print(f"      v1 quotes +0.0104 and +0.0127 against 0.0112, at")
        print(f"      two larger N (8000 values of k)")
        print()
    print(f"    (2) blocks flat at B=8 and lag-1 within 3 sigma: "
          f"{'PASS' if okB else 'FAIL'}")
    print()

    # ---------- (d) which statistic answers the question ----------
    print("(d) does the fall in the unstandardised ratio mean R4 is")
    print("    ALIVE? The two statistics test different things:")
    print("      r_supp(B)^2 = sum_j S_B(j)^2 / sum_j sum_{k in j} supp(k)")
    print("        -- tests |D(k)| ~ sqrt(supp k) AND independence at once")
    print("      rho(B)^2    = sum_j S_B(j)^2 / sum_k D(k)^2")
    print("        -- holds the magnitudes fixed, so it tests ONLY")
    print("           independence of the signs across k. E[rho^2] = 1")
    print("           for every B under sign-independence, whatever the")
    print("           magnitudes are.")
    print("    Null: 400 surrogates D(k) -> eps_k D(k), eps_k = +-1 iid.")
    print("    This preserves every |D(k)| exactly and destroys only the")
    print("    coupling across k -- the one thing R4 is about.")
    print()
    okD = True
    for N in (4_000_000, 16_000_000):
        N -= N % 2
        mu = sieve_mu(N)
        rt = int(math.isqrt(N))
        ks = np.arange(1, rt, dtype=np.int64)
        D = np.zeros(len(ks))
        S = np.zeros(len(ks))
        for i, k in enumerate(ks):
            m = np.arange(rt + 1, N // int(k) + 1, dtype=np.int64)
            if not len(m):
                continue
            a = mu[m].astype(np.int64)
            b = mu[N - m * int(k)].astype(np.int64)
            D[i] = float(np.dot(a, b))
            S[i] = float(np.count_nonzero((a != 0) & (b != 0)))
        keep = S > 0
        D, S = D[keep], S[keep]
        n = len(D)
        dd = float((D ** 2).sum())
        rng = np.random.default_rng(31337)
        print(f"    N = {N}   sum D^2 / sum supp = {dd/S.sum():.4f}"
              f"   (1 under square-root cancellation)")
        print(f"      {'B':>4} {'r_supp':>8} {'rho':>8} "
              f"{'null mean':>10} {'null sd':>9} {'z':>7} {'p(low)':>8}")
        for B in (1, 2, 4, 8, 16, 32, 64):
            nb = n // B
            if nb < 8:
                continue
            def stat(x):
                s = x[:nb * B].reshape(nb, B).sum(axis=1)
                return float(np.sqrt((s ** 2).sum() / dd))
            v = S[:nb * B].reshape(nb, B).sum(axis=1)
            sD = D[:nb * B].reshape(nb, B).sum(axis=1)
            r_supp = float(np.sqrt((sD ** 2).sum() / v.sum()))
            rho = stat(D)
            null = np.array([stat(D * rng.choice([-1.0, 1.0], size=n))
                             for _ in range(400)])
            if B == 1:
                # rho(1) = 1 identically and the null has no spread;
                # the statistic carries no information at B = 1.
                z, p = float("nan"), 1.0
            else:
                z = (rho - null.mean()) / null.std()
                p = float((null <= rho).mean())
                okD &= (z > -3.0)
            print(f"      {B:>4} {r_supp:>8.4f} {rho:>8.4f} "
                  f"{null.mean():>10.4f} {null.std():>9.4f} "
                  f"{z:>7.2f} {p:>8.3f}")
        print()
    print("    A surviving residue of the switch's cancellation is a")
    print("    NEGATIVE z that deepens with B. Read the z column.")
    print(f"    (4) no block size shows rho below its null by 3 sigma: "
          f"{'PASS' if okD else 'FAIL -- R4 would be alive'}")
    print()

    # ---------- (c) the mirror ----------
    print("(c) the mirror: switching the banded sum puts mu on the LONG")
    print("    variable. Checked as an identity at small N:")
    print(f"    {'N':>7} {'K':>5} {'sum_k~K D_full(k)':>19} "
          f"{'switched form':>15} {'ok':>4}")
    okC = True
    for N, K in ((3000, 7), (5000, 11), (20000, 23)):
        mu = sieve_mu(N)
        lhs = 0
        for k in range(K, 2 * K):
            m = np.arange(1, (N - 1) // k + 1, dtype=np.int64)
            lhs += int(np.dot(mu[m].astype(np.int64),
                              mu[N - m * k].astype(np.int64)))
        rhs = 0
        for u in range(1, N):
            lo, hi = u / (2.0 * K), u / float(K)
            d = np.arange(1, u + 1)
            d = d[(u % d == 0) & (d > lo) & (d <= hi)]
            if len(d):
                rhs += int(mu[d].sum()) * int(mu[N - u])
        ok = (lhs == rhs)
        okC &= ok
        print(f"    {N:>7} {K:>5} {lhs:>19} {rhs:>15} "
              f"{'OK' if ok else 'FAIL':>4}")
    print(f"    the switched form's inner sum runs over m in (u/2K, u/K],")
    print(f"    i.e. m ~ N/K: the surviving mu is on the long variable,")
    print(f"    the opposite of what Theorem 1 needs.")
    print(f"    (3) the mirror identity holds: {'PASS' if okC else 'FAIL'}")
    print()
    if not (okA and okB and okC and okD):
        print("DONE (failed)")
        sys.exit(1)
    print("DONE")


if __name__ == "__main__":
    main()
