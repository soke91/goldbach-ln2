# -*- coding: utf-8 -*-
"""
Re-verification of Proposition 8 (`prop:Dpp`, polynomial weights) of
v1/paper/wall_v1.tex and its measured table.

THE STATEMENTS UNDER TEST (wall_v1.tex §2.3 and theorem_A.tex §6.3).

    Measured at N = 1e6, 4e6, 1.6e7: the two pieces of CP_2 are the
    same order, with ratio 0.771, 0.790, 0.810 drifting towards 1, and
    CP_2/(N log N) = 2.886, 2.949, 2.997.

      N                                  1e6     4e6    1.6e7
      r=1: sum_p Lambda(N-p) log^2 p /N  22.51   25.04  27.46
      r=2: 2 sum_{pq} Lambda(N-pq)
                       log p log q /N    17.36   19.79  22.25
      ratio r2/r1                        0.771   0.790  0.810
      CP_2/(N log N)                     2.886   2.949  2.997

    Calibration: sum_p Lambda(N-p) log p = 1.7565N, 1.7633N, 1.7614N
    against S(1e6) = 1.7604.

    The canonical tuning f(x) = x^2 - 2 gamma x moves CP from
    39.88N, 44.84N, 49.71N to 37.85N, 42.80N, 47.67N -- about five
    percent.

METHOD HERE. Written from the statement. b = mu * log^D = Lambda_D, and
on a squarefree u the value is elementary: Lambda_2(p) = log^2 p,
Lambda_2(pq) = 2 log p log q, and 0 for omega(u) >= 3. So CP_2 is built
directly from a smallest-prime-factor sieve, with no convolution and no
appeal to v1's implementation. The tuned weight is
b = Lambda_2 - 2 gamma Lambda, so CP_tuned = CP_2 - 2 gamma * G where
G = sum_p Lambda(N-p) log p is the calibration row itself -- an
internal consistency relation the table must satisfy.

PRE-REGISTRATION (written before the run).

  (1) RULE. Every quoted entry must reproduce to the digits printed.
      A disagreement is the finding.
  (2) The internal relation CP_tuned = CP_2 - 2 gamma G must hold to
      the printed digits; if the table satisfies it, the tuning row is
      not an independent measurement but arithmetic on the other two.
  (3) SEPARATE POINT, on the argument rather than the numbers.
      Proposition 8(ii) says "by classical sieve bounds
      CP_D ~ N (log N)^{D-1}". For D = 2 that is CP_2 ~ N log N, so
      CP_2/(N log N) should approach a CONSTANT. The measured column
      2.886, 2.949, 2.997 rises monotonically over a factor 16 in N.
      RULE: report the fitted drift. If the column is still climbing,
      the claim CP_2 ~ N log N is a statement the data have not
      reached, and the closure rests on nonnegativity alone -- which is
      what the paper says anyway, so this is a check on the wording of
      (ii), not on the closure.
"""
import sys
import math

import numpy as np

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

GAMMA = 0.5772156649015329


def sieve(X):
    """smallest prime factor, mu, Lambda."""
    spf = np.zeros(X + 1, dtype=np.int32)
    for i in range(2, int(X ** 0.5) + 1):
        if spf[i] == 0:
            sl = spf[i * i::i]
            sl[sl == 0] = i
    idx = np.arange(X + 1, dtype=np.int32)
    spf[spf == 0] = idx[spf == 0]
    spf[0] = spf[1] = 0
    is_p = np.zeros(X + 1, dtype=bool)
    is_p[2:] = spf[2:] == idx[2:]
    lam = np.zeros(X + 1, dtype=np.float64)
    for p in range(2, X + 1):
        if is_p[p]:
            lg = math.log(p)
            q = p
            while q <= X:
                lam[q] = lg
                q *= p
    return spf, is_p, lam


def lambda2_squarefree(X, spf, is_p):
    """Lambda_2 on squarefree u <= X: log^2 p at u=p, 2 log p log q at
    u = pq with p != q, 0 elsewhere. Also returns Lambda restricted to
    squarefree u (= log p at u = p)."""
    L2 = np.zeros(X + 1, dtype=np.float64)
    L1 = np.zeros(X + 1, dtype=np.float64)
    logs = np.zeros(X + 1, dtype=np.float64)
    logs[2:] = np.log(np.arange(2, X + 1, dtype=np.float64))
    for p in range(2, X + 1):
        if not is_p[p]:
            continue
        lp = math.log(p)
        L2[p] = lp * lp
        L1[p] = lp
        # u = p*q with q > p prime, q <= X/p
        if p * p > X:
            continue
        qs = np.arange(p + 1, X // p + 1)
        qs = qs[is_p[qs]]
        if len(qs):
            L2[p * qs] = 2.0 * lp * np.log(qs.astype(np.float64))
    return L1, L2


def main():
    print("Re-verification of Proposition 8's measured table")
    print()
    NS = [1_000_000, 4_000_000, 16_000_000]
    X = max(NS)
    spf, is_p, lam = sieve(X)
    L1, L2 = lambda2_squarefree(X, spf, is_p)

    hdr = (f"{'N':>10} {'r=1 /N':>9} {'r=2 /N':>9} {'r2/r1':>7} "
           f"{'CP_2/N':>9} {'CP_2/(NlogN)':>13} {'G/N':>8} "
           f"{'CPtuned/N':>10}")
    print(hdr)
    print("-" * len(hdr))
    quoted = {1_000_000: (22.51, 17.36, 0.771, 2.886, 1.7565, 37.85),
              4_000_000: (25.04, 19.79, 0.790, 2.949, 1.7633, 42.80),
              16_000_000: (27.46, 22.25, 0.810, 2.997, 1.7614, 47.67)}
    rows = []
    for N in NS:
        u = np.arange(1, N)
        w = lam[N - u]
        r1 = float((w * np.where(is_p[u], L2[u], 0.0)).sum())
        r2 = float((w * np.where(~is_p[u], L2[u], 0.0)).sum())
        G = float((w * L1[u]).sum())
        cp = r1 + r2
        cpt = cp - 2.0 * GAMMA * G
        rows.append((N, cp, G))
        print(f"{N:>10} {r1/N:>9.2f} {r2/N:>9.2f} {r2/r1:>7.3f} "
              f"{cp/N:>9.2f} {cp/(N*math.log(N)):>13.3f} {G/N:>8.4f} "
              f"{cpt/N:>10.2f}")
        q = quoted[N]
        print(f"{'v1 quotes':>10} {q[0]:>9.2f} {q[1]:>9.2f} {q[2]:>7.3f} "
              f"{'':>9} {q[3]:>13.3f} {q[4]:>8.4f} {q[5]:>10.2f}")

    print()
    S = 2.0 * 0.6601618158 * (1.0 + 1.0 / (5 - 2))
    print(f"(1) calibration target: S(N) for N = 2^a 5^b is "
          f"2*C_2*(1+1/3) = {S:.4f}; v1 quotes 1.7604")
    print()
    print("(2) internal relation CP_tuned = CP_2 - 2*gamma*G:")
    for N, cp, G in rows:
        q = quoted[N]
        pred = q[3] * math.log(N) - 2 * GAMMA * q[4]
        print(f"    N={N:>9}: from v1's own quoted CP_2 and G  -> "
              f"{pred:.2f}   v1 quotes CP_tuned/N = {q[5]:.2f}")
    print()
    print("(3) is CP_2/(N log N) settling?")
    ln = np.log([r[0] for r in rows])
    lv = np.log([r[1] / (r[0] * math.log(r[0])) for r in rows])
    d = np.polyfit(ln, lv, 1)[0]
    print(f"    fitted CP_2/(N log N) ~ N^{d:+.4f} over a factor 16 in N")
    print(f"    a constant would give 0. The column rises by "
          f"{100*(rows[-1][1]/(rows[-1][0]*math.log(rows[-1][0])))/(rows[0][1]/(rows[0][0]*math.log(rows[0][0])))-100:.1f}%"
          f" over that range.")
    print("DONE")


if __name__ == "__main__":
    main()
