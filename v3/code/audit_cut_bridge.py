# -*- coding: utf-8 -*-
r"""Does the fixed-cut theorem reach the moving-cut object by summation?

Supports {#rem:cutbridge}.

WHAT IS AT STAKE

Two formulations of the same double sum are in play.  The fixed cut
gives every k the same endpoint t, and the theorem of this note bounds

    T1(t; K) = sum_{k<K, (k,N)=1} mu(k) Emu(t;k)     uniformly in t.

The corrected formulation gives k its own endpoint Y_k = ceil(N-alpha k)-1
and asks about

    T1* = sum_{k<K, (k,N)=1} mu(k) Emu(Y_k; k).

Neither is formally a consequence of the other: the supremum in the
first sits OUTSIDE the k-sum, so using it on the second requires
absolute values inside, which destroys the cancellation in k that the
whole statement is about.  There is, however, an exact bridge.  Since
Y_k decreases in k,

    T1* = T1(Y_J; K) + sum_{j<J} [ T1(Y_j; j+1) - T1(Y_{j+1}; j+1) ],

J the largest index -- a telescoping identity in the two-parameter
family T1(t; K').  So the fixed-cut statement reaches the moving-cut
object if (a) it is uniform in the outer truncation K' as well as in t,
and (b) the telescoping sum can be bounded term by term.  This run
measures whether either is true, and which of the two objects is the
larger.

The reason to ask: if the fixed cut is the harder object and the bridge
is cheap, then the fixed-cut statement is the stronger one and the
corrected formulation is a corollary of it.  If the bridge is not
cheap, the direct re-derivation under the moving cut is not a
convenience but a necessity, and the two statements are separate
results.

WHAT IS MEASURED

For N even, alpha = N^theta, K = (N-1)/alpha, a_n = Lambda(n) mu(N-n),
C(t) = sum_{n<=t} a_n, and

    d_{K'}(n) = sum_{k | N-n, (k,N)=1, k < K'} mu(k),
    W_{K'}    = sum_{k<K', (k,N)=1} mu(k)/phi(k),
    T1(t; K') = sum_{n<=t} a_n d_{K'}(n) - W_{K'} C(t),

which is the same T1 regrouped by divisor count instead of by residue
class -- a different decomposition, computed here from scratch.

  Y0  control.  |T1*|/N against the value an independently written
      script produced for the same N and alpha.

  Y1  which object is larger: |T1*| against sup_{1<=t<N} |T1(t;K)|,
      the supremum taken over every t, not a grid.

  Y2  what the outer truncation costs inside the admissible band:
      max over K' in [N^(1/2), K] and over t of |T1(t;K')|, against
      sup_t |T1(t;K)|.  Values of K' below N^(1/2) are printed but not
      judged -- there the completion has no proof at any N, and a
      finite computation cannot see that.

  Y3  the price of the bridge: sum_j |T1(Y_j;j+1) - T1(Y_{j+1};j+1)|
      against N.  The telescoping identity is itself checked first:
      the decomposition must reproduce T1*.

  Y4  the share carried by the range where the completion has no
      proof: T1* split at R_n = ceil((N-n)/alpha)-1 = N^(1/2), the
      cofactor threshold, into bulk (R_n >= N^(1/2)) and boundary.

  Y5  the rate at which the bridge's cost falls: the least-squares
      slope of log(cost/N) against log N.  See the disclosure below --
      this check was added after Y3's threshold was seen to be the
      wrong question.

DISCLOSURE

A code smoke test was run at N = 20000 and 40000 -- outside the field
registered here -- to check that the program runs.  Its numbers were
visible, and they showed that Y3's registered threshold asks the wrong
question: the target is not N but N(log N)^(-A) for every A, so a cost
below N settles nothing.  Y3's rule is left exactly as registered and
whatever it returns is reported as it stands.  Y5 was added afterward,
and its threshold is a prediction about the mechanism rather than about
the size: if the cost falls because each increment enjoys square-root
cancellation over its interval of length alpha, then summing K = N/alpha
of them gives cost/N ~ N^(-theta/2), a slope of -theta/2 = -0.22.  The
five N below are a different field from the two the smoke test used, so
the exponent is being tested out of sample.

FALSIFICATION, registered before the run

  Y0  REFUTED if |T1*|/N differs from the earlier value by more than
      1e-8, which is the precision that value was printed to.  Then the
      two implementations disagree and nothing below is interpretable.
  Y1  REFUTED if |T1*| exceeds sup_t |T1(t;K)| at any N.  Then the
      moving-cut object is the larger one, the fixed-cut statement is
      weaker than what the corrected formulation needs, and no amount
      of uniformity recovers it.
  Y2  REFUTED if the two-parameter maximum over the admissible band
      exceeds 3 times sup_t |T1(t;K)|.  Then uniformity in the outer
      truncation is not free even where the mechanism is supposed to
      work.
  Y3  REFUTED if sum_j |increment_j| / N falls below 1 at any N, or if
      the telescoping identity fails to reproduce T1* to 1e-12
      relative to N.  The first would mean the bridge is cheap -- the
      fixed-cut family reaches the moving-cut object term by term, and
      re-deriving under the moving cut was unnecessary.  The second
      would mean the identity above is misstated.
  Y4  REFUTED if |T1*_boundary|/N at the largest N is not below its
      value at the smallest N.
  Y5  REFUTED if the fitted slope differs from -theta/2 by more than
      0.10.  Then the cost is not falling at the rate square-root
      cancellation over an interval of length alpha would give, and
      that explanation of Y3's outcome is wrong.

  PREDICTION.  Y0, Y1, Y2, Y4, Y5 hold.  Y3 was predicted to hold when
  it was registered; the smoke test disclosed above makes that
  prediction obsolete rather than confirmed, and it is left standing so
  the record shows what was asked before the numbers were seen.

  Expect the admissible band R_n in [N^(1/2), N^(theta')] to be thin at
  these N, so Y4's decrease should be slow.

  WHAT A CHEAP BRIDGE WOULD AND WOULD NOT SETTLE.  If Y5 holds, the
  cost falls like a power, far below the target N(log N)^(-A) -- and
  that is exactly the point at which the computation stops being
  evidence.  A power saving here IS square-root cancellation for
  Lambda over an interval of length alpha = N^(1-theta'), which is a
  short-interval prime estimate no unconditional argument supplies.
  The numerics would say the bridge is cheap and no proof would follow.
  That is the reason the moving cut was re-derived directly rather than
  bridged, and Y5 is the check that says so in numbers.

NULL.  None applies: every quantity is a deterministic finite sum with
no sampling and no sign input.  The controls are Y0, which re-derives a
quantity from an independently written script, and the decomposition
itself, which groups the same sum by divisor count rather than by
residue class.
"""
import io
import math
import os
import sys

import numpy as np

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(os.path.dirname(HERE), "results")

NS = (200_000, 400_000, 800_000, 1_600_000, 3_200_000)
THETA = 0.44           # alpha = N^theta, as in audit_moving_switch
TOL_ID = 1e-12         # telescoping identity, relative to N
TOL_CTL = 1e-8         # control: the prior value is printed to 8 decimals
Y2_CAP = 3.0           # allowed factor for the two-parameter maximum
N_CHECKPOINTS = 24     # K' values sampled, log-spaced over [1, K]

# |T1*|/N as produced by code/audit_moving_switch.py at the same N and
# alpha, by a different decomposition (residue classes, k-loop).
PRIOR = {
    200_000: 0.02284282,
    400_000: 0.00655656,
    800_000: 0.01205837,
    1_600_000: 0.00871321,
    3_200_000: 0.00612307,
}


def sieve(n):
    """Lambda, mu, and the smallest-prime-factor table, to n."""
    spf = np.zeros(n + 1, dtype=np.int64)
    for i in range(2, int(n ** 0.5) + 1):
        if spf[i] == 0:
            spf[i * i::i] = np.where(spf[i * i::i] == 0, i, spf[i * i::i])
    lam = np.zeros(n + 1, dtype=np.float64)
    mu = np.zeros(n + 1, dtype=np.int64)
    mu[1] = 1
    for v in range(2, n + 1):
        p = int(spf[v]) or v
        w = v // p
        mu[v] = 0 if w % p == 0 else -mu[w]
        t, e = v, 0
        while t % p == 0:
            t //= p
            e += 1
        lam[v] = math.log(p) if t == 1 else 0.0
    return lam, mu, spf


def factor_set(v, spf):
    out = set()
    while v > 1:
        p = int(spf[v]) or v
        out.add(p)
        while v % p == 0:
            v //= p
    return out


def main():
    lines = []

    def say(s=""):
        print(s)
        lines.append(s)

    say("STATISTIC: T1(t;K') = sum_{n<=t} a_n d_{K'}(n) - W_{K'} C(t),")
    say("           a_n = Lambda(n) mu(N-n), regrouped by divisor count;")
    say("           its exact supremum over t, its two-parameter maximum")
    say("           over K', the moving-cut value T1*, the telescoping")
    say("           cost of passing from one to the other, and the split")
    say("           of T1* at the cofactor threshold R_n = N^(1/2).")
    say("FIELD: N even in %s; alpha = N^theta; k < K = (N-1)/alpha with"
        % ",".join(str(v) for v in NS))
    say("       (k,N) = 1 and k squarefree; n < N; R_n = ceil((N-n)/alpha)-1;")
    say("       Y_k = ceil(N - alpha k) - 1. Lambda, mu from one sieve.")
    say("CONSTANTS: THETA = %.2f, TOL_ID = %.0e, TOL_CTL = %.0e,"
        % (THETA, TOL_ID, TOL_CTL))
    say("           Y2_CAP = %.1f, N_CHECKPOINTS = %d, PRIOR from"
        % (Y2_CAP, N_CHECKPOINTS))
    say("           code/audit_moving_switch.py (8 decimals as printed).")
    say("NULL: none applies -- deterministic finite sums, no sampling, no")
    say("      sign input. Controls: Y0 against an independently written")
    say("      script; the decomposition groups by divisor count, not by")
    say("      residue class.")
    say("DENOM: every ratio printed as .../N is divided by N; the Y1 and")
    say("      Y2 ratios are divided by sup_t |T1(t;K)| at the same N.")
    say()
    say(__doc__.strip())
    say()
    say("=" * 72)
    say("sieving to %d ..." % max(NS))
    lam, mu, spf = sieve(max(NS))
    say()

    hdr = ("  %-9s %-6s %-6s %-11s %-11s %-9s %-9s"
           % ("N", "alpha", "K", "sup|T1|/N", "|T1*|/N", "Y1 ratio", "Y2 ratio"))
    say(hdr)
    say("  " + "-" * (len(hdr) - 2))

    ok0 = ok1 = ok2 = ok3 = True
    rows, bdry_share, band = [], [], []
    for N in NS:
        alpha = N ** THETA
        K = (N - 1) / alpha
        jmax = int(math.ceil(K)) - 1
        PN = factor_set(N, spf)

        # phi on the outer range, and the admissible k
        phi = np.arange(jmax + 2, dtype=np.float64)
        for p in range(2, jmax + 2):
            if phi[p] == p:                     # p prime
                phi[p::p] *= (1.0 - 1.0 / p)
        phi[1] = 1.0

        nn = np.arange(0, N, dtype=np.int64)
        a = lam[:N] * mu[N - nn].astype(np.float64)
        a[0] = 0.0
        C = np.cumsum(a)

        good = [k for k in range(1, jmax + 1)
                if mu[k] != 0 and not (factor_set(k, spf) & PN)]
        goodset = set(good)

        # checkpoints for Y2, log-spaced in K'
        cps = sorted(set(
            int(round(x)) for x in
            np.exp(np.linspace(0.0, math.log(max(jmax, 2)), N_CHECKPOINTS))
            if 1 <= x <= jmax))
        cpset = set(cps)

        Y = np.zeros(jmax + 2, dtype=np.int64)
        for j in range(1, jmax + 1):
            Y[j] = math.ceil(N - alpha * j) - 1

        # one increasing sweep in k: builds d, snapshots the two-parameter
        # maximum at the checkpoints, and reads the telescoping increments
        d = np.zeros(N, dtype=np.int64)
        W = 0.0
        inc = []
        band_max = 0.0
        small_max = 0.0
        root = math.sqrt(N)
        for j in range(1, jmax + 1):
            if j in goodset:
                m = np.arange(1, (N - 1) // j + 1, dtype=np.int64)
                d[N - j * m] += int(mu[j])
                W += float(mu[j]) / phi[j]
            if j <= jmax - 1:
                lo, hi = int(Y[j + 1]), int(Y[j])
                if hi > lo >= 0:
                    seg = math.fsum((a[lo + 1:hi + 1]
                                     * d[lo + 1:hi + 1]).tolist())
                    inc.append(seg - W * (float(C[hi]) - float(C[lo])))
            if j in cpset:
                S = np.cumsum(a * d) - W * C
                v = float(np.max(np.abs(S)))
                if j + 1 >= root:
                    band_max = max(band_max, v)
                else:
                    small_max = max(small_max, v)

        S = np.cumsum(a * d) - W * C
        supT1 = float(np.max(np.abs(S)))
        band_max = max(band_max, supT1)
        lead = float(S[int(Y[jmax])]) if Y[jmax] >= 1 else 0.0

        # the moving-cut object, built from the cofactor cut
        ds = np.zeros(N, dtype=np.int64)
        for k in good:
            m0 = int(math.floor(alpha)) + 1
            if m0 * k >= N:
                continue
            m = np.arange(m0, (N - 1) // k + 1, dtype=np.int64)
            ds[N - k * m] += int(mu[k])
        Wpref = np.zeros(jmax + 2, dtype=np.float64)
        run = 0.0
        for k in range(1, jmax + 1):
            if k in goodset:
                run += float(mu[k]) / phi[k]
            Wpref[k] = run
        Rn = np.ceil((N - nn[1:]) / alpha).astype(np.int64) - 1
        np.clip(Rn, 0, jmax, out=Rn)
        g = a[1:] * (ds[1:].astype(np.float64) - Wpref[Rn])
        T1s = math.fsum(g.tolist())

        # Y3: does the telescoping decomposition reproduce it
        tele = lead + math.fsum(inc)
        r_id = abs(tele - T1s) / N
        cost = math.fsum(abs(v) for v in inc) / N

        # Y4: split at the cofactor threshold
        bulk = math.fsum(g[Rn >= root].tolist())
        bdry = math.fsum(g[Rn < root].tolist())
        bdry_share.append(abs(bdry) / N)
        band.append((root, float(jmax)))

        r0 = abs(abs(T1s) / N - PRIOR[N])
        y1 = abs(T1s) / supT1
        y2 = band_max / supT1
        ok0 &= r0 <= TOL_CTL
        ok1 &= abs(T1s) <= supT1
        ok2 &= y2 <= Y2_CAP
        ok3 &= (r_id <= TOL_ID) and (cost >= 1.0)

        say("  %-9d %-6.0f %-6d %-11.6f %-11.6f %-9.4f %-9.4f"
            % (N, alpha, jmax, supT1 / N, abs(T1s) / N, y1, y2))
        rows.append((N, r0, r_id, cost, small_max / supT1,
                     abs(bulk) / N, abs(bdry) / N))

    say()
    say("  %-9s %-11s %-11s %-11s %-11s %-11s %-11s"
        % ("N", "Y0 diff", "Y3 id", "Y3 cost/N", "small K'", "bulk/N",
           "bdry/N"))
    say("  " + "-" * 82)
    for N, r0, r_id, cost, sm, bu, bd in rows:
        say("  %-9d %-11.3e %-11.3e %-11.2f %-11.4f %-11.6f %-11.6f"
            % (N, r0, r_id, cost, sm, bu, bd))

    say()
    say("  admissible band for the outer truncation, K' in [N^(1/2), K]:")
    for N, (lo, hi) in zip(NS, band):
        say("    N = %-9d  [%.0f, %.0f]   width in log: %.3f"
            % (N, lo, hi, math.log(hi / lo) if hi > lo else float("nan")))

    ok4 = bdry_share[-1] < bdry_share[0]

    xs = np.log(np.array([r[0] for r in rows], dtype=np.float64))
    ys = np.log(np.array([r[3] for r in rows], dtype=np.float64))
    slope, icpt = np.polyfit(xs, ys, 1)
    ok5 = abs(float(slope) - (-THETA / 2.0)) <= 0.10

    say()
    say("  Y5  cost/N against N, least squares in logs:")
    say("      fitted slope   %+.4f" % slope)
    say("      square-root cancellation over alpha predicts  %+.4f"
        % (-THETA / 2.0))
    say("      residuals:")
    for x, y, r in zip(xs, ys, rows):
        say("        N = %-9d  log(cost/N) = %+.4f   fitted %+.4f"
            % (r[0], y, slope * x + icpt))

    say()
    say("Y0  the moving-cut value reproduces the earlier one")
    say("    Y0 %s" % ("hold" if ok0 else "REFUTED"))
    say("Y1  the fixed-cut supremum is the larger object")
    say("    Y1 %s" % ("hold" if ok1 else "REFUTED"))
    say("Y2  uniformity in the outer truncation is cheap in the band")
    say("    Y2 %s" % ("hold" if ok2 else "REFUTED"))
    say("Y3  the telescoping bridge is exact, and costs more than the target")
    say("    Y3 %s" % ("hold" if ok3 else "REFUTED"))
    say("Y4  the unproved boundary carries a shrinking share")
    say("    Y4 %s" % ("hold" if ok4 else "REFUTED"))
    say("Y5  the cost falls at the rate short-interval cancellation gives")
    say("    Y5 %s" % ("hold" if ok5 else "REFUTED"))
    say()
    say("  Y3 does not enter the exit status. Its threshold was shown to")
    say("  ask the wrong question before this field was run (see")
    say("  DISCLOSURE), so passing or failing it decides nothing. The")
    say("  rule and its verdict are printed unchanged; Y5 carries the")
    say("  question Y3 meant to ask.")
    say()
    say("  K' below N^(1/2) is printed, not judged: there the completion")
    say("  leaves a cofactor that is not short, and no computation at these")
    say("  N can show that. The column is there so the claim is checkable")
    say("  in the only sense it can be -- as a size, not as a proof.")
    say()
    say("=" * 72)
    say("Y0 %s  Y1 %s  Y2 %s  Y3 %s  Y4 %s  Y5 %s"
        % tuple("hold" if v else "REFUTED"
                for v in (ok0, ok1, ok2, ok3, ok4, ok5)))

    io.open(os.path.join(RES, "audit_cut_bridge.txt"), "w",
            encoding="utf-8", newline="\n").write("\n".join(lines) + "\n")
    return 0 if (ok0 and ok1 and ok2 and ok4 and ok5) else 1


if __name__ == "__main__":
    raise SystemExit(main())
