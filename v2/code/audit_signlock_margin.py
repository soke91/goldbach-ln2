# -*- coding: utf-8 -*-
r"""How close is the sign lock to breaking, and where would it break?

Supports {#rem:signlockmargin}.

WHAT IS AT STAKE

{#rem:signlock} established that on a family whose radical contains
2,3,5,7,11,13 every H(N;k) is negative across a factor 128 in N, and
that the mechanism is forced rather than statistical: primality of
N - mk requires m to be coprime to rad(N), and among rough m the ones
with a single prime factor outnumber those with two.  It also noted
that this bias thins like 1/log, so the lock should be a transient.

"Should be" is not a measurement.  Pushing N until the lock breaks is
the direct test and it is expensive; the margin is the cheap one.  For
each k split H(N;k) by the sign of mu(m):

    P(k) = sum_{m: mu(m)=+1} Lambda(N-mk),
    M(k) = sum_{m: mu(m)=-1} Lambda(N-mk),
    r_k  = P(k) / M(k),

so H(N;k) < 0 exactly when r_k < 1.  The lock is the statement
max_k r_k < 1.  How max_k r_k moves with N says whether the lock is
about to go, and where.

WHAT IS MEASURED

  T1  the lock, restated as a margin: max_k r_k on the many-prime
      family at ten N, now reaching 1.5e7.

  T2  its trend: max_k r_k against log N, and the N at which the fitted
      line reaches 1.

  T3  the control: max_k r_k on 2^a 5^b, where the lock does not hold.

  T4  reported, not judged: the median r_k, which says whether the
      typical term is near the edge or far from it.

FALSIFICATION, registered before the run

  T1  REFUTED if max_k r_k >= 1 at any N of the many-prime family.
      That is the lock breaking, and it would mean {#rem:signlock}
      measured a transient that ends inside the computable range.
  T2  REFUTED if max_k r_k at the largest N is not above its value at
      the smallest.  Then the margin is not closing, the 1/log
      reasoning does not describe what the maximum does, and the lock
      is not shown to be a transient by anything measured here.
  T3  REFUTED if max_k r_k < 1 on the control family at any N -- the
      control is supposed to have terms of both signs, and if it does
      not, this statistic is not separating the families.
  T4  reported, not judged.

  PREDICTION.  T1 and T3 hold.  T2 holds, and the extrapolated crossing
  sits far above anything computable -- which would leave the lock a
  transient in principle and a fact in practice, and would say so with
  a number instead of an argument.

NULL.  None applies: deterministic sums over one sieve, no sampling and
no sign input.  The control is the second family, run through identical
code.
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

THETA = 0.56
CORE = 30030                       # 2*3*5*7*11*13
FAM_MANY = [CORE * (2 ** j) for j in range(10)]        # 3.0e4 .. 1.5e7
FAM_CTRL = [200_000 * (2 ** j) for j in range(6)]      # 2.0e5 .. 6.4e6


def sieve(n):
    """mu, Lambda and the prime flag to n, vectorised."""
    mu = np.ones(n + 1, dtype=np.int8)
    lam = np.zeros(n + 1, dtype=np.float32)
    comp = np.zeros(n + 1, dtype=bool)
    comp[:2] = True
    mu[0] = 0
    for p in range(2, n + 1):
        if comp[p]:
            continue
        comp[p * p::p] = True
        mu[p::p] = -mu[p::p]
        sq = p * p
        if sq <= n:
            mu[sq::sq] = 0
        lp = math.log(p)
        q = p
        while q <= n:
            lam[q] = lp
            if q > n // p:
                break
            q *= p
    isp = ~comp
    isp[:2] = False
    return mu, lam, isp


def factor_set(v):
    out, d = set(), 2
    while d * d <= v:
        if v % d == 0:
            out.add(d)
            while v % d == 0:
                v //= d
        d += 1 if d == 2 else 2
    if v > 1:
        out.add(v)
    return out


def margins(N, mu, lam):
    """r_k for every admissible k, and the count of nonzero terms."""
    PN = factor_set(N)
    K = int(N ** THETA)
    rs = []
    for k in range(2, K):
        if mu[k] == 0:
            continue
        fk = factor_set(k)
        if fk & PN:
            continue
        ms = np.arange(1, (N - 1) // k + 1, dtype=np.int64)
        for q in fk:
            ms = ms[ms % q != 0]
        w = lam[N - ms * k]
        nz = w != 0
        if not nz.any():
            continue
        s = mu[ms[nz]]
        ww = w[nz].astype(np.float64)
        P = float(ww[s > 0].sum())
        M = float(ww[s < 0].sum())
        if M <= 0.0:
            rs.append(float("inf") if P > 0 else float("nan"))
        else:
            rs.append(P / M)
    return np.array(rs, dtype=np.float64)


def main():
    lines = []

    def say(s=""):
        print(s)
        lines.append(s)

    say("STATISTIC: r_k = P(k)/M(k) with P and M the Lambda-mass of the")
    say("           mu(m) = +1 and mu(m) = -1 terms of H(N;k); its")
    say("           maximum and median over the admissible k; and the")
    say("           fitted crossing of max_k r_k with 1 against log N.")
    say("FIELD: N = 30030*2^j, j = 0..9 (radical contains 2,3,5,7,11,13")
    say("       at every point), and the control N = 2e5*2^j, j = 0..5;")
    say("       k squarefree with 2 <= k < N^THETA and (k,N) = 1; m < N/k")
    say("       with (m,k) = 1; mu, Lambda and the prime flag from one")
    say("       vectorised sieve to max N.")
    say("CONSTANTS: THETA = %.2f, CORE = %d" % (THETA, CORE))
    say("NULL: none applies -- deterministic sums over one sieve, no")
    say("      sampling and no sign input. The control is the second")
    say("      family, run through identical code.")
    say("DENOM: r_k is the positive Lambda-mass over the negative one;")
    say("      the lock is max_k r_k < 1.")
    say()
    say(__doc__.strip())
    say()
    say("=" * 74)
    NMAX = max(max(FAM_MANY), max(FAM_CTRL))
    say("sieving to %d ..." % NMAX)
    mu, lam, _ = sieve(NMAX)

    out = {}
    for tag, fam in (("many", FAM_MANY), ("ctrl", FAM_CTRL)):
        say()
        say("  family %s" % ("radical contains 2,3,5,7,11,13"
                             if tag == "many" else "N = 2^a 5^b"))
        hdr = "  %-11s %-7s %-12s %-12s %-10s" % (
            "N", "#k", "max r_k", "median r_k", "locked")
        say(hdr)
        say("  " + "-" * (len(hdr) - 2))
        rows = []
        for N in fam:
            r = margins(N, mu, lam)
            r = r[np.isfinite(r)]
            mx = float(r.max()) if r.size else float("nan")
            md = float(np.median(r)) if r.size else float("nan")
            rows.append((N, r.size, mx, md))
            say("  %-11d %-7d %-12.6f %-12.6f %-10s"
                % (N, r.size, mx, md, "yes" if mx < 1.0 else "NO"))
        out[tag] = rows

    many, ctrl = out["many"], out["ctrl"]
    ok1 = all(r[2] < 1.0 for r in many)
    ok2 = many[-1][2] > many[0][2]
    ok3 = all(r[2] >= 1.0 for r in ctrl)

    xs = np.log(np.array([r[0] for r in many], dtype=np.float64))
    ys = np.array([r[2] for r in many], dtype=np.float64)
    a, b = np.polyfit(xs, ys, 1)
    say()
    say("  max_k r_k against log N: slope %+.6f, intercept %+.6f" % (a, b))
    if a > 0:
        cross = math.exp((1.0 - b) / a)
        say("  the fitted line reaches 1 at log N = %.3f, i.e. N = 10^%.2f"
            % ((1.0 - b) / a, math.log10(cross)))
    else:
        say("  the fitted line does not rise; no crossing to report.")
    say("  a straight line in log N is the simplest reading and not a")
    say("  justified model; it is printed so the reading is checkable.")

    say()
    say("T1  the lock holds at every N of the many-prime family")
    say("    T1 %s" % ("hold" if ok1 else "REFUTED"))
    say("T2  the margin closes as N grows")
    say("    T2 %s" % ("hold" if ok2 else "REFUTED"))
    say("T3  the control family is not locked")
    say("    T3 %s" % ("hold" if ok3 else "REFUTED"))
    say("T4  median r_k is in the table, reported not judged.")
    say()
    say("=" * 74)
    say("T1 %s  T2 %s  T3 %s"
        % tuple("hold" if v else "REFUTED" for v in (ok1, ok2, ok3)))

    io.open(os.path.join(RES, "audit_signlock_margin.txt"), "w",
            encoding="utf-8", newline="\n").write("\n".join(lines) + "\n")
    return 0 if (ok1 and ok3) else 1


if __name__ == "__main__":
    raise SystemExit(main())
