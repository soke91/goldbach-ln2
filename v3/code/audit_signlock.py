# -*- coding: utf-8 -*-
r"""Does the sign of H(N;k) stay locked as N grows, or is it a small-N thing?

Supports {#rem:signlock}.

WHAT IS AT STAKE

Two measurements in this repository, taken for different reasons, found
the same thing: at N carrying several small odd primes the terms of
E_3 do not cancel at all.  {#rem:thetalawarith} reports
|sum a| / l1 = 1.00000000 at every swept theta' on a primorial-like
family, and {#rem:threshfam} reports that at
N = 2*3^2*5*7*11*13*17 all 513 terms are negative and |E_3| = B(N)
exactly.  Neither asked whether this survives N growing.

It matters because of where it sits.  Proposition onesided reaches its
threshold through |C(N)| <= A(N) N, and at N with many small prime
factors 1 - A(N) collapses -- the threshold is thinnest exactly there.
If the sign of H(N;k) is locked on that family for all large N, then
|E_3| = B(N) there, the one-sided condition has no cancellation left to
save it, and the failure is structural rather than numerical.  That
would be a no-go about the THRESHOLD rather than about the weight,
which is a different statement from the one this program already has.

But the repository has been burned by exactly this shape once already:
a crossing quoted as decisive turned out to be one arithmetic family's
number.  So the phenomenon is probed before it is believed, and the
probe is registered before it is run.

WHAT IS MEASURED

For squarefree k < N^theta' coprime to N, with
H(N;k) = sum_{m < N/k, (m,k)=1} Lambda(N-mk) mu(m):

  S1  the locked fraction: #{k : H(N;k) < 0} / #{k : H(N;k) != 0}, on a
      family whose radical contains 2,3,5,7,11,13 and on the swept
      family 2^a 5^b, at matched sizes.

  S2  the cancellation ratio rho = |sum_k (log k) H| / sum_k |(log k) H|
      on both families.  rho = 1 means no cancellation whatsoever.

  S3  the stated mechanism: the fraction of the m that carry a term
      (that is, N - mk prime) which are themselves prime.  The
      explanation on record is that those m "are almost all primes, so
      mu(m) = -1 dominates".

  S4  reported, not judged: |E_3|/B(N) on both families, which is 1
      exactly when the sign is locked.

FALSIFICATION, registered before the run

  S1  REFUTED if the locked fraction on the many-prime family is below
      1 at any N tested.  Then the sign is not locked and the two
      earlier measurements were about the sizes they were taken at.
  S2  REFUTED if rho on the many-prime family is below 1 - 1e-9 at any
      N, or if rho on the control family 2^a 5^b is above 1 - 1e-3 at
      every N.  The second half is the control: if the swept family
      also fails to cancel, this statistic is not measuring what it is
      being read as.
  S3  REFUTED if the prime fraction among contributing m falls below
      1/2 at the largest N while S1 still holds.  That would leave the
      phenomenon standing and its stated explanation dead, which is the
      outcome worth knowing.
  S4  reported, not judged.

  PREDICTION.  S1 and S2 hold across the tested range.  S3 is refuted:
  the contributing m run up to N/k, and the density of primes there is
  1/log(N/k), so "almost all primes" cannot survive N growing even if
  the sign lock does.  If that is how it comes out, the phenomenon is
  real and the explanation on record is wrong, and the mechanism has to
  be found before anything is claimed.

NULL.  None applies: deterministic sums over one integer sieve, no
sampling and no sign input.  The control is the second family, which is
run through identical code.
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

THETA = 0.56              # theta', as everywhere in this program
CORE = 30030              # 2*3*5*7*11*13
LOCK_TOL = 1e-9           # rho below 1 by more than this refutes S1/S2
CTRL_TOL = 1e-3           # the control family must cancel by at least this
PRIME_FLOOR = 0.5         # S3's threshold on the prime fraction

# many-prime family: radical contains 2,3,5,7,11,13 at every point
FAM_MANY = [CORE * j for j in (1, 2, 4, 8, 16, 32, 64, 128)]
# control: the family every sweep in this program used, over the same
# range of sizes as FAM_MANY so the comparison is like for like
FAM_CTRL = [200_000 * (2 ** j) for j in range(5)]


def sieve(n):
    spf = np.zeros(n + 1, dtype=np.int32)
    for i in range(2, int(n ** 0.5) + 1):
        if spf[i] == 0:
            spf[i * i::i] = np.where(spf[i * i::i] == 0, i, spf[i * i::i])
    lam = np.zeros(n + 1, dtype=np.float64)
    mu = np.zeros(n + 1, dtype=np.int8)
    isp = np.zeros(n + 1, dtype=bool)
    mu[1] = 1
    for v in range(2, n + 1):
        p = int(spf[v]) or v
        w = v // p
        mu[v] = 0 if w % p == 0 else -mu[w]
        t = v
        while t % p == 0:
            t //= p
        if t == 1:
            lam[v] = math.log(p)
            isp[v] = (v == p)
    return lam, mu, isp, spf


def factor_set(v, spf):
    out = set()
    while v > 1:
        p = int(spf[v]) or v
        out.add(p)
        while v % p == 0:
            v //= p
    return out


def scan(N, lam, mu, isp, spf):
    """H(N;k) for every admissible k, plus the prime fraction of the m."""
    PN = factor_set(N, spf)
    K = int(N ** THETA)
    Hs, ks = [], []
    m_tot = m_prime = 0
    for k in range(2, K):
        if mu[k] == 0:
            continue
        fk = factor_set(k, spf)
        if fk & PN:
            continue
        ms = np.arange(1, (N - 1) // k + 1, dtype=np.int64)
        for q in fk:
            ms = ms[ms % q != 0]
        vals = N - ms * k
        hit = lam[vals] != 0.0
        if hit.any():
            mm = ms[hit]
            m_tot += mm.size
            m_prime += int(isp[mm].sum())
        h = float((lam[vals] * mu[ms].astype(np.float64)).sum())
        ks.append(k)
        Hs.append(h)
    ks = np.array(ks, dtype=np.float64)
    Hs = np.array(Hs, dtype=np.float64)
    nz = Hs != 0.0
    locked = float((Hs[nz] < 0).sum()) / max(int(nz.sum()), 1)
    a = np.log(ks) * Hs
    l1 = float(np.abs(a).sum())
    rho = abs(float(a.sum())) / l1 if l1 > 0 else float("nan")
    E3 = float(a.sum())
    B = l1
    pf = m_prime / m_tot if m_tot else float("nan")
    return locked, rho, abs(E3) / B if B > 0 else float("nan"), pf, int(nz.sum())


def main():
    lines = []

    def say(s=""):
        print(s)
        lines.append(s)

    say("STATISTIC: the locked fraction #{k: H(N;k)<0}/#{k: H(N;k)!=0};")
    say("           the cancellation ratio rho = |sum_k (log k)H| /")
    say("           sum_k |(log k)H|; |E_3|/B(N); and the fraction of the")
    say("           m carrying a term that are themselves prime.")
    say("FIELD: two families at matched sizes -- N = 30030*j with")
    say("       j = 1,2,4,...,128, whose radical contains 2,3,5,7,11,13")
    say("       at every point, and N = 2e5*2^j, j = 0..7, the family")
    say("       every sweep in this program used; k squarefree,")
    say("       2 <= k < N^THETA, (k,N) = 1; m < N/k with (m,k) = 1;")
    say("       Lambda, mu, and the primality flag from one sieve.")
    say("CONSTANTS: THETA = %.2f, CORE = %d, LOCK_TOL = %.0e,"
        % (THETA, CORE, LOCK_TOL))
    say("           CTRL_TOL = %.0e, PRIME_FLOOR = %.2f"
        % (CTRL_TOL, PRIME_FLOOR))
    say("NULL: none applies -- deterministic sums over one sieve, no")
    say("      sampling and no sign input. The control is the second")
    say("      family, run through identical code.")
    say("DENOM: rho and |E_3|/B are ratios of the same sum's signed and")
    say("      absolute totals; the prime fraction is over the m that")
    say("      carry a term.")
    say()
    say(__doc__.strip())
    say()
    say("=" * 74)
    NMAX = max(max(FAM_MANY), max(FAM_CTRL))
    say("sieving to %d ..." % NMAX)
    lam, mu, isp, spf = sieve(NMAX)

    out = {}
    for tag, fam in (("many", FAM_MANY), ("ctrl", FAM_CTRL)):
        say()
        say("  family %s" % ("2,3,5,7,11,13 | N" if tag == "many"
                             else "N = 2^a 5^b (the swept one)"))
        hdr = ("  %-11s %-7s %-10s %-13s %-10s %-10s"
               % ("N", "#k", "locked", "rho", "|E3|/B", "prime frac"))
        say(hdr)
        say("  " + "-" * (len(hdr) - 2))
        rows = []
        for N in fam:
            locked, rho, eb, pf, nk = scan(N, lam, mu, isp, spf)
            rows.append((N, nk, locked, rho, eb, pf))
            say("  %-11d %-7d %-10.6f %-13.10f %-10.6f %-10.6f"
                % (N, nk, locked, rho, eb, pf))
        out[tag] = rows

    many, ctrl = out["many"], out["ctrl"]
    ok1 = all(r[2] >= 1.0 for r in many)
    ok2 = (all(r[3] >= 1.0 - LOCK_TOL for r in many)
           and any(r[3] <= 1.0 - CTRL_TOL for r in ctrl))
    pf_last = many[-1][5]
    ok3 = not (ok1 and pf_last < PRIME_FLOOR)

    say()
    say("  the control family's rho: %s"
        % ", ".join("%.6f" % r[3] for r in ctrl))
    say("  the many-prime family's prime fraction, first to last:")
    say("    %s" % ", ".join("%.6f" % r[5] for r in many))
    say()
    say("S1  the sign of H(N;k) is locked on the many-prime family")
    say("    S1 %s" % ("hold" if ok1 else "REFUTED"))
    say("S2  no cancellation there, and the control family does cancel")
    say("    S2 %s" % ("hold" if ok2 else "REFUTED"))
    say("S3  the stated mechanism -- the contributing m are almost all prime")
    say("    S3 %s" % ("hold" if ok3 else "REFUTED"))
    say("S4  |E_3|/B(N) is in the tables above, reported not judged.")
    say()
    say("=" * 74)
    say("S1 %s  S2 %s  S3 %s"
        % tuple("hold" if v else "REFUTED" for v in (ok1, ok2, ok3)))

    io.open(os.path.join(RES, "audit_signlock.txt"), "w",
            encoding="utf-8", newline="\n").write("\n".join(lines) + "\n")
    return 0 if (ok1 and ok2) else 1


if __name__ == "__main__":
    raise SystemExit(main())
