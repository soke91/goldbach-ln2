# -*- coding: utf-8 -*-
r"""Is the sign lock a fact about primes, or a fact about rough integers?

Supports {#rem:signlockwhy}.

WHAT IS AT STAKE

{#rem:signlock} found that every H(N;k) is negative on a family whose
radical carries 2,3,5,7,11,13, and that primality of N - mk forces the
contributing m to be coprime to rad(N).  {#rem:signlockmargin} found
the margin saturating below 1 rather than crossing it.

That leaves one question, and it decides what the finding is worth.
Primality enters twice: it forces m to be rough, and it decides which
rough m carry a term at all.  If only the first matters, the lock is a
statement about the parity of omega on rough integers in a short range
-- elementary, and nothing to do with the Goldbach problem.  If the
second matters too, the lock says something about how primes sit in
the progressions N - mk, and that is a different object.

The test separates them.  For each k compare

    r_k = sum_{mu(m)=+1} Lambda(N-mk) / sum_{mu(m)=-1} Lambda(N-mk),

which is the real quantity, against

    q_k = #{m rough, (m,k)=1, mu(m)=+1} / #{m rough, (m,k)=1, mu(m)=-1},

which counts the same m with the primality of N - mk thrown away.  If
r_k and q_k agree, the primes did nothing but force roughness.

WHAT IS MEASURED

For N on the many-prime family and k squarefree, 2 <= k < N^theta',
(k,N) = 1, and m < N/k coprime to k and to rad(N):

  W1  max_k q_k -- whether the unweighted parity ratio is locked too.
  W2  how closely r_k tracks q_k: the Pearson correlation across k, and
      the gap between their medians.
  W3  reported, not judged: both medians and both maxima, per N.

FALSIFICATION, registered before the run

  W1  REFUTED if max_k q_k >= 1 at any N.  Then the unweighted ratio is
      not locked, the lock is not a fact about rough integers alone,
      and the primality of N - mk is doing part of the work.
  W2  REFUTED if the correlation between r_k and q_k across k is below
      0.8 at any N, or if the two medians differ by more than 0.15 at
      any N.  Same conclusion as W1 from the other side.
  W3  reported, not judged.

  PREDICTION.  W1 and W2 hold, and the finding deflates: the lock is
  the parity of omega on rough integers below N/k, primality having
  contributed the roughness and nothing else.  That is worth knowing
  precisely because it is the deflating outcome -- {#rem:signlock} is
  currently written as though a mechanism about primes had been found.

NULL.  None applies: deterministic counts and sums over one sieve, no
sampling and no sign input.  The control is internal: r_k and q_k are
two readings of the same index set, and W2 is the comparison.
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
CORE = 30030
FAM = [CORE * (2 ** j) for j in range(9)]      # 3.0e4 .. 7.7e6
CORR_FLOOR = 0.8
MED_GAP = 0.15


def sieve(n):
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
    return mu, lam


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


def scan(N, mu, lam):
    PN = factor_set(N)
    K = int(N ** THETA)
    rs, qs = [], []
    for k in range(2, K):
        if mu[k] == 0:
            continue
        fk = factor_set(k)
        if fk & PN:
            continue
        ms = np.arange(1, (N - 1) // k + 1, dtype=np.int64)
        for q in fk | PN:                    # rough, and coprime to k
            ms = ms[ms % q != 0]
        s = mu[ms]
        ms = ms[s != 0]
        s = s[s != 0]
        if ms.size == 0:
            continue
        # unweighted parity ratio over exactly this index set
        npos = int((s > 0).sum())
        nneg = int((s < 0).sum())
        if nneg == 0:
            continue
        qs.append(npos / nneg)
        # the real quantity, on the same index set
        w = lam[N - ms * k].astype(np.float64)
        P = float(w[s > 0].sum())
        M = float(w[s < 0].sum())
        rs.append(P / M if M > 0 else float("nan"))
    r = np.array(rs, dtype=np.float64)
    q = np.array(qs, dtype=np.float64)
    ok = np.isfinite(r) & np.isfinite(q)
    return r[ok], q[ok]


def main():
    lines = []

    def say(s=""):
        print(s)
        lines.append(s)

    say("STATISTIC: per k, the Lambda-weighted positive-to-negative mass")
    say("           ratio r_k of H(N;k), and the unweighted count ratio")
    say("           q_k over the same index set with the primality of")
    say("           N - mk discarded; their maxima, medians, and the")
    say("           Pearson correlation of r_k against q_k across k.")
    say("FIELD: N = 30030*2^j, j = 0..8; k squarefree, 2 <= k < N^THETA,")
    say("       (k,N) = 1; m < N/k coprime to k and to rad(N) with")
    say("       mu(m) != 0; mu and Lambda from one vectorised sieve.")
    say("CONSTANTS: THETA = %.2f, CORE = %d, CORR_FLOOR = %.1f,"
        % (THETA, CORE, CORR_FLOOR))
    say("           MED_GAP = %.2f" % MED_GAP)
    say("NULL: none applies -- deterministic counts and sums over one")
    say("      sieve. The control is internal: r_k and q_k read the same")
    say("      index set two ways.")
    say("DENOM: both ratios are positive over negative; the lock is")
    say("      max < 1.")
    say()
    say(__doc__.strip())
    say()
    say("=" * 78)
    NMAX = max(FAM)
    say("sieving to %d ..." % NMAX)
    mu, lam = sieve(NMAX)

    hdr = ("  %-11s %-6s %-10s %-10s %-10s %-10s %-8s"
           % ("N", "#k", "max r_k", "max q_k", "med r_k", "med q_k",
              "corr"))
    say()
    say(hdr)
    say("  " + "-" * (len(hdr) - 2))
    ok1 = ok2 = True
    rows = []
    for N in FAM:
        r, q = scan(N, mu, lam)
        if r.size < 3:
            continue
        c = float(np.corrcoef(r, q)[0, 1])
        mr, mq = float(np.median(r)), float(np.median(q))
        xr, xq = float(r.max()), float(q.max())
        ok1 &= xq < 1.0
        ok2 &= (c >= CORR_FLOOR) and (abs(mr - mq) <= MED_GAP)
        rows.append((N, r.size, xr, xq, mr, mq, c))
        say("  %-11d %-6d %-10.6f %-10.6f %-10.6f %-10.6f %-8.4f"
            % (N, r.size, xr, xq, mr, mq, c))

    say()
    say("W1  the unweighted parity ratio is locked too")
    say("    W1 %s" % ("hold" if ok1 else "REFUTED"))
    say("W2  the weighted ratio tracks the unweighted one")
    say("    W2 %s" % ("hold" if ok2 else "REFUTED"))
    say("W3  both maxima and both medians are in the table, reported")
    say("    not judged.")
    say()
    say("  If W1 and W2 hold, the primality of N - mk contributed the")
    say("  roughness of m and nothing further, and the lock is a")
    say("  statement about the parity of omega on rough integers below")
    say("  N/k. If either fails, it is not.")
    say()
    say("=" * 78)
    say("W1 %s  W2 %s" % tuple("hold" if v else "REFUTED"
                               for v in (ok1, ok2)))

    io.open(os.path.join(RES, "audit_signlock_why.txt"), "w",
            encoding="utf-8", newline="\n").write("\n".join(lines) + "\n")
    return 0 if (ok1 and ok2) else 1


if __name__ == "__main__":
    raise SystemExit(main())
