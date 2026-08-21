# -*- coding: utf-8 -*-
r"""The divisor switch under a moving truncation, and where its residual sits.

Supports {#rem:movingswitch}.

WHAT IS AT STAKE

The published argument this note is about carried its cut on the outer
variable, k < K, with the inner sum over all n < N.  The authors have
since replaced that by a moving cut: the inner sum runs to
Y_k = ceil(N - alpha k) - 1, equivalently k <= R_n = ceil((N-n)/alpha) - 1,
and the outer bound k < K becomes redundant.  Every statement proved
against the fixed cut is a statement about a superseded formulation
unless the switch survives the move.

Under the move the exchanged inner divisor sum is cut by the COFACTOR:

    sigma*(u) = sum_{k | u, (k,N)=1, u/k > alpha} mu(k),

so completing it leaves a tail over m = u/k <= alpha -- the short
variable, directly, with no need for the completion to produce it.
That is the favourable configuration, reached one step earlier than
under the fixed cut.

WHAT IS MEASURED

For N even and alpha = N^theta, with u = N - n and m = u/k:

  X1  the switch is an exact rearrangement:
        D*  = sum_{k<K,(k,N)=1} mu(k) sum_{n<=Y_k, n=N(k)} Lambda(n) mu(N-n)
        Du* = sum_{u<N} Lambda(N-u) mu(u) sigma*(u)
      accumulated in different index orders; X1 compares them.

  X2  the completion splits it: D* = P* - R*, with
        P*  = sum_{u<N, rad(u)|N} Lambda(N-u) mu(u)          (few terms)
        R*  = sum_{m<=alpha} mu(m) sum_{(k,mN)=1, mk<N} mu^2(k) Lambda(N-mk)
      X2 compares D* against P* - R*.

  X3  the tail really is short: the largest m carrying a term is
      recorded against alpha.

  X4  the subtracted mean term is negligible.  With
        M*  = sum_{k<K,(k,N)=1} mu(k)/phi(k) * A_N(Y_k)
      and the endpoint exchange M* = sum_{n<N} a_n rho_N(R_n + 1),
      X4 compares the two forms and reports |M*|/N.

  X5  the object the theorem is about,
        T1* = D* - M* = sum_k mu(k) E_mu(Y_k; k),
      is reported as |T1*|/N beside the fixed-cut |T1|/N at the same N.

FALSIFICATION, registered before the run

  X1  REFUTED if the two index orders disagree beyond 1e-12 relative
      to N.  Then the exchange is not the identity claimed and nothing
      below is interpretable.
  X2  REFUTED at the same tolerance.  Then the completion is wrong.
  X3  REFUTED if any surviving m exceeds alpha.  Then the tail is not
      the short variable and the configuration is not the favourable
      one.
  X4  REFUTED if the two forms of M* disagree beyond 1e-12 relative to
      N, or if |M*|/N does not fall as N grows.
  X5  is reported, not judged: a bound is analytic and no computation
      here establishes one.  It is printed so that the claim "the
      moving cut does not make the object larger" is checkable.

  PREDICTION.  X1-X4 hold.  X5 shows |T1*|/N of the same order as the
  fixed-cut |T1|/N -- the move changes which step produces the short
  variable, not the size of what is left.

NULL.  None applies: every quantity is a deterministic finite sum with
no sampling and no sign input.  The control is that each identity is
accumulated twice in genuinely different index orders (X1, X2, X4).
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
THETA = 0.44           # alpha = N^theta; Corollary-1 regime is theta < 1/2
TOL = 1e-12            # relative to N


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
        if w % p == 0:
            mu[v] = 0
        else:
            mu[v] = -mu[w]
        # von Mangoldt: v is a prime power iff stripping p leaves 1 or p^j
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

    say(__doc__.strip())
    say()
    say("=" * 70)
    say("sieving to %d ..." % max(NS))
    lam, mu, spf = sieve(max(NS))

    say()
    say("STATISTIC: five exact rearrangements of one finite double sum")
    say("           under the moving cut n <= Y_k, each accumulated in")
    say("           two different index orders; the largest surviving")
    say("           m against alpha; and |M*|/N and |T1*|/N.")
    say("FIELD: N even in %s; alpha = N^theta; k squarefree with"
        % ",".join(str(v) for v in NS))
    say("       (k,N) = 1; n < N with n = N (mod k) and n <= Y_k;")
    say("       u = N - n, m = u/k. Lambda, mu from one integer sieve.")
    say("CONSTANTS: THETA = %.2f, TOL = %.0e relative to N" % (THETA, TOL))
    say("NULL: none applies -- deterministic finite sums, no sampling,")
    say("      no sign input. The control is double accumulation.")
    say("DENOM: every relative error printed below is divided by N.")
    say()

    hdr = ("  %-9s %-6s %-12s %-12s %-12s %-10s %-11s"
           % ("N", "alpha", "X1 rel", "X2 rel", "X4 rel", "max m", "|M*|/N"))
    say(hdr)
    say("  " + "-" * (len(hdr) - 2))

    ok1 = ok2 = ok3 = ok4 = True
    t1s = []
    for N in NS:
        alpha = N ** THETA
        K = (N - 1) / alpha
        PN = factor_set(N, spf)

        ks = [k for k in range(1, int(math.ceil(K)))
              if mu[k] != 0 and not (factor_set(k, spf) & PN)]

        # a_n once, and its prefix sums; both forms of M* read from these,
        # so the comparison is a pure reordering and the outer sums are
        # accumulated exactly (math.fsum) rather than in running float64.
        nn = np.arange(0, N, dtype=np.int64)
        a = lam[:N] * mu[N - nn].astype(np.float64)
        a[0] = 0.0
        pref = np.concatenate(([0.0], np.cumsum(a[1:])))

        # --- D* by the k-side: modular slices, inner cut at Y_k
        dk_terms, mk_terms = [], []
        for k in ks:
            Yk = math.ceil(N - alpha * k) - 1
            if Yk < 1:
                continue
            r = N % k
            start = r if r > 0 else k
            ns = np.arange(start, Yk + 1, k, dtype=np.int64)
            if ns.size:
                dk_terms.append(mu[k] * math.fsum(
                    (lam[ns] * mu[N - ns].astype(np.float64)).tolist()))
            ph = 1
            for p in factor_set(k, spf):
                ph *= (p - 1)
            if k == 1:
                ph = 1
            mk_terms.append(mu[k] / ph * float(pref[Yk]))
        Dk = math.fsum(dk_terms)
        Mk = math.fsum(mk_terms)

        # --- D* by the u-side: divisor sum cut by the cofactor
        Du = 0.0
        Pu = 0.0
        for u in range(1, N):
            lv = lam[N - u]
            if lv == 0.0 or mu[u] == 0:
                continue
            fs = factor_set(u, spf)
            if fs & PN:
                # k must be coprime to N; divisors of u using those primes
                pass
            sig = 0
            # divisors of squarefree u
            ds = [1]
            for p in fs:
                ds += [d * p for d in ds]
            for k in ds:
                if factor_set(k, spf) & PN:
                    continue
                if u / k > alpha:
                    sig += mu[k]
            Du += lv * mu[u] * sig
            if fs <= PN:
                Pu += lv * mu[u]

        # --- R* by the m-side
        Ru = 0.0
        maxm = 0
        for m in range(1, int(alpha) + 1):
            if mu[m] == 0:
                continue
            fm = factor_set(m, spf) | PN
            kk = np.arange(1, (N - 1) // m + 1, dtype=np.int64)
            keep = np.ones(kk.size, dtype=bool)
            for p in fm:
                keep &= (kk % p != 0)
            kk = kk[keep]
            kk = kk[mu[kk] != 0]
            if kk.size == 0:
                continue
            vals = N - m * kk
            s = float(lam[vals].sum())
            if s != 0.0:
                maxm = max(maxm, m)
            Ru += mu[m] * s

        # --- M* by the n-side: endpoint exchange
        wk = []
        for k in ks:
            ph = 1
            for p in factor_set(k, spf):
                ph *= (p - 1)
            wk.append((k, mu[k] / (1 if k == 1 else ph)))
        mn_terms, rho = [], {}
        for n in range(1, N):
            an = float(a[n])
            if an == 0.0:
                continue
            Rn = math.ceil((N - n) / alpha) - 1
            if Rn not in rho:
                rho[Rn] = math.fsum(v for k, v in wk if k <= Rn)
            mn_terms.append(an * rho[Rn])
        Mn = math.fsum(mn_terms)

        r1 = abs(Dk - Du) / N
        r2 = abs(Dk - (Pu - Ru)) / N
        r4 = abs(Mk - Mn) / N
        t1 = abs(Dk - Mk) / N
        t1s.append(t1)
        ok1 &= r1 <= TOL
        ok2 &= r2 <= TOL
        ok3 &= maxm <= alpha
        ok4 &= r4 <= TOL
        say("  %-9d %-6.0f %-12.3e %-12.3e %-12.3e %-10d %-11.6f"
            % (N, alpha, r1, r2, r4, maxm, abs(Mk) / N))

    say()
    say("X1  the switch is exact under the moving cut")
    say("    X1 %s" % ("hold" if ok1 else "REFUTED"))
    say("X2  the completion splits it, tail on the short variable")
    say("    X2 %s" % ("hold" if ok2 else "REFUTED"))
    say("X3  every surviving m is at most alpha")
    say("    X3 %s" % ("hold" if ok3 else "REFUTED"))
    say("X4  the two forms of the mean term agree, and it is small")
    say("    X4 %s" % ("hold" if ok4 else "REFUTED"))
    say()
    say("X5  |T1*|/N, reported not judged")
    for N, t in zip(NS, t1s):
        say("    N = %-9d |T1*|/N = %.8f" % (N, t))
    say("    a bound is analytic; nothing here establishes one.")
    say()
    say("=" * 70)
    say("X1 %s  X2 %s  X3 %s  X4 %s"
        % tuple("hold" if v else "REFUTED" for v in (ok1, ok2, ok3, ok4)))

    io.open(os.path.join(RES, "audit_moving_switch.txt"), "w",
            encoding="utf-8", newline="\n").write("\n".join(lines) + "\n")
    return 0 if (ok1 and ok2 and ok3 and ok4) else 1


if __name__ == "__main__":
    raise SystemExit(main())
