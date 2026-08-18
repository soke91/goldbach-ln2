# -*- coding: utf-8 -*-
r"""
What the split is worth, at each N, against what the signs are worth.

WHAT IS AT STAKE

Remark {#rem:residuesigned} closes with a comparison: "Remark
{#rem:splitbudget} measured the whole elementary/residue division at
about 0.06 in theta', and this one step at 0.21 to 0.29." A single
number on one side, a range on the other. Remark {#rem:signedgain}
has since shown the 0.21-0.29 is a declining series, so the shape of
the other side matters.

It is a series too. results/lab_split_budget.txt prints the two
exponents at five N -- log K*_H/log N of 0.6552 to 0.6716 and
log K*_R/log N of 0.7477 to 0.7382 -- and their difference, which is
what the split is worth, is not one number. **Comparing the smallest
value of one series with the whole range of another is not a
comparison.** Both sides have to be taken at the same N.

The published table was read before this script was written, so its
monotone appearance is not a prediction here; V3 states it as
something to check rather than to discover. What is not visible in
that table, and is what this measures, is whether either trend clears
its own noise and whether the RATIO of the two losses moves at all.

The implementation is the bitmask one of audit_level_slope_reach.py,
independent of lab_split_budget.py's, so V1 and V2 are cross-checks
of the published split and not reruns.

BACKS: Remark {#rem:splitvalue} in paper/theorem_A.md.

PRE-REGISTERED PREDICTIONS (written before this script was run)

  V1  The control: the three crossings reproduce the published
      K*_H, K*_P, K*_R at every N exactly.
  V2  The control: the three exponents reproduce the published values
      to within 0.0001.
  V3  The split's value in theta' -- the K*_R exponent minus the
      K*_H exponent -- falls with N, and the fall clears two standard
      errors. (The fall is visible in the published table; whether it
      clears its own noise is not.)
  V4  The ratio of the two losses does not move: the least-squares
      slope of (signed gain)/(split value) against log N stays below
      two standard errors, so "the signs are worth about three times
      the split" is a statement about the range and not a snapshot.

REFUTATION RULE (fixed before the run)

  V1  REFUTED on any mismatch -- not the same statistic, and nothing
      below may be compared with {#rem:splitbudget}.
  V2  REFUTED at 0.0001 at any N, likewise.
  V3  REFUTED if the slope fails to reach two standard errors. The
      split's value would then be a level, and since {#rem:signedgain}
      showed the signs' value declining, the ratio between them would
      have to be the moving quantity.
  V4  REFUTED if the ratio's slope reaches two standard errors. Then
      the two losses are being worked off at different rates and no
      single factor between them may be quoted -- the sentence in
      {#rem:residuesigned} would need an N attached.

  All four gate.

  NO NULL IS RUN and none applies. Crossings of deterministic sums
  against computed thresholds are located and differenced; there is no
  background to detect against. The null for the split itself is
  lab_split_budget.py's own permutation of each half's magnitudes
  across k, and for the signs lab_residue_signed.py's sixteen-draw
  randomisation; neither is repeated here.
"""

import io
import math
import os
import re
import sys

import numpy as np

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
RES = os.path.join(ROOT, "results")
OUT = os.path.join(RES, "audit_split_value.txt")

NS = [200_000, 400_000, 800_000, 1_600_000, 3_200_000]
KCAP = 100_000
QSIEVE = 30
CLIM = 4_000_000


def primes_upto(n):
    s = np.ones(n + 1, dtype=bool)
    s[:2] = False
    for p in range(2, int(n ** 0.5) + 1):
        if s[p]:
            s[p * p::p] = False
    return np.flatnonzero(s).astype(np.int64)


def lambda_and_mu(n):
    """von Mangoldt and Moebius, the cofactor kept in int32"""
    pr = primes_upto(n)
    lgp = np.log(pr.astype(np.float64))
    lam = np.zeros(n + 1, dtype=np.float64)
    lam[pr] = lgp
    for i, p in enumerate(pr):
        p = int(p)
        if p * p > n:
            break
        q = p * p
        while q <= n:
            lam[q] = lgp[i]
            if q > n // p:
                break
            q *= p
    del pr, lgp
    mu = np.ones(n + 1, dtype=np.int8)
    cof = np.arange(n + 1, dtype=np.int32)
    for p in primes_upto(int(math.isqrt(n))):
        p = int(p)
        mu[p::p] = -mu[p::p]
        if p * p <= n:
            mu[p * p::p * p] = 0
        cof[p::p] //= p
        pk = p * p
        while pk <= n:
            cof[pk::pk] //= p
            if pk > n // p:
                break
            pk *= p
    big = cof > 1
    del cof
    mu[big] = -mu[big]
    del big
    mu[0] = 0
    return lam, mu


def residue_mask(n, qs):
    """bit i of mask[v] is set exactly when qs[i] divides v"""
    m = np.zeros(n + 1, dtype=np.uint16)
    for i, q in enumerate(qs):
        m[0::q] |= np.uint16(1 << i)
    return m


def factor_set(n):
    v, out, d = n, set(), 2
    while d * d <= v:
        if v % d == 0:
            out.add(d)
            while v % d == 0:
                v //= d
        d += 1
    if v > 1:
        out.add(v)
    return out


def read_split():
    """the published crossings and exponents of the two halves"""
    src = io.open(os.path.join(RES, "lab_split_budget.txt"),
                  encoding="utf-8").read()
    i = src.index("N            K*_H      K*_P      K*_R")
    ks = {}
    for ln in src[i:].splitlines()[1:]:
        f = ln.split()
        if len(f) < 4 or not f[0].isdigit():
            break
        ks[int(f[0])] = (int(f[1]), int(f[2]), int(f[3]))
    j = src.index("log K*_H/log N   log K*_P/log N   log K*_R/log N")
    ex = {}
    for ln in src[j:].splitlines()[1:]:
        f = ln.split()
        if len(f) < 4 or not f[0].isdigit():
            break
        ex[int(f[0])] = (float(f[1]), float(f[2]), float(f[3]))
    return ks, ex


def read_gain():
    """the signed gain at each N, read from the results file"""
    src = io.open(os.path.join(RES, "audit_signed_gain.txt"),
                  encoding="utf-8").read()
    i = src.index("N            log K*/log N abs   signed    gain")
    g = {}
    for ln in src[i:].splitlines()[1:]:
        f = ln.split()
        if len(f) < 4 or not f[0].isdigit():
            break
        g[int(f[0])] = float(f[3])
    return g


def halves(N, lam, mu, sqf, vmask, qs, artin, twin):
    """the crossings of S(N)N by B_H, beta B_P and B_R"""
    PN = factor_set(N)
    S_ = twin
    for q in sorted(PN):
        if q > 2:
            S_ *= (1.0 + 1.0 / (q - 2.0))

    ks, Hs, Ps = [], [], []
    for k in range(2, KCAP):
        if not sqf[k]:
            continue
        if any(k % q == 0 for q in PN):
            continue
        M = (N - 1) // k
        if M < 2:
            continue
        ms = np.arange(1, M + 1, 2, dtype=np.int64)
        ms = ms[sqf[ms]]
        kb, ck = 0, 1.0
        for i, q in enumerate(qs):
            if k % q == 0:
                kb |= 1 << i
            else:
                ck *= q / (q - 1.0)
        for q in factor_set(k):
            if q > 2:
                ms = ms[ms % q != 0]
        if ms.size == 0:
            continue
        vals = N - ms * k
        g = mu[ms].astype(np.float64)
        keep = (vmask[vals] & np.uint16(~kb & 0xFFFF)) == 0
        ks.append(k)
        Hs.append(float((lam[vals] * g).sum()))
        Ps.append(ck * float(g[keep].sum()))
    ks = np.array(ks, dtype=np.int64)
    H = np.array(Hs)
    P = np.array(Ps)
    beta = float((H * P).sum() / (P * P).sum())
    R = H - beta * P
    w = np.log(ks.astype(np.float64))
    thr = S_ * N

    def cross(v):
        c = np.cumsum(w * np.abs(v))
        j = int(np.searchsorted(c, thr))
        return int(ks[j]) if j < ks.size else None

    return cross(H), cross(beta * P), cross(R), beta, S_


def fit(x, y):
    a, b = np.polyfit(x, y, 1)
    r = y - (a * x + b)
    n = x.size
    rms = float(np.sqrt((r ** 2).mean()))
    se = math.sqrt(float((r ** 2).sum() / (n - 2))
                   / float(((x - x.mean()) ** 2).sum()))
    return float(a), rms, se, abs(float(a)) / se


def main():
    lines = []

    def say(s=""):
        print(s)
        lines.append(s)

    pubk, pubex = read_split()
    gain = read_gain()
    say("read %d published crossings and %d exponent triples from "
        "results/lab_split_budget.txt," % (len(pubk), len(pubex)))
    say("  and %d signed gains from results/audit_signed_gain.txt"
        % len(gain))

    NMAX = max(NS)
    qs = [int(q) for q in primes_upto(QSIEVE) if q > 2]
    say("sieving to %d, sieve weight over the odd primes %s"
        % (NMAX, ", ".join(map(str, qs))))
    lam, mu = lambda_and_mu(NMAX)
    sqf = mu != 0
    vmask = residue_mask(NMAX, qs)

    artin, twin = 1.0, 2.0
    for p in primes_upto(CLIM):
        p = int(p)
        artin *= 1.0 - 1.0 / (p * (p - 1.0))
        if p > 2:
            twin *= 1.0 - 1.0 / (p - 1.0) ** 2

    got = []
    for N in NS:
        kh, kp, kr, beta, S_ = halves(
            N, lam, mu, sqf, vmask, qs, artin, twin)
        got.append((N, kh, kp, kr, beta))
        say("  N = %-10d beta %.4f  budget S(N) = %.4f  "
            "K*_H %-7s K*_P %-7s K*_R %s"
            % (N, beta, S_, kh, kp, kr))
        say("BUDGET kstar_splitvalue_SN_N%d %.6f" % (N, S_))
    say("RADICALS %d"
        % len(set(tuple(sorted(q for q in factor_set(g[0]) if q > 2))
                  for g in got)))

    # ------------------------------------------------------------- V1
    say()
    say("V1  the control: the three crossings")
    say("  N            K*_H here / published   K*_P            K*_R")
    v1 = True
    for N, kh, kp, kr, beta in got:
        p3 = pubk.get(N)
        if (kh, kp, kr) != p3:
            v1 = False
        say("  %-12d %-7s %-15s %-7s %-7s %-7s %s"
            % (N, kh, p3[0], kp, p3[1], kr, p3[2]))
    say("  V1 %s" % ("hold" if v1 else "REFUTED"))

    # ------------------------------------------------------------- V2
    say()
    say("V2  the control: the three exponents")
    say("  N            log K*_H/log N  log K*_P/log N  "
        "log K*_R/log N  worst diff")
    v2 = True
    exps = {}
    for N, kh, kp, kr, beta in got:
        e = tuple(math.log(k) / math.log(N) for k in (kh, kp, kr))
        exps[N] = e
        d = max(abs(a - b) for a, b in zip(e, pubex[N]))
        if not (d < 0.0001):
            v2 = False
        say("  %-12d %-15.4f %-15.4f %-15.4f %.6f"
            % (N, e[0], e[1], e[2], d))
    say("  V2 %s   (cap 0.0001)" % ("hold" if v2 else "REFUTED"))

    # ------------------------------------------------------------- V3
    say()
    say("V3  what the split is worth, at each N")
    say("  N            split value   signed gain   ratio")
    xs, sp, ra = [], [], []
    for N, kh, kp, kr, beta in got:
        e = exps[N]
        val = e[2] - e[0]
        xs.append(math.log(N))
        sp.append(val)
        ra.append(gain[N] / val)
        say("  %-12d %-13.4f %-13.4f %.4f"
            % (N, val, gain[N], gain[N] / val))
    x = np.array(xs)
    ysp = np.array(sp)
    a, rms, se, t = fit(x, ysp)
    say("  the split's least-squares slope against log N = %+.6f" % a)
    say("  r.m.s. residual %.4f, standard error %.6f, t = %.2f"
        % (rms, se, t))
    say("SCATTER slope_audit_split_value %.4f" % rms)
    say("TSTAT slope_audit_split_value %.2f" % t)
    say("SPREAD slope_audit_split_value %.4f" % float(x.max() - x.min()))
    if t < 2.0:
        say("UNRESOLVED SIGN slope_audit_split_value")
    v3 = (a < 0.0) and (t >= 2.0)
    say("  V3 the split's value falls and clears two s.e.   %s"
        % ("hold" if v3 else "REFUTED"))
    say("  the range it takes over this sweep, which is what may be")
    say("  quoted for it: %.4f to %.4f" % (min(sp), max(sp)))
    say("PERN split_value %d %.4f %.4f" % (len(sp), min(sp), max(sp)))
    say("PERN signed_gain %d %.4f %.4f"
        % (len(gain), min(gain.values()), max(gain.values())))

    # ------------------------------------------------------------- V4
    say()
    say("V4  and the ratio between the two losses")
    yra = np.array(ra)
    a2, rms2, se2, t2 = fit(x, yra)
    say("  ratio by N: " + ", ".join("%.4f" % v for v in ra))
    say("  slope %+.6f, r.m.s. residual %.4f, standard error %.6f, "
        "t = %.2f" % (a2, rms2, se2, t2))
    say("SCATTER slope_audit_split_value_ratio %.4f" % rms2)
    say("TSTAT slope_audit_split_value_ratio %.2f" % t2)
    say("SPREAD slope_audit_split_value_ratio %.4f"
        % float(x.max() - x.min()))
    if t2 < 2.0:
        say("UNRESOLVED SIGN slope_audit_split_value_ratio")
    v4 = t2 < 2.0
    say("RATIO signed_gain split_value %.4f %.4f" % (min(ra), max(ra)))
    say("  V4 the ratio does not move (t = %.2f)   %s"
        % (t2, "hold" if v4 else "REFUTED"))

    say()
    say("  DIAGNOSTIC (post hoc). What the comparison in")
    say("  {#rem:residuesigned} looks like when both sides are taken")
    say("  at the same N. Quoting the split at one end of its range")
    say("  and the signs across all of theirs gives")
    say("  %.4f / %.4f = %.2f; the per-N ratios above run %.2f to"
        % (max(gain.values()), min(sp), max(gain.values()) / min(sp),
           min(ra)))
    say("  %.2f. The comparison is real either way, and it is %.2f"
        % (max(ra), max(gain.values()) / min(sp) - max(ra)))
    say("  larger when the two sides are taken at different N.")

    say()
    say("=" * 70)
    ok = v1 and v2 and v3 and v4
    say("both losses decline together and their ratio is flat"
        if ok else "REFUTED")

    head = [
        "STATISTIC: the crossings of S(N)N by sum(log k)|H|, by",
        "           sum(log k)|beta P| and by sum(log k)|R|, their",
        "           exponents log K*/log N, the difference between the",
        "           residue's and H's -- which is what the",
        "           elementary/residue split is worth in theta' -- and",
        "           the ratio of that to the signed gain of Remark",
        "           {#rem:signedgain} at the same N, each with the",
        "           least-squares slope against log N, its r.m.s.",
        "           residual and its standard error.",
        "NULL: none is run and none applies. Crossings of deterministic",
        "      sums against computed thresholds are located and",
        "      differenced; there is no background to detect against.",
        "      The null for the split is lab_split_budget.py's",
        "      permutation of each half's magnitudes across k, and for",
        "      the signs lab_residue_signed.py's sixteen draws.",
        "FIELD: N = 2e5 through 3.2e6 by doubling; k squarefree and",
        "       coprime to N with 2 <= k < " + str(KCAP) + "; m odd,",
        "       squarefree and coprime to k, m <= (N-1)/k; the sieve",
        "       weight over the odd primes below " + str(QSIEVE) + ";",
        "       S(N) from an Euler product at the fixed bound "
        + str(CLIM) + ".",
        "       Every N is 2^a 5^b, one odd radical, as RADICALS says.",
        "       The published crossings are read from",
        "       results/lab_split_budget.txt and the signed gains from",
        "       results/audit_signed_gain.txt.",
        "",
    ]
    io.open(OUT, "w", encoding="utf-8", newline="\n").write(
        "\n".join(head + lines) + "\n")
    print("\nwrote %s" % OUT)
    if not ok:
        raise SystemExit(1)
    return 0


if __name__ == "__main__":
    sys.exit(main())
