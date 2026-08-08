# -*- coding: utf-8 -*-
r"""
What the signs across k are worth once the elementary half is gone.

WHAT IS AT STAKE

[eq:direct] contains the signed sum over k; [eq:directcond] throws the
signs away because a bound on |H(N;k)| is what an estimate supplies.
Remark {#rem:signedlevel} priced that discard for H: keeping the signs
moves K* by a factor 1.793 to 2.223, about 0.053 in theta', and the
control says the gain is not mu's doing -- holding every |H(N;k)|
fixed and redrawing the signs, all 16 draws fail to cross anywhere in
the walk. **mu is worse than random signs there.**

But Remark {#rem:residue} located that lean in the elementary half.
P's mass-weighted f+ is 0.0796 to 0.1836, an order below the sign-draw
band; R's is 0.5516 to 0.4832, nearly centred. So the object whose
signs across k are pathological is P, and the object the conditional
reduction of {#rem:provablehalf} leaves behind is R.

Nobody has priced the signs for R. If R's signs across k are
random-like, the signed residue sum is a walk rather than a drift, and
the level it permits is larger than the 0.5654 to 0.5799 that
{#rem:residuelevel} measured from absolute values -- which is the
knife-edge the programme currently sits on.

BACKS: Remark {#rem:residuesigned} in paper/theorem_A.md.

PRE-REGISTERED PREDICTIONS (written before this script was run)

  S1  The control: the absolute crossing reproduces the operative
      K*_R of results/audit_residue_level.txt -- 993, 1447, 2019,
      3319, 5923 -- to within 1 per cent at every N.
  S2  Keeping the signs helps: the signed crossing is later than the
      absolute one at every N.
  S3  And unlike H, R's signs are not worse than random: mu's signed
      crossing is not earlier than the median of 16 draws that hold
      every |R(N;k)| fixed and redraw the sign of each term.
  S4  The gain clears the barrier the absolute sum missed: the signed
      exponent exceeds 0.56 at every N, including N = 8e5 where the
      absolute one reads 0.5599.

REFUTATION RULE (fixed before the run)

  S1  REFUTED at 1 per cent at any N.
  S2  REFUTED if the signed crossing is at or before the absolute one
      at any N.
  S3  REFUTED if mu crosses earlier than the draws' median at any N.
      That is the outcome worth having either way: a refutation would
      say the lean {#rem:signedlevel} found survives the split, and
      the elementary half is not where it lives.
  S4  REFUTED if the signed exponent fails to exceed 0.56 at any N.

  All four gate.

  NULL: the sign randomisation of S3 -- every |R(N;k)| held fixed and
  the sign of each term redrawn, 16 draws. A coin on supp(mu^2) is NOT
  used: it changes the magnitudes as well, and by {#rem:whycoinwins}
  it beats mu for reasons unrelated to the signs across k. Permuting
  signs alone isolates the one thing at issue, and is the same control
  lab_signed_level.py ran for H so that the two answers are
  comparable.
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
OUT = os.path.join(ROOT, "results", "lab_residue_signed.txt")

NS = [200_000, 400_000, 800_000, 1_600_000, 3_200_000]
KCAP = 100_000
QSIEVE = 30
CLIM = 4_000_000
THETA = 0.56
DRAWS = 16
SEED = 20260808


def primes_upto(n):
    s = np.ones(n + 1, dtype=bool)
    s[:2] = False
    for p in range(2, int(n ** 0.5) + 1):
        if s[p]:
            s[p * p::p] = False
    return np.flatnonzero(s).astype(np.int64)


def sieves(n):
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
    mu = np.ones(n + 1, dtype=np.int8)
    rem = np.arange(n + 1, dtype=np.int32)
    for p in primes_upto(int(math.isqrt(n))):
        p = int(p)
        mu[p::p] = -mu[p::p]
        if p * p <= n:
            mu[p * p::p * p] = 0
        q = p
        while q <= n:
            rem[q::q] //= p
            if q > n // p:
                break
            q *= p
    big = rem > 1
    del rem
    mu[big] = -mu[big]
    del big
    mu[0] = 0
    return lam, mu


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


def read_H_ratios():
    """what the signs bought for H -- read, not copied"""
    p = os.path.join(ROOT, "results", "lab_signed_level.txt")
    src = io.open(p, encoding="utf-8").read()
    m = re.search(r"Y4  K\*_signed/K\*_H: ([\d., ]+?)\s+\(", src)
    return [float(v) for v in m.group(1).split(",")]


def read_published():
    """the operative absolute crossings -- read, not copied"""
    p = os.path.join(ROOT, "results", "audit_residue_level.txt")
    src = io.open(p, encoding="utf-8").read()
    i = src.index("log K*_R/log N")
    out = {}
    for ln in src[i:].splitlines()[1:]:
        f = ln.split()
        if len(f) < 4 or not f[0].isdigit():
            break
        out[int(f[0])] = int(f[2])
    return out


def main():
    lines = []

    def say(s=""):
        print(s)
        lines.append(s)

    pub = read_published()
    say("read %d operative K*_R from results/audit_residue_level.txt"
        % len(pub))

    NMAX = max(NS)
    say("sieving to %d ..." % NMAX)
    lam, mu = sieves(NMAX)
    sqf = mu != 0
    QS = [int(q) for q in primes_upto(QSIEVE) if q > 2]
    say("  sieve weight over the odd primes %s"
        % ", ".join(map(str, QS)))

    artin, twin = 1.0, 2.0
    for p in primes_upto(CLIM):
        p = int(p)
        artin *= 1.0 - 1.0 / (p * (p - 1.0))
        if p > 2:
            twin *= 1.0 - 1.0 / (p - 1.0) ** 2

    rng = np.random.default_rng(SEED)
    res = []
    for N in NS:
        PN = factor_set(N)
        A_, S_ = artin, twin
        for q in sorted(PN):
            A_ /= (1.0 - 1.0 / (q * (q - 1.0)))
            if q > 2:
                S_ *= (1.0 + 1.0 / (q - 2.0))

        ks, Hs, Ps = [], [], []
        for k in range(2, KCAP):
            if not sqf[k] or any(k % q == 0 for q in PN):
                continue
            M = (N - 1) // k
            if M < 2:
                continue
            ms = np.arange(1, M + 1, 2, dtype=np.int64)
            ms = ms[sqf[ms]]
            for q in factor_set(k):
                if q > 2:
                    ms = ms[ms % q != 0]
            if ms.size == 0:
                continue
            vals = N - ms * k
            g = mu[ms].astype(np.float64)
            w = np.ones(ms.size, dtype=np.float64)
            for q in QS:
                if k % q == 0:
                    continue
                w *= np.where(vals % q == 0, 0.0, q / (q - 1.0))
            ks.append(k)
            Hs.append(float((lam[vals] * g).sum()))
            Ps.append(float((g * w).sum()))
        ks = np.array(ks, dtype=np.int64)
        H = np.array(Hs)
        P = np.array(Ps)
        lw = np.log(ks.astype(float))
        beta = float((H * P).sum() / (P * P).sum())
        R = H - beta * P
        thr = S_ * (1.0 - A_) * N
        res.append((N, ks, lw, R, thr, beta, S_, A_))
        say("  N = %-10d  #k = %-7d beta = %.6f" % (N, ks.size, beta))

    def cross_abs(ks, lw, R, thr):
        cum = np.cumsum(lw * np.abs(R))
        j = int(np.searchsorted(cum, thr))
        return int(ks[j]) if j < ks.size else None

    def cross_signed(ks, lw, R, thr):
        """first k at which the signed walk leaves [-thr, thr]"""
        cum = np.cumsum(lw * R)
        bad = np.flatnonzero(np.abs(cum) > thr)
        return int(ks[bad[0]]) if bad.size else None

    def expo(N, k):
        return math.log(k) / math.log(N) if k else float("nan")

    # ------------------------------------------------------------- S1
    say()
    say("S1  the control: the absolute crossing")
    say("  N            K*_R (here)   published    ratio")
    s1 = True
    for N, ks, lw, R, thr, beta, S_, A_ in res:
        k = cross_abs(ks, lw, R, thr)
        r = k / pub[N] if k else float("nan")
        if not (abs(r - 1.0) < 0.01):
            s1 = False
        say("  %-12d %-13s %-12d %.4f" % (N, str(k), pub[N], r))
    say("  S1 %s" % ("hold" if s1 else "REFUTED"))

    # ------------------------------------------------------------- S2
    say()
    say("S2/S4  keeping the signs across k")
    say("  N            K*_R abs   K*_R signed  factor   exponent  "
        "clears")
    s2 = s4 = True
    sg = []
    for N, ks, lw, R, thr, beta, S_, A_ in res:
        ka = cross_abs(ks, lw, R, thr)
        kv = cross_signed(ks, lw, R, thr)
        sg.append((ka, kv))
        if kv is None:
            say("  %-12d %-10d %-12s -        -         yes (no "
                "crossing below k = %d)" % (N, ka, "none", KCAP))
            continue
        if kv <= ka:
            s2 = False
        e = expo(N, kv)
        if e <= THETA:
            s4 = False
        say("  %-12d %-10d %-12d %-8.3f %-9.4f %s"
            % (N, ka, kv, kv / ka, e, "yes" if e > THETA else "NO"))
    say("  S2 signed later than absolute at every N   %s"
        % ("hold" if s2 else "REFUTED"))
    say("  S4 signed exponent above %.2f at every N   %s"
        % (THETA, "hold" if s4 else "REFUTED"))

    # ------------------------------------------------------------- S3
    say()
    say("S3  the control: same |R(N;k)|, signs redrawn, %d draws"
        % DRAWS)
    say("  N            mu         draws min  median     max        "
        "verdict")
    s3 = True
    for N, ks, lw, R, thr, beta, S_, A_ in res:
        aR = np.abs(R)
        kv = cross_signed(ks, lw, R, thr)
        got = []
        for _ in range(DRAWS):
            s = rng.integers(0, 2, size=aR.size) * 2.0 - 1.0
            k2 = cross_signed(ks, lw, s * aR, thr)
            got.append(k2 if k2 is not None else KCAP)
        got.sort()
        med = float(np.median(got))
        mine = kv if kv is not None else KCAP
        if mine < med:
            s3 = False
        say("  %-12d %-10s %-10d %-10.0f %-10d %s"
            % (N, str(kv) if kv else ">%d" % KCAP, got[0], med, got[-1],
               "not earlier" if mine >= med else "EARLIER"))
    say("  S3 %s" % ("hold" if s3 else "REFUTED"))

    say()
    say("  the budget constant crossed throughout, declared:")
    for N, ks, lw, R, thr, beta, S_, A_ in res:
        say("BUDGET kstar_R_signed_S1AN_N%d %.6f" % (N, S_ * (1.0 - A_)))

    say()
    say("  DIAGNOSTIC (post hoc). The same comparison for H, which")
    say("  lab_signed_level.py ran against the OTHER budget, is not")
    say("  directly comparable; what is comparable is the ratio the")
    hr = read_H_ratios()
    say("  signs buy. For H it was %.3f to %.3f in K*, read from"
        % (min(hr), max(hr)))
    say("  results/lab_signed_level.txt. Here:")
    say("  N            signed/absolute in K*   in the exponent")
    for i, (N, ks, lw, R, thr, beta, S_, A_) in enumerate(res):
        ka, kv = sg[i]
        if kv is None:
            say("  %-12d %-23s %s" % (N, ">%.1f" % (KCAP / ka), "-"))
            continue
        say("  %-12d %-23.3f %+.4f"
            % (N, kv / ka, expo(N, kv) - expo(N, ka)))

    say()
    say("=" * 70)
    ok = s1 and s2 and s3 and s4
    say("the residue's signs across k are random-like and buy the "
        "level the absolute sum could not" if ok else "REFUTED")

    head = [
        "STATISTIC: the first crossing of sum_{k<K}(log k)|R(N;k)| above",
        "           S(N)(1-A(N))N, the first k at which the SIGNED walk",
        "           sum_{k<K}(log k)R(N;k) leaves [-S(N)(1-A(N))N,",
        "           +S(N)(1-A(N))N], their ratio and exponents, and the",
        "           same signed crossing with the sign of each term",
        "           redrawn at random.",
        "NULL: the sign randomisation of S3 -- every |R(N;k)| held fixed",
        "      and the sign of each term redrawn, 16 draws. A coin on",
        "      supp(mu^2) is not used: it changes the magnitudes as well",
        "      and beats mu for reasons unrelated to the signs across k.",
        "      This is the same control lab_signed_level.py ran for H,",
        "      so that the two answers are comparable.",
        "FIELD: N = 2e5 to 3.2e6 by doubling; k squarefree and coprime",
        "       to N with 2 <= k < 100000; m odd, squarefree, coprime to",
        "       the odd part of k, m < N/k; the sieve weight uses the",
        "       odd primes up to 30; beta refitted as sum(H P)/sum(P^2)",
        "       on the same k-range; S(N) and A(N) from Euler products",
        "       at the fixed bound 4000000; the published absolute",
        "       crossings are read from results/audit_residue_level.txt;",
        "       numpy default_rng seed 20260808.",
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
