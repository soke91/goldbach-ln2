# -*- coding: utf-8 -*-
r"""
Forcing the configuration the mechanism says is best

WHAT IS AT STAKE

rem:identityfloor measured the mechanism behind the multiplicative
ensemble's low end and found it dominant -- the ratio |T|/(S(N)N)
falls with m, the number of the ten smallest primes a draw sends to
-1, at t = -17.08 -- and then found something neither earlier remark
had suspected: **no draw of 512 reached m = 10.**  Ten independent
signs give that one time in 1024.  So the observed minimum 1.9149 is
the best of m <= 9, not the ensemble's infimum, and what the
mechanism's best configuration gives has never been seen.

The algebra says what to expect and it is not comfortable.  W(v)
vanishes as soon as two primes dividing v are negative, so forcing the
first J primes to -1 empties the support of W wherever v has two of
them; as J grows almost every v is killed and |T| falls.  **If it
falls past mu, then rem:identitydeep's "0 of 512" is a statement about
how rare such draws are, not about whether the ensemble contains
them** -- and the honest form of that remark becomes a quantile rather
than an exclusion.

This run forces the first J primes negative for J = 0, 5, 10, 15, 20,
25, 30, leaves the rest iid, and takes 64 draws at each.  The
prediction registered below is the one the algebra points at, not the
one that would be convenient.

BACKS: Remark {#rem:identityforced} in paper/wall_v3.md.

PRE-REGISTERED PREDICTIONS (written before this script was run)

  Q1  THE GATE.  mu's own ratios reproduce rem:identitynull's, 1.0039
      at N = 2e5 and 0.9865 at N = 4e5, to four decimals.
  Q2  The ratio falls with J: the median over the 64 draws is
      monotone decreasing across the seven values of J.
  Q3  **Forcing reaches mu.**  At some J <= 30, at least one draw's
      ratio is at or below mu's, at either N.  The multiplicative
      ensemble therefore does contain configurations reaching mu and
      rem:identitydeep's exclusion is about their rarity.
  Q4  And the rarity is quantifiable: for the smallest J that reaches
      mu, an unforced draw realises that configuration with
      probability 2^-J, so at least 2^J draws would be needed in
      expectation -- a number to set beside that remark's 512.

REFUTATION RULE (fixed before the run)

  Q1  REFUTED outside four decimals; nothing below is reported.
  Q2  REFUTED by any rise.  A rise would mean forcing more primes
      negative does not empty W monotonically, which the derivation
      does not allow, so it would point at this script rather than at
      the mathematics.
  Q3  **REFUTED if no J reaches mu.**  That is the stronger outcome
      for the paper and the one to state plainly if it comes: even
      the mechanism's most favourable configurations, forced rather
      than waited for, do not bring a multiplicative sign function to
      where mu is, and rem:identitydeep's claim is about the ensemble
      and not about depth.  **The unresolved case is named**: forcing
      only the first thirty primes is a choice, so a failure at
      J = 30 bounds the question at thirty and does not settle it,
      and the remark must say "not at J <= 30" rather than "not at
      all".
  Q4  Not independently refutable -- it is arithmetic on Q3's
      outcome, reported only if Q3 holds, and omitted if it does not.

  WHAT THIS CANNOT DO.  Two N and 64 draws per J.  Forcing changes
  the ensemble, so a ratio reached under forcing is not a ratio the
  unforced ensemble was ever seen to reach; the bridge between them
  is the probability in Q4 and nothing stronger.
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
OUT = os.path.join(ROOT, "results", "audit_identity_forced.txt")
SRC = os.path.join(ROOT, "results", "audit_directidentity_null.txt")

NS = [200_000, 400_000]
NMAX = max(NS)
CLIM = 4_000_000
JS = [0, 5, 10, 15, 20, 25, 30]
PERJ = 64
SEED = 20260910
DEC = 4


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
    return pr, lam, mu


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


def read_pub():
    src = io.open(SRC, encoding="utf-8").read()
    blk = src[src.index("V2  T / (S(N)N)"):]
    out = []
    for n in NS:
        m = re.search(r"^\s*%d\s+([\d.]+)\s" % n, blk, re.M)
        if not m:
            raise SystemExit("no V2 row for N = %d" % n)
        out.append(float(m.group(1)))
    return out


HEAD = [
    "STATISTIC: |T|/(S(N)N) for mu and for ensembles of %d random"
    % PERJ,
    "           multiplicative sign patterns in which the first J",
    "           primes are forced to -1 and the rest are iid, for",
    "           J = %s; the median and minimum at each J and" % (JS,),
    "           whether any draw reaches mu.",
    "FIELD: N = %s; k over every squarefree 2 <= k < N; S(N) from"
    % (NS,),
    "       the twin product at %d; field, k-range and S(N) taken"
    % CLIM,
    "       unchanged from code/audit_directidentity_null.py. mu's",
    "       published ratios are READ from",
    "       results/audit_directidentity_null.txt and re-measured here",
    "       as the gate.",
    "SEED: numpy default_rng at %d; without it the file does not"
    % SEED,
    "      reproduce its own control.",
    "",
]


def main():
    lines = []

    def say(t=""):
        print(t)
        sys.stdout.flush()
        lines.append(t)

    pub = read_pub()
    for n, v in zip(NS, pub):
        say("READ audit_directidentity_null.txt %d %.4f" % (n, v))
    say("  mu's T/(S(N)N) at the two N, read from that file")
    say("PRINTBOUND audit_identity_forced %d %.8f"
        % (DEC, 0.5 * 10.0 ** (-DEC)))
    say("  J over %s, %d draws at each" % (JS, PERJ))

    say("sieving to %d" % NMAX)
    pr, lam, mu = sieves(NMAX)
    sqf = mu != 0
    muf = mu.astype(np.float64)
    del mu
    twin = 2.0
    for p in primes_upto(CLIM):
        p = int(p)
        if p > 2:
            twin *= 1.0 - 1.0 / (p - 1.0) ** 2

    ks, lgs, idxs, Ss = {}, {}, {}, {}
    for N in NS:
        kk = np.flatnonzero(sqf[2:N]).astype(np.int64) + 2
        ks[N] = kk
        lgs[N] = np.log(kk.astype(np.float64))
        idxs[N] = np.arange(1, N, dtype=np.int64)
        S = twin
        for q in sorted(factor_set(N)):
            if q > 2:
                S *= (1.0 + 1.0 / (q - 2.0))
        Ss[N] = S

    def ratio_of(sg):
        out = []
        for N in NS:
            f0 = np.zeros(N, dtype=np.float64)
            f0[1:] = lam[1:N] * sg[N - idxs[N]]
            kk = ks[N]
            lg = lgs[N]
            s = 0.0
            for i, k in enumerate(kk):
                k = int(k)
                r = N % k
                a = f0[r::k].sum() if r else f0[k::k].sum()
                s += lg[i] * sg[k] * a
            del f0
            out.append(abs(s) / (Ss[N] * N))
        return out

    real = ratio_of(muf)

    # -------------------------------------------------------------- Q1
    say()
    say("Q1  does this construction reproduce mu's ratios?")
    ok = True
    for N, r, v in zip(NS, real, pub):
        g = abs(round(r, DEC) - round(v, DEC)) < 10.0 ** (-DEC) / 2
        ok &= g
        say("  N = %-8d here %.4f against its %.4f  %s"
            % (N, r, v, "ok" if g else "MISMATCH"))
    say("  Q1 %s   (cap: %d decimals at both)"
        % ("hold" if ok else "REFUTED", DEC))
    if not ok:
        io.open(OUT, "w", encoding="utf-8", newline="\n").write(
            "\n".join(HEAD + lines) + "\n")
        raise SystemExit(1)

    rng = np.random.default_rng(SEED)
    res = {}
    for J in JS:
        pool = np.zeros((PERJ, len(NS)))
        for d in range(PERJ):
            t = np.zeros(NMAX + 1, dtype=np.uint64)
            hi = rng.integers(0, 1 << 32, size=len(pr),
                              dtype=np.uint64)
            lo = rng.integers(0, 1 << 32, size=len(pr),
                              dtype=np.uint64)
            pm = (hi << np.uint64(32)) | lo
            bits = (pm & np.uint64(1)).astype(np.int64)
            bits[:J] = 1                     # forced negative
            del hi, lo, pm
            sg = np.zeros(NMAX + 1, dtype=np.float64)
            par = np.zeros(NMAX + 1, dtype=np.int8)
            for i, p in enumerate(pr):
                if bits[i]:
                    par[int(p)::int(p)] ^= 1
            sg[sqf] = 1.0 - 2.0 * par[sqf].astype(np.float64)
            pool[d] = ratio_of(sg)
            del sg, par, t, bits
        res[J] = pool
        say("  J = %-3d done, %d draws" % (J, PERJ))

    say()
    say("        J   forced primes to        N=%d min   median"
        "      N=%d min   median" % (NS[0], NS[1]))
    meds, mins = [], []
    for J in JS:
        pool = res[J]
        top = int(pr[J - 1]) if J else 0
        meds.append(float(np.median(pool[:, 0])))
        mins.append([float(pool[:, i].min()) for i in range(len(NS))])
        say("      %3d   %16d   %9.4f  %8.4f   %9.4f  %8.4f"
            % (J, top, mins[-1][0], float(np.median(pool[:, 0])),
               mins[-1][1], float(np.median(pool[:, 1]))))
        say("POINT forcedmin_%d %.5f" % (J, mins[-1][0]))
    say("SCALES %d" % len(JS))
    say("  mu is at %.4f and %.4f" % (real[0], real[1]))

    # -------------------------------------------------------------- Q2
    say()
    say("Q2  does the ratio fall with J?")
    rises = [(JS[i], JS[i + 1]) for i in range(len(JS) - 1)
             if meds[i + 1] > meds[i]]
    q2 = not rises
    say("  medians at N=%d: %s"
        % (NS[0], ", ".join("%.4f" % m for m in meds)))
    if rises:
        say("  rises at %s" % ", ".join("J %d to %d" % r
                                        for r in rises))
    say("  Q2 %s   (cap: no rise)" % ("hold" if q2 else "REFUTED"))

    # -------------------------------------------------------------- Q3
    say()
    say("Q3  does forcing reach mu?")
    hit = None
    for J in JS:
        pool = res[J]
        if any(float(pool[:, i].min()) <= real[i]
               for i in range(len(NS))):
            hit = J
            break
    q3 = hit is not None
    if q3:
        say("  the smallest J that reaches mu is %d" % hit)
    else:
        say("  no draw at any J reaches mu; the closest is %.4f "
            "against mu's %.4f"
            % (min(m[0] for m in mins), real[0]))
    say("COUNT forced_hitJ %d" % (hit if q3 else -1))
    say("  Q3 %s   (cap: some J <= %d)"
        % ("hold" if q3 else "REFUTED", max(JS)))

    # -------------------------------------------------------------- Q4
    say()
    say("Q4  how rare is that configuration?")
    if q3:
        say("  an unforced draw realises J = %d negatives among the "
            "first %d primes" % (hit, hit))
        say("  with probability 2^-%d, so about %.3e draws would be "
            "needed" % (hit, 2.0 ** hit))
        say("POINT forced_rarity %.5e" % (2.0 ** hit))
        say("  against rem:identitydeep's 512")
    else:
        say("  Q3 did not hold, so the smallest J that reaches mu -- "
            "the object")
        say("  Q4 prices -- does not exist. Q4 asserted it did, so Q4 "
            "fails with it.")
        say("  Q4 REFUTED   (its subject is empty)")
        say("  NOTE: the rule said Q4 would be 'omitted' if Q3 fell. "
            "A prediction")
        say("  that can be omitted is not one. The rule is not "
            "changed here; it is")
        say("  applied -- an existence claim whose subject is empty "
            "is refuted --")
        say("  and its misspecification is recorded as the finding it "
            "is.")

    say()
    say("=" * 70)
    say("Q1 %s  Q2 %s  Q3 %s"
        % tuple("hold" if v else "REFUTED" for v in (ok, q2, q3)))
    say()
    if q3:
        say("the multiplicative ensemble does contain configurations "
            "that reach mu.")
        say("rem:identitydeep's exclusion is therefore about how rare "
            "they are and")
        say("not about whether they exist, and the honest form of "
            "that remark is a")
        say("quantile: mu sits where a draw with the first %d primes "
            "negative sits," % hit)
        say("which 512 draws had no business finding.")
    else:
        say("even forced, the mechanism's best configurations do not "
            "bring a")
        say("multiplicative sign function to mu. rem:identitydeep's "
            "claim is about")
        say("the ensemble and not about depth -- bounded at J <= %d, "
            "which is a" % max(JS))
        say("choice this run made and not a limit it established.")

    io.open(OUT, "w", encoding="utf-8", newline="\n").write(
        "\n".join(HEAD + lines) + "\n")
    print("\nwrote %s" % os.path.normpath(OUT))
    return 0


if __name__ == "__main__":
    sys.exit(main())
