# -*- coding: utf-8 -*-
r"""
The floor under the multiplicative ensemble, and where it comes from

WHAT IS AT STAKE

Two runs have now seen the same thing from different directions.
rem:identityseeds found the ensemble's minimum moving by at most
0.1702 across three seeds while its median moved by 1.1706;
rem:identitybig found the minimum identical to four decimals over the
first 64, 128 and 256 draws, a drift of 0.0000.  **The low end of this
ensemble behaves like a floor rather than a tail**, and both remarks
offered the same reading and both called it an observation: draws in
which many f(p) = -1 make the divisor sums vanish wholesale.

The reading is not a guess.  For multiplicative s, writing
log k = sum over p | k of log p and collecting,

    W(v) = sum_{k | v} s(k) log k
         = sum_{p | v} log(p) s(p) prod_{q | v, q != p} (1 + s(q)),

and T = sum_v Lambda(N-v) s(v) W(v).  **If two or more primes dividing
v carry s(q) = -1 then every term of that outer sum has a zero factor
and W(v) = 0.**  So W is supported on the v with at most one negative
prime, and the more primes a draw sends to -1 the sparser that support
becomes.  That is a derivation, not a hypothesis, and it makes the
floor predictable rather than merely visible.

What is not derived is whether the effect is the one actually driving
the measured minimum, or whether something else dominates and this is
a spectator.  That is what is measured here.

Every draw is scored by m, the number of the ten smallest primes it
sends to -1, and the ratio |T|/(S(N)N) is put against it.  The seed,
sieve top, pass structure and draw order are those of
audit_directidentity_deep.py, so draw d here is draw d there and the
minima must reproduce.

BACKS: Remark {#rem:identityfloor} in paper/wall_v3.md.

PRE-REGISTERED PREDICTIONS (written before this script was run)

  P1  THE GATE.  The minima reproduce rem:identitydeep's, 1.9149 at
      N = 2e5 and 2.0030 at N = 4e5, to four decimals.  Same seed,
      same sieve top, same pass structure: if the draws are not the
      same draws, nothing below is about that run's floor.
  P2  **The ratio falls with m**: regressing log of the ratio on m
      over the 512 draws gives a negative slope resolved at |t| > 2.
  P3  And the floor is where m is largest: the draw achieving the
      minimum sits in the top decile of m.
  P4  **It is a floor and not a tail**: among the draws in the top
      decile of m the ratio's standard deviation is at most half what
      it is among the bottom decile.

REFUTATION RULE (fixed before the run)

  P1  REFUTED outside four decimals; nothing below is reported.
  P2  REFUTED if the slope is positive or unresolved.  **Unresolved
      is the case to watch** -- m over ten primes takes eleven values
      and 512 draws concentrate near the middle, so a real dependence
      could still fail to resolve, and that would mean this design
      cannot see it rather than that it is absent.
  P3  REFUTED if the minimising draw is outside the top decile.  Then
      the derived mechanism is a spectator: it is still true that
      many negative primes empty W, but something else is choosing
      the minimum, and what the minimising draw's m actually is would
      be the finding.
  P4  REFUTED above half.  Then the low end is a tail after all and
      the two earlier remarks' reading was wrong, which would matter
      more than P2 or P3 -- the stability those runs measured would
      need another explanation.

  WHAT THIS CANNOT DO.  Ten primes is a choice, fixed above; a
  different count could sort the draws differently and this run does
  not vary it.  Nothing here measures the support of W directly, only
  the ratio it produces.

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
OUT = os.path.join(ROOT, "results", "audit_identity_floor.txt")
SRC = os.path.join(ROOT, "results", "audit_directidentity_deep.txt")

NS = [200_000, 400_000]
NMAX = 800_000            # the sieve top of the run being reproduced
CLIM = 4_000_000
PASSES = 8
PERPASS = 64
DRAWS = PASSES * PERPASS
SEED = 20260906           # the seed of the run being reproduced
NSMALL = 10
DEC = 4
HALF = 0.5


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


def ols(x, y):
    A = np.column_stack([np.ones(len(y)), x])
    c, *_ = np.linalg.lstsq(A, y, rcond=None)
    r = y - A.dot(c)
    s2 = float((r ** 2).sum()) / (len(y) - 2)
    cov = s2 * np.linalg.inv(A.T.dot(A))
    return c, np.sqrt(np.diag(cov))


def read_mins():
    src = io.open(SRC, encoding="utf-8").read()
    out = []
    for n in NS:
        m = re.search(r"^\s*multiplicative\s+%d\s+[\d.]+\s+"
                      r"([\d.]+)\s" % n, src, re.M)
        if not m:
            raise SystemExit("no multiplicative row for N = %d" % n)
        out.append(float(m.group(1)))
    return out


HEAD = [
    "STATISTIC: for each of %d random multiplicative draws, the" % DRAWS,
    "           ratio |T|/(S(N)N) and m, the number of the %d"
    % NSMALL,
    "           smallest primes the draw sends to -1; the slope of",
    "           log(ratio) on m, where the minimising draw sits in m,",
    "           and the ratio's spread in the top and bottom deciles",
    "           of m.",
    "FIELD: N = %s; k over every squarefree 2 <= k < N; S(N) from"
    % (NS,),
    "       the twin product at %d. Seed, sieve top (%d), pass"
    % (CLIM, NMAX),
    "       structure and draw order are those of",
    "       code/audit_directidentity_deep.py, so draw d here is draw",
    "       d there; that run's minima are READ from its result file",
    "       and re-measured here as the gate.",
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

    pub = read_mins()
    for n, v in zip(NS, pub):
        say("READ audit_directidentity_deep.txt multiplicative %.4f"
            % v)
    say("  the minima this run has to reproduce, at N = %s" % (NS,))
    say("PRINTBOUND audit_identity_floor %d %.8f"
        % (DEC, 0.5 * 10.0 ** (-DEC)))
    say("  %d draws, scored by the %d smallest primes"
        % (DRAWS, NSMALL))

    say("sieving to %d" % NMAX)
    pr, lam, mu = sieves(NMAX)
    sqf = mu != 0
    del mu
    small = [int(p) for p in pr[:NSMALL]]
    say("  the primes scored: %s" % ", ".join(str(p) for p in small))
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

    rng = np.random.default_rng(SEED)
    say()
    say("  NOTE, disclosed: the first execution failed P1 -- the "
        "minima came out")
    say("  higher than the run's own at both N, so the draws were "
        "not its draws.")
    say("  Their values are not quoted here: that execution's result "
        "file has")
    say("  been overwritten and a number no file carries is not one "
        "to print.")
    say("  The cause is that")
    say("  audit_directidentity_deep.py draws both arms from one "
        "generator, iid")
    say("  first, so the multiplicative masks there come after the "
        "iid arm's")
    say("  consumption. That consumption is replayed here and "
        "nothing else is")
    say("  changed; P1 is what catches it either way.")
    for _ in range(PASSES):
        rng.integers(0, 1 << 32, size=NMAX + 1, dtype=np.uint64)
        rng.integers(0, 1 << 32, size=NMAX + 1, dtype=np.uint64)
    pool = np.zeros((DRAWS, len(NS)))
    ms = np.zeros(DRAWS, dtype=np.int64)
    row = 0
    for k in range(PASSES):
        t = np.zeros(NMAX + 1, dtype=np.uint64)
        hi = rng.integers(0, 1 << 32, size=len(pr), dtype=np.uint64)
        lo = rng.integers(0, 1 << 32, size=len(pr), dtype=np.uint64)
        pm = (hi << np.uint64(32)) | lo
        for i, p in enumerate(pr):
            t[int(p)::int(p)] ^= pm[i]
        smallmask = pm[:NSMALL]
        del hi, lo, pm
        for d in range(PERPASS):
            bit = (t >> np.uint64(d)) & np.uint64(1)
            sg = np.where(sqf, 1.0 - 2.0 * bit.astype(np.float64), 0.0)
            pool[row] = ratio_of(sg)
            ms[row] = int(sum(int((smallmask[j] >> np.uint64(d))
                                  & np.uint64(1))
                              for j in range(NSMALL)))
            row += 1
            del sg, bit
        del t
        say("  pass %d of %d, %d draws" % (k + 1, PASSES, row))

    mins = [float(pool[:, i].min()) for i in range(len(NS))]

    # -------------------------------------------------------------- P1
    say()
    say("P1  are these the same draws?")
    ok = True
    for N, a, b in zip(NS, mins, pub):
        g = abs(round(a, DEC) - round(b, DEC)) < 10.0 ** (-DEC) / 2
        ok &= g
        say("  N = %-8d minimum here %.4f against its %.4f  %s"
            % (N, a, b, "ok" if g else "MISMATCH"))
    say("  P1 %s   (cap: %d decimals at both)"
        % ("hold" if ok else "REFUTED", DEC))
    if not ok:
        io.open(OUT, "w", encoding="utf-8", newline="\n").write(
            "\n".join(HEAD + lines) + "\n")
        raise SystemExit(1)

    say()
    say("      m   draws     ratio min   ratio median   ratio max")
    for m in range(NSMALL + 1):
        sel = ms == m
        if not sel.any():
            continue
        v = pool[sel, 0]
        say("     %2d  %6d  %11.4f  %13.4f  %10.4f"
            % (m, int(sel.sum()), float(v.min()),
               float(np.median(v)), float(v.max())))
        say("POINT floorm_%d %.5f" % (m, float(np.median(v))))
    say("SCALES 1")

    # -------------------------------------------------------------- P2
    say()
    say("P2  does the ratio fall with m?")
    c, se = ols(ms.astype(np.float64), np.log(pool[:, 0]))
    t = c[1] / se[1]
    say("  log(ratio) on m: slope %+.5f +- %.5f, t = %.2f"
        % (c[1], se[1], t))
    say("TSTAT floor_slope %.2f" % t)
    say("SPREAD floor_slope %.5f" % float(ms.max() - ms.min()))
    if abs(t) < 2.0:
        say("UNRESOLVED SIGN floor_slope")
    p2 = c[1] < 0 and abs(t) > 2.0
    say("  P2 %s   (cap: negative and |t| > 2)"
        % ("hold" if p2 else "REFUTED"))

    # -------------------------------------------------------------- P3
    say()
    say("P3  is the floor where m is largest?")
    j = int(np.argmin(pool[:, 0]))
    thr = float(np.quantile(ms.astype(np.float64), 0.9))
    p3 = ms[j] >= thr
    say("  the minimising draw has m = %d; the top decile of m starts "
        "at %.1f" % (ms[j], thr))
    say("  the largest m any draw reaches is %d" % int(ms.max()))
    say("  P3 %s   (cap: in the top decile)"
        % ("hold" if p3 else "REFUTED"))

    # -------------------------------------------------------------- P4
    say()
    say("P4  is it a floor or a tail?")
    hi_ = pool[ms >= thr, 0]
    lo_ = pool[ms <= float(np.quantile(ms.astype(np.float64), 0.1)), 0]
    sh, sl = float(hi_.std(ddof=1)), float(lo_.std(ddof=1))
    say("  top decile of m: %d draws, ratio sd %.4f" % (len(hi_), sh))
    say("  bottom decile:   %d draws, ratio sd %.4f" % (len(lo_), sl))
    say("  ratio of spreads %.4f" % (sh / sl if sl else float("inf")))
    say("POINT floor_spreadratio %.5f" % (sh / sl if sl else 0.0))
    p4 = sl > 0 and sh <= HALF * sl
    say("  P4 %s   (cap: at most %.1f of the bottom decile's)"
        % ("hold" if p4 else "REFUTED", HALF))

    say()
    say("=" * 70)
    say("P1 %s  P2 %s  P3 %s  P4 %s"
        % tuple("hold" if v else "REFUTED"
                for v in (ok, p2, p3, p4)))
    say()
    if p2 and p3 and p4:
        say("the floor is what the algebra says it is. W(v) is "
            "supported on the v")
        say("with at most one negative prime, so a draw that sends "
            "many primes to")
        say("-1 empties it, and those draws are both the lowest and "
            "the tightest.")
        say("What the two earlier remarks called an observation is "
            "now measured.")
    elif not p4:
        say("the low end is a tail and not a floor, so the reading "
            "rem:identityseeds")
        say("and rem:identitybig offered is wrong and the stability "
            "they measured")
        say("needs another explanation. The spreads are above.")
    else:
        say("the derived mechanism is a spectator: many negative "
            "primes do empty W,")
        say("and something else is choosing the minimum. The "
            "minimising draw's m")
        say("is printed above and is where to look.")

    io.open(OUT, "w", encoding="utf-8", newline="\n").write(
        "\n".join(HEAD + lines) + "\n")
    print("\nwrote %s" % os.path.normpath(OUT))
    return 0


if __name__ == "__main__":
    sys.exit(main())
