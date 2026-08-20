# -*- coding: utf-8 -*-
r"""
The object item 5 says has no name is a truncated Möbius sieve weight

WHAT IS AT STAKE

OPEN item 5 ends its description of what the sign axis is left with by
saying, in as many words, that the object **has no name**: "a
correlation of Lambda against a convolution supported on indices with
two or more distinct prime factors, restricted to the tenth of the
dilations whose sums are largest."  rem:deficitlog then measured that
this field cannot separate a logarithm from a line (corr 0.996737), so
the shape of G cannot be settled by fitting here and what is left is
an unconditional statement about |sum a| itself.  An unconditional
statement needs the object to be one that unconditional statements
exist about.

The algebra gives it a name.  rem:denominator has
sum a = sum_j Lambda(N-j) LK(j) with

    LK(j) = sum_{k | j, k squarefree, k < N^theta, (k,N)=1,
                 (j/k,k)=1} mu(j/k) log k,

and rem:support splits it: j prime below the level gives a Goldbach
count, j a higher prime power gives exactly zero, everything else has
two or more prime factors.  For squarefree j coprime to N the
conditions on k are automatic and the *untruncated* sum is mu*log =
Lambda(j) = 0 as soon as omega(j) >= 2.  **So the composite part is
nothing but the divisors the truncation threw away**, and rewriting
those by their cofactors d = j/k, which satisfy d <= j/N^theta,

    LK(j) = - sum_{d | j, d <= j/N^theta} mu(d) log(j/d).

That is a truncated Möbius sum over the divisors of j below a level --
the Eratosthenes-Legendre sieve weight, cut at D_j = j/N^theta.  Two
things follow if it is right.  The composite part is supported on
j >= N^theta, since a smaller j has no divisor above the level at all
and the sum is empty.  And |sum a| is a correlation of Lambda against
a classical sieve weight, which is a family that unconditional bounds
exist for, rather than an object with no name.

This run does not fit anything.  It checks an identity.

BACKS: Remark {#rem:sieveweight} in paper/wall_v3.md.

PRE-REGISTERED PREDICTIONS (written before this script was run)

  H1  THE GATE.  Built from the (k,m) definition, sum a reproduces
      the published |sum a| at N = 200000, the POINT marker of
      results/audit_deficit_direct.txt, to four decimals.

      *Disclosed.*  H1 as first written named rem:denominator's T4
      row "trunc only" -- -8657.0, -17462.3, -23198.0 -- as the
      published restricted sum.  It is not: that row applies the
      truncation *without* the other three restrictions, so it is a
      different quantity from sum a.  The error was found when the
      script failed to parse that table and exited, before any
      quantity had been computed; no verdict existed to be changed.
      The gate now reads a marker that is sum a itself.
  H2  **THE IDENTITY.**  For every squarefree j < N coprime to N with
      omega(j) >= 2, LK(j) equals -sum_{d|j, d <= j/N^theta} mu(d)
      log(j/d) to a relative 1e-9.
  H3  **THE SUPPORT.**  No j < N^theta contributes to the composite
      part: the composite sum restricted to j < N^theta is exactly
      zero.
  H4  Rebuilding sum a as (prime part) + (sieve-weight part) matches
      the direct sum a to a relative 1e-12.
  H5  The sieve weight is not a small correction to something else:
      the composite part carries at least 0.99 of |sum a| at every N
      here, as rem:support's exponents (+0.716454 against +0.717916)
      imply.

REFUTATION RULE (fixed before the run)

  H1  REFUTED outside four decimals; nothing below is reported, since
      a construction that does not reproduce the published sum is not
      the published object.  The substitution of the marker for the
      T4 row is disclosed above and changes what is read, not what is
      required.
  H2  **REFUTED by any j failing.**  Then the cofactor rewriting is
      wrong and the object is not the sieve weight this remark would
      name it.  The failure is to be reported with the smallest
      failing j, because that names which condition was dropped --
      the (k,N)=1 restriction is the one this derivation waves at, and
      a failure confined to j sharing a factor with N would say the
      identity holds on the coprime part and needs a correction term
      elsewhere.  **That case is named here so it cannot be presented
      afterwards as agreement.**
  H3  REFUTED by any nonzero contribution below the level.
  H4  REFUTED outside 1e-12 relative.
  H5  REFUTED below 0.99. Then the composite part is not the object
      and naming it names the wrong thing.

  A SECOND BLOCK, REGISTERED AFTER H2-H4 WERE REFUTED AND SAYING SO

  H2, H3 and H4 are refuted above and stay refuted; the run that
  refuted them is not re-labelled.  Two causes were found afterwards
  and each is a defect in what this script asked, not in what it
  measured:

    (a) **The level is an integer.**  audit_lean_identity.py takes
        K = int(N**theta) and k in range(2, K), so the truncation is
        at int(N^theta) and not at N^theta.  At N = 50000 those are
        427 and 427.979732, and j = 16653 = 3*7*13*61 has the divisor
        427 sitting in the gap.  H2 was written with the real level
        and is refuted with it.
    (b) **H3 and H4 tested the wrong set.**  The derivation is for
        *squarefree* j; H3 and H4 ranged over every j with Lambda(j)
        = 0, which includes j like 12 = 2^2*3 where mu*log is not
        Lambda and the cofactor rewriting was never claimed.

  So, pre-registered before the second run and after seeing the first:

  H6  With the level taken as int(N^theta), the identity holds on
      squarefree j coprime to N with omega(j) >= 2, to a relative
      1e-9 at every such j.
  H7  On that same set, no j below int(N^theta) contributes: the sum
      of |LK| over squarefree composite j below the level is exactly
      zero.
  H8  Rebuilding sum a as the prime part, plus the sieve weight on
      squarefree composite j, plus LK itself on the j the derivation
      does not cover, matches sum a to a relative 1e-12.
  H9  **And the name has to cover the object.**  The squarefree
      composite j carry at least 0.90 of the composite part; if the j
      the derivation does not cover carry more than a tenth, calling
      the object a sieve weight names a piece and not the thing.

  REFUTATION for the second block.  H6 refuted by any failing j, and
  then the cofactor rewriting is wrong for reasons that are not the
  level and the naming fails outright.  H7 by any nonzero sum.  H8
  outside 1e-12.  **H9 refuted below 0.90, and that is the one that
  would matter**: the identity could hold on squarefree j and still
  describe a minority of the object, in which case this remark may
  not say the object has a name -- only that a part of it does.

  WHAT THIS CANNOT DO.  An identity is not a bound.  Naming the
  object as a sieve weight does not supply an unconditional estimate
  for the correlation, and nothing here says the classical bounds are
  strong enough to decide item 5 -- only that the object is in a
  family where such bounds are the right thing to look for.  No
  exponent is measured here and no forecast is made.
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
OUT = os.path.join(ROOT, "results", "audit_sieve_weight.txt")
SRC = os.path.join(ROOT, "results",
                   "audit_deficit_direct.txt")

THETA = 0.56
NS = [20_000, 50_000, 100_000]
NGATE = 200_000
DEC = 4
RELID = 1e-9
RELSUM = 1e-12
SHARE = 0.99
SHARE2 = 0.90


def primes_upto(n):
    s = np.ones(n + 1, dtype=bool)
    s[:2] = False
    for p in range(2, int(n ** 0.5) + 1):
        if s[p]:
            s[p * p::p] = False
    return np.flatnonzero(s).astype(np.int64)


def lambda_and_mu(n):
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
    rem = np.arange(n + 1, dtype=np.int64)
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


def lk_direct(N, mu, sqf):
    """LK(j) for every j < N, from the (k,m) definition itself"""
    PN = factor_set(N)
    K = int(N ** THETA)
    out = np.zeros(N, dtype=np.float64)
    for k in range(2, K):
        if not sqf[k] or any(k % q == 0 for q in PN):
            continue
        lg = math.log(k)
        ms = np.arange(1, (N - 1) // k + 1, dtype=np.int64)
        for q in factor_set(k):
            ms = ms[ms % q != 0]
        out[ms * k] += lg * mu[ms].astype(np.float64)
        del ms
    return out


def lk_sieve(N, mu, lev):
    """- sum_{d | j, d <= j/lev} mu(d) log(j/d), for every j < N"""
    out = np.zeros(N, dtype=np.float64)
    dmax = int(N / lev) + 2
    for d in range(1, dmax + 1):
        m = int(mu[d])
        if m == 0:
            continue
        lo = int(math.ceil(d * lev))
        lo = max(lo, d)
        if lo >= N:
            continue
        js = np.arange(((lo + d - 1) // d) * d, N, d, dtype=np.int64)
        js = js[js >= d * lev]
        if js.size:
            out[js] -= m * np.log(js.astype(np.float64) / d)
        del js
    return out


def read_pub():
    """|sum a| at NGATE, from the POINT marker that is sum a itself"""
    src = io.open(SRC, encoding="utf-8").read()
    m = re.search(r"^POINT deficitdirect_%d ([\d.eE+-]+) " % NGATE,
                  src, re.M)
    if not m:
        raise SystemExit("no POINT marker for N = %d" % NGATE)
    return float(m.group(1))


HEAD = [
    "STATISTIC: the identity LK(j) = - sum_{d|j, d <= j/N^theta} mu(d)",
    "           log(j/d) on the composite part of sum a, its support,",
    "           and the rebuild of sum a from the prime part and the",
    "           truncated Mobius sieve weight.",
    "FIELD: N = %s, the three N of rem:denominator's T4; j over every"
    % NS,
    "       index below N; k over the squarefree k < N^%.2f coprime to"
    % THETA,
    "       N, the k-range of code/audit_lean_identity.py. The",
    "       published truncation-alone values are READ from",
    "       results/audit_denominator.txt.",
    "DERIVED: for squarefree j coprime to N with omega(j) >= 2 the",
    "         untruncated mu*log is Lambda(j) = 0, so LK(j) is exactly",
    "         minus what the truncation removed, rewritten by",
    "         cofactors d = j/k <= j/N^theta.",
    "",
]


def main():
    lines = []

    def say(t=""):
        print(t)
        sys.stdout.flush()
        lines.append(t)

    pub = read_pub()
    say("READ audit_deficit_direct.txt %d %.4f" % (NGATE, pub))
    say("  |sum a| at the gate N, from the marker that is sum a "
        "itself")
    say("  NOTE, disclosed: H1 first named rem:denominator's T4 "
        "'trunc only'")
    say("  row as the published restricted sum. It is not -- that "
        "row drops the")
    say("  other three restrictions. Found when the parse failed, "
        "before any")
    say("  quantity had been computed; no verdict existed to change.")
    say("PRINTBOUND audit_sieve_weight %d %.8f"
        % (DEC, 0.5 * 10.0 ** (-DEC)))
    say("  theta %.2f, identity tolerance %.0e, rebuild %.0e, share "
        "%.2f" % (THETA, RELID, RELSUM, SHARE))

    say("sieving to %d" % max(NS + [NGATE]))
    lam, mu = lambda_and_mu(max(NS + [NGATE]))
    sqf = mu != 0

    rows = []
    for N in NS:
        lk = lk_direct(N, mu, sqf)
        sa = float((lam[N - np.arange(1, N)] * lk[1:]).sum())
        rows.append((N, lk, sa))
        say("  N = %-8d sum a = %.4f" % (N, sa))

    # -------------------------------------------------------------- H1
    say()
    say("H1  the gate: does this reproduce rem:denominator's sum a?")
    lkg = lk_direct(NGATE, mu, sqf)
    sag = float((lam[NGATE - np.arange(1, NGATE)] * lkg[1:]).sum())
    del lkg
    h1 = abs(round(abs(sag), DEC) - round(pub, DEC)) < 10.0 ** (-DEC)
    say("  N = %-8d here %.4f against its %.4f  %s"
        % (NGATE, abs(sag), pub, "ok" if h1 else "MISMATCH"))
    say("  H1 %s   (cap: the published decimal)"
        % ("hold" if h1 else "REFUTED"))
    if not h1:
        io.open(OUT, "w", encoding="utf-8", newline="\n").write(
            "\n".join(HEAD + lines) + "\n")
        raise SystemExit(1)

    # -------------------------------------------------------------- H2
    say()
    say("H2  the identity, on squarefree j coprime to N with omega>=2")
    h2 = True
    worst = 0.0
    badj = None
    for N, lk, sa in rows:
        sv = lk_sieve(N, mu, N ** THETA)
        j = np.arange(1, N, dtype=np.int64)
        keep = sqf[1:N].copy()
        for q in sorted(factor_set(N)):
            keep &= (j % int(q)) != 0
        keep &= lam[1:N] == 0.0                    # not 1, not a prime
        keep &= j > 1
        jj = j[keep]
        a = lk[jj]
        b = sv[jj]
        den = np.maximum(np.abs(a), 1.0)
        rel = np.abs(a - b) / den
        w = float(rel.max()) if rel.size else 0.0
        if rel.size and w > worst:
            worst = w
            badj = int(jj[int(np.argmax(rel))])
        say("  N = %-8d %-7d such j, worst relative %.3e"
            % (N, jj.size, w))
        del sv, j, keep, jj, a, b, den, rel
    h2 = worst <= RELID
    say("POINT identityworst %.6e" % worst)
    if not h2:
        say("  the smallest failing j is %s, which names the dropped "
            "condition" % badj)
    say("  H2 %s   (cap: %.0e relative)"
        % ("hold" if h2 else "REFUTED", RELID))

    # -------------------------------------------------------------- H3
    say()
    say("H3  is the composite part empty below the level?")
    h3 = True
    for N, lk, sa in rows:
        K = N ** THETA
        j = np.arange(1, N, dtype=np.int64)
        comp = (lam[1:N] == 0.0) & (j > 1)
        below = comp & (j < K)
        s = float(np.abs(lk[j[below]]).sum())
        h3 &= s == 0.0
        say("  N = %-8d %-6d composite j below %.0f, |LK| summing "
            "to %.1f" % (N, int(below.sum()), K, s))
        del j, comp, below
    say("  H3 %s   (cap: exactly zero)"
        % ("hold" if h3 else "REFUTED"))

    # -------------------------------------------------------------- H4
    say()
    say("H4  rebuild: prime part plus sieve-weight part")
    h4 = True
    for N, lk, sa in rows:
        sv = lk_sieve(N, mu, N ** THETA)
        j = np.arange(1, N, dtype=np.int64)
        isp = lam[1:N] > 0.0
        pr = float((lam[N - j[isp]] * lk[j[isp]]).sum())
        cm = ~isp & (j > 1)
        co = float((lam[N - j[cm]] * sv[j[cm]]).sum())
        rel = abs((pr + co) - sa) / max(abs(sa), 1.0)
        h4 &= rel <= RELSUM
        say("  N = %-8d prime %+.4f  sieve %+.4f  total %+.4f  "
            "relative %.2e" % (N, pr, co, pr + co, rel))
        say("POINT rebuild_%d %.6e" % (N, rel))
        del sv, j, isp, cm
    say("  H4 %s   (cap: %.0e relative)"
        % ("hold" if h4 else "REFUTED", RELSUM))

    # -------------------------------------------------------------- H5
    say()
    say("H5  does the composite part carry the object?")
    h5 = True
    for N, lk, sa in rows:
        j = np.arange(1, N, dtype=np.int64)
        isp = lam[1:N] > 0.0
        pr = float((lam[N - j[isp]] * lk[j[isp]]).sum())
        cm = ~isp & (j > 1)
        co = float((lam[N - j[cm]] * lk[j[cm]]).sum())
        sh = abs(co) / max(abs(sa), 1e-30)
        h5 &= sh >= SHARE
        say("  N = %-8d composite %+.4f is %.4f of |sum a|, prime "
            "%+.4f" % (N, co, sh, pr))
        say("POINT compshare_%d %.6f" % (N, sh))
        del j, isp, cm
    say("SCALES %d" % len(rows))
    say("  H5 %s   (cap: %.2f)"
        % ("hold" if h5 else "REFUTED", SHARE))

    # ---------------------------------------- the second block
    say()
    say("  the two causes found after H2-H4, and the corrected checks")
    ngap = 50_000
    say("  (a) the level is int(N^theta), not N^theta: at N = %d"
        % ngap)
    say("      those are %d and %.6f, and j = 16653 = 3*7*13*61 has"
        % (int(ngap ** THETA), ngap ** THETA))
    say("      the divisor 427 sitting in the gap")
    say("  (b) H3 and H4 ranged over every j with Lambda(j) = 0, "
        "which includes")
    say("      non-squarefree j the derivation never covered")

    say()
    say("H6  the identity with the integer level, on squarefree j")
    h6 = True
    w6 = 0.0
    bad6 = None
    sets = {}
    for N, lk, sa in rows:
        lev = float(int(N ** THETA))
        sv = lk_sieve(N, mu, lev)
        j = np.arange(1, N, dtype=np.int64)
        keep = sqf[1:N].copy()
        for q in sorted(factor_set(N)):
            keep &= (j % int(q)) != 0
        keep &= lam[1:N] == 0.0
        keep &= j > 1
        jj = j[keep]
        sets[N] = (jj, sv, lev)
        rel = np.abs(lk[jj] - sv[jj]) / np.maximum(np.abs(lk[jj]), 1.0)
        w = float(rel.max()) if rel.size else 0.0
        if w > w6:
            w6, bad6 = w, int(jj[int(np.argmax(rel))])
        say("  N = %-8d level %d, %-7d such j, worst relative %.3e"
            % (N, int(lev), jj.size, w))
        del j, keep, rel
    h6 = w6 <= RELID
    say("POINT identityint %.6e" % w6)
    if not h6:
        say("  the smallest failing j is %s" % bad6)
    say("  H6 %s   (cap: %.0e relative)"
        % ("hold" if h6 else "REFUTED", RELID))

    say()
    say("H7  is that set empty below the integer level?")
    h7 = True
    for N, lk, sa in rows:
        jj, sv, lev = sets[N]
        below = jj[jj < lev]
        t = float(np.abs(lk[below]).sum())
        h7 &= t == 0.0
        say("  N = %-8d %-5d squarefree composite j below %d, |LK| "
            "summing to %.4e" % (N, below.size, int(lev), t))
        say("POINT belowlevel_%d %.6e" % (N, t))
    say("  H7 %s   (cap: exactly zero)"
        % ("hold" if h7 else "REFUTED"))
    if not h7:
        say("  NOTE, disclosed: the cap was written as an exact zero "
            "on a sum of")
        say("  floating-point magnitudes, which no computation can "
            "meet. The sums")
        say("  above are the number the cap should have been "
            "compared with. The")
        say("  rule is not rewritten; the cap was unmeetable as "
            "written and that")
        say("  is recorded as the defect it is.")

    say()
    say("H8  rebuild with the sieve weight where it applies")
    h8 = True
    for N, lk, sa in rows:
        jj, sv, lev = sets[N]
        j = np.arange(1, N, dtype=np.int64)
        isp = lam[1:N] > 0.0
        pr = float((lam[N - j[isp]] * lk[j[isp]]).sum())
        cov = np.zeros(N, dtype=bool)
        cov[jj] = True
        rest = (~isp) & (j > 1) & (~cov[1:N])
        sw = float((lam[N - jj] * sv[jj]).sum())
        rs = float((lam[N - j[rest]] * lk[j[rest]]).sum())
        rel = abs((pr + sw + rs) - sa) / max(abs(sa), 1.0)
        h8 &= rel <= RELSUM
        say("  N = %-8d prime %+.4f  sieve %+.4f  uncovered %+.4f  "
            "relative %.2e" % (N, pr, sw, rs, rel))
        del j, isp, cov, rest
    say("  H8 %s   (cap: %.0e relative)"
        % ("hold" if h8 else "REFUTED", RELSUM))

    say()
    say("H9  does the name cover the object?")
    h9 = True
    for N, lk, sa in rows:
        jj, sv, lev = sets[N]
        j = np.arange(1, N, dtype=np.int64)
        isp = lam[1:N] > 0.0
        cov = np.zeros(N, dtype=bool)
        cov[jj] = True
        comp = (~isp) & (j > 1)
        co = float((lam[N - j[comp]] * lk[j[comp]]).sum())
        sw = float((lam[N - jj] * lk[jj]).sum())
        sh = abs(sw) / max(abs(co), 1e-30)
        h9 &= sh >= SHARE2
        say("  N = %-8d squarefree composite carry %+.4f of the "
            "composite %+.4f, share %.4f" % (N, sw, co, sh))
        say("POINT sqfshare_%d %.6f" % (N, sh))
        del j, isp, cov, comp
    say("  H9 %s   (cap: %.2f)"
        % ("hold" if h9 else "REFUTED", SHARE2))

    say()
    say("=" * 70)
    say("H6 %s  H7 %s  H8 %s  H9 %s"
        % tuple("hold" if v else "REFUTED" for v in (h6, h7, h8, h9)))
    say("H1 %s  H2 %s  H3 %s  H4 %s  H5 %s"
        % tuple("hold" if v else "REFUTED"
                for v in (h1, h2, h3, h4, h5)))
    say()
    if h6 and h7 and h8 and h9:
        say("the object has a name. the composite part of sum a -- "
            "which is the")
        say("whole of it -- is the correlation of Lambda against the "
            "truncated")
        say("Mobius sieve weight of Eratosthenes and Legendre, cut at "
            "D_j = j/N^theta,")
        say("and it is supported on j >= N^theta because a smaller j "
            "has no divisor")
        say("above the level to throw away. an unconditional bound "
            "for item 5 is")
        say("a bound on a sieve-weighted Goldbach correlation, which "
            "is a thing")
        say("that has a literature. that is what this identity buys "
            "and it is all")
        say("it buys: an identity is not a bound.")
    else:
        say("the rewriting is not the identity this run claimed, so "
            "the object is")
        say("not named here and item 5 keeps the description it had.")

    io.open(OUT, "w", encoding="utf-8", newline="\n").write(
        "\n".join(HEAD + lines) + "\n")
    print("\nwrote %s" % os.path.normpath(OUT))
    return 0


if __name__ == "__main__":
    sys.exit(main())
