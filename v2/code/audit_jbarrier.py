# -*- coding: utf-8 -*-
r"""
The same sum has two square-root barriers, and nobody measured the second

WHAT IS AT STAKE

rem:sieveweight named the object: sum a = sum_j Lambda(N-j) LK(j) with
LK the Eratosthenes-Legendre sieve weight cut at D_j = j/K,
K = floor(N^theta), supported on j >= K.  Item 5's demand is
|sum a| <~ l2, and l2 is a norm **over k** -- the dilations -- while
sum a is a sum **over j**.  The same number has two groupings and
therefore two square-root barriers, and only one of them has ever
been computed.

    l2 = sqrt( sum_k (log k)^2 H(N;k)^2 )        over k, published
    D  = sqrt( sum_j Lambda(N-j)^2 LK(j)^2 )     over j, never computed

Each is what the sum would be if its terms carried random signs in
that grouping.  Which one is the target matters, because the demand
is a comparison against l2 and l2 is not obviously the relevant floor
for a sum indexed by j.

**If e(D) = e(l2), the demand is exactly "achieve square-root
cancellation"** -- hard by every known method, not known impossible,
and the honest description of item 5 changes from a gap between two
fitted exponents to a named barrier.  **If e(D) is resolvedly above
e(l2), the demand asks for less than square-root cancellation over
j**, which no method supplies, and item 5 would be closed against the
sign axis rather than left open.  Both readings are consequences and
neither is assumed here.

rem:denominator wrote that "sum_k (log k) H(N;k) is already as small
as chance allows".  That is a statement about the k grouping.  Nothing
in this repository has said what chance allows in the j grouping, and
the two need not agree.

BACKS: Remark {#rem:jbarrier} in paper/wall_v3.md.

PRE-REGISTERED PREDICTIONS (written before this script was run)

  J1  THE GATE.  At N = 200000 this reproduces the published
      |sum a| = 87895.3236, the POINT marker of
      results/audit_deficit_direct.txt, and the published
      l1/l2 = 11.9596 of results/audit_lean_identity.txt, both to four
      decimals.
  J2  **THE ONE THAT DECIDES THE READING.**  e(D) and e(l2) agree:
      their difference is smaller than 0.03 in exponent.  The two
      groupings of one sum have the same square-root scale.
  J3  The sum stands above both by the same amount: the exponent of
      |sum a|/D is positive and within 0.03 of the exponent of
      |sum a|/l2, which the field publishes as +0.134019.
  J4  And it stands above them everywhere, not on average: |sum a| > D
      and |sum a| > l2 at every N in this field.

REFUTATION RULE (fixed before the run)

  J1  REFUTED outside four decimals on either; nothing below is
      reported, because a construction that does not reproduce the
      published sum is not the published object.
  J2  REFUTED outside 0.03.  **If it is refuted with D the larger,
      that is the outcome that matters most and must be stated as
      such**: the demand would be asking |sum a| below what random
      signs give in its own index, and no method known to this
      repository or to the literature it cites delivers below
      square-root cancellation.  If it is refuted with l2 the larger,
      the demand is weaker than square-root and the barrier is not
      the obstruction.
      **The unresolved case is named**: nine N over 2.41 decades give
      exponents with standard errors this run prints, and a gap
      smaller than the larger of the two errors is not a gap.  In
      that case J2 is to be read as "not resolved either way" and
      neither consequence above may be drawn -- the verdict word is
      still whatever the 0.03 cap says, but the reading is barred.
  J3  REFUTED outside 0.03, or by a negative exponent.
  J4  REFUTED by any N where either inequality fails.

  A SECOND BLOCK, REGISTERED AFTER J2 AND SAYING SO

  J2 holds at its 0.03 cap and its unresolved clause fires: the gap
  -0.000261 is inside the larger individual standard error 0.005519,
  so by the rule fixed above **no reading is drawn from J2 and that
  stands**.  But the rule used the wrong error bar, and this is the
  M9 family again.  The two exponents are fitted on the *same* nine
  N, so their errors are strongly correlated and the error of the
  gap is not the larger of the two -- it is the error of a fit to the
  ratio, which this run did not print.  A rule that compares a
  difference to an error bar the difference does not have cannot
  detect agreement, only fail to detect disagreement.

  And the pointwise ratios say more than the exponents did: D/l2 runs
  0.485 to 0.511 across the field.  So, pre-registered before the
  second run and after seeing the first:

  K1  Fitting log(D/l2) directly on log N, the slope's standard error
      is smaller than either individual fit's, and the slope is
      unresolved: |t| < 2.  That is the statistic J2 should have
      used.
  K2  The ratio is flat, not merely equal in exponent: its range
      across the field is under 10 per cent of its mean.
  K3  **And its value is one half**: the mean of D/l2 is within 0.02
      of 0.5.

  REFUTATION for the second block.  K1 refuted at |t| >= 2, and then
  the two barriers do drift apart and the sign of the drift decides
  which of J2's two consequences applies -- D rising means the demand
  asks below square-root cancellation in j, l2 rising means it asks
  less.  K2 refuted above 10 per cent, and then "flat" is the wrong
  word and only the exponent statement survives.  K3 refuted outside
  0.02, and then the constant is not one half and this run reports
  the number it is without proposing where it comes from -- **no
  mechanism is offered for a constant this run only measured**, and
  a near miss is a miss.

  WHAT THIS CANNOT DO.  D is a heuristic floor, not a theorem: it is
  what the sum would be under random signs in the j grouping, and
  measuring it does not prove any lower bound on |sum a|.  Nine N on
  one radical over 2.41 decades.  No forecast is made and no closure
  N is computed; rem:shapepower and rem:deficitlog stand.
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
OUT = os.path.join(ROOT, "results", "audit_jbarrier.txt")
SRCA = os.path.join(ROOT, "results", "audit_deficit_direct.txt")
SRCB = os.path.join(ROOT, "results", "audit_lean_identity.txt")

THETA = 0.56
NS = [25_000 * (1 << j) for j in range(9)]
NGATE = 200_000
DEC = 4
GAPCAP = 0.03
PUBDEF = 0.134019
TFLAT = 2.0
RANGEPC = 10.0
HALFCAP = 0.02
HALF = 0.5


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


def both_groupings(N, lam, mu, sqf):
    """a_k over k and LK(j) over j, in one pass over the k-range"""
    PN = factor_set(N)
    K = int(N ** THETA)
    lk = np.zeros(N, dtype=np.float64)
    l1 = l2sq = 0.0
    for k in range(2, K):
        if not sqf[k] or any(k % q == 0 for q in PN):
            continue
        lg = math.log(k)
        ms = np.arange(1, (N - 1) // k + 1, dtype=np.int64)
        for q in factor_set(k):
            ms = ms[ms % q != 0]
        mm = mu[ms].astype(np.float64)
        a = lg * float((lam[N - ms * k] * mm).sum())
        l1 += abs(a)
        l2sq += a * a
        lk[ms * k] += lg * mm
        del ms, mm
    j = np.arange(1, N, dtype=np.int64)
    w = lam[N - j] * lk[1:]
    sa = float(w.sum())
    dsq = float((w * w).sum())
    del j, w, lk
    return l1, math.sqrt(l2sq), abs(sa), math.sqrt(dsq)


def fit(x, y):
    x = np.asarray(x, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)
    b, a0 = np.polyfit(x, y, 1)
    r = y - (b * x + a0)
    se = math.sqrt(float((r ** 2).sum() / (len(x) - 2))
                   / float(((x - x.mean()) ** 2).sum()))
    return float(b), se


def read_pub():
    sa = re.search(r"^POINT deficitdirect_%d ([\d.eE+-]+) " % NGATE,
                   io.open(SRCA, encoding="utf-8").read(), re.M)
    ll = re.search(r"N = %d\s+#k = \d+\s+G = [\d.]+\s+l1/l2 = ([\d.]+)"
                   % NGATE, io.open(SRCB, encoding="utf-8").read())
    if not sa or not ll:
        raise SystemExit("missing a published value for N = %d" % NGATE)
    return float(sa.group(1)), float(ll.group(1))


HEAD = [
    "STATISTIC: the two square-root barriers of one sum -- l2 over the",
    "           dilations k and D = sqrt(sum_j Lambda(N-j)^2 LK(j)^2)",
    "           over the indices j -- their exponents, and where",
    "           |sum a| stands against each.",
    "FIELD: N = 25000*2^j for j < %d, one radical throughout, 2.41"
    % len(NS),
    "       decades; k over the squarefree k < N^%.2f coprime to N, the"
    % THETA,
    "       k-range of code/audit_lean_identity.py; j over every index",
    "       below N. |sum a| and l1/l2 at N = %d are READ from" % NGATE,
    "       results/audit_deficit_direct.txt and",
    "       results/audit_lean_identity.txt as the gate.",
    "",
]


def main():
    lines = []

    def say(t=""):
        print(t)
        sys.stdout.flush()
        lines.append(t)

    psa, pll = read_pub()
    say("READ audit_deficit_direct.txt %d %.4f" % (NGATE, psa))
    say("READ audit_lean_identity.txt %d %.4f" % (NGATE, pll))
    say("  |sum a| and l1/l2 at the gate N")
    say("PRINTBOUND audit_jbarrier %d %.8f"
        % (DEC, 0.5 * 10.0 ** (-DEC)))
    say("  theta %.2f, gap cap %.2f, published deficit %.6f"
        % (THETA, GAPCAP, PUBDEF))

    NMAX = max(NS + [NGATE])
    say("sieving to %d" % NMAX)
    lam, mu = lambda_and_mu(NMAX)
    sqf = mu != 0

    rows = []
    for N in NS:
        l1, l2, sa, d = both_groupings(N, lam, mu, sqf)
        rows.append((N, l1, l2, sa, d))
        say("  N = %-9d |sum a| %12.2f   l2 %11.2f   D %11.2f"
            % (N, sa, l2, d))
        say("POINT jbarrier_%d %.6e %.6e" % (N, d, l2))
    say("SCALES %d" % len(rows))

    # -------------------------------------------------------------- J1
    say()
    say("J1  the gate: does this reproduce the published pair?")
    g = [r for r in rows if r[0] == NGATE]
    if not g:
        l1g, l2g, sag, dg = both_groupings(NGATE, lam, mu, sqf)
    else:
        _, l1g, l2g, sag, dg = g[0]
    a = abs(round(sag, DEC) - round(psa, DEC)) < 10.0 ** (-DEC)
    b = abs(round(l1g / l2g, DEC) - round(pll, DEC)) < 10.0 ** (-DEC)
    j1 = a and b
    say("  |sum a| here %.4f against its %.4f  %s"
        % (sag, psa, "ok" if a else "MISMATCH"))
    say("  l1/l2   here %.4f against its %.4f  %s"
        % (l1g / l2g, pll, "ok" if b else "MISMATCH"))
    say("  J1 %s   (cap: %d decimals on both)"
        % ("hold" if j1 else "REFUTED", DEC))
    if not j1:
        io.open(OUT, "w", encoding="utf-8", newline="\n").write(
            "\n".join(HEAD + lines) + "\n")
        raise SystemExit(1)

    x = np.array([math.log(r[0]) for r in rows])
    esa, ssa = fit(x, [math.log(r[3]) for r in rows])
    el2, sl2 = fit(x, [math.log(r[2]) for r in rows])
    ed, sd = fit(x, [math.log(r[4]) for r in rows])
    say()
    say("        quantity        exponent        s.e.")
    for nm, e, s in (("|sum a|", esa, ssa), ("l2 over k", el2, sl2),
                     ("D over j", ed, sd)):
        say("  %-14s   %+.6f     %.6f" % (nm, e, s))
    say("TSTAT jbarrier_gap %.2f"
        % ((ed - el2) / max(sd, sl2)))
    say("SPREAD jbarrier_gap %.6f" % max(sd, sl2))

    # -------------------------------------------------------------- J2
    say()
    say("J2  do the two barriers have the same exponent?")
    gap = ed - el2
    big = max(sd, sl2)
    j2 = abs(gap) < GAPCAP
    say("  e(D) - e(l2) = %+.6f, against the larger error %.6f"
        % (gap, big))
    say("POINT barriergap %.6f" % gap)
    say("  J2 %s   (cap: %.2f)"
        % ("hold" if j2 else "REFUTED", GAPCAP))
    if abs(gap) < big:
        say("  UNRESOLVED: the gap is inside the larger standard "
            "error, so neither")
        say("  reading may be drawn, as the rule says")

    # -------------------------------------------------------------- J3
    say()
    say("J3  does the sum stand above both by the same amount?")
    gd = esa - ed
    gl = esa - el2
    j3 = gd > 0 and abs(gd - gl) < GAPCAP
    say("  e(|sum a|/D)  = %+.6f" % gd)
    say("  e(|sum a|/l2) = %+.6f, published %.6f" % (gl, PUBDEF))
    say("  difference %+.6f" % (gd - gl))
    say("POINT abovebarrier %.6f" % gd)
    say("  J3 %s   (cap: positive and within %.2f)"
        % ("hold" if j3 else "REFUTED", GAPCAP))

    # -------------------------------------------------------------- J4
    say()
    say("J4  does it stand above them at every N?")
    badd = [r[0] for r in rows if not r[3] > r[4]]
    badl = [r[0] for r in rows if not r[3] > r[2]]
    j4 = not badd and not badl
    say("  |sum a| <= D at %s" % (badd if badd else "no N"))
    say("  |sum a| <= l2 at %s" % (badl if badl else "no N"))
    say("  J4 %s   (cap: every N)" % ("hold" if j4 else "REFUTED"))

    # ---------------------------------------- the second block
    say()
    say("  registered after J2: its rule compared a difference to an")
    say("  error bar the difference does not have, since both fits "
        "use the")
    say("  same nine N and their errors are correlated")

    say()
    say("K1  the statistic J2 should have used")
    ratio = np.array([r[4] / r[2] for r in rows])
    er, ser = fit(x, np.log(ratio))
    t1 = er / ser
    k1 = abs(t1) < TFLAT and ser < min(sd, sl2)
    say("  slope of log(D/l2) on log N: %+.6f +- %.6f, t = %+.2f"
        % (er, ser, t1))
    say("  against the individual errors %.6f and %.6f" % (sd, sl2))
    say("TSTAT jbarrier_ratio %.2f" % t1)
    say("SPREAD jbarrier_ratio %.6f" % ser)
    say("  K1 %s   (cap: |t| < %.1f and the smaller error)"
        % ("hold" if k1 else "REFUTED", TFLAT))

    say()
    say("K2  is the ratio flat?")
    mn = float(ratio.mean())
    rng = 100.0 * (float(ratio.max()) - float(ratio.min())) / mn
    k2 = rng < RANGEPC
    say("    N            D/l2")
    for r, v in zip(rows, ratio):
        say("  %-11d %.6f" % (r[0], v))
    say("  mean %.6f, range %.2f per cent" % (mn, rng))
    say("POINT ratiomean %.6f" % mn)
    say("SERIES jbarrier_ratio %s"
        % " ".join("%.4f" % v for v in ratio))
    say("FLAT jbarrier_ratio" if k2 else "DRIFTS jbarrier_ratio")
    say("  K2 %s   (cap: %.1f per cent)"
        % ("hold" if k2 else "REFUTED", RANGEPC))

    say()
    say("K3  is it one half?")
    k3 = abs(mn - HALF) <= HALFCAP
    say("  mean %.6f against %.1f, off by %+.6f" % (mn, HALF, mn - HALF))
    say("  K3 %s   (cap: %.2f)"
        % ("hold" if k3 else "REFUTED", HALFCAP))
    if k3:
        say("  no mechanism is offered for this constant; it is "
            "measured here and")
        say("  not derived, as the rule for K3 says")

    say()
    say("=" * 70)
    say("K1 %s  K2 %s  K3 %s"
        % tuple("hold" if v else "REFUTED" for v in (k1, k2, k3)))
    say("J1 %s  J2 %s  J3 %s  J4 %s"
        % tuple("hold" if v else "REFUTED" for v in (j1, j2, j3, j4)))
    say()
    if j2 and abs(gap) >= big:
        say("the two groupings of one sum have the same square-root "
            "scale, so the")
        say("demand of item 5 is exactly 'achieve square-root "
            "cancellation' -- a named")
        say("barrier rather than a gap between two fitted exponents. "
            "the sum stands")
        say("above it and has never been below it in this field.")
    elif not j2 and gap > 0:
        say("the j-side barrier is above l2. the demand asks "
            "|sum a| to fall below")
        say("what random signs give in its own index, and nothing "
            "this repository")
        say("cites delivers below square-root cancellation. that "
            "closes the sign")
        say("axis against item 5 rather than leaving it open.")
    elif not j2:
        say("l2 sits above the j-side barrier, so the demand is "
            "weaker than")
        say("square-root cancellation and the barrier is not what "
            "obstructs it.")
    else:
        say("the gap is inside the errors. the two barriers are not "
            "shown to")
        say("differ and are not shown to agree, and neither reading "
            "is drawn --")
        say("this field does not decide which floor item 5 is "
            "measured against.")

    io.open(OUT, "w", encoding="utf-8", newline="\n").write(
        "\n".join(HEAD + lines) + "\n")
    print("\nwrote %s" % os.path.normpath(OUT))
    return 0


if __name__ == "__main__":
    sys.exit(main())
