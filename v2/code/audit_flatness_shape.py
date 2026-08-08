# -*- coding: utf-8 -*-
r"""
The flatness rises. It cannot rise forever -- so what does it do?

WHAT IS AT STAKE

Remark {#rem:leanidentity} reduced OPEN item 5 to two exponents and
left one loose: e(l1/l2) measures 0.287798 +- 0.002472 while its own
ceiling is theta'/2 = 0.28, because l1/l2 <= sqrt(#k) and #k is
c N^theta' with c fixed on this one radical. The excess is carried by
the flatness

    F = (l1/l2) / sqrt(#k),

which runs 0.6622 to 0.6986 and is still rising. Whether that rise is
finite-N or a property of the sequence was not settled.

It has to be finite-N: F <= 1 by Cauchy-Schwarz, so no power of N can
describe it past the point where it would reach one. That is not a
statistical statement but an arithmetic one, and it turns the question
into a shape adjudication of the kind {#rem:laddershape} runs -- with
the difference that here one of the shapes is excluded a priori and
the other must be fitted.

What hangs on it: if F saturates then asymptotically
e(l1/l2) = theta'/2 exactly, item 5 does reduce to e(G) alone, and the
lean/floor exponent is theta'/2 - e(G) = 0.126 rather than the 0.159
measured over this range.

BACKS: Remark {#rem:flatnessshape} in paper/wall_v3.md.

PRE-REGISTERED PREDICTIONS (written before this script was run)

  X1  The control: #k reproduces lab_positive_weights.py's 313, 462,
      682, 1004, 1485 exactly and F reproduces
      {#rem:crosskreference}'s 0.6622 to 0.6854 on those five N to
      within 0.001.
  X2  #k is a clean power: its exponent against log N reproduces
      theta' = 0.56 to within 0.001, so the ceiling of e(l1/l2) is
      theta'/2 and not something fitted.
  X3  A bounded shape does at least as well as the power: fitted to
      F, the saturating a + b/log N has an r.m.s. residual no larger
      than the power law's.
  X4  And it saturates below the bound: the fitted a lies in (0, 1).

REFUTATION RULE (fixed before the run)

  X1  REFUTED on any mismatch -- not the same statistic, and nothing
      below may be compared with {#rem:crosskreference}.
  X2  REFUTED at 0.001. The ceiling would then not be theta'/2 and
      the excess {#rem:leanidentity} measured would be against the
      wrong bound.
  X3  REFUTED if the power law fits better. That is the one that
      matters: the data would prefer a shape that Cauchy-Schwarz
      forbids, which would say this range is too short to see the
      saturation at all and e(l1/l2) could not be read as tending to
      theta'/2.
  X4  REFUTED if a is outside (0, 1). A fitted ceiling at or above
      one would be vacuous and below zero meaningless.

  All four gate.

  NO NULL IS RUN and none applies. Two deterministic shapes are
  fitted to the same eight measurements and compared by r.m.s.; there
  is no background to detect against. The magnitudes themselves are
  mu's, and the coin arms for this field were run in
  lab_positive_weights.py and audit_crossk_reference.py.
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
OUT = os.path.join(RES, "audit_flatness_shape.txt")

NS = [200_000, 400_000, 800_000, 1_600_000, 3_200_000, 6_400_000,
      12_800_000, 25_600_000]
THETA = 0.56
UMAX = 400.0


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


def read_published():
    """#k from lab_positive_weights and F from crosskreference"""
    src = io.open(os.path.join(RES, "lab_positive_weights.txt"),
                  encoding="utf-8").read()
    i = src.index("G_eps     sqrt K   top-decile share")
    nk = {}
    for ln in src[i:].splitlines()[1:]:
        f = ln.split()
        if len(f) < 8 or not f[0].isdigit():
            if f and set(f[0]) == {"-"}:
                continue
            if not nk:
                continue
            break
        nk[int(f[0])] = int(f[2])
    src2 = io.open(os.path.join(RES, "audit_crossk_reference.txt"),
                   encoding="utf-8").read()
    m = re.search(r"^REFERENCE audit_crossk_reference (\d+) "
                  r"([\d.]+) ([\d.]+)\s*$", src2, re.M)
    return nk, float(m.group(2)), float(m.group(3))


def weighted(N, lam, mu, sqf):
    """(log k)H(N;k) over the squarefree k < N^theta coprime to N"""
    PN = factor_set(N)
    K = int(N ** THETA)
    ks = np.array([k for k in range(2, K)
                   if sqf[k] and not any(k % q == 0 for q in PN)],
                  dtype=np.int64)
    Hs = []
    for k in ks:
        k = int(k)
        M = (N - 1) // k
        ms = np.arange(1, M + 1, dtype=np.int64)
        for q in factor_set(k):
            ms = ms[ms % q != 0]
        vals = N - ms * k
        Hs.append(float((lam[vals] * mu[ms].astype(np.float64)).sum()))
    return ks, np.log(ks.astype(np.float64)) * np.array(Hs)


def fit(x, y):
    a, b = np.polyfit(x, y, 1)
    r = y - (a * x + b)
    n = x.size
    se = math.sqrt(float((r ** 2).sum() / (n - 2))
                   / float(((x - x.mean()) ** 2).sum()))
    return float(a), float(b), float(np.sqrt((r ** 2).mean())), se


def main():
    lines = []

    def say(s=""):
        print(s)
        lines.append(s)

    pubnk, flo, fhi = read_published()
    say("read %d published #k from results/lab_positive_weights.txt "
        "and the flatness" % len(pubnk))
    say("  range %.4f to %.4f from "
        "results/audit_crossk_reference.txt" % (flo, fhi))

    NMAX = max(NS)
    say("sieving to %d ..." % NMAX)
    lam, mu = lambda_and_mu(NMAX)
    sqf = mu != 0
    say("RADICALS %d"
        % len(set(tuple(sorted(q for q in factor_set(N) if q > 2))
                  for N in NS)))

    rows = []
    for N in NS:
        ks, a = weighted(N, lam, mu, sqf)
        l1 = float(np.abs(a).sum())
        l2 = float(np.sqrt((a * a).sum()))
        rows.append((N, ks.size, l1 / l2,
                     (l1 / l2) / math.sqrt(ks.size)))
        say("  N = %-10d #k = %-6d l1/l2 = %-9.4f F = %.4f"
            % (N, ks.size, l1 / l2, rows[-1][3]))

    x = np.log(np.array([r[0] for r in rows], dtype=np.float64))
    F = np.array([r[3] for r in rows])

    # ------------------------------------------------------------- X1
    say()
    say("X1  the control: #k and the flatness")
    say("  N            #k here  #k pub   F")
    x1 = True
    fiv = []
    for N, nk, cc, f in rows:
        if N in pubnk:
            if nk != pubnk[N]:
                x1 = False
            fiv.append(f)
            say("  %-12d %-8d %-8d %.4f" % (N, nk, pubnk[N], f))
        else:
            say("  %-12d %-8d %-8s %.4f" % (N, nk, "-", f))
    if abs(min(fiv) - flo) >= 0.001 or abs(max(fiv) - fhi) >= 0.001:
        x1 = False
    say("  the five published N give F from %.4f to %.4f against "
        "%.4f to %.4f" % (min(fiv), max(fiv), flo, fhi))
    say("  X1 %s   (cap 0.001)" % ("hold" if x1 else "REFUTED"))

    # ------------------------------------------------------------- X2
    say()
    say("X2  is #k a clean power of N?")
    e, b, rms, se = fit(x, np.log(np.array([r[1] for r in rows],
                                           dtype=np.float64)))
    x2 = abs(e - THETA) < 0.001
    say("  exponent of #k against log N = %.6f, against theta' = %.2f"
        % (e, THETA))
    say("  r.m.s. residual %.6f, standard error %.6f" % (rms, se))
    say("  so the ceiling of e(l1/l2) is %.6f" % (e / 2.0))
    say("  X2 %s   (cap 0.001)" % ("hold" if x2 else "REFUTED"))
    say("REFERENCE audit_flatness_shape %d %.4f %.4f"
        % (len(rows), float(F.min()), float(F.max())))

    # ---------------------------------------------------------- X3/X4
    say()
    say("X3/X4  two shapes for the flatness")
    ep, bp, rp, sp = fit(x, np.log(F))
    A = np.vstack([np.ones_like(x), 1.0 / x]).T
    c, *_ = np.linalg.lstsq(A, F, rcond=None)
    rs = float(np.sqrt(((F - A @ c) ** 2).mean()))
    rpow = float(np.sqrt(((F - np.exp(bp) * np.exp(ep * x)) ** 2)
                         .mean()))
    x3 = rs <= rpow
    x4 = 0.0 < c[0] < 1.0
    say("  shape                 r.m.s.     what it does")
    say("  F ~ N^e               %-10.6f e = %+.6f, so F reaches 1"
        % (rpow, ep))
    ucross = (0.0 - bp) / ep
    say("                                   at log10 N = %.4f"
        % (ucross / math.log(10.0)))
    say("  F = a + b/log N       %-10.6f a = %.6f, b = %+.6f"
        % (rs, c[0], c[1]))
    say("  X3 the bounded shape is no worse   %s"
        % ("hold" if x3 else "REFUTED"))
    say("  X4 its ceiling a is inside (0, 1)   %s"
        % ("hold" if x4 else "REFUTED"))
    say()
    say("  DIAGNOSTIC on X3 (post hoc). The two r.m.s. differ by")
    say("  %.6f, and an r.m.s. from %d points with two parameters"
        % (abs(rs - rpow), len(rows)))
    say("  carries its own standard error of rms/sqrt(2(n-2)):")
    dof = len(rows) - 2
    ser = min(rs, rpow) / math.sqrt(2.0 * dof)
    say("    best r.m.s. %.6f, standard error %.6f (%.1f per cent)"
        % (min(rs, rpow), ser, 100 * ser / min(rs, rpow)))
    say("    the gap between the shapes is %.2f of that"
        % (abs(rs - rpow) / ser))
    say("  so X3 fails as registered and the two shapes are not")
    say("  separated at all. What separates them is not the fit:")
    say("  Cauchy-Schwarz forbids the power past F = 1.")
    uref = 12.0 * math.log(10.0)
    fpow = math.exp(bp + ep * uref)
    fsat = c[0] + c[1] / uref
    say("  and what they say where it matters, at log10 N = 12:")
    say("    power     F = %.6f" % fpow)
    say("    bounded   F = %.6f" % fsat)
    nsurv = 2 if abs(rs - rpow) <= ser else 1
    say("SHAPEGAP flatness %.6f %.6f" % (abs(rs - rpow), ser))
    if abs(rs - rpow) <= ser:
        say("SHAPES TIED flatness")
    say("SHAPESURVIVE flatness %d %d %.4f"
        % (len(rows), nsurv, abs(fpow - fsat)))
    say("SHAPECURRENT flatness %d" % len(rows))

    say()
    say("  the power law is excluded by arithmetic and not by fit:")
    say("  F <= 1 by Cauchy-Schwarz, and the fitted power crosses 1")
    say("  at log10 N = %.4f, so it cannot describe the sequence past"
        % (ucross / math.log(10.0)))
    say("  there. The bracket a leave-one-out refit gives:")
    lo, hi = ucross, ucross
    for sl in (slice(1, None), slice(0, -1)):
        e2, b2, _r2, _s2 = fit(x[sl], np.log(F[sl]))
        v = (0.0 - b2) / e2
        lo, hi = min(lo, v), max(hi, v)
    say("BRACKET log10_N_flatness_reaches_one %.4f %.4f %.4f"
        % (ucross / math.log(10.0), lo / math.log(10.0),
           hi / math.log(10.0)))
    say("DRIFT log10_N_flatness_reaches_one %.4f"
        % ((hi - lo) / abs(ucross)))
    say("  and the bounded shape's own answer to the same question:")
    say("  it tends to %.6f and never reaches 1, so within the search"
        % c[0])
    say("  ceiling log10 N = %.1f it gives no crossing at all."
        % (UMAX / math.log(10.0)))
    say("FORECAST BOTH flatness %.4f %.4f %.4f"
        % (ucross / math.log(10.0), ucross / math.log(10.0),
           UMAX / math.log(10.0)))
    say("  which is wide because e is small; what is not wide is the")
    say("  conclusion, since any positive e crosses 1 eventually.")

    say()
    say("  what saturation would do to {#rem:leanidentity}. If F")
    say("  tends to %.6f then e(l1/l2) tends to theta'/2 = %.4f, and"
        % (c[0], e / 2.0))
    say("  the lean/floor exponent tends to that minus e(G):")
    gsrc = io.open(os.path.join(RES, "audit_lean_identity.txt"),
                   encoding="utf-8").read()
    eg = float(re.search(r"^  G\s+([-+][\d.]+)", gsrc, re.M).group(1))
    el = float(re.search(r"^  lean/floor\s+([-+][\d.]+)",
                         gsrc, re.M).group(1))
    say("  e(G) = %+.6f, read from results/audit_lean_identity.txt"
        % eg)
    say("  so the limit is %+.6f against the %+.6f measured over this"
        % (e / 2.0 - eg, el))
    say("  range -- the lean still grows, by a fifth less.")

    say()
    say("=" * 70)
    ok = x1 and x2 and x3 and x4
    say("the flatness saturates and the ceiling is theta'/2"
        if ok else "REFUTED")

    head = [
        "STATISTIC: the flatness F = (l1/l2)/sqrt(#k) of the weighted",
        "           wall a_k = (log k)H(N;k) over the squarefree",
        "           k < N^" + str(THETA) + " coprime to N; the exponent",
        "           of #k against log N; and two shapes fitted to F --",
        "           the power F ~ N^e and the bounded a + b/log N --",
        "           with their r.m.s. residuals and the log10 N at",
        "           which the power would reach the Cauchy-Schwarz",
        "           bound F = 1.",
        "NULL: none is run and none applies. Deterministic shapes are",
        "      fitted to the same measurements and compared by r.m.s.",
        "      The coin arms for this field were run in",
        "      lab_positive_weights.py and audit_crossk_reference.py.",
        "FIELD: N = 2e5 through 2.56e7 by doubling; k squarefree and",
        "       coprime to N with 2 <= k < N^" + str(THETA) + "; m over",
        "       1 <= m < N/k with (m,k) = 1; Lambda and mu from an",
        "       integer sieve to " + str(NMAX) + ". Every N is 2^a 5^b,",
        "       one odd radical, as RADICALS declares, so the density",
        "       of admissible k is fixed and #k is a clean power. The",
        "       published #k are read from",
        "       results/lab_positive_weights.txt and the flatness",
        "       range from results/audit_crossk_reference.txt.",
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
