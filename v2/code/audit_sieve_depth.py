# -*- coding: utf-8 -*-
r"""
Is the missing piece sieve depth? The sweep ends at primality.

WHAT IS AT STAKE

OPEN item 4: nothing elementary carries the sign lean on the k the
lean is measured over. Remark {#rem:logweightpredictor} ruled out the
Lambda weight as the missing piece -- putting log(N-mk) on the
sieve-weighted P lowers the sign agreement at every N. What has never
been varied is the other half of P's definition: the sieve runs over
the odd primes below 30, and that number was inherited, not chosen.

The depth matters more on long inner sums than on short ones, which is
exactly where P fails. A term survives P when N - mk has no odd prime
factor below Q; it contributes to H when N - mk is prime. Those two
sets converge as Q grows, and they COINCIDE at Q = ceil(sqrt(N)),
since a value below N with no prime factor at or below its square root
is prime.

So the sweep is not a search over models. It is a path from the
published predictor to H's own support, and its far end isolates one
thing: at the top rung the only difference left between the predictor
and H is the weight log p. Remark {#rem:logweightpredictor} says the
weight does not help the agreement; if the top rung also fails to
carry the lean, no unweighted elementary predictor can, and if it
succeeds the missing piece was depth all along.

One simplification makes the sweep cheap. The published weight skips
the primes q dividing k. That skip is vacuous here: k is coprime to N,
so q | k gives N - mk = N mod q with q not dividing N, and such a q
never divides a term. The survivor set is therefore "no odd prime
factor below Q" outright, one boolean array per Q.

BACKS: Remark {#rem:sievedepth} in paper/wall_v3.md.

PRE-REGISTERED PREDICTIONS (written before this script was run)

  R1  The control: at Q = 29 the sign agreement reproduces
      {#rem:survivorrange}'s 0.8129, 0.7632, 0.7446, 0.7367, 0.7759,
      0.7579 and the lean ratios its 0.6571 to 0.8676, both to within
      0.001.
  R2  Depth helps: the agreement rises with Q at every N.
  R3  And the far end carries the lean: at Q = ceil(sqrt(N)), where
      the survivors are exactly the primes, the predicted |0.5 - f| is
      within a factor 1.05 of mu's at every N.
  R4  And its trend: at that Q the decay slope is within two standard
      errors of mu's.

REFUTATION RULE (fixed before the run)

  R1  REFUTED at 0.001 anywhere -- not the same statistic, and nothing
      below may be compared with {#rem:survivorrange}.
  R2  REFUTED if the agreement fails to rise at any N. Depth would
      then not be the axis, and the sweep would say so before its end
      point does.
  R3  REFUTED outside a factor 1.05 at any N. That is the one that
      matters: the survivors ARE the primes there, so the only thing
      left is the weight log p, and {#rem:logweightpredictor} has
      already ruled the weight out as a help -- OPEN item 4 would then
      be closed in the negative, with no unweighted elementary
      predictor able to carry the lean.
  R4  REFUTED beyond two standard errors, the same conclusion for the
      trend rather than the level.

  All four gate.

  THE NULL is the 256 global sign vectors of audit_lean_floor.py on
  the identical magnitudes, so every rung's lean is read on the same
  scale as mu's. The permutation control for this predictor family was
  run in lab_survivor_selection.py at 0.5372 to 0.5414.
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
OUT = os.path.join(RES, "audit_sieve_depth.txt")

NS = [200_000, 400_000, 800_000, 1_600_000, 3_200_000, 6_400_000]
THETA = 0.56
QS = [29, 53, 101, 211, 503, 1009]
DRAWS = 256
SEED = 20260808


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


def nosmall(n, q):
    """not divisible by any odd prime at or below q -- the sieve's
    own convention, which strikes a small prime along with its
    multiples because a sieve does not know it is prime"""
    s = np.ones(n + 1, dtype=bool)
    s[0] = False
    for p in primes_upto(q):
        p = int(p)
        if p == 2:
            continue
        s[p::p] = False
    return s


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
    """the agreements and lean ratios of {#rem:survivorrange}"""
    src = io.open(os.path.join(RES, "audit_survivor_range.txt"),
                  encoding="utf-8").read()
    i = src.index("N            #k     agreement  mu lean   P lean"
                  "    ratio")
    agr, rat = {}, {}
    for ln in src[i:].splitlines()[1:]:
        f = ln.split()
        if len(f) < 6 or not f[0].isdigit():
            break
        agr[int(f[0])] = float(f[2])
        rat[int(f[0])] = float(f[5])
    return agr, rat


def fit(x, y):
    a, b = np.polyfit(x, y, 1)
    r = y - (a * x + b)
    n = x.size
    se = math.sqrt(float((r ** 2).sum() / (n - 2))
                   / float(((x - x.mean()) ** 2).sum()))
    return float(a), float(np.sqrt((r ** 2).mean())), se, abs(a) / se


def main():
    lines = []

    def say(s=""):
        print(s)
        lines.append(s)

    MAJS = []          # the predictor's own majority share

    pagr, prat = read_published()
    say("read %d agreements and lean ratios from "
        "results/audit_survivor_range.txt" % len(pagr))

    NMAX = max(NS)
    say("sieving to %d ..." % NMAX)
    lam, mu = lambda_and_mu(NMAX)
    sqf = mu != 0
    oddsqf = sqf.copy()
    oddsqf[::2] = False
    rng = np.random.default_rng(SEED)
    say("RADICALS %d"
        % len(set(tuple(sorted(q for q in factor_set(N) if q > 2))
                  for N in NS)))

    qtop = {N: int(math.isqrt(N)) + 1 for N in NS}
    say("the sweep runs over Q = %s and, at each N, the top rung "
        "Q = ceil(sqrt(N))" % ", ".join(map(str, QS)))
    say("  where the survivors are exactly the primes: %s"
        % ", ".join("%d->%d" % (N, qtop[N]) for N in NS))
    surv = {}
    for q in sorted(set(QS) | set(qtop.values())):
        surv[q] = nosmall(NMAX, q)

    agr = {}          # (N, q) -> sign agreement
    lean = {}         # (N, q) -> |0.5 - f|
    muln, flo = {}, {}
    for N in NS:
        PN = factor_set(N)
        qs_here = QS + [qtop[N]]
        ks, Hs = [], []
        Ps = dict((q, []) for q in qs_here)
        for k in range(2, int(N ** THETA)):
            if not sqf[k]:
                continue
            if any(k % q == 0 for q in PN):
                continue
            M = N // k
            if M < 1:
                continue
            ms = np.arange(1, M + 1, 2, dtype=np.int64)
            ms = ms[oddsqf[ms]]
            for q in factor_set(k):
                if q > 2:
                    ms = ms[ms % q != 0]
            if ms.size == 0:
                continue
            vals = N - ms * k
            g = mu[ms].astype(np.float64)
            ks.append(k)
            Hs.append(float((lam[vals] * g).sum()))
            for q in qs_here:
                Ps[q].append(float(g[surv[q][vals]].sum()))
        ks = np.array(ks, dtype=np.int64)
        H = np.array(Hs)
        a = np.log(ks.astype(np.float64)) * H
        l1 = float(np.abs(a).sum())
        w = np.abs(a)
        sh = np.sign(H)
        muln[N] = abs(0.5 - float(a[a > 0].sum() / l1))
        eps = (rng.integers(0, 2, size=(DRAWS, ks.size))
               .astype(np.int8) * 2 - 1)
        flo[N] = float(np.median(np.abs((eps @ w) / (2.0 * l1))))
        for q in qs_here:
            P = np.array(Ps[q])
            sp = np.sign(P)
            MAJS.append(max(float((sp > 0).mean()),
                            float((sp < 0).mean())))
            ok = (sh != 0) & (sp != 0)
            agr[(N, q)] = float((sh[ok] == sp[ok]).mean())
            spz = np.where(sp == 0, 1.0, sp)
            lean[(N, q)] = abs(0.5 - float(w[spz > 0].sum() / l1))
        say("  N = %-10d #k = %-6d done" % (N, ks.size))

    # ---------------------------------------------------------- R1/R2
    say()
    say("R1/R2  the sign agreement as the sieve deepens")
    say("  N            " + "".join("Q=%-8d" % q for q in QS)
        + "Q=sqrt(N)")
    r1 = r2 = True
    sat, fell = set(), []
    for N in NS:
        row = [agr[(N, q)] for q in QS] + [agr[(N, qtop[N])]]
        if abs(row[0] - pagr[N]) >= 0.001:
            r1 = False
        for i in range(len(row) - 1):
            if row[i + 1] < row[i]:
                r2 = False
                fell.append((N, i, row[i] - row[i + 1]))
            elif row[i + 1] == row[i]:
                r2 = False
                sat.add(N)
        say("  %-12d %s" % (N, "".join("%-10.4f" % v for v in row)))
    say("  R1 Q = %d reproduces the published agreement   %s"
        "   (cap 0.001)" % (QS[0], "hold" if r1 else "REFUTED"))
    say("  R2 the agreement rises with Q at every N   %s"
        % ("hold" if r2 else "REFUTED"))
    if sat or fell:
        say("  DIAGNOSTIC on R2 (post hoc). Every failure is past")
        say("  ceil(sqrt(N)), and both kinds are arithmetic rather")
        say("  than statistical.")
        if sat:
            say("  Ties at N = %s: beyond the square root no composite"
                % ", ".join(str(N) for N in sorted(sat)))
            say("  below N has a factor left to remove, so the columns")
            say("  saturate. The sweep's ladder overshoots sqrt(N)")
            say("  there -- %s against the rungs %s."
                % (", ".join("%d" % qtop[N] for N in sorted(sat)),
                   ", ".join(map(str, QS))))
        if fell:
            say("  Falls, %d of them, the largest %.4f:"
                % (len(fell), max(d for _n, _i, d in fell)))
            for N, i, d in fell:
                say("    N = %-10d Q = %d to %d, by %.4f"
                    % (N, QS[i], QS[i + 1] if i + 1 < len(QS)
                       else qtop[N], d))
            say("  A sieve strikes a small prime along with its")
            say("  multiples, so raising Q past sqrt(N) removes true")
            say("  contributors -- the primes at or below Q -- without")
            say("  removing any composite. That is why the top rung is")
            say("  ceil(sqrt(N)) and not the largest rung on the")
            say("  ladder, and it is where R3 and R4 are judged.")

    # ---------------------------------------------------------- R3/R4
    say()
    say("R3  the lean each depth gives, as a ratio to mu's")
    say("  N            " + "".join("Q=%-8d" % q for q in QS)
        + "Q=sqrt(N)")
    r3 = True
    for N in NS:
        row = ([lean[(N, q)] / muln[N] for q in QS]
               + [lean[(N, qtop[N])] / muln[N]])
        if abs(row[0] - prat[N]) >= 0.001:
            r1 = False
        if not (1.0 / 1.05 <= row[-1] <= 1.05):
            r3 = False
        say("  %-12d %s" % (N, "".join("%-10.4f" % v for v in row)))
    say("  R1 also covers the published lean ratios   %s"
        % ("hold" if r1 else "REFUTED"))
    say("  R3 the top rung is within a factor 1.05   %s"
        % ("hold" if r3 else "REFUTED"))
    say("PERN sievedepth_top_over_mu %d %.4f %.4f"
        % (len(NS),
           min(lean[(N, qtop[N])] / muln[N] for N in NS),
           max(lean[(N, qtop[N])] / muln[N] for N in NS)))

    say()
    say("R4  and the trends")
    x = np.log(np.array(NS, dtype=np.float64))
    bm, rm, sem, tm = fit(x, np.log(np.array([muln[N] for N in NS])))
    bt, rt, set_, tt = fit(
        x, np.log(np.array([lean[(N, qtop[N])] for N in NS])))
    sd = math.sqrt(sem ** 2 + set_ ** 2)
    r4 = abs(bm - bt) <= 2.0 * sd
    say("  mu            slope %+.6f, standard error %.6f"
        % (bm, sem))
    say("  Q = sqrt(N)   slope %+.6f, standard error %.6f"
        % (bt, set_))
    say("  difference    %+.6f against %.6f = 2 s.e., i.e. %.2f s.e."
        % (bm - bt, 2.0 * sd, abs(bm - bt) / sd))
    say("SCATTER slope_audit_sieve_depth %.4f" % rt)
    say("TSTAT slope_audit_sieve_depth %.2f" % tt)
    say("SPREAD slope_audit_sieve_depth %.4f"
        % float(x.max() - x.min()))
    if tt < 2.0:
        say("UNRESOLVED SIGN slope_audit_sieve_depth")
    say("  R4 %s" % ("hold" if r4 else "REFUTED"))

    say()
    say("  the ledger, on the criterion these predictors were built")
    say("  for. Modd and the log-weighted variant are read from")
    say("  results/audit_oddmertens_range.txt and")
    say("  results/audit_logweight_predictor.txt:")
    osrc = io.open(os.path.join(RES, "audit_oddmertens_range.txt"),
                   encoding="utf-8").read()
    oagr = [float(v) for v in re.findall(
        r"^  \d+\s+\d+\s+\d+\s+([\d.]+)\s+[\d.]+\s+[\d.]+"
        r"\s+[\d.]+\s*$", osrc, re.M)]
    orat = [float(v) for v in re.findall(
        r"^  \d+\s+\d+\s+\d+\s+[\d.]+\s+[\d.]+\s+[\d.]+"
        r"\s+([\d.]+)\s*$", osrc, re.M)]
    oslope = float(re.search(r"Modd      slope ([-+][\d.]+)",
                             osrc).group(1))
    lsrc = io.open(os.path.join(RES, "audit_logweight_predictor.txt"),
                   encoding="utf-8").read()
    lrow = re.search(r"^PREDICTOR ladder_lean sieve_P_log "
                     r"([\d.]+) ([\d.]+) ([-+][\d.]+)\s*$",
                     lsrc, re.M)
    rows = [("modd", min(oagr), max(orat), oslope),
            ("sieve_P", min(pagr.values()), max(prat.values()),
             float(re.search(r"^PREDICTOR ladder_lean sieve_P "
                             r"[\d.]+ [\d.]+ ([-+][\d.]+)\s*$",
                             lsrc, re.M).group(1))),
            ("sieve_P_log", float(lrow.group(1)),
             float(lrow.group(2)), float(lrow.group(3))),
            ("sieve_to_sqrtN",
             min(agr[(N, qtop[N])] for N in NS),
             max(lean[(N, qtop[N])] / muln[N] for N in NS), bt)]
    say("  predictor        worst agreement  worst lean ratio  slope")
    for nm, a0, r0, sl in rows:
        say("  %-16s %-16.4f %-17.4f %+.6f" % (nm, a0, r0, sl))
        say("PREDICTOR ladder_lean %s %.4f %.4f %+.6f"
            % (nm, a0, r0, sl))
    say("PREDICTOR CRITERION ladder_lean agreement")
    best = max(rows, key=lambda t: t[1])[0]
    say("PREDICTOR BEST ladder_lean %s" % best)
    say()
    say("  and the sieve level each one uses, which is what decides")
    say("  whether it is elementary in the sense {#rem:provablehalf}")
    say("  needs -- every condition to a BOUNDED modulus:")
    say("  predictor        level")
    for nm, lv in (("modd", "0"), ("sieve_P", str(QS[0])),
                   ("sieve_P_log", str(QS[0])),
                   ("sieve_to_sqrtN", "sqrt")):
        say("  %-16s %s" % (nm, lv))
        say("LEVEL %s %s" % (nm, lv))
    say("UNBOUNDED LEVEL sieve_to_sqrtN")
    say("  so the winner on the declared criterion is the one whose")
    say("  level grows with N. At fixed level 29 the agreement is")
    say("  %.4f to %.4f and the lean ratio %.4f to %.4f, and both"
        % (min(pagr.values()), max(pagr.values()),
           min(prat.values()), max(prat.values())))
    say("  get worse as N grows.")

    say()
    say("  and against the floor of {#rem:leanfloor}:")
    say("  N            mu/floor   Q=sqrt(N)/floor")
    for N in NS:
        say("  %-12d %-10.4f %.4f"
            % (N, muln[N] / flo[N], lean[(N, qtop[N])] / flo[N]))

    say()
    say("=" * 70)
    ok = r1 and r2 and r3 and r4
    say("the missing piece is sieve depth" if ok else "REFUTED")

    mj = max(MAJS)
    say()
    say("  the predictor's own majority sign share, at its worst over "
        "everything")
    say("  reported above: %.4f. An agreement is only a measurement "
        "where the" % mj)
    say("  predictor has variance; where it takes one sign almost "
        "everywhere,")
    say("  the agreement is the other side's marginal rate read back.")
    say("MARGINAL %s %.4f" % ("audit_sieve_depth", mj))
    if mj >= 0.9:
        say("DEGENERATE %s" % "audit_sieve_depth")

    head = [
        "STATISTIC: on the squarefree k < N^" + str(THETA)
        + " coprime to N, the",
        "           sign agreement of H(N;k) with",
        "           P_Q = sum_m mu(m) [N - mk has no odd prime factor",
        "           at or below Q], for Q = " + ", ".join(map(str, QS)),
        "           and for Q = ceil(sqrt(N)), where the survivors are",
        "           exactly the primes; the mass-weighted lean each",
        "           gives on mu's own (log k)|H| magnitudes; the top",
        "           rung's slope against log N with its standard",
        "           error; and all of it against the median lean of "
        + str(DRAWS),
        "           sign vectors on those magnitudes.",
        "NULL: " + str(DRAWS) + " global sign vectors on the identical",
        "      magnitudes, the convention of audit_lean_floor.py. The",
        "      permutation control for this predictor family was run",
        "      in lab_survivor_selection.py at 0.5372 to 0.5414.",
        "FIELD: N = 2e5 through 6.4e6 by doubling; k squarefree and",
        "       coprime to N with 2 <= k < N^" + str(THETA) + "; m odd,",
        "       squarefree and coprime to k, m <= N//k, the convention",
        "       of lab_survivor_selection.py; Lambda and mu from an",
        "       integer sieve to " + str(NMAX) + "; numpy default_rng",
        "       seed " + str(SEED) + ". The published weight's skip of",
        "       primes dividing k is vacuous here and is dropped: k is",
        "       coprime to N, so such a prime never divides N - mk.",
        "       Every N is 2^a 5^b, one odd radical, as RADICALS says.",
        "       The published agreements and ratios are read from",
        "       results/audit_survivor_range.txt.",
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
