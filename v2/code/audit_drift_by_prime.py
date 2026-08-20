# -*- coding: utf-8 -*-
r"""
How the corrections' drift depends on the prime

WHAT IS AT STAKE

rem:localcorrections measured that the weight corrections move with N
and that they move most at the smallest prime: c(3) drifts at
-0.011551 +- 0.000876 with t = -13.19, c(7) at -0.005152 +- 0.000543
with t = -9.49, and c(11) at -0.001589 +- 0.001154 with t = -1.38,
unresolved.  Three primes is not a shape, and that remark said so.

This run measures the drift for every prime the field supports at all
eight N, and asks what it is a function of.

**And it names what it cannot ask.**  The natural candidates -- 1/p,
1/(p-1), 1/(p+1) -- agree to within a fifth for every p >= 11 and
differ only at 3 and 7.  Two points cannot choose among three
functions, so **no candidate is selected here**; what can be asked is
the exponent of the decay, and that is what GG4 asks.  G69 exists for
exactly this and the correlation is printed.

BACKS: Remark {#rem:driftbyprime} in paper/wall_v3.md.

PRE-REGISTERED PREDICTIONS (written before this script was run)

  GG1 THE GATE.  A at N = 200000 reproduces rem:maintermremoval's
      POINT mainA marker to a relative 1e-12, and c(3)'s drift
      reproduces rem:localcorrections' -0.011551 to six decimals.
  GG2 **The drift is real at the small primes.**  At least four
      primes have a drift resolved at |t| above 3.
  GG3 **And it dies away with p.**  Over the resolved primes, |s(p)|
      is largest at the smallest p and smallest at the largest.
  GG4 **The exponent.**  Fitting log|s(p)| on log p over the resolved
      primes, the slope is within 0.3 of -1.

REFUTATION RULE (fixed before the run)

  GG1 REFUTED outside either tolerance; nothing below is reported.
  GG2 REFUTED below four resolved primes.  Then the drift is a fact
      about c(3) and c(7) and not about the corrections as a family,
      and GG3 and GG4 have nothing to run on -- their verdicts would
      stand without a reading.
  GG3 REFUTED if the ordering fails at either end.  Then the drift is
      not a decreasing function of p and the picture
      rem:localcorrections drew from three primes was the three
      primes and not the pattern.
  GG4 **REFUTED outside -1 +- 0.3.**  Then the drift is not of order
      1/p, and since 1/p is what a density in a progression to
      modulus p would give, the correction's N-dependence does not
      come from that.  **What it comes from would then be unmeasured
      and must be left so** -- this branch has fitted enough shapes to
      quantities it could not derive.

  **A CRASH FIXED BEFORE ANY VERDICT EXISTED, DISCLOSED.**  The first
  execution died before GG1 with a TypeError: the noise threshold was
  stored in the same dict as the (c, w) pairs, and a set comprehension
  unpacks its target before it applies its filter, so the guard
  written for that key never ran.  The threshold now lives in its own
  dict; no quantity and no rule is changed and no verdict had been
  produced.

  **THE UNRESOLVED CASE, NAMED, WITH A NUMBER.**  1/p, 1/(p-1) and
  1/(p+1) are within a fifth of each other for p >= 11 here, so
  **this design cannot choose among them and will not be read as
  choosing**; GG4 is a statement about the exponent alone.  The
  correlation of log p with log(p-1) over the resolved set is printed
  and, if it exceeds 0.99, COEFF NOT SEPARABLE is emitted -- the
  exponent stays readable but no functional form does.  And if fewer
  than four primes resolve, GG4's verdict stands without a reading,
  as GG2's rule says.

  WHAT THIS CANNOT DO.  One radical family, eight N, and drifts
  fitted on eight points each.  A decay exponent measured over the
  primes this field reaches says nothing about larger p, and nothing
  here is derived -- the whole of it is description, and
  rem:residuemodel's closure of the reduction stands.
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
OUT = os.path.join(ROOT, "results", "audit_drift_by_prime.txt")
SRCM = os.path.join(ROOT, "results", "audit_mainterm_removal.txt")
SRCL = os.path.join(ROOT, "results", "audit_local_corrections.txt")

THETA = 0.56
NS = [25_000, 50_000, 100_000, 200_000, 400_000, 800_000,
      1_600_000, 3_200_000]
NGATE = 200_000
RELID = 1e-12
DEC = 6
TCAP = 3.0
MINRES = 4
EXPTARGET = -1.0
EXPTOL = 0.3
CORRCAP = 0.99
NOISEFRAC = 0.001


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


def phi(n):
    r = n
    for p in factor_set(n):
        r = r // p * (p - 1)
    return r


def pieces(N, lam, mu, sqf):
    PN = factor_set(N)
    K = int(N ** THETA)
    D = (N - 1) // K
    out = {}
    for d in range(1, D + 1):
        md = int(mu[d])
        if md == 0 or any(d % q == 0 for q in PN):
            continue
        ms = np.arange(K, (N - 1) // d + 1, dtype=np.int64)
        if ms.size == 0:
            continue
        keep = sqf[ms]
        for q in factor_set(d) | PN:
            keep &= (ms % int(q)) != 0
        if d == 1:
            keep &= lam[ms] == 0.0
        ms = ms[keep]
        if ms.size == 0:
            continue
        out[d] = (md, float((lam[N - d * ms]
                             * np.log(ms.astype(np.float64))).sum()))
        del ms, keep
    return out, D


def model(d, D):
    return (1.0 / phi(d)) * (1.0 - d / D) / (1.0 - 1.0 / D)


def fit(x, y):
    x = np.asarray(x, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)
    b, a0 = np.polyfit(x, y, 1)
    r = y - (b * x + a0)
    se = math.sqrt(float((r ** 2).sum() / (len(x) - 2))
                   / float(((x - x.mean()) ** 2).sum()))
    return float(b), se


def read_pub():
    m = re.search(r"^POINT mainA_%d ([-+]?[\d.eE+-]+)\s*$" % NGATE,
                  io.open(SRCM, encoding="utf-8").read(), re.M)
    c3 = re.search(r"^POINT cslope_3 ([-+]?[\d.]+)\s*$",
                   io.open(SRCL, encoding="utf-8").read(), re.M)
    if not m or not c3:
        raise SystemExit("a published value is missing")
    return float(m.group(1)), float(c3.group(1))


HEAD = [
    "STATISTIC: the drift of each prime's weight correction c(p, N)",
    "           against log N, one slope per prime, and how those",
    "           slopes depend on p.",
    "FIELD: N = %s; p over the primes coprime to N that contribute"
    % NS,
    "       at every one of those N with |w(p)| above the noise",
    "       threshold, inside D = floor((N-1)/K), K = floor(N^%.2f)."
    % THETA,
    "       A at N = %d and c(3)'s drift are READ from" % NGATE,
    "       results/audit_mainterm_removal.txt and",
    "       results/audit_local_corrections.txt.",
    "NOTE: 1/p, 1/(p-1) and 1/(p+1) are within a fifth of each other",
    "      for p >= 11 here, so no functional form is chosen; only the",
    "      exponent of the decay is read.",
    "",
]


def main():
    lines = []

    def say(t=""):
        print(t)
        sys.stdout.flush()
        lines.append(t)

    pubA, pubc3 = read_pub()
    say("READ audit_mainterm_removal.txt %d %.17e" % (NGATE, pubA))
    say("READ audit_local_corrections.txt cslope_3 %.6f" % pubc3)
    say("PRINTBOUND audit_drift_by_prime %d %.10f"
        % (DEC, 0.5 * 10.0 ** (-DEC)))
    say("  theta %.2f, |t| cap %.1f, minimum resolved %d, exponent "
        "%.1f +- %.1f," % (THETA, TCAP, MINRES, EXPTARGET, EXPTOL))
    say("  CORR cap %.2f, noise %.4f of |w(3)|" % (CORRCAP, NOISEFRAC))

    NMAX = max(NS)
    say("sieving to %d" % NMAX)
    lam, mu = lambda_and_mu(NMAX)
    sqf = mu != 0

    cs_by_N, noise_by_N, A200 = {}, {}, None
    for N in NS:
        ps, D = pieces(N, lam, mu, sqf)
        A = ps[1][1]
        if N == NGATE:
            A200 = A
        ws = {d: v / A for d, (md, v) in ps.items() if d != 1}
        noise = NOISEFRAC * abs(ws[3])
        cs_by_N[N] = {d: (ws[d] / model(d, D), abs(ws[d]))
                      for d in ws if d != D}
        noise_by_N[N] = noise
        say("  N = %-9d D = %-5d %d contributing d" % (N, D, len(ws)))

    common = None
    for N in NS:
        here = {p for p, cw in cs_by_N[N].items()
                if len(factor_set(p)) == 1
                and cw[1] >= noise_by_N[N]}
        common = here if common is None else (common & here)
    ps_common = sorted(common)
    say("SCALES %d" % len(NS))
    say("  %d primes contribute above the threshold at all eight N: %s"
        % (len(ps_common), ps_common))

    # ------------------------------------------------------------- GG1
    x = np.array([math.log(N) for N in NS])
    s3, e3 = fit(x, [cs_by_N[N][3][0] for N in NS])
    say()
    say("GG1  the gate")
    ra = abs(A200 - pubA) / max(abs(pubA), 1.0)
    okc = abs(round(s3, DEC) - round(pubc3, DEC)) < 10.0 ** (-DEC)
    gg1 = ra <= RELID and okc
    say("  A relative %.2e against %.0e   %s"
        % (ra, RELID, "ok" if ra <= RELID else "MISMATCH"))
    say("  c(3) drift here %+.6f against its %+.6f  %s"
        % (s3, pubc3, "ok" if okc else "MISMATCH"))
    say("  GG1 %s" % ("hold" if gg1 else "REFUTED"))
    if not gg1:
        io.open(OUT, "w", encoding="utf-8", newline="\n").write(
            "\n".join(HEAD + lines) + "\n")
        raise SystemExit(1)

    say()
    say("      p     drift        s.e.        t")
    sl = {}
    for p in ps_common:
        b, se = fit(x, [cs_by_N[N][p][0] for N in NS])
        sl[p] = (b, se)
        say("  %5d  %+.6f   %.6f   %+7.2f" % (p, b, se, b / se))
        say("POINT drift_%d %.6f" % (p, b))
        say("SPREAD drift_%d %.6f" % (p, se))

    # ------------------------------------------------------------- GG2
    say()
    say("GG2  how many primes have a resolved drift?")
    res = [p for p in ps_common if abs(sl[p][0] / sl[p][1]) > TCAP]
    gg2 = len(res) >= MINRES
    say("  resolved at |t| above %.1f: %s" % (TCAP, res))
    say("COUNT resolved %d" % len(res))
    say("  GG2 %s   (cap: at least %d)"
        % ("hold" if gg2 else "REFUTED", MINRES))

    # ------------------------------------------------------------- GG3
    say()
    say("GG3  does the drift die away with p?")
    if res:
        mags = [abs(sl[p][0]) for p in res]
        gg3 = (mags[0] == max(mags)) and (mags[-1] == min(mags))
        say("  |s(p)| over the resolved: %s"
            % " ".join("%.6f" % m for m in mags))
    else:
        gg3 = False
        say("  no resolved primes")
    say("  GG3 %s   (cap: largest at the smallest p, smallest at the "
        "largest)" % ("hold" if gg3 else "REFUTED"))

    # ------------------------------------------------------------- GG4
    say()
    say("GG4  the exponent of the decay")
    if len(res) >= 2:
        lx = np.array([math.log(p) for p in res])
        ly = np.array([math.log(abs(sl[p][0])) for p in res])
        b4, se4 = fit(lx, ly)
        gg4 = abs(b4 - EXPTARGET) <= EXPTOL
        if math.isfinite(se4):
            say("  slope of log|s| on log p: %+.6f +- %.6f"
                % (b4, se4))
            say("TSTAT driftexp %.2f" % ((b4 - EXPTARGET) / se4))
            say("SPREAD driftexp %.6f" % se4)
        else:
            say("  slope of log|s| on log p: %+.6f, standard error "
                "not finite" % b4)
            say("  %d resolved primes give %d degrees of freedom, so "
                "the line passes"
                % (len(res), max(len(res) - 2, 0)))
            say("  through its points exactly and no t is emitted -- "
                "a t against an")
            say("  infinite error is not a t")
        say("POINT driftexp %.6f" % b4)
        c = float(np.corrcoef(lx, np.log(np.array(res) - 1.0))[0, 1])
        say("CORR driftbyprime_regressors %.6f" % abs(c))
        if abs(c) >= CORRCAP:
            say("COEFF NOT SEPARABLE driftbyprime")
            say("  log p and log(p-1) correlate at %.6f over the "
                "resolved set, so" % c)
            say("  no functional form is chosen and only the exponent "
                "is read")
    else:
        gg4 = False
        b4 = float("nan")
        say("  fewer than two resolved primes; no fit is made")
    say("  GG4 %s   (cap: %.1f +- %.1f)"
        % ("hold" if gg4 else "REFUTED", EXPTARGET, EXPTOL))
    if len(res) < MINRES:
        say("  UNREADABLE: fewer than %d resolved primes, so GG4's "
            "verdict stands" % MINRES)
        say("  without a reading, as GG2's rule says")

    say()
    say("=" * 70)
    say("GG1 %s  GG2 %s  GG3 %s  GG4 %s"
        % tuple("hold" if v else "REFUTED"
                for v in (gg1, gg2, gg3, gg4)))
    say()
    if gg2 and gg3 and gg4:
        say("the corrections' drift dies away with p at an exponent "
            "this field")
        say("puts near minus one, which is what a density in a "
            "progression to")
        say("modulus p would give. no functional form is chosen -- "
            "the candidates")
        say("are not separable here -- and nothing is derived. it is "
            "a description")
        say("of how the description moves.")
    elif not gg2:
        say("too few primes resolve. the drift is a fact about the "
            "smallest primes")
        say("and not about the corrections as a family, and the shape "
            "questions")
        say("have nothing to run on.")
    elif not gg3:
        say("the drift is not a decreasing function of p, so what "
            "rem:localcorrections")
        say("drew from three primes was the three primes and not the "
            "pattern.")
    else:
        say("the drift dies away but not at the exponent a modulus-p "
            "density would")
        say("give. where it comes from is unmeasured and is left so.")

    io.open(OUT, "w", encoding="utf-8", newline="\n").write(
        "\n".join(HEAD + lines) + "\n")
    print("\nwrote %s" % os.path.normpath(OUT))
    return 0


if __name__ == "__main__":
    sys.exit(main())
