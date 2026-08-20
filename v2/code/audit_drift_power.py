# -*- coding: utf-8 -*-
r"""
Buying power: do 11 and 13 resolve when the field is three N longer

WHAT IS AT STAKE

rem:driftbyprime measured the drift of every weight correction c(p, N)
over twenty primes and found two resolved: p = 3 at t = -13.19 and
p = 7 at t = -9.49, with nothing beyond p = 7 reaching |t| = 2 and the
signs mixed.  It refused to say whether that means 3 and 7 are special
or the other eighteen are too noisy, because **the errors at p >= 11
are larger than the effect at p = 3** -- so the absence is an absence
of power as much as an absence of drift.

Power is buyable here and the reason is arithmetic.  The errors at
p >= 11 are large because |w(p)| is small, and |w(p)| is small partly
because p sits far up a d-range that ends at D = floor((N-1)/K).
Raising N raises D, so the same p moves to a relatively lower place in
the range and carries more weight.  Three more N -- 6.4e6, 1.28e7,
2.56e7 -- take D from 728 to about 1815.

**And the buying is itself tested.**  HH3 asks whether the error at
p = 11 actually falls; if it does not, the extension bought nothing
and HH2's verdict would be about the extension rather than about the
primes.  That guard is the same shape as the fitted-N leg that made
rem:omegafactor's transport failure readable, and it is written into
the rule rather than discovered afterwards.

BACKS: Remark {#rem:driftpower} in paper/wall_v3.md.

PRE-REGISTERED PREDICTIONS (written before this script was run)

  HH1 THE GATE.  A at N = 200000 reproduces rem:maintermremoval's
      POINT mainA marker to a relative 1e-12, and the drift of c(3)
      refitted on rem:driftbyprime's original eight N alone
      reproduces its -0.011551 to six decimals.
  HH2 **THE QUESTION.**  On the eleven N, p = 11 and p = 13 both have
      a drift resolved at |t| above 3.
  HH3 **THE GUARD.**  The standard error of the p = 11 drift falls by
      at least 30 per cent against its eight-N value of 0.001154.
  HH4 And the two resolved primes stay put: the eleven-N drifts of
      c(3) and c(7) are within 0.002 of their eight-N values.

REFUTATION RULE (fixed before the run)

  HH1 REFUTED outside either tolerance; nothing below is reported.
  HH2 **REFUTED if either fails to resolve.**  Then more power did
      not bring them in, and -- provided HH3 holds so that power was
      really bought -- 3 and 7 are special rather than merely
      loudest.  That is the reading, and it is only available with
      HH3.
  HH3 **REFUTED below a 30 per cent fall, and that voids HH2's
      reading.**  The extension would have bought no power, so a
      failure at 11 and 13 would say nothing about them; the verdict
      stands and the reading is barred.
  HH4 REFUTED outside 0.002.  Then the drift is not the stable linear
      trend both remarks have fitted, the longer field moves it, and
      the slope model itself needs revisiting before either reading
      is used.

  **THE UNRESOLVED CASE, NAMED, WITH A NUMBER.**  A drift resolves or
  not against its own error, and this run prints both for every prime
  at both field lengths so the change is checkable rather than
  asserted.  **If the p = 11 error falls but its drift also falls, so
  that |t| stays under 3 while both shrink, that is not a failure to
  resolve -- it is a measurement that the drift is smaller than the
  eight-N fit suggested**, and the remark must report the drift's own
  movement beside its t rather than only the verdict.

  WHAT THIS CANNOT DO.  One radical family; eleven N over 3.0
  decades.  A prime that resolves here is resolved on this field and
  not shown to drift asymptotically, and one that does not is not
  shown to be still.  Nothing here bounds anything, and
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
OUT = os.path.join(ROOT, "results", "audit_drift_power.txt")
SRCM = os.path.join(ROOT, "results", "audit_mainterm_removal.txt")
SRCD = os.path.join(ROOT, "results", "audit_drift_by_prime.txt")

THETA = 0.56
OLD = [25_000, 50_000, 100_000, 200_000, 400_000, 800_000,
       1_600_000, 3_200_000]
NEW = [6_400_000, 12_800_000, 25_600_000]
NS = OLD + NEW
NGATE = 200_000
RELID = 1e-12
DEC = 6
TCAP = 3.0
ERRFALL = 30.0
STABLE = 0.002
NOISEFRAC = 0.001
WATCH = (11, 13)


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
    src = io.open(SRCD, encoding="utf-8").read()
    out = {}
    for p in (3, 7, 11, 13):
        d = re.search(r"^POINT drift_%d ([-+]?[\d.]+)\s*$" % p,
                      src, re.M)
        e = re.search(r"^SPREAD drift_%d ([\d.]+)\s*$" % p, src, re.M)
        if not d or not e:
            raise SystemExit("no drift markers for p = %d" % p)
        out[p] = (float(d.group(1)), float(e.group(1)))
    if not m:
        raise SystemExit("no mainA marker")
    return float(m.group(1)), out


HEAD = [
    "STATISTIC: the drift of each prime's weight correction c(p, N)",
    "           fitted on eight N and again on eleven, the change in",
    "           its standard error, and whether p = 11 and p = 13",
    "           resolve when the field is longer.",
    "FIELD: N = %s; the first eight are rem:driftbyprime's and the"
    % NS,
    "       last three are new. d over the squarefree d <= D coprime",
    "       to N, D = floor((N-1)/K), K = floor(N^%.2f). A at N = %d"
    % (THETA, NGATE),
    "       and the eight-N drifts are READ from",
    "       results/audit_mainterm_removal.txt and",
    "       results/audit_drift_by_prime.txt.",
    "",
]


def main():
    lines = []

    def say(t=""):
        print(t)
        sys.stdout.flush()
        lines.append(t)

    pubA, pubd = read_pub()
    say("READ audit_mainterm_removal.txt %d %.17e" % (NGATE, pubA))
    for p in sorted(pubd):
        say("READ audit_drift_by_prime.txt drift_%d %.6f"
            % (p, pubd[p][0]))
        say("READ audit_drift_by_prime.txt sderr_%d %.6f"
            % (p, pubd[p][1]))
        say("  the eight-N drift and its standard error for p = %d"
            % p)
    say("PRINTBOUND audit_drift_power %d %.10f"
        % (DEC, 0.5 * 10.0 ** (-DEC)))
    say("  theta %.2f, |t| cap %.1f, error fall %.0f per cent, "
        "stability %.3f" % (THETA, TCAP, ERRFALL, STABLE))
    say("  watching p = %s" % (WATCH,))

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
        noise_by_N[N] = NOISEFRAC * abs(ws[3])
        cs_by_N[N] = {d: (ws[d] / model(d, D), abs(ws[d]))
                      for d in ws if d != D}
        say("  N = %-10d D = %-5d %d contributing d" % (N, D, len(ws)))
    say("SCALES %d" % len(NS))

    def common(sub):
        c = None
        for N in sub:
            here = {p for p, cw in cs_by_N[N].items()
                    if len(factor_set(p)) == 1
                    and cw[1] >= noise_by_N[N]}
            c = here if c is None else (c & here)
        return sorted(c)

    def drifts(sub):
        x = np.array([math.log(N) for N in sub])
        return {p: fit(x, [cs_by_N[N][p][0] for N in sub])
                for p in common(sub)}

    d8 = drifts(OLD)
    d11 = drifts(NS)

    # ------------------------------------------------------------- HH1
    say()
    say("HH1  the gate")
    ra = abs(A200 - pubA) / max(abs(pubA), 1.0)
    ok3 = abs(round(d8[3][0], DEC)
              - round(pubd[3][0], DEC)) < 10.0 ** (-DEC)
    hh1 = ra <= RELID and ok3
    say("  A relative %.2e against %.0e   %s"
        % (ra, RELID, "ok" if ra <= RELID else "MISMATCH"))
    say("  c(3) drift on the original eight %+.6f against its %+.6f"
        "  %s" % (d8[3][0], pubd[3][0], "ok" if ok3 else "MISMATCH"))
    say("  HH1 %s" % ("hold" if hh1 else "REFUTED"))
    if not hh1:
        io.open(OUT, "w", encoding="utf-8", newline="\n").write(
            "\n".join(HEAD + lines) + "\n")
        raise SystemExit(1)

    say()
    say("      p    eight-N drift    s.e.      t       "
        "eleven-N drift   s.e.      t")
    for p in sorted(set(d8) & set(d11)):
        b8, e8 = d8[p]
        b11, e11 = d11[p]
        say("  %5d   %+.6f  %.6f  %+7.2f    %+.6f  %.6f  %+7.2f"
            % (p, b8, e8, b8 / e8, b11, e11, b11 / e11))
        say("POINT drift11_%d %.6f" % (p, b11))
        say("SPREAD drift11_%d %.6f" % (p, e11))

    # ------------------------------------------------------------- HH3
    say()
    say("HH3  the guard: did the extension buy power?")
    e11 = d11[11][1] if 11 in d11 else float("nan")
    fall = 100.0 * (pubd[11][1] - e11) / pubd[11][1]
    hh3 = fall >= ERRFALL
    say("  p = 11 standard error %.6f -> %.6f, a fall of %.2f per cent"
        % (pubd[11][1], e11, fall))
    say("POINT errfall %.4f" % fall)
    say("  HH3 %s   (cap: at least %.0f per cent)"
        % ("hold" if hh3 else "REFUTED", ERRFALL))

    # ------------------------------------------------------------- HH2
    say()
    say("HH2  do 11 and 13 resolve?")
    hh2 = True
    for p in WATCH:
        if p not in d11:
            hh2 = False
            say("  p = %-3d does not contribute above the threshold at "
                "every N" % p)
            continue
        b, e = d11[p]
        ok = abs(b / e) > TCAP
        hh2 &= ok
        say("  p = %-3d drift %+.6f +- %.6f, t %+.2f   %s"
            % (p, b, e, b / e, "resolved" if ok else "not resolved"))
        say("TSTAT drift11_%d %.2f" % (p, b / e))
        b8 = pubd[p][0]
        say("        its eight-N drift was %+.6f, so the drift itself "
            "moved %+.6f" % (b8, b - b8))
    say("  HH2 %s   (cap: |t| above %.1f on both)"
        % ("hold" if hh2 else "REFUTED", TCAP))
    if not hh3:
        say("  VOID: HH3 failed, so the extension bought no power and "
            "HH2's reading")
        say("  is barred, as the rule says")

    # ------------------------------------------------------------- HH4
    say()
    say("HH4  do 3 and 7 stay put?")
    hh4 = True
    for p in (3, 7):
        mv = abs(d11[p][0] - pubd[p][0])
        ok = mv <= STABLE
        hh4 &= ok
        say("  p = %-3d %+.6f -> %+.6f, moved %.6f   %s"
            % (p, pubd[p][0], d11[p][0], mv, "ok" if ok else "OUT"))
        say("POINT move_%d %.6f" % (p, mv))
    say("  HH4 %s   (cap: %.3f)"
        % ("hold" if hh4 else "REFUTED", STABLE))

    # a diagnostic, after the verdicts and predicted by nothing:
    # HH4 capped the absolute movement; the relative one is the
    # informative quantity and it is the same for both resolved primes
    say()
    say("  a diagnostic, after the verdicts and predicted by nothing")
    say("  HH4 capped how far the drifts moved, not by what factor. "
        "Relatively:")
    for p in (3, 7, 11):
        if p not in d11:
            continue
        r = d11[p][0] / pubd[p][0]
        say("    p = %-3d %+.6f -> %+.6f, a factor of %.5f"
            % (p, pubd[p][0], d11[p][0], r))
        say("POINT shrink_%d %.6f" % (p, r))
    r3 = d11[3][0] / pubd[3][0]
    r7 = d11[7][0] / pubd[7][0]
    say("    the two resolved primes shrink by %.5f and %.5f, "
        "differing by %.5f" % (r3, r7, abs(r3 - r7)))
    say("POINT shrinkgap %.6f" % abs(r3 - r7))
    say("  a slope that shrinks when the field lengthens is not a "
        "slope; that the")
    say("  two shrink together says the linear model is failing the "
        "same way at")
    say("  both, and this run measures it without proposing what "
        "replaces it")

    say()
    say("=" * 70)
    say("HH1 %s  HH2 %s  HH3 %s  HH4 %s"
        % tuple("hold" if v else "REFUTED"
                for v in (hh1, hh2, hh3, hh4)))
    say()
    if hh3 and hh2:
        say("power brought them in. the drift is a property of the "
            "corrections as a")
        say("family and rem:driftbyprime's absence at p >= 11 was an "
            "absence of")
        say("power, which that remark said it could not tell from an "
            "absence of")
        say("drift. now it can.")
    elif hh3 and not hh2:
        say("power was bought and they still do not resolve. 3 and 7 "
            "are special")
        say("rather than merely loudest, and the drift is not a "
            "property of the")
        say("family.")
    else:
        say("the extension bought no power, so nothing is learned "
            "about 11 and 13")
        say("and the verdict above is about this extension and not "
            "about them.")

    io.open(OUT, "w", encoding="utf-8", newline="\n").write(
        "\n".join(HEAD + lines) + "\n")
    print("\nwrote %s" % os.path.normpath(OUT))
    return 0


if __name__ == "__main__":
    sys.exit(main())
