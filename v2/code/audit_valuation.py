# -*- coding: utf-8 -*-
r"""
The variable is the odd prime's exponent, not the radical and not the size

WHAT IS AT STAKE

rem:radicalblind refuted the radical as the variable: five bases
sharing {2,19} disagreed at chi-square 14.233 against a cap of 9.49,
and three remarks lost their readings.  It reported a v2 slope at
t = +2.70 and then barred reading it, because a base of fixed size
cannot raise v2 without lowering the odd primes' exponents -- the two
were correlated at -0.993884.

Writing this run made the reason plainer than that.  **In a doubling
family v2 is not a property of the family at all.**  For base
2^a p^b the family is N_j = 2^(a+j) p^b, so v2(N) runs a, a+1, ...,
a+9 within the family while v_p(N) stays at b.  The label "v2" in
rem:radicalblind was the *starting* v2 and every family covered
overlapping ranges of it.  The invariant that actually differs
between two families of one radical is b.

Read that way its own numbers group exactly:

    {2,19}   v19 = 3   +0.066885  +0.069391
             v19 = 2   +0.106752
             v19 = 1   +0.104068  +0.098167
    {2,3,7}  v3 = 1..5  all +0.2489 to +0.2532, errors near 0.002

Size does not explain it there -- 19456 with v19 = 1 sits between
13718 and 27436 with v19 = 3, and its drift is with the other v19 = 1.
But that is a reading of data already seen, so it is registered here
as a prediction and tested on bases chosen to separate v_p from size
on purpose:

    p = 5   20480 (b1)  25600 (b2)  16000 (b3)  20000 (b4)  25000 (b5)
    p = 7   14336 (b1)  28672 (b1)  25088 (b2)  21952 (b3)  19208 (b4)

Sizes are held inside a factor of two while b runs one to five, so the
two are nearly uncorrelated by construction rather than by luck.

BACKS: Remark {#rem:valuation} in paper/wall_v3.md.

PRE-REGISTERED PREDICTIONS (written before this script was run)

  R1  THE GATE.  Bases 25000 and 16000 reproduce rem:basecontrol's
      ten-N drifts +0.138572 and +0.144190 to six decimals.
  R2  **The exponent is the variable.**  For p = 5 the five drifts
      regressed on v_5 give a slope resolved at |t| above 3.
  R3  The same for p = 7, on its own five.
  R4  **And the size is not.**  Regressed on log base instead, both
      give |t| below 2.

REFUTATION RULE (fixed before the run)

  R1  REFUTED outside six decimals; nothing below is reported.
  R2  **REFUTED at |t| of 3 or below.**  Then v_p is not the variable
      either, and this branch has refuted the radical and the
      exponent and has no candidate left -- which is worth more than
      another fit, and must be said in those words.
  R3  REFUTED at |t| of 3 or below.  R2 and R3 are separate: a
      dependence at one prime and not the other would say the effect
      is not a property of the exponent as such, and the remark may
      then report only the prime where it appears.
  R4  **REFUTED at |t| of 2 or above on either prime.**  Then size
      moves the drift too and this design has not separated what it
      was built to separate; no reading of R2 or R3 may be given
      without that qualification.

  **THE UNRESOLVED CASE, NAMED, WITH A NUMBER.**  Five bases give
  three degrees of freedom on a slope, which is little.  **This run
  prints the smallest slope in v_p it could have resolved at |t| = 3**
  -- three times the standard error the fit returns -- and a refuted
  R2 licenses only "no exponent dependence larger than that", never
  "the exponent does not matter".  The same number governs R4 in the
  other direction: an unresolved size slope excludes only sizes
  effects above twice its own error, and that number is printed too.
  A verdict quoted without its power is the defect this branch has
  now made eight times; the number goes beside the word.

  WHAT THIS CANNOT DO.  Two primes, five bases each, one decade and
  a half of base size.  The exponent b and the starting v2 still move
  together within each prime, since the base size is held fixed --
  what is separated here is b from size, not b from a.  Nothing here
  measures |sum a| against any floor or moves item 5's demand.
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
OUT = os.path.join(ROOT, "results", "audit_valuation.txt")
SRC = os.path.join(ROOT, "results", "audit_base_control.txt")
SRCB = os.path.join(ROOT, "results",
                    "audit_radical_blind.txt")

THETA = 0.56
GROUPS = [(5, [20_480, 25_600, 16_000, 20_000, 25_000]),
          (7, [14_336, 28_672, 25_088, 21_952, 19_208])]
GATES = [25_000, 16_000]
NPER = 10
DEC = 6
TVAL = 3.0
TSIZE = 2.0


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


def vp(n, p):
    e = 0
    while n % p == 0:
        n //= p
        e += 1
    return e


def pairs(N, lam, mu, sqf):
    PN = factor_set(N)
    K = int(N ** THETA)
    lk = np.zeros(N, dtype=np.float64)
    l2sq = 0.0
    for k in range(2, K):
        if not sqf[k] or any(k % q == 0 for q in PN):
            continue
        lg = math.log(k)
        ms = np.arange(1, (N - 1) // k + 1, dtype=np.int64)
        for q in factor_set(k):
            ms = ms[ms % q != 0]
        mm = mu[ms].astype(np.float64)
        l2sq += (lg * float((lam[N - ms * k] * mm).sum())) ** 2
        lk[ms * k] += lg * mm
        del ms, mm
    j = np.arange(1, N, dtype=np.int64)
    sa = abs(float((lam[N - j] * lk[1:]).sum()))
    del j, lk
    return sa, math.sqrt(l2sq)


def fit(x, y):
    x = np.asarray(x, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)
    b, a0 = np.polyfit(x, y, 1)
    r = y - (b * x + a0)
    se = math.sqrt(float((r ** 2).sum() / (len(x) - 2))
                   / float(((x - x.mean()) ** 2).sum()))
    return float(b), se


def drift(base, lam, mu, sqf, keep=None):
    xs, ys = [], []
    for jj in range(NPER):
        N = base * (1 << jj)
        sa, l2 = pairs(N, lam, mu, sqf)
        xs.append(math.log(N))
        ys.append(math.log(sa / l2))
    if keep is not None:
        keep.append((xs, ys))
    return fit(xs, ys)


def read_pub():
    src = io.open(SRC, encoding="utf-8").read()
    out = {}
    for b in GATES:
        m = re.search(r"^POINT basedrift_%d ([-+]?[\d.]+)\s*$" % b,
                      src, re.M)
        if not m:
            raise SystemExit("no basedrift marker for %d" % b)
        out[b] = float(m.group(1))
    return out


HEAD = [
    "STATISTIC: the drift of log(|sum a|/l2) against log N for five",
    "           bases at each of two odd primes, with the prime's",
    "           exponent b running 1 to 5 while the base size is held",
    "           inside a factor of two, and that drift regressed on b",
    "           and on log base separately.",
    "FIELD: base * 2^j for j < %d; bases %s at p = %d and %s at"
    % (NPER, GROUPS[0][1], GROUPS[0][0], GROUPS[1][1]),
    "       p = %d. k over the squarefree k < N^%.2f coprime to N;"
    % (GROUPS[1][0], THETA),
    "       j over every index below N. Two drifts are READ from",
    "       results/audit_base_control.txt as the gate.",
    "NOTE: in a doubling family v2(N) = a + j is not fixed along the",
    "      family; the odd prime's exponent b is.",
    "",
]


def main():
    lines = []

    def say(t=""):
        print(t)
        sys.stdout.flush()
        lines.append(t)

    pub = read_pub()
    for b in GATES:
        say("READ audit_base_control.txt %d %.6f" % (b, pub[b]))
    say("  the two drifts this run shares with rem:basecontrol")
    say("PRINTBOUND audit_valuation %d %.10f"
        % (DEC, 0.5 * 10.0 ** (-DEC)))
    say("  theta %.2f, |t| cap on the exponent %.1f, on the size %.1f"
        % (THETA, TVAL, TSIZE))
    say("RADICALS %d" % len(GROUPS))

    NMAX = max(b for _, g in GROUPS for b in g) * (1 << (NPER - 1))
    say("sieving to %d" % NMAX)
    lam, mu = lambda_and_mu(NMAX)
    sqf = mu != 0

    res = {}
    for p, bases in GROUPS:
        say()
        say("p = %d" % p)
        rows = []
        for base in bases:
            kp = []
            b, se = drift(base, lam, mu, sqf, kp)
            rows.append((base, vp(base, p), b, se, kp[0]))
            say("  base %-6d v%d %-3d drift %+.6f +- %.6f"
                % (base, p, vp(base, p), b, se))
            say("POINT valdrift_%d %.6f" % (base, b))
            say("SPREAD valdrift_%d %.6f" % (base, se))
        res[p] = rows
    say("SCALES %d" % (sum(len(g) for _, g in GROUPS) * NPER))

    # -------------------------------------------------------------- R1
    say()
    say("R1  the gate")
    r1 = True
    for p, rows in res.items():
        for base, _, b, _, _ in rows:
            if base in pub:
                g = abs(b - pub[base]) < 10.0 ** (-DEC)
                r1 &= g
                say("  base %-6d here %+.6f against its %+.6f  %s"
                    % (base, b, pub[base], "ok" if g else "MISMATCH"))
    say("  R1 %s   (cap: %d decimals)"
        % ("hold" if r1 else "REFUTED", DEC))
    if not r1:
        io.open(OUT, "w", encoding="utf-8", newline="\n").write(
            "\n".join(HEAD + lines) + "\n")
        raise SystemExit(1)

    # ---------------------------------------------------- R2, R3 and R4
    verd = {}
    for name, (p, _) in zip(("R2", "R3"), GROUPS):
        rows = res[p]
        bb = np.array([float(r[1]) for r in rows])
        dd = np.array([r[2] for r in rows])
        ss = np.array([math.log(float(r[0])) for r in rows])
        sl, se = fit(bb, dd)
        t = sl / se
        ok = abs(t) > TVAL
        verd[name] = ok
        slz, sez = fit(ss, dd)
        tz = slz / sez
        c = float(np.corrcoef(bb, ss)[0, 1])
        say()
        say("%s  does the drift depend on the exponent at p = %d?"
            % (name, p))
        say("  on v%d      slope %+.6f +- %.6f, t %+.2f"
            % (p, sl, se, t))
        say("  on log base slope %+.6f +- %.6f, t %+.2f"
            % (slz, sez, tz))
        say("  corr(v%d, log base) %+.6f" % (p, c))
        say("CORR valuation%d_regressors %.6f" % (p, abs(c)))
        if abs(c) >= 0.99:
            say("COEFF NOT SEPARABLE valuation%d" % p)
        say("TSTAT valexp_%d %.2f" % (p, t))
        say("SPREAD valexp_%d %.6f" % (p, se))
        say("TSTAT valsize_%d %.2f" % (p, tz))
        say("SPREAD valsize_%d %.6f" % (p, sez))
        say("  smallest exponent slope this could resolve at |t| = "
            "%.1f: %.6f" % (TVAL, TVAL * se))
        say("POINT valpower_%d %.6f" % (p, TVAL * se))
        say("  %s %s   (cap: |t| above %.1f)"
            % (name, "hold" if ok else "REFUTED", TVAL))
        verd[("size", p)] = abs(tz) < TSIZE

    say()
    say("R4  and is the size not the variable?")
    r4 = all(verd[("size", p)] for p, _ in GROUPS)
    for p, _ in GROUPS:
        say("  p = %d: size slope %s at the %.1f cap"
            % (p, "unresolved" if verd[("size", p)] else "RESOLVED",
               TSIZE))
    say("  R4 %s   (cap: |t| below %.1f on both)"
        % ("hold" if r4 else "REFUTED", TSIZE))

    # a diagnostic, after the verdicts and predicted by nothing:
    # log(|sum a|/l2) is not a straight line -- rem:deficitdirect put
    # its curvature at -0.007380 and rem:deficitshape found the drift
    # steepening -- so a ten-point slope depends on where its window
    # sits. if one family's own slope moves across its own sub-windows
    # by as much as the bases differ, the base scatter is window noise
    # and not arithmetic.
    say()
    say("  a diagnostic, after the verdicts and predicted by nothing")
    say("  the fitted quantity is a slope of a curved function, so a "
        "window that")
    say("  sits elsewhere returns a different slope. each family's "
        "own first six")
    say("  and last six points, against the spread between bases:")
    for p, _ in GROUPS:
        rows = res[p]
        spans = []
        for base, b, d, se, (xs, ys) in rows:
            lo, _ = fit(xs[:6], ys[:6])
            hi, _ = fit(xs[-6:], ys[-6:])
            spans.append(abs(hi - lo))
            say("    p %d base %-6d full %+.6f   low %+.6f   high "
                "%+.6f   span %.6f" % (p, base, d, lo, hi,
                                       abs(hi - lo)))
        ds = [r[2] for r in rows]
        between = max(ds) - min(ds)
        within = float(np.mean(spans))
        say("    p %d: mean within-family window span %.6f, spread "
            "between bases %.6f" % (p, within, between))
        say("POINT windowspan_%d %.6f" % (p, within))
        say("POINT betweenbase_%d %.6f" % (p, between))
    say("  and the size span of each group of bases, here and in "
        "rem:radicalblind,")
    say("  beside how far that group's drifts disagreed:")
    for p, bases in GROUPS:
        rows = res[p]
        ds = [r[2] for r in rows]
        say("    p %-8d bases span %.3fx   drifts spread %.6f"
            % (p, max(bases) / float(min(bases)), max(ds) - min(ds)))
        say("POINT sizespan_%d %.6f" % (p, max(bases) / float(min(bases))))
    bs = io.open(SRCB, encoding="utf-8").read()
    for tag, key in (("{2,19}", "2_19"), ("{2,3,7}", "2_3_7")):
        got = [int(m) for m in re.findall(
            r"^POINT blinddrift_(\d+) ", bs, re.M)]
        cm = re.search(r"^POINT blindchi_%s ([\d.]+)" % key, bs, re.M)
        grp = [b for b in got if tag == "{2,19}"
               and b in (13718, 19456, 23104, 27436, 38912)
               or tag == "{2,3,7}"
               and b in (15876, 16128, 21168, 21504, 23814)]
        if grp and cm:
            sp = max(grp) / float(min(grp))
            say("    %-10s bases span %.3fx   chi-square %s"
                % (tag, sp, cm.group(1)))
            say("READ audit_radical_blind.txt chi_%s %s"
                % (key, cm.group(1)))
            say("POINT sizespan_%s %.6f" % (key, sp))
    say("  a within-family span of the same size as the between-base "
        "spread")
    say("  would mean the scatter this branch has been reading as "
        "arithmetic is")
    say("  the slope of a curve moving under its own window")

    say()
    say("=" * 70)
    say("R1 %s  R2 %s  R3 %s  R4 %s"
        % (("hold" if r1 else "REFUTED"),
           ("hold" if verd["R2"] else "REFUTED"),
           ("hold" if verd["R3"] else "REFUTED"),
           ("hold" if r4 else "REFUTED")))
    say()
    if verd["R2"] and verd["R3"] and r4:
        say("the variable is the odd prime's exponent. it resolves at "
            "both primes")
        say("with the base size held inside a factor of two and its "
            "own slope")
        say("unresolved, so this is the first variable in this branch "
            "that")
        say("survives a control. what it is not is a mechanism: "
            "nothing here says")
        say("why the exponent of a prime should move a quantity whose "
            "construction")
        say("depends on N only through its radical.")
    elif not r4:
        say("the size moves the drift as well, so this design has not "
            "separated")
        say("what it was built to separate and no exponent reading is "
            "given.")
    elif verd["R2"] or verd["R3"]:
        say("the exponent resolves at one prime and not the other, so "
            "it is not a")
        say("property of the exponent as such. only the prime where "
            "it appears is")
        say("reported and the branch still has no variable that "
            "works at both.")
    else:
        say("the exponent is not the variable either. this branch has "
            "now refuted")
        say("the radical and the exponent and has no candidate left, "
            "which is")
        say("worth more than another fit and is said in those words.")

    io.open(OUT, "w", encoding="utf-8", newline="\n").write(
        "\n".join(HEAD + lines) + "\n")
    print("\nwrote %s" % os.path.normpath(OUT))
    return 0


if __name__ == "__main__":
    sys.exit(main())
