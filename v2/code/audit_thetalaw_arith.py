# -*- coding: utf-8 -*-
r"""Is the theta'-law's crossing a fact about N, or about one family of N?

Supports {#rem:thetalawarith}; audits {#rem:thetalaw}'s U4 and U5.

WHAT IS AT STAKE

{#rem:thetalaw}'s U4/U5 carry item 5's verdict. They fit e(G) against
theta' and find the demand line e(G) = theta'/2 met only at
theta' = 0.4041, below the 1/2 the reduction needs -- with a bracket
[0.3017, 0.5351] whose upper end is inside the admissible region, so
the branch is unresolved rather than closed.

Every N behind that estimate is 2^a 5^b. Both S(N) and A(N) are then
constant across the field, and so is the prime support that fixes the
local densities. A ratio measured on a constant-arithmetic field has
been shown once already in this repository to be a property of the
family and not of the range: the demand ratio read 1.6 there and 40 at
primorial N in the same range.

This run asks the same question of the number that decides item 5.

WHAT IS MEASURED

For each arithmetic family F and each theta' in THETAS: over the N in F
inside a fixed range, with a_k = (log k) H(N;k) over squarefree
k < N^theta' coprime to N,

    l1(N)  = sum_k |a_k|,       |sum a|(N) = |sum_k a_k|,
    e(l1), e(|sum a|)           fitted in log N,
    e(G)   = e(l1) - e(|sum a|),

then e(G) is fitted against theta' and crossed with the derived ceiling
theta'/2.  The reported quantity is the crossing point per family.

FIELDS.  Four, all inside the same [LO, HI]:
    A  2^a 5^b            -- the published field, reproduced
    B  2^a 3^b            -- one small prime swapped
    C  2 * squarefree odd -- prime support varies point to point
    D  primorial-like     -- N divisible by 2,3,5,7,...

FALSIFICATION, registered before the run

  T1  REFUTED if family A's crossing does not reproduce the published
      0.4041 to within its own published bracket. Then this
      implementation is not measuring the published quantity and
      nothing below is interpretable.

  T2  the claim under test: the crossing is a property of the range.
      REFUTED if the crossings across the four families span more than
      0.05 in theta'. Then the published 0.4041 is one family's number.

  T3  REFUTED if any family's crossing lands at or above 0.50. Then the
      branch is open at the admissible end for that arithmetic, and
      item 5 is not closed by this axis.

  PREDICTION.  T2 is refuted and T3 is not: the crossings move with the
  arithmetic -- the local factors that set l1's size differ by family --
  but all stay below 1/2, so the direction of {#rem:thetalaw}'s verdict
  survives while its precision does not.

  A crossing is an extrapolation past the swept theta' whenever it lands
  outside THETAS; it is reported with the fit's own bracket and is not
  read as a measurement.

NULL.  None applies: every quantity is a deterministic sum over a fixed
finite set. There is no sampling and no sign input, so there is nothing
for a null to preserve. The control here is the family comparison
itself -- four fields measured the same way, which is what M4 asks for.
"""
import io
import importlib.util
import math
import os
import sys

import numpy as np

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(os.path.dirname(HERE), "results")

LO, HI = 200_000, 12_800_000
THETAS = (0.40, 0.46, 0.52, 0.56, 0.60, 0.64)
CEILING = 0.5           # the derived ceiling is theta'/2
PUBLISHED_CROSS = 0.4041
PUBLISHED_BRACKET = (0.3017, 0.5351)
SPAN_CAP = 0.05         # T2
ADMISSIBLE = 0.50       # T3


def module(name):
    p = os.path.join(HERE, name + ".py")
    spec = importlib.util.spec_from_file_location(name, p)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


SPL = module("audit_gain_split")


def fam_2a5b(lo, hi):
    out = []
    a = 1
    while 2 ** a <= hi:
        b = 1
        while 2 ** a * 5 ** b <= hi:
            v = 2 ** a * 5 ** b
            if v >= lo:
                out.append(v)
            b += 1
        a += 1
    return sorted(set(out))


def fam_2a3b(lo, hi):
    out = []
    a = 1
    while 2 ** a <= hi:
        b = 1
        while 2 ** a * 3 ** b <= hi:
            v = 2 ** a * 3 ** b
            if v >= lo:
                out.append(v)
            b += 1
        a += 1
    return sorted(set(out))


def fam_twice_squarefree(lo, hi, mu, want=14):
    """N = 2*q, q odd squarefree -- prime support varies point to point.

    Spread over log N so the fit has the same lever arm as the others.
    """
    out, seen = [], set()
    lg0, lg1 = math.log(lo), math.log(hi)
    for i in range(want):
        target = int(math.exp(lg0 + (lg1 - lg0) * i / (want - 1.0)))
        t = target if target % 2 == 0 else target + 1
        for step in range(0, 40000, 2):
            for v in (t + step, t - step):
                if not (lo <= v <= hi) or v in seen:
                    continue
                q = v // 2
                if q % 2 == 1 and mu[q] != 0:
                    out.append(v)
                    seen.add(v)
                    break
            else:
                continue
            break
    return sorted(set(out))


def fam_primorial_like(lo, hi, want=14):
    """N carrying many small primes: multiples of a primorial."""
    prims = [2, 6, 30, 210, 2310, 30030]
    out = set()
    lg0, lg1 = math.log(lo), math.log(hi)
    for i in range(want):
        target = math.exp(lg0 + (lg1 - lg0) * i / (want - 1.0))
        best = None
        for P in prims:
            j = max(1, int(round(target / P)))
            v = P * j
            if not (lo <= v <= hi):
                continue
            score = (-len(SPL.factor_set(v)), abs(math.log(v) - math.log(target)))
            if best is None or score < best[0]:
                best = (score, v)
        if best:
            out.add(best[1])
    return sorted(out)


def fit(x, y):
    a, b = np.polyfit(x, y, 1)
    r = y - (a * x + b)
    se = math.sqrt(float((r ** 2).sum() / (x.size - 2))
                   / float(((x - x.mean()) ** 2).sum())) \
        if x.size > 2 else float("inf")
    return float(a), float(b), float(np.sqrt((r ** 2).mean())), float(se)


def crossing(slope, icept):
    """where e(G) = slope*t + icept meets t/2"""
    den = CEILING - slope
    return float("nan") if abs(den) < 1e-12 else icept / den


def main():
    lines = []

    def say(s=""):
        print(s)
        lines.append(s)

    say(__doc__.strip())
    say()
    say("=" * 70)

    base = SPL.THETA
    say("sieving to %d ..." % HI)
    lam, mu = SPL.lambda_and_mu(HI)
    sqf = mu != 0

    fams = [
        ("A  2^a 5^b", fam_2a5b(LO, HI)),
        ("B  2^a 3^b", fam_2a3b(LO, HI)),
        ("C  2*sqfree", fam_twice_squarefree(LO, HI, mu)),
        ("D  primorial", fam_primorial_like(LO, HI)),
    ]

    say()
    say("STATISTIC: for a_k = (log k) H(N;k) over squarefree k < N^theta'")
    say("           coprime to N, the exponents e(l1) and e(|sum a|)")
    say("           fitted in log N, their difference e(G), then e(G)")
    say("           fitted against theta' and crossed with theta'/2.")
    say("           The reported number is that crossing, per family.")
    say("FIELD: four arithmetic families inside [%d, %d]; k squarefree"
        % (LO, HI))
    say("       with 2 <= k < N^theta', (k,N) = 1; m over 1 <= m < N/k")
    say("       with (m,k) = 1. Lambda and mu from one integer sieve.")
    say("CONSTANTS: LO = %d, HI = %d, THETAS = %s, ceiling slope = %.1f,"
        % (LO, HI, ",".join("%.2f" % t for t in THETAS), CEILING))
    say("           T2 span cap = %.2f, T3 admissible edge = %.2f"
        % (SPAN_CAP, ADMISSIBLE))
    say("NULL: none applies -- deterministic sums over fixed finite sets,")
    say("      no sampling, no sign input. The control is the family")
    say("      comparison itself (M4).")
    say("DENOM: no relative errors are printed.")
    say()

    for nm, NS in fams:
        say("  %-13s %2d points   P(N) classes: %d   N: %d .. %d"
            % (nm, len(NS),
               len(set(tuple(sorted(SPL.factor_set(v))) for v in NS)),
               NS[0], NS[-1]))
    say()

    results = {}
    for nm, NS in fams:
        x = np.log(np.array(NS, dtype=float))
        egs, ses, ratios = [], [], []
        for th in THETAS:
            SPL.THETA = th
            l1s, sas = [], []
            for N in NS:
                ks, a = SPL.weighted(N, lam, mu, sqf)
                l1s.append(float(np.abs(a).sum()))
                sas.append(abs(float(a.sum())))
            el = fit(x, np.log(np.array(l1s)))
            ea = fit(x, np.log(np.array(sas)))
            eg = fit(x, np.log(np.array(l1s) / np.array(sas)))
            egs.append(eg[0])
            ses.append(eg[3])
            # the model-free number: 1 means no cancellation at all
            ratios.append(float(np.mean(np.array(sas) / np.array(l1s))))
        sl = fit(np.array(THETAS, dtype=float), np.array(egs))
        cr = crossing(sl[0], sl[1])
        lo3 = fit(np.array(THETAS[:3], dtype=float), np.array(egs[:3]))
        hi3 = fit(np.array(THETAS[3:], dtype=float), np.array(egs[3:]))
        b = sorted([crossing(lo3[0], lo3[1]), crossing(hi3[0], hi3[1])])
        # degenerate: no cancellation across dilations at all, so
        # |sum a| = l1 identically and e(G) is 0 at every theta'.
        # A crossing is then meaningless -- there is no line to cross
        # the ceiling with.
        deg = max(abs(v) for v in egs) < 1e-9
        results[nm] = (egs, ses, sl, cr, b, deg, ratios[:])

    say("e(G) at each theta', by family")
    say("  theta'   " + "".join("%-14s" % nm.split()[0] for nm, _ in fams))
    for i, th in enumerate(THETAS):
        say("  %-8.2f " % th
            + "".join("%+-13.6f " % results[nm][0][i] for nm, _ in fams))
    say()
    say("  gap e(G) - theta'/2")
    say("  theta'   " + "".join("%-14s" % nm.split()[0] for nm, _ in fams))
    for i, th in enumerate(THETAS):
        say("  %-8.2f " % th
            + "".join("%+-13.6f " % (results[nm][0][i] - th / 2.0)
                      for nm, _ in fams))
    say()

    say("  |sum a| / l1 -- 1 means no cancellation across dilations")
    say("  theta'   " + "".join("%-14s" % nm.split()[0] for nm, _ in fams))
    for i, th in enumerate(THETAS):
        say("  %-8.2f " % th
            + "".join("%-13.8f " % results[nm][6][i] for nm, _ in fams))
    say()

    say("crossing with the derived ceiling theta'/2")
    say("  family          slope       intercept   crossing   bracket")
    crs, degs = [], []
    for nm, _ in fams:
        egs, ses, sl, cr, b, deg, rat = results[nm]
        if deg:
            degs.append(nm)
            say("  %-13s  e(G) = 0 at every theta' -- DEGENERATE, no"
                " crossing exists" % nm)
            say("DEGENERATE %s |sum a|/l1 = %.8f at theta' = 0.56"
                % (nm.split()[0], rat[THETAS.index(0.56)]))
            continue
        crs.append(cr)
        flag = "" if b[0] <= cr <= b[1] else "   (crossing outside its own bracket)"
        say("  %-13s %+-11.6f %+-11.6f %-10.4f [%.4f, %.4f]%s"
            % (nm, sl[0], sl[1], cr, b[0], b[1], flag))
        say("BRACKET crossing_%s %.4f %.4f %.4f"
            % (nm.split()[0], cr, b[0], b[1]))
        say("DRIFT crossing_%s %.6f" % (nm.split()[0], abs(b[1] - b[0])))
    say()

    span = max(crs) - min(crs)
    a_cr = crs[0]
    t1 = PUBLISHED_BRACKET[0] <= a_cr <= PUBLISHED_BRACKET[1]
    t2 = span <= SPAN_CAP
    live = [nm for nm, _ in fams if not results[nm][5]]
    t3_hits = [nm for nm, c in zip(live, crs) if c >= ADMISSIBLE]

    say("T1  family A reproduces the published crossing")
    say("    published %.4f, bracket [%.4f, %.4f]; here %.4f"
        % (PUBLISHED_CROSS, PUBLISHED_BRACKET[0], PUBLISHED_BRACKET[1],
           a_cr))
    say("    T1 %s" % ("hold" if t1 else "REFUTED"))
    say()
    say("T2  the crossing is a property of the range, not the family")
    say("    crossings %s" % ", ".join("%.4f" % c for c in crs))
    say("    span %.4f against cap %.4f" % (span, SPAN_CAP))
    say("    T2 %s" % ("hold" if t2 else "REFUTED"))
    say()
    say("T3  no family crosses at or above %.2f" % ADMISSIBLE)
    if t3_hits:
        say("    at or above: %s" % ", ".join(t3_hits))
    else:
        say("    none reaches %.2f" % ADMISSIBLE)
    say("    T3 %s" % ("hold" if not t3_hits else "REFUTED"))
    say()
    if degs:
        say("DEGENERATE FAMILIES  %s" % ", ".join(degs))
        say("    For these the terms all carry one sign, so")
        say("    |sum a| = l1 exactly and e(G) is 0 at every theta'.")
        say("    There is no cancellation across dilations to measure,")
        say("    and no theta' brings the ceiling down to meet a flat")
        say("    zero. The gap there is the whole ceiling, not a margin.")
        say()
    say("SPREAD crossing_across_families %.4f" % span)
    say("COPRIME thetalaw_arith 4")
    say("FIELD SPLIT four arithmetic families, measured identically")
    say()
    say("=" * 70)
    say("T1 %s  T2 %s  T3 %s"
        % ("hold" if t1 else "REFUTED",
           "hold" if t2 else "REFUTED",
           "hold" if not t3_hits else "REFUTED"))

    SPL.THETA = base
    io.open(os.path.join(RES, "audit_thetalaw_arith.txt"), "w",
            encoding="utf-8", newline="\n").write("\n".join(lines) + "\n")
    return 0 if t1 else 1


if __name__ == "__main__":
    raise SystemExit(main())
