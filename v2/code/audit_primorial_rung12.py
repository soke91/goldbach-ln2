# -*- coding: utf-8 -*-
r"""
The thirteenth rung, and the ladder's first out-of-sample test.

WHAT IS AT STAKE

OPEN item 1 is blocked on one thing: the primorial ladder's level
exponent is trusted only to log10 N = 8.3256, and {#rem:densenoise}
showed the noise there is the field's rather than the statistic's, so
"the only way to push it is to raise N". Sixteen cycles of the sign
axis have since shown that N = 1.024e8 costs seconds on this machine.
The ladder's rungs are 30030 * 2^j, so the next one is
30030 * 2^12 = 123002880, at log10 N = 8.0899 -- three tenths of a
decade past the top published rung and the first rung ever computed
after the line was fitted.

That makes it the ladder's first genuine out-of-sample point.
{#rem:primorialrung11} refitted twelve rungs and said that refit is
"what any later extrapolation must use"; this rung was not in it. Two
things can be read from it and one cannot. It can say whether the
level exponent is still above 1/2 by more than the ladder's own
scatter, and whether the ladder has bent. It cannot say anything
about 10^10.6 -- the forecast for theta' -- and none is attempted.

BACKS: Remark {#rem:rung12} in paper/theorem_A.md.

PRE-REGISTERED PREDICTIONS (written before this script was run)

  B1  The control. Rung 11 recomputed here reproduces
      results/audit_primorial_rung11.txt's exponent inside the bound
      its printing forces.
  B2  The new rung clears the barrier by more than the floor: its
      exponent exceeds 1/2 by more than the ladder's r.m.s. scatter.
  B3  The rise continues: the new exponent exceeds rung 11's.
  B4  And the ladder has not bent: the new rung's residual from the
      published twelve-rung line is inside that line's r.m.s.
      residual.

REFUTATION RULE (fixed before the run)

  B1  REFUTED outside the printing bound. Then this is not the
      statistic that remark measured. THIS ONE GATES.
  B2  REFUTED if the margin is at or below the scatter. Then the
      barrier is not cleared at this rung and {#rem:primorialrung11}'s
      reading rests on rungs 10 and 11 alone.
  B3  REFUTED if the exponent does not rise. One rung falling is not
      a bend -- B4 is the test for that -- but a fall would mean the
      rise is not monotone and no rung may be quoted as "the top".
  B4  REFUTED if the residual exceeds the line's r.m.s. Then the
      ladder bends at the first point outside the fit, and every
      extrapolation from it, including {#rem:primorialdense}'s
      surviving shape, loses its basis.

  B1 gates. B2 to B4 are the measurement and do not gate.

  NO NULL IS RUN and none applies. A deterministic curve is located
  against a computed threshold. The coin arms for this statistic were
  run in lab_primorial_ladder.py and lab_primorial_share.py, and the
  scatter they left is the floor B2 is judged against.
"""

import importlib.util
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
CODE = os.path.join(ROOT, "code")
RES = os.path.join(ROOT, "results")
OUT = os.path.join(RES, "audit_primorial_rung12.txt")

CLIM = 4_000_000                    # the fixed Euler bound (G20)


def module(name):
    p = os.path.join(CODE, name + ".py")
    spec = importlib.util.spec_from_file_location(name, p)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


R11 = module("audit_primorial_rung11")
# bound to a bare name so the Euler-product bound is visible to G20
primes_upto = R11.primes_upto
BASE = R11.BASE
CONTROL = BASE * (1 << 11)          # 61501440, the top published rung
NEW = BASE * (1 << 12)              # 123002880


def read_published():
    """rung 11's exponent, the ladder scatter and the twelve-rung line"""
    src = io.open(os.path.join(RES, "audit_primorial_rung11.txt"),
                  encoding="utf-8").read()
    m = re.search(r"^  N = " + str(CONTROL) +
                  r"\s+thr [\d.]+\s+#k \d+\s+beta [\d.]+\s+K\*_R "
                  r"(\d+)\s+exp ([\d.]+)\s*$", src, re.M)
    dec = len(m.group(2).split(".")[1])
    m2 = re.search(r"slope ([+-][\d.]+), r\.m\.s\. residual "
                   r"([\d.]+), standard error ([\d.]+)", src)
    m3 = re.search(r"^FLOOR primorial_rung11 ([\d.]+)\s*$", src, re.M)
    return (int(m.group(1)), float(m.group(2)), dec,
            float(m2.group(1)), float(m2.group(2)), float(m3.group(1)))


def main():
    lines = []

    def say(s=""):
        print(s)
        lines.append(s)

    kpub, epub, dec, slope12, rms12, floor = read_published()
    say("read rung 11 from results/audit_primorial_rung11.txt: "
        "K*_R %d, exponent %.4f" % (kpub, epub))
    say("  with the twelve-rung line's slope %+.6f and r.m.s. "
        "residual %.4f, and the ladder floor %.4f"
        % (slope12, rms12, floor))
    say("  the statistic, the sieve, the k-cap and the Euler bound "
        "are imported from code/audit_primorial_rung11.py")

    ns, ex, scatter = R11.read_ladder()
    say("  and %d published rungs from "
        "results/audit_primorial_rung10.txt, scatter %.4f"
        % (len(ns), scatter))

    say()
    say("the rungs: the control %d at log10 N = %.4f and the new "
        "%d at %.4f" % (CONTROL, math.log10(CONTROL), NEW,
                        math.log10(NEW)))
    say("RADICALS %d"
        % len(set(tuple(sorted(q for q in R11.factor_set(N) if q > 2))
                  for N in (CONTROL, NEW))))

    qs = [int(q) for q in primes_upto(R11.QSIEVE) if q > 2]
    say()
    say("sieving to %d, sieve weight over the odd primes %s"
        % (NEW, ", ".join(map(str, qs))))
    lam, mu = R11.lambda_and_mu(NEW)
    sqf = mu != 0
    vmask = R11.residue_mask(NEW, qs)
    # the Euler products at the fixed bound, built exactly as
    # code/audit_primorial_rung11.py builds them (G20)
    artin, twin = 1.0, 2.0
    assert CLIM == R11.CLIM
    for p in primes_upto(CLIM):
        p = int(p)
        artin *= 1.0 - 1.0 / (p * (p - 1.0))
        if p > 2:
            twin *= 1.0 - 1.0 / (p - 1.0) ** 2
    say("  the Euler products at the fixed bound %d: Artin %.9f, "
        "twin %.9f" % (CLIM, artin, twin))

    got = {}
    say()
    for N in (CONTROL, NEW):
        out = R11.measure(N, lam, mu, sqf, vmask, qs, artin, twin)
        if out is None:
            say("  N = %-12d no crossing below k = %d" % (N, R11.KCAP))
            continue
        kstar, e, bpn, beta, nk = out
        got[N] = (kstar, e)
        say("  N = %-12d thr %.6f  #k %-7d beta %.6f  K*_R %-8d "
            "exp %.4f" % (N, bpn, nk, beta, kstar, e))
        say("BUDGET kstar_R_S1AN_N%d %.6f" % (N, bpn))
    say("RADICALS 1")

    e11, e12 = got[CONTROL][1], got[NEW][1]

    # -------------------------------------------------------------- B1
    say()
    say("B1  the control at rung 11")
    rnd = 0.5 * 10.0 ** (-dec)
    d = abs(e11 - epub)
    b1 = d <= rnd
    say("  exponent here %.4f against the published %.4f, departure "
        "%.6f; the table prints %d decimals, so the bound is %.8f"
        % (e11, epub, d, dec, rnd))
    say("  and K*_R here %d against %d" % (got[CONTROL][0], kpub))
    say("PRINTBOUND audit_primorial_rung12 %d %.8f" % (dec, rnd))
    say("  B1 %s   (cap: the printing bound)"
        % ("hold" if b1 else "REFUTED"))

    # -------------------------------------------------------------- B2
    say()
    say("B2  does the new rung clear the barrier?")
    marg = e12 - 0.5
    b2 = marg > scatter
    say("  the new exponent is %.4f, a margin of %.4f over one half, "
        "against the ladder's scatter %.4f" % (e12, marg, scatter))
    say("MARGIN audit_primorial_rung12 %.4f %.4f" % (marg, scatter))
    if marg <= scatter:
        say("INSIDE FLOOR audit_primorial_rung12")
    say("FLOOR primorial_rung12 %.4f" % scatter)
    say("  B2 %s   (cap: the scatter)" % ("hold" if b2 else "REFUTED"))

    # -------------------------------------------------------------- B3
    say()
    say("B3  does the rise continue?")
    b3 = e12 > e11
    say("  rung 11 %.4f, rung 12 %.4f, difference %+.4f"
        % (e11, e12, e12 - e11))
    say("  B3 %s" % ("hold" if b3 else "REFUTED"))

    # -------------------------------------------------------------- B4
    say()
    say("B4  has the ladder bent?")
    # the published line is fitted against ln N, not log10 N -- its
    # SPREAD says so -- and the same convention is used here so the
    # slopes are comparable. Predictions do not depend on it.
    x = np.log(np.array(ns + [CONTROL], dtype=np.float64))
    y = np.array(ex + [e11])
    a, b = np.polyfit(x, y, 1)
    x0 = math.log(NEW)
    pred = a * x0 + b
    resid = e12 - pred
    b4 = abs(resid) <= rms12
    n = x.size
    rss = float(((y - (a * x + b)) ** 2).sum())
    s2 = rss / (n - 2)
    sxx = float(((x - x.mean()) ** 2).sum())
    sepred = math.sqrt(s2 * (1.0 + 1.0 / n + (x0 - x.mean()) ** 2
                             / sxx))
    say("  the twelve-rung line refitted here: slope %+.6f against "
        "the published %+.6f" % (a, slope12))
    say("  at log N = %.4f it predicts %.4f; measured %.4f, "
        "residual %+.4f" % (x0, pred, e12, resid))
    say("  against the line's own in-sample r.m.s. residual %.4f"
        % rms12)
    say("  B4 %s   (cap: the line's r.m.s.)"
        % ("hold" if b4 else "REFUTED"))
    say()
    say("  and against the width an out-of-sample point should be "
        "judged by,")
    say("  which is what {#rem:primorialrung11} had already diagnosed "
        "for its own P4")
    say("  and which this script's B4 repeated the error of: the "
        "prediction")
    say("  standard error at a new abscissa, s*sqrt(1 + 1/n + "
        "(x0-xbar)^2/Sxx).")
    say("  in-sample r.m.s.        %.4f" % rms12)
    say("  prediction s.e. at x0   %.4f" % sepred)
    say("  measured departure      %.4f   %s"
        % (abs(resid), "inside" if abs(resid) <= sepred
           else "OUTSIDE"))
    say("  that is %.2f prediction standard errors."
        % (abs(resid) / sepred))

    say()
    say("  the thirteen rungs refitted, which is what any later "
        "extrapolation must use:")
    x13 = np.append(x, x0)
    y13 = np.append(y, e12)
    a13, b13 = np.polyfit(x13, y13, 1)
    r13 = y13 - (a13 * x13 + b13)
    rms13 = float(np.sqrt((r13 ** 2).mean()))
    se13 = math.sqrt(float((r13 ** 2).sum() / (x13.size - 2))
                     / float(((x13 - x13.mean()) ** 2).sum()))
    say("  slope %+.6f, r.m.s. residual %.4f, standard error %.6f, "
        "t = %.2f" % (a13, rms13, se13, abs(a13) / se13))
    say("SCATTER slope_audit_primorial_rung12 %.4f" % rms13)
    say("TSTAT slope_audit_primorial_rung12 %.2f" % (abs(a13) / se13))
    say("SPREAD slope_audit_primorial_rung12 %.4f"
        % float(x13.max() - x13.min()))
    if abs(a13) / se13 < 2.0:
        say("UNRESOLVED SIGN slope_audit_primorial_rung12")
    say("  and the r.m.s. against the twelve-rung %.4f: %+.4f"
        % (rms12, rms13 - rms12))
    say("  no crossing is forecast from this. The barrier is behind "
        "the ladder, and the")
    say("  published theta' forecast of {#rem:primorialdense} lives "
        "two and a half decades")
    say("  past this rung; {#rem:shapetrust} and {#rem:shapepower} "
        "are why none is published.")

    say()
    say("=" * 70)
    say("B1 %s  B2 %s  B3 %s  B4 %s"
        % tuple("hold" if v else "REFUTED"
                for v in (b1, b2, b3, b4)))

    head = [
        "STATISTIC: the truncation K*_R at which",
        "           sum_{k<K}(log k)|R(N;k)| first reaches",
        "           S(N)(1-A(N))N, and its exponent log K*_R / log N,",
        "           at N = 30030*2^12 = 123002880 and, as a control,",
        "           at N = 30030*2^11; the margin of the new exponent",
        "           over 1/2 against the ladder's own scatter; the",
        "           residual of the new rung from the twelve-rung",
        "           line published in",
        "           results/audit_primorial_rung11.txt; and the",
        "           thirteen-rung refit.",
        "NULL: none is run and none applies. A deterministic curve is",
        "      located against a computed threshold; there is no",
        "      background to detect against. The coin arms for this",
        "      statistic were run in lab_primorial_ladder.py and",
        "      lab_primorial_share.py, and the scatter they left is",
        "      the floor the margin is judged against here.",
        "FIELD: N = 61501440 and 123002880, the odd radical",
        "       3*5*7*11*13 fixed so the threshold is constant; k",
        "       squarefree and coprime to N with 2 <= k < 100000; m",
        "       odd, squarefree and coprime to k, m < N/k; the sieve",
        "       weight over the odd primes below 30; the Euler",
        "       products at the fixed bound 4000000. One odd radical,",
        "       as the RADICALS line declares. The statistic, the",
        "       sieve, the k-cap and the Euler bound are imported",
        "       from code/audit_primorial_rung11.py; the published",
        "       rungs come from results/audit_primorial_rung10.txt",
        "       and rung 11 with its line from",
        "       results/audit_primorial_rung11.txt.",
        "",
    ]
    io.open(OUT, "w", encoding="utf-8", newline="\n").write(
        "\n".join(head + lines) + "\n")
    print("\nwrote %s" % OUT)
    if not b1:
        raise SystemExit(1)
    return 0


if __name__ == "__main__":
    sys.exit(main())
