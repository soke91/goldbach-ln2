# -*- coding: utf-8 -*-
r"""
The curvature's out-of-sample errors have a sign. Is that chance?

WHAT IS AT STAKE

[rem:rung18] recorded six out-of-sample departures of the quadratic
-- +0.0046, +0.0023, -0.0012, -0.0023, -0.0033, -0.0041 -- and noted
that the last four are negative and growing, the last one at 0.95 of
its own prediction error. Each is individually inside its error, so
no registered test has failed. But a one-signed run that grows is
what a shape looks like when it begins to fail from above, and if
that is what this is, then the ladder is flattening against the
quadratic and every level above 0.56 is crossed later than the
quadratic places it -- including the 0.56 bracket's own upper end,
which [rem:rung18] left unresolved.

Nothing has measured it, and the six numbers as published cannot:
**they are not one series.** The first four were computed on the
ladder at the published cap and the last two on the uniform ladder of
[rem:laddercap]. Reading a run off a list that changes definition in
the middle is exactly the error [rem:rung16] found in the cap.

So the series is rebuilt homogeneously here, on the uniform ladder
alone, by walking forward: for each j, fit the quadratic on rungs
0..j-1 only and predict rung j. That is the same construction the
rung remarks used at the top, applied at every rung where it can be,
and it yields a series long enough to test rather than six numbers
from two ladders.

The null is the one that matters for the question asked. If the
quadratic is the right shape and the scatter is noise, then walk-
forward departures are a mean-zero sequence and long one-signed runs
are rare. So ladders are simulated from the fitted quadratic with
Gaussian noise at the ladder's own r.m.s., walked forward the same
way, and the terminal run measured on each. The question is what
fraction of those reach the observed run.

No new sieving is done. Every exponent is read from a result file.

BACKS: Remark {#rem:signrun} in paper/theorem_A.md.

PRE-REGISTERED PREDICTIONS (written before this script was run)

  S1  The control. Walking forward reproduces what the rung remarks
      published at the two rungs where they did it: fitted on rungs
      0..16 the quadratic departs from rung 17 by -0.0033 with
      prediction error 0.0044, and fitted on 0..17 it departs from
      rung 18 by -0.0041 with prediction error 0.0043, to the
      decimals those remarks print.
  S2  The run is real and at least four. The walk-forward departures
      end in a run of at least four consecutive negatives.
  S3  The run is not chance. Under the null -- the fitted quadratic
      plus Gaussian noise at the ladder's own r.m.s., walked forward
      identically -- fewer than 5 per cent of simulated ladders end
      in a same-signed run as long as the observed one.
  S4  The ladder is flattening. Regressing the walk-forward
      departures on log N gives a negative slope resolved at |t| > 2.

REFUTATION RULE (fixed before the run)

  S1  REFUTED if either departure or either prediction error differs
      in the decimals published. Then this is not the construction
      the remarks used and nothing below is about their numbers.
      THIS ONE GATES.
  S2  REFUTED if the terminal run is shorter than four. Then the run
      read off the published six was an artefact of splicing two
      ladders, and {#rem:rung18}'s closing paragraph is withdrawn.
  S3  REFUTED at 5 per cent or above. That is the outcome that
      keeps the quadratic: runs of this length would be ordinary for
      a correct shape with this much noise and this few points, and
      the departures would carry no information about the shape
      failing. **This is the prediction most likely to be refuted
      and the reason the null is run rather than argued.**
  S4  REFUTED if the slope is positive, or if |t| <= 2. A run can be
      one-signed without drifting; if the slope is unresolved then
      "growing" is not established and only the sign is, which is a
      weaker statement than {#rem:rung18} made.

  S1 gates. S2 to S4 are the measurement and do not gate.

  THE NULL IS RUN, and it is S3. The arm is ladders drawn from the
  fitted quadratic at the observed abscissae with i.i.d. Gaussian
  errors of the ladder's r.m.s. residual, walked forward by the same
  code path that produces the observed series.
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
OUT = os.path.join(RES, "audit_ladder_signrun.txt")

SEED = 20260823
DRAWS = 20000
START = 6                           # the first rung predicted forward
BASE = 30030
DEC = 4                             # the decimals the rung remarks print


def read_ladder():
    """the nineteen uniform rungs at cap 10^6, from three result files"""
    src = io.open(os.path.join(RES, "audit_ladder_cap.txt"),
                  encoding="utf-8").read()
    ns, ex = [], []
    for ln in src.splitlines():
        m = re.match(r"^  (\d+)\s+(\d+)\s+(.*)$", ln)
        if not m:
            continue
        f = [t for t in m.group(3).split() if t != "cap-invariant"]
        if len(f) < 8 or f[6] == "none":
            continue
        ns.append(int(m.group(2)))
        ex.append(float(f[7]))
    pub = {}
    for j, name in ((17, "audit_primorial_rung17.txt"),
                    (18, "audit_primorial_rung18.txt")):
        s = io.open(os.path.join(RES, name), encoding="utf-8").read()
        N = BASE * (1 << j)
        m = re.search(r"^  N = " + str(N) +
                      r"\s+thr [\d.]+\s+#k \d+\s+beta [\d.]+\s+"
                      r"K\*_R \d+\s+exp ([\d.]+)\s*$", s, re.M)
        ns.append(N)
        ex.append(float(m.group(1)))
        d = re.search(r"^  quadratic\s+[\d.]+\s+[\d.]+\s+"
                      r"([-+][\d.]+)\s+([\d.]+)", s, re.M)
        pub[j] = (float(d.group(1)), float(d.group(2)))
    return ns, ex, pub


def quadfit(x, y):
    A = np.column_stack([np.ones_like(x), x, x * x])
    c, *_ = np.linalg.lstsq(A, y, rcond=None)
    r = y - A.dot(c)
    s2 = float((r ** 2).sum()) / (x.size - 3)
    cov = s2 * np.linalg.inv(A.T.dot(A))
    return c, cov, s2, float(np.sqrt((r ** 2).mean()))


def walk(x, y, start=START):
    """fit on 0..j-1, predict j: the departure and its own error"""
    dep, se, at = [], [], []
    for j in range(start, x.size):
        c, cov, s2, _ = quadfit(x[:j], y[:j])
        v = np.array([1.0, x[j], x[j] * x[j]])
        p = float(v.dot(c))
        dep.append(y[j] - p)
        se.append(math.sqrt(s2 + float(v.dot(cov).dot(v))))
        at.append(x[j])
    return np.array(dep), np.array(se), np.array(at)


def tailrun(d):
    """how many of the last departures share one sign"""
    if d.size == 0:
        return 0
    s = np.sign(d[-1])
    n = 0
    for v in d[::-1]:
        if np.sign(v) != s:
            break
        n += 1
    return n


def main():
    lines = []

    def say(t=""):
        print(t)
        sys.stdout.flush()
        lines.append(t)

    ns, ex, pub = read_ladder()
    x = np.log(np.array(ns, dtype=np.float64))
    y = np.array(ex)
    say("read %d uniform rungs at cap 1000000 from "
        "results/audit_ladder_cap.txt," % len(ns))
    say("  results/audit_primorial_rung17.txt and rung18.txt")
    say("SEED %d" % SEED)
    say("DRAWS %d" % DRAWS)

    dep, se, at = walk(x, y)
    run = tailrun(dep)

    # -------------------------------------------------------------- S1
    say()
    say("S1  the control at the two rungs the remarks walked forward")
    s1 = True
    rnd = 0.5 * 10.0 ** (-DEC)
    for j in (17, 18):
        i = j - START
        pd, ps = pub[j]
        ok = abs(dep[i] - pd) <= rnd and abs(se[i] - ps) <= rnd
        s1 = s1 and ok
        say("  rung %d  departure %+.4f against the published %+.4f, "
            "error %.4f against %.4f, %s"
            % (j, dep[i], pd, se[i], ps, "equal" if ok else "DIFFERENT"))
    say("PRINTBOUND audit_ladder_signrun %d %.8f" % (DEC, rnd))
    say("  S1 %s   (cap: the decimals those remarks print)"
        % ("hold" if s1 else "REFUTED"))
    if not s1:
        io.open(OUT, "w", encoding="utf-8", newline="\n").write(
            "\n".join(lines) + "\n")
        raise SystemExit(1)

    # -------------------------------------------------------------- S2
    say()
    say("S2  the walk-forward series, on one ladder")
    say("  rung  log10 N   measured   predicted  departure  own error"
        "  ratio")
    for i in range(dep.size):
        j = START + i
        say("  %-5d %-9.4f %-10.6f %-10.6f %+-10.4f %-10.4f %.2f"
            % (j, at[i] / math.log(10.0), y[j], y[j] - dep[i],
               dep[i], se[i], abs(dep[i]) / se[i]))
    say("SIGNRUN ladder_walkforward %d" % run)
    say("  the departures end in a run of %d with sign %+d, out of "
        "%d points" % (run, int(np.sign(dep[-1])), dep.size))
    s2 = run >= 4
    say("  S2 %s   (cap: a run of four)"
        % ("hold" if s2 else "REFUTED"))

    # -------------------------------------------------------------- S3
    say()
    say("S3  the null: is a run that long ordinary for a right shape?")
    c, cov, s2v, rms = quadfit(x, y)
    A = np.column_stack([np.ones_like(x), x, x * x])
    fit = A.dot(c)
    rng = np.random.default_rng(SEED)
    hits = 0
    for _ in range(DRAWS):
        yy = fit + rng.normal(0.0, rms, size=x.size)
        d2, _s, _a = walk(x, yy)
        if tailrun(d2) >= run:
            hits += 1
    frac = hits / float(DRAWS)
    s3 = frac < 0.05
    say("  the fitted quadratic is %+.8f in (log N)^2 with r.m.s. "
        "residual %.4f" % (c[2], rms))
    say("  %d of %d null ladders end in a run of %d or longer: "
        "%.4f" % (hits, DRAWS, run, frac))
    say("NULL ladder_signrun %.4f" % frac)
    say("  S3 %s   (cap: 0.05)" % ("hold" if s3 else "REFUTED"))

    # -------------------------------------------------------------- S4
    say()
    say("S4  is the run drifting, or only one-signed?")
    sl, ic = np.polyfit(at, dep, 1)
    n = at.size
    resid = dep - (sl * at + ic)
    sse = float((resid ** 2).sum()) / (n - 2)
    sxx = float(((at - at.mean()) ** 2).sum())
    sse_sl = math.sqrt(sse / sxx)
    t = sl / sse_sl
    s4 = sl < 0.0 and abs(t) > 2.0
    say("  the departure slope against log N is %+.8f +- %.8f, "
        "t = %.2f" % (sl, sse_sl, t))
    say("TSTAT ladder_signrun_slope %.2f" % t)
    if abs(t) < 2.0:
        say("UNRESOLVED SIGN ladder_signrun_slope")
    say("SPREAD ladder_signrun_slope %.4f" % (at.max() - at.min()))
    say("SCATTER slope_audit_ladder_signrun %.4f"
        % float(np.sqrt((resid ** 2).mean())))
    say("  S4 %s   (cap: negative and |t| > 2)"
        % ("hold" if s4 else "REFUTED"))

    say()
    say("=" * 70)
    say("S1 %s  S2 %s  S3 %s  S4 %s"
        % tuple("hold" if v else "REFUTED" for v in (s1, s2, s3, s4)))

    head = [
        "STATISTIC: the out-of-sample departure of the quadratic at",
        "           each rung of the uniform ladder, fitted on the",
        "           rungs below it only; the length of the run of",
        "           one sign the series ends in; the fraction of null",
        "           ladders reaching that run; and the slope of the",
        "           departures against log N with its t.",
        "NULL: RUN, and it is S3. Ladders are drawn from the",
        "      quadratic fitted on all nineteen rungs, at the same",
        "      abscissae, with i.i.d. Gaussian errors of the ladder's",
        "      own r.m.s. residual, and walked forward by the same",
        "      code path as the observed series. The arm asks how",
        "      often a correct shape with this noise and this many",
        "      points produces a terminal run as long as observed.",
        "FIELD: N = 30030*2^j for j = 0..18, the odd radical",
        "       3*5*7*11*13 fixed along the ladder; the level",
        "       exponent log K*_R / log N at the uniform cap 1000000",
        "       as printed in results/audit_ladder_cap.txt for",
        "       j = 0..16, results/audit_primorial_rung17.txt for",
        "       j = 17 and results/audit_primorial_rung18.txt for",
        "       j = 18. No sieving is done here; every exponent is",
        "       read. The walk forward starts at j = 6, the first",
        "       rung with enough points below it for a quadratic and",
        "       a residual variance.",
        "",
    ]
    io.open(OUT, "w", encoding="utf-8", newline="\n").write(
        "\n".join(head + lines) + "\n")
    print("\nwrote %s" % OUT)
    if not s1:
        raise SystemExit(1)
    return 0


if __name__ == "__main__":
    sys.exit(main())
