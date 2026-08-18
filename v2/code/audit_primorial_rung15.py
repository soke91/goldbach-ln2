# -*- coding: utf-8 -*-
r"""
The sixteenth rung: the curvature's second out-of-sample test.

WHAT IS AT STAKE

{#rem:curvebound} found the quadratic predicting rungs it was not
fitted to -- +0.0046 and +0.0023 inside its own errors, where the line
gave +0.0098 and +0.0097 outside them -- and put theta' a decade
earlier than the line at log10 N = 9.6068. But its drift, 0.5688 from
refitting on the lower rungs, is as wide as its bracket, so the
extrapolation is not stable and one more point is the only thing that
narrows it.

30030 * 2^15 = 984023040 at log10 N = 8.9930 is that point. It asks
the quadratic the same question a second time, and it says whether
adding a rung settles the crossing or moves it again.

BACKS: Remark {#rem:rung15} in paper/theorem_A.md.

PRE-REGISTERED PREDICTIONS (written before this script was run)

  H1  The control. Rung 14 recomputed here reproduces
      results/audit_primorial_rung14.txt inside the bound its
      printing forces.
  H2  The margin keeps growing: the new exponent's margin over 1/2
      exceeds rung 14's 0.0333.
  H3  The curvature predicts again: the quadratic fitted on the
      fifteen published rungs puts this one inside its own prediction
      standard error.
  H4  And beats the line again: its departure is smaller than the
      line's on the same fifteen.
  H5  The extrapolation settles: refitting with this rung moves the
      0.56 crossing by less than the 0.5688 drift
      {#rem:curvebound} declared.

REFUTATION RULE (fixed before the run)

  H1  REFUTED outside the printing bound. THIS ONE GATES.
  H2  REFUTED if the margin does not grow. Six rungs have grown in a
      row; a seventh that does not would end it.
  H3  REFUTED if the rung falls outside. Then the quadratic's two
      successes were the two points nearest the fit and it does not
      predict at reach -- which would leave item 1 with no shape that
      predicts, not with a worse one.
  H4  REFUTED if the line does better. Same reading as H3 by the
      comparison that matters.
  H5  REFUTED if the crossing moves further than the declared drift.
      Then each new rung relocates theta' by as much as the last, the
      extrapolation is not converging, and no crossing may be quoted
      at any reach -- which is a stronger statement than
      {#rem:shapepower}'s and would close item 1's forecast for good.

  H1 gates. H2 to H5 are the measurement and do not gate. "Inside"
  in H3 includes "too noisy to tell" and the ratio printed is what
  says which (M9).

  NO NULL IS RUN and none applies. A deterministic curve is located
  against a computed threshold. The coin arms for this statistic were
  run in lab_primorial_ladder.py and lab_primorial_share.py.
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
OUT = os.path.join(RES, "audit_primorial_rung15.txt")

CLIM = 4_000_000                    # the fixed Euler bound (G20)
SEED = 20260823
DRAWS = 4000
TARGET = 0.56


def module(name):
    p = os.path.join(CODE, name + ".py")
    spec = importlib.util.spec_from_file_location(name, p)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


R11 = module("audit_primorial_rung11")
primes_upto = R11.primes_upto
BASE = R11.BASE
CONTROL = BASE * (1 << 14)          # 492011520, the rung 14 point
NEW = BASE * (1 << 15)              # 984023040


def read_all():
    """the fifteen rung exponents, rung 14's margin, and the crossing"""
    src = io.open(os.path.join(RES, "audit_primorial_rung10.txt"),
                  encoding="utf-8").read()
    i = src.index("N            log10 N   exponent   fitted     "
                  "residual")
    ns, ex, dec = [], [], 0
    for ln in src[i:].splitlines()[1:]:
        f = ln.split()
        if len(f) < 3 or not f[0].isdigit():
            break
        ns.append(int(f[0]))
        ex.append(float(f[2]))
        dec = max(dec, len(f[2].split(".")[1]))
    for j in (11, 12, 13, 14):
        s = io.open(os.path.join(RES, "audit_primorial_rung%d.txt" % j),
                    encoding="utf-8").read()
        N = BASE * (1 << j)
        m = re.search(r"^  N = " + str(N) +
                      r"\s+thr [\d.]+\s+#k \d+\s+beta [\d.]+\s+"
                      r"K\*_R \d+\s+exp ([\d.]+)\s*$", s, re.M)
        ns.append(N)
        ex.append(float(m.group(1)))
    s14 = io.open(os.path.join(RES, "audit_primorial_rung14.txt"),
                  encoding="utf-8").read()
    marg14 = float(re.search(r"the new exponent is [\d.]+, margin "
                             r"([\d.]+),", s14).group(1))
    scat = float(re.search(r"^FLOOR primorial_rung14 ([\d.]+)\s*$",
                           s14, re.M).group(1))
    sb = io.open(os.path.join(RES, "audit_curve_bound.txt"),
                 encoding="utf-8").read()
    m = re.search(r"^BRACKET ladder_quadratic_theta_prime ([\d.]+) "
                  r"([\d.]+) ([\d.]+)\s*$", sb, re.M)
    dr = float(re.search(r"^DRIFT ladder_quadratic_theta_prime "
                         r"([\d.]+)\s*$", sb, re.M).group(1))
    return (ns, ex, dec, marg14, scat, float(m.group(1)),
            float(m.group(2)), float(m.group(3)), dr)


def quadfit(x, y):
    A = np.column_stack([np.ones_like(x), x, x * x])
    c, *_ = np.linalg.lstsq(A, y, rcond=None)
    r = y - A.dot(c)
    n = x.size
    s2 = float((r ** 2).sum()) / (n - 3)
    cov = s2 * np.linalg.inv(A.T.dot(A))
    return c, cov, s2, float(np.sqrt((r ** 2).mean()))


def cross(c, level):
    a2, b2, c2 = c[0] - level, c[1], c[2]
    if abs(c2) < 1e-18:
        return None if abs(b2) < 1e-18 else -a2 / b2
    disc = b2 * b2 - 4.0 * c2 * a2
    if disc < 0:
        return None
    rs = [r for r in ((-b2 + math.sqrt(disc)) / (2.0 * c2),
                      (-b2 - math.sqrt(disc)) / (2.0 * c2)) if r > 0]
    return min(rs) if rs else None


def main():
    lines = []

    def say(s=""):
        print(s)
        lines.append(s)

    (ns, ex, dec, marg14, scat, pt56, lo56, hi56,
     drift56) = read_all()
    say("read %d rung exponents, rung 14's margin %.4f, the scatter "
        "%.4f," % (len(ns), marg14, scat))
    say("  and the quadratic's 0.56 crossing %.4f [%.4f, %.4f] with "
        "drift %.4f" % (pt56, lo56, hi56, drift56))
    say("  from results/audit_primorial_rung10 through rung14 and "
        "results/audit_curve_bound.txt")
    say("SEED %d" % SEED)
    say("DRAWS %d" % DRAWS)
    say()
    say("the rungs: the control %d at log10 N = %.4f and the new %d "
        "at %.4f" % (CONTROL, math.log10(CONTROL), NEW,
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
    e14, e15 = got[CONTROL][1], got[NEW][1]

    # -------------------------------------------------------------- H1
    say()
    say("H1  the control at rung 14")
    rnd = 0.5 * 10.0 ** (-dec)
    d = abs(e14 - ex[-1])
    h1 = d <= rnd
    say("  exponent here %.4f against the published %.4f, departure "
        "%.6f; the bound from %d decimals is %.8f"
        % (e14, ex[-1], d, dec, rnd))
    say("PRINTBOUND audit_primorial_rung15 %d %.8f" % (dec, rnd))
    say("  H1 %s   (cap: the printing bound)"
        % ("hold" if h1 else "REFUTED"))

    # -------------------------------------------------------------- H2
    say()
    say("H2  does the margin keep growing?")
    marg = e15 - 0.5
    h2 = marg > marg14
    say("  the new exponent is %.4f, margin %.4f, against rung 14's "
        "%.4f and the scatter %.4f" % (e15, marg, marg14, scat))
    say("MARGIN audit_primorial_rung15 %.4f %.4f" % (marg, scat))
    if marg <= scat:
        say("INSIDE FLOOR audit_primorial_rung15")
    say("FLOOR primorial_rung15 %.4f" % scat)
    say("  H2 %s   (cap: rung 14's margin)"
        % ("hold" if h2 else "REFUTED"))

    # -------------------------------------------------------- H3, H4
    x = np.log(np.array(ns, dtype=np.float64))
    y = np.array(ex)
    c, cov, s2, rms = quadfit(x, y)
    a, b = np.polyfit(x, y, 1)
    xn = math.log(NEW)
    v = np.array([1.0, xn, xn * xn])
    pq = float(v.dot(c))
    sp = math.sqrt(s2 + float(v.dot(cov).dot(v)))
    pl = a * xn + b
    dq, dl = e15 - pq, e15 - pl
    h3 = abs(dq) <= sp
    h4 = abs(dq) < abs(dl)
    say()
    say("H3/H4  does the curvature predict this one too?")
    say("  fitted on the %d published rungs:" % len(ns))
    say("  shape        predicts   measured   departure   pred s.e.  "
        " ratio")
    say("  quadratic    %-10.4f %-10.4f %+-11.4f %-11.4f %.2f"
        % (pq, e15, dq, sp, abs(dq) / sp))
    say("  line         %-10.4f %-10.4f %+-11.4f" % (pl, e15, dl))
    say("  H3 %s   (cap: the prediction standard error)"
        % ("hold" if h3 else "REFUTED"))
    say("  H4 %s   (cap: the line's departure)"
        % ("hold" if h4 else "REFUTED"))

    # -------------------------------------------------------------- H5
    say()
    say("H5  does the crossing settle?")
    x16 = np.append(x, xn)
    y16 = np.append(y, e15)
    c16, cov16, s216, rms16 = quadfit(x16, y16)
    p16 = cross(c16, TARGET)
    moved = abs(p16 / math.log(10.0) - pt56)
    h5 = moved < drift56
    rng = np.random.default_rng(SEED)
    draws = rng.multivariate_normal(c16, cov16, size=DRAWS)
    vals = [cross(dd, TARGET) for dd in draws]
    vals = [w / math.log(10.0) for w in vals
            if w is not None and w > x16.max()]
    lo = float(np.percentile(vals, 2.5))
    hi = float(np.percentile(vals, 97.5))
    say("  the sixteen-rung quadratic is %+.8f in (log N)^2, r.m.s. "
        "%.4f" % (c16[2], rms16))
    say("  it reaches 0.56 at log10 N = %.4f, bracket [%.4f, %.4f] "
        "from %d of %d draws"
        % (p16 / math.log(10.0), lo, hi, len(vals), DRAWS))
    say("BRACKET ladder_quadratic16_theta_prime %.4f %.4f %.4f"
        % (p16 / math.log(10.0), lo, hi))
    say("DRIFT ladder_quadratic16_theta_prime %.4f" % moved)
    say("  the fifteen-rung value was %.4f, so it moved %.4f against "
        "the declared drift %.4f" % (pt56, moved, drift56))
    say("SHAPES 2")
    say("SCATTER slope_audit_primorial_rung15 %.4f" % rms16)
    say("  H5 %s   (cap: the declared drift)"
        % ("hold" if h5 else "REFUTED"))
    say("  no forecast is made from this; {#rem:shapepower} is why.")

    say()
    say("=" * 70)
    say("H1 %s  H2 %s  H3 %s  H4 %s  H5 %s"
        % tuple("hold" if v_ else "REFUTED"
                for v_ in (h1, h2, h3, h4, h5)))

    head = [
        "STATISTIC: the truncation K*_R at which",
        "           sum_{k<K}(log k)|R(N;k)| first reaches",
        "           S(N)(1-A(N))N, and its exponent log K*_R / log N,",
        "           at N = 30030*2^15 = 984023040 and, as a control,",
        "           at N = 30030*2^14; the margin over 1/2 against",
        "           rung 14's; the quadratic and the line fitted on",
        "           the fifteen published rungs and asked for this",
        "           one, with the prediction standard error at its",
        "           abscissa; and the sixteen-rung quadratic's 0.56",
        "           crossing with a bracket from its own parameter",
        "           covariance, against the fifteen-rung value.",
        "NULL: none is run and none applies. A deterministic curve is",
        "      located against a computed threshold; there is no",
        "      background to detect against. The coin arms for this",
        "      statistic were run in lab_primorial_ladder.py and",
        "      lab_primorial_share.py. The bracket is drawn from the",
        "      fit's own parameter covariance with the fixed SEED.",
        "FIELD: N = 492011520 and 984023040, the odd radical",
        "       3*5*7*11*13 fixed so the threshold is constant; k",
        "       squarefree and coprime to N with 2 <= k < 100000; m",
        "       odd, squarefree and coprime to k, m < N/k; the sieve",
        "       weight over the odd primes below 30; the Euler",
        "       products at the fixed bound 4000000. One odd radical,",
        "       as the RADICALS line declares. The statistic, the",
        "       sieve and the k-cap are imported from",
        "       code/audit_primorial_rung11.py; the fifteen published",
        "       rungs come from results/audit_primorial_rung10.txt",
        "       and rung11 through rung14, and the published crossing",
        "       from results/audit_curve_bound.txt.",
        "",
    ]
    io.open(OUT, "w", encoding="utf-8", newline="\n").write(
        "\n".join(head + lines) + "\n")
    print("\nwrote %s" % OUT)
    if not h1:
        raise SystemExit(1)
    return 0


if __name__ == "__main__":
    sys.exit(main())
