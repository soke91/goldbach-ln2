# -*- coding: utf-8 -*-
r"""
The family the deficit belongs to, if the norms are square-root times logs

WHAT IS AT STAKE

rem:deficitregion closed the sign axis' computational branch on a
specific pathology: fitting the deficit with polynomials in x = log N,
"the degrees say this field will resolve a new coefficient in whatever
family is offered without the residual conceding that the shape was
found" -- degree eight at t = +7.41 while the r.m.s. residual moved
only 0.032663 -> 0.030092.  It named the cure: "a derivation that says
which family the deficit belongs to."

**A quantity that is a power of log N is exactly what does that to a
polynomial fit.**  log x is not a polynomial in x, so a polynomial
family chases it forever, resolving coefficients and buying nothing.
That is the signature rem:deficitregion measured.

And there is a reason to expect logs here.  Writing the second moment
of H over k with its diagonal,

    l2^2 = sum_k (log k)^2 sum_{m,m'} mu(m)mu(m') L(N-mk) L(N-m'k),

the diagonal m = m' is sum_k (log k)^2 sum_m mu^2(m) L(N-mk)^2, the
same shape prop:V evaluates exactly in the C(N) branch.  If that term
carries l2, then l2 is of order sqrt(N) times a power of log N: the
inner sum is of order (N/k) log N, and sum over squarefree k coprime
to N of (log k)^2/k is of order (theta x)^3/3 with the same constant
(6/pi^2) prod p/(p+1) that rem:targetderived measured to 0.370 per
cent in #k.  **The exponent 0.583897 that six remarks quote for l2
would then be a logarithm absorbed into a power.**

If |sum a| is sqrt(N) times a power of log N as well, the sqrt(N)
cancels in the ratio and

    y(x) = log(|sum a| / l2) = const + C log x

exactly -- **two parameters, no x-term at all** -- and the deficit,
which is y'(x), is C/x rather than a constant plus polynomial drift.
That form makes a check available before any fitting: it forces
y''/y' = -1/x.  rem:deficitdirect published -0.007380 and 0.134019,
whose ratio is -0.05507, giving x = 18.16 -- inside the field's
[12.2061, 22.7030].  This run asks whether that is the shape or the
coincidence it could be.

Nothing is measured here.  All 156 points are read from the POINT
markers of results/audit_deficit_direct.txt and
results/audit_alpha_reach.txt, as rem:deficitshape read them.

BACKS: Remark {#rem:deficitlog} in paper/wall_v3.md.

PRE-REGISTERED PREDICTIONS (written before this script was run)

  F1  THE GATE.  The read points reproduce rem:deficitdirect: a line
      on the 81 published points returns the whole-field deficit
      0.134019, and the quadratic on all 156 returns the beta of that
      file's BETA marker, to the decimals printed there.
  F2  **The derived form fits.**  y = log(|sum a|/l2) fitted as
      c + C log x, two parameters, has r.m.s. residual no worse than
      the cubic's 0.031762, which has four.
  F3  Its implied deficit matches: C divided by the field's mean x
      reproduces 0.134019 within 10 per cent.
  F4  **And it does not keep resolving.**  Adding a linear x term to
      the log family gives |t| < 3, against the polynomial family's
      degree-eight t = +7.41.  A family that contains the shape stops
      buying coefficients; that is the whole difference.
  F5  Out of sample, walking forward from the fortieth point, the log
      family's departure is smaller than the cubic's 0.028534.

REFUTATION RULE (fixed before the run)

  F1  REFUTED outside the printed decimals; nothing below is reported.
  F2  REFUTED if the two-parameter log form's residual is the larger.
      Then the logs do not describe y and the derivation above is
      wrong about these norms.
  F3  REFUTED outside 10 per cent.
  F4  REFUTED at |t| >= 3.  Then the log family buys coefficients too
      and is no better than the polynomials -- the pathology is not
      about logarithms and rem:deficitregion's verdict stands
      untouched.
  F5  REFUTED if the departure is the larger.

  **THE UNRESOLVED CASE, NAMED.**  Over this field x runs 12.2061 to
  22.7030 and log x runs 2.5019 to 3.1224, so log x is nearly linear
  in x here and the two families are nearly the same family.  G69
  exists for this.  **If CORR of the two regressors is 0.99 or above,
  F2, F4 and F5 cannot separate the families and none of them may be
  read as a win for either**: the verdict to print is that the field
  cannot tell them apart, which is a statement about reach and not
  about shape -- the same thing rem:curvereach found on the level
  axis.  This is the outcome to expect and it must not be dressed as
  anything else.

  WHAT THIS CANNOT DO.  It does not measure the diagonal, so the
  derivation that motivates the log form is unverified here and is
  offered as motivation only.  Nothing here says where the deficit
  reaches zero, and a log family that fits inside the field licenses
  no more extrapolation than a polynomial that fits inside it --
  rem:shapepower is not repealed by finding a better-motivated shape.
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
OUT = os.path.join(RES, "audit_deficit_log.txt")

SRCS = (("audit_deficit_direct.txt", "deficitdirect"),
        ("audit_alpha_reach.txt", "alphareach"))
NPUB = 81
START = 40
DEC = 4
TOLPC = 10.0
TCAP = 3.0
CORRCAP = 0.99
CUBICRMS = 0.031762
CUBICOOS = 0.028534
DEG8T = 7.41


def read_points():
    out = {}
    for fn, tag in SRCS:
        src = io.open(os.path.join(RES, fn), encoding="utf-8").read()
        for m in re.finditer(r"^POINT %s_(\d+) ([\d.eE+-]+) "
                             r"([\d.eE+-]+)\s*$" % tag, src, re.M):
            out[int(m.group(1))] = (float(m.group(2)),
                                    float(m.group(3)))
    return out


def read_beta():
    src = io.open(os.path.join(RES, "audit_deficit_direct.txt"),
                  encoding="utf-8").read()
    m = re.search(r"^BETA deficit_direct ([-+][\d.]+) ([\d.]+)\s*$",
                  src, re.M)
    if not m:
        raise SystemExit("no BETA marker in audit_deficit_direct.txt")
    return float(m.group(1)), float(m.group(2))


def lsq(cols, y):
    """least squares on the given design columns; coefs, se, rms"""
    A = np.column_stack(cols)
    coef, *_ = np.linalg.lstsq(A, y, rcond=None)
    r = y - A @ coef
    n, p = A.shape
    dof = max(n - p, 1)
    s2 = float((r ** 2).sum()) / dof
    cov = s2 * np.linalg.pinv(A.T @ A)
    se = np.sqrt(np.abs(np.diag(cov)))
    return coef, se, math.sqrt(float((r ** 2).mean()))


def walkfwd(cols_of, y, start):
    """mean |departure| predicting each point from the ones before"""
    d = []
    for i in range(start, len(y)):
        cols = [c[:i] for c in cols_of]
        A = np.column_stack(cols)
        coef, *_ = np.linalg.lstsq(A, y[:i], rcond=None)
        row = np.array([c[i] for c in cols_of])
        d.append(abs(float(row @ coef) - float(y[i])))
    return float(np.mean(d))


HEAD = [
    "STATISTIC: y = log(|sum a| / l2) over the sign axis' field, fitted",
    "           by the derived two-parameter form c + C log x against",
    "           the polynomial family of rem:deficitshape; the",
    "           regressors' correlation, the log family's next",
    "           coefficient, and walk-forward departures.",
    "FIELD: the %d N of the sign axis, x = log N in [12.2061, 22.7030]."
    % (NPUB + 75),
    "       Nothing is measured here: every |sum a| and l2 is READ from",
    "       a POINT marker, %d from results/audit_deficit_direct.txt"
    % NPUB,
    "       and the rest from results/audit_alpha_reach.txt, and the",
    "       quadratic's beta is READ from the BETA marker of the first.",
    "",
]


def main():
    lines = []

    def say(t=""):
        print(t)
        sys.stdout.flush()
        lines.append(t)

    pts = read_points()
    pbeta, pse = read_beta()
    ns = sorted(pts)
    say("READ audit_deficit_direct.txt BETA %.6f" % pbeta)
    say("READ audit_deficit_direct.txt BETASE %.6f" % pse)
    say("  the quadratic's drift and its error, this run's gate")
    say("PRINTBOUND audit_deficit_log %d %.8f"
        % (DEC, 0.5 * 10.0 ** (-DEC)))
    say("  read %d POINT markers, %d of them published below 1.024e8"
        % (len(ns), NPUB))
    say("  caps: %.1f per cent, |t| %.1f, CORR %.2f"
        % (TOLPC, TCAP, CORRCAP))
    say("  compared against the cubic's %.6f r.m.s. and %.6f "
        "out of sample," % (CUBICRMS, CUBICOOS))
    say("  and against the polynomial family's degree-eight t = %.2f"
        % DEG8T)

    x = np.array([math.log(n) for n in ns])
    y = np.array([math.log(pts[n][0] / pts[n][1]) for n in ns])
    lx = np.log(x)
    one = np.ones_like(x)
    say("SCALES %d" % len(ns))
    say("  x runs %.4f to %.4f, log x runs %.4f to %.4f"
        % (x.min(), x.max(), lx.min(), lx.max()))

    # -------------------------------------------------------------- F1
    say()
    say("F1  the gate: do the read points reproduce rem:deficitdirect?")
    xp, yp = x[:NPUB], y[:NPUB]
    c1, _, _ = lsq([np.ones_like(xp), xp], yp)
    cq, seq, _ = lsq([one, x, 0.5 * x * x], y)
    g1a = abs(round(float(c1[1]), DEC) - 0.1340) < 10.0 ** (-DEC)
    g1b = abs(float(cq[2]) - pbeta) < 10.0 ** (-6)
    wfd = float(c1[1])
    say("  whole-field deficit on the %d published: %.6f" % (NPUB, wfd))
    say("  quadratic beta on all %d: %.6f against its %.6f"
        % (len(ns), cq[2], pbeta))
    f1 = g1a and g1b
    say("  F1 %s   (cap: %d decimals)"
        % ("hold" if f1 else "REFUTED", DEC))
    if not f1:
        io.open(OUT, "w", encoding="utf-8", newline="\n").write(
            "\n".join(HEAD + lines) + "\n")
        raise SystemExit(1)

    # ------------------------------------------------------- the CORR
    say()
    say("the separability question, before any verdict is read")
    corr = float(np.corrcoef(x, lx)[0, 1])
    say("CORR deficitlog_regressors %.6f" % corr)
    sep = corr < CORRCAP
    say("  corr(x, log x) over this field is %.6f" % corr)
    if not sep:
        say("COEFF NOT SEPARABLE deficitlog")
        say("  at or above %.2f, so F2, F4 and F5 cannot separate the"
            % CORRCAP)
        say("  families and none is read as a win for either, as the "
            "rule says")

    # -------------------------------------------------------------- F2
    say()
    say("F2  does the derived two-parameter form fit?")
    cl, sel, rmsl = lsq([one, lx], y)
    cline, _, rmsline = lsq([one, x], y)
    say("  log form   c %+.6f  C %+.6f (s.e. %.6f)  r.m.s. %.6f"
        % (cl[0], cl[1], sel[1], rmsl))
    say("  line       r.m.s. %.6f" % rmsline)
    say("  cubic (4 parameters, read)  r.m.s. %.6f" % CUBICRMS)
    f2 = rmsl <= CUBICRMS
    say("POINT logrms %.6f" % rmsl)
    say("  F2 %s   (cap: %.6f)"
        % ("hold" if f2 else "REFUTED", CUBICRMS))

    # -------------------------------------------------------------- F3
    say()
    say("F3  does its implied deficit match the published one?")
    xm = float(x.mean())
    imp = float(cl[1]) / xm
    pc = 100.0 * (imp - wfd) / wfd
    f3 = abs(pc) <= TOLPC
    say("  C / mean x = %.6f / %.4f = %.6f against F1's own %.6f "
        "(%+.2f per cent)" % (cl[1], xm, imp, wfd, pc))
    say("POINT impliedeficit %.6f" % imp)
    say("  F3 %s   (cap: %.1f per cent)"
        % ("hold" if f3 else "REFUTED", TOLPC))

    # -------------------------------------------------------------- F4
    say()
    say("F4  does the log family keep buying coefficients?")
    c3, se3, rms3 = lsq([one, lx, x], y)
    t3 = float(c3[2]) / float(se3[2]) if se3[2] else float("inf")
    f4 = abs(t3) < TCAP
    say("  adding an x term: %+.8f +- %.8f, t = %+.2f"
        % (c3[2], se3[2], t3))
    say("  r.m.s. %.6f against the two-parameter %.6f" % (rms3, rmsl))
    say("TSTAT deficitlog_next %.2f" % t3)
    say("SPREAD deficitlog_next %.8f" % float(se3[2]))
    say("  the polynomial family's degree eight read %.2f" % DEG8T)
    say("  F4 %s   (cap: |t| < %.1f)"
        % ("hold" if f4 else "REFUTED", TCAP))

    # -------------------------------------------------------------- F5
    say()
    say("F5  which family predicts forward better?")
    dl = walkfwd([one, lx], y, START)
    dc = walkfwd([one, x, x * x, x ** 3], y, START)
    f5 = dl < CUBICOOS
    say("  from point %d: log form %.6f, cubic here %.6f, cubic read "
        "%.6f" % (START, dl, dc, CUBICOOS))
    say("POINT logoos %.6f" % dl)
    say("  F5 %s   (cap: %.6f)"
        % ("hold" if f5 else "REFUTED", CUBICOOS))

    say()
    say("=" * 70)
    say("F1 %s  F2 %s  F3 %s  F4 %s  F5 %s"
        % tuple("hold" if v else "REFUTED"
                for v in (f1, f2, f3, f4, f5)))
    say()
    if not sep:
        say("the field cannot tell the two families apart. corr(x, "
            "log x) is %.6f" % corr)
        say("over this range, so whatever F2, F4 and F5 read, they "
            "read it about")
        say("reach and not about shape. That is the rule this run "
            "registered and")
        say("it is the same limit rem:curvereach found on the level "
            "axis: the")
        say("discrimination is not in the data, and no fit performed "
            "on it can")
        say("put it there.")
    elif f2 and f4:
        say("the deficit's family is the logarithmic one. two "
            "parameters do what")
        say("four did, and the family stops buying coefficients where "
            "the")
        say("polynomials never did -- which is what rem:deficitregion "
            "said would")
        say("distinguish a family that contains the shape from one "
            "that chases it.")
        say("this licenses no extrapolation. rem:shapepower is not "
            "repealed by a")
        say("better-motivated shape, and nothing here says where the "
            "deficit")
        say("reaches zero.")
    else:
        say("the logarithmic family is not it either. the pathology "
            "rem:deficitregion")
        say("measured is not explained by a power of log N in the "
            "norms, and its")
        say("verdict stands where it stood.")

    io.open(OUT, "w", encoding="utf-8", newline="\n").write(
        "\n".join(HEAD + lines) + "\n")
    print("\nwrote %s" % os.path.normpath(OUT))
    return 0


if __name__ == "__main__":
    sys.exit(main())
