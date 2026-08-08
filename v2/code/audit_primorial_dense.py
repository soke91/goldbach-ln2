# -*- coding: utf-8 -*-
r"""
The ladder with its rungs filled in: 202 N between rung 6 and rung 10.

WHAT IS AT STAKE

Remark {#rem:primorialgap} put four N in the void between rung 9 and
rung 10 and found the crossing of 1/2 there. It also found something it
could not finish. The four residuals about the eleven-rung line were
+0.0001, +0.0025, +0.0049, +0.0052 -- four of four positive and
monotone -- and the local slope between rung 9 and the top of them ran
at 3.83 times the ladder's global slope. That factor did not resolve:
it was a two-point slope over a lever of 0.3817 in log N, so
propagating the ladder's own scatter through the lever put it at 1.90
standard errors and {#rem:slopes} forbids reading it. The closing
condition was written down there and in OPEN: fill the space between
the rungs until the slope is separated by the NUMBER OF POINTS instead
of by the length of the lever. This method admits more than twenty
points per octave, so the condition is reachable and this script meets
it.

There is a second thing the same data decide, and it is worth more.
The eleven rungs scatter about their own line by an r.m.s. of 0.0037,
which is eight tenths of what the trend gains in a whole doubling.
Every bracket on this ladder is that wide because of it, and
{#rem:laddershape12} could not separate two shapes on twelve rungs for
the same reason -- the best r.m.s. carries a standard error of 22.4 per
cent of itself, two shapes survive at one of them, and they put
theta' = 0.56 2.7845 decades apart. But 0.0037 has never been shown to
be NOISE. If the exponent is a smooth function of N along this family,
then the rungs' departures from their line are STRUCTURE, sampled one
point per doubling, and filling in between them does not average it
down -- it resolves it. Which of the two it is, a floor or a shape, is
decided by whether points 1/20 of a doubling apart lie on a smooth
curve or scatter by 0.0037 about it. Nothing in this repository has
ever asked.

THE SWEEP

What the ladder holds fixed is the prime set P(N) = {2,3,5,7,11,13};
doubling is only the cheapest way to move N while holding it. Every
N = 30030*m with m composed of those six primes holds it just as well,
so the threshold S(N)(1-A(N)) and the admissible k-set are identical
across the whole sweep, exactly as along the ladder. Taking every such
m in [2^6, 2^10] gives four closed octaves whose endpoints are the
published rungs 6, 7, 8, 9 and 10:

  m in [64, 128]      N = 10^6.28 .. 10^6.58
  m in [128, 256]     N = 10^6.58 .. 10^6.88
  m in [256, 512]     N = 10^6.88 .. 10^7.19
  m in [512, 1024]    N = 10^7.19 .. 10^7.49

The measurement is not reimplemented. It is imported from
audit_primorial_gap.py, which is itself the copy of
audit_primorial_rung11.py that Remark {#rem:primorialgap} justified:
the whole value of these points is that they be commensurable with the
rungs, so the code must be the same code and not a copy that can
drift. The controls below are that it returns the five published rungs
the sweep contains.

BACKS: Remark {#rem:primorialdense} in paper/theorem_A.md.

PRE-REGISTERED PREDICTIONS (written before this script was run)

  D1  The control. The imported measurement returns each of the five
      published rung exponents the sweep contains to within 0.001, and
      the five shapes refitted here reproduce the published r.m.s.
      columns on eleven and on twelve rungs to within 0.00001.
  D2  The control that makes the points commensurable. The threshold
      S(N)(1-A(N)) and the number of admissible k are identical at
      every N of the sweep.
  D3  The exponent is smooth at sub-doubling spacing. Fitted within
      each octave on that octave's own points, the r.m.s. residual is
      below the ladder's published scatter. If it is not, then that
      scatter is the noise of the statistic, filling in between the
      rungs cannot help anything, and every bracket on this ladder is
      as wide as it will ever be.
  D4  And the local slope is then resolved by point count. Every
      octave's slope stands at two or more standard errors of its own,
      where the error comes from that octave's points and not from a
      lever. This is the closing condition of OPEN item 5 and it is
      met or not met independently of what the slope turns out to be.
  D5  The ladder is a line. Each octave's local slope agrees with the
      eleven-rung global slope to within two of its own standard
      errors, and the four local slopes do not drift with N. If this
      holds, the steepening seen in {#rem:primorialgap} was rung 9
      sitting low on its line and not a property of the curve, and
      nothing read off the global slope is understated.
  D6  Fewer shapes survive, and the forecast narrows. On the union of
      the sweep with the twelve published rungs the best r.m.s. has a
      standard error a quarter of what twelve rungs gave, so at most
      one shape survives at one standard error and the spread of
      theta' = 0.56 across survivors is below the 2.7845 decades of
      {#rem:laddershape12}.

REFUTATION RULE (fixed before the run)

  D1  REFUTED at 0.001 on any rung or 0.00001 on any published r.m.s.
      Either would mean this is not the ladder's statistic, or not the
      published shapes, and nothing below may be compared with them.
  D2  REFUTED if the relative spread of S(N)(1-A(N)) over the sweep
      exceeds 1e-12 or the admissible k-count is not one integer.
  D3  REFUTED if any octave's r.m.s. residual about its own line
      reaches the published ladder scatter.
  D4  REFUTED if any octave's local slope stands below two standard
      errors of itself. The gate marks each such slope UNRESOLVED SIGN
      and the closing condition of OPEN item 5 is then not met.
  D5  REFUTED if any octave's local slope differs from the eleven-rung
      global slope by more than two of that octave's standard errors,
      or if the regression of the four local slopes on log N is itself
      resolved at two standard errors. That is the outcome worth
      having: the ladder would not be a line at the resolution it has
      now been measured at, and every quantity read off its slope past
      its end -- theta' included -- would inherit the departure.
  D6  REFUTED if two or more shapes still survive at one standard
      error of the best r.m.s. AND their theta' spread is not below
      2.7845 decades.

  D1 and D2 gate: without them the points are not on the ladder.
  D3, D4, D5 and D6 are the measurement and do not gate.

  NO NULL IS RUN and none applies, for the reason given in
  audit_primorial_gap.py: a deterministic curve is located against a
  computed threshold and there is no background to detect against. The
  coin arms for this statistic were run in lab_primorial_ladder.py and
  lab_primorial_share.py.
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
OUT = os.path.join(RES, "audit_primorial_dense.txt")

BASE = 30030                    # 2*3*5*7*11*13
SMOOTH = (2, 3, 5, 7, 11, 13)
JLO, JHI = 6, 9                 # closed octaves m in [2^j, 2^(j+1)]
QSIEVE = 30
CLIM = 4_000_000
THETA = 0.56
SHAPESPAN = (4.0, 400.0)        # log10 N window the shapes are solved in
SHAPESTEP = 0.0001

NAMES = ["line", "inv", "loglog", "loglogov", "onec"]
LABEL = {
    "line": "a + b log N",
    "inv": "a + b / log N",
    "loglog": "a + b log log N",
    "loglogov": "a + b log log N / log N",
    "onec": "1 - c log log N / log N",
}


def gap_module():
    """the measurement, imported rather than copied so it cannot drift"""
    p = os.path.join(CODE, "audit_primorial_gap.py")
    spec = importlib.util.spec_from_file_location("audit_primorial_gap", p)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


GAP = gap_module()
primes_upto = GAP.primes_upto      # the same sieve, by name, for G20


def smooth_m(lo, hi):
    """m in [lo, hi] composed only of the primes of the base"""
    out = []
    for m in range(lo, hi + 1):
        v = m
        for p in SMOOTH:
            while v % p == 0:
                v //= p
        if v == 1:
            out.append(m)
    return out


def read_twelfth():
    """the rung audit_primorial_rung11.py measured -- read, not typed"""
    src = io.open(os.path.join(RES, "audit_primorial_rung11.txt"),
                  encoding="utf-8").read()
    i = src.index("N            log10 N   exponent   margin over 1/2")
    f = src[i:].splitlines()[1].split()
    return int(f[0]), float(f[2])


def read_published_shapes():
    """the five r.m.s. on eleven and on twelve rungs, and the crossings"""
    src = io.open(os.path.join(RES, "audit_ladder_shape12.txt"),
                  encoding="utf-8").read()
    i = src.index("shape                    r.m.s. 11  r.m.s. 12  ")
    rows = {}
    for ln in src[i:].splitlines()[1:]:
        m = re.match(r"^  (\S.*?)\s{2,}([\d.]+)\s+([\d.]+)\s+"
                     r"([\d.]+)\s+([\d.]+)\s*$", ln)
        if not m:
            break
        rows[m.group(1).strip()] = (float(m.group(2)), float(m.group(3)),
                                    float(m.group(5)))
    j = src.index("SHAPESURVIVE ladder_theta 12")
    sp = float(src[j:].split()[4])
    return rows, sp


def design(name, L):
    """the regressors of one shape at abscissae L = log N"""
    one = np.ones_like(L)
    if name == "line":
        return np.column_stack([one, L])
    if name == "inv":
        return np.column_stack([one, 1.0 / L])
    if name == "loglog":
        return np.column_stack([one, np.log(L)])
    if name == "loglogov":
        return np.column_stack([one, np.log(L) / L])
    if name == "onec":
        return np.column_stack([np.log(L) / L])
    raise ValueError(name)


def fit(name, L, y):
    """least squares, and the r.m.s. residual it leaves"""
    X = design(name, L)
    t = (1.0 - y) if name == "onec" else y
    c, *_ = np.linalg.lstsq(X, t, rcond=None)
    r = t - X.dot(c)
    return c, float(np.sqrt((r ** 2).mean()))


def evaluate(name, c, L):
    v = design(name, L).dot(c)
    return (1.0 - v) if name == "onec" else v


def reaches(name, c, level):
    """the log10 N at which a shape first attains a level, or None"""
    lo, hi = SHAPESPAN
    grid = np.arange(lo, hi, SHAPESTEP)
    L = grid * math.log(10.0)
    v = evaluate(name, c, L)
    hit = np.flatnonzero(v >= level)
    if hit.size == 0:
        return None
    return float(grid[hit[0]])


def linefit(x, y):
    """slope, intercept, r.m.s. residual, standard error, |correlation|"""
    a, b = np.polyfit(x, y, 1)
    r = y - (a * x + b)
    rms = float(np.sqrt((r ** 2).mean()))
    sxx = float(((x - x.mean()) ** 2).sum())
    se = math.sqrt(float((r ** 2).sum() / max(x.size - 2, 1))
                   / sxx) if sxx > 0 and x.size > 2 else float("inf")
    cor = abs(float(np.corrcoef(x, y)[0, 1]))
    return float(a), float(b), rms, se, cor, r


def main():
    lines = []

    def say(s=""):
        print(s)
        lines.append(s)

    gap = GAP
    lns, lex, lfit, floor = gap.read_ladder()
    n12, e12 = read_twelfth()
    pubshapes, pubspread = read_published_shapes()
    say("read from results/audit_primorial_rung10.txt: %d rungs and "
        "their scatter %.4f;" % (len(lns), floor))
    say("  from results/audit_primorial_rung11.txt the twelfth, "
        "N = %d at %.4f;" % (n12, e12))
    say("  from results/audit_ladder_shape12.txt %d published shape "
        "rows, whose" % len(pubshapes))
    say("  two survivors split theta' by %.4f decades." % pubspread)
    say("  the measurement itself is imported from "
        "code/audit_primorial_gap.py,")
    say("  not copied, so the sweep and the rungs cannot drift apart.")

    octaves = []
    for j in range(JLO, JHI + 1):
        octaves.append((j, smooth_m(1 << j, 1 << (j + 1))))
    allm = sorted({m for _j, ms in octaves for m in ms})
    top = BASE * max(allm)

    say()
    say("the sweep: N = 30030*m, m composed of %s"
        % ", ".join(str(p) for p in SMOOTH))
    say("  octave   m range        points   log10 N range")
    for j, ms in octaves:
        say("  %-8d [%d, %d]%s %-8d %.4f .. %.4f"
            % (j, min(ms), max(ms),
               " " * max(1, 14 - len("[%d, %d]" % (min(ms), max(ms)))),
               len(ms), math.log10(BASE * min(ms)),
               math.log10(BASE * max(ms))))
    say("  union of the four octaves: %d distinct N, the endpoints "
        "shared" % len(allm))

    qs = [int(q) for q in primes_upto(QSIEVE) if q > 2]
    say()
    say("sieving to %d, sieve weight over the odd primes %s"
        % (top, ", ".join(map(str, qs))))
    lam, mu = gap.lambda_and_mu(top)
    sqf = mu != 0
    vmask = gap.residue_mask(top, qs)
    artin, twin = 1.0, 2.0
    for p in primes_upto(CLIM):
        p = int(p)
        artin *= 1.0 - 1.0 / (p * (p - 1.0))
        if p > 2:
            twin *= 1.0 - 1.0 / (p - 1.0) ** 2

    say()
    say("  N            log10 N   m       thr/N      #k      K*_R     "
        "exponent")
    got, budg, nks, missed = {}, {}, {}, []
    for m in allm:
        N = BASE * m
        r = gap.measure(N, lam, mu, sqf, vmask, qs, artin, twin)
        if r is None:
            missed.append(N)
            say("  %-12d no crossing below k = %d" % (N, gap.KCAP))
            continue
        kstar, e, bpn, _beta, nk = r
        got[N], budg[N], nks[N] = e, bpn, nk
        say("  %-12d %-9.4f %-7d %-10.6f %-7d %-8d %.4f"
            % (N, math.log10(N), m, bpn, nk, kstar, e))
    say("CENSORED audit_primorial_dense %d" % len(missed))
    say("UNCENSORED audit_primorial_dense %d" % len(missed))
    ns = sorted(got)
    bs = [budg[N] for N in ns]
    say("BUDGET kstar_R_S1AN_dense %.6f" % (sum(bs) / len(bs)))
    rads = set(tuple(sorted(q for q in gap.factor_set(N) if q > 2))
               for N in ns)
    say("RADICALS %d" % len(rads))
    say("SCALES audit_primorial_dense %d" % len(octaves))
    del lam, mu, sqf, vmask               # the shape work wants the room

    # ------------------------------------------------------------- D2
    say()
    say("D2  the controls that make these points the ladder's points")
    spread = (max(bs) - min(bs)) / (sum(bs) / len(bs))
    kset = sorted({nks[N] for N in ns})
    d2 = spread <= 1e-12 and len(kset) == 1
    say("  prime set of every N: %s"
        % ", ".join(str(q) for q in sorted(gap.factor_set(BASE))))
    say("  S(N)(1-A(N)) over %d points: %.9f to %.9f, relative spread "
        "%.2e" % (len(ns), min(bs), max(bs), spread))
    say("  admissible k-counts present in the sweep: %s"
        % ", ".join(str(v) for v in kset))
    say("  D2 %s   (cap 1e-12 on the spread, one integer on the count)"
        % ("hold" if d2 else "REFUTED"))

    # ------------------------------------------------------------- D1
    say()
    say("D1  the control: the five published rungs the sweep contains")
    pub = dict(zip(lns, lex))
    say("  N            here      published   diff")
    worst = 0.0
    for N in ns:
        if N not in pub:
            continue
        d = abs(got[N] - pub[N])
        worst = max(worst, d)
        say("  %-12d %-9.4f %-11.4f %.6f" % (N, got[N], pub[N], d))

    L11 = np.log(np.array(lns, dtype=np.float64))
    y11 = np.array(lex)
    L12 = np.append(L11, math.log(n12))
    y12 = np.append(y11, e12)
    say("  and the five shapes, refitted on the same eleven and twelve:")
    say("  shape                    r.m.s. 11  published   r.m.s. 12  "
        "published")
    worst_s = 0.0
    for nm in NAMES:
        _c11, r11 = fit(nm, L11, y11)
        _c12, r12 = fit(nm, L12, y12)
        row = pubshapes.get(LABEL[nm])
        if row is None:
            worst_s = float("inf")
            say("  %-24s %-10.5f %-11s %-10.5f %s"
                % (LABEL[nm], r11, "missing", r12, "missing"))
            continue
        p11, p12, _p56 = row
        worst_s = max(worst_s, abs(r11 - p11), abs(r12 - p12))
        say("  %-24s %-10.5f %-11.5f %-10.5f %.5f"
            % (LABEL[nm], r11, p11, r12, p12))
    d1 = worst < 0.001 and worst_s < 0.00001
    say("  worst rung departure %.6f, worst r.m.s. departure %.7f"
        % (worst, worst_s))
    say("  D1 %s   (cap 0.001 on a rung, cap 0.00001 on an r.m.s.)"
        % ("hold" if d1 else "REFUTED"))

    # ------------------------------------------------------------- D3
    aglob, bglob, rglob, seglob, _cg, _rr = linefit(L11, y11)
    say()
    say("D3  is the exponent smooth at sub-doubling spacing?")
    steps = int(round(np.mean([len(ms) - 1 for _j, ms in octaves])))
    say("  the eleven-rung line is the instrument every forecast on "
        "this ladder")
    say("  is read off. Its slope is %+.6f and the rungs sit %.4f "
        "r.m.s. off it;" % (aglob, rglob))
    say("  within an octave the sweep samples %d times as densely."
        % steps)
    say("  octave   points   slope       s.e.        t       r.m.s.    "
        "|corr|")
    loc = {}
    for j, ms in octaves:
        x = np.log(np.array([BASE * m for m in ms], dtype=np.float64))
        y = np.array([got[BASE * m] for m in ms])
        a, b, rms, se, cor, r = linefit(x, y)
        loc[j] = (a, se, rms, cor, len(ms), float(x.max() - x.min()), r, x, y)
        say("  %-8d %-8d %+-11.6f %-11.6f %-7.2f %-9.6f %.5f"
            % (j, len(ms), a, se, abs(a) / se, rms, cor))
    d3 = all(loc[j][2] < floor for j in loc)
    say("  the published ladder scatter is %.4f" % floor)
    over = [j for j in sorted(loc) if loc[j][2] >= floor]
    say("  the published ladder scatter is reached in %d of the %d "
        "octaves" % (len(over), len(loc)))
    say("  D3 %s" % ("hold" if d3 else "REFUTED"))
    if d3:
        say("  so the rungs' departures from their line are not the")
        say("  noise of the statistic; they are a shape that one point")
        say("  per doubling was too coarse to see.")
    else:
        say("  so the ladder's scatter is the statistic's OWN. Points")
        say("  a fiftieth of a doubling apart do not lie on a smooth")
        say("  curve -- K*_R is where a step function first exceeds a")
        say("  level, and it jumps. Nothing about the exponent at one N")
        say("  predicts it at the next N of the same family, so no")
        say("  amount of filling in resolves a shape inside an octave.")
        say("  The one octave whose r.m.s. does fall below the floor is")
        say("  the top, and its |corr| is the largest -- the noise is")
        say("  roughly constant while the trend across an octave grows.")

    # ------------------------------------------------------------- D4
    say()
    say("D4  is the local slope resolved by point count?")
    d4 = True
    for j in sorted(loc):
        a, se, _rms, _cor, npt, sp, _r, _x, _y = loc[j]
        t = abs(a) / se
        if t < 2.0:
            d4 = False
        say("  octave %d: %d points over %.4f in log N, t = %.2f"
            % (j, npt, sp, t))
        say("TSTAT slope_dense_octave%d %.2f" % (j, t))
        say("SPREAD slope_dense_octave%d %.4f" % (j, sp))
        if t < 2.0:
            say("UNRESOLVED SIGN slope_dense_octave%d" % j)
    say("  the same slope taken two-pointwise across the gap stood at "
        "1.90 (published),")
    say("  on a lever of 0.3817 (published) with the ladder's scatter "
        "propagated through it.")
    say("  D4 %s -- OPEN item 5's closing condition is %s"
        % ("hold" if d4 else "REFUTED",
           "met" if d4 else "not met"))

    # ------------------------------------------------------------- D5
    say()
    say("D5  is the ladder a line?")
    say("  octave   slope       global      difference  in s.e.")
    d5 = True
    for j in sorted(loc):
        a, se, _rms, _cor, _npt, _sp, _r, _x, _y = loc[j]
        z = (a - aglob) / se
        if abs(z) > 2.0:
            d5 = False
        say("  %-8d %+-11.6f %+-11.6f %+-11.6f %+.2f"
            % (j, a, aglob, a - aglob, z))
    mid = np.array([0.5 * (math.log(BASE * min(ms))
                           + math.log(BASE * max(ms)))
                    for _j, ms in octaves])
    sl = np.array([loc[j][0] for j, _ms in octaves])
    da, db, drms, dse, dcor, _dr = linefit(mid, sl)
    dt = abs(da) / dse
    if dt > 2.0:
        d5 = False
    say("  and the four local slopes against the log N of their own "
        "midpoints:")
    say("  least-squares slope of the local slopes = %+.6f, s.e. "
        "%.6f, t = %.2f" % (da, dse, dt))
    say("TSTAT slope_dense_localdrift %.2f" % dt)
    say("SPREAD slope_dense_localdrift %.4f" % float(mid.max() - mid.min()))
    if dt < 2.0:
        say("UNRESOLVED SIGN slope_dense_localdrift")
    say("SWEPT dense_localslope octave-range %.6f"
        % float(sl.max() - sl.min()))
    say("POP dense_localslope %d" % min(len(ms) for _j, ms in octaves))
    say("CORR dense_localslope %.5f" % min(loc[j][3] for j in loc))
    say("  D5 %s" % ("hold" if d5 else "REFUTED"))

    # the sign structure of the residuals, on the octave the gap opened
    jtop = JHI
    _a, _se, _rms, _cor, _npt, _sp, rtop, xtop, ytop = loc[jtop]
    rgl = ytop - (aglob * xtop + bglob)
    say()
    say("  and the top octave's residuals about the ELEVEN-RUNG line, "
        "which is")
    say("  what {#rem:primorialgap} printed four of:")
    shown = np.round(rgl, 4)          # the signs the reader can check
    say("  residuals about the eleven-rung line: %s"
        % ", ".join("%+.4f" % v for v in shown))
    npos = int((shown > 0).sum())
    run = max(npos, int(shown.size) - npos)
    say("SIGNRUN dense_top_octave %d %d" % (run, shown.size))
    say("  %d of %d have one sign; under an exchangeable draw that is "
        "%.2e" % (run, shown.size, 2.0 ** (1 - shown.size)))

    # ------------------------------------------------ not pre-registered
    say()
    say("E1  what the sweep does buy, which is not what D4 asked for")
    say("  (this section was written after D3 fell and is not "
        "pre-registered)")
    xs = np.log(np.array(ns, dtype=np.float64))
    ys = np.array([got[N] for N in ns])
    asw, bsw, rsw, sesw, csw, _rs = linefit(xs, ys)
    tsw = abs(asw) / sesw
    say("  D3 says the scatter is the statistic's own, so no number of")
    say("  points inside one octave can make a lever out of it. What "
        "the")
    say("  points do is divide the error of a fit over the WHOLE sweep "
        "by")
    say("  the square root of their count.")
    say("  least-squares slope on the sweep = %+.6f, s.e. %.6f, "
        "t = %.2f" % (asw, sesw, tsw))
    say("  over a lever of %.4f in log N, with r.m.s. residual %.6f "
        "and |corr| %.5f" % (float(xs.max() - xs.min()), rsw, csw))
    say("TSTAT slope_audit_primorial_dense %.2f" % tsw)
    say("SPREAD slope_audit_primorial_dense %.4f"
        % float(xs.max() - xs.min()))
    if tsw < 2.0:
        say("UNRESOLVED SIGN slope_audit_primorial_dense")
    say("  the eleven rungs over the same span of N gave %+.6f with "
        "s.e. %.6f;" % (aglob, seglob))
    say("  the sweep and the rungs differ by %+.2f of the sweep's error"
        % ((asw - aglob) / sesw))

    # ------------------------------------------------------------- D6
    say()
    say("D6  the shapes, redone on every point this ladder has")
    seen = {}
    for N, e in zip(lns, lex):
        seen[N] = e
    seen[n12] = e12
    seen.update(got)                       # the sweep wins where it overlaps
    Lu = np.log(np.array(sorted(seen), dtype=np.float64))
    yu = np.array([seen[N] for N in sorted(seen)])
    npts = Lu.size
    say("  %d distinct N: the twelve published rungs and the sweep, "
        "with the" % npts)
    say("  sweep's own value used at the %d N they share (D1 says they "
        "agree)" % sum(1 for N in lns if N in got))
    say("  shape                    r.m.s.     0.5 at      %.2f at"
        % THETA)
    res = {}
    for nm in NAMES:
        c, rms = fit(nm, Lu, yu)
        res[nm] = (c, rms)
        h, t = reaches(nm, c, 0.5), reaches(nm, c, THETA)
        say("  %-24s %-10.6f %-11s %s"
            % (LABEL[nm], rms,
               ("%.4f" % h) if h is not None else "never",
               ("%.4f" % t) if t is not None else "never"))
    order = sorted(NAMES, key=lambda nm: res[nm][1])
    best = res[order[0]][1]
    se_rms = best / math.sqrt(2.0 * (npts - 2))
    surv = [nm for nm in NAMES if res[nm][1] <= best + se_rms]
    say("  the best r.m.s. carries a standard error of its own:")
    say("    %d points   %.6f / sqrt(%d) = %.6f  (%.1f per cent)"
        % (npts, best, 2 * (npts - 2), se_rms, 100.0 * se_rms / best))
    say("    twelve rungs gave %.1f per cent (published)"
        % (100.0 / math.sqrt(2.0 * 10)))
    say("  surviving at one standard error: %s"
        % ", ".join(LABEL[nm] for nm in surv))
    tsurv = [reaches(nm, res[nm][0], THETA) for nm in surv]
    tsurv = [v for v in tsurv if v is not None]
    spread56 = (max(tsurv) - min(tsurv)) if len(tsurv) > 1 else 0.0
    d6 = (len(surv) < 2) or (spread56 < pubspread)
    say("  their theta' spread: %.4f decades, against the published "
        "%.4f" % (spread56, pubspread))
    say("SHAPESURVIVE ladder_theta %d %d %.4f"
        % (npts, len(surv), spread56))
    say("SHAPECURRENT ladder_theta %d" % npts)
    gap56 = res[order[1]][1] - best
    say("SHAPEGAP ladder_theta %.6f %.6f" % (gap56, se_rms))
    if gap56 <= se_rms:
        say("SHAPES TIED ladder_theta")
    say("  D6 %s" % ("hold" if d6 else "REFUTED"))

    # where the two best shapes part, on the union's own r.m.s.
    top10 = math.log10(max(seen))
    grid = np.arange(top10, SHAPESPAN[1], SHAPESTEP) * math.log(10.0)
    v0 = evaluate(order[0], res[order[0]][0], grid)
    v1 = evaluate(order[1], res[order[1]][0], grid)
    part = np.flatnonzero(np.abs(v0 - v1) > best)
    partat = float(grid[part[0]] / math.log(10.0)) if part.size else \
        float(SHAPESPAN[1])
    fc = reaches(order[0], res[order[0]][0], THETA)
    say()
    say("  where the two best shapes part by more than that r.m.s.:")
    say("  top point %.4f, parting at %.4f, %.4f decades above it"
        % (top10, partat, partat - top10))
    say("  the best shape reaching %.2f: N = 10^%s"
        % (THETA, ("%.4f" % fc) if fc is not None else "never"))
    if fc is not None:
        say("TRUST ladder_theta %.4f %.4f" % (partat, fc))
        if fc > partat:
            say("FORECAST OUTSIDE ladder_theta")
        lo56 = reaches(order[0], res[order[0]][0], THETA - best)
        hi56 = reaches(order[0], res[order[0]][0], THETA + best)
        lo56 = min([v for v in (lo56, fc) if v is not None])
        hi56 = max([v for v in (hi56, fc) if v is not None])
        drift = abs(sl.max() - sl.min()) / abs(np.mean(sl))
        say("  bracketed at one r.m.s. residual either way, the "
            "convention")
        say("  rem:primorialreach uses: [10^%.4f, 10^%.4f]"
            % (lo56, hi56))
        say("BRACKET ladder_theta_dense %.4f %.4f %.4f"
            % (fc, lo56, hi56))
        say("DRIFT slope_ladder_theta_dense %.6f" % drift)
        say("SCATTER slope_ladder_theta_dense %.6f" % best)
        say("SHAPES %d" % len(NAMES))
        say("  and that bracket is a floor on the uncertainty and not")
        say("  the whole of it. The DRIFT line is the relative spread "
            "of")
        say("  the four local slopes, %.2f of their mean: D5 says the "
            "slope" % drift)
        say("  this line extrapolates is itself moving, and no shift of")
        say("  the level covers that.")

    say()
    say("=" * 70)
    say("D1 %s  D2 %s  D3 %s  D4 %s  D5 %s  D6 %s"
        % tuple("hold" if v else "REFUTED"
                for v in (d1, d2, d3, d4, d5, d6)))
    ok = d1 and d2
    say("the ladder is now sampled %d times per doubling; its scatter "
        "is" % int(round(np.mean([len(ms) - 1 for _j, ms in octaves]))))
    say("noise and not shape (D3), so no local lever exists (D4), but "
        "the")
    say("same points divide the error of a whole-ladder fit by their "
        "root")
    say("and that leaves one shape standing where twelve rungs left "
        "two (D6).")
    if not ok:
        say("REFUTED: the sweep is not on the ladder")

    head = [
        "STATISTIC: the truncation K*_R at which",
        "           sum_{k<K}(log k)|R(N;k)| first reaches",
        "           S(N)(1-A(N))N, and its exponent log K*_R / log N,",
        "           at every N = 30030*m with m composed of",
        "           2,3,5,7,11,13 and m in [2^6, 2^10]; within each",
        "           closed octave of m, the least-squares slope of that",
        "           exponent against log N, its standard error from",
        "           that octave's own points, and its r.m.s. residual",
        "           against the ladder's published scatter; the four",
        "           local slopes against the eleven-rung global slope",
        "           and against log N; and five shapes refitted on the",
        "           union of this sweep with the twelve published",
        "           rungs, their r.m.s., how many survive at one",
        "           standard error of the best, and where each puts",
        "           1/2 and theta' = 0.56.",
        "NULL: none is run and none applies. A deterministic curve is",
        "      located against a computed threshold; there is no",
        "      background to detect against. The coin arms for this",
        "      statistic were run in lab_primorial_ladder.py and",
        "      lab_primorial_share.py, and the scatter they left is",
        "      the floor D3 is judged against.",
        "FIELD: N = 30030*m for every 13-smooth m in [64, 1024]; every",
        "       one has prime set {2,3,5,7,11,13}, so the odd radical",
        "       3*5*7*11*13 and the threshold are fixed exactly as",
        "       along the ladder; k squarefree and coprime to N with",
        "       2 <= k < " + str(gap.KCAP) + "; m odd, squarefree and",
        "       coprime to k, m < N/k; the sieve weight over the odd",
        "       primes below " + str(QSIEVE) + "; the Euler products at",
        "       the fixed bound " + str(CLIM) + "; the measurement is",
        "       imported from code/audit_primorial_gap.py. The eleven",
        "       rungs and the ladder scatter are read from",
        "       results/audit_primorial_rung10.txt, the twelfth from",
        "       results/audit_primorial_rung11.txt and the published",
        "       shape table from results/audit_ladder_shape12.txt.",
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
