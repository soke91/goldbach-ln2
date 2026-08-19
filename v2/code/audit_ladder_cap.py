# -*- coding: utf-8 -*-
r"""
The whole ladder recomputed at one cap: does the branch survive it?

WHAT IS AT STAKE

Remark {#rem:rung16} found that the k-cap is inside the statistic and
not beside it.  K*_R is the truncation at which
sum_{k<K}(log k)|R(N;k)| reaches S(N)(1-A(N))N with R = H - beta*P,
and beta is a least-squares fit over the same k-range, so the cap the
search runs to is also the cap beta is estimated on.  Widening it
moves beta, moves R, and moves K*_R at every rung the ladder has
already printed.

Rung 16 handled that by freezing beta's window at the published
k < 100000 and widening only the search, which leaves every published
integer alone -- and that is the right convention for adding a rung,
because it is the only one under which the new rung and the old ones
are the same statistic.  It is not an answer to the other question:
**would the branch's conclusions be the same on a ladder whose cap
had been generous from the start?**

Sixteen rungs, a resolved curvature, three out-of-sample hits and a
settled extrapolation all rest on exponents computed at one arbitrary
cap.  Nothing has ever asked what they would have been at another.
This script asks.

The sweep is cheap for a reason worth stating.  H(N;k) and P(N;k) do
not depend on the cap at all -- the cap enters only through which k
are fitted and how far the cumulative sum is read.  So one pass to
k < 1000000 at each N yields every cap's answer as a subset
operation, and four caps cost one pass rather than four.

BACKS: Remark {#rem:laddercap} in paper/theorem_A.md.

PRE-REGISTERED PREDICTIONS (written before this script was run)

  U1  The control.  At the published cap this code returns the
      published K*_R of results/audit_primorial_rung10.txt and
      rung11 through rung15 exactly, at every one of those rungs.
  U2  The cap converges.  At each of the top five rungs the exponent
      moves less between caps 400000 and 1000000 than between caps
      100000 and 400000.  The statistic has a cap-free limit and the
      published cap is an approximation to it, not an arbitrary point
      on a drift.
  U3  And the approximation is inside the noise.  At every rung the
      exponent moves less than the ladder's published scatter 0.0037
      between the published cap and 1000000.
  U4  The escalation survives.  On the uniform ladder at cap 1000000
      the margin over 1/2 grows at each of the top six rungs, as it
      does on the published one.
  U5  The curvature survives.  On the uniform ladder the coefficient
      of (log N)^2 is positive with t > 2, as {#rem:laddercurve}
      found at the published cap.
  U6  The crossing survives.  The uniform ladder's 0.56 crossing
      falls inside the bracket [9.4216, 9.9340] that
      {#rem:rung15} published.

REFUTATION RULE (fixed before the run)

  U1  REFUTED by a single rung whose K*_R differs.  Then this is not
      the ladder's statistic and nothing below is comparable with
      it.  THIS ONE GATES.
  U2  REFUTED if any of the five moves more on the second step than
      the first.  That would say the cap is not converging and the
      exponent is a function of a parameter with no natural value --
      the ladder would then have to publish the cap as a declared
      level, and no single number for the exponent would be
      defensible.
  U3  REFUTED if any rung moves by the scatter or more.  Then the
      published exponents are not robust to the cap and every number
      in this branch carries an undeclared systematic larger than
      the noise it is judged against.
  U4  REFUTED by one rung of the top six that does not grow.  The
      six growing margins are what the branch's "the barrier stays
      crossed" rests on; if they are an artefact of the cap the
      reading changes.
  U5  REFUTED if the coefficient is not positive, or if t <= 2.
      {#rem:laddercurve} corrected two earlier readings to land on a
      resolved upward curve; a curvature that dissolves when the cap
      moves would withdraw that correction.
  U6  REFUTED if the crossing falls outside the published bracket.
      The bracket is the branch's most-quoted interval and it has
      never been tested against anything but its own draws.

  U1 gates.  U2 to U6 are the measurement and do not gate: a
  refutation among them is a finding about the ladder, not a fault
  in this script.

  NO NULL IS RUN and none applies.  A deterministic curve is located
  against a computed threshold at four caps; there is no background
  to detect against.  The coin arms for this statistic were run in
  lab_primorial_ladder.py and lab_primorial_share.py.
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
OUT = os.path.join(RES, "audit_ladder_cap.txt")

CLIM = 4_000_000                    # the fixed Euler bound (G20)
SEED = 20260823
DRAWS = 4000
TARGET = 0.56
CAPS = (100_000, 200_000, 400_000, 1_000_000)
TOP = 16                            # the highest rung recomputed here


def module(name):
    p = os.path.join(CODE, name + ".py")
    spec = importlib.util.spec_from_file_location(name, p)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


R11 = module("audit_primorial_rung11")
R16 = module("audit_primorial_rung16")
primes_upto = R11.primes_upto
BASE = R11.BASE
BLOCK = R16.BLOCK


def hp(N, lp, mu, vmask, qs, kmax, block=BLOCK):
    """H(N;k) and P(N;k) for every admissible k < kmax

    Neither depends on the cap; the cap only chooses which of these
    are fitted and how far the cumulative sum is read.  That is why
    one pass serves every cap.
    """
    PN = R11.factor_set(N)
    ks, Hs, Ps = [], [], []
    for k in range(2, kmax):
        if mu[k] == 0:
            continue
        if any(k % q == 0 for q in PN):
            continue
        M = (N - 1) // k
        if M < 2:
            continue
        kb, ck = 0, 1.0
        for i, q in enumerate(qs):
            if k % q == 0:
                kb |= 1 << i
            else:
                ck *= q / (q - 1.0)
        keepbits = np.uint16(~kb & 0xFFFF)
        drop = [q for q in R11.factor_set(k) if q > 2]
        h, pp, seen = 0.0, 0.0, 0
        for lo in range(1, M + 1, block):
            hi = min(lo + block, M + 1)
            ms = np.arange(lo if lo % 2 else lo + 1, hi, 2,
                           dtype=np.int64)
            if ms.size == 0:
                continue
            ms = ms[mu[ms] != 0]
            for q in drop:
                ms = ms[ms % q != 0]
            if ms.size == 0:
                continue
            seen += ms.size
            vals = N - ms * k
            g = mu[ms].astype(np.float64)
            pv = lp[vals]
            nz = pv != 0
            if nz.any():
                h += float((np.log(pv[nz].astype(np.float64))
                            * g[nz]).sum())
            keep = (vmask[vals >> 1] & keepbits) == 0
            pp += float(g[keep].sum())
        if seen == 0:
            continue
        ks.append(k)
        Hs.append(h)
        Ps.append(ck * pp)
    return (np.array(ks, dtype=np.int64), np.array(Hs),
            np.array(Ps))


def kstar_at(N, ks, H, P, artin, twin, cap):
    """the ladder's K*_R with both the fit and the search at one cap"""
    PN = R11.factor_set(N)
    A_, S_ = artin, twin
    for q in sorted(PN):
        A_ /= (1.0 - 1.0 / (q * (q - 1.0)))
        if q > 2:
            S_ *= (1.0 + 1.0 / (q - 2.0))
    m = ks < cap
    if not m.any():
        return None
    k_, H_, P_ = ks[m], H[m], P[m]
    beta = float((H_ * P_).sum() / (P_ * P_).sum())
    cum = np.cumsum(np.log(k_.astype(np.float64))
                    * np.abs(H_ - beta * P_))
    thr = S_ * (1.0 - A_) * N
    j = int(np.searchsorted(cum, thr))
    if j >= k_.size:
        return None
    kk = int(k_[j])
    return kk, math.log(kk) / math.log(N), beta, int(m.sum())


def read_published():
    """the published rungs, their K*_R where printed, and the bracket"""
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
    star = {}
    for j in (11, 12, 13, 14, 15):
        s = io.open(os.path.join(RES, "audit_primorial_rung%d.txt" % j),
                    encoding="utf-8").read()
        N = BASE * (1 << j)
        m = re.search(r"^  N = " + str(N) +
                      r"\s+thr [\d.]+\s+#k \d+\s+beta [\d.]+\s+"
                      r"K\*_R (\d+)\s+exp ([\d.]+)\s*$", s, re.M)
        ns.append(N)
        ex.append(float(m.group(2)))
        star[N] = int(m.group(1))
    s15 = io.open(os.path.join(RES, "audit_primorial_rung15.txt"),
                  encoding="utf-8").read()
    scat = float(re.search(r"^FLOOR primorial_rung15 ([\d.]+)\s*$",
                           s15, re.M).group(1))
    m = re.search(r"^BRACKET ladder_quadratic16_theta_prime ([\d.]+) "
                  r"([\d.]+) ([\d.]+)\s*$", s15, re.M)
    return (ns, ex, dec, star, scat, float(m.group(2)),
            float(m.group(3)))


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
        sys.stdout.flush()
        lines.append(s)

    ns, ex, dec, star, scat, blo, bhi = read_published()
    say("read %d published rung exponents, the scatter %.4f, and the "
        "published" % (len(ns), scat))
    say("  0.56 bracket [%.4f, %.4f] from "
        "results/audit_primorial_rung15.txt" % (blo, bhi))
    say("SEED %d" % SEED)
    say("DRAWS %d" % DRAWS)
    say("CAPS %s" % " ".join(str(c) for c in CAPS))

    rungs = [BASE * (1 << j) for j in range(TOP + 1)]
    assert rungs[:len(ns)] == ns, "the published rungs are not this ladder"
    NEW = rungs[TOP]
    qs = [int(q) for q in primes_upto(R11.QSIEVE) if q > 2]
    say()
    say("sieving to %d once; every rung and every cap is read off it"
        % NEW)
    lp, mu = R16.prime_and_mu_block(NEW)
    vmask = R16.residue_mask_odd(NEW, qs)
    say("BYTES resident_arrays %d"
        % (lp.nbytes + mu.nbytes + vmask.nbytes))
    artin, twin = 1.0, 2.0
    assert CLIM == R11.CLIM
    for p in primes_upto(CLIM):
        p = int(p)
        artin *= 1.0 - 1.0 / (p * (p - 1.0))
        if p > 2:
            twin *= 1.0 - 1.0 / (p - 1.0) ** 2
    say("  the Euler products at the fixed bound %d: Artin %.9f, "
        "twin %.9f" % (CLIM, artin, twin))

    kmax = max(CAPS)
    tab = {}
    say()
    say("  rung  N             "
        + "".join("%-18s" % ("cap %d" % c) for c in CAPS))
    say("  " + " " * 19
        + "".join("%-18s" % "K*_R     exponent" for c in CAPS))
    for j, N in enumerate(rungs):
        ks, H, P = hp(N, lp, mu, vmask, qs, kmax)
        row = {}
        for c in CAPS:
            row[c] = kstar_at(N, ks, H, P, artin, twin, c)
        tab[j] = row
        say("  %-5d %-13d %s"
            % (j, N, "".join(
                ("%-8d %-9.6f " % (row[c][0], row[c][1]))
                if row[c] else "%-8s %-9s " % ("none", "---")
                for c in CAPS)))
        flat = len(set(row[c][0] for c in CAPS if row[c])) == 1
        if flat and all(row[c] for c in CAPS):
            lines[-1] = lines[-1].rstrip() + "   cap-invariant"
            print("    (cap-invariant)")

    # -------------------------------------------------------------- U1
    say()
    say("U1  the control at the published cap")
    u1 = True
    bad = []
    for j, N in enumerate(rungs):
        if N not in star:
            continue
        r = tab[j][CAPS[0]]
        if r is None or r[0] != star[N]:
            u1 = False
            bad.append(j)
    say("  the %d rungs whose K*_R is printed reproduce %s"
        % (len(star), "exactly" if u1 else "with %d differing" % len(bad)))
    for j, N in enumerate(rungs):
        if N in star:
            r = tab[j][CAPS[0]]
            say("    rung %-3d K*_R %-8s against the published %d"
                % (j, r[0] if r else "none", star[N]))
    say("  U1 %s   (cap: exact on K*_R)" % ("hold" if u1 else "REFUTED"))
    if not u1:
        io.open(OUT, "w", encoding="utf-8", newline="\n").write(
            "\n".join(lines) + "\n")
        raise SystemExit(1)

    # -------------------------------------------------------------- U2
    say()
    say("U2  does the cap converge?")
    u2 = True
    for j in range(TOP - 4, TOP + 1):
        a = tab[j][CAPS[0]]
        b = tab[j][CAPS[2]]
        c_ = tab[j][CAPS[3]]
        if not (a and b and c_):
            u2 = False
            say("    rung %-3d incomplete" % j)
            continue
        s1 = abs(b[1] - a[1])
        s2_ = abs(c_[1] - b[1])
        if not s2_ < s1:
            u2 = False
        say("    rung %-3d %d->%d moves %.6f, %d->%d moves %.6f  %s"
            % (j, CAPS[0], CAPS[2], s1, CAPS[2], CAPS[3], s2_,
               "shrinks" if s2_ < s1 else "GROWS"))
    say("  U2 %s   (cap: the first step's size)"
        % ("hold" if u2 else "REFUTED"))

    # -------------------------------------------------------------- U3
    say()
    say("U3  is the cap's cost inside the noise?")
    worst, wj = 0.0, -1
    for j in range(TOP + 1):
        a, c_ = tab[j][CAPS[0]], tab[j][CAPS[3]]
        if not (a and c_):
            continue
        d = abs(c_[1] - a[1])
        if d > worst:
            worst, wj = d, j
    u3 = worst < scat
    say("  the largest move between cap %d and cap %d is %.6f, at "
        "rung %d" % (CAPS[0], CAPS[3], worst, wj))
    say("CAPCOST audit_ladder_cap %.6f %.4f" % (worst, scat))
    say("  against the ladder's scatter %.4f, a ratio of %.2f"
        % (scat, worst / scat))
    say("  U3 %s   (cap: the ladder's scatter)"
        % ("hold" if u3 else "REFUTED"))

    # -------------------------------------------------------------- U4
    say()
    say("U4  does the escalation survive a uniform cap?")
    big = CAPS[3]
    margins = []
    for j in range(TOP - 5, TOP + 1):
        r = tab[j][big]
        margins.append(None if r is None else r[1] - 0.5)
    u4 = all(m is not None for m in margins) and all(
        margins[i] > margins[i - 1] for i in range(1, len(margins)))
    say("  margins over 1/2 on the top six rungs at cap %d: %s"
        % (big, ", ".join("%.4f" % m if m is not None else "none"
                          for m in margins)))
    say("  U4 %s   (cap: each rung against the one below)"
        % ("hold" if u4 else "REFUTED"))

    # -------------------------------------------------------------- U5
    say()
    say("U5  does the curvature survive?")
    xs, ys = [], []
    for j in range(TOP + 1):
        r = tab[j][big]
        if r is None:
            continue
        xs.append(math.log(rungs[j]))
        ys.append(r[1])
    x = np.array(xs)
    y = np.array(ys)
    c, cov, s2v, rms = quadfit(x, y)
    se2 = math.sqrt(float(cov[2, 2]))
    t2 = c[2] / se2
    u5 = c[2] > 0.0 and t2 > 2.0
    say("  on %d rungs at cap %d the (log N)^2 coefficient is "
        "%+.8f +- %.8f, t = %.2f" % (x.size, big, c[2], se2, t2))
    say("  r.m.s. residual %.4f" % rms)
    say("TSTAT ladder_cap_quadratic %.2f" % t2)
    if abs(t2) < 2.0:
        say("UNRESOLVED SIGN ladder_cap_quadratic")
    say("SPREAD ladder_cap_quadratic %.4f" % (x.max() - x.min()))
    say("  U5 %s   (cap: positive and t > 2)"
        % ("hold" if u5 else "REFUTED"))

    # -------------------------------------------------------------- U6
    say()
    say("U6  does the crossing survive?")
    p = cross(c, TARGET)
    rng = np.random.default_rng(SEED)
    draws = rng.multivariate_normal(c, cov, size=DRAWS)
    vv = [cross(dd, TARGET) for dd in draws]
    vv = [w / math.log(10.0) for w in vv
          if w is not None and w > x.max()]
    lo = float(np.percentile(vv, 2.5))
    hi = float(np.percentile(vv, 97.5))
    here = p / math.log(10.0) if p else None
    u6 = here is not None and blo <= here <= bhi
    say("  the uniform ladder reaches 0.56 at log10 N = %.4f, bracket "
        "[%.4f, %.4f]" % (here, lo, hi))
    say("  the published bracket is [%.4f, %.4f]" % (blo, bhi))
    say("BRACKET ladder_cap_theta_prime %.4f %.4f %.4f"
        % (here, lo, hi))
    say("DRIFT ladder_cap_theta_prime %.4f" % abs(here - (blo + bhi) / 2.0))
    say("SHAPES 1")
    say("SCATTER slope_audit_ladder_cap %.4f" % rms)
    say("  U6 %s   (cap: the published bracket)"
        % ("hold" if u6 else "REFUTED"))
    say("  no forecast is made from this; {#rem:shapepower} is why.")

    say()
    say("=" * 70)
    say("U1 %s  U2 %s  U3 %s  U4 %s  U5 %s  U6 %s"
        % tuple("hold" if v_ else "REFUTED"
                for v_ in (u1, u2, u3, u4, u5, u6)))

    head = [
        "STATISTIC: the ladder's level exponent log K*_R / log N at",
        "           N = 30030*2^j for j = 0..16, recomputed at four",
        "           k-caps 100000, 200000, 400000 and 1000000 with",
        "           both the beta fit and the truncation search at",
        "           the same cap; the exponent's movement between",
        "           caps against the ladder's published scatter; and",
        "           the margin, the (log N)^2 coefficient and the",
        "           0.56 crossing of the uniform ladder at the",
        "           largest cap, against what the published ladder",
        "           reports for each.",
        "NULL: none is run and none applies. A deterministic curve is",
        "      located against a computed threshold at four caps;",
        "      there is no background to detect against. The coin",
        "      arms for this statistic were run in",
        "      lab_primorial_ladder.py and lab_primorial_share.py.",
        "      The bracket is drawn from the fit's own parameter",
        "      covariance with the fixed SEED.",
        "FIELD: N = 30030*2^j for j = 0..16, the odd radical",
        "       3*5*7*11*13 fixed so the threshold is constant along",
        "       the ladder; k squarefree and coprime to N with",
        "       2 <= k < cap, beta fitted on the same range; m odd,",
        "       squarefree and coprime to k, m < N/k; the sieve",
        "       weight over the odd primes below 30; the Euler",
        "       products at the fixed bound 4000000. The sieve, the",
        "       block statistic and the packed arrays are imported",
        "       from code/audit_primorial_rung16.py, whose C1",
        "       compares them against code/audit_primorial_rung11.py",
        "       elementwise; the published exponents come from",
        "       results/audit_primorial_rung10.txt and rung11",
        "       through rung15.",
        "",
    ]
    io.open(OUT, "w", encoding="utf-8", newline="\n").write(
        "\n".join(head + lines) + "\n")
    print("\nwrote %s" % OUT)
    if not u1:
        raise SystemExit(1)
    return 0


if __name__ == "__main__":
    sys.exit(main())
