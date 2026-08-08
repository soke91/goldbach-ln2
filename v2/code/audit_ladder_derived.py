# -*- coding: utf-8 -*-
r"""
The model failed on its constant. Put the measured constant back in.

WHAT IS AT STAKE

OPEN item 1 is the location of theta' = 0.56, and after
{#rem:laddershape12} and {#rem:shapetrust} the bottleneck is named:
the two surviving shapes are FITTED, they are separated by 2.7845
decades, and the ladder can only be read 0.34 decades past its top
rung. More rungs narrow the field and not the forecast. What would
resolve it is a DERIVED shape.

{#rem:laddermodel} tried one and it failed -- but it failed at a
single identified place. The heuristic assumes |R| ~ c_R sqrt(N/k)
with c_R ~ gamma sqrt(log N); on the primorial ladder c_R/sqrt(log N)
is not constant but climbs 0.3724 to 0.6139, and c_R is fitted
directly as (log N)^1.3838 with correlation 0.98535. So the model was
refuted with its replacement in hand and nobody put it back.

Doing so is not a fit. With |R| = gamma (log N)^rho sqrt(N/k) the
crossing sum(log k)|R| = S(1-A)N has a closed form, and since
sum_{k<K}(log k)k^{-1/2} ~ 2 sqrt(K)(log K - 2), writing K = N^e gives

    e = 1 - 2(1+rho) loglog N / log N + d / log N

with the leading coefficient DERIVED from rho and only d free. That
is a one-parameter shape against every rival's two, and 1 - c loglog
N/log N was already one of the five candidates -- but without the
d/log N term, which is why it fitted at r.m.s. 0.03934.

BACKS: Remark {#rem:ladderderived} in paper/theorem_A.md.

PRE-REGISTERED PREDICTIONS (written before this script was run)

  M1  The control: this script recomputes the ladder from its own
      sieve and reproduces the eleven published exponents to 1e-4,
      the twelfth (0.5099) to 1e-4, and the eleven published c_R to
      0.001.
  M2  The law: fitting log c_R on loglog N reproduces the published
      exponent 1.3838 and correlation 0.98535 to 0.001 on eleven
      rungs, and on twelve the exponent stays within 0.05 of it.
  M3  Correcting the constant fixes the crossings: with gamma and rho
      fitted once to c_R and nothing fitted to the exponents, the
      model's K* is within 5 per cent of the measured at every rung.
      That is {#rem:laddermodel}'s own J2 cap, so the two are
      comparable.
  M4  And it fixes the shape: the derived shape above, with the
      coefficient forced to 2(1+rho) and only d free, fits the twelve
      exponents with an r.m.s. no worse than 1.5 times the best
      fitted shape's 0.00370.

REFUTATION RULE (fixed before the run)

  M1  REFUTED beyond the caps -- not the same ladder, and nothing
      below may be compared with the published rungs.
  M2  REFUTED beyond 0.001 (control) or 0.05 (twelve rungs).
  M3  REFUTED beyond 5 per cent at any rung. That would say the
      failure of {#rem:laddermodel}'s J2 is not the constant's law
      but the sqrt(N/k) profile underneath it.
  M4  REFUTED beyond 1.5 times 0.00370. That is the one that
      matters: the derived coefficient would then be measurably
      wrong, and OPEN item 1 would have no derived shape after the
      one repair its own diagnosis pointed at.

  M1 and M2 gate M3 and M4.

  NO NULL IS RUN and none applies, for the reason
  {#rem:laddermodel} gives: a deterministic model with no free
  parameter is compared with a measured crossing of the same measured
  sum at twelve N, and there is no background to detect against. The
  sign controls for R on this ladder were run in
  lab_primorial_ladder.py.
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
OUT = os.path.join(RES, "audit_ladder_derived.txt")

BASE = 30030                       # 2*3*5*7*11*13
LADDER = [BASE * (1 << j) for j in range(12)]
KCAP = 100_000
QSIEVE = 30
CLIM = 4_000_000


def primes_upto(n):
    s = np.ones(n + 1, dtype=bool)
    s[:2] = False
    for p in range(2, int(n ** 0.5) + 1):
        if s[p]:
            s[p * p::p] = False
    return np.flatnonzero(s).astype(np.int64)


def sieves(n, qs):
    """Lambda in float32, mu in int8, and a bitmask of the small qs

    Independent of audit_ladder_model.py in the weight: there the
    weight is a product over q of 0 or q/(q-1), which is a constant
    times the indicator that no admissible q divides N-mk. Here the
    indicator is read off one bitmask instead of nine remainders.
    """
    pr = primes_upto(n)
    lgp = np.log(pr.astype(np.float64))
    lam = np.zeros(n + 1, dtype=np.float32)
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
    del pr, lgp
    mu = np.ones(n + 1, dtype=np.int8)
    cof = np.arange(n + 1, dtype=np.int32)
    for p in primes_upto(int(math.isqrt(n))):
        p = int(p)
        mu[p::p] = -mu[p::p]
        if p * p <= n:
            mu[p * p::p * p] = 0
        q = p
        while q <= n:
            cof[q::q] //= p
            if q > n // p:
                break
            q *= p
    big = cof > 1
    del cof
    mu[big] = -mu[big]
    del big
    mu[0] = 0
    rmask = np.zeros(n + 1, dtype=np.uint16)
    for i, q in enumerate(qs):
        rmask[::q] |= np.uint16(1 << i)
    return lam, mu, rmask


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


def read_published():
    """the eleven rungs' exponents and c_R, and the twelfth exponent"""
    src = io.open(os.path.join(RES, "audit_ladder_model.txt"),
                  encoding="utf-8").read()
    i = src.index("  N            c_R        sqrt(log N)   "
                  "c_R/sqrt(log N)")
    cR = {}
    for ln in src[i:].splitlines()[1:]:
        f = ln.split()
        if len(f) < 4 or not f[0].isdigit():
            break
        cR[int(f[0])] = float(f[1])
    i = src.index("J1  the control")
    ex = {}
    for ln in src[i:].splitlines()[2:]:
        f = ln.split()
        if len(f) < 4 or not f[0].isdigit():
            break
        ex[int(f[0])] = float(f[2])
    law = re.search(r"c_R ~ \(log N\)\^([\d.]+), correlation ([\d.]+)",
                    src)
    src2 = io.open(os.path.join(RES, "audit_primorial_rung11.txt"),
                   encoding="utf-8").read()
    m = re.search(r"^  (\d+) +([\d.]+) +([\d.]+) +", src2[
        src2.index("N            log10 N   exponent"):], re.M)
    ex[int(m.group(1))] = float(m.group(3))
    src3 = io.open(os.path.join(RES, "audit_ladder_shape12.txt"),
                   encoding="utf-8").read()
    best = min(float(m[1]) for m in re.findall(
        r"^  \S.*?\s(\d\.\d{5})\s+(\d\.\d{5})\s+"
        r"(\d+\.\d{2})\s+(\d+\.\d{4})\s*$",
        src3[src3.index("J2/J3"):], re.M))
    i4 = src3.index("J4  and what")
    ans = {}
    for nm, v in re.findall(
            r"^  (a \+ b log(?: log)? N) +([\d.]+)\s*$",
            src3[i4:src3.index("spread", i4)], re.M):
        ans[nm] = float(v)
    return ex, cR, float(law.group(1)), float(law.group(2)), best, ans


def rung(N, lam, mu, rmask, sqf, qs, thrc):
    """H, P, beta, |R|, the crossing and c_R at one rung"""
    PN = factor_set(N)
    ks, Hs, Ps = [], [], []
    for k in range(2, KCAP):
        if not sqf[k] or any(k % q == 0 for q in PN):
            continue
        M = (N - 1) // k
        if M < 2:
            continue
        ms = np.arange(1, M + 1, 2, dtype=np.int64)
        ms = ms[sqf[ms]]
        fk = factor_set(k)
        for q in fk:
            if q > 2:
                ms = ms[ms % q != 0]
        if ms.size == 0:
            continue
        vals = N - ms * k
        g = mu[ms].astype(np.float64)
        keep = np.uint16(0)
        C = 1.0
        for i, q in enumerate(qs):
            if k % q:
                keep |= np.uint16(1 << i)
                C *= q / (q - 1.0)
        ks.append(k)
        Hs.append(float(np.dot(lam[vals], g)))
        Ps.append(C * float(g[(rmask[vals] & keep) == 0].sum()))
    ks = np.array(ks, dtype=np.int64)
    H = np.array(Hs)
    P = np.array(Ps)
    lw = np.log(ks.astype(float))
    beta = float((H * P).sum() / (P * P).sum())
    aR = np.abs(H - beta * P)
    j = int(np.searchsorted(np.cumsum(lw * aR), thrc * N))
    kstar = int(ks[min(j, ks.size - 1)])
    e = math.log(kstar) / math.log(N)
    scale = np.sqrt(N / ks.astype(float))
    sel = ks <= kstar
    cR = float((aR[sel] / scale[sel]).mean())
    return ks, lw, scale, kstar, e, cR


def fitline(x, y):
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    b = float(((x - x.mean()) * (y - y.mean())).sum()
              / ((x - x.mean()) ** 2).sum())
    a = float(y.mean() - b * x.mean())
    r = float(((x - x.mean()) * (y - y.mean())).sum()
              / math.sqrt(((x - x.mean()) ** 2).sum()
                          * ((y - y.mean()) ** 2).sum()))
    return a, b, r


def main():
    lines = []

    def say(s=""):
        print(s)
        lines.append(s)

    pubex, pubcR, publaw, pubcorr, bestrms, ans = read_published()
    say("read %d published exponents, %d published c_R, the law "
        "(log N)^%.4f" % (len(pubex), len(pubcR), publaw))
    say("  correlation %.5f, the best fitted r.m.s. on twelve %.5f, "
        "and %d shape answers" % (pubcorr, bestrms, len(ans)))

    qs = [int(q) for q in primes_upto(QSIEVE) if q > 2]
    NMAX = max(LADDER)
    say("sieving to %d over %d small primes ..." % (NMAX, len(qs)))
    lam, mu, rmask = sieves(NMAX, qs)
    sqf = mu != 0

    artin, twin = 1.0, 2.0
    for p in primes_upto(CLIM):
        p = int(p)
        artin *= 1.0 - 1.0 / (p * (p - 1.0))
        if p > 2:
            twin *= 1.0 - 1.0 / (p - 1.0) ** 2

    rows = []
    for N in LADDER:
        A_, S_ = artin, twin
        for q in sorted(factor_set(N)):
            A_ /= (1.0 - 1.0 / (q * (q - 1.0)))
            if q > 2:
                S_ *= (1.0 + 1.0 / (q - 2.0))
        thrc = S_ * (1.0 - A_)
        ks, lw, scale, kstar, e, cR = rung(
            N, lam, mu, rmask, sqf, qs, thrc)
        rows.append((N, thrc, ks, lw, scale, kstar, e, cR))
        say("  N = %-9d K*_R %-6d exp %.4f  c_R %.4f  thr %.6f"
            % (N, kstar, e, cR, thrc))

    # ------------------------------------------------------------- M1
    say()
    say("M1  the control: the ladder recomputed")
    say("  N            exp        published   diff      c_R      "
        "published   diff")
    m1 = True
    for N, thrc, ks, lw, scale, kstar, e, cR in rows:
        d1 = abs(e - pubex[N])
        d2 = abs(cR - pubcR[N]) if N in pubcR else 0.0
        if d1 >= 1e-4 or d2 >= 1e-3:
            m1 = False
        say("  %-12d %-10.4f %-11.4f %-9.6f %-8.4f %-11s %.6f"
            % (N, e, pubex[N], d1, cR,
               ("%.4f" % pubcR[N]) if N in pubcR else "-", d2))
    say("  M1 %s   (cap 0.0001 in the exponent, cap 0.001 in c_R)"
        % ("hold" if m1 else "REFUTED"))

    # ------------------------------------------------------------- M2
    say()
    say("M2  the law the constant actually follows")
    L1 = np.array([math.log(r[0]) for r in rows])
    L2 = np.log(L1)
    cRs = np.array([r[7] for r in rows])
    g11, rho11, r11 = fitline(L2[:11], np.log(cRs[:11]))
    g12, rho12, r12 = fitline(L2, np.log(cRs))
    say("  eleven rungs   c_R ~ (log N)^%.4f   correlation %.5f"
        % (rho11, r11))
    say("  published      c_R ~ (log N)^%.4f   correlation %.5f"
        % (publaw, pubcorr))
    say("  twelve rungs   c_R ~ (log N)^%.4f   correlation %.5f"
        % (rho12, r12))
    m2 = (abs(rho11 - publaw) < 1e-3 and abs(r11 - pubcorr) < 1e-3
          and abs(rho12 - publaw) < 0.05)
    say("  M2 %s   (cap 0.001 on eleven, cap 0.05 on twelve)"
        % ("hold" if m2 else "REFUTED"))
    gam = math.exp(g12)
    say("  so gamma = %.6f and rho = %.6f, both fitted to c_R alone"
        % (gam, rho12))

    # ------------------------------------------------------------- M3
    say()
    say("M3  the crossings with the corrected constant")
    say("  N            measured   corrected  ratio    "
        "published ratio")
    m3 = True
    rat = []
    for i, (N, thrc, ks, lw, scale, kstar, e, cR) in enumerate(rows):
        cmod = gam * L1[i] ** rho12
        cum = np.cumsum(lw * cmod * scale)
        jm = int(np.searchsorted(cum, thrc * N))
        km = int(ks[min(jm, ks.size - 1)])
        r = km / kstar
        rat.append(r)
        if abs(r - 1.0) >= 0.05:
            m3 = False
        say("  %-12d %-10d %-10d %-8.4f %s"
            % (N, kstar, km, r, "-"))
    say("  the ratios span %.4f to %.4f" % (min(rat), max(rat)))
    say("  M3 %s   (cap 5 per cent, the same cap the uncorrected "
        "model was read at)" % ("hold" if m3 else "REFUTED"))
    say("PERN ladder_derived_ratio %d %.4f %.4f"
        % (len(rat), min(rat), max(rat)))

    # ------------------------------------------------------------- M4
    say()
    say("M4  the derived shape against the fitted ones")
    c = 2.0 * (1.0 + rho12)
    ex = np.array([r[6] for r in rows])
    u, x = L2 / L1, 1.0 / L1
    y = ex - 1.0 + c * u
    d = float((x * y).sum() / (x * x).sum())
    fit = 1.0 - c * u + d * x
    rms = float(np.sqrt(((ex - fit) ** 2).mean()))
    say("  the derived coefficient is 2(1+rho) = %.4f and the one "
        "free constant" % c)
    say("  fits at d = %.4f" % d)
    say("  N            exponent   derived    residual")
    for i, N in enumerate(LADDER):
        say("  %-12d %-10.4f %-10.4f %+.4f" % (N, ex[i], fit[i],
                                               ex[i] - fit[i]))
    say("  r.m.s. %.5f against the best fitted shape's %.5f, a factor "
        "of %.2f" % (rms, bestrms, rms / bestrms))
    m4 = rms <= 1.5 * bestrms
    say("  M4 %s   (cap 1.5 times %.5f)"
        % ("hold" if m4 else "REFUTED", bestrms))

    # what the data would have wanted, and whether it can say
    say()
    say("  post hoc, not pre-registered: freeing the coefficient")
    X = np.vstack([-u, x]).T
    sol, *_ = np.linalg.lstsq(X, ex - 1.0, rcond=None)
    cf, df = float(sol[0]), float(sol[1])
    fit2 = 1.0 - cf * u + df * x
    rms2 = float(np.sqrt(((ex - fit2) ** 2).mean()))
    say("  free coefficient %.4f against the derived %.4f, r.m.s. "
        "%.5f" % (cf, c, rms2))
    ru = float(((u - u.mean()) * (x - x.mean())).sum()
               / math.sqrt(((u - u.mean()) ** 2).sum()
                           * ((x - x.mean()) ** 2).sum()))
    say("  but the two regressors correlate at %.5f, so the "
        "coefficient and" % ru)
    say("  the constant are not separable on this ladder: the free "
        "fit buys")
    say("  %.5f of r.m.s. for one more parameter."
        % (rms - rms2))
    say("CORR ladder_derived_regressors %.5f" % ru)
    if abs(ru) >= 0.99:
        say("COEFF NOT SEPARABLE ladder_derived")

    # where it puts 1/2 and 0.56
    say()
    say("  and where the derived shape reaches the two levels")
    say("  level    log10 N   fitted shapes say")
    for lev in (0.50, 0.56):
        lo, hi = 5.0, 200.0
        for _ in range(200):
            mid = 0.5 * (lo + hi)
            v = 1.0 - c * math.log(mid) / mid + d / mid
            if v < lev:
                lo = mid
            else:
                hi = mid
        say("  %.2f     %-9.4f %s"
            % (lev, 0.5 * (lo + hi) / math.log(10.0),
               ", ".join("%.4f" % v for v in sorted(ans.values()))
               if lev > 0.55 else "-"))

    say()
    say("  the arithmetic and the budget, declared:")
    say("  %d N, %d distinct odd radical: %d"
        % (len(LADDER),
           len(set(tuple(sorted(q for q in factor_set(N) if q > 2))
                   for N in LADDER)),
           BASE // 2))

    ok = m1 and m2 and m3 and m4
    say()
    say("=" * 70)
    say("REFUTED" if not ok else "all four hold")

    hdr = [
        "STATISTIC: at each rung of N = 30030*2^j, j = 0..11, the",
        "           crossing K*_R of sum(log k)|H-beta P| against",
        "           S(N)(1-A(N))N, its exponent, and the constant",
        "           c_R = mean |R|/sqrt(N/k) below it; then c_R's own",
        "           law gamma (log N)^rho, the crossings that law",
        "           predicts with nothing fitted to the exponents,",
        "           and the shape e = 1 - 2(1+rho) loglog N/log N +",
        "           d/log N it derives, read against the fitted",
        "           shapes' r.m.s.",
        "NULL: none is run and none applies. A deterministic model",
        "      with no freedom in the crossing is compared with a",
        "      measured crossing of the same measured sum at twelve",
        "      N; there is no background to detect against. The sign",
        "      controls for R on this ladder were run in",
        "      lab_primorial_ladder.py.",
        "FIELD: N = 30030*2^j, j = 0..11, one odd radical throughout",
        "       so the budget is constant; k squarefree and coprime",
        "       to N with 2 <= k < 100000; m odd, squarefree, coprime",
        "       to the odd part of k, m < N/k; the sieve weight uses",
        "       the odd primes up to 30, read off one bitmask; beta",
        "       refitted as sum(HP)/sum(P^2) on that k-range; S(N)",
        "       and A(N) from Euler products at the fixed bound",
        "       4000000. The published exponents, c_R and law are",
        "       read from results/audit_ladder_model.txt, the twelfth",
        "       exponent from results/audit_primorial_rung11.txt and",
        "       the fitted shapes from",
        "       results/audit_ladder_shape12.txt.",
        "",
    ]
    io.open(OUT, "w", encoding="utf-8", newline="\n").write(
        "\n".join(hdr + lines) + "\n")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
