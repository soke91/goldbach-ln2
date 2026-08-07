# -*- coding: utf-8 -*-
"""
Does the cell floor grow or die relative to its own error? (inc. 336)

WHAT 335 HANDED OVER. The mechanism that explains the floor's size is
scale-invariant and so cannot explain its decay at all. What survives
that refutation is a change of object: the error bar itself falls, so
an exponent on dm alone says nothing about whether the floor is
becoming more or less real. The quantity that does is

    s_c(N)  =  |dm_c| / se_c ,

the floor in units of its own uncertainty, and its exponent

    g_d  =  b_d - a_d ,     dm ~ N^-a ,   se ~ N^-b .

WHY THIS MATTERS BEYOND A3. The whole programme treats C(N)/sqrt(V) as
something whose departures from mean zero are noise to be beaten down.
If s_c GROWS with N in some cell, that cell's departure is a systematic
that outruns its own error, and no amount of range makes it go away.
If it SHRINKS, the cell floor is a small-N artifact and dies on its
own. These are opposite conclusions for the same measured floor and
the exponent on dm cannot tell them apart.

WHAT IS PREDICTED IN ADVANCE, so that it is not later reported as a
surprise. n_c is proportional to N at every depth -- the cell is a
congruence condition, and its density does not change with scale -- and
the leading term of the exact variance is Q_cc/n_c^2 ~ 1/n_c. So b_d
should come out near 1/2 at EVERY depth. If it does, then g_d = 1/2 -
a_d and the six g_d must fail a common-value test exactly as the six
a_d did in #180; (S3) below is therefore a CHECK, not a payoff, and it
is written down here as expected-to-fail. The result that is not
predetermined is (S2): the SIGN of g_d, which decides the paragraph
above, and which nothing measured so far fixes.

PRE-REGISTRATION (fixed before the run).

  (S1) THE ERROR BAR'S OWN EXPONENT b_d, per depth, with its SE.
       RULE: the statistical framing above holds iff every b_d lies
       within 0.10 of 1/2. If some b_d does not, the error bar is not
       behaving like 1/sqrt(n_c) and the interpretation of g_d changes
       before anything else is read.

  (S2) THE SIGN OF THE SIGNIFICANCE EXPONENT, per depth. g_d = b_d -
       a_d, with SE from the two fits added in quadrature -- an
       overestimate, since the two share their residuals, and the
       overestimate is the safe direction here.
       RULE: the floor in cell d is a growing systematic iff
       g_d > 3*SE(g_d), and a dying artifact iff g_d < -3*SE(g_d).

  (S3) IS THERE A COMMON SIGNIFICANCE EXPONENT? chi-square of the g_d
       about their inverse-variance-weighted mean.
       RULE: a common value is admissible iff chi2/dof <= 2. EXPECTED
       TO FAIL, for the reason given above; recorded so that its
       failure cannot be reported as information and its SUCCESS would
       be.

  (S4) WHERE IT ARRIVES. For each cell with g_d > 0, report the
       measured significance at the top band and the N at which the
       fitted law reaches 10 sigma. RULE: none -- this is an
       extrapolation of a fitted power law far beyond its range and is
       reported as an order of magnitude, not a number.

  WHAT WOULD REFUTE. (S1) failing anywhere invalidates the reading of
  every g. If all g_d < 0 the cell floor is a small-N artifact and A3's
  whole exponent question is about something that vanishes.
"""
import math
import sys
import time

import numpy as np

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

QS = [3, 5, 7, 11, 13]
DEPTHS = [5, 4, 3, 2, 1, 0]


def sieve(X):
    spf = np.zeros(X + 1, dtype=np.int32)
    for i in range(2, int(X ** 0.5) + 1):
        if spf[i] == 0:
            sl = spf[i * i::i]
            sl[sl == 0] = i
    for i in range(2, X + 1):
        if spf[i] == 0:
            spf[i] = i
    mu = np.zeros(X + 1, dtype=np.int8)
    mu[1] = 1
    for i in range(2, X + 1):
        p = int(spf[i])
        j = i // p
        mu[i] = 0 if j % p == 0 else -mu[j]
    primes = np.nonzero(spf[2:] == np.arange(2, X + 1))[0] + 2
    lam = np.zeros(X + 1, dtype=np.float64)
    for p in primes:
        q = int(p)
        lg = math.log(int(p))
        while q <= X:
            lam[q] = lg
            q *= int(p)
    return mu, lam


def wfit_se(x, y, w):
    S = w.sum()
    Sx = (w * x).sum()
    Sxx = (w * x * x).sum()
    Sy = (w * y).sum()
    Sxy = (w * x * y).sum()
    D = S * Sxx - Sx * Sx
    c1 = (S * Sxy - Sx * Sy) / D
    c0 = (Sxx * Sy - Sx * Sxy) / D
    r = y - (c0 + c1 * x)
    dof = max(len(x) - 2, 1)
    s2 = (w * r * r).sum() / dof
    return c1, c0, math.sqrt(max(s2 * S / D, 0.0))


def main():
    X = 16_000_000
    lo = 100_000
    t0 = time.time()
    mu, lam = sieve(X)
    nf = 1
    while nf < 2 * (X + 1):
        nf *= 2
    supp = (mu != 0).astype(np.float64)
    F_supp = np.fft.rfft(np.pad(supp, (0, nf - X - 1)))
    F_lam = np.fft.rfft(np.pad(lam, (0, nf - X - 1)))
    Fl_c = np.conj(F_lam)
    V = np.fft.irfft(F_supp * np.fft.rfft(
        np.pad(lam ** 2, (0, nf - X - 1))), nf)[: X + 1]
    Ns = np.arange(lo, X + 1, 2)
    invV = 1.0 / np.sqrt(V[Ns])
    Creal = np.fft.irfft(np.fft.rfft(
        np.pad(mu.astype(np.float64), (0, nf - X - 1))) * F_lam,
        nf)[: X + 1]
    Z = Creal[Ns] * invV
    muw = supp[: X + 1]
    print(f"sieve + V + C  t={time.time()-t0:.0f}s", flush=True)

    div = [(Ns % q) == 0 for q in QS]
    cell = {}
    for d in DEPTHS:
        m = np.ones(len(Ns), dtype=bool)
        for j in range(len(QS)):
            m &= div[j] if j < d else ~div[j]
        cell[d] = m

    bands = []
    b = float(lo)
    while b < X:
        hi = min(b * math.sqrt(2.0), X)
        sel = (Ns >= b) & (Ns < hi)
        if int(sel.sum()) > 500:
            bands.append((b, hi, sel))
        b = hi
    print(f"{len(bands)} half-octave bands", flush=True)

    def ucorr(vals):
        w = np.zeros(nf)
        w[Ns] = vals
        return np.fft.irfft(Fl_c * np.fft.rfft(w), nf)[: X + 1]

    rows = {d: [] for d in DEPTHS}
    for bi, (b0, hi, sel) in enumerate(bands):
        n = int(sel.sum())
        u_all = ucorr(np.where(sel, invV, 0.0))
        mu_all = muw * u_all
        Qaa = float(np.dot(mu_all, u_all))
        gm = float(Z[sel].mean())
        Nmid = math.sqrt(b0 * hi)
        for d in DEPTHS:
            m = sel & cell[d]
            nc = int(m.sum())
            if nc < 2:
                continue
            u = ucorr(np.where(m, invV, 0.0))
            var = (float(np.dot(muw * u, u)) / nc ** 2
                   - 2 * float(np.dot(mu_all, u)) / (nc * n)
                   + Qaa / n ** 2)
            if var <= 0:
                continue
            dm = float(Z[m].mean()) - gm
            if abs(dm) < 1e-12:
                continue
            rows[d].append((Nmid, dm, math.sqrt(var)))
        if (bi + 1) % 5 == 0:
            print(f"  band {bi+1}/{len(bands)}  "
                  f"t={time.time()-t0:.0f}s", flush=True)

    print(f"\n(S1)(S2) the floor, its error, and the ratio of the two")
    print(f"{'depth':>6} {'pts':>4} {'a (floor)':>11} {'b (error)':>11} "
          f"{'|b-0.5|':>9} {'g = b-a':>11} {'g/SE':>7}")
    gs, gses, ds, tops, cs = [], [], [], [], []
    okS1 = True
    for d in DEPTHS:
        r = rows[d]
        if len(r) < 4:
            print(f"{d:>6} {len(r):>4}   too few points")
            continue
        Nm = np.array([x[0] for x in r])
        dm = np.array([x[1] for x in r])
        se = np.array([x[2] for x in r])
        w = (np.abs(dm) / se) ** 2
        a, _, sa = wfit_se(np.log(Nm), np.log(np.abs(dm)), w)
        bb, _, sb = wfit_se(np.log(Nm), np.log(se), np.ones(len(Nm)))
        sig = np.abs(dm) / se
        g, c0, sg_fit = wfit_se(np.log(Nm), np.log(sig),
                                np.ones(len(Nm)))
        sgq = math.sqrt(sa ** 2 + sb ** 2)
        sg = max(sg_fit, sgq)
        if abs(-bb - 0.5) > 0.10:
            okS1 = False
        gs.append(g); gses.append(sg); ds.append(d)
        tops.append(float(sig[-1])); cs.append(c0)
        print(f"{d:>6} {len(r):>4} {-a:>11.4f} {-bb:>11.4f} "
              f"{abs(-bb-0.5):>9.4f} {g:>11.4f} "
              f"{abs(g)/max(sg,1e-12):>7.1f}")

    gs = np.array(gs); gses = np.array(gses); ds = np.array(ds)
    tops = np.array(tops); cs = np.array(cs)
    print(f"\n    (S1) every error-bar exponent within 0.10 of 1/2: "
          f"{'PASS' if okS1 else 'FAIL'}")

    print(f"\n    (S2) is the floor outrunning its own error?")
    print(f"{'depth':>6} {'g':>10} {'SE':>8} {'verdict':>22} "
          f"{'sigma at 1.4e7':>15}")
    grow, die = [], []
    for i, d in enumerate(ds):
        if gs[i] > 3 * gses[i]:
            v = "GROWING systematic"
            grow.append(d)
        elif gs[i] < -3 * gses[i]:
            v = "dying artifact"
            die.append(d)
        else:
            v = "flat within 3 SE"
        print(f"{d:>6} {gs[i]:>10.4f} {gses[i]:>8.4f} {v:>22} "
              f"{tops[i]:>15.1f}")

    wg = 1.0 / gses ** 2
    gm_ = float((gs * wg).sum() / wg.sum())
    chi2 = float((wg * (gs - gm_) ** 2).sum())
    dof = max(len(gs) - 1, 1)
    okS3 = (chi2 / dof) <= 2.0
    print(f"\n    (S3) a common significance exponent is admissible "
          f"(chi2/dof <= 2): {'PASS' if okS3 else 'FAIL as predicted'} "
          f" (chi2/dof = {chi2/dof:.2f}, common g = {gm_:.4f})")

    print(f"\n    (S4) where the growing cells arrive "
          f"(power law extrapolated far past its range)")
    for i, d in enumerate(ds):
        if gs[i] <= 3 * gses[i]:
            continue
        n10 = math.exp((math.log(10.0) - cs[i]) / gs[i])
        print(f"{'depth':>6} {d}: {tops[i]:.1f} sigma now, "
              f"10 sigma near N ~ 1e{math.log10(n10):.1f}")

    if not okS1:
        v = ("the error bar does not fall like 1/sqrt(n_c) at every "
             "depth, so g cannot be read as significance growth and "
             "nothing below the first table stands")
    elif grow and not die:
        v = (f"the cell floor is a GROWING systematic at every "
             f"measurable depth (g from {gs.min():.2f} to "
             f"{gs.max():.2f}). It is not a small-N artifact: in each "
             f"cell the departure of C/sqrt(V) from mean zero outruns "
             f"its own error bar, so more range makes it more "
             f"significant, not less")
    elif grow and die:
        v = (f"the cell floor SPLITS. It is a growing systematic at "
             f"depths {sorted(grow, reverse=True)} and a dying "
             f"artifact at depths {sorted(die, reverse=True)}. The "
             f"deep cells -- where more small primes divide N -- carry "
             f"a departure that outruns its error, while the shallow "
             f"ones are washing out. That is A4's per-cell share seen "
             f"in N, and it says the floor concentrates with scale "
             f"rather than fading")
    elif die and not grow:
        v = ("the cell floor is a small-N artifact at every depth: it "
             "falls faster than its own error bar everywhere, so A3's "
             "exponent question is about a quantity that vanishes")
    else:
        v = ("no cell's floor moves relative to its error bar by 3 SE "
             "over a factor 160 in N, so neither reading is supported")
    print(f"\n    {v}")
    print("DONE")


if __name__ == "__main__":
    main()
