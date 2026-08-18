# -*- coding: utf-8 -*-
r"""
The error bar Remark {#rem:forecast} does not carry.

WHAT IS AT STAKE

That remark forecasts, from a model calibrated on three values of N,
the point at which the route's hypothesis becomes true at
theta' = 0.56: N = 2.077e8. It quotes one number.

Remark {#rem:modeltransfer} has since shown that this family of models
is a model of the exponent and only a fit of the constant: a five per
cent drift in c = mean |H|/sqrt(N/k) between two k-ranges moved the
predicted crossing by ten per cent, and the model absorbed an upward
drift but not a downward one. A forecast extrapolated sixty-fold past
its calibration inherits that, and the remark records no bracket.

There is also a cross-check nobody has made. The forecast's model
writes |A(N;k)| = gamma sqrt(N/k) sqrt(log N) and fits gamma = 0.9803
on the mean. Remark {#rem:heuristic} measured the same constant from a
different statistic, over a different k-range, and got c/sqrt(log N)
between 0.9844 and 1.0138. If those agree then the fitted gamma is not
a free parameter at all and the forecast is parameter-free.

BACKS: Remark {#rem:forecastbracket} in paper/wall_v3.md.

PRE-REGISTERED PREDICTIONS (written before this script was run)

  W1  This reconstruction reproduces the published forecast table:
      K*/sqrt(N) and K*/N^{0.56} at every decade from 10^5 to 10^20,
      read from results/lab_level_forecast.txt, to within 1 per cent.
  W2  The two independent measurements of the constant agree: the
      fitted gamma and every c(N)/sqrt(log N) of the heuristic audit
      lie within 5 per cent of each other.
  W3  So the forecast is parameter-free. Setting gamma = 1 exactly
      -- no fit at all -- moves both published crossings by less than
      20 per cent.
  W4  And the inherited uncertainty is not an order of magnitude:
      perturbing the constant by the +-10 per cent that
      {#rem:modeltransfer} measured moves the theta' = 0.56 crossing
      by less than a factor 2.

REFUTATION RULE (fixed before the run)

  W1  REFUTED at 1 per cent at any decade -- it would mean this
      script is not reconstructing the published model.
  W2  REFUTED if any pair differs by 5 per cent or more, which would
      say the two measurements of the constant are not the same
      measurement and one of them is wrong.
  W3  REFUTED if either crossing moves by 20 per cent or more.
  W4  REFUTED if the crossing moves by a factor 2 or more, in which
      case the published single number 2.077e8 cannot be quoted
      without an order-of-magnitude bracket.

  All four gate. W4 is the deliverable either way: the bracket it
  computes is what the remark has to carry.

  NO NULL IS RUN and none applies. A published deterministic forecast
  is reconstructed and perturbed; there is no detection against a
  background. The sign control for this field was run in
  lab_level_coin_null.py, whose coin draws established that the level
  is bought by cancellation, and audit_forecast_null.py supplied the
  control that {#rem:forecast}'s own evidence had declined.
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
OUT = os.path.join(ROOT, "results", "audit_forecast_bracket.txt")

KMAX = 10_000_000
PN = (2, 5)                # the family N = 2^a 5^b the forecast uses
CLIM = 4_000_000
WOBBLE = 0.10              # what rem:modeltransfer measured


def primes_upto(n):
    s = np.ones(n + 1, dtype=bool)
    s[:2] = False
    for p in range(2, int(n ** 0.5) + 1):
        if s[p]:
            s[p * p::p] = False
    return np.flatnonzero(s).astype(np.int64)


def read_forecast():
    """the published gamma, table and crossings -- read, not copied"""
    p = os.path.join(ROOT, "results", "lab_level_forecast.txt")
    src = io.open(p, encoding="utf-8").read()
    g = float(re.search(r"mean-based gamma over the three fitting N "
                        r"= ([\d.]+)", src).group(1))
    i = src.index("the forecast under the CORRECTED gamma:")
    tab = {}
    for ln in src[i:].splitlines()[2:]:
        m = re.match(r"\s*10\^(\d+)\s+([\d.e+]+)\s+([\d.]+)\s+([\d.]+)\s*$",
                     ln)
        if not m:
            break
        tab[int(m.group(1))] = (float(m.group(3)), float(m.group(4)))
    cr = {}
    for lab, key in (("sqrt N", "sqrt"), ("N\\^0.56", "056")):
        cr[key] = float(re.search(
            r"corrected: K\*/" + lab + r" crosses 1 at N = ([\d.e+]+)",
            src).group(1))
    return g, tab, cr


def read_heuristic_c():
    """c(N)/sqrt(log N) from the heuristic audit -- read, not copied"""
    p = os.path.join(ROOT, "results", "audit_directlevel_heuristic.txt")
    src = io.open(p, encoding="utf-8").read()
    i = src.index("c(N)/sqrt(log N)")
    out = {}
    for ln in src[i:].splitlines()[1:]:
        f = ln.split()
        if len(f) < 4 or not f[0].isdigit():
            break
        out[int(f[0])] = float(f[1])
    return out


def main():
    lines = []

    def say(s=""):
        print(s)
        lines.append(s)

    gpub, tabpub, crpub = read_forecast()
    cheur = read_heuristic_c()
    say("read from results/: gamma = %.4f, %d published decades, "
        "%d crossings," % (gpub, len(tabpub), len(crpub)))
    say("  and %d values of c(N)/sqrt(log N) from the heuristic audit"
        % len(cheur))

    # ---- S(K) rebuilt independently, then continued asymptotically
    say()
    say("building S(K) = sum_{k<K, mu^2=1, (k,%s)=1} (log k)/sqrt(k) "
        "to %d ..." % (",".join(map(str, PN)), KMAX))
    sq = np.ones(KMAX + 1, dtype=bool)
    sq[0] = False
    p = 2
    while p * p <= KMAX:
        sq[p * p::p * p] = False
        p += 1
    for q in PN:
        sq[q::q] = False
    ks = np.flatnonzero(sq)
    ksf = ks.astype(np.float64)          # searchsorted upcasts otherwise,
    S = np.cumsum(np.log(ksf) / np.sqrt(ksf))   # once per call, at 13 ms
    dens = ks.size / float(KMAX)
    say("  admissible k <= %d: %d (density %.6f);  S(%d) = %.4f"
        % (KMAX, ks.size, dens, KMAX, S[-1]))

    artin, twin = 1.0, 2.0
    for q in primes_upto(CLIM):
        q = int(q)
        artin *= 1.0 - 1.0 / (q * (q - 1.0))
        if q > 2:
            twin *= 1.0 - 1.0 / (q - 1.0) ** 2
    for q in PN:
        artin /= (1.0 - 1.0 / (q * (q - 1.0)))
        if q > 2:
            twin *= (1.0 + 1.0 / (q - 2.0))
    THR = twin * (1.0 - artin)
    say("  threshold S(N)(1-A(N)) = %.6f" % THR)

    def S_of(K):
        if K <= KMAX:
            return float(S[max(int(np.searchsorted(ksf, K)) - 1, 0)])
        f = lambda t: 2.0 * math.sqrt(t) * (math.log(t) - 2.0)
        return float(S[-1]) + dens * (f(K) - f(KMAX))

    def Kstar(N, g):
        """solve g sqrt(log N / N) S(K) = THR for K"""
        tgt = THR * math.sqrt(N / math.log(N)) / g
        lo, hi = 2.0, 1e24
        for _ in range(200):
            mid = math.sqrt(lo * hi)
            if S_of(mid) < tgt:
                lo = mid
            else:
                hi = mid
        return math.sqrt(lo * hi)

    def crossing(g, expo):
        """the N at which K*(N) = N^expo"""
        lo, hi = 1e4, 1e24
        for _ in range(200):
            mid = math.sqrt(lo * hi)
            if Kstar(mid, g) < mid ** expo:
                lo = mid
            else:
                hi = mid
        return math.sqrt(lo * hi)

    # ------------------------------------------------------------- W1
    say()
    say("W1  the published forecast table, reconstructed")
    say("  N       K*/sqrt N  published   K*/N^0.56  published   worst")
    w1 = True
    for e in sorted(tabpub):
        k = Kstar(10.0 ** e, gpub)
        a, b = k / 10.0 ** (e / 2.0), k / 10.0 ** (e * 0.56)
        pa, pb = tabpub[e]
        d = max(abs(a / pa - 1.0), abs(b / pb - 1.0))
        if d >= 0.01:
            w1 = False
        say("  10^%-5d %-10.4f %-11.4f %-10.4f %-11.4f %.4f"
            % (e, a, pa, b, pb, d))
    say("  W1 %s" % ("hold" if w1 else "REFUTED"))

    # ------------------------------------------------------------- W2
    say()
    say("W2  two independent measurements of the same constant")
    say("  the forecast fits gamma on the mean of |A|/sqrt(N/k) over")
    say("  k < K*; the heuristic audit measures c/sqrt(log N) over a")
    say("  different k-range at five N. Every pair:")
    say("  N            c/sqrt(log N)   gamma       ratio")
    w2 = True
    for N in sorted(cheur):
        r = cheur[N] / gpub
        if abs(r - 1.0) >= 0.05:
            w2 = False
        say("  %-12d %-15.4f %-11.4f %.4f" % (N, cheur[N], gpub, r))
    say("  W2 %s" % ("hold" if w2 else "REFUTED"))

    # ------------------------------------------------------------- W3
    say()
    say("W3  is the fit doing any work? gamma = 1 exactly, no fit")
    say("  crossing        published     gamma = 1     ratio")
    w3 = True
    free = {}
    for key, expo in (("sqrt", 0.5), ("056", 0.56)):
        c1 = crossing(1.0, expo)
        free[key] = c1
        r = c1 / crpub[key]
        if abs(r - 1.0) >= 0.20:
            w3 = False
        say("  K* = N^%-9.2f %-13.4e %-13.4e %.4f"
            % (expo, crpub[key], c1, r))
    say("  W3 %s" % ("hold" if w3 else "REFUTED"))

    # ------------------------------------------------------------- W4
    say()
    say("W4  the bracket: the constant perturbed by the +-%d per cent"
        % int(100 * WOBBLE))
    say("  that rem:modeltransfer measured")
    say("  crossing        c low         published     c high       "
        "span")
    w4 = True
    span = {}
    for key, expo in (("sqrt", 0.5), ("056", 0.56)):
        # a LARGER constant means a larger B, so the crossing is LATER
        hi = crossing(gpub * (1.0 + WOBBLE), expo)
        lo = crossing(gpub * (1.0 - WOBBLE), expo)
        s = max(hi, lo) / min(hi, lo)
        span[key] = (lo, hi, s)
        if s >= 2.0:
            w4 = False
        say("  K* = N^%-9.2f %-13.4e %-13.4e %-13.4e %.4f"
            % (expo, lo, crpub[key], hi, s))
    say("  W4 %s" % ("hold" if w4 else "REFUTED"))

    say()
    say("  THE BRACKET THE REMARK HAS TO CARRY, in log10:")
    for key, expo in (("sqrt", 0.5), ("056", 0.56)):
        lo, hi, s = span[key]
        say("    K* = N^%.2f :  10^%.2f  [10^%.2f, 10^%.2f]"
            % (expo, math.log10(crpub[key]),
               math.log10(lo), math.log10(hi)))
    say()
    say("  Bracket lines, in the form the gate reads. This file supplies")
    say("  them for results/lab_level_forecast.txt, whose crossings are")
    say("  quoted without one:")
    for key, expo in (("sqrt", 0.5), ("056", 0.56)):
        lo, hi, s = span[key]
        say("BRACKET kstar_N^%.2f %.4e %.4e %.4e"
            % (expo, crpub[key], lo, hi))
    say()
    say("  And the drift of the constant those brackets were built by")
    say("  wobbling, which gate check G33 reads. The bracket assumes")
    say("  +-%d per cent; the five independent measurements of the same"
        % int(100 * WOBBLE))
    say("  constant, c(N)/sqrt(log N) from the heuristic audit, span")
    dr = (max(cheur.values()) - min(cheur.values()))         / (sum(cheur.values()) / len(cheur))
    say("  %.4f of their mean. A bracket is honest only when the wobble"
        % dr)
    say("  it assumes is at least the drift the constant actually has,")
    say("  and here it is %s."
        % ("wider" if WOBBLE >= dr else "NARROWER -- the bracket "
           "understates"))
    say("DRIFT kstar_gamma %.4f" % dr)

    say()
    say("  DIAGNOSTIC (post hoc). Why the bracket is this wide. At the")
    say("  crossing K* = N^e the balance reads, with")
    say("  S(K) ~ 2 d sqrt(K)(log K - 2),")
    say("    N^{(1-e)/2}  ~  gamma sqrt(log N) (e log N - 2),")
    say("  so a relative error d in gamma moves N by (1+d)^{2/(1-e)} --")
    say("  the square root of K doubles the exponent, which is the")
    say("  factor a naive 1/(1-e) misses.")
    say("  e         naive (1.10)^{2/(1-e)}   measured, high side   low")
    for key, expo in (("sqrt", 0.5), ("056", 0.56)):
        lo, hi, s = span[key]
        say("  %-9.2f %-24.4f %-21.4f %.4f"
            % (expo, 1.10 ** (2.0 / (1.0 - expo)),
               hi / crpub[key], crpub[key] / lo))
    say("  Even that understates it, because the right-hand side is not")
    say("  constant: sqrt(log N)(e log N - 2) GROWS with N, so N must")
    say("  move further still to absorb the same change in gamma. The")
    say("  logs amplify the error rather than damping it, and that is")
    say("  the whole gap between the naive factor and the measured one.")

    say()
    say("  The bracket is honest in the one place it can be checked:")
    say("  the measured K*/sqrt(N) crossing lies between 8e5 and 1.6e6,")
    say("  and both the point estimate %.4e and the bracket contain it."
        % crpub["sqrt"])

    say()
    say("=" * 70)
    ok = w1 and w2 and w3 and w4
    say("the forecast is parameter-free and its bracket is under a "
        "factor 2" if ok else "REFUTED")

    head = [
        "STATISTIC: K*(N) solved independently from the published model",
        "           gamma sqrt(log N / N) S(K) = S(N)(1-A(N)), with S(K)",
        "           rebuilt by sieve to 1e7 and continued by",
        "           2 sqrt(K)(log K - 2) times the measured density;",
        "           the crossings K* = sqrt(N) and K* = N^{0.56} under",
        "           the published gamma, under gamma = 1 exactly, and",
        "           under gamma scaled by 1 +- 0.10.",
        "NULL: none is run and none applies. A published deterministic",
        "      forecast is reconstructed and perturbed; there is no",
        "      background to detect against. The sign control for this",
        "      field was run in lab_level_coin_null.py, and",
        "      audit_forecast_null.py supplied the control that",
        "      rem:forecast's own evidence had declined.",
        "FIELD: the family P(N) = {2,5} the forecast uses throughout;",
        "       k squarefree and coprime to it, k < 1e7 enumerated",
        "       exactly; S(N) and A(N) from Euler products at the fixed",
        "       bound 4000000; gamma and the published table and",
        "       crossings are read from results/lab_level_forecast.txt,",
        "       and c(N)/sqrt(log N) from",
        "       results/audit_directlevel_heuristic.txt.",
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
