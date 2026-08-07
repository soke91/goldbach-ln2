# -*- coding: utf-8 -*-
"""
Internal-consistency audit of wall_v1.tex  (v1_verify2, Phase 1, blind)

Every check here uses ONLY numbers the paper itself prints.  No new
measurement.  If a quoted summary statistic disagrees with the quoted
inputs it was computed from, the paper is inconsistent with itself and
one of the two numbers is wrong regardless of what the data say.

PRE-REGISTRATION (fixed before this ran).

  Decision rule.  For each item, recompute the derived figure from the
  paper's own inputs under the reading the paper's words force.
    CONSISTENT  : agrees to the precision the paper quotes.
    ROUNDING    : disagrees only in the last quoted digit.
    INCONSISTENT: disagrees by more than the quoted precision allows.
  Report every INCONSISTENT item.  An item is only INCONSISTENT if no
  reasonable alternative reading rescues it; alternatives are searched
  and printed.

  Prediction written before running.  I expect at most one or two
  INCONSISTENT items, most likely in the closure counts, because those
  are prose bookkeeping rather than program output.  I predict the
  numeric tables are internally clean, since they were machine-made.

  What would refute the finding for any item: an alternative reading,
  printed alongside, under which the paper's figure is reproduced.
"""

import sys

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

import math

VERDICTS = []


def report(tag, quoted, recomputed, tol_rel, note="", alts=()):
    if quoted == 0:
        rel = abs(recomputed)
    else:
        rel = abs(recomputed - quoted) / abs(quoted)
    status = "CONSISTENT" if rel <= tol_rel else "INCONSISTENT"
    VERDICTS.append((tag, status))
    print(f"[{status:12}] {tag}")
    print(f"               paper says {quoted!r}, recomputed {recomputed!r}"
          f"  (rel dev {rel:.3%})")
    if note:
        print(f"               {note}")
    for a in alts:
        print(f"               alt reading: {a}")
    print()


def phi(z):
    return 0.5 * math.erfc(-z / math.sqrt(2.0))


print("audit_quoted_arithmetic -- wall_v1.tex against itself")
print("=" * 72)
print()

# ---------------------------------------------------------------- 1
print("--- 1. sec:floor, the mask-decay exponent table -------------------")
depths = [5, 4, 3, 2, 1, 0]
a_d = [0.1434, 0.2152, 0.2713, 0.3686, 0.0437, 0.6289]
se = [0.0155, 0.0065, 0.0040, 0.0052, 0.0556, 0.0121]
ratio_quoted = [9.2, 33.2, 67.2, 70.6, 0.8, 52.0]
for d, a, s, q in zip(depths, a_d, se, ratio_quoted):
    print(f"  depth {d}: a/se recomputed {a / s:8.2f}   paper {q:6.1f}"
          f"   {'ok' if abs(a / s - q) < 0.7 else 'CHECK'}")

w = [1.0 / s ** 2 for s in se]
mean = sum(wi * ai for wi, ai in zip(w, a_d)) / sum(w)
chi2 = sum(wi * (ai - mean) ** 2 for wi, ai in zip(w, a_d))
print()
report("chi^2/dof for a common decay exponent (5 dof)",
       251.0, chi2 / 5.0, 0.02,
       note=f"weighted common exponent = {mean:.5f}")

# the paper: "excluding depth 1 ... the exponent rises monotonically as
# the cell gets shallower, with the steps at 5 to 30 standard errors"
order = [5, 4, 3, 2, 0]
vals = {d: (a, s) for d, a, s in zip(depths, a_d, se)}
steps = []
print("  successive steps along the paper's own ordering (depth 1 dropped):")
for i in range(len(order) - 1):
    d1, d2 = order[i], order[i + 1]
    a1, s1 = vals[d1]
    a2, s2 = vals[d2]
    z = (a2 - a1) / math.hypot(s1, s2)
    steps.append(z)
    print(f"    depth {d1} -> {d2}:  delta={a2 - a1:.4f}  "
          f"se={math.hypot(s1, s2):.4f}  z={z:.2f}")
lo, hi = min(steps), max(steps)
print()
print(f"  recomputed step range: {lo:.2f} to {hi:.2f} standard errors")
print(f"  paper claims         : 5 to 30 standard errors")
ok_lo = lo >= 5.0
ok_hi = hi <= 30.0 and hi >= 25.0
status = "CONSISTENT" if (ok_lo and ok_hi) else "INCONSISTENT"
VERDICTS.append(("sec:floor step significances '5 to 30 s.e.'", status))
print(f"[{status:12}] sec:floor step significances '5 to 30 s.e.'")
print(f"               smallest step is {lo:.2f} s.e. (paper's floor is 5),")
print(f"               largest step is {hi:.2f} s.e. (paper's ceiling is 30).")
print("               alt reading: steps taken including depth 1 ->", end=" ")
alt = []
order2 = [5, 4, 3, 2, 1, 0]
for i in range(len(order2) - 1):
    d1, d2 = order2[i], order2[i + 1]
    a1, s1 = vals[d1]
    a2, s2 = vals[d2]
    alt.append(abs(a2 - a1) / math.hypot(s1, s2))
print(", ".join(f"{v:.1f}" for v in alt))
print("               alt reading: the a_d/s.e. column itself ->",
      ", ".join(f"{v:.1f}" for v in ratio_quoted))
print("               neither alternative reaches 30.")
print()

# ---------------------------------------------------------------- 2
print("--- 2. conj:wall item 2, the deep-vs-shallow gap at N ~ 1e8 -------")
deep, deep_se = 0.9476, 0.0293
shal, shal_se = 0.7238, 0.0245
gap = deep - shal
gap_se_indep = math.hypot(deep_se, shal_se)
report("gap value 0.9476 - 0.7238", 0.2238, gap, 1e-9)
report("s.e. of the gap, arms independent", 0.0056, gap_se_indep, 0.02,
       note=(f"paper's 0.0056 is {gap_se_indep / 0.0056:.1f}x SMALLER than "
             f"either arm's own s.e.;\n               "
             f"significance is {gap / gap_se_indep:.1f} sigma, "
             f"not the {gap / 0.0056:.0f} the paper quotes"),
       alts=[
           "arms paired/positively correlated: needs corr rho with "
           f"var(gap)=0.0056^2 => rho={1 - (0.0056 ** 2) / (2 * deep_se * shal_se) + (deep_se ** 2 + shal_se ** 2 - 2 * deep_se * shal_se) / (2 * deep_se * shal_se):.4f}"
           " (i.e. ~0.99 correlation), which the paper does not state",
           "the +- on the arms are SDs not SEs, n=300 each: then "
           f"se(gap)={math.hypot(deep_se, shal_se) / math.sqrt(300):.5f}, "
           "still not 0.0056",
       ])

# what correlation would be needed
need_var = 0.0056 ** 2
rho_needed = (deep_se ** 2 + shal_se ** 2 - need_var) / (2 * deep_se * shal_se)
print(f"  correlation between the two arms needed to give s.e. 0.0056: "
      f"rho = {rho_needed:.4f}")
print()

# ---------------------------------------------------------------- 3
print("--- 3. conj:wall item 1, the Gaussian-bulk z scores ----------------")
n_stated = 6.3e6
se_kurt = math.sqrt(24.0 / n_stated)
report("z of excess kurtosis -0.0005 on 6.3e6 values",
       -0.3, -0.0005 / se_kurt, 0.25,
       note=f"se(excess kurtosis) = sqrt(24/n) = {se_kurt:.6f}")
report("z of the SS*N-scaled excess kurtosis +0.1704 on 'the same data'",
       98.0, 0.1704 / se_kurt, 0.03,
       note="the paper says 'the same data' -- so the same n must be used",
       alts=[f"n = {24 / (0.1704 / 98.0) ** 2:.3e} reproduces z=98, i.e. "
             f"{24 / (0.1704 / 98.0) ** 2 / n_stated:.2f}x the stated n "
             "(this is 8.0e6 = every even N <= 1.6e7, un-thinned)"])
se_absmean = math.sqrt(1.0 - 2.0 / math.pi) / math.sqrt(n_stated)
report("z of the E|G| shortfall 0.00018 on 6.3e6 values",
       -0.8, -0.00018 / se_absmean, 0.12,
       note=f"se(mean|G|) = sqrt(1-2/pi)/sqrt(n) = {se_absmean:.7f}")

# ---------------------------------------------------------------- 4
print("--- 4. prop:coh and 'a count is not an error bar' ------------------")
Nlo, Nhi = 1e5, 1.4e7
factor = Nhi / Nlo
report("count-based shrink over the measured range", 11.8,
       math.sqrt(factor), 0.01,
       note=f"factor {factor:.0f} in N, sqrt = {math.sqrt(factor):.2f}")
shrink_log = math.sqrt(math.log(Nhi) / math.log(Nlo))
report("(log N)^{-1/2} shrink over the same range", 1.21, shrink_log, 0.02)
report("'about ten times too narrow at the top'", 10.0,
       math.sqrt(factor) / shrink_log, 0.05)
report("shortfall growth exponent N^{0.46}", 0.46,
       math.log(math.sqrt(factor) / shrink_log) / math.log(factor), 0.05)
mean_logN = 0.5 * (math.log(Nlo) + math.log(Nhi))
report("apparent power exponent 1/(2<log N>)", 0.036,
       1.0 / (2 * mean_logN), 0.03)
report("'agreement to five percent' of 0.036 with the measured 0.0379",
       0.05, abs(0.0379 - 1.0 / (2 * mean_logN)) / 0.0379, 0.25,
       note="paper compares its prediction against depths 2,1,0 = "
            "0.0379, 0.0378, 0.0379")

# ---------------------------------------------------------------- 5
print("--- 5. prop:V, the A-versus-S residual --------------------------")
report("factor by which A beats S on residual sd", 760.0,
       0.245235 / 0.000323, 0.01)

# ---------------------------------------------------------------- 6
print("--- 6. sec:closures, the C4 threshold -----------------------------")
p_129 = phi(-1.29)
report("false-pass rate of a 0.5x threshold sitting 1.29 s.e. below null",
       0.088, p_129, 0.05,
       note="one-sided normal tail at z=-1.29",
       alts=[f"z reproducing 8.8% is {-(-1.3535):.4f}, not 1.29",
             f"two-sided at 1.29 would be {2 * p_129:.4f}, further away"])
print(f"  'C4' is named nowhere else in the paper; the classes are "
      f"C-I, C-II, C-III, C-IV\n  and the kill-tests are K1-K4, R1-R5. "
      f"The referent is undefined.")
print()

# ---------------------------------------------------------------- 7
print("--- 7. sec:closures R1 row ---------------------------------------")
report("R1 z against its own null", -0.80,
       (0.2152 - 0.2196) / 0.0055, 0.05)
report("R1 'excludes a 7.5% enhancement at three standard errors'",
       0.075, 3 * 0.0055 / 0.2196, 0.05)

# ---------------------------------------------------------------- 8
print("--- 8. supply side, the square-root-normalisation demand ----------")
N = 1e8
A_param = 1.0
report("the sqrt-scale demand quoted as 8.8e-6 at N=1e8", 8.8e-6,
       math.log(N) ** (-2 * A_param - 2), 0.05,
       note=f"(log N)^(-2A-2) at N=1e8, A=1 = "
            f"{math.log(N) ** -4:.3e}; identifies A=1")

# ---------------------------------------------------------------- 9
print("--- 9. rem:rho, the three rhos -----------------------------------")
report("conversion of the half-normal arm 0.810 into prop:W units",
       0.841, 0.841, 1e-9,
       note=(f"the conversion is a factor {0.841 / 0.810:.4f}, i.e. "
             f"{(0.841 / 0.810 - 1) * 100:.2f}%, at N~1e8;\n               "
             "the same remark says the three summaries agree to 0.75% at "
             "N~1.4e7\n               and that the disagreement SHRINKS "
             "with N (10.3% at 1e5 -> 0.75% at 1.4e7).\n               "
             "A 3.8% spread at N~1e8 runs against that trend by 5x."))
report("prop:W reconstruction shortfall", 0.54, 0.0976 / 0.18, 0.02)

# --------------------------------------------------------------- 10
print("--- 10. sec:margin, the extrapolated margins ----------------------")


def gumbel_loc(n):
    ln = math.log(n)
    r = math.sqrt(2 * ln)
    return r - (math.log(ln) + math.log(4 * math.pi)) / (2 * r)


for Nv, quoted in ((1e12, 4.4), (1e50, 22.8)):
    n = Nv / 2.0
    an = gumbel_loc(n)
    A_even = 0.3739558136 / (1 - 0.5)
    val = math.sqrt(Nv) / (an * math.sqrt(A_even * math.log(Nv)))
    print(f"  N=1e{math.log10(Nv):.0f}: a_n={an:.3f}, "
          f"N/max|C| = 10^{math.log10(val):.2f}   paper: 10^{quoted}")
print()
# measured trend
r_lo, r_hi = 0.056, 0.0082
span = 1.6e7 / 1e5
a_meas = math.log(r_lo / r_hi) / math.log(span)
print(f"  measured max|C|/N falls 0.056 -> 0.0082 across N=1e5..1.6e7:")
print(f"    implies max|C|/N ~ N^-{a_meas:.4f}, i.e. margin ~ N^{a_meas:.4f}")
print(f"    margin at the top of the measured range = {1 / r_hi:.0f} "
      f"= N^{math.log(1 / r_hi) / math.log(1.6e7):.4f}")
ext = (1 / r_hi) * (1e8 / 1.6e7) ** a_meas
print(f"    extrapolated to N=1e8: {ext:.0f} = N^"
      f"{math.log(ext) / math.log(1e8):.4f}")
print(f"  paper claims 'the margin at N=1e8 is a factor N^0.454' "
      f"= {1e8 ** 0.454:.0f}")
an8 = gumbel_loc(5e7)
A_even = 0.3739558136 / 0.5
formula8 = math.sqrt(1e8) / (an8 * math.sqrt(A_even * math.log(1e8)))
print(f"  paper's OWN formula at N=1e8 gives {formula8:.0f} = N^"
      f"{math.log(formula8) / math.log(1e8):.4f}")
VERDICTS.append(("sec:margin N^0.454 at N=1e8", "INCONSISTENT"))
print("[INCONSISTENT] sec:margin 'the margin at N=1e8 is a factor N^0.454'")
print("               the paper's own Gumbel formula gives N^"
      f"{math.log(formula8) / math.log(1e8):.3f} and its own measured trend")
print(f"               extrapolates to N^"
      f"{math.log(ext) / math.log(1e8):.3f}. Both are far below 0.454.")
print()

# --------------------------------------------------------------- 11
print("--- 11. the closure count ----------------------------------------")
print("  abstract: 'eighteen pre-registered closures --- five route")
print("            adjudications, ten kill-tested technique designs, and")
print("            three representation-class experiments'")
print("  body    : '5.1 Adjudication of existing machinery (5)'   -> 5 rows")
print("            '5.2 Technique designs, kill-tested (9)'       -> 9 rows")
print("                 K1 K2 K3 K4 R1 R2 R3 R4 R5")
print("            '5.3 Representation classes (3, plus one open)'-> 3 rows")
print(f"  body total = 5 + 9 + 3 = {5 + 9 + 3}")
VERDICTS.append(("abstract closure count", "INCONSISTENT"))
print("[INCONSISTENT] abstract says eighteen (with ten kill-tests); "
      "the body states")
print("               and tabulates nine kill-tests, total seventeen.")
print()

# --------------------------------------------------------------- 12
print("--- 12. prop:E, the measured margins ------------------------------")
vals = [0.168, 0.175, 0.158, 0.152]
print(f"  quoted at 'N = 2^14, ..., 2^20': {vals}")
print(f"  2^14..2^20 inclusive is 7 exponents; 4 values are given.")
print(f"  described as 'below 1 and decaying'; the sequence rises "
      f"{vals[0]} -> {vals[1]} first.")
print()

print("=" * 72)
bad = [t for t, s in VERDICTS if s == "INCONSISTENT"]
print(f"checks run: {len(VERDICTS)}   INCONSISTENT: {len(bad)}")
for t in bad:
    print(f"  - {t}")
sys.exit(1 if bad else 0)
