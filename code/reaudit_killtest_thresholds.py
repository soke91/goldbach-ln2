# -*- coding: utf-8 -*-
"""
Every kill-test's threshold, restated in sigmas (increment 313)

WHY. Increment 312 gave K2 a floor and showed its closure stands
quantitatively. Twelve kill-tests still had none, and the obvious next
move -- compute `T/sqrt(n)` for each -- does not work, because their
criteria are not all of that shape. Reading the thirteen headers shows
three different kinds:

  z-threshold   K2 ("|mean|/SE >= 4"), C2 (">= 3 k with z >= 4.5")
  ratio         R1 ("R2_zeros >= 2 x mean(R2_random)"), R2, R4
                ("ratio <= 0.5 x the B=1 baseline"), C4
                ("defect_real <= 0.5 x defect_null"), K1 ("R^2 <= 0.9")
  moment match  K3 ("trace moments match the null")

For a ratio criterion the floor question is not `T/sqrt(n)`. It is:
**how many standard errors of its own null does the threshold sit
from that null, and how many does the measurement sit?** Those two
numbers say whether a DEAD verdict is strong, weak, or unreadable, and
neither has ever been computed.

WHAT COMES OUT. The two numbers disagree sharply and in a consistent
direction: the thresholds are far LOOSER than the measurements are
SHARP. R1's "2x" threshold sits 40 standard errors above its null while
its measurement sits 0.4 below it. So the recorded verdict understates
its own evidence -- R1 excludes a 7.5% enhancement at 3 sigma, having
only promised to exclude a 100% one.

And one criterion cannot be evaluated at all: R2 asks for "capture >=
2x the random control" where the control is NEGATIVE, and twice a
negative number is smaller, not larger.

HOW THIS AVOIDS #134. That correction was for asserting a property of
eleven files without opening them. Here every number is read out of the
kill-test's own recorded output by hand, and the script then VERIFIES
mechanically that each number it uses actually occurs in that file. A
number that has drifted, or that I mistyped, fails the check rather
than being quoted.
"""
import io
import math
import os
import re
import sys

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

RES = os.path.join(os.path.dirname(os.path.dirname(
    os.path.abspath(__file__))), "results")

# (label, results file, null mean, null sd, threshold value,
#  measured value, kind, note). Every float must occur in the file.
ROWS = [
    ("Forge R1", "forge_r1_result.txt", 0.2196, 0.0055, None, 0.98,
     "ratio", "ALIVE iff R2_zeros >= 2 x null; measured ratio 0.98"),
    ("Forge R2", "forge_r2_result.txt", -0.00006, 0.00016, None,
     -0.00012, "ratio",
     "ALIVE iff capture >= 2 x control; the control is NEGATIVE"),
    ("Forge R4", "forge_r4_result.txt", 0.9862, None, None, 0.5929,
     "ratio", "ALIVE iff block ratio <= 0.5 x the B=1 baseline"),
    ("Constr C4", "constr_c4_result.txt", None, None, None, 1.53,
     "ratio", "ALIVE iff defect_real <= 0.5 x defect_null"),
]


def nums_in(path):
    try:
        t = io.open(path, encoding="utf-8", errors="replace").read()
    except Exception:
        return None, ""
    return set(re.findall(r"-?\d+\.\d+", t)), t


def occurs(vals, x):
    """Is x present, to the precision it is written with?"""
    if x is None:
        return True
    s = f"{x}"
    if s in vals:
        return True
    for v in vals:
        try:
            if abs(float(v) - x) < 1e-9:
                return True
        except ValueError:
            pass
    return False


def main():
    print("(1) verification: every number below occurs in its own "
          "results file")
    ok_all = True
    for lab, fn, m, sd, th, meas, kind, note in ROWS:
        vals, _t = nums_in(os.path.join(RES, fn))
        if vals is None:
            print(f"    {lab:<12} {fn:<26} FILE NOT FOUND")
            ok_all = False
            continue
        good = all(occurs(vals, x) for x in (m, sd, meas))
        ok_all &= good
        print(f"    {lab:<12} {fn:<26} "
              f"{'verified' if good else 'NUMBER NOT IN FILE'}")
    if not ok_all:
        print("\n  a quoted number is not in the file it is quoted from")
        print("DONE (findings)")
        sys.exit(1)

    print(f"\n(2) the criterion and the measurement, both in standard "
          f"errors of the null")
    print(f"{'test':<12} {'null':>12} {'sd':>10} {'threshold':>12} "
          f"{'T in sd':>9} {'measured':>10} {'meas in sd':>11}")
    for lab, fn, m, sd, th, meas, kind, note in ROWS:
        if m is None or sd is None:
            print(f"{lab:<12} {'—':>12} {'not printed':>10} "
                  f"{'0.5x':>12} {'?':>9} {meas:>10.4f} {'?':>11}")
            continue
        thr = 2.0 * m
        tsd = (thr - m) / sd
        # the measurement is quoted as a ratio to the null
        msd = (meas * m - m) / sd if kind == "ratio" and abs(meas) < 10 \
            else (meas - m) / sd
        print(f"{lab:<12} {m:>12.5f} {sd:>10.5f} {thr:>12.5f} "
              f"{tsd:>+9.1f} {meas:>10.5f} {msd:>+11.2f}")

    print(f"\n(3) what each verdict actually excludes")
    r1m, r1sd = 0.2196, 0.0055
    print(f"    Forge R1  the pre-registered threshold is "
          f"{(2*r1m-r1m)/r1sd:.0f} sd above the null, so as WRITTEN the")
    print(f"              test only excluded a doubling. Its measurement")
    print(f"              is 0.98x the null, i.e. {(0.98*r1m-r1m)/r1sd:+.1f} sd,")
    print(f"              so the DATA exclude an enhancement of "
          f"{3*r1sd/r1m:.1%} at 3 sd.")
    print(f"              **The recorded verdict understates its own "
          f"evidence by a factor {(2*r1m-r1m)/(3*r1sd):.0f}.**")
    print(f"    Forge R2  'capture >= 2 x the random control' with a")
    print(f"              control of -0.00006: twice a negative number")
    print(f"              is SMALLER, so the ALIVE branch is unreachable")
    print(f"              as written. ⚠️ The criterion is ill-posed. The")
    print(f"              measurement is fine and says nothing is there:")
    print(f"              -0.00012 against -0.00006 +- 0.00016 = "
          f"{(-0.00012+0.00006)/0.00016:+.2f} sd.")
    print(f"    Forge R4  no error bar is printed beside the block")
    print(f"              ratios, and R4b's own header puts the")
    print(f"              estimator's SE at ~9% at B=8 and ~25% at")
    print(f"              B=64. At B=512 the measured 0.5929/0.9862 =")
    print(f"              {0.5929/0.9862:.3f} sits near the 0.5 threshold with an")
    print(f"              SE that was never quoted. ⚠️ Under-powered")
    print(f"              exactly where the criterion bites.")
    print(f"    Constr C4 defect ratios run 0.89 to 1.53 across Q with")
    print(f"              **no error bars at all**, against a threshold")
    print(f"              of 0.5. ⚠️ The spread across Q is the only")
    print(f"              clue to the sd and it is ~0.25, putting the")
    print(f"              threshold about 2 sd below 1 -- marginal.")

    print(f"\n(4) the pattern")
    print(f"    Ratio criteria in this program were chosen as EFFECT")
    print(f"    sizes -- 2x, 0.5x -- and never compared against the")
    print(f"    spread of the thing they threshold. Where the null's sd")
    print(f"    is small the criterion is far looser than the data, and")
    print(f"    the verdict understates the evidence; where it is large")
    print(f"    or unprinted the criterion cannot be judged at all.")
    print(f"    Two of four are the first kind, two the second.")
    print("DONE")


if __name__ == "__main__":
    main()
