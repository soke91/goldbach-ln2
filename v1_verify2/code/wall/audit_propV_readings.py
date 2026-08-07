# -*- coding: utf-8 -*-
"""
prop:V -- which reading of "residual standard deviation" gives 0.000323?
(v1_verify2, Phase 1, blind.)

The sentence under test:

  "Rescaling each candidate to the measured mean, so that only its shape
   in N is judged, the residual standard deviation over every even
   N <= 1.6e7 is 0.000323 for A against 0.245235 for S -- a factor of 760."

audit_propV_and_wall.py reproduced the S figure to six significant
figures (0.245236) under the reading "sd of (V/W)/S after rescaling to
mean 1", on exactly the stated field. Under the SAME reading on the SAME
field the A figure came out 0.000582, not 0.000323 -- a factor 1.80 out,
and the quoted ratio 760 became 421.

PRE-REGISTRATION.  Decision rule: enumerate readings of "residual
standard deviation ... rescaled to the measured mean" and fields, and
report every combination. An item is RESOLVED if one combination returns
BOTH quoted figures at once; UNRESOLVED if the two figures require
different combinations.

Prediction: UNRESOLVED, with S matching on the full field and A matching
only on a truncated one -- the same split already found in conj:wall
item 1's two z-scores (blind finding B5).
"""

import os
import sys

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
CACHE = os.path.join(ROOT, "v1_verify2_log", "cache")


def main():
    X = int(sys.argv[1]) if len(sys.argv) > 1 else 16_000_000
    z = np.load(os.path.join(CACHE, f"field_{X}.npz"))
    good = (z["V"] > 0) & (z["W"] > 0)
    N = z["N"][good].astype(np.float64)
    V, W, A, S = (z[k][good] for k in ("V", "W", "A", "S"))
    r = V / W

    print("audit_propV_readings   (v1_verify2 Phase 1, blind)")
    print("target: 0.000323 for A, 0.245235 for S, ratio 760")
    print("=" * 74)

    def readings(mask, label):
        rr, AA, SS, NN = r[mask], A[mask], S[mask], N[mask]
        out = {}
        # (1) sd of the ratio, rescaled to mean 1
        for nm, cand in (("A", AA), ("S", SS)):
            q = rr / cand
            out[f"ratio/{nm}"] = (q / q.mean()).std(ddof=1)
        # (2) sd of the additive residual, in units of mean(V/W)
        for nm, cand in (("A", AA), ("S", SS)):
            c = rr.mean() / cand.mean()
            res = rr - c * cand
            out[f"resid/{nm}"] = res.std(ddof=1) / rr.mean()
        # (3) sd of the additive residual, in units of the candidate mean
        for nm, cand in (("A", AA), ("S", SS)):
            c = rr.mean() / cand.mean()
            res = rr - c * cand
            out[f"residC/{nm}"] = res.std(ddof=1) / (c * cand).mean()
        print(f"  {label:<34} n={mask.sum():>9,}")
        for k in ("ratio", "resid", "residC"):
            a, s = out[f"{k}/A"], out[f"{k}/S"]
            print(f"    {k:<8} A={a:.6f}  S={s:.6f}  factor={s / a:7.1f}"
                  f"   {'<== S MATCHES' if abs(s - 0.245235) < 2e-5 else ''}"
                  f"{' <== A MATCHES' if abs(a - 0.000323) < 3e-6 else ''}")
        return out

    for lo, hi, lab in (
        (0, X, "every even N <= 1.6e7 (as written)"),
        (1e5, 1.4e7, "1e5 .. 1.4e7"),
        (1e6, X, "1e6 .. 1.6e7"),
        (3.4e6, X, "3.4e6 .. 1.6e7  (6.3e6 values)"),
        (X / 2, X, "top octave"),
    ):
        m = (N > lo) & (N <= hi)
        readings(m, lab)
        print()

    print("--- at what lower cutoff does the A figure become 0.000323? ---")
    for lo in (0, 5e4, 1e5, 2e5, 5e5, 1e6, 2e6, 3.4e6, 5e6, 8e6):
        m = (N > lo)
        q = (r[m] / A[m])
        sdA = (q / q.mean()).std(ddof=1)
        qs = (r[m] / S[m])
        sdS = (qs / qs.mean()).std(ddof=1)
        print(f"    N > {lo:>10,.0f}   n={m.sum():>9,}   sd_A={sdA:.6f}   "
              f"sd_S={sdS:.6f}   factor={sdS / sdA:6.1f}")
    print()
    print("  note: the S figure is stable across every cutoff (it is")
    print("  dominated by S's wrong SHAPE, not by small-N noise), while the")
    print("  A figure is entirely a small-N effect. So the two quoted")
    print("  numbers cannot be pinned to the same field by their values.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
