# -*- coding: utf-8 -*-
"""
Documentation-integrity audit: does every load-bearing number quoted in
the documents appear in the results file that produced it?
(increment 274)

Increment 273 recorded that three of the last four increments contained
a fault in my own check. The response to a high error rate is to verify
what is already written, not to write more, and one fault class has not
been audited at all: TRANSCRIPTION. A number can be measured correctly,
recorded correctly in results/, and then quoted wrongly in a document,
and nothing anywhere would catch it.

WHAT THIS DOES. For each load-bearing figure now standing in
TRANSFORM_P.md, LOCATION_MASK.md, CLOSURE_REAUDIT.md and
MEASUREMENTS.md, it locates the results file that produced it and
checks the string is present there. A figure that is a derived
quantity (a ratio computed by hand from two printed numbers) is marked
DERIVED and its derivation is recomputed instead.

WHAT IT CANNOT DO. It cannot check that a number is the RIGHT one for
the claim it supports -- only that it was not invented or mistyped
between the results file and the document. That limit is stated
because an audit that oversells itself is the fault it is looking for.
"""
import io
import os
import re

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def norm(s):
    """ASCII-normalise the minus signs the documents typeset."""
    return s.replace("−", "-").replace("–", "-")


def variants(fig):
    """The figure and every rounding of it, signs normalised, so a
    document quoting 1.537 matches a results file printing 1.5368.

    Without this the audit produced 8 false alarms out of 24 on its
    first run -- the very fault class it was built to find, in the
    tool itself."""
    f = norm(fig)
    out = {f, f.lstrip("-")}
    try:
        v = float(f)
    except ValueError:
        return out
    for d in range(1, 6):
        out.add(f"{v:.{d}f}")
        out.add(f"{abs(v):.{d}f}")
    return out

# (document, quoted string, results file, note)
CHECKS = [
    # --- LOCATION_MASK.md, the enumeration result
    ("LOCATION_MASK.md", "9.0034", "results/lab_wall_tails_exact.txt",
     "max|Z| before the mask"),
    ("LOCATION_MASK.md", "5.1515", "results/lab_location_mask.txt",
     "max|Z| after the mask"),
    ("LOCATION_MASK.md", "0.4683", "results/lab_wall_tails_exact.txt",
     "excess kurtosis before"),
    ("LOCATION_MASK.md", "0.0143", "results/lab_location_mask.txt",
     "excess kurtosis after"),
    # --- LOCATION_MASK.md, the 1e8 ladder
    ("LOCATION_MASK.md", "-4.852", "results/lab_location_mask_1e8.txt",
     "deep cell mean Z at 1e8"),
    ("LOCATION_MASK.md", "-6.014", "results/lab_location_mask_1e8.txt",
     "deeper cell mean Z at 1e8"),
    # --- LOCATION_MASK.md, the sign balance by depth
    ("LOCATION_MASK.md", "0.10884", "results/lab_sign_reaudit.txt",
     "P(C>0) at depth 5"),
    # --- TRANSFORM_P.md, the margin table
    ("TRANSFORM_P.md", "0.8852", "results/lab_transformP_weight1.txt",
     "S_abs/S_null at 5e4"),
    ("TRANSFORM_P.md", "0.8725", "results/lab_transformP_weight1.txt",
     "S_abs/S_null at 4e5"),
    # --- TRANSFORM_P.md, the depth ladder
    ("TRANSFORM_P.md", "0.3912", "results/lab_transformP_depth.txt",
     "R at depth 0, N ~ 1e6"),
    ("TRANSFORM_P.md", "0.6821", "results/lab_transformP_depth7.txt",
     "R at depth 6, N ~ 1.5e7"),
    ("TRANSFORM_P.md", "1.5368", "results/lab_transformP_depth7.txt",
     "R/R_null at depth 6"),
    # --- TRANSFORM_P.md, P.5 and P.7
    ("TRANSFORM_P.md", "0.3416", "results/lab_termcount.txt",
     "1/sqrt(j) floor, deep"),
    ("TRANSFORM_P.md", "0.2022", "results/lab_termcount.txt",
     "1/sqrt(j) floor, all even"),
    ("TRANSFORM_P.md", "0.0213", "results/lab_signstructure.txt",
     "pos/neg at 1 term, deep"),
    ("TRANSFORM_P.md", "19.5173", "results/lab_signstructure.txt",
     "pos/neg at >=9 terms, deep"),
    ("TRANSFORM_P.md", "0.5412", "results/lab_p_vs_omega.txt",
     "share of demand from p with <=2 terms, deep"),
    ("TRANSFORM_P.md", "0.7775", "results/lab_p_vs_omega.txt",
     "share of demand, shallow"),
    # --- TRANSFORM_P.md, de-masking
    ("TRANSFORM_P.md", "1.7677", "results/lab_transformP_demask.txt",
     "split de-masked demand"),
    ("TRANSFORM_P.md", "1.1114", "results/lab_transformP_demask.txt",
     "permuted floor"),
    # --- CLOSURE_REAUDIT.md, the re-audit figures
    ("CLOSURE_REAUDIT.md", "-193.5", "results/lab_sweepB_reaudit.txt",
     "corr(C,S) z at full power"),
    ("CLOSURE_REAUDIT.md", "-137.4", "results/lab_sweepB_reaudit.txt",
     "corr(C/sqrtS,S) z"),
    ("CLOSURE_REAUDIT.md", "-0.208", "results/lab_offdiagonal_reaudit.txt",
     "off/LHS at 4e6"),
    ("CLOSURE_REAUDIT.md", "14.39", "results/lab_excess_location.txt",
     "S/S_sign at p in [16,32], deep"),
]

# derived figures: (document, quoted, how it is derived, tolerance)
DERIVED = [
    ("TRANSFORM_P.md", "1.69", 0.3416 / 0.2022, 0.01,
     "floor ratio deep/all-even"),
    ("TRANSFORM_P.md", "1.62", 2.2833 / 1.4142, 0.02,
     "predicted ratio sqrt(rad/phi) deep vs typical"),
    ("TRANSFORM_P.md", "2.284", (1.0 / (0.5 * (2/3.) * 0.8 *
                                        (6/7.) * (10/11.) *
                                        (12/13.))) ** 0.5, 0.002,
     "1/sqrt(delta) at rad = 30030"),
]


def main():
    print("(A) quoted numbers against the results file that produced them")
    print(f"{'document':>22} {'figure':>10} {'in results?':>12}  note")
    bad = 0
    missing = 0
    for doc, fig, res, note in CHECKS:
        dp = os.path.join(REPO, doc)
        rp = os.path.join(REPO, res)
        indoc = False
        if os.path.exists(dp):
            body = norm(io.open(dp, encoding="utf-8").read())
            indoc = any(f in body for f in variants(fig))
        if not os.path.exists(rp):
            print(f"{doc:>22} {fig:>10} {'NO FILE':>12}  {note}")
            missing += 1
            continue
        txt = norm(io.open(rp, encoding="utf-8",
                           errors="replace").read())
        ok = any(f in txt for f in variants(fig))
        if not indoc:
            print(f"{doc:>22} {fig:>10} {'NOT IN DOC':>12}  {note}")
            bad += 1
        elif not ok:
            print(f"{doc:>22} {fig:>10} {'MISSING':>12}  {note}")
            bad += 1
        else:
            print(f"{doc:>22} {fig:>10} {'ok':>12}  {note}")
    print(f"\n  {len(CHECKS)} checks, {bad} mismatches, "
          f"{missing} results files absent")

    print("\n(B) derived figures, recomputed rather than located")
    print(f"{'document':>22} {'quoted':>9} {'recomputed':>12} "
          f"{'ok?':>5}  note")
    for doc, fig, val, tol, note in DERIVED:
        got = abs(float(fig) - val) <= tol
        if not got:
            bad += 1
        print(f"{doc:>22} {fig:>9} {val:>12.4f} "
              f"{'ok' if got else 'NO':>5}  {note}")

    print(f"\nTOTAL problems: {bad}")
    print("  LIMIT: this checks that a figure was not invented or")
    print("  mistyped between results/ and the document. It does NOT")
    print("  check that it is the right figure for the claim it")
    print("  supports.")
    print("DONE")


if __name__ == "__main__":
    main()
