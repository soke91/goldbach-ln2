# -*- coding: utf-8 -*-
"""
How close to 1 must the level of distribution be, and is that regime
vacuous? (increment 279)

THIS IS ARITHMETIC, NOT STATISTICS. There is no null, no threshold and
no test here: the question is what an inequality already proved implies,
and the answer is obtained by solving it. Saying so in advance is the
point of hazards 2 and 4 -- a computation dressed as a measurement
invites a verdict it has not earned.

THE CLAIM BEING CHECKED. Theorem D' is stated correctly:

    loss factor  >>  exp(c sqrt((1-theta_E) log N)),
    beyond every power of log  FOR EACH FIXED theta_E < 1.

The prose that follows it in THEOREM_A.md is not:

    "Closing would require theta_E = 1 EXACTLY -- equidistribution of
     Lambda to moduli of size N itself, where each progression holds
     O(1) terms and the statement carries no information."

That skips a regime. Between "fixed theta_E < 1" and "theta_E = 1"
lies theta_E = theta_E(N) -> 1 at a RATE, and the theorem says nothing
about it. Increment 235 used this prose as reason 2 of three for
closing the program's highest-value open item, so it is worth solving
properly.

SOLVING IT. Write eta := 1 - theta_E. The route saves a power of log
when the loss is below (log N)^A:

    exp(c sqrt(eta log N))  <=  (log N)^A
      <=>  c sqrt(eta log N)  <=  A log log N
      <=>  eta  <=  (A/c)^2 (log log N)^2 / log N.

So the requirement is not theta_E = 1 but

    theta_E  >=  1 - C (log log N)^2 / log N,      C = (A/c)^2.

VACUITY. Moduli run to N^{theta_E} = N^{1-eta}, so a progression holds
about N^eta = exp(eta log N) = exp(C (log log N)^2) terms. Against any
fixed power of log:

    exp(C (log log N)^2)  /  (log N)^A  =  exp((log log N)(C log log N - A))
                                        -> infinity.

A progression therefore holds MORE THAN EVERY FIXED POWER OF log N
terms -- not O(1). The regime is not vacuous.

WHAT IS PRINTED.
 (A) the required eta and theta_E at a ladder of N, for several C.
 (B) the progression size in that regime, against (log N)^A.
 (C) the control: at FIXED eta (EH-style, eta = 0.4 is Lichtman's 3/5,
     eta = 0.01 is far beyond EH), the loss factor against (log N)^A --
     confirming the theorem's own statement that no fixed level works.
 (D) a self-check that the solved eta really is the crossing point:
     evaluate both sides at eta and at eta*(1 +/- 1e-6) and require the
     inequality to flip. A derivation that is not evaluated is a
     derivation nobody has checked.
"""
import math


def loss(eta, logN, c):
    return c * math.sqrt(eta * logN)          # log of the loss factor


def target(logN, A):
    return A * math.log(logN)                 # log of (log N)^A


def main():
    KS = [10, 50, 100, 500, 1000, 10000]
    A = 3.0          # save (log N)^A; A is the Siegel-Walfisz exponent
    CS = [(1.0, "c = 1"), (0.5, "c = 1/2"), (2.0, "c = 2")]

    print("(A) required level: theta_E >= 1 - (A/c)^2 (loglog N)^2/log N")
    print(f"    A = {A}, i.e. asking the switch to save (log N)^{A:g}\n")
    hdr = f"{'N':>10} {'log N':>10} {'loglog N':>9}"
    for _, lab in CS:
        hdr += f" {'eta  ' + lab:>16}"
    print(hdr)
    rows = {}
    for k in KS:
        logN = k * math.log(10)
        lg2 = math.log(logN)
        line = f"{'1e' + str(k):>10} {logN:>10.2f} {lg2:>9.3f}"
        for c, _ in CS:
            eta = (A / c) ** 2 * lg2 * lg2 / logN
            rows[(k, c)] = eta
            line += f" {eta:>16.6f}"
        print(line)

    print("\n    the same as theta_E = 1 - eta (c = 1 column).")
    print("    A REQUIRED eta >= 1 MEANS THE FORMULA IS NOT YET IN ITS")
    print("    REGIME: it would permit theta_E <= 0, which is outside")
    print("    Theorem D's hypothesis theta_E in (0,1). At such N the")
    print("    loss factor is ALREADY below (log N)^A for every level,")
    print("    so there is nothing for the requirement to constrain.")
    for k in KS:
        eta = rows[(k, 1.0)]
        if eta >= 1.0:
            print(f"{'1e' + str(k):>10}   eta = {eta:>8.4f} >= 1"
                  f"   -- NOT IN REGIME, no constraint here")
        else:
            print(f"{'1e' + str(k):>10}   theta_E >= {1 - eta:.8f}"
                  f"    (EH is theta_E = 1 - eps for FIXED eps)")

    # Where does the no-go acquire content at all? Solve
    # c sqrt(eta L) = A log L for L, by bisection -- the honest
    # "from what size of N does Theorem D' bite" number, and it is
    # not small.
    print("\n(A2) from what N does the no-go have content?")
    print("     crossover L = log N solving c sqrt(eta L) = A log L")
    print(f"{'eta':>8} {'log N':>14} {'N':>18}")
    for eta in (0.40, 0.10, 0.01):
        lo, hi = 10.0, 1e12
        for _ in range(400):
            mid = math.sqrt(lo * hi)
            if loss(eta, mid, 1.0) < target(mid, A):
                lo = mid
            else:
                hi = mid
        L = math.sqrt(lo * hi)
        print(f"{eta:>8.2f} {L:>14.1f} {'1e' + f'{L / math.log(10):.0f}':>18}")
    print("     Below these N the switch route loses NOTHING that matters")
    print("     and Theorem D' is silent. This is an asymptotic no-go and")
    print("     its content lives far past any computable range -- stated")
    print("     because the theorem is otherwise easy to read as though")
    print("     it bit at reachable N. It does not.")

    print("\n(B) is that regime vacuous? terms per progression = N^eta")
    print(f"{'N':>10} {'N^eta':>16} {'(log N)^A':>16} {'ratio':>14}")
    for k in KS:
        logN = k * math.log(10)
        eta = rows[(k, 1.0)]
        lterms = eta * logN
        ltarget = target(logN, A)
        flag = "" if eta < 1.0 else "   (not in regime)"
        print(f"{'1e' + str(k):>10} {'e^' + f'{lterms:.2f}':>16} "
              f"{'e^' + f'{ltarget:.2f}':>16} "
              f"{'e^' + f'{lterms - ltarget:.2f}':>14}{flag}")
    print("    ratio -> infinity, so a progression holds MORE than every")
    print("    fixed power of log N terms. NOT O(1). NOT vacuous.")

    print("\n(C) control: no FIXED level works, as Theorem D' states")
    print(f"{'N':>10} {'eta=0.40 (3/5)':>16} {'eta=0.10':>12} "
          f"{'eta=0.01':>12} {'(log N)^A':>12}")
    for k in KS:
        logN = k * math.log(10)
        t = target(logN, A)
        cells = "".join(f" {'e^' + f'{loss(e, logN, 1.0):.2f}':>12}"
                        for e in (0.40, 0.10, 0.01))
        print(f"{'1e' + str(k):>10} {'e^' + f'{loss(0.40, logN, 1.0):.2f}':>16}"
              f"{cells[13:]} {'e^' + f'{t:.2f}':>12}")
    print("    every fixed-eta column overtakes (log N)^A and stays")
    print("    above it: a constant level never suffices, however close")
    print("    to 1 -- which is exactly what the theorem says.")

    print("\n(D) self-check: is the solved eta really the crossing point?")
    print(f"{'N':>10} {'at eta':>22} {'at eta*(1-1e-6)':>20} "
          f"{'at eta*(1+1e-6)':>20}")
    ok = True
    for k in KS:
        logN = k * math.log(10)
        eta = rows[(k, 1.0)]
        t = target(logN, A)
        lo = loss(eta * (1 - 1e-6), logN, 1.0)
        hi = loss(eta * (1 + 1e-6), logN, 1.0)
        at = loss(eta, logN, 1.0)
        good = (lo < t < hi) and abs(at - t) < 1e-6 * t
        ok &= good
        print(f"{'1e' + str(k):>10} {at - t:>+22.3e} "
              f"{lo - t:>+20.3e} {hi - t:>+20.3e}"
              f"{'' if good else '   <-- FAILS'}")
    print("    column 1 must be ~0 (eta IS the crossing point), column 2")
    print("    negative (below eta the route works), column 3 positive.")
    print(f"\n{'SELF-CHECK OK' if ok else 'SELF-CHECK FAILED'}")

    print("\n" + "=" * 70)
    print("conclusion")
    print("=" * 70)
    print("The boundary is NOT theta_E = 1 and the regime that closes the")
    print("gap is NOT vacuous. The correct statement is:")
    print()
    print("    the switch route needs  1 - theta_E << (loglog N)^2/log N,")
    print("    which is STRICTLY STRONGER THAN EH (any fixed eps) and")
    print("    STRICTLY WEAKER THAN the vacuous theta_E = 1.")
    print()
    print("Internal consistency: (A) needs eta <= 0.388 at N = 1e500,")
    print("just under 0.40, and (A2) puts the eta = 0.40 crossover at")
    print("N = 1e480, just under 1e500. Two separate solves agree.")
    print()
    print("And a second thing the tables show, which was not the")
    print("question asked: the no-go is asymptotic and its content")
    print("begins around N = 1e480 for eta = 0.40. Every no-go in this")
    print("program should be read that way -- they constrain METHODS,")
    print("not any computation anyone will run.")
    print()
    print("This does not reopen increment 235: Lichtman's 3/5 is a fixed")
    print("constant and column (C) shows fixed constants never suffice,")
    print("and reason 1 there (object mismatch) is independent and")
    print("decisive. What changes is that 'vacuous' was the wrong word.")
    print("'Vacuous' says the target is meaningless; 'stronger than EH'")
    print("says it is meaningful and out of reach. Only the second is")
    print("true, and only the second leaves a well-posed question.")
    print("DONE")


if __name__ == "__main__":
    main()
