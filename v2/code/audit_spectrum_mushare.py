# -*- coding: utf-8 -*-
r"""
Is there a mu-share in the wall's spectrum beyond the principal arcs?

WHAT IS AT STAKE

OPEN.md's wall item 6 has stood untouched: the spectral measure is
atomic at the rationals j/q with weights mu^2(q)/phi^2(q), and that
is Lambda's structure.  **mu's own contribution has been measured
only through the principal-arc deficit** -- v1 recorded 8.399 at
q = 3 and 15.163 at q = 5 -- and whether anything of mu lies beyond
that was never asked.

The table those two numbers come from already carries the answer's
ingredients.  v1/results/wall/lab_atoms_perq.txt lists, per modulus,
the share of the spectrum the real field puts there and the share the
coin puts there.  The coin replaces mu by random signs and keeps
everything else, so **the difference between the two shares is mu's
contribution by construction**, not by inference.

And the difference is not flat.  The prime moduli carry less real
share than coin share and the composite ones carry more, which is a
dependence on omega(q) -- the number of prime factors -- and omega
is not in mu^2(q)/phi^2(q).  If that dependence is real, mu's share
is structured and structured by something Lambda's weights do not
contain.

Nothing is measured here.  Every value is read from v1's table.

BACKS: Remark {#rem:spectrummushare} in paper/wall_v3.md.

PRE-REGISTERED PREDICTIONS (written before this script was run)

  K1  The control.  The two moduli OPEN.md names reproduce from the
      table: the exponential-sum ratio is 8.399 at q = 3 and 15.163
      at q = 5, to the three decimals printed.
  K2  **mu's share depends on omega(q)**: regressing the log of
      (real share / coin share) on omega gives a positive slope
      resolved at |t| > 2.
  K3  And omega is not standing in for size: adding log q to the
      regression leaves omega's coefficient positive and resolved.
  K4  The two moduli OPEN.md quotes are not special: both sit in the
      omega = 1 group, so what was measured there is the omega = 1
      case rather than a fact about 3 and 5.

REFUTATION RULE (fixed before the run)

  K1  REFUTED outside the printed decimals; then this is not the
      table those numbers came from.  THIS ONE GATES.
  K2  REFUTED if the slope is negative or unresolved.  **Unresolved
      is the likely failure and means only that sixteen moduli
      cannot see it** -- not that mu's share is unstructured; the
      printed table is what a reader should judge that on.  A
      resolved negative slope would be a different finding again and
      would say mu concentrates on the prime moduli.
  K3  REFUTED if omega loses its sign or its resolution once log q
      is in.  Then the dependence is on the size of the modulus and
      omega was a proxy, which is a weaker claim and not the one
      predicted -- and with sixteen points and two regressors,
      collinearity is the reason to expect it.
  K4  REFUTED if either modulus has omega other than 1.  Arithmetic,
      not measurement; it is here because the reading of K2 depends
      on it.

  K1 gates.  K2 to K4 are the measurement and do not gate.

  NO NULL IS RUN and none applies -- the coin arm is already in the
  table and is the comparison, not a background to detect against.
  v1 ran it in lab_spectrum_coin_control.py.
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
V1 = os.path.abspath(os.path.join(ROOT, "..", "v1", "results", "wall"))
OUT = os.path.join(RES, "audit_spectrum_mushare.txt")

DEC = 3
ROW = (r"^\s+(\d+)\s+(\d+)\s+([\d.]+)\s+([\d.]+)\s+(yes|NO)\s+"
       r"([\d.]+)\s+([\d.]+)\s*$")


def read_table():
    src = io.open(os.path.join(V1, "lab_atoms_perq.txt"),
                  encoding="utf-8").read()
    out = []
    for m in re.finditer(ROW, src, re.M):
        out.append((int(m.group(1)), int(m.group(2)),
                    float(m.group(3)), float(m.group(4)),
                    m.group(5), float(m.group(6)),
                    float(m.group(7))))
    return out


def omega(q):
    n, k = q, 0
    d = 2
    while d * d <= n:
        if n % d == 0:
            k += 1
            while n % d == 0:
                n //= d
        d += 1
    return k + (1 if n > 1 else 0)


def ols(X, y):
    A = np.column_stack([np.ones(len(y))] + X)
    c, *_ = np.linalg.lstsq(A, y, rcond=None)
    r = y - A.dot(c)
    n = len(y)
    s2 = float((r ** 2).sum()) / (n - A.shape[1])
    cov = s2 * np.linalg.inv(A.T.dot(A))
    return c, np.sqrt(np.diag(cov)), float(np.sqrt((r ** 2).mean()))


def main():
    lines = []

    def say(t=""):
        print(t)
        sys.stdout.flush()
        lines.append(t)

    tab = read_table()
    rnd = 0.5 * 10.0 ** (-DEC)
    say("read %d moduli from v1/results/wall/lab_atoms_perq.txt; "
        "nothing is measured here" % len(tab))
    say("PRINTBOUND audit_spectrum_mushare %d %.8f" % (DEC, rnd))

    # -------------------------------------------------------------- K1
    say()
    say("K1  the control on the two moduli OPEN.md names")
    got = {q: E for q, _f, _P, E, _a, _sr, _sc in tab}
    k1 = (abs(got.get(3, 0.0) - 8.399) <= rnd
          and abs(got.get(5, 0.0) - 15.163) <= rnd)
    say("  q = 3 gives %.3f against v1's published 8.399; q = 5 "
        "gives %.3f against its published 15.163"
        % (got.get(3, float("nan")),
                    got.get(5, float("nan"))))
    say("  K1 %s   (cap: the three decimals v1 prints)"
        % ("hold" if k1 else "REFUTED"))
    if not k1:
        io.open(OUT, "w", encoding="utf-8", newline="\n").write(
            "\n".join(lines) + "\n")
        raise SystemExit(1)

    say()
    say("  the shares, and what mu does to them")
    say("     q  omega   share real   share coin   real/coin")
    qs, om, rat = [], [], []
    for q, _f, _P, _E, _a, sr, sc in tab:
        w = omega(q)
        r = sr / sc
        qs.append(q)
        om.append(w)
        rat.append(r)
        say("  %5d  %-6d %-12.4f %-12.4f %.3f" % (q, w, sr, sc, r))
    say("SCALES 1")
    qs = np.array(qs, dtype=np.float64)
    om = np.array(om, dtype=np.float64)
    y = np.log(np.array(rat))

    # -------------------------------------------------------------- K2
    say()
    say("K2  does mu's share depend on omega(q)?")
    c, se, rms = ols([om], y)
    t = c[1] / se[1]
    k2 = c[1] > 0.0 and abs(t) > 2.0
    say("  log(real/coin) on omega: slope %+.4f +- %.4f, t = %.2f, "
        "r.m.s. %.4f" % (c[1], se[1], t, rms))
    say("TSTAT mushare_omega %.2f" % t)
    if abs(t) < 2.0:
        say("UNRESOLVED SIGN mushare_omega")
    say("SPREAD mushare_omega %.4f" % (om.max() - om.min()))
    say("SCATTER slope_audit_spectrum_mushare %.4f" % rms)
    say("  K2 %s   (cap: positive and |t| > 2)"
        % ("hold" if k2 else "REFUTED"))

    # -------------------------------------------------------------- K3
    say()
    say("K3  is omega standing in for the size of the modulus?")
    c2, se2, rms2 = ols([om, np.log(qs)], y)
    t2 = c2[1] / se2[1]
    tq = c2[2] / se2[2]
    k3 = c2[1] > 0.0 and abs(t2) > 2.0
    say("  with log q added: omega %+.4f +- %.4f (t = %.2f), log q "
        "%+.4f +- %.4f (t = %.2f)"
        % (c2[1], se2[1], t2, c2[2], se2[2], tq))
    say("TSTAT mushare_omega_adj %.2f" % t2)
    say("SPREAD mushare_omega_adj %.4f" % (om.max() - om.min()))
    if abs(t2) < 2.0:
        say("UNRESOLVED SIGN mushare_omega_adj")
    say("CORR mushare_regressors %.5f"
        % abs(float(np.corrcoef(om, np.log(qs))[0, 1])))
    say("  r.m.s. %.4f against %.4f without it" % (rms2, rms))
    say("  K3 %s   (cap: omega positive and resolved with log q in)"
        % ("hold" if k3 else "REFUTED"))

    # -------------------------------------------------------------- K4
    say()
    say("  NOT PRE-REGISTERED, reported because honesty needs it:")
    say("  q = 15 sits at %.3f, far above every other modulus, and a "
        "single" % max(rat))
    say("  point that size could carry a slope on its own. Leaving "
        "each modulus")
    say("  out in turn gives")
    ts = []
    for i in range(len(y)):
        m = np.ones(len(y), dtype=bool)
        m[i] = False
        ci, sei, _ = ols([om[m]], y[m])
        ts.append(ci[1] / sei[1])
    ts = np.array(ts)
    say("  slope t from %.2f to %.2f over the %d leave-one-out fits, "
        "and dropping" % (ts.min(), ts.max(), len(ts)))
    say("  q = 15 alone leaves t = %.2f"
        % ts[[q == 15 for q in qs].index(True)])
    say("TSTAT mushare_omega_loo_min %.2f" % ts.min())
    say("SPREAD mushare_omega_loo_min %.4f" % (om.max() - om.min()))
    say("  so the dependence is not one modulus' doing"
        if abs(ts).min() > 2.0 else
        "  so at least one modulus carries the resolution on its own")

    say()
    say("K4  are the two quoted moduli special?")
    k4 = omega(3) == 1 and omega(5) == 1
    grp = [r for q, w, r in zip(qs, om, rat) if w == 1]
    say("  omega(3) = %d and omega(5) = %d; the omega = 1 group has "
        "%d moduli with" % (omega(3), omega(5), len(grp)))
    say("  real/coin from %.3f to %.3f" % (min(grp), max(grp)))
    say("  K4 %s   (cap: both have omega = 1)"
        % ("hold" if k4 else "REFUTED"))

    say()
    say("what this settles")
    if k2 and k3:
        say("  mu's share is structured, and structured by omega, "
            "which mu^2(q)/phi^2(q)")
        say("  does not contain -- so there is something of mu in "
            "the spectrum beyond")
        say("  Lambda's weights and beyond the two principal arcs "
            "OPEN.md quotes")
    elif k2 and not k3:
        say("  the dependence is on the modulus's size and omega was "
            "a proxy; mu's share")
        say("  is still structured but by something the table cannot "
            "separate from q")
    else:
        say("  %d moduli cannot resolve the dependence; the "
            "table is what a reader" % len(tab))
        say("  should judge it on and nothing here says the share is "
            "unstructured")

    say()
    say("=" * 70)
    say("K1 %s  K2 %s  K3 %s  K4 %s"
        % tuple("hold" if v else "REFUTED"
                for v in (k1, k2, k3, k4)))

    head = [
        "STATISTIC: per modulus q, the ratio of the spectral share",
        "           the real field puts at j/q to the share the coin",
        "           puts there -- mu's contribution by construction,",
        "           since the coin replaces mu by random signs and",
        "           keeps the rest -- regressed on omega(q) alone and",
        "           on omega(q) with log q.",
        "NULL: none is run and none applies. The coin arm is already",
        "      in the table and is the comparison itself, not a",
        "      background to detect against; v1 ran it in",
        "      lab_spectrum_coin_control.py.",
        "FIELD: the %d moduli of" % len(tab),
        "       v1/results/wall/lab_atoms_perq.txt at n = 3873870,",
        "       with the periodogram and exponential-sum ratios and",
        "       the real and coin shares v1 recorded. Nothing is",
        "       measured here; omega(q) is computed by trial",
        "       division and the two quoted moduli are checked",
        "       against the 8.399 and 15.163 that OPEN.md's wall",
        "       item 6 carries.",
        "",
    ]
    io.open(OUT, "w", encoding="utf-8", newline="\n").write(
        "\n".join(head + lines) + "\n")
    print("\nwrote %s" % OUT)
    if not k1:
        raise SystemExit(1)
    return 0


if __name__ == "__main__":
    sys.exit(main())
