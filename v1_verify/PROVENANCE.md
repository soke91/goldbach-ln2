# Provenance — every check in this tree, and what backs it

One row per check. `Code` is the script; `Results` is the file its
output was read from, under `v1_verify/results/<same subdirectory>/`.
Statement numbers are `wall_v1.tex`'s, resolved by
`code/verify/lint_numbering.py` and **not** copied from
`v1/PROVENANCE.md`, which is off by one from `conj:L` onward.

One-shot gate: `python v1_verify/code/verify/verify_all.py` — twelve
rows, each against a pre-registered interval, nonzero exit on failure.

## Reproduction stamps and lints (`code/verify/`)

| Check | What it does | Code | Results |
|---|---|---|---|
| the stamp | 12 rows, every finding's load-bearing number at reduced size | `verify_all.py` | `verify_all.txt` |
| the corrected paper | environments, braces, `\ref` resolution, astral characters, refuted forms, `v1` still read-only | `lint_corrected_paper.py` | `lint_corrected_paper.txt` |
| statement numbers | resolves the shared counter from source; flags every citation that disagrees | `lint_numbering.py` | `lint_numbering.txt` |
| renumbering | one-shot fix of the numbers this tree inherited | `_renumber.py` | — (applied, not a measurement) |

## The wall (`code/wall/`)

| Statement | Claim under test | Code | Results |
|---|---|---|---|
| Lemma 13 (`lem:MP`) | `sum_{N<=X} C^2 = sum_h M(h)P(h)` | `audit_lemMP_identity.py` | `audit_lemMP_identity.txt` |
| Proposition 15 (`prop:W`) | "under that input `rho -> 1`" | `audit_propW_chowla_gap.py` | `audit_propW_chowla_gap.txt` |
| — the reconstruction | `-0.0976`; the shift-mass table | `audit_propW_reconstruction.py` | `audit_propW_reconstruction.txt` |
| Proposition 20 (`prop:coh`) | the derivation, term by term | `audit_propcoh_cancellation.py` | `audit_propcoh_cancellation.txt` |
| `conj:wall` item 4 | the surrogate null | `audit_zeta_regression_null.py` | `audit_zeta_regression_null.txt` |
| — attribution, local | ordinates against their own neighbourhood | `audit_zeta_local_background.py` | `audit_zeta_local_background.txt` |
| — attribution, coin | Lemma 17 applied to the lines | `audit_zeta_coin_local.py` | `audit_zeta_coin_local.txt` |
| Lemma 17 (`lem:coin`) | applied to `rho` and to the zeta `R^2` | `audit_coin_control_v1claims.py` | `audit_coin_control_v1claims.txt` |
| `conj:wall` item 1 | "removing cell means alone" | `audit_conjwall_scalemask.py` | `audit_conjwall_scalemask.txt` |
| `conj:wall` items 1, 3 | the error bars, and the `t=5` tail | `audit_gaussian_errorbars.py` | `audit_gaussian_errorbars.txt` |

## The demand side (`code/demand/`)

| Statement | Claim under test | Code | Results |
|---|---|---|---|
| Theorem 1 (`thm:A`) | `sup_t |T_1(t)|`, computed exactly for every `t` | `audit_thmA_direct_sup.py` | `audit_thmA_direct_sup.txt` |
| Proposition 7 (`prop:E`) | the margin table, and whether it is grid-converged | `audit_propE_grid.py` | `audit_propE_grid.txt` |
| Proposition 8 (`prop:Dpp`) | the `CP_2` table and the tuned column | `audit_propDpp_table.py` | `audit_propDpp_table.txt` |

## The supply side (`code/supply/`)

| Statement | Claim under test | Code | Results |
|---|---|---|---|
| E1, §3.1 | the normalisation, and the margin at accessible `N` | `audit_E1_norm_and_gauss.py` | `audit_E1_norm_and_gauss.txt` |
| Conjecture 10 (`conj:L`) | band ratios, kurtosis, half-normal | `audit_E1_norm_and_gauss.py` | `audit_E1_norm_and_gauss.txt` |
| R4, §`sec:R4` | the switch identity, the block ratio and its null, the mirror | `audit_R4_switch.py` | `audit_R4_switch.txt` |
| C-III (2), §`sec:c3` | the Heath–Brown weight table, and the definition behind it | `audit_c3_hb_weight.py` | `audit_c3_hb_weight.txt` |
| K1, §7.2 | the orbit regression, with three nulls | `audit_K1_orbit.py` | `audit_K1_orbit.txt` |
| — its reach | how much of the 63-divisor orbit is live | `audit_K1_orbit_reach.py` | `audit_K1_orbit_reach.txt` |
| K3, §7.2 | the pair matrix, its trace moments, and the coin null | `audit_K3_wishart.py` | `audit_K3_wishart.txt` |
| R2, §7.2 | both criteria, their nulls' spreads, and the pair count | `audit_R2_criterion.py` | `audit_R2_criterion.txt` |
| C-I, §7.3 | the rational-peak excess at 8, 64 and 400 null draws | `audit_C1_nulldraws.py` | `audit_C1_nulldraws.txt` |
| R1, §7.2 | the chance floor, the six-draw null, and what the precision is a fraction of | `audit_R1_power.py` | `audit_R1_power.txt` |

## Not covered by this tree

Stated so the gaps are visible rather than implied:

- the five route adjudications of §7.1 — these are readings of the
  source papers' lemma hypotheses, not computations
- kill-test R3 (analytic), and the representation classes C-II and
  C-IV, whose known defects are already stated in the paper's own
  Remark on what a null verdict costs.
- three of the five route adjudications of §7.1; the MRT/Lichtman and
  Tao 2016 rows were checked against the sources and hold.
  K2 and K4 were checked at the level of their thresholds and counts
  (both calibrated) but not re-implemented.
- the reproduction stamps' own pre-registered intervals in
  `v1/code/verify/`

## Conventions

Every script carries its pre-registration in its docstring, including
the prediction it was written to test, so that a prediction that failed
cannot be reported afterwards as a surprise. Four of them did fail and
are recorded in `paper/ADVERSARIAL_FINDINGS.md` under "Refuted
objections".
