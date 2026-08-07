# Provenance — every statement in `paper/wall_v1.tex`, and what backs it

One row per numbered statement. `Code` is the script that produces the
figure; `Results` is the file the figure was read from, under
`v1/results/<same subdirectory>/`. Where a statement is analytic and has
no computation behind it, the row says so rather than pointing at
something adjacent.

A figure quoted in a document and absent from the result file that
produced it is a defect this program has made and now checks for
mechanically: `v1_log/code/audit_quoted_numbers.py`.

## §2 The demand side

| Statement | Claim | Code (`v1/code/demand/`) | Results |
|---|---|---|---|
| Theorem 1 (`thm:A`) | `sum_k mu(k) E_mu(t;k) << N(log N)^-A` | `thmA_E3.py`, `thmA_audit.py`, `thmA_fix.py`, `thmA_scale.py` | `thmA_fix_err.txt`, `thmA_fix_result.txt`, `thmA_scale_err.txt` (+1 more) |
| — the density identity | exact rational arithmetic, all 243 squarefree `m<400` | `thmA_density_check.py` | — **no saved output** |
| — the `log k` branch | where the `(q,N)>1` trap lives | `thmA_logw.py`, `thmA_mtlog.py` | — **no saved output** |
| Corollary 2 (`cor:B`) | `E_4 << N(log N)^-A`, via `E_4 = int T_1(t)dt/(N-t)` | `complement_identity.py` | — **no saved output** |
| Theorem 3 (`thm:C`) | `E_3` equivalent to binary Goldbach; Huang–Li (22) | `thmC_alpha_scan.py`, `hl_assembly.py`, `hl_rederive.py`, `hl_S1_check.py` | `hl_assembly_err.txt`, `hl_assembly_result.txt`, `hl_rederive_err.txt` (+1 more) |
| Theorem 5 (`thm:D`) | no weight extracts `C(N)`: `‖b‖_1/|B_w| >> exp(c sqrt(log N/2))` | `thmD_tradeoff.py` | `thmD_tradeoff.txt` |
| Theorem 6 (`thm:Dprime`) | the no-go survives `EH`; closing needs `theta_E = 1` | `thmD_tradeoff.py` | — **no saved output** |
| Proposition 7 (`prop:E`) | the circle method has zero margin; `0.168, 0.175, 0.158, 0.152` | `circle_margin.py`, `reaudit_prop_E.py` | `circle_margin_err.txt`, `circle_margin_result.txt`, `reaudit_prop_E.txt` |
| Proposition 8 (`prop:Dpp`) | polynomial weights; `CP_2` ratio `0.771, 0.790, 0.810` | `thmD2_polyweight.py` | — **no saved output** |
| §2.4 the (18) defect | the omitted term `Delta`, and its repair | analytic; `residual_identity.py`, `vaughan_pieces.py` support it | — **no saved output** |

## §3 The supply side

| Statement | Claim | Code (`v1/code/supply/`) | Results |
|---|---|---|---|
| E1 (weak form) | normalisation is `M_k^2`, not `M_k`; `0.305, 0.310, 0.319` at `1e8` | `e1_corrected_norm.py`, `e1_target_audit.py`, `e1_target_1e9.py` | `e1_target_audit_err.txt`, `e1_target_audit_result.txt` |
| Remark 10 | the mandatory target-versus-measurement comparison | `e1_target_multiN.py` | — **no saved output** |
| Conjecture 11 (`conj:L`) | field = mask × exactly Gaussian fluctuation | `e1_mask_model.py`, `e1_local_model.py`, `e1_exact_cells.py`, `e1_blind_1e9.py` | — **no saved output** |
| — pair spectrum | `lambda_max` at the Wishart null, `z = -0.19` | `e1_wishart_null.py`, `e1_pair_local.py` | — **no saved output** |
| — E1 at `1e9` | band ratios `0.966 / 0.950 / 0.922` | `e1_law_1e9.py`, `e1_law_1e9b.py`, `e1_cells_1e9.py` | — **no saved output** |
| — definitive grid | 4 `N` × 2 bands × 300 `k`, mean `z = +0.02` | `e1_grid_final.py`, `e1_settle_800.py`, `e1_settle_900.py` | — **no saved output** |

## §7 The negative map

| Statement | Claim | Code (`v1/code/supply/`) | Results |
|---|---|---|---|
| K1 | multiplicative Fejér kernel dead, `R^2 = 0.466` | `e1_forge_kt1.py` | — **no saved output** |
| K2 | manufactured pair congruence dead, 0/10 | `e1_forge_kt2.py` | — **no saved output** |
| — its floor | `4/sqrt(n) ≈ 0.163` sd; needs `14,000×` at `A=1` | `reaudit_killtest_power.py`, `reaudit_killtest_economics.py` | `reaudit_killtest_economics.txt`, `reaudit_killtest_power.txt` |
| K3 | Wishart moments dead-centre, `z = -0.48/+0.42/-0.03` | `e1_forge_kt3.py` | — **no saved output** |
| K4 | `N`-average descent dead, 0/8 | `e1_forge_kt4.py` | — **no saved output** |
| R1 | zero spectrum `0.2152` vs `0.2196 ± 0.0055` | `e1_forge_r1.py`, `reaudit_killtest_thresholds.py` | `reaudit_killtest_thresholds.txt` |
| R2 | determinant phase blind, `R^2 = -0.0001/+0.0004` | `e1_forge_r2.py`, `e1_forge_r2b.py` | — **no saved output** |
| R3 | character transform dead by analysis | analytic; `e1_thin_closure.py` measures the thin-progression cost | — **no saved output** |
| R4 | divisor switch does not localize | `e1_forge_r4.py`, `e1_forge_r4b.py` | — **no saved output** |
| — its error bars | 5σ at `B=8`, no information at `B=512` | `reaudit_r4_errorbars.py` | `reaudit_r4_errorbars.txt` |
| R5 | circle method, zero margin | see Proposition 7 | — |
| C-I | abelian spectrum mask-exact, 0/6 | `e1_constr_c1.py` | — **no saved output** |
| C-II | fired at `z = +10.97`, collapsed under 64-draw nulls | `e1_constr_c2.py`, `e1_constr_c2b.py` | — **no saved output** |
| — its error bars | the `0.5×` threshold sits `-1.29` sd below its null | `reaudit_c4_errorbars.py` | `reaudit_c4_errorbars.txt` |
| C-IV | manufactured modularity, 0/6 levels | `e1_constr_c4.py` | — **no saved output** |
| C-III (1) | the pencil/parallel obstruction to any change of variable | `c3_pencil_check.py`, `c3_transform_probe.py` | — **no saved output** |
| C-III (2) | Heath–Brown weight: `0.939…0.969` outside the covered region | `c3_hb_mass.py` | `c3_hb_mass_1e6.txt`, `c3_hb_mass_err.txt`, `c3_hb_mass_result.txt` |
| — kill-test nulls | which kill-tests carry a coin control | `audit_killtest_nulls.py` | `audit_killtest_nulls.txt` |

## §4–§6 The wall

| Statement | Claim | Code (`v1/code/wall/`) | Results |
|---|---|---|---|
| Proposition 12 (`prop:V`) | `V = W·A(N)(1+o(1))`; `A` not `S`, factor `760` | `lab_wall_variance_law.py`, `lab_variance_law_reaudit.py`, `lab_V_asymptotic.py` | `lab_V_asymptotic.txt`, `lab_variance_law_reaudit.txt`, `lab_wall_variance_law.txt` |
| — `V/W` vs `A(N)` | mean `1.000000`, sd `0.000145` top octave | `lab_variance_ratio.py` | `lab_variance_ratio.txt` |
| — unambiguous as a summary | ratio-of-means equals mean-of-ratios to `0.00%` | `audit_aggregates.py` | `audit_aggregates.txt` |
| Lemma 13 (`lem:MP`) | `sum_N C(N)^2 = sum_h M(h)P(h)` | `wall_secondmoment.py`, `attack_wall_identity.py` | `attack_wall_identity.txt`, `wall_secondmoment_err.txt`, `wall_secondmoment_result.txt` |
| Conjecture 14 (`conj:wall`) | `C = m(N) + sqrt(V)·G` | — (the four measurements below) | — |
| — bulk Gaussian | kurtosis `-0.0005` (`z=-0.3`); `+0.1704` at `z=98` under `S·N` | `lab_gaussian_half_audit.py` | `lab_gaussian_half_audit.txt` |
| — tail Gaussian | Gumbel deviation `+0.54 ± 0.45`; ratios `0.999/0.997/0.878` | `lab_wall_tails.py`, `lab_wall_tails_exact.py`, `lab_wall_extremes.py` | `lab_wall_extremes.txt`, `lab_wall_tails.txt`, `lab_wall_tails_exact.txt` |
| — zeta phases | `R^2 = 3.90e-3` vs surrogate max `5.09e-6` | `lab_E_zeta_spectrum.py`, `lab_wall_spectral_share.py` | `lab_E_zeta_spectrum.txt`, `lab_wall_spectral_share.txt` |
| Proposition 16 (`prop:W`) | `rho - 1 = (1/V) sum_h c(h) S(h)` | `lab_offdiag_chowla.py`, `lab_offdiagonal_reaudit.py` | `lab_offdiag_chowla.txt`, `lab_offdiagonal_reaudit.txt` |
| Remark 17 | `rho` names three quantities, differing `10.3%` at `1e5` | `audit_rho_definitions.py` | `audit_rho_definitions.txt` |
| — mass by shift | `1.1 / 3.0 / 23.1 / 48.9 / 23.8` percent | `lab_wall_process_mass.py` | `lab_wall_process_mass.txt` |
| — major arcs | `8.40×` at `q=3`, `15.16×` at `q=5` | `lab_atoms_perq.py`, `lab_atoms_normalisation.py`, `lab_atoms_localise.py` | `lab_atoms_localise.txt`, `lab_atoms_normalisation.txt`, `lab_atoms_perq.txt` |
| — the spectral process | atoms on `j/q` with weights `mu^2(q)/phi^2(q)` | `lab_wall_process.py`, `lab_mask_vs_atoms.py` | `lab_mask_vs_atoms.txt`, `lab_wall_process.txt`, `lab_wall_process_mass.txt` |
| Lemma 18 (`lem:coin`) | `eps^2 = mu^2`, so `V` is unchanged | `lab_rho_coin_control.py`, `lab_mask_coin_control.py`, `lab_spectrum_coin_control.py` | `lab_mask_coin_control.txt`, `lab_rho_coin_control.txt`, `lab_spectrum_coin_control.txt` |
| — what it invalidates | the de-trended `rho` fit and its rate | `lab_rho_decay.py`, `lab_rho_rate.py`, `lab_rho_uncentred.py` | `lab_rho_decay.txt`, `lab_rho_rate.txt`, `lab_rho_uncentred.txt` |
| Lemma 19 (`lem:placebo`) | permuting cell labels leaves `Z` byte-identical | `lab_cell_floor_mechanism.py` | `lab_cell_floor_mechanism.txt` |
| Lemma 20 (`lem:cellmom`) | `Var(m_c - m) = Q_cc/n_c^2 - 2Q_ca/(n_c n) + Q_aa/n^2` | `lab_mask_exact_floor.py`, `lab_cell_floor_exact.py`, `audit_cellvar_bias.py` | `audit_cellvar_bias.txt`, `lab_cell_floor_exact.txt`, `lab_mask_exact_floor.txt` |
| Proposition 21 (`prop:coh`) | `se ∝ (log N)^-1/2`; measured `b = 0.0379` vs predicted `0.036` | `lab_mask_significance.py` | `lab_mask_significance.txt` |
| §6.3 shape | corr with the singular series of the shift, `0.9997–1.0000` | `lab_cell_floor_shift.py`, `lab_cell_floor_hresolved.py` | `lab_cell_floor_hresolved.txt`, `lab_cell_floor_shift.txt` |
| §6.3 size | accounted to `~6%`; the cap inflates `1.19 → 1.69` | `lab_cellfloor_weighting.py` | `lab_cellfloor_weighting.txt` |
| §6.3 rarity | per-cell share `0.94` deepest vs `0.018` pooled | `lab_mask_share_percell.py` | `lab_mask_share_percell.txt` |
| §6.4 exponents | the six `a_d` with s.e.; `chi^2/dof = 251` | `lab_mask_exponent_se.py` | `lab_mask_exponent_se.txt` |
| Proposition 22 (`prop:scaleinv`) | `D_c` scale-invariant; predicted exponents `0.000 ± 0.0007` | `lab_mask_exponent_predict.py` | `lab_mask_exponent_predict.txt` |
| §6 decay, undetermined | `N^-a` vs `(log N)^-b` not separated over a factor `160` | `lab_mask_amplitude_law.py`, `lab_wall_exponent_reaudit.py` | `lab_mask_amplitude_law.txt`, `lab_wall_exponent_reaudit.txt` |
| §9 the margin | `max|C|/N` falls `0.056 → 0.0082`; margin `N^0.454` at `1e8` | `lab_margin_extrapolate.py`, `attack_wall_units.py` | `lab_margin_extrapolate.txt`, `attack_wall_units.txt` |
| §9 target units | trivial `(log N)^A` above, C–S `(log N)^{A+1/2}` above | `attack_wall_units.py`, `attack_uncertainty_slack.py` | `attack_uncertainty_slack.txt`, `attack_wall_units.txt` |

## Reproduction stamps

| Script (`v1/code/verify/`) | What it gates |
|---|---|
| `verify_all.py` | the core corpus, every stamp against a pre-registered interval, exits nonzero on failure. Intervals are widths measured from the rows' own spread (`centre ± 4σ`), not chosen. |
| `verify_deep.py` | the deep-`N` arm: 300 deep (`30030 \| N`) versus 300 shallow at `N ≈ 1e8`, gating **the gap** `+0.2238 ± 0.0056`. |
| `verify_propositions.py` | the propositions' identities. |

## What is not here

`v1_log/` holds the process record — the correction log, the open-question
register, the notation register, and the exploratory code and results the
paper does not cite. It is not distributed. Nothing in `v1/` reads from
it, so this directory is complete on its own.

## The gap this table exposes

Regenerating the Results column from disk, rather than writing it from
the script names, showed that **a substantial number of the cited
computations have code but no archived output**: they were run with
their output read from the console and never redirected to a file. The
rows above say `no saved output` where that is the case, and they say it
rather than naming a plausible file that does not exist.

This does not put any number in the paper in doubt — the figures were
read from the runs when they happened, and the scripts are
deterministic — but it does mean those rows are reproducible rather than
*reproduced*. Closing the gap is mechanical:

```
python <script> > v1/results/<subdir>/<script stem>.txt 2>&1
```

Two naming conventions are present in the archived files, `<stem>.txt`
and the older `<stem>_result.txt` / `<stem>_err.txt` pair. New runs
should use the first.
