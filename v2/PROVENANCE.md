# Provenance — every numbered statement in `paper/`, and what backs it

One row per numbered statement. `Code` is the script that produces the figure;
`Results` is the file it was read from. Analytic statements say so rather than
pointing at something adjacent.

The gate checks this correspondence mechanically: a cited script or result that
does not exist fails the commit, and so does a printed figure that appears in no
result file. This table is generated from the same markers, so it cannot drift
from what the gate sees.

## `theorem_A.md — the demand side`

| Statement | Claim | Code (`v2/code/`) | Results |
|---|---|---|---|
| Theorem (`thm:A`) | the flat-weighted sum is $\ll_A N(\log N)^{-A}$, unconditionally | analytic | — |
| Corollary (`cor:B`) | hence $E_4$ needs no hypothesis; the demand collapses to $E_3$ | analytic | — |
| Theorem (`thm:C`) | $E_3$ is unconditionally equivalent to Huang--Li's equation (22) | `audit_E3_constant.py` | `results/audit_E3_constant.txt` |
| Proposition (`prop:onesided`) | the demand is one-sided, and its threshold is not a log power | `lab_onesided_margin.py` | `results/lab_onesided_margin.txt` |
| Proposition (`prop:nolog`) | the demand needs no saving in $\log N$ | `lab_onesided_demand.py` | `results/lab_onesided_demand.txt` |
| Proposition (`prop:direct`) | the wall cancels out of the count | `lab_direct_route.py` | `results/lab_direct_route.txt` |
| Proposition (`prop:untrunc`) | the untruncated sum is the count | `lab_direct_identity.py` | `results/lab_direct_identity.txt` |
| Proposition (`prop:layers`) | the count as signed layers | `lab_layer_decomposition.py` | `results/lab_layer_decomposition.txt` |
| Proposition (`prop:combined`) | the count over a combined modulus | `lab_combined_modulus.py` | `results/lab_combined_modulus.txt` |
| Lemma (`lem:complete`) | the complete divisor sum is $\mathbf 1_{\rad(u)\mid N}$ | `audit_switch_identity.py` | `results/audit_switch_identity.txt` |
| Lemma (`lem:degen`) | terms with $(k,N)>1$ contribute $O(N^{o(1)})$ | analytic | — |
| Lemma (`lem:BV`) | $\tau_3$-weighted Bombieri--Vinogradov | analytic | — |
| Lemma (`lem:density`) | the local factor is exactly $p^{-1}$; $1/\zeta$ to the first power | `audit_density_identity.py` | `results/audit_density_identity.txt` |
| Lemma (`lem:mu`) | $G(x)\ll e^{-c\sqrt{\log x}}$ in the relevant range | analytic | — |
| Proposition (`prop:MT`) | the main term dies uniformly in $t$ | analytic | — |
| Proposition (`prop:movingcut`) | the bound survives the corrected cut, unconditionally for $\theta'>1/2$ | analytic | — |
| Lemma (`lem:completelog`) | the $\log k$ branch returns $-\Lambda(u')$ | analytic | — |
| Lemma (`lem:Gb`) | $b=\mu*w$ is exactly the complete divisor transform | analytic | — |
| Proposition (`prop:flatsum`) | both ends are the same sum of dilated walls | `lab_weight_gap.py` | `results/lab_weight_gap.txt` |
| Lemma (`lem:extract`) | extraction | `audit_extraction_tradeoff.py` | `results/audit_extraction_tradeoff.txt` |
| Lemma (`lem:bv`) | BV-accessibility | analytic | — |
| Theorem (`thm:D`) | no weight extracts $C(N)$; loss $\exp(c_1\sqrt{\log N/2})$ | analytic | — |
| Theorem (`thm:Dprime`) | the no-go survives $EH(N^{\theta_E})$ for every $\theta_E<1$ | analytic | — |
| Proposition (`prop:Dpp`) | polynomial weights | analytic | — |
| Proposition (`prop:E`) | zero margin | `audit_circle_margin.py` | `results/audit_circle_margin.txt` |

## `wall_v3.md — the wall`

| Statement | Claim | Code (`v2/code/`) | Results |
|---|---|---|---|
| Conjecture (`conj:L`) | the factorization law | `audit_support_density.py` | `results/audit_support_density.txt` |
| Proposition (`prop:V`) | the exact scale: $V\sim\AAA(N)N\log N$ | `lab_second_moment.py` | `results/lab_second_moment.txt` |
| Lemma (`lem:MP`) | the exact aggregate second-moment identity | analytic | — |
| Proposition (`prop:W`) | the amplification | `audit_amplification.py` | `results/audit_amplification.txt` |
| Lemma (`lem:coin`) | the coin control | analytic | — |
| Lemma (`lem:placebo`) | the placebo key | `lab_mask_placebo.py` | `results/lab_mask_placebo.txt` |
| Proposition (`prop:coindisc`) | the coin obstruction is arithmetic, not informational | `lab_coin_discriminator.py` | `results/lab_coin_discriminator.txt` |
| Proposition (`prop:dilate`) | the demand-side discrepancy is a dilate | `lab_dilate_identity.py` | `results/lab_dilate_identity.txt` |
| Proposition (`prop:posweights`) | the weights are nonnegative | `lab_positive_weights.py` | `results/lab_positive_weights.txt` |
| Lemma (`lem:cellmom`) | exact cell moments | `lab_cellmom_montecarlo.py` | `results/lab_cellmom_montecarlo.txt` |
| Proposition (`prop:coh`) | coherent sums | `lab_cell_floor.py` | `results/lab_cell_floor.txt` |
| Proposition (`prop:placebo`) | the mask survives its placebo | `lab_mask_placebo.py` | `results/lab_mask_placebo.txt` |
| Proposition (`prop:scaleinv`) | the size mechanism is scale-invariant | `lab_cell_singular.py` | `results/lab_cell_singular.txt` |

---

38 numbered statements. Remarks carry their own evidence markers and are
reachable the same way; they are not listed here because they are not claims the
papers stand on.
