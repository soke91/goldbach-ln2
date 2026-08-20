# MANIFEST — 논문 ↔ 진술 ↔ 코드 ↔ 결과

배포본의 모든 진술과 인쇄된 수치가 어디서 왔는지의 대응표.
저장소의 `<!-- evidence: -->` 마커를 대체한다 (SPEC §3.3-5).

- 코드: `deploy/code/<이름>.py` (원본 `v2/code/`)
- 결과: `deploy/results/<이름>.txt` (원본 `v2/results/`)
- 공용 모듈: `deploy/lib/goldbach/`
- 실행: `python deploy/code/<이름>.py` — 각 스크립트는 독립 실행되고,
  자기 사전등록 판정을 출력하며, 실패 시 비-0 종료한다.

---

## P1 — `P1-mobius-fixed-class.tex`

| # | 진술 | 등급 | 근거 |
|---|---|---|---|
| 1 | Theorem (thm:A) | 증명 (§3, 7단계) | 해석적 |
| 2 | Corollary (cor:B) | 증명 (§4) | 해석적 |
| 3 | Theorem (thm:C) | 증명 (§5) | 해석적 |
| 4 | Proposition (prop:onesided) | 증명 (§6) | 해석적 |
| 5 | Proposition (prop:nolog) | 증명 (§6) | 해석적 |
| 6 | Lemma (lem:complete) | 증명 | 해석적 |
| 7 | Lemma (lem:degen) | 증명 | 해석적 |
| 8 | Lemma (lem:BV) | 증명 | 해석적 |
| 9 | Lemma (lem:density) | 증명 | `audit_density_identity` |
| 10 | Lemma (lem:mu) | 증명 | 해석적 |
| 11 | Proposition (prop:MT) | 증명 | 해석적 |
| 12 | Lemma (lem:completelog) | 증명 | 해석적 |

| 수치 | 스크립트 | 결과 |
|---|---|---|
| 교환 항등식 `eq:switch` 이 $10^{-16}$ 이내 | `audit_switch_identity.py` | `.txt` |
| 밀도 항등식, 유리수 산술, 제곱없는 $m<400$ 무불일치 | `audit_density_identity.py` | `.txt` |
| $\mathfrak{A}\widetilde G(1)=-\mathfrak{S}$ 의 부호와 값 | `audit_E3_constant.py` | `.txt` |
| $\theta'$ 스윕, 유한-$N$ 잔차 | `lab_theta_sweep.py` | `.txt` |
| 문턱 상수의 중앙값·최소·원시근사 argmin | `lab_onesided_margin.py` | `.txt` |
| $B(N)/N$ 과 문턱 비 | `lab_onesided_demand.py` | `.txt` |

---

## P2 — `P2-no-go-divisor-switch.tex`

| # | 진술 | 등급 | 근거 |
|---|---|---|---|
| 1 | Lemma (lem:Gb) | 증명 | 해석적 |
| 2 | Proposition (prop:dilate) | 증명 | `lab_dilate_identity` |
| 3 | Proposition (prop:posweights) | 증명 | `lab_positive_weights` |
| 4 | Proposition (prop:flatsum) | 증명 | `lab_weight_gap` |
| 5 | Proposition (prop:untrunc) | 증명 | `lab_direct_identity` |
| 6 | Proposition (prop:layers) | 증명 | `lab_layer_decomposition` |
| 7 | Proposition (prop:combined) | 증명 | `lab_combined_modulus` |
| 8 | Proposition (prop:direct) | 증명 | `lab_direct_route` |
| 9 | Lemma (lem:extract) | 증명 | `audit_extraction_tradeoff` |
| 10 | Lemma (lem:bv) | 증명 | 해석적 |
| 11 | **Theorem (thm:D)** | 증명 | 해석적 |
| 12 | **Theorem (thm:Dprime)** | 증명 | 해석적 |
| 13 | Proposition (prop:Dpp) | 증명 | `audit_polyweight` |
| 14 | Proposition (prop:E) | 증명 | `audit_circle_margin` |

| 수치 | 스크립트 |
|---|---|
| `eq:dilate` 최악 상대오차 $\sim10^{-14}$ | `lab_dilate_identity.py` |
| `eq:posweights` 상대오차 $\le1.8\cdot10^{-16}$ | `lab_positive_weights.py` |
| 평탄합/로그가중합 비, 감쇠 지수 | `lab_weight_gap.py` |
| 새 문턱 대 옛 문턱, 잔차 감쇠 | `lab_direct_route.py` |
| `lem:extract` 항등식 검증, $10\,262$ 개 $d_0$ | `audit_extraction_tradeoff.py` |
| $\mathrm{CP}_2$ 두 조각의 크기 | `audit_polyweight.py` |
| 원 방법 여유 표 (FFT, $N=2^{14..20}$) | `audit_circle_margin.py` |

---

## P3 — `P3-wall-second-moment.tex`

| # | 진술 | 등급 | 근거 |
|---|---|---|---|
| 1 | **Proposition (prop:V)** | 증명 | `lab_second_moment` |
| 2 | **Lemma (lem:MP)** | 증명 | 해석적 |
| 3 | **Proposition (prop:W)** | 증명 | `audit_amplification` |
| 4 | Lemma (lem:coin) | 증명 | 해석적 |
| 5 | Proposition (prop:coindisc) | 증명 | `lab_coin_discriminator` |

| 수치 | 스크립트 |
|---|---|
| $\mathfrak{A}$ 대 $\mathfrak{S}$ 판별, 셀별 $V/(W\mathfrak{A})$, 2차 형태 | `lab_second_moment.py` |
| 증폭 $\Gamma$, 절대 예산 $\mathcal B(X)$ | `audit_amplification.py` |
| $T(x)\to-4/\pi^2$, 역필터 증가 $1.92^m$ | `lab_coin_discriminator.py` |
| iid 널이 곱셈적 널보다 좁음 ($0.6455$) | `audit_cn_coin_deep.py` |
| 극단 여유 $10^{4.466}$, $10^{22.842}$ | `audit_margin.py` |
| 상수 $\mathfrak{A}=0.787275$, $\mathfrak{S}$ | `audit_constants.py` |

---

## P4 — `P4-coherent-cell-floor.tex`

| # | 진술 | 등급 | 근거 |
|---|---|---|---|
| 1 | **Lemma (lem:cellmom)** | 증명 | `lab_cellmom_montecarlo` |
| 2 | **Observation (prop:coh)** | 유도 + 측정 | `lab_cell_floor` |
| 3 | Lemma (lem:placebo) | 증명 | 해석적 |
| 4 | Measurement (prop:placebo) | 순수 측정 + 널 | `lab_mask_placebo` |
| 5 | Proposition (prop:scaleinv) | 증명 + 측정 | `lab_cell_singular` |

| 수치 | 스크립트 |
|---|---|
| 닫힌 형태 대 2000 뽑기 MC, 깊이 0–5 | `lab_cellmom_montecarlo.py` |
| $\mathrm{se}\propto N^{-b}$, $b\approx0.039$ | `lab_cell_floor.py` |
| 플라시보 10 순열, 바닥 붕괴, $D_c$ 대 $\mathrm{se}_c$ | `lab_mask_placebo.py` |
| $D_c$ 규모 지수 $\approx0$ | `lab_cell_singular.py` |

> **등급 주의.** 2번과 4번은 저장소에서 `Proposition` 이었다. 배포본은
> 증명 유무로 이름을 가른다 (SPEC §2, 규칙 T1).

---

## P5 — `P5-negative-map.tex`

| # | 진술 | 등급 | 근거 |
|---|---|---|---|
| 1 | Measurement (지지집합 밀도) | 측정 | `audit_support_density` |
| 2 | **Conjecture (conj:L)** | 추측 — 증거 등급 §Note 2 | `audit_support_density` + 재검증 pass2 |
| 3 | Proposition (prop:R4) | 증명 | `audit_r4_blocks` |

| 수치 | 스크립트 |
|---|---|
| 대역 밀도 $0.3303,0.3298,0.3320$, 소멸집합 $28.0\%$ | `audit_support_density.py` |
| Heath–Brown 무게 몫 $0.772590,0.833180,0.886081$ | `audit_hb_weight.py` |
| R4 블록합 | `audit_r4_blocks.py` |
| 체 변종 대조 목록 | `audit_sieve.py` |

폐쇄 판정 17건(적합 5 + 킬테스트 9 + 표현계급 4)의 통계량은 논문에
싣지 않는다. 각 설계의 판정규칙·널·결과 파일은 원 저장소 `v2/`에 있다.

---

## 전체 파일 목록

```
deploy/
  README.md      배포 안내
  SPEC.md        논문 작성 규약 (양식·규칙·수록/삭제 기준)
  INVENTORY.md   프로젝트 분석 (진술 개수, 논문 주제 수)
  MANIFEST.md    이 파일
  papers/        P1..P5 .tex  (5편)
  code/          인용된 스크립트 30개
  results/       그 결과 파일 30개
  lib/goldbach/  공용 모듈
```
