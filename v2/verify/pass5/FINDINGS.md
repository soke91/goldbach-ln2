# pass5 — 코드 감사 (BRIEF C, 블라인드)

대상 sha256: `A-blind/SEAL.txt` 그대로.
종료 시 재확인: **논문 5편 일치**. `code/` 30개 · `results/` 30개 · `lib/` 8개는
`SEAL.txt`가 트리 해시를 어떤 나열 형식으로 잡았는지 명시하지 않아 그 값을
재현하지 못했다(세 가지 관례를 시도했고 모두 불일치). 대신 감사 시작 시점에
68개 파일의 개별 해시를 직접 기록해 두고 종료 시 재확인했으며, **68개 전부
변경 없음**. 이 패스는 `A-blind/` 안의 어떤 파일도 쓰지 않았고 저장소의
`deploy/`·`v2/code/`도 건드리지 않았다.

> **게이트 제안 0.** `SEAL.txt`가 트리 해시의 나열 형식(구분자·정렬 키·개행)을
> 적지 않으면 그 해시는 재확인 불가능하다. 봉인 자체가 검증 불가능한 봉인이다.

---

## 사전등록

무엇을 어느 깊이로 볼지, 무엇이 나오면 결함으로 볼지를 먼저 고정했다.

1. **체 계층 전수 대조** — `code/`의 모든 μ·Λ·φ 구현을 시행나눗셈 인수분해와
   `n ≤ 10^4` 전 구간에서 대조. 한 개라도 불일치면 결함. (표본 아님)
2. **수치 계층** — FFT 합성곱 경로(9개 스크립트)를 직접합·`math.fsum`과 대조.
   인쇄된 자릿수를 float64가 지탱하지 못하면 결함.
3. **식 충실도** — 논문 문장에서 식을 먼저 재구성하고 코드를 나중에 연다.
   합의 시작·끝 지수, 정규화 분모, μ vs μ², φ vs q, 절댓값 위치가 다르면 결함.
4. **널** — 부호가 μ²의 지지집합 위에서만 뽑히는가, 순열이 실제 순열인가,
   뽑기 수가 인쇄 정밀도를 감당하는가.
5. **죽은 검사** — 항상 참인 조건, 항등적으로 성립해 아무것도 판정하지 않는
   비교, 인쇄 자릿수보다 좁은 tol.
6. 결함마다 **독립 구현으로 값 차이를 실측**하고 판정이 뒤집히는지 계산한다.

---

## 판정 여유표 (작업 순서의 근거)

여유가 작은 순. 결함의 값은 이 여유에 반비례한다.

| 판정 | 근거 스크립트 | 임계 | 측정 | 여유 |
|---|---|---|---|---|
| MC/closed − 1 ≤ 0.05 | lab_cellmom_montecarlo | 0.05 | 0.0320 | **1.6×** |
| placebo max\|z_c\| < 5 | lab_mask_placebo | 5.0 | 3.2040 | **1.6×** |
| B(N) ≤ (1−ε)S(1−𝔄)N | lab_onesided_demand | 1.0 | 1.580 | **미달 (1.6×)** |
| B_H ≤ S·N | lab_direct_route | 1.0 | 0.3338 | 3.0× |
| Γ log N/N → 1/𝔄 | audit_amplification | 1.2702 | 1.3590 | 7% |
| r̃ ≥ S(1−𝔄)N | lab_onesided_margin | 1.0 | 3.7038 | 3.7× |
| max\|z_c\| > Bonferroni | lab_cell_floor | ≈4.5 | 9.06 | 2.0× |
| V/(W·𝔄) = 1 | lab_second_moment | 5e−7 | 2e−7 | 2.5× |
| 예산 ℬ(X) > 1 | audit_amplification | 1.0 | 13.3–30.5 | 13–30× |
| 부호 판정 𝔄·G̃ = −S | audit_E3_constant | 2S·N | 3.52 대 0.00018 | 10^4× |
| switch 항등식 | audit_switch_identity | 5e−10 | 1.0e−16 | 5·10^6× |

여유가 1.6배인 세 자리가 이 패스의 첫 작업 대상이었고, F1·F3이 거기서 나왔다.

---

## 확정 (식이 논문과 다르고, 값 차이를 계산했다)

### F1. [X3] 평탄합 Σ_k H(N;k)에서 k = 1 항이 빠져 있다 — P2 Measurement 3.4의 여덟 숫자와 지수 하나가 바뀐다

- **스크립트 / 줄**: `lab_weight_gap.py:183` — `ks = np.array([k for k in range(2, K) ...])`
- **논문이 정의하는 양**: P2 [eq:Tw]·[eq:Bw]와 Proposition {#prop:flatsum}은
  모든 합을 `k < K, (k,N) = 1`로 색인한다. **k = 1은 이 색인 집합 안에 있다**
  (1 < K, (1,N) = 1, μ(1) = 1). 그러므로
  `T_1 = Σ_{k<K} H(N;k) − B_1 C(N)`의 평탄합은 k = 1부터다.
- **코드가 계산하는 양**: `range(2, K)` — k ≥ 2.
- **차이**: `H(N;1) = Σ_{m<N} Λ(N−m)μ(m) = C(N)`. 즉 평탄합에서 **정확히 C(N)이
  빠진다**. 같은 절단이 B(K)에도 적용되어(k=1 항 = 1) T_1 자체는 불변이다 —
  T_1의 k=1 항은 `A(N;1) − C(N)/φ(1) = 0`이므로 원래도 0이다.
- **값 차이 (실측, `code/c02_flatsum_k1.py` → `results/c02_flatsum_k1.txt`)**:
  내 독립 구현이 대상의 숫자를 먼저 그대로 재현한 뒤(|Σ H|/N = 0.079409 …
  0.014925, 비 0.1807 … 0.1188, 지수 −0.3620 — 전부 일치), k = 1을 복원하면

  | N | 2e5 | 4e5 | 8e5 | 1.6e6 | 3.2e6 | 6.4e6 | 1.28e7 | 2.56e7 |
  |---|---|---|---|---|---|---|---|---|
  | 인쇄된 비 | 0.1807 | 0.1740 | 0.1624 | 0.1389 | 0.1456 | 0.1216 | 0.1258 | 0.1188 |
  | k=1 복원 | 0.1785 | **0.1559** | **0.1484** | 0.1483 | 0.1367 | 0.1334 | 0.1241 | 0.1135 |

  최대 편차 10.4% (N = 4e5). 적합 지수는 **−0.3620 → −0.3505**.
  `|T_1|/N` 열은 자릿수 끝까지 불변임을 확인했다(차이 ≤ 2.3e−18).
- **파급**: 인쇄된 수 **9개** — P2 Measurement {#meas:flatsum}의 비 8개
  (0.1807 … 0.1188)와 적합 지수 −0.3620. 기대는 진술 **1개** —
  「the gap widens, as it must if one side is ≪_A N(log N)^{−A} and the other
  is ≍ N」.
- **판정이 뒤집히는가**: **아니오.** Z2(비 < 0.2)의 여유는 임계의 1.12배이고
  값은 최대 −1.2%p 움직이므로 0.1785 < 0.2로 산다. Z3(평탄합 지수가 더 가파름)은
  −0.3505 대 −0.2658로 여전히 가파르다. **인쇄된 수는 틀렸고 판정은 산다.**
- **근거**: `pass5/code/c02_flatsum_k1.py` → `pass5/results/c02_flatsum_k1.txt`

`range(2, K)`는 8개 스크립트의 공통 관례다(`audit_E3_constant`, `audit_sieve`,
`lab_dilate_identity`, `lab_direct_route`, `lab_onesided_demand`,
`lab_positive_weights`, `lab_theta_sweep`, `lab_weight_gap`). 나머지 7곳은
전부 log k 가중(log 1 = 0)이거나 T_w 형태(k=1 항이 항등적으로 0)라서 **무해함을
개별 확인**했다. 해를 입는 것은 가중이 없는 평탄합 하나뿐이다.

### F2. [X3] P5 §c3 (2)의 세 숫자 중 가운데가, 같은 문장이 배제한 관례에서 나온 값이다

- **스크립트**: `audit_hb_weight.py` (결과 `results/audit_hb_weight.txt`, E4 행)
- **논문이 정의하는 양**: P5는 「only the rounding `z = ⌈M^{1/J}⌉` is admissible
  at all three M, since rounding down violates `z^J ≥ x`」라고 **명시**한 뒤,
  같은 문장에서 J = 8의 j ∈ {6,7,8} 몫을
  `0.772590, 0.833180, 0.886081` (M = 10^4, 10^5, 10^6)로 인쇄한다.
- **코드가 계산하는 양**: 결과 파일은 세 관례(floor/round/ceil)를 모두 인쇄하고,
  `z^J ≥ x` 성립 여부를 열로 표시한다. M = 10^5, J = 8에서
  round는 z = 4 → `z^8 = 65536 < 10^5` → **NO**, ceil은 z = 5 → yes.
- **차이**: 인쇄된 세 숫자는 **round 행**(E4)이다. 논문이 유일하게 허용한
  ceil 행은 M = 10^5에서 **0.840039**다.
- **값 차이**: 0.833180 → **0.840039** (+0.82%). M = 10^4·10^6은 세 관례가
  같은 값이라 불변. (실측 불필요 — 두 값이 대상 결과 파일 안에 나란히 있다.)
- **파급**: 인쇄된 수 **1개**. 기대는 진술 1개 — 「the fraction increases with M」.
- **판정이 뒤집히는가**: **아니오.** 0.7726 → 0.8400 → 0.8861은 여전히 증가하고,
  「majority sits outside a ≤ M^η」도 유지된다. 그러나 **논문이 한 문장 안에서
  배제한 관례의 수를 인용한 것**이고, 이는 게이트가 구조적으로 못 본다 —
  0.833180은 결과 파일에 실재하기 때문이다.
- **근거**: `A-blind/results/audit_hb_weight.txt` (J = 8 표의 ceil M=1e5 행 →
  0.840039) 대 `A-blind/papers/P5-negative-map.tex:411`.

### F3. [X3] P1 Measurement {#meas:relocate}의 「임계의 1.6배 이내」는 N = 2^a·5^b 만으로 된 필드의 성질이다 — 같은 범위의 원시근사 N에서는 32–40배다

- **스크립트**: `lab_onesided_demand.py` (필드 `NS = [2e5, 4e5, 8e5, 1.6e6, 3.2e6]`)
- **논문이 정의하는 양**: P1 [eq:Bsum]의
  `B(N) = Σ_{k<K,(k,N)=1}(log k)|E_μ(N;k)|`를 임계 `S(N)(1−𝔄(N))N`에 대고 재는 것.
  Proposition {#prop:nolog}가 걸리는 대상 그 자체다.
- **코드가 계산하는 양**: 식은 논문 그대로다(확인함). 문제는 **필드**다.
  2e5 = 2^6·5^5, 4e5 = 2^7·5^5, …, 3.2e6 = 2^10·5^5 — **다섯 N 전부 P(N) = {2,5}**.
  그래서 `1−𝔄(N) = 0.2127`과 `S(N) = 1.7604`, 따라서 임계 `0.3745N`이
  다섯 점에서 **같은 수**다(결과 파일의 해당 두 열은 상수 열이다).
  이 sweep은 크기만 바꾸고 임계를 정하는 산술은 한 번도 바꾸지 않는다.
- **값 차이 (실측, `code/c03_demand_field.py` → `results/c03_demand_field.txt`)**:
  내 독립 구현이 대상의 다섯 비를 먼저 그대로 재현했다
  (2.159, 1.975, 1.950, 1.748, 1.580 — 논문 인쇄값과 일치). 같은 범위 안에서
  N의 산술만 바꾸면

  | N | P(N) | 임계/N | B/N | **B/임계** |
  |---|---|---|---|---|
  | 3 200 000 | 2,5 | 0.3745 | 0.5916 | **1.58** |
  | 223 092 | 2,3,6197 | 0.2707 | 0.8956 | 3.31 |
  | 210 210 | 2,3,5,7,11,13 | 0.0873 | 2.7617 | **31.6** |
  | 510 510 | 2,3,5,7,11,13,17 | 0.0733 | 2.9270 | **39.9** |
  | 2 042 040 | 2,3,5,7,11,13,17 | 0.0733 | 2.7208 | 37.1 |

  대상 필드의 최악 2.16 대 여기의 최악 39.9 — **18.5배**. 그리고 움직이는 것은
  주로 임계다(B/N도 0.59 → 2.93으로 오르지만 임계는 0.3745 → 0.0733으로 내린다).
- **덤으로 나온 것**: 원시근사 N에서 `|E_3| = B`가 **정확히** 성립한다
  (N = 210210에서 177개 항 전부 같은 부호, |E_3|/B = 1.00000000). 대상 필드에서
  결과 파일이 「the signs buy a factor 2 to 3」(|E_3|/B = 0.35–0.54)이라고 적은
  그 이득이 여기서는 **정확히 1배, 즉 0**이다.
- **파급**: 인쇄된 수 **10개** — P1 meas:relocate의 B(N)/N 5개와 비 5개.
  기대는 진술 1개 — 「already within a factor 1.6 of the threshold at the top of
  the accessible range, and falling」.
- **판정이 뒤집히는가**: **판정 자체는 아니오** — 논문은 Prop {#prop:nolog}의
  조건이 성립한다고 주장하지 않고 「a constant-factor bound valid for all large
  N is what is needed, and no computation supplies that」로 닫는다.
  그러나 **문장이 전달하는 상태는 뒤집힌다**: 「accessible range」에서 1.6배가
  아니라, 접근 가능한 같은 범위 안에서 40배인 N이 있다.
- **이 패킷 자체가 반증을 담고 있다**: `results/lab_direct_route.txt`의 Y3 표는
  산술적으로 다양한 7개 N에서 같은 옛 비가 1.63 … **37.87**로, 스팬 **23.25**임을
  이미 인쇄해 두었다(P2 Measurement {#meas:direct}도 그 23.25를 인용한다).
  즉 P2는 그 민감도를 측정해 보고했고, P1은 같은 양을 상수-산술 필드에서 재서
  1.6배로 인쇄했다. **두 논문이 같은 양에 대해 서로 다른 인상을 준다.**
- **근거**: `pass5/code/c03_demand_field.py` → `pass5/results/c03_demand_field.txt`

---

## 확정 (죽은 검사 — 아무것도 판정하지 않는다)

### F4. [X2] `audit_switch_identity.py`의 S4는 항등적으로 참이다 — [eq:PR]은 검증된 적이 없다

`audit_switch_identity.py`가 `R = Σ g·(rad_ok − sig)`로 R을 **정의**하고
`P = Σ g·rad_ok`, `D_u = Σ g·sig`를 쓴다. 그러므로 `P − R ≡ D_u`가 대수적으로
성립하고, S4가 재는 `|D_k − (P−R)|/N`은 S1이 재는 `|D_k − D_u|/N`과 **같은 수**다.
결과 파일이 그대로 보여준다 — 25행 중 24행에서 두 열이 자릿수까지 동일하고,
남은 한 행만 float 반올림으로 2.73e−17 대 3.64e−17이다.

파급: S4는 Lemma {#lem:complete}가 거짓이어도 통과한다(완전합을 무엇으로 넣든
P와 R에 같은 배열이 들어가 상쇄되므로). [eq:PR]의 내용을 실제로 지탱하는 것은
S3(모든 무제곱 u < N 전수 대조)뿐이고, S3는 통과한다. **명제는 살지만
「[eq:PR] 검증됨」은 S3의 재진술이지 독립 증거가 아니다.** X2로 등급한 이유는
이 자리가 「기계 정밀도로 검증됨」을 주장하는 항등식 검사이기 때문이다.

부수: 이 스크립트의 필드 N = 2.5e4 … 4e5는 **다섯 개 모두 rad(N) = 10, ω(N) = 2**다.
Lemma {#lem:complete}의 `(k,N)=1` 취급은 P(N) = {2,5}에서만 시험됐고,
`|supp P|`는 다섯 N에서 전부 4이며, `max_t |P(t)|`는 세 N에서 0.0000이다 —
[eq:P]의 상한 검사가 세 점에서 공집합을 상한과 비교한다.

### F5. [X4] `audit_density_identity.py`의 Q4는 항등식이 아니라 상쇄를 적어 놓은 것이다

```python
v = (Fraction(1, p-1) - Fraction(1, p*(p-1))
     - Fraction(1, p*p*(p-1)) + Fraction(1, p*p*(p-1)))
```

뒤의 두 항은 **기호적으로 상쇄**한다. 실제로 검사되는 것은
`1/(p−1) − 1/(p(p−1)) = 1/p`뿐이고, 이는 Q1이 이미 전수로 확인하는 국소인자다.
Lemma {#lem:density} 증명이 요구하는 것 — `(v_p(d), v_p(e)) ∈ {0,1}^2`의 네 항이
정의로부터 **그 네 항이 맞는지** — 은 검사되지 않는다. 정의에서 직접 국소인자를
세우면 (0,0)→1/(p−1), (1,0)→−1/φ(p³), (0,1)→−1/φ(p²), (1,1)→+1/φ(p³)이고
논문의 표시와 일치함을 손으로 확인했다(따라서 **논문은 옳다**). 그러나
스크립트는 그것을 확인하지 않는다.

### F6. [X4] `audit_E3_constant.py`의 수렴 checkpoint가 조용히 사라진다

`G_tilde`/`B_log`가 `if m in checkpoints`로 부분합을 기록하는데, 루프가
`mu[m] == 0`이나 `(m,N) > 1`인 m을 먼저 건너뛴다. checkpoint로 고른
`x//64, x//16, x//4`는 N = 10^6·4·10^6에서 전부 2·5의 거듭제곱이라 μ = 0이다.
결과: **N = 10^6에서도 4·10^6에서도 G̃의 수렴표가 한 줄뿐**이다.
논문이 인용하는 `𝔄(N)G̃(x) = −1.760250 at x = 4·10^6`은 그 한 점이고,
「수렴한다」는 주장은 이 출력으로 확인할 수 없다. FIELD 줄은
「truncations x, K in geometric checkpoints」라고 적고 있어 사실과 다르다.

---

## 확정 (근거 포인터가 가리키는 파일에 그 수가 없다)

### F7. [X4] P4의 두 Measurement가 잘못된 스크립트를 지목한다

- **Measurement {#meas:mc}** 끝의 `(lab_cellmom_montecarlo.py)`가 덮는 마지막
  문단 — 「at the top octave of 1.6·10^7, Var/(Q_cc/n_c²) is
  0.113119, 0.016302, 0.288567, 0.551516 … 0.113118, 0.016303, 0.288571 …
  0.557654 … 266 … 1 687 911」 — 의 숫자는 **`results/lab_cell_floor.txt`**에 있다.
  `lab_cellmom_montecarlo.py`는 밴드 (10^5, 2·10^5]에서만 돌고 X = HI = 2·10^5이라
  1.6·10^7을 계산하지 않는다(grep 확인: 해당 값 0건).
- **Measurement {#meas:Dc}** 끝의 `(lab_mask_placebo.py)`가 덮는 D_c 행
  (0.302138, 0.046766, …)과 상관 0.9805는 **`results/lab_cell_singular.txt`**에 있다.
  `lab_mask_placebo.py`는 S₂(h)도 D_c도 계산하지 않는다(se_c 행만 그 파일에 있다).

「인용된 값이 결과 파일에 실재하는가」를 파일 단위로 묻는 게이트라면 이 둘은
빨간불이고, 저장소 전체를 grep하는 게이트라면 초록불이다. 지금은 초록불이다.

---

## 확정 (수치·인쇄 — 값은 바뀌지만 판정 여유가 크다)

### F8. [X4] `audit_amplification.py`의 `mobius_upto`가 μ(0) = 1을 준다

```python
mu = np.ones(n+1, dtype=np.int64)   # mu[0]이 1로 남는다
```

전수 대조(c01)에서 `code/`의 **유일한** 체 오류다. 사용처 두 곳을 확인했다:
`mu[N−q]`는 q < N이라 인덱스 0에 닿지 않고, `am[1:] = mu[1:Xv+1]`은 0을 배제한다.
**모든 인쇄된 수에 영향 없음**을 확인했다. 다른 스크립트의 μ 배열은 전부
μ(0) = 0이다.

### F9. [X4] `lab_theta_sweep.py`의 FIELD 줄이 코드와 다르다

결과 헤더: 「S(N) as an Euler product over p < 1e6」.
코드: `CLIM = 4_000_000`, `for p in primes_upto(CLIM)`.
값에는 영향이 없다(4·10^6 절단의 S(N) 오차 ≈ 2·10^−8, 인쇄는 6자리).
재현 지시로서는 틀렸다.

### F10. [X4] 인쇄된 밴드가 측정값을 담지 않는다

`lab_cell_floor.txt` C5: 「top 16000000 max|z| = 9.064 **in [9.1,13.0]**」.
9.064는 [9.1, 13.0] 안에 없다. 판정 규칙이 0.05의 여유를 허용하므로 `hold`는
맞지만 인쇄된 문구는 사실과 다르고, P4 본문의 「runs 9.1 to 13.0」은
실측 9.064–12.968을 바깥쪽으로 반올림한 것이다.

### F11. [X4] P1 {#note:thetasweep}의 「≍10^−3 N」이 측정값과 한 자릿수 다르다

논문: 「it reads 0.4558, 0.3729, 0.3108, while the right-hand side it is
compared against is ≍10^{−3}N」.
결과 파일의 그 우변 `|r̃ − S(N−C)|/N`은 **0.01810311, 0.01081673, 0.00633949**다.
가장 작은 것도 6.3·10^−3이고 가장 큰 것은 1.8·10^−2 — 「≍10^−3」의 18배.
비(15.19, 18.38, 26.84)는 논문이 따로 정확히 인용하므로 결론은 영향 없다.

---

## 의심 (식은 논문대로인데 필드·보고 구조가 결론을 지탱하지 못한다)

### S1. N-sweep이 전부 2^a·5^b다 — P1·P2의 적합 지수 전부에 걸린다

`lab_weight_gap`(2e5·2^j, 8점), `lab_theta_sweep`(2e5,4e5,8e5),
`lab_onesided_demand`(2e5 … 3.2e6), `lab_dilate_identity`(2e5,4e5,8e5),
`audit_switch_identity`(2.5e4·2^j), `audit_E3_constant`(2e5 … 1.6e6),
`audit_polyweight`(1e6, 4e6, 1.6e7) — **N이 전부 2^a·5^b**다. 그래서 S(N)과 𝔄(N)이
필드 전체에서 상수이고, 이 필드에서 적합한 지수
(−0.3620, −0.2658, −0.2794, −1.4526 …)는 「크기의 함수로서의 감쇠」가 아니라
「하나의 산술 유형 안에서의 감쇠」다. `audit_polyweight`는 이를 스스로 밝힌다
(「which is the same at all three N because N has the same prime support {2,5}」).
나머지는 밝히지 않는다. F3이 이 관찰의 정량판이다.

### S2. 30개 중 17개의 결과 파일이 자기 사전등록 판정에서 REFUTED로 끝난다

`audit_amplification, audit_constants, audit_margin, audit_polyweight,
audit_r4_blocks, audit_support_density, lab_cell_floor, lab_cell_singular,
lab_coin_discriminator, lab_combined_modulus, lab_direct_identity,
lab_direct_route, lab_layer_decomposition, lab_mask_placebo,
lab_onesided_demand, lab_onesided_margin, lab_positive_weights,
lab_second_moment, lab_theta_sweep` — 그리고 이 중 다수가 `raise SystemExit(1)`로
끝난다. 논문은 이 파일들의 **숫자**를 인용하고, 대부분은 post-hoc DIAGNOSTIC이
왜 예측이 틀렸는지 정직하게 적어 두었으며 배포본은 그 수정을 반영했다
(예: `audit_amplification`의 Γ는 소수 관례 열 1.5128e3이 인쇄되어 있고
Λ 관례 열 1.5489e3이 아니다 — 대상 논문이 옳은 쪽을 골랐다).
그러나 **「스크립트가 통과한다」는 이 corpus에 대해 참인 적이 없고**, 게이트는
값의 실재만 보므로 이 상태를 구조적으로 못 본다. 인용값이 사전등록 통과분인지
post-hoc 진단분인지가 문장에서 구분되지 않는다.

### S3. 절단 상수의 경계가 스크립트마다 다르다

Artin 상수·2C₂를 `4·10^6`(audit_E3_constant, audit_amplification,
audit_density_identity, lab_theta_sweep), `10^7`(lab_second_moment의 ARTIN_LIM),
`측정 범위 자체`(lab_onesided_margin은 1.6·10^7, lab_second_moment의 twin)에서
자른다. 차이는 상대 10^−8 수준이라 인쇄 자릿수 아래지만, `audit_constants.py`가
이 문제를 이미 스스로 감사해 「the last printed digit of the threshold depends on
the bound the script happened to sieve to」라고 결론짓는다. 그 관찰을 재확인했고
동의한다.

### S4. 「기계 정밀도」 주장의 분모가 선언되어 있지 않다

P1 §Step 1은 switch 항등식이 「equal to 10^{−16} **relative to N**」이라고
분모를 명시한다(좋다). 그러나 P2 {#meas:flatsum}의
「holds to a worst relative error of 1.875·10^{−16}」, {#meas:extract}의
「agree to 6.5·10^{−18}」, {#meas:direct}의 잔차들은 분모가 무엇인지 문장에서
결정되지 않는다. switch의 경우 실측 D는 N의 6·10^−2배 크기이므로 N 분모와
D 분모는 자릿수가 다르다(1.0e−16 대 1.6e−15).

---

## 검사했고 통과한 것

깊이를 두 등급으로 구분한다 — **[재구현]** 은 독립 구현으로 값을 재서 대조한 것,
**[정독]** 은 논문 식을 먼저 재구성하고 코드와 문자 단위로 대조했으나 값을 새로
재지는 않은 것.

### 체 계층 — 전수 대조 [재구현]

`pass5/code/c01_sieve_bruteforce.py` → `pass5/results/c01_sieve_bruteforce.txt`.
`code/`의 μ·Λ·φ 진입점 **22개 전부**를 시행나눗셈 인수분해와 `0 ≤ n ≤ 10^4`
**전 구간**에서 대조. 불일치 총계 **1개**(F8의 μ(0)).

- Λ(p^k) = log p 정확(Λ(4) = Λ(8) = log 2 확인), Λ(1) = 0.
- μ의 제곱인수 0, μ(1) = 1, 경계 n = 0·1 확인.
- φ: `audit_extraction_tradeoff.sieve_mu_phi` 전수 일치, φ(1) = 1.
- `audit_sieve.py`의 네 독립 구현(stride / 인수분해 / 점화 / 비트마스크)이
  서로 그리고 내 기준과 일치.
- 스크립트별 μ 구현이 두 계열뿐임을 확인했고(spf 점화형, 부호뒤집기+잔여
  인수형) 두 계열 모두 전수 통과. 체 상한도 확인: 조회 인덱스 `N−mk`, `N−q`가
  전부 체 범위 안이고 int64 곱셈 한계에서 멀다.

### 수치 계층 — FFT 대 직접합 [재구현]

`pass5/code/c04_fft_vs_direct.py` → `pass5/results/c04_fft_vs_direct.txt`.

- 패딩: 대상들이 쓰는 `2^25`는 길이 1.6·10^7 두 배열의 선형 합성곱이 요구하는
  `2X−1 = 31 999 999`를 넘는다 — **되감김 없음**. `lab_cellmom_montecarlo`(2^19),
  `lab_mask_placebo`·`lab_cell_floor`(2^23), `audit_amplification`의
  `autocorr`(≥ 2n+2)도 각각 충분함을 확인.
- 정확도: V(N)은 상대 3.7e−16 … 1.4e−15, r̃(N)은 0 … 2.0e−15,
  C(N)은 절대 2.3e−13 … 1.5e−11(|C| 대비 최악 2.0e−13 — 가장 상쇄가 심한
  C(10^6) = 71.22 지점에서). `math.fsum`(정확 반올림)과의 대조도 같은 자릿수.
  **논문이 인쇄하는 6–7자리를 float64 FFT가 지탱한다.**
- 교차: `audit_amplification`(직접합, 소수거듭제곱 열거)의
  V(4·10^6) = 44684177.8625와 `lab_second_moment`(FFT)의 W/V = 1.270800이
  내 세 경로(직접·fsum·FFT) 전부와 일치.
- 셀 바닥의 3항 상쇄(순열 라벨에서 5자릿수)도 각 항의 float 오차보다
  9자릿수 위임을 확인.

### 정규화 — 논문이 스스로 함정이라고 지목한 자리 [재구현]

P3 {#note:AvsS}가 「One reconstruction error of exactly this kind divides by
W(N) where the definition calls for V(N)」이라고 경고한다.
`lab_second_moment.py`는 `W = np.cumsum(lam2)` 후 **`W[ev−1]`**을 쓴다 —
`W(N) = Σ_{w<N}`의 올바른 오프셋이고 `V[ev]`는 `Σ_{v<N}`. 둘 다 맞다.
직접합으로 V(4e6)·W(4e6)을 다시 재서 W/V = 1.270800을 재현했다.
(주의: `W[N]`을 쓰는 오류 버전도 이 스케일에선 6자리까지 같다 —
**이 수치 검사로는 그 오류를 못 잡는다.** 잡는 것은 코드 정독뿐이다.)

### 구간 끝 관례 [정독]

- `audit_support_density.py`: dyadic 밴드를 `K < k ≤ 2K`로 구현
  (`range(59,117)` = 58개). P5 {#meas:supp}가 **두 관례를 다 인쇄하고**
  차이(0.3303 대 0.3345)가 재려는 효과보다 크다고 명시한다 — 코드와 문장 일치.
  m-범위도 `arange(root+1, N//k+1)` = `√N < m ≤ N/k`로 [eq:E1]과 일치.
- 진행 슬라이스 `f[r::k]`(r = N mod k, 0이면 k)가 `n ≡ N (mod k), 1 ≤ n < N`을
  정확히 훑음을 6개 스크립트에서 확인. `m < N/k`는 `arange(1, (N−1)//k + 1)`로 정확.
- `V(N) = Σ_{v<N}`, `C(N) = Σ_{n<N}`, `r̃(N) = Σ_{n<N}`, `U(N) = Σ_{n<N}`의
  끝점이 FFT 인덱스에서 정확(양끝 항이 Λ(0) = μ(0) = 0으로 죽는다).

### 널 [정독]

- **부호가 μ²의 지지집합 위에서만 뽑히는가**: `lab_dilate_identity`,
  `lab_positive_weights`, `lab_cellmom_montecarlo`, `lab_mask_placebo` 전부
  `supp = flatnonzero(mu != 0)` 위에서만 ±1을 뽑고 나머지는 0으로 둔다.
  따라서 V(N)이 보존된다 — **널이 널이다**.
- **순열이 실제 순열인가**: `lab_mask_placebo`는 `rng.permutation(true_lab)`으로
  라벨만 섞는다. 셀 크기가 정확히 보존되고 Z는 점별로 불변 — Lemma
  {#lem:placebo}가 요구하는 그대로. 복원추출 아님.
- **씨앗**: 난수를 쓰는 모든 스크립트가 `default_rng(20260807/20260808)`로
  고정하고 결과 헤더에 적는다 — 결과 파일만으로 재현 가능.
- **뽑기 수 대 인쇄 정밀도**: `lab_cellmom_montecarlo`의 2000회 뽑기는
  sd 정밀도 `1/√(2·1999) = 0.0158`이다. P4 {#meas:mc}가 인쇄하는
  비 0.9920 … 1.0320은 넷째 자리까지 적혀 있으나 둘째 자리 아래는 잡음이며,
  최악 편차 0.0320은 그 정밀도의 **2σ**다. 논문이 정밀도 수치와
  「their deviations are strongly correlated」를 함께 적어 두었으므로
  결함이 아니라 관찰로 남긴다.
- `lab_second_moment`·`lab_onesided_margin`은 「NULL: none applies」를 선언하고,
  실제로 결정론적 산술 비교뿐임을 확인했다.

### 식 충실도 [정독] — 개별 확인

- `audit_amplification.py`: `Σ_{h≠0}c(h) = θ(N)² − Σ(log p)²`가 [eq:cdef]의
  **소수** 관례와 일치(prime powers 아님). `S(h) = Ah[h]/(X−h)`의 항 수 X−|h| 일치,
  `am[0] = 0`이라 u > |h| 일치. `B = 2Σ_{h=1}^{X−1}c(h)|S(h)|/V`의 인자 2(±h)와
  범위 `0 < |h| < X` 일치. 두 정규화(X−|h| 대 X)를 둘 다 인쇄하고 논문이 둘 다
  인용한다.
- `audit_E3_constant.py`: `λ(m) = Π_{p|m}(1−1/(p(p−1)))^{−1}`,
  `G̃(x) = Σ_{m≤x}μ(m)λ(m)1_{(m,N)=1}log m/m`, `𝔄(N) = Π_{p∤N}(1−1/(p(p−1)))`
  전부 논문과 일치(짝수 N이라 p = 2 배제가 무해함을 확인). 𝔄(N) = 0.787275가
  Artin/((1−1/2)(1−1/20))와 일치. 순차 누적의 오차 상한 ≈ 1.6·10^−9로
  판정 여유(3.52)의 10^9배 아래.
- 셀 3종(`lab_cellmom_montecarlo`, `lab_cell_floor`, `lab_mask_placebo`):
  `u_c(v) = Σ_{N∈c}Λ(N−v)/√V(N)`를 `IFFT(conj(F_Λ)·F_b)[v]`로 얻는 교차상관의
  **방향**을 직접 확인(`Σ_t Λ(t)b(t+v)`, t = N−v). Lemma {#lem:cellmom}의 세 항
  `Q_cc/n_c² − 2Q_ca/(n_c n) + Q_aa/n²`을 정의에서 손으로 재유도해 코드와 일치.
- `lab_theta_sweep.py`: `R(N,θ') = |E_3 − (r̃ − S(N)(N−C))|/N`이 논문 정의와
  문자 단위 일치. `E_3 = Σ μ(k)log k(A(N;k) − C/φ(k))` — 분모가 φ(k)이지 k가
  아님을 확인.
- `lab_second_moment.py`: `V`, `W`, `𝔄(N)`, `S(N)`의 스트라이드 오일러 곱
  (`AN[p::p] /= (1−1/(p(p−1)))`, `SN[p::p] *= (1+1/(p−2))`, p > 2 조건) 전부 일치.
- `lab_onesided_margin.py`: `U(N) = Σ_{n<N}Λ(n)μ²(N−n)`, `C(N)`, `r̃(N)`을
  Λ\*μ²·Λ\*μ·Λ\*Λ로 만드는 세 합성곱의 인덱스 전부 일치. 임계
  `S(N)(1−𝔄(N))`가 Prop {#prop:onesided}와 일치.
- `audit_density_identity.py`: `c_{D,E}(m)`의 `m·lcm(d²,e)`와 `1_{(d,N)=1}`,
  `Σ_{g|m}μ(g)/(φ(m/g)·g·φ(g))` 전부 일치. `Fraction`으로 **진짜 유리수 산술**
  (float+tol 아님) — BRIEF §3.3의 「exact 주장」 항목 통과.
- `lab_dilate_identity.py` / `lab_positive_weights.py`: [eq:dilate]와
  [eq:posweights]가 코드와 일치. A(N;k)를 n-슬라이스로, H를 m-열거로 얻어
  **순환 논증이 아님**을 확인(같은 배열을 두 번 읽지 않는다).
- `audit_switch_identity.py`: [eq:switch] 양변의 인덱스 순서가 실제로 다르고
  (k는 모듈러 슬라이스, u는 divisor-sum 배열), σ_K를 stride와 약수열거 두
  경로로 만들어 대조 — 순환 아님. S3는 모든 무제곱 u < N 전수.

### 스크립트 사이의 일관성 [정독]

- S(N) = 1.760432, 𝔄(N) = 0.787275, 1−𝔄 = 0.2127이 `audit_E3_constant`,
  `lab_theta_sweep`, `lab_onesided_demand`, `lab_onesided_margin`,
  `lab_second_moment`, 그리고 내 독립 구현에서 전부 일치.
- 1/𝔄(4·10^6) = 1.270204가 `audit_amplification`과 `lab_second_moment`에서 일치.
- `results/lab_direct_identity.txt`가 `AGREE` 줄로 `lab_signed_level.py`·
  `lab_layer_decomposition.py`와의 교차 대조를 인쇄한다(다섯 N에서 전부 AGREE).
  이 corpus는 교차 검증이 실제로 걸려 있는 편이다.
- 내 `c02`·`c03`의 독립 구현이 `lab_weight_gap`과 `lab_onesided_demand`의
  인쇄값을 **모든 자릿수에서 재현**했다(0.079409 … 0.014925 / 2.159 … 1.580 /
  지수 −0.3620). 두 스크립트의 산술 자체에 오류는 없다 — 틀린 것은 색인 집합과
  필드다.

---

## 안 본 것

정직하게 적는다. 30개 중 등급별로 다음과 같이 봤다.

- **[재구현 + 정독] 12개** — audit_switch_identity, audit_density_identity,
  audit_E3_constant, audit_amplification, lab_weight_gap, lab_onesided_demand,
  lab_onesided_margin, lab_second_moment, lab_theta_sweep,
  lab_cellmom_montecarlo, lab_mask_placebo, lab_cell_floor.
- **[정독] 5개** — lab_dilate_identity, lab_positive_weights,
  audit_support_density, lab_coin_discriminator(H1–H3만), audit_constants.
- **[결과 파일 + 핵심 루프만] 6개** — audit_hb_weight(F2가 여기서 나왔다),
  lab_direct_route, lab_direct_identity, audit_polyweight,
  audit_extraction_tradeoff, audit_sieve.
- **[결과 파일만 — 식을 코드와 대조하지 않았다] 7개** —
  `audit_circle_margin`, `audit_margin`, `audit_r4_blocks`, `audit_cn_coin_deep`,
  `lab_cell_singular`, `lab_combined_modulus`, `lab_layer_decomposition`.
  **이 일곱은 감사되지 않았다.** 특히 `audit_cn_coin_deep`(널이 붙은 측정 —
  BRIEF §4 우선순위 4), `audit_margin`(M1·M2 REFUTED),
  `audit_r4_blocks`(D3·D6 REFUTED), `lab_combined_modulus`(REFUTED)가
  다음 패스가 먼저 볼 자리다.
- **`lib/goldbach/` 8개 파일** — 형식적으로 훑었으나, **30개 스크립트 중 어느
  것도 이 라이브러리를 import하지 않는다**(`sys.path`·`import goldbach` 둘 다
  0건). 29개가 각자 체를 다시 구현한다. 패킷에 들어 있으나 감사 대상의 실행
  경로에 없다.
- 패킷에 없는데 참조되는 스크립트: `lab_signed_level.py`, `lab_sign_structure.py`,
  `lab_lean_decay.py`, `audit_truncation_exponent.py`,
  `audit_threshold_arithmetic.py`. 이들의 값이 `AGREE` 줄과 Measurement
  {#meas:direct}의 7점 표에 들어와 있으나 **이 패스는 그 코드를 볼 수 없었다.**
  F3의 보강 근거 중 하나(스팬 23.25)가 그중 하나에서 온다.
- P1 §4–§7의 해석적 논증, P2 no-go의 증명, P5의 17개 방향 서술 —
  **수학 재검증은 이 패스의 일이 아니다**(A의 일). 코드가 계산하는 양이
  논문이 정의한 양인지만 봤다.

---

## 게이트로 걸 제안 (수정은 이 패스의 일이 아니다)

1. `SEAL.txt`는 트리 해시의 나열 형식(구분자·정렬 키·개행)을 명시한다. (F0)
2. 결과 파일의 최종 판정 줄을 게이트가 읽고, REFUTED로 끝나는 파일을 인용하는
   논문 문장은 인용값이 **사전등록 통과분인지 post-hoc 진단분인지** 선언한다. (S2)
3. 상대오차를 인쇄하는 줄은 `DENOM:`을 함께 선언한다. (S4)
4. 오일러 곱 상수는 한 곳에서 고정 경계로 계산하고, 결과 파일이 `CONST-BOUND:`로
   자기 경계를 선언한다. (S3)
5. Measurement의 스크립트 인용을 **파일 단위로** 검증한다 — 인용된 값이 인용된
   그 파일 안에 있는지. (F7)
6. N-sweep을 쓰는 결과 파일은 `FIELD-ARITH:` 줄로 필드에 등장하는 P(N)의 집합을
   인쇄한다. 다섯 점이 전부 같은 P(N)이면 그 자리에서 보인다. (F3, S1)
7. `k < K, (k,N) = 1`로 색인된 합을 계산하는 스크립트는 **k = 1 항의 값을 별도
   줄로 인쇄**한다. 0이면 무해가 증명되고, 아니면 F1이 보인다. (F1)
8. 「convergence checkpoint」를 인쇄하는 루프는 요청된 checkpoint 수와 실제
   인쇄된 수가 같은지 검사한다. (F6)

---

## 요약

- **X1(닫힌 경로를 잘못 닫음) 0건.** 음성 판정을 떠받치는 계산 —
  audit_amplification의 예산 ℬ(X) = 13.3–30.5, lab_theta_sweep의 오차/신호 =
  15–27, audit_E3_constant의 부호, lab_second_moment의 V 대 W — 을 식 단위로
  대조했고 전부 논문대로였다. 여유도 크다(13배 – 10^4배).
- **X2 1건.** F4 — [eq:PR]의 검사가 항등적으로 참이라 아무것도 판정하지 않는다.
  같은 스크립트의 S3가 실질을 지탱하므로 명제 자체는 산다.
- **X3(인쇄된 수가 바뀌지만 판정은 삼) 3건.** F1(9개 수), F2(1개 수),
  F3(10개 수 + 한 문장의 인상).
- **X4 6건**, 의심 4건, 게이트 제안 8건, 미감사 7개 스크립트.

이 corpus는 방어가 잘 되어 있다. 체 계층은 전수 대조에서 오류 1개(무해)뿐이고,
FFT 계층은 인쇄 자릿수를 지탱하며, 널은 지지집합을 보존하고 순열은 진짜 순열이다.
논문이 스스로 지목한 함정(V 대 W, dyadic 관례, 상수 절단)은 전부 올바르게 처리되어
있다. 찾은 것은 그보다 한 겹 바깥 — **색인 집합의 끝점 하나, 같은 문장이 배제한
관례에서 인용한 수 하나, 그리고 산술이 상수인 필드에서 잰 비를 「접근 가능한 범위」의
성질로 적은 자리 하나**다.

가장 비싼 것은 F3이다. 결함이 아니라 필드 선택이므로 어떤 자동 검사도 볼 수 없고,
이 패킷 자체가 그 반증(스팬 23.25)을 다른 파일에 인쇄해 두고 있는데도 두 논문이
같은 양에 대해 서로 다른 인상을 남긴다.
