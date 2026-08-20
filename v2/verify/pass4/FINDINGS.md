# pass4 — 블라인드 수학 재검증

대상 sha256 (SEAL.txt 그대로):

```
9dc9b2aa2aa6b83159425ebf8a717c2bb1feba04634ceab952e61085e674d8d9  papers/P1-mobius-fixed-class.tex
00035519620a527bb8f90d33c076483fdd731a89d728a32e8c31c51f241a0c29  papers/P2-no-go-divisor-switch.tex
ef3df86975554dd11d5b3ff6e1c957fa7a4c164ed76dee1013dc6277636a4bbb  papers/P3-wall-second-moment.tex
e7f6c4eed65231b6a2e5d3c08704d0d7ee35579c53dc84c61d4d87df3a75bffd  papers/P4-coherent-cell-floor.tex
c4e84934c5b75ba7029c36677fc49b32297980d6545ad0fcad9469c35778ac1f  papers/P5-negative-map.tex

supporting trees (sha256 over the sorted name:hash listing):
0d85bf6f53c8b66b0505e3dc1c253dfb2952baacdcec64d58ed11e0737a874ad  code/  (30 files)
c6eb846f71482a8710c19e0994a8b09b7089a61d7134d2c734e1f8111b358411  results/  (30 files)
```

**종료 시 재확인: 논문 다섯 편 일치.** 시작·종료 두 번 계산했고 다섯 개
모두 SEAL과 같다.

두 개의 supporting tree 는 **확인 못 함**. SEAL이 적은 조리법("sha256 over
the sorted name:hash listing")으로는 어떤 인코딩을 써도 두 값이 나오지
않는다 (`name:hash\n`, 마지막 개행 유무, `hash  name`, `hash name`,
`name hash`, `dir/name:hash` 등 일곱 가지를 시도했고 전부 불일치). 파일
개수 30/30은 맞다. 대신 내가 계산한 지문을 시작·종료에 각각 남겼고
**둘 다 변하지 않았다**:

```
code/ (*.py 30개)  8ce53a1e181602a8a67e4876d28d1fce0512a62de182903771caca998fd9f97e  (시작=종료)
results/ (30개)    32bc47f35036391685441032f688ac5eb68e55c73b0b4d4c441516dc518c414b  (시작=종료)
```

리뷰 도중 `A-blind/code/__pycache__/` 가 새로 생겼다 (30개 `.pyc`,
타임스탬프 21:23). `.py` 원본 30개의 mtime과 해시는 그대로이므로 표적은
움직이지 않았다. 나는 그 스크립트를 실행하지도 import 하지도 않았다 —
같은 패킷을 쓰는 다른 패스(C)의 흔적으로 보인다. 기록만 해 둔다.

---

## 0. 이 패스가 한 일

논문 다섯 편을 처음 보는 심사자로 읽었다. 저장소는 열지 않았고 패킷 밖
파일은 읽지 않았다. 인쇄된 수는 `results/` 의 해당 파일에 대조했고, 그중
핵심적인 것은 **논문의 진술만 보고 독립 구현해서 다시 계산했다**
(`pass4/code/a1`–`a6`). 남의 스크립트는 읽기만 했고 다시 돌리지 않았다.

**판정 요약: 다섯 편의 중심 판정 중 뒤집히는 것은 없다.** 가장 무거운 두
결함(F1, F2)은 P2의 no-go 안에 있고, 둘 다 **no-go를 더 강하게 만드는
방향으로** 틀렸다. 나머지는 문장·정의·수치 라벨의 결함이고, 그중 F3·F4·F5·
F7은 문장이 말하는 것과 계산된 것이 다른 경우다.

---

## 확정 (재현 가능한 결함)

### F1. P2 (5)식은 항등식이 아니다 — 이 주제 전체의 대상인 $T_w$ 를 빠뜨렸다

- **위치**: P2 §1.2, eq.(5) (PDF 4쪽). Theorem 14의 증명이 이 식을 먹는다
  ("Feeding it into (5) ... gives the bound on $|C(N)|$").
- **주장**: "Running the divisor switch of [6] with the weight $w$ gives the
  identity
  $B_w\,C(N) = \sum_{u<N}\Lambda(N-u)\mu^2(u)b_u - \mathcal R_w + O(N^{o(1)})$."
- **문제**: 정확한 재배열은
  $B_w C(N) = \mathrm{CP} - \mathcal R_w - T_w$
  이다 ($T_w$ 는 같은 논문 (2)식). 인쇄된 식은 $T_w$ 를 통째로 버리고 그
  자리를 $O(N^{o(1)})$ 이라고 적었다. $T_w$ 는 $w=\log$ 에서 $E_3(\alpha)$
  자신이고, $N^{o(1)}$ 이 아니다.
- **근거**: `pass4/code/a4_extract_identity.py` → `pass4/results/a4_extract_identity.txt`.
  $N=2\cdot10^5,4\cdot10^5,8\cdot10^5$, $\theta'=0.56$ 에서 다섯 양을 각자
  정의대로 따로 계산했다.

  | $N$ | $w$ | (5)식 잔차 / $N$ | $T_w$ 복원 시 잔차 / $N$ |
  |---|---|---|---|
  | 200000 | $1$ | 0.07843 | 2.6e-18 |
  | 200000 | $\log$ | **0.43769** | 1.7e-16 |
  | 400000 | $\log$ | **0.38372** | 3.0e-17 |
  | 800000 | $\log$ | **0.31716** | 2.2e-17 |

  $T_w$ 를 넣으면 기계 정밀도로 항등식이 되고, 인쇄된 대로면 $\asymp N$
  만큼 어긋난다. 빠진 항의 크기 0.43769, 0.38372, 0.31716 은 패킷 자신의
  `lab_weight_gap.txt` 가 인쇄한 $|E_3|/N=0.4377,\,0.3837,\,0.3172$ 와
  같은 수다 — 버려진 것이 정확히 $E_3$ 임을 확인해 준다.
- **영향**: Theorem 14는 **음성** 주장이므로 결론은 살아남고 오히려 강해진다
  ($T_w$ 를 되돌리면 이 경로가 줄 수 있는 상한은 더 나빠진다). 다만
  (a) (5)식은 항등식이 아니고, (b) $T_w$ 를 복원하면
  $|C(N)|\ll\exp(c_1\sqrt{\tfrac12\log N})N(\log N)^{-A}$ 의 도출에
  $T_w$ 의 무조건 상한이 필요해지는데, 이 방법이 주는 유일한 무조건 상한은
  $|T_w|\le \mathrm{CP}+\mathcal R_w+N|B_w|$ 이고 마지막 항이 **자명한
  $+N$** 을 도로 집어넣는다. Theorem 14는 $T_w$ 를 명시한 채로 다시
  써야 한다.

### F2. P2 Proposition 21(ii)는 $D=1$ 에서 이항 골드바흐를 "고전 체 상계"라고 부른다

- **위치**: P2 Proposition 21(ii)와 그 증명 (PDF 13–14쪽).
- **주장**: "If $f$ is a single monomial $x^D$, $D\ge1$, then every term of
  $\mathrm{CP}_D$ is nonnegative ... **by classical sieve bounds**
  $\mathrm{CP}_D\asymp N(\log N)^{D-1}$, with a fixed sign. **It is never
  $o(N)$.**" 증명은 "the order of magnitude is that of the sieve upper and
  lower bounds for prime-plus-almost-prime sums".
- **문제**: $D=1$ 이면 $\mathrm{CP}_1=\sum_{p<N}\Lambda(N-p)\log p$ 이고,
  하한 $\mathrm{CP}_1\gg N$ 은 바로 이항 골드바흐의 하한이다. $r=1$ 조각은
  prime-plus-almost-prime이 아니라 prime-plus-**prime** 이고, 거기에 체
  하한은 없다 — 같은 논문 Note 15가 인용하는 패리티 장벽이 그것이다.
  게다가 같은 명제의 (iii)이 이 사실을 명시한다: "the $r=1$ piece with
  $D=1$ is $\sum_p\Lambda(N-p)\log p$, **whose asymptotic is the assertion
  to be proved**". (ii)와 (iii)이 서로 모순이다.
- **근거**: 논문 내부 대조. 외부 계산 불필요.
- **영향**: no-go가 실제로 필요로 하는 것은 "$\mathrm{CP}_D$ 가 $o(N)$ 임이
  **알려져 있지 않다**" 뿐이고 그것은 참이므로, (ii)를 그렇게 약화하면
  결론은 그대로다. $D=1$ 은 P1 Theorem 3이 이미 닫은 골드바흐 가지이기도
  하다. 그러나 인쇄된 대로면 골드바흐 등가 하한을 고전 문헌에 귀속시킨
  것이 된다 — BRIEF §3.5가 더 위험하다고 적은 쪽의 실패다.

### F3. P2 Measurement 11은 두 비율 목록을 같은 분자로 적었는데, 두 번째 목록의 분자는 다른 합이다

- **위치**: P2 Measurement 11 (PDF 12쪽).
- **주장**: "with $B_H(N)=\sum_{k<K}(\log k)|H(N;k)|$, the new ratio
  $B_H(N)/(\SS(N)N)$ reads $0.4578,\,0.4064,\,0.4079,\,0.3769,\,0.3338$ ...
  against $2.1591,\,1.9747,\,1.9500,\,1.7483,\,1.5798$ for the old ratio
  $B_H(N)/(\SS(N)(1-\AAA(N))N)$".
- **문제**: 두 번째 목록은 $B_H$ 로 만들어지지 않는다. 그것은
  $B(N)=\sum(\log k)|E_\mu(N;k)|$ — P1 (9)식의 양 — 으로 만들어진다.
  $E_\mu(N;k)=\mu(k)H(N;k)-C(N)/\varphi(k)$ 이므로 절댓값 안에서 평균항이
  살아남아 $B\neq B_H$ 다.
- **근거**: `pass4/code/a1_bh_vs_b.py` → `pass4/results/a1_bh_vs_b.txt`
  (네 가지 분자를 전부 독립 계산).

  | $N$ | $B_H/(\SS(1-\AAA)N)$ (논문이 적은 식) | $B/(\SS(1-\AAA)N)$ | 인쇄된 값 |
  |---|---|---|---|
  | 200000 | 2.1519 | **2.1591** | 2.1591 |
  | 400000 | 1.9105 | **1.9747** | 1.9747 |
  | 800000 | 1.9175 | **1.9500** | 1.9500 |
  | 1600000 | 1.7720 | **1.7483** | 1.7483 |
  | 3200000 | 1.5691 | **1.5798** | 1.5798 |

  "new" 목록은 $B_H/(\SS N)$ 으로 다섯 자리 전부 재현된다 (0.4578, 0.4064,
  0.4079, 0.3769, 0.3338). 패킷 자신의 `lab_direct_route.txt` Y3 열 이름은
  `old B/(S(1-A)N)` 으로 **맞게** 적혀 있다 — 논문이 옮기면서 $B$ 를
  $B_H$ 로 바꿨다.
- **영향**: 이 측정의 논지 전체가 Note 10의 "**What has changed is the
  threshold constant**, from $\SS(1-\AAA)$ to $\SS$" 이다 — 즉 분자를
  고정한 채 분모만 바뀌었다는 것. 실제로는 분자도 바뀌었다. 정성적 결론
  (새 비율은 내내 1 미만, 옛 비율은 내내 1 초과, 일곱 점 스프레드 2.00 대
  23.25)은 살아남는다. 살아남지 못하는 것은 적힌 대로의 주장이다.

### F4. P1 (9)식의 $B(N)$ 정의는 자기 수치가 쓰는 제곱무관 제한을 빠뜨렸다

- **위치**: P1 Proposition 6, eq.(9); Measurement 18.
- **주장**: $B(N):=\sum_{k<K,(k,N)=1}(\log k)|\Emu(N;k)|$, "the quantity
  Huang--Li reach by the triangle inequality"; Measurement 18은
  $B(N)/N=0.8086,\,0.7395,\,0.7303,\,0.6547,\,0.5916$ 을 인쇄한다.
- **문제**: $E_3=\sum_k\mu(k)(\log k)\Emu$ 에 삼각부등식을 쓰면 나오는 것은
  $\sum_k\mu^2(k)(\log k)|\Emu|$ — **제곱무관 $k$ 만** — 이다. 인쇄된
  수는 그 합이다. 정의대로 (모든 $k$) 계산하면 0.8106, 0.7564, 0.7422,
  0.6617, 0.5976 이 나온다.
- **근거**: `pass4/results/a1_bh_vs_b.txt` 의 `B_all/N` 과 `B_sqf/N` 열.
  `B_sqf` 가 인쇄된 다섯 수와 자릿수까지 일치한다.
- **영향**: Proposition 6은 적힌 대로도 **참이다** (모든-$k$ 합이 제곱무관
  합보다 크므로 가설이 더 셀 뿐이다). 깨지는 것은 정의와 측정의 대응이다.
  비율은 1.5798 → 1.5958 로 움직이므로 "already within a factor 1.6" 은
  어느 쪽으로도 성립한다.

### F5. P3 §5의 두 수치는 §1.3이 명시적으로 주장하지 않는다고 한 법칙에서 나온다

- **위치**: P3 §1.3 "What is not claimed" 대 §5 "The margin the problem
  demands" (PDF 9쪽).
- **주장(§1.3)**: "**No law for $C(N)$ is conjectured here.** In particular
  no distributional statement of the form $C(N)=m(N)+\sqrt{V(N)}\,G(N)$ is
  asserted." 그리고 "**No level or rate for $\rho=\mathrm{Var}\,C/V$ is
  quoted.** By Lemma 9 the centred estimator cannot distinguish $\mu$ from a
  sign pattern, so any such figure calibrates nothing."
- **주장(§5)**: "With $\max|C|\approx a_n\sqrt{V(N)}$ and $a_n$ the Gumbel
  location ... which is $10^{4.466}$ at $N=10^{12}$ and $10^{22.842}$ at
  $N=10^{50}$".
- **문제**: $\max|C|\approx a_n\sqrt{V}$ 에 검벨 위치를 쓰는 것은
  $C(N)=\sqrt{V(N)}G(N)$ 을 독립 가우스 장으로 놓는 것, 즉 $\rho=1$ 을
  가정하는 것이고, 정확히 §1.3이 주장하지 않는다고 한 형태다. 더구나 그것은
  Lemma 9의 **동전 모형** 자체이므로, Lemma 9에 따르면 그 두 수는 $\mu$ 에
  대한 진술이 아니라 동전에 대한 진술이다. 그런데도 유효숫자 네 자리로
  인쇄되고, 이 논문이 문제 크기에 대해 내놓는 유일한 정량 진술이다.
- **근거**: 논문 내부 대조. 산술 자체는 재현된다 —
  $\sqrt{10^{12}}/(\sqrt{2\log(10^{12}/2)}\sqrt{0.787275\cdot\log10^{12}})
  =2.92\cdot10^4=10^{4.466}$ ✓, $10^{50}$ 에서 $10^{22.842}$ ✓
  (`audit_margin.txt` M4와도 일치).
- **영향**: §5의 결론("expected to be small by a wide margin and provably
  small by none")은 이 두 수 없이도 성립한다. 두 수는 삭제하거나, "이것은
  동전 모형의 여유이고 Lemma 9에 의해 $\mu$ 에 대한 측정이 아니다"라고
  명시해야 한다. 지금은 자기 §1.3과 자기 Lemma 9를 동시에 위반한다.

### F6. P5의 경로 수가 맞지 않고, 초록과 요약이 서로 다른 수를 말한다

- **위치**: P5 제목, 초록(2회), §3 소절 제목, Summary.
- **문제**:
  - 표의 행을 기계적으로 센 결과 **5 + 9 + 4 = 18**이다 (adjudications 5:
    MRT15/Li20, Tao16, Li23, Dirichlet-polynomial, Partial slices;
    designs 9: K1–K4, R1–R5; classes 4: C-I, C-II, C-IV, C-III). 소절
    제목들이 붙인 (5)/(9)/(3 closed, 1 open) 도 합이 18이다. 그런데 본문은
    네 곳에서 "**seventeen**"이라고 적는다 (제목 line 36, 초록 line 52와
    line 80, Summary line 498).
  - 초록: "**Sixteen close**; one is reopened here"; Summary: "**Fifteen
    close**, K1 is reopened, and C-III is open". 표의 판정으로 세면
    4 Blocked + 8 dead + 3 closed = **15** 이고 K1·C-III가 open, "Partial
    slices" 는 Partial이다. Summary가 맞고 초록이 틀렸다.
  - 제목 "seventeen pre-registered **closures**" 는 닫힌 것이 15개인 본문과
    두 겹으로 어긋난다.
- **근거**: `pass4/results/a5_grades.txt` 및 표 행 grep (본 문서 §근거 명령).
- **영향**: 수학은 아니지만 이 논문의 산출물이 곧 "지도"이므로, 지도의
  항목 수가 맞지 않으면 지도로 쓸 수 없다. 세 곳(제목·초록·요약)을 18과
  15로 맞춰야 하고, "Partial slices" 를 닫힌 것으로 셀지 말지를 정해야
  한다.

### F7. P4 Measurement 12는 표 하나에 두 개의 표본 해상도를 섞었고, 표제 주장이 그 섞음에 의존한다

- **위치**: P4 Measurement 12 (PDF 8쪽).
- **주장**: "Measured over $(10^6,2\cdot10^6]$, $(2\cdot10^6,4\cdot10^6]$
  and $(4\cdot10^6,8\cdot10^6]$, the fitted exponents are"
  $-0.000879,\,-0.016226,\,+0.000936,\,+0.002904,\,+0.000241,\,-0.008428$
  "--- **zero to three decimals at every depth but one**."
- **문제**: `lab_cell_singular.txt` 의 M2/M4 표를 보면, 선언된 field
  (셀당 400000 표본쌍)에서의 지수는
  $+0.012633,\,+0.065546,\,+0.000936,\,+0.002904,\,+0.000241,\,-0.008428$
  이다. 논문 표의 depth 0과 depth 1 항목은 **열 배 표본**의 값
  ($-0.000879$, $-0.016226$) 이다. 논문은 depth 1에 대해서만 재표본을
  공개한다 ("its fitted exponent flips sign from $+0.065546$ to $-0.016226$
  when the sample is taken ten times larger"); **depth 0의 $+0.012633\to
  -0.000879$ 은 공개하지 않는다.**
- **영향**: 선언된 field 그대로면 "zero to three decimals at every depth but
  one" 의 예외는 하나가 아니라 **둘**이다 ($+0.013$ 과 $+0.066$).
  Proposition 11 자체는 증명되어 있으므로 명제는 무사하고, 무너지는 것은
  측정이 명제를 뒷받침하는 강도다. 표에 어느 항목이 어느 해상도인지 적거나
  전 depth를 같은 해상도로 다시 재야 한다.

### F8. P4의 "conservative by about two orders of magnitude" 는 신호가 없는 셀의 수치다

- **위치**: P4 초록, Note 9, Summary (세 곳 모두 같은 표현).
- **주장**: "the floor **collapses** by factors from $3.8$ to $105$: it is
  itself a second measurement of that correspondence, and quoting against it
  is conservative by **about two orders of magnitude**."
- **문제**: 붕괴 배수는 depth별로 105 (depth 0), 70 (depth 2), 43.6
  (depth 3), 17.0 (depth 4), **3.79 (depth 5)** 이다. 그런데 효과를
  나르는 셀은 depth 3·4·5 이고 (Measurement 8: $z=-2.40,\,-6.00,\,-11.03$),
  그중 가장 센 depth 5에서 보수성은 **3.8배 — 반 자릿수** 다. Note 9의
  바로 다음 문장이 그것을 그대로 보여 준다: "with the permuted floor in the
  denominator the depth-5 cell would read $z\approx-42$ rather than $-11$"
  ($=3.8$배). "두 자릿수"는 효과가 $z=+0.11$ 인 depth 0의 수다.
- **근거**: `lab_mask_placebo.txt` L3 (se_true/se_perm) 및
  `pass4/results/a6_cellfloor.txt` (내가 독립 재계산한 se_c·z_c).
- **영향**: "보수적이다"라는 방향은 맞고 결론은 살아남는다. 배수만
  셀별로 적으면 된다.

### F9. P4 Note 6의 "$5.8$ to $160$" 은 어떤 정규화로도 얻어지지 않는다

- **위치**: P4 Note 6 (PDF 6쪽) 및 Summary.
- **주장**: "at the top octave the exact floor exceeds
  $\mathrm{sd}(Z)/\sqrt{n_c}$ by factors of $5.8$ to $160$, **growing with
  the cell**".
- **근거**: `pass4/code/a6_cellfloor.py` → `pass4/results/a6_cellfloor_top.txt`.
  Lemma 1을 진술만 보고 독립 구현했고, 최상위 옥타브 $(8\cdot10^6,1.6\cdot10^7]$
  에서 인쇄된 $n_c$, $se_c$, $z_c$ 를 **전부 다섯 자리까지 재현했다**
  ($se_c=1.18233$e-1 … $4.11529$e-1; $z_c=+0.0440,+0.6905,-0.0561,-1.5529,
  -4.4589,-9.0642$). 그 위에서 문제의 비율은

  | depth | 0 | 1 | 2 | 3 | 4 | 5 |
  |---|---|---|---|---|---|---|
  | $n_c$ | 1534466 | 1687911 | 654281 | 114017 | 9059 | 266 |
  | $se_c/(\mathrm{sd}(Z)/\sqrt{n_c})$ | 158.39 | 62.46 | 123.73 | 84.18 | 32.62 | **7.26** |

  ($\mathrm{sd}(Z)=0.924687$ 실측. 널값 $\mathrm{sd}(Z)=1$ 로 잡으면
  146.46 … 6.71.)
- **문제**: (i) 하한이 5.8이 아니라 7.26 (널값 규약이면 6.71). (ii) 두
  끝점의 **비**는 $\mathrm{sd}(Z)$ 규약에 무관한 불변량인데, 실측은
  $158.39/7.26=21.8$ 이고 논문의 $160/5.8=27.6$ 이다 — 어떤 하나의
  $\mathrm{sd}(Z)$ 로도 두 끝점을 동시에 낼 수 없다. (iii) "growing with
  the cell" 은 depth 1에서 깨진다: depth 1이 가장 큰 셀
  ($n_c=1687911$)인데 비율은 62.46으로 depth 0의 158.39보다 작다. 이는
  같은 논문이 다른 곳에서 지적하는 $Q_{cc}/n_c^2$ 의 비단조성과 같은
  현상이다.
- **영향**: 정성적 주장("절대적으로도 한두 자릿수 좁다")은 살아남는다.
  끝점 두 개와 "growing with the cell" 을 고쳐야 한다.

### F10. P4 Note 6의 "factor $140$" 은 그것이 귀속된 구간이 아니다

- **위치**: P4 Note 6.
- **주장**: "Over a factor $140$ in $N$ --- the sub-range of
  $10^5<N\le1.6\cdot10^7$ on which the exact floor was fitted ---
  $n_c^{-1/2}$ says the error bar shrinks by $11.8$; it shrinks by $1.21$,
  against the $1.19$ that $(\log N)^{-1/2}$ predicts."
- **문제**: $10^5$ 부터 $1.6\cdot10^7$ 은 배수 **160** 이지 140이 아니다.
  그리고 Measurement 5의 적합은 "eight octaves from $(6.25\cdot10^4,
  1.25\cdot10^5]$ up to $(8\cdot10^6,1.6\cdot10^7]$" — $N$ 으로는 배수
  **256**, 옥타브 중점끼리는 **128** — 위에서 이뤄졌다.
  인쇄된 세 수는 오직 $R=140$ 과만 맞는다:
  $\sqrt{140}=11.83$, $140^{0.039451}=1.2153$,
  $\sqrt{\log(1.6\cdot10^7)/\log(1.6\cdot10^7/140)}=1.1934$.
  $R=160$ 이면 12.65 / 1.2217 / 1.2003, $R=128$ 이면 11.31 / 1.2110 /
  1.1889 이다.
- **영향**: 표제 결론 "about ten times too narrow" 은
  $11.8/1.21=9.75$ (140), $12.65/1.22=10.4$ (160) 로 어느 쪽이든 산다.
  구간과 배수를 일치시키면 된다.

### F11. P3 Measurement 12의 "about $430$ standard deviations" 는 표준편차가 아니라 20회 최댓값이다

- **위치**: P3 Measurement 12 (PDF 8쪽).
- **주장**: "At the top the separation is about $430$ **standard
  deviations** of the sign ensemble."
- **근거**: `pass4/code/a2_identities.py` (I6) →
  `pass4/results/a2_identities.txt`. $x=2\cdot10^6$ 에서 독립 부호 200회
  추첨의 표준편차는 $4.61\cdot10^{-4}$ 이고, $|T(x)|=0.405295$ 이므로
  $|T|/\mathrm{sd}=\mathbf{879.0}$ 이다. 인쇄된 430은
  $0.405295/0.000935=433.5$ 이고, $0.000935$ 는 패킷의
  `lab_coin_discriminator.txt` 가 "coin: max $|T_\varepsilon(x)|$"
  (20회 추첨의 **최댓값**)으로 인쇄한 수다.
- **영향**: 방향은 보수적이므로(실제 분리는 두 배 크다) 결론은 무사하다.
  다만 "표준편차"라는 이름이 붙은 수가 표준편차가 아니고, 이 논문의 §4는
  통째로 "널을 먼저 재고 나서 임계를 정하라"는 절이다. 패킷의 결과 파일에는
  표준편차가 인쇄되어 있지 않다.

### F12. 동반 논문 상호참조가 내부 라벨 이름으로 되어 있고, 하나는 다른 진술을 가리킨다

- **위치와 사실**(PDF 인쇄 번호로 확인):
  | 인용 | 가리키는 곳 | 실제 인쇄 번호 |
  |---|---|---|
  | P2 Prop. 9: `[6, Thm. C]` | P1의 Theorem C | P1에 "Theorem C"는 없다 — **Theorem 3** |
  | P2 Lemma 13 증명: `[6, eq. (R1)]` | P1의 (R1)식 | P1에 "(R1)"은 없다 — **(16)** |
  | P2 Note 10: `[6, Prop. 5]` | "the constant-factor condition" | 그것은 P1 **Proposition 6**; P1 Proposition 5는 **부호 있는 $E_3$ 의 한쪽 조건** |
  | P5 R5 행: `[P2, Prop. E]` | P2의 Prop. E | P2에 "Proposition E"는 없다 — **Proposition 21** |
  | P4: `[P3, Prop. 1]` | P3 Proposition 1 | **맞음** (유일하게 인쇄 번호로 인용) |
- **영향**: 앞의 세 개는 독자가 찾지 못하는 포인터다. **`[6, Prop. 5]` 는
  그보다 나쁘다**: 존재하는 다른 명제로 조용히 보낸다. eq.(13)은
  절댓값 합에 대한 상수배 조건이므로 P1 Prop. 6과 같은 모양이고, P1
  Prop. 5는 부호 있는 한쪽 조건이다 — P1이 Note 4에서 가장 애써 갈라놓은
  구분이 바로 그 둘이다.

### F13. P3의 `\texttt` 두 곳이 탭 문자로 깨져 있고, 그대로 PDF에 인쇄된다

- **위치**: `P3-wall-second-moment.tex` 408행, 497행. 백슬래시 자리에
  리터럴 **TAB** 이 들어가 `<TAB>exttt{...}` 가 되었다.
- **PDF 결과**: `P3.pdf` 에 " extttaudit cn coin deep.py" (Note 10 끝) 와
  " extttaudit margin.py, extttaudit constants.py" (§5 끝) 로 인쇄된다.
- **근거**: `grep -n "exttt" P3*.tex | grep -v '\texttt'` 가 두 행을
  집어내고, `pdftotext` 출력의 382행·452행이 깨진 문자열을 보여 준다.
- **영향**: LaTeX는 오류도 경고도 내지 않는다 (그래서 "warnings 0"이
  잡지 못했다). §5는 자기 수치의 스크립트 포인터를 이 한 줄에만 가지고
  있고, 그 두 스크립트(`audit_margin.py`, `audit_constants.py`)는 P3의
  Reproducibility 표에도 없다 — 즉 §5의 두 수치는 렌더된 문서에서
  근거로 가는 길이 완전히 끊긴다.

### F14. P3 Note 10의 $0.6455$ 는 패킷에 없는 파일에서 온 절반을 쓴다

- **위치**: P3 Note 10.
- **주장**: "resolved to $512$ draws on one statistic, the i.i.d. spread was
  $0.6455$ of the multiplicative".
- **문제**: `code/audit_cn_coin_deep.py` 는 이 비율의 분모(0.63553)를
  `results/audit_cn_multnull_deep.txt` 에서 읽고, 자기 교정점은
  `results/audit_cn_class.txt` 에서 읽는다 (`read_mark(SRC, "POINT
  classm2_%d")`). **두 결과 파일 모두, 그리고 그것들을 만드는 스크립트도
  패킷에 없다.** 그래서 이 스크립트는 이 패킷 안에서 실행되지 않고,
  0.6455 는 재현되지 않는다.
- **부수**: `audit_cn_coin_deep.py`, `audit_margin.py`,
  `audit_constants.py` 셋 다 P3 본문이 인용하는데 P3의 Reproducibility
  표에는 없다 (표에는 세 개만 있다).
- **영향**: Note 10이 하는 일은 "i.i.d. 널은 옳은 널보다 좁다"는 경고이고,
  그 경고 자체는 방향이 보수적이라 논문의 다른 주장을 위협하지 않는다.
  다만 그 유일한 정량 근거가 패킷 안에서 확인 불가다.

### F15. P5 §6 (2)의 가운데 수치는 그 문장이 배제한 반올림에서 나온 값이다

- **위치**: P5 §6, item (2).
- **주장**: "the $j\in\{6,7,8\}$ share at $J=8$ reads
  $0.772590,\,0.833180,\,0.886081$ at $M=10^4,10^5,10^6$ ... and **only the
  rounding $z=\lceil M^{1/J}\rceil$ is admissible at all three $M$**, since
  rounding down violates $z^J\ge x$".
- **문제**: `audit_hb_weight.txt` 의 $J=8$ 표에서 $M=10^5$ 은
  floor/round $\to z=4$, share $0.833180$, `z^J >= x` **NO**
  ($4^8=65536<10^5$); ceil $\to z=5$, share $0.840039$, **yes**. 즉
  인쇄된 $0.833180$ 은 같은 문장이 부적격이라고 선언한 반올림의 값이다.
  적격한 세 값은 $0.772590,\,\mathbf{0.840039},\,0.886081$ 이다.
- **영향**: 정성적 주장(다수가 $a\le M^\eta$ 밖에 있고 $M$ 과 함께 증가)은
  $0.7726,\,0.8400,\,0.8861$ 로도 그대로 성립한다.

### F16. P1의 문자 $B$ 가 네 가지 뜻을 나른다

- $B(K)=\sum_{k<K,(k,N)=1}\mu(k)/\varphi(k)$ (Step 0, 369행)
- $B(N)=\sum_{k<K,(k,N)=1}(\log k)|\Emu(N;k)|$ (eq.(9))
- $B_{\log}(K)$ (§5 (ii), 776행)
- Lemma 10의 자유 지수 $B$ ($\tau_3(q)^B$, 이어서 "with $B=3$")
- **문제**: $B(K)$ 와 $B(N)$ 은 인수만 다르고 $K<N$ 이므로, 지면에서 둘을
  가르는 것이 아무것도 없다. 이 논문 §1.1은 스스로 기준을 세운다:
  Huang--Li의 $A$ 를 $\AAA$ 로 고쳐 쓰며 "**one symbol must not carry two
  meanings**". (P2는 여기에 $B_w$, $B_H(N)$, $B(s)=W(s)/\zeta(s)$ 를
  더한다.)
- **근거**: `pass4/results/a5_grades.txt` G2.

### F17. P1 Note 15의 "$\asymp10^{-3}N$" 은 한 자릿수 틀렸다

- **주장**: "Over $N=2\cdot10^5,\,4\cdot10^5,\,8\cdot10^5$ at
  $\theta'=0.56$ it reads $0.4558,\,0.3729,\,0.3108$, while the right-hand
  side it is compared against is $\asymp10^{-3}N$."
- **문제**: `lab_theta_sweep.txt` 는 그 세 $N$ 에서 우변을
  $0.01810311,\,0.01081673,\,0.00633949$ 배 $N$ 으로 인쇄한다 — 즉
  $\asymp10^{-2}N$. $10^{-3}N$ 에 이르는 것은 $N=1.6\cdot10^6$
  ($1.07\cdot10^{-3}$, `audit_E3_constant.txt`) 이다.
- **영향**: 없음에 가깝다. 두 문장 뒤에 같은 비교가 비율
  $15.19,\,18.38,\,26.84$ 로 정확히 다시 나온다.

### F18. P3 §5의 $\AAA=0.787275$ 는 "generic even $N$ 의 값"이 아니다

- $0.787275$ 는 소인수 지지집합이 $\{2,5\}$ 인 $N$ (예: $10^6$) 의
  $\AAA(N)$ 이다 (`audit_density_identity.txt` 의 `N = 1000000 P(N) =
  [2, 5] A(N) = 0.78727541`). 작은 소수 중 2만 나누는 짝수라면
  $\AAA=\text{Artin}/(1-\tfrac12)=0.747912$ 다.
- **영향**: 인쇄된 지수에 미치는 영향은 $\log_{10}$ 으로 0.01 남짓이라
  무시할 만하다. 틀린 것은 수가 아니라 라벨이다.

### F19. 그 밖의 확정된 소소한 것들

- **P4 Measurement 3**: "the three shallow ratios agree to **six digits** ---
  $0.113118,\,0.016303,\,0.288571$" — 상대인 최상위 옥타브 값은
  $0.113119,\,0.016302,\,0.288567$ 이므로 세 번째 쌍은 소수 넷째 자리에서
  갈린다. 여섯 자리 일치가 아니다.
- **P4 Measurement 12**: "the sampling error of the pair average is the
  **same $0.0013$ throughout**" — `lab_cell_singular.txt` 의 se(D_c)는
  $0.001358,\,0.001318,\,0.001359,\,0.001370,\,0.001183,\,\mathbf{0.000162}$
  로 depth 5에서 한 자릿수 작다. 이 문장은 depth 1이 유일한 이상점이라는
  논증에 쓰인다.
- **P2 Measurement 11**: "reaches $2\%$ of $\SS$ only near $N=10^{10.16}$" —
  `lab_direct_route.txt` 는 "0.020 at N = 10^10.16" 을 인쇄하는데 그 양은
  residual/$N$ 이다. $\SS(N)=1.760$ 의 2%라면 residual/$N=0.0352$ 이고
  더 일찍 도달한다.
- **P2 Measurement 22**: "The first row reproduces $\sqrt{6/\pi^2}=0.7797$
  **exactly**" — 첫 행 첫 항목은 0.7798이다.
- **P4 §4 무번호 문단** ("Over every octave from $6.25\cdot10^4$ to
  $1.6\cdot10^7$, $\max_c|z_c|$ runs $9.1$ to $13.0$ ... depths $3,4,5$ at
  $z=-1.6,-4.5,-9.1$"): 수치는 `lab_cell_floor.txt` C5/C6과 맞고 내가
  독립 재계산으로도 확인했다. 문제는 이 문단에 스크립트 인용이 없고,
  Reproducibility 표가 `lab_cell_floor.py` 를 Measurement 5에만 걸어
  두었다는 것 — 독자에게 포인터가 없다.
- **P1 Step 4의 두 꼬리 상계**: $d>D_0$ 꼬리는 윗줄이
  $\ll\log N\sum_m\sum_d(N/(md^2)+1)$ 이므로 두 번째 항이
  $\sqrt{NM}\log N=N^{1-\theta'/2}\log N$ 인데 $N^{1-\theta'/2}$ 로
  인쇄된다. $e>E_0$ 꼬리의 두 번째 항도 $M\log^2N$ 인데 $M$ 으로
  인쇄된다. $\theta'$ 고정·$A$ 임의이므로 둘 다 무해하지만, 표시된 식이
  바로 위 줄이 주는 것과 다르다.

---

## 의심 (근거가 아직 하나)

- **P1 §5 (i)의 "$\sum f(m)/m$ 과 $\sum f(m)$ 이 둘 다 죽는다"**.
  Lemma 12는 $G(x)=\sum_{m\le x}f(m)/m$ 만 다룬다. $\sum_{m\le x}f(m)$ 은
  $G$ 에서 부분합으로 나오고 $\ll xe^{-c\sqrt{\log x}}$ 가 맞을 것이나,
  논문은 "by Lemma 12" 라고만 적는다. 표준적인 한 줄이 빠진 것으로 보이고,
  결론에는 영향이 없다고 판단하지만 독립 확인은 못 했다.
- **P1 §5 (i)의 절단 예산**. $\log k$ 가 붙으면 Step 4의 두 꼬리에
  $\log N$ 이 하나 더 붙는데 $D_0,E_0$ 은 (10)식 그대로 쓰인다. $A$ 가
  임의이므로 무해하지만 명시가 없다.
- **P3 Lemma 5의 $h<0$**. $M(h)=\sum_{v,\,v+h\le X}\mu(v)\mu(v+h)$ 는
  음의 $h$ 에서 $v+h\ge1$ 조건이 진술에 없다. $M(-h)=M(h)$ 이므로
  실질 문제는 아니지만 진술은 불완전하다.

---

## 확인 못 함

**여기가 이 패스가 멈추는 자리다.** 아래는 "통과"가 아니라 "못 봄"이다.

1. **[HL] arXiv:2005.03811v2 의 본문 전부.** 패킷에 사본이 없고 블라인드
   조건상 구할 수 없다. 따라서 다음이 전부 미확인이다 — 그들의 (7), (10),
   (18), (22)식, Lemma 1, Lemma 4, Theorem 1, Corollary 1의 실제 진술;
   $E_3,E_4$ 가 정말 그들의 §3에서 그 형태로 나오는지; **P1 §7 전체**
   (그들의 (18)이 $n$ 의존 제약을 떨어뜨렸다는 결함 보고와 $\Delta$ 의
   형태). P1 §7은 다른 논문에 대한 **결함 주장**이므로 이 미확인은
   가볍지 않다. 이 패스는 $\Delta$ 의 내부 대수만 확인했다
   ($m\le\alpha$ 항의 $\mu(k)\mu(N-n)=\mu(m)\mu^2(k)\mathbf 1_{(k,m)=1}$,
   $\log(k/(N-n))=-\log m$ 은 맞다).
2. **[Vau88]** $\|S_\Lambda\|_1\gg N^{1/2}$. P2 Proposition 21(ii)가
   여기에 기댄다. 측정치 $\|S_\Lambda\|_1/\sqrt N=1.946\to2.346$ 은
   부합하지만 정리 자체는 확인 못 했다.
3. **[Mir49]** 미르스키 정리를 P1 Prop. 5 / P3 Prop. 1이 쓰는 정확한
   형태로. 수치($U/(\AAA N)$ 평균 0.999996)는 상수를 뒷받침한다.
4. **[Bom76], [Tao16], [MRT15], [Li20], [Li23]** 의 보조정리 가정문.
   P1 §8이 Li20의 $\mu$ 에 대한 $3/5$ 와 Li23의 소수에 대한 $66/107$ 을
   갈라 놓는 문단, P2 Note 15가 봄비에리 점근체에 귀속시키는 것, P5
   adjudication 표의 다섯 판정 — 전부 원문 대조가 필요하고 못 했다.
5. **P5의 열여덟(열일곱) 판정 중 열여섯 개.** 논문 자신이 "The verdicts of
   Section 3 are stated without their supporting statistics ... recorded in
   the repository accompanying this report" 라고 적고, 그 저장소는 이
   패킷에 없다. 패킷으로 확인 가능한 것은 Measurement 2, §6 (2)의 HB
   가중치, §5의 R4 블록 뿐이다.
6. **P4 Measurement 3의 시뮬레이션.** 닫힌 형식은 독립 구현으로
   재현했지만 2000회 추첨은 다시 돌리지 않았다.
7. **SEAL의 supporting tree 해시 두 개** (위 §0).
8. **P3 Note 10이 인용하는 다중곱 널** (F14).

---

## 검사했고 통과한 것

무엇을 봤는데 문제가 없었는지. **"통과"이고 "안 봄"이 아니다.**

### 등급 (BRIEF §3.1)
- 다섯 편의 모든 번호 붙은 진술을 열거하고 각각 증명 블록의 유무를 봤다
  (`pass4/results/a5_grades.txt` G1). **Theorem / Corollary / Proposition /
  Lemma 중 증명 없는 것은 하나도 없다.** P1의 Theorem 1·3, Corollary 2,
  Proposition 5·6은 진술 자리에 증명이 없지만 각각 §3, §5, §4, §6에
  전용 증명(또는 `\begin{proof}[Proof of ...]`)이 있다.
- `Observation` / `Measurement` 는 전부 측정 범위를 진술 안에 적고 있다
  (P1 Meas 17 "Over even $N\le1.6\cdot10^7$", P4 Meas 8 "Run on the octave
  $(2\cdot10^6,4\cdot10^6]$" 등).
- "표준적이다 / 잘 알려져 있다" 로 닫히는 증명은 하나 — P2 Prop. 21(ii)의
  "by classical sieve bounds" — 이고 그것이 F2다.

### P1
- **Theorem 1의 증명 Step 0–7을 줄 단위로 읽었다.** 완성(completion)
  단계와 그것이 없으면 $m<N^{1-\theta'}$ 가 거짓이 된다는 지적, Lemma 8의
  퇴화 논증($p\mid(k,N)$ 이면 $N-mk=p^\ell$), Lemma 10의 코시–슈바르츠
  ($\tau_3^{2B}$ 자명 상계 $\times$ BV 지수 $2A+C_1$ → $N(\log N)^{-A}$),
  $(m,d,e)\mapsto q$ 의 중복도 $\le\tau(q)^3\le\tau_3(q)^3$ 과 $B=3$,
  레벨 $q\le MD_0^2E_0\le N^{1/2-\delta/2}$, $T_m(t)=\min(t,N-mK)$ 가
  $n\equiv N\ (q)$ 와 정확히 동치인 점, Proposition 13의 아벨 합
  ($T_M(t)=0$, 도함수 $-K$) — 전부 맞는다.
- **Lemma 11 (밀도)** — 유리수 정확 산술로 독립 재계산.
  $[1,400)$ 의 제곱무관 $m$ 243개, 불일치 **0**
  (`a2_identities.py` I1). 국소 인자가 정확히 $p^{-1}$ 이고 따라서
  $1/\zeta$ 가 1차로 온다는 지적도 맞다.
- **§5의 상수** — 기호적으로 확인. $\AAA H(1)=\SS(N)$ 의 국소 짝짓기
  $(1-\tfrac1{p(p-1)})(1-\tfrac1{(p^2-p-1)(p-1)})=1-\tfrac1{(p-1)^2}$ 는
  항등식이다 ($(p^2-p-1)(p-1)-1=p^3-2p^2=p^2(p-2)$). $B_{\log}(\infty)
  =-f'(0)=-h(0)=-\SS(N)$ 도 확인했다 ($h(0)=\prod_{p\mid N}\frac{p}{p-1}
  \prod_{p\nmid N}(1-\frac1{(p-1)^2})=\SS(N)$; $p=2$ 는 $N$ 이 짝수라
  첫 곱에 들어가 계수 2를 준다).
- **(30)식 수치** — 독립 재계산: $\AAA(N)\widetilde G(4\cdot10^6)
  =-1.760250$, $\SS(N)=1.760432$, 부호 음수 (`a2_identities.py` I2).
  인쇄된 값과 소수 여섯째 자리까지 일치.
  (부기: `results/audit_E3_constant.txt` 의 마지막 줄은
  "finding stands: ... has A(N)*Gtilde(1) = S(N)" 로 **부호가 빠져 있다**.
  자기 표는 $-1.760250$ 을 인쇄한다. 논문은 부호를 맞게 적는다.)
- **Lemma 14와 골드바흐 합으로의 이행.** $u_N>1$ 항이
  $O(N^{o(1)}\log^2N)$ 이라는 주장은 처음에 의심스러워 보이나 옳다:
  $p\mid u_N$, $p\mid N$ 이면 $p\mid N-u$ 이므로 $\Lambda(N-u)\ne0$ 이
  $N-u=p^\ell$ 을 강제하고, 그런 $n$ 은 $\ll(\log N)^2$ 개뿐이며 각각
  $u$ 를 결정한다.
- **Corollary 2** 의 아벨 합
  ($\sum a_n\log(N-n)=\int_1^{N-1}\Emu(t;k)/(N-t)\,dt$) 과 $\log N$ 손실.
- **Proposition 5의 임계 크기.** primorial에서
  $1-\AAA\asymp\sum_{q>p}\frac1{q(q-1)}\asymp\frac1{p\log p}$, $p\asymp\log N$
  → $\gg1/(\log N\log\log N)$; $\AAA(N)<1$ 은 모든 $N$ 에 대해 참.
  $U(N)\sim\AAA(N)N$ 의 국소밀도 논증도 맞다.
- **Measurement 17 — 여덟 수치와 네 argmin 전부 정의에서 독립 재계산**
  (`a3_margin_extract.py` J1): median 0.333459, 0.001-분위 0.119639,
  최소 0.060890 at $N=9699690$, $\times\log N\log\log N$ 최소 2.482019 at
  510510 / 1.737799 at 2310 / 0.810202 at 30, 그리고 가장 작은 열 개
  전부 일치.
- Measurement 18의 다섯 수치는 `lab_onesided_demand.txt` 와 일치
  (정의 문제는 F4).
- Note 15의 세 잔차 0.4558, 0.3729, 0.3108 과 $\theta'$ 스윕
  (0.170167 at 0.51, 4.991178 at 0.90), 비율 15.19/18.38/26.84 — 전부
  `lab_theta_sweep.txt` 와 일치.
- **§1.4 "What is not claimed" 대 본문**: 다섯 항목 모두 본문과 부합.
  특히 Note 4가 세 진술(가설 $E_3\ll_AN(\log N)^{-A}$ / (22) / 점근)을
  갈라 두는 것은 정확하고, 초록이 본문보다 세게 말하는 자리를 찾지 못했다.

### P2
- **Proposition 2 (dilate), 6 (untruncated), 7 (layers), 8 (combined
  modulus), 9 (wall cancels)** — 전부 재유도. 특히 Prop. 9의 부호 계산이
  P1 Theorem 3과 맞물리는 것을 확인했다:
  $\tilde r=\SS N+\sum_k(\log k)H-C(B_{\log}+\SS)+O_A$.
- **Lemma 12 (extraction)** — 독립 구현으로 정확히 확인.
  $N=99999998$, $K=30199$, $\lfloor K/d_0\rfloor=7$ 인 $d_0$ **184개**
  (논문과 일치), 무차별 $B_w$ 대 인수분해형의 최대 차 **0.000e+00**,
  $\gcd(d_0,15)$ 별 값 $0.25/0.75/0.5/1.0$ (논문과 일치)
  (`a3_margin_extract.py` J2). 네 값이 나오는 이유도 확인했다:
  $7\mid N$ 이므로 $j=7$ 항이 $(j,d_0N)=1$ 에 걸려 죽고 $j\in\{1,3,5\}$ 만
  남는다. 표가 $\lfloor K/d_0\rfloor$ 만으로 색인되지 않는다는 논문의
  경고도 옳다.
- **Theorem 14 / Theorem 16의 지수 산술**: $K/D=N^{1/2+\delta}$,
  $K/D=N^{\theta'}/N^{\theta_E-1+\theta'}=N^{1-\theta_E}$, 그리고
  $\theta_E=1$ 에서 진술이 공허해진다는 지적 — 맞다.
- **Note 19**: $W=\zeta''-c\zeta'$, $B=W/\zeta$ 는 모든 $c$ 에 대해
  $s=1$ 에서 이중극을 남긴다 ($\zeta',\zeta''$ 에 $(s-1)^{-1}$ 항이 없고
  $1/\zeta(s)=(s-1)-\gamma(s-1)^2+\cdots$).
- **Proposition 21(i)** 과 $\sqrt{6/\pi^2}=0.7797$; Measurement 22의
  여섯 행 전부 `audit_circle_margin.txt` 와 일치.
- Measurement 5의 모든 수치 (최악 상대오차 $1.875\cdot10^{-16}$, 비율
  0.1807…0.1188, $|T_1|/N=0.01425$ 대 $|T_{\log}|/N=0.1245$, 지수
  $-0.3620$ 대 $-0.2658$, dilate 항등식 오차 세 개, posweights $1.8\cdot
  10^{-16}$) — 전부 결과 파일에 실재. 비율 열은 $|\sum H|/|\sum(\log k)H|$
  이고 논문 문장이 그렇게 말한다 (결과 파일 헤더가 이것을 $|E_3|$ 로 잘못
  적었으나 논문은 맞다).
- Measurement 20의 열두 수치 전부 `audit_polyweight.txt` 와 일치. 세 $N$
  이 같은 소인수 지지 $\{2,5\}$ 를 가진다는 지적도 맞다
  ($10^6=2^65^6$, $4\cdot10^6=2^85^6$, $1.6\cdot10^7=2^{10}5^6$).
- 초록의 "Seven finite rearrangements (Propositions 2--9)" — 실제로
  일곱 개다.

### P3
- **Proposition 1** 의 국소밀도 논증: $q\mid N$ 에서 1, $q\nmid N$ 에서
  $1-1/\varphi(q^2)=1-1/(q(q-1))$ — 맞다. 따라서 $\AAA$ 이고 $\SS$ 가
  아니라는 Note 2의 요지도 맞다.
- **Lemma 5 (aggregate second moment)** — 독립 구현으로 $X=400,800,1600$
  에서 상대오차 $\le5.7\cdot10^{-16}$ (`a2_identities.py` I3).
- **Note 8 (절단이 하중을 받는다)** — 독립 확인. 상자/단체 비는
  $X=800,1600,3200$ 에서 $1.5382,\,1.5658,\,1.5656$ ("near 1.57" ✓),
  $X=6400$ 에서 1.6838 — 1로 가지 않는다 ✓ (`a2_identities.py` I4).
- **Proposition 3**: $\sum_{h\ne0}c(h)=\theta(N)^2-\sum_p(\log p)^2
  =N^2(1+o(1))$, $\Gamma\sim N/(\AAA\log N)$, 그리고 $c(h)\ge0$ 이라
  삼각부등식이 sharp라는 논증 — 맞다. (15)식이 $c$ 를 소수에 대한 합으로
  못박은 덕분에 논문의 $\Gamma$ 수치(1.5128e3, 1.8412e4, 3.5759e5)가
  자기 정의와 일관된다. 결과 파일의 "pub" 열은 소수멱 변형이라 다르다 —
  **논문 쪽이 맞다.**
- **Proposition 11(i)** $T(x)\to-4/\pi^2$ (홀수 제곱무관 밀도) 와
  **(ii)** $\widehat\Lambda(0)=\log2\ne0$ 에 의한 가역성 — 맞다.
  $|a(200)|^{1/200}=1.9164$ ("about $1.92^m$") 도 일치.
- **Measurement 3 — $V$ 와 $W$ 를 독립 FFT로 다시 계산해 확인**
  (`a3_margin_extract.py` J3). $\AAA$ 열은 컷오프 $10^5,5\cdot10^4,10^4,
  10^3,10^2$ 에서 $0.000323,\,0.000346,\,0.000398,\,0.000470,\,0.000529$ 로
  인쇄된 값과 정확히 같고, **$\SS$ 열은 다섯 컷오프 모두 정확히
  $0.245235$** — 즉 "while the figure for $\SS$ is unaffected" 는 참이다.
  이 절은 패킷의 결과 파일이 인쇄하지 않는 유일한 항목이었다.
- Measurement 4의 여섯 수치와 이차형식 $\AAA(N\log N-N)$ 의 근거
  ($1-1/\log N=0.939716$) — 일치.
- Lemma 9와 Note 10의 논리 방향("reproduced $\Rightarrow$ not measuring
  $\mu$", 역은 아님)은 정확하고, 역을 쓰는 자리를 본문에서 찾지 못했다.
  — 단 §5가 그 동전 모형으로 수치를 낸다 (F5).

### P4
- **Lemma 1 — 진술만 보고 독립 구현했고, 두 옥타브에서 인쇄된
  $n_c$, $se_c$, $z_c$ 를 전부 다섯 자리까지 재현했다**
  (`a6_cellfloor.py`; `a6_cellfloor.txt`, `a6_cellfloor_top.txt`).
  $(2\cdot10^6,4\cdot10^6]$: $se=1.24001$e-1 … $4.36854$e-1,
  $z=+0.1069,+1.1133,-0.2314,-2.3990,-5.9997,-11.0258$.
  $(8\cdot10^6,1.6\cdot10^7]$: $z=+0.0440,\dots,-9.0642$.
  세 항이 모두 필요하다는 Note 2의 지적도 옳다.
- **Lemma 7 (placebo key)** 의 증명과 Measurement 8의 열 번 치환
  (최대 3.2040, 평균 1.6741, 상관 $-0.9106$ 대 $+0.0042$), 바닥 붕괴
  105 (depth 0) / 3.8 (depth 5), Note 9의 $z\approx-42$ 산술
  ($-11.0258\times0.43685/0.11518=-41.8$) — 전부 확인.
- **Proposition 11 (scale invariance)** 의 증명 논리(셀이 $3,5,7,11,13$
  법 잉여를 고정하므로 $h$ 의 국소법칙이 스케일 무관) — 맞다.
- Measurement 10의 $D_c$ 행과 상관 0.9805, Observation 4의
  $Q_{cc}/n_c^2$ 행과 휴리스틱 0.049540 (느슨함 배수 2.494 ≈ "about 2.5"),
  Measurement 5의 지수 $0.039451,\,0.039671,\,0.039388$ 과
  $1/(2\langle\log N\rangle)=0.036038$ — 전부 결과 파일과 일치.
- §1.3의 네 항목 중 "No decay exponent for the effect is quoted" 와
  "No distributional law for $Z$ is asserted" 는 본문과 부합.

### P5
- **Proposition 6** — $N=200,1000,5000,20000$ 에서 이중합이 정확히
  $\mu(N-1)$ (`a2_identities.py` I5).
- **Measurement 2** — Z0–Z4 전부 결과 파일과 대조. 두 가지 "dyadic" 규약을
  명시하고 인쇄된 밀도가 어느 쪽인지 밝힌 것(0.3303/0.3298/0.3320 이
  $K<k\le2K$)이 맞다. $\prod_p(1-2/p^2)=0.3226341$ 대비 2.2–2.9% 초과도
  맞다 (0.3303/0.3226341=1.0238 등). $4\mid k$ 또는 $25\mid k$ 인
  2799/9999 = 28.0% 가 정확히 소멸집합이라는 것도 결과 파일과 일치.
- Note 1의 $(\log N)^{-2A-2}$ at $N=10^8,A=1$ $\approx8.7\cdot10^{-6}$
  (논문은 $8.8\cdot10^{-6}$ — 반올림).
- §6 (2)의 $J=3$ 값 0.848668 과 $J=8$ 양끝 값, "rounding down violates
  $z^J\ge x$ at $J=3$ and $J=8$ alike" — `audit_hb_weight.txt` 와 일치
  (가운데 값 문제는 F15).
- Note 4가 증거를 등급 매기고 두 개를 single-witness로 못박은 것,
  Note 5가 "a threshold chosen as an effect size is not a threshold" 를
  적은 것 — 방법론상 맞고, 본문에서 그것을 어기는 자리를 찾지 못했다.

### 교차 일관성
- 다섯 편이 공유하는 상수 정의 ($\AAA$, $\SS$) 는 P1·P2·P3에서 자구까지
  동일하다.
- P1 §6의 "The no-go of [7] costs a factor $\exp(c_1\sqrt{\tfrac12\log N})$"
  는 P2 Theorem 14의 결론과 맞고, 그것이 $\log N\log\log N$ 을 넘는다는
  추론도 맞다 ($\sqrt{\log N}$ 이 $\log\log N$ 을 압도).
- P2·P3·P4·P5가 서로에게 귀속시키는 내용은 (F12의 번호 문제를 빼면)
  가리켜진 진술과 실제로 일치한다. 특히 P5의 "Dividing by $W(N)$ where
  the definition calls for $V(N)$ changes the answer by exactly $\AAA(N)$
  [P3]" 는 P3 Note 2와 맞는다 (엄밀히는 $\AAA(N)(1+o(1))$).
- 기호 충돌은 F16 하나뿐. 논문 사이에서 같은 기호가 다른 뜻인 자리는
  찾지 못했다 ($H$, $C$, $V$, $Z$, $D_c$, $b$, $T_w$ 전부 일관).

---

## 안 봄

통과와 구별한다.

- **`code/` 의 30개 스크립트를 코드로서** 읽지 않았다. 결과 파일을 읽고
  중요한 양은 독립 구현으로 다시 계산했다. 코드 감사는 pass5의 일이다.
  (예외: F14의 의존성 확인을 위해 `audit_cn_coin_deep.py` 의 입출력
  부분만 봤다.)
- `lib/goldbach/` 전체.
- PDF 조판 — F13의 두 곳과 상호참조 번호 외에는 보지 않았다.
- P5 Conjecture 3을 수학 진술로서. (추측이고 증거 등급이 Note 4에
  명시되어 있으며, 그 증거의 두 팔이 패킷 밖이다.)
- P4 Measurement 3의 몬테카를로 재실행.
- 각 논문의 참고문헌 서지사항(권/쪽/연도)의 정확성.

---

## 다음 패스가 볼 자리 (우선순위)

1. **[HL] 원문 대조.** 이 패스가 구조적으로 볼 수 없었던 가장 큰 덩어리이고,
   P1 §7이 남의 논문에 대한 결함 주장이므로 가장 비싸다.
2. **P2 Theorem 14를 $T_w$ 를 명시한 채 다시 도출** (F1). 결론은 살
   것으로 보이나, 도출이 순환하지 않는지는 다시 써 봐야 안다.
3. **인쇄된 수의 분자·정의 대응 전수 점검** (F3, F4가 같은 종류다).
   "이름이 같고 정의가 다른" 실패가 이 패킷에서 두 번 나왔다.
4. **P5의 열여섯 개 미확인 판정**. 지금 상태로는 지도의 대부분이
   패킷 밖 근거에 걸려 있다.
