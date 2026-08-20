# INVENTORY — 전체 프로젝트 분석

`z:\업무\goldbach-ln2-real` 전체를 훑어 **무엇이 있는지**, **논문 몇 편이
나오는지**, **정리·증명이 몇 개인지** 정리한다.

- 대상: `v1/`(동결), `v1_verify/`(구조만), `v2/`(현행 산출물),
  `lib/goldbach/`(공용 모듈)
- 산출물의 유일한 집: `v2/paper/theorem_A.md`(7,065줄),
  `v2/paper/wall_v3.md`(7,064줄)
- 게이트 상태: `v2/gate/gate.txt` → **failures: 0** (77개 검사)

---

## 1. 규모

| 항목 | 수 |
|---|---|
| 산출물 문서 | 2 (합 14,129줄) |
| 계산 스크립트 (`v2/code/`) | 213 |
| 결과 파일 (`v2/results/`) | 209 |
| v1 스크립트 / 결과 | 137 / 91 |
| 독립 재검증 패스 | 3 (`v2/verify/pass1–3`) |
| 게이트 검사 | 77 (G1–G77) |
| 증거 마커(`evidence:`) | 243 (고유 스크립트 207) |

---

## 2. 진술 개수 — 정확한 집계

번호 붙은 진술(Remark 제외) **총 37개**. 두 문서가 같은 라벨을 재진술하는
곳(`wall_v3` §Results imported)은 중복으로 세지 않았다.

| 종류 | 개수 | 라벨 |
|---|---|---|
| **Theorem** | 4 | `thm:A`, `thm:C`, `thm:D`, `thm:Dprime` |
| **Corollary** | 1 | `cor:B` |
| **Proposition** | 18 | `prop:onesided` `prop:nolog` `prop:direct` `prop:untrunc` `prop:layers` `prop:combined` `prop:flatsum` `prop:MT` `prop:Dpp` `prop:E` / `prop:V` `prop:W` `prop:coindisc` `prop:dilate` `prop:posweights` `prop:coh` `prop:placebo` `prop:scaleinv` |
| **Lemma** | 13 | `lem:complete` `lem:degen` `lem:BV` `lem:density` `lem:mu` `lem:completelog` `lem:Gb` `lem:extract` `lem:bv` / `lem:MP` `lem:coin` `lem:placebo` `lem:cellmom` |
| **Conjecture** | 1 | `conj:L` |
| **Remark** | 218 | (105 + 113) — 실험실 노트 계층, 배포본에서 제거 |

### 2.1 증명 현황

| 구분 | 개수 |
|---|---|
| 명시적 `**Proof.**` 블록 | 21 |
| 절 단위 증명 (`thm:A` 7단계, `cor:B`, `thm:C`) | 3 |
| 본문 내 유도로 증명 완결 | 10 |
| **증명된 진술 합계** | **34** |
| Observation 등급 (유도는 있으나 상수 미통제) | 1 — `prop:coh` |
| Measurement 등급 (순수 측정 + 사전등록 널) | 1 — `prop:placebo` |
| Conjecture | 1 — `conj:L` |

> **배포본에서 반드시 고칠 것**: 저장소는 `prop:coh`·`prop:placebo` 를
> 다른 34개와 같은 `Proposition` 으로 부른다. 등급이 다르므로 이름도
> 달라야 한다 (SPEC §2, 규칙 T1).

### 2.2 무조건적인 것 / 조건부인 것

- **무조건적 (가정 없음)**: `thm:A`, `cor:B`, `thm:C`, `thm:D`,
  `prop:E`, `prop:V`, `lem:MP`, `prop:dilate`, `prop:posweights`,
  `prop:untrunc`, `prop:layers`, `prop:combined`, `prop:flatsum`,
  `lem:cellmom`, `lem:coin`, `lem:placebo`, `prop:coindisc`,
  그리고 §Δ의 결함 보고 — **18개**
- **조건부**: `thm:Dprime`($EH(N^{\theta_E})$ 가정),
  `prop:onesided`·`prop:nolog`·`prop:direct`(골드바흐 충분조건 형태)
- **음성(no-go)**: `thm:D`, `thm:Dprime`, `prop:Dpp`, `prop:E`
- **골드바흐를 향한 순진전**: **0**. `thm:A` 는 요구의
  골드바흐-중립 절반만 제거한다.

---

## 3. 논문 주제 — 몇 편이 나오는가

**핵심 답: 자립 가능한 논문 5편.** 세분하면 7편까지 가능하나, 4·5편을
합치는 편이 낫다고 판단했다. 아래 표의 P1–P5 가 `deploy/papers/` 에 있다.

| # | 논문 | 성격 | 주재료 | 자립도 |
|---|---|---|---|---|
| **P1** | 고정 잉여류 위 뫼비우스 가중 상관합의 무조건 유계와 Huang–Li 환원의 두 귀결 | 해석적, **무조건 정리** | `thm:A` `cor:B` `thm:C` `prop:onesided` `prop:nolog` §Δ + 보조정리 7 | ★★★ 단독 투고 가능 |
| **P2** | 나눗수 교환 경로의 no-go: 어떤 가중치도 $C(N)$ 을 추출하지 못한다 | 해석적, **음성 정리** | `thm:D` `thm:Dprime` `prop:Dpp` `prop:E` `lem:Gb` `lem:extract` `lem:bv` + 항등식 5 | ★★★ 단독 투고 가능 |
| **P3** | 벽의 정확한 2차 모멘트, 그리고 Chowla 가 그것을 통제하지 못하는 이유 | 해석적 + 계산 | `prop:V` `lem:MP` `prop:W` `lem:coin` `prop:coindisc` | ★★★ 단독 투고 가능 |
| **P4** | 개수는 오차 막대가 아니다 — 산술 필드 셀 평균의 정확한 요동 바닥 | **방법론**, 타 분야 전용 가능 | `lem:cellmom` `prop:coh` `lem:placebo` `prop:placebo` `prop:scaleinv` | ★★★ 단독 투고 가능 |
| **P5** | 음성 지도: $\mu$-쌍 필드에 이름 있는 결합면이 없다 (실험적 보고) | **실험수학 보고** | 17개 사전등록 폐쇄 + `conj:L` | ★★☆ Experimental Math 류 |

### 세분하면 나오는 추가 2편 (권장하지 않음)

| # | 주제 | 왜 안 쪼개는가 |
|---|---|---|
| P2a | 골드바흐 개수의 확대 항등식 (`prop:untrunc` `prop:layers` `prop:combined` `prop:direct` `prop:flatsum` `prop:dilate` `prop:posweights`) | 전부 유한 교환으로 얻는 항등식이라 그 자체로는 얇다. P2의 §2(설계공간)로 넣으면 no-go 의 동기가 된다. |
| P3a | Huang–Li 식 (18)의 결함과 보정항 $\Delta$ | 한 쪽짜리 정오표. P1 §6 으로 넣는 것이 정상적인 관행. |

### 3.1 논문 간 의존 관계

```
P1 (thm:A, thm:C)  ─── 근거 제공 ──▶  P2 (설계공간이 비었음)
     │
     └── C(N)=o(N) 로 환원 ──▶  P3 (벽의 정확한 사실)
                                    │
                                    └── 오차막대 문제 ──▶ P4 (방법론)
                                    │
                                    └── 측정 한계 ──▶ P5 (음성 지도)
```

P1 은 다른 어느 것도 인용하지 않고 성립한다. P2 는 P1 의 `thm:C` 를
동기로 쓰되 논리적으로는 독립(자체 보조정리로 닫힌다). P3·P4 는
P1 의 결론(`C(N)=o(N)` 로의 환원)을 배경으로만 쓴다. P5 는 P3 의
`lem:coin` 이 필요하다.

---

## 4. 무엇이 새롭고 무엇이 고전인가 (과대주장 방지)

| 진술 | 새로움 등급 |
|---|---|
| `thm:A` | **적용이 새롭다.** 메커니즘(나눗수 교환 → 짧은 변수 위 μ → BV)은 고전. 기여는 이 소비가 불필요함을 보인 것. |
| `thm:C` | **이 프로그램의 최대 결과.** $E_3$ 이 Huang–Li 식 (22)와 **항등적으로 동치**임. 추정이 아니라 항등식 수준의 폐쇄. |
| `cor:B` | `thm:A` 의 직접 귀결 (Abel 합). |
| `thm:D`/`thm:Dprime` | **장르는 고전** (Bombieri 점근 체). 새로운 것은 정량적 형태: 두 문턱이 $N^{1/2}$ 만큼 벌어진다는 것, 그리고 완전 EH 를 줘도 닫히지 않는다는 것. |
| `prop:E` | 관찰. Parseval 이 아래에서 막는다는 사실의 정확한 기록. |
| `prop:V` | Mirsky 정리 + 부분합. **국소인자가 $\mathfrak{S}$ 가 아니라 $\mathfrak{A}$** 라는 것이 요점. |
| `lem:MP` | 정확·무조건. 절단이 load-bearing (안 하면 $1.57$ 배 어긋남). |
| `prop:W` | **조건부 경로가 왜 안 닫히는지의 정량화.** 증폭 계수 $\sim N/(\mathfrak{A}\log N)$, 부호가 양이라 삼각부등식이 날카롭고, 실측 예산은 이미 초과. |
| `lem:cellmom`/`prop:coh` | **다른 분야로 전용 가능.** 산술 입력을 공유하는 필드의 셀 평균은 자기평균하지 않는다. |
| `conj:L` | 4개 스탬프 중 1개만 독립 재검증. 그렇게 적는다. |

---

## 5. 배포에서 제거되는 것 (오류·정정 계층)

저장소가 의도적으로 보존한 것들이며, 배포본에서는 전부 뺀다.
자세한 기준은 `SPEC.md` §3.3.

| 계층 | 대략 개수 | 예 |
|---|---|---|
| 버전 간 정정 서술 | 약 40문단 | "Version 3 said this agrees to six decimals … Five of the six do" |
| 실패한 사전등록 규칙 | 30+ (X5 X6 X7 C2 C3 C4 I1–I3 L3 M1 M2 M4 N4 Z1–Z6 H4 U4 B5 Y1 Y2 …) | "That was the audit's rule X6" |
| 철회된 수치 | 6 | `prop:V`의 `0.000582`, R1·R2 정밀도, E1 행의 "exact square-root cancellation", `rem:levelmeas` 레벨값, `rem:artifact` 예측, v3 증폭 행 |
| 반증된 자기 예측 | 5 | `rem:toprdom` `rem:cap` `rem:band` `rem:artifact` `rem:thetasweep`(I1–I3) |
| 취소된 추측 | 1 | $C(N)=m(N)+\sqrt{V}G(N)$ 가우시안 법칙 (증거 부실로 철회) |
| 운영 메타 | 전부 | `evidence:` 마커, G1–G77, M1–M9, 접두사 규약, `OPEN.md`/`DECISIONS.md` |
| 1인칭 실패 고백 | 10+ | "두 번 물렸다", "내 출력이 걸렸다" |

**단, 다음 7개는 정정처럼 보이지만 지우면 논문이 틀린다** (SPEC §3.4):
`rem:bound`(오차항 $N(\log N)^{-A}$), `rem:sign`(부호 $-\mathfrak{S}$),
`rem:threeway`(세 갈래), `rem:trap`($(q,N)=1$), `prop:onesided`가
`thm:D`를 재개방하지 않음, `lem:coin`의 단방향성,
`rem:secondorder`(2차 형태 $\mathfrak{A}(N\log N-N)$).

---

## 6. 열린 문제 (배포본 §Open 에 그대로 싣는다)

1. $C(N)=\sum_{n<N}\Lambda(n)\mu(N-n)=o(N)$ — **벽 자체**. 미해결.
2. $h$ 를 가로지르는 **부호 상쇄**를 공급하는 대상의 이름
   (`prop:W` 가 남기는 유일한 경로).
3. Huang–Li 소비가 요구하는 가중치가 μ 에 대해 레벨 $3/5$ 가 알려진
   triply well-factorable 계급에 드는가.
4. `conj:L` 의 나머지 3개 스탬프 독립 재검증.
5. 마스크의 감쇠 지수 (얕은 세 깊이에서는 바닥을 못 넘어 적합 불가).
6. K1 (곱셈적 Fejér 커널) — 설계가 돌린 필드에서 계산 불가였으므로 재개방.
7. C-III (Motohashi 형) — 유한 시험이 없다.

---

## 7. 판정

- **정리·명제·보조정리 37개, 그중 증명된 것 34개, 무조건적인 것 18개.**
- **배포 가능한 논문 5편** (해석 3 + 방법론 1 + 실험보고 1).
- 가장 강한 단독 산출물은 **P1의 Theorem C** — 추정이 아니라 **항등식**으로
  demand 측을 닫는다.
- 골드바흐를 향한 순진전은 **0**이며, 배포본은 이를 초록과 §1.3에
  명시한다.
