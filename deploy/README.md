# deploy — 논문 다섯 편

`v2/paper/` 두 문서를 투고 형태의 논문 다섯 편으로 사영한 것. 저장소가
보존하는 정정 이력·실패한 사전등록·철회된 수치는 싣지 않는다.

```
papers/   P1..P5 (.tex, 각각 독립 컴파일)
pdf/      컴파일 산출물
code/     논문이 인용하는 스크립트 30개
results/  그 결과 파일 30개
lib/      스크립트가 쓰는 공용 모듈
```

> **상태.** 투고도 발표도 하지 않았다. 재검증 세 패스 중 둘이 끝났다 —
> 수학(`pass4`)과 코드(`pass5`). **셋째(`pass6`, 사영 감사: 줄이면서
> 떠받치던 것이 빠졌는지)는 아직 안 돌았다.** 무엇이 나왔고 무엇이
> 바뀌었는지는 `../v2_verify/README.md`.

---

## 다섯 편

| 파일 | 무엇 | 성격 | 진술 |
|---|---|---|---|
| `P1-mobius-fixed-class.tex` | 고정 잉여류 위 뫼비우스 가중 상관합의 무조건 유계 | 무조건 정리 둘 + 항등식 | 12 |
| `P2-no-go-divisor-switch.tex` | 어떤 가중치도 $C(N)$ 을 추출하지 못한다 | 음성 정리 | 14 |
| `P3-wall-second-moment.tex` | 벽의 정확한 2차 모멘트, 그리고 Chowla 가 안 되는 이유 | 정확한 사실 + 장애물 | 5 |
| `P4-coherent-cell-floor.tex` | 개수는 오차 막대가 아니다 | 방법론 | 5 |
| `P5-negative-map.tex` | 음성 지도 | 실험적 보고 | 3 + 추측 1 |

**핵심 셋**

1. **P1 Theorem 1** — $w_k=1$ 가지는 무조건 $\ll_A N(\log N)^{-A}$.
   Huang–Li 의 $E_4$ 소비가 불필요해지고, 요구 전체가 스칼라 $E_3$ 하나로
   붕괴한다.
2. **P1 Theorem 3** — 그 스칼라는 Huang–Li 식 (22) 와 **항등적으로 동치**다.
   추정이 아니라 항등식 수준의 폐쇄이고, 원인은 $\mu*\log=\Lambda$.
3. **P2 Theorem 11–12** — 그 사이 설계공간은 비어 있다. 하드 입력을 공짜로
   줘도, 완전 Elliott–Halberstam 을 줘도 비어 있다. 손실
   $\exp(c_1\sqrt{\log N/2})$ 의 $1/2$ 은 문자 그대로 $\sqrt N$ 장벽의
   지수다.

**골드바흐를 향한 순진전: 0.** 다섯 편 모두 초록과 "What is not claimed"
에서 이를 명시한다.

또한 P1 §7 이 [HL] 식 (18) 의 결함 하나를 보고한다 — $n$-의존 제약이
$n$-무관 상한으로 바뀌면서 누락된 항 $\Delta$. arXiv v1·v2 원문 대조로
확인했고, 그들의 Theorem 1 과 Corollary 1 은 그대로 성립한다.

---

## 컴파일

```bash
cd papers
pdflatex P1-mobius-fixed-class.tex && pdflatex P1-mobius-fixed-class.tex
```

TeX Live(scheme-small) 에서 다섯 편 전부 2회 통과, **오류 0 · 미해결 참조 0 ·
경고 0 · overfull 0**. 쪽수 P1 14 · P2 12 · P3 7 · P4 6 · P5 8.

외부 스타일 파일 없음. `amsmath`·`amssymb`·`amsthm`·`geometry`·`hyperref`
(+ P5 는 `longtable`) 뿐이고, 참고문헌은 `thebibliography` 내장이라 BibTeX
불필요하다. 남은 것은 underfull hbox 뿐이며 지면에 드러나지 않는다.

## 재현

```bash
python code/<이름>.py > results/<이름>.txt 2>&1; echo $?
```

각 스크립트는 독립 실행되고, 자기 사전등록 판정을 출력하며, 실패 시 비-0
으로 종료한다. Python 과 numpy, 노트북 한 대가 전부다.

---

## 원본과의 관계

배포본은 요약이 아니라 **한 축을 없앤 사영**이다. 없앤 축은 역사 —
어느 판본이 무엇을 틀렸고 어떤 사전등록이 깨졌는지. 그 축은 `v2/paper/`
에 그대로 있다.

| | `v2/paper/` | `deploy/papers/` |
|---|---|---|
| 형식 | 마크다운 | LaTeX |
| Remark | 218 | 0 — 필요한 것은 본문과 Note 로 흡수 |
| 정정 이력 | 보존 | 싣지 않음 |
| 철회된 수치 | 철회 사실과 함께 보존 | 인쇄하지 않음 |
| 진술 등급 | 전부 `Proposition` | 증명 유무로 `Proposition` / `Observation` / `Measurement` 분리 |
| 기호 | $A$ 가 두 뜻 | 국소인자를 $\mathfrak{A}$ 로 분리 |

그 사영이 옳았는지 — 줄이면서 문장을 떠받치던 것을 두고 오지 않았는지 —
는 `pass6` 이 판정한다. 아직 안 돌았다.
