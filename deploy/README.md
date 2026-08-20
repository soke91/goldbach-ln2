# deploy — 배포용 최정본

`v2/paper/` 두 문서(14,129줄)를 **투고 가능한 논문 5편**으로 정리한 것.
저장소가 보존하는 정정 이력·실패한 사전등록·철회된 수치는 전부 제거했다.

```
README.md      이 파일
SPEC.md        논문 작성 규약 — 양식, 규칙, 넣을 것, 뺄 것
INVENTORY.md   전체 프로젝트 분석 — 진술 37개, 논문 주제 5편
MANIFEST.md    논문 ↔ 진술 ↔ 코드 ↔ 결과 대응표
papers/        P1..P5 (.tex, 각각 독립 컴파일)
code/          논문이 인용하는 스크립트 30개
results/       그 결과 파일 30개
lib/goldbach/  스크립트가 쓰는 공용 모듈
```

---

## 논문 5편

| 파일 | 제목 | 성격 | 진술 |
|---|---|---|---|
| `P1-mobius-fixed-class.tex` | 고정 잉여류 위 뫼비우스 가중 상관합의 무조건 유계 | 무조건 정리 2 + 항등식 | 12 |
| `P2-no-go-divisor-switch.tex` | 어떤 가중치도 $C(N)$ 을 추출하지 못한다 | 음성 정리 | 14 |
| `P3-wall-second-moment.tex` | 벽의 정확한 2차 모멘트, Chowla 가 안 되는 이유 | 정확한 사실 + 장애물 | 5 |
| `P4-coherent-cell-floor.tex` | 개수는 오차 막대가 아니다 | 방법론 (타 분야 전용 가능) | 5 |
| `P5-negative-map.tex` | 음성 지도 — 17건 사전등록 폐쇄 | 실험수학 보고 | 3 + 추측 1 |

**핵심 결과 셋**

1. **Theorem 1 (P1)** — $w_k=1$ 가지의 상관합은 무조건 $\ll_A N(\log N)^{-A}$.
   Huang–Li 의 $E_4$ 소비는 불필요하고, 요구 전체가 스칼라 $E_3$ 하나로 붕괴한다.
2. **Theorem 3 (P1)** — 그 스칼라는 Huang–Li 식 (22) 와 **항등적으로 동치**다.
   추정이 아니라 항등식 수준의 폐쇄이며, 원인은 $\mu*\log=\Lambda$.
3. **Theorem 11–12 (P2)** — 그 사이의 설계공간은 비어 있고, 완전
   Elliott–Halberstam 을 줘도 비어 있다. 손실 $\exp(c_1\sqrt{\log N/2})$ 의
   $1/2$ 은 문자 그대로 $\sqrt N$ 장벽의 지수다.

**골드바흐를 향한 순진전: 0.** 모든 논문이 초록과 §"What is not claimed"
에서 이를 명시한다.

---

## 컴파일

```bash
cd papers
pdflatex P1-mobius-fixed-class.tex && pdflatex P1-mobius-fixed-class.tex
```

외부 스타일 파일 없음. `amsmath`·`amssymb`·`amsthm`·`geometry`·`hyperref`
(+ P5 는 `longtable`) 만 쓴다. 참고문헌은 `thebibliography` 내장이라
BibTeX 불필요. arXiv 업로드에 그대로 올라간다.

## 재현

```bash
python code/<이름>.py > results/<이름>.txt 2>&1; echo $?
```

각 스크립트는 독립 실행되고, 자기 사전등록 판정을 출력하며, 실패 시
비-0 으로 종료한다. Python + numpy, 노트북 한 대면 된다.

---

## 배포본이 원본과 다른 점

| | 저장소 `v2/paper/` | 배포본 `deploy/papers/` |
|---|---|---|
| 형식 | 마크다운 | LaTeX |
| Remark | 218 | 0 (필요한 것은 본문·Note 로 흡수) |
| 정정 이력 | 보존 | **전부 삭제** |
| 실패한 사전등록 규칙 | 규칙 번호까지 기록 | 싣지 않음 |
| 철회된 수치 | 철회 사실과 함께 기록 | 인쇄하지 않음 |
| 진술 등급 | 전부 `Proposition` | 증명 유무로 `Proposition`/`Observation`/`Measurement` 분리 |
| 기호 | $A$ 가 두 뜻 (자유 지수 / 국소인자) | 국소인자를 $\mathfrak{A}$ 로 분리 |

삭제 기준의 전문은 `SPEC.md` §3. **정정처럼 보이지만 지우면 논문이
틀리는 7개 항목**은 `SPEC.md` §3.4 에 열거돼 있고, 전부 배포본에
살아 있다.

## 원본과의 관계

배포본은 원본의 요약이 아니라 **한 축(역사)을 없앤 사영**이다.
"어떻게 알게 됐는가"는 `v2/` 에 그대로 있고, `git log` 와
`v2/gate/gate.py` (77개 검사, 현재 failures: 0) 가 그것을 지킨다.
