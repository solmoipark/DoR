# 데모: 가상 소성점토 CC-Gwangju-2026 (실제 xGEMS, TINN_v4, 2026-09-03)

입력: `templates/example_calcined_clay.yaml` (값은 가정). 배합 OPC 70 / 소성점토 30, w/b 0.50, 20 °C.

## A. 예측 → 열역학 봉투 (`dorgems envelope --input … --real`)
- 역할 calcined_clay → 슬롯 metakaolin. OOD: 학습 표본 4편(sparse_role), 마할라노비스 p99 초과, 유사배합 1편(5배합).
  베이지안은 이 역할을 `other`로 풀링 → 구간이 넓다.
- DoR (q05/q50/q95 %): 1 d 5/14/23 · 7 d 10/28/47 · 28 d 14/39/66 · 90 d 16/44/76 · 365 d 17/46/79.
- 실행 3회 모두 정식 재료 주입(native), 자기검증 통과. 28 d 봉투: pH 13.4–13.5, CSHQ 43–51 g,
  **CH 8.6 g (q05) / 0 (q50, q95)**, 스트라틀링자이트 0 / 7.7 / 28.9 g, 공극률(q50) 0.362.
  → 중앙값 반응도(39 %)만 돼도 열역학적으로 CH가 소진되고 스트라틀링자이트가 생긴다.

## B. 관측 대조 (`dorgems compare --input … --real`)
| 물리량 | 재령 | 관측 | 모델(보정 후) | z |
|---|---|---|---|---|
| 결합수 | 7 / 28 / 90 d | 13.0 / 16.5 / 18.0 | 10.2 / 15.5 / 17.3 | −0.8 / −0.3 / −0.2 |
| 화학수축 | 28 d | 0.058 | 0.062 | +0.35 |
| CH | 7 / 28 / 90 d | 9.5 / 6.0 / 4.5 | 0 / 0 / 0 | −1.8 / −1.2 / −0.9 |
결합수·화학수축은 보정(×0.60) 후 일치; CH는 모델이 완전 소진(보조 물리량). 물리량당 n<5라 판정은 insufficient_data.

## C. 역추정 (`dorgems infer --input … --real --alpha-grid 11`)
- xGEMS 93회, ESS 1680, KL 0.11 (관측이 사전을 약간만 좁힘).
- 사후 DoR %: 7 d 15/30/44 · 28 d 21/42/61 · 90 d 24/49/70 (사전 28 d 14/39/66).
- 결합수·화학수축은 사후예측 안에 들어옴; CH는 어떤 α에서도 0이라 맞출 수 없음(가중 0.25).
- 단조성 플래그: 7 d에서 결합수가 α에 대해 감소(C-S-H 물이 저수화물로 대체되는 구간) — 해석 시 주의.
- 상태: 시나리오 C는 `consistency_check`, 슬롯 metakaolin은 `consistency_check`.

원자료: 이 폴더의 prediction.json, envelope.csv, comparison.csv/json, inference.json, ppc.csv, forward_map_report.md.
