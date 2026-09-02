# 커널 조사: OPC 페이스트의 포틀란다이트 과대와 황산염 상 부재 (2026-09-03)

G2-3(실기 OPC 참조 검사)에서 모델 CH가 문헌보다 크게 높았던 원인을 InverseGems 커널
(`e84d7a9`) + TINN_v4 GEMS3K 시스템에서 직접 추적한 기록. 재현: `scripts/`는 없고 아래
계산은 `dorgems.gems.forward.run_forward(capture_species=True)`로 수행(캡처 JSON에
`phase_elements_amounts` 포함).

## 관찰 (OPC 100 %, w/b 0.5, 20 °C, 기본 kinetics)

| 재령 d | PK α C3S/C2S/C3A/C4AF | CH g/100 g | CSHQ g (Ca/Si) | 기타 상 | S in system |
|---|---|---|---|---|---|
| 1 | 0.13/0.07/0.11/0.05 | 3.8 | 7.8 (1.63) | C3AH6 1.6, C3(AF)S0.84H 0.7 | **0** |
| 7 | 0.48/0.26/0.41/0.21 | 14.5 | 29.5 (1.63) | C3AH6 6.2, 하이드로가넷 2.7 | 0 |
| 28 | 0.82/0.60/0.81/0.53 | **25.3** | 51.6 (1.63) | C3AH6 12.2, 하이드로가넷 6.8 | 0 |
| 90 | 0.90/0.72/0.88/0.78 | 28.1 | 56.3 (1.63) | C3AH6 13.2, 하이드로가넷 10.0 | 0 |
| 365 | 0.94/0.78/0.91/0.84 | 29.4 | 58.9 (1.63) | C3AH6 13.8, 하이드로가넷 10.7 | 0 |

xGEMS 입력(28 d): `C3S 54.27 g, C2S 5.01 g, C3A 8.70 g, C4AF 3.87 g, H2O 50 g` — **그 외 아무것도 없음**.

## 원인

1. **OPC의 소량 산화물이 입력에서 사라진다.** `xgems_input_builder.py:159-176`은 OPC를 Bogue 4상으로만
   바꾸고, 산화물 조성의 SO3 2.6 %, MgO 2.0 %, Na2O 0.5 %, K2O 0.7 %(합 5.8 %)와 Bogue 합(92.6 %)
   밖의 7.4 %를 어디에도 넣지 않는다. 결과:
   - 시스템에 S가 0 → 에트린자이트·모노설페이트가 생길 수 없고 C3A는 C3AH6로 간다.
   - Na/K가 0 → 공극수 pH가 포틀란다이트 완충값 **12.662로 고정**(모든 배합·재령에서 동일). 실제 OPC
     공극수는 13.2–13.6. GemsPilot의 pH 앵커 12.661534도 같은 결함의 산물이다.
   - 질량 수지: 결합재 100 g 중 7.4 g이 반응도 미반응도 아닌 채 소실(`unreacted_masses` 20.75 g은
     Bogue 4상만 집계).
2. **기본 OPC 조성이 C3S 과다.** `configs/materials.yaml`의 OPC(CaO 65.9, SiO2 20.2)는 Bogue로
   C3S 66.2 % / C2S 8.3 %. 통상 CEM I는 C3S 55–65 % / C2S 12–20 %.
3. **CSHQ의 Ca/Si 상한 1.63.** Cemdata18 CSHQ 고용체는 Ca/Si ≈ 1.67이 상한이라 잉여 Ca는 전부 CH로
   간다. 포틀란다이트 공존 시 실측 C-S-H Ca/Si는 1.7–1.85.

CH 수지 검산(28 d): C3S 54.27 g = 0.238 mol × (3 − 1.63) = 0.326 mol + C2S 0.029 mol × 0.37 =
0.011 mol → 0.337 mol × 74.09 = **25.0 g** (계산 25.3 g). 즉 CH 값은 위 세 가정의 산술적 귀결이다.

## 대조 계산 (실제 xGEMS, 28 d / 90 d)

| 케이스 | CH 28 d | CH 90 d | 황산염 상 28 d | pH |
|---|---|---|---|---|
| 기본 OPC 100 | 25.3 | 28.1 | 없음 | 12.662 |
| 기본 OPC 95 + 석고 5 (별도 성분) | 23.8 | 26.4 | AFm(OH_SO4/SO4_OH) 9.1 g, AFt 1.2 g | 12.662 |
| CEM I형 조성(CaO 63.5, SiO2 20.5, Al2O3 5.0, Fe2O3 3.0, SO3 3.0) 95 + 석고 5 | **21.8** | 24.3 | AFm 6.2 g, AFt 4.7 g | 12.662 |

문헌(OPC, w/b 0.4–0.5, 28 d, TGA): CH 15–22 g/100 g 시멘트가 흔함. 모델은 현실적 조성에서도
상단(21.8)에 있고, 나머지 ~3 g은 CSHQ Ca/Si 상한(1.63 vs 1.75–1.85 → 0.24 mol × 0.15 × 74 ≈ 2.7 g)으로
설명된다. pH는 조성과 무관하게 고정되므로 알칼리 결손이 확정적이다.

## DoR 에이전트에 미치는 영향

- 시나리오 B의 CH 대조는 커널 결함(1·2·3)과 DB basis 문제가 겹쳐 있어 현재로선 판정에 쓸 수 없다.
  결합수·화학수축 대조는 이 결함의 영향이 작다(결합수 24.0 g은 문헌 범위).
- 시나리오 C에서 CH를 우도에 넣으면 α가 과대 추정되는 방향으로 편향된다(모델 CH가 높을수록 같은
  관측 CH를 맞추려면 더 큰 α가 필요). G3-3 실측 결과 해석 시 반드시 고려.
- twin(G2-5)에서 측정 DoR을 고정했을 때 CH가 0으로 소진되는 현상은 이 결함과 방향이 반대(모델 CH
  과대)이므로, SCM 배합에서는 별도 원인(측정 DoR의 열역학적 과대, CH 가용성)이 더 크다.

## 커널 패치 제안 (P-IG-6, 별도 PR)

1. `build_xgems_input`: OPC 산화물 조성의 SO3를 `Gp`(또는 `Anh`) 종으로, 반응분율 = 1(황산염은
   초기 용해)로 추가. Bogue 식이 이미 C3S에서 2.85·SO3를 빼므로 Ca 수지는 일관됨.
2. Na2O·K2O를 알칼리 종(시스템에 `Na2O`, `K2SO4`, `Na2SO4` 존재)으로 추가, 반응분율 1 — 공극수 pH·
   C-S-H 알칼리 흡수 복원.
3. MgO를 `Brc`(brucite) 전구체 또는 periclase로 추가(반응분율은 클링커 평균 α로 근사).
4. `materials.yaml` OPC 기본값을 CEM I 대표 조성(C3S ≈ 58 %)으로 교체하거나 `opc_phase_mass_percent`
   (이미 인자 존재)를 recipe/forward_query에서 받을 수 있게 노출.
5. 회귀 앵커(GemsPilot `agent_qa_generated.yaml`, DoRGems `docs/real_anchors_TINN_v4.json`)는 패치 후
   재생성해야 한다 — 현재 앵커는 결함이 있는 화학을 고정하고 있다.

CSHQ Ca/Si 상한은 열역학 DB의 성질이라 커널 패치 대상이 아니다. 비교 시 σ_model에 +2–3 g의 계통
편향으로 반영하거나, CH 대신 결합수·화학수축을 1차 검증량으로 쓰는 것이 현실적이다.

## 패치 P-IG-6 적용 결과 (InverseGems `dorgems/p-ig-6-opc-minor-oxides`, 커밋 `535f90e`)

구현: `reaction_parameters.py`에 `opc_minor_oxides` 정책(기본 enabled; SO3 1.0 → CaSO4로, 동반 CaO
재투입; Na2O·K2O 1.0; MgO = 클링커 평균 α), `xgems_input_builder._add_opc_minor_oxides`, source ledger
행 추가, 시그니처 payload 포함, 테스트 5개(`tests/test_opc_minor_oxides.py`). InverseGems 스위트
226 passed / 1 failed(기존 `test_feature_table`).

실제 xGEMS(TINN_v4), OPC 100 w/b 0.5:

| 재령 d | pH (전 → 후) | CH g (전 → 후) | 황산염 상(후) | 기타(후) |
|---|---|---|---|---|
| 1 | 12.66 → **13.10** | 3.8 → 4.7 | 에트린자이트 6.3 g, 석고 0.9 g | 브루사이트 0.3 |
| 7 | 12.66 → **13.67** | 14.5 → 15.1 | AFt 10.5, AFm 2.5 g | 브루사이트 1.2 |
| 28 | 12.66 → **13.60** | 25.3 → 26.1 | AFm 9.4, AFt 2.3 g | 브루사이트 2.3, 하이드로가넷 6.8 |
| 90 | 12.66 → **13.58** | 28.1 → 28.9 | AFm 10.6, AFt 0.8 g | 브루사이트 2.5 |

OPC 60 / 슬래그 40 (w/b 0.45, 28 d): pH 13.39, CH 9.6 g, AFm 9.0 g, MgAl-OH-LDH 5.2 g(전에는 없던
하이드로탈사이트가 생김). 질량 수지 100.2 g(전 92.6 g). porosity 0.412 → 0.327(OPC), 0.386 → 0.326(슬래그).
결합수 24.0 → 29.1 g(AFt/AFm의 결정수), 화학수축 0.076 → 0.058 mL/g.

**CH는 패치로 줄지 않는다(+0.8 g).** CaSO4의 Ca가 추가되고 알칼리 흡수로 C-S-H Ca/Si가 1.63 → 1.58로
낮아지기 때문. G2-3(등급 A 39건/11편) median r: 12.4 → 13.3 g. 따라서 CH 과대의 남은 원인은
(2) 기본 OPC 조성의 C3S 66 %와 (3) CSHQ Ca/Si 상한이며, 이 둘은 각각 `materials.yaml` 기본값 교체
(또는 `opc_phase_mass_percent` 노출)와 열역학 DB 선택의 문제다. 패치의 실효는 pH·황산염 상·Mg 상·
질량 수지의 정상화이고, 이는 QXRD 상 비율 비교(§8.3)와 공극률에 직접 영향을 준다.

앵커 갱신: DoRGems `docs/real_anchors_TINN_v4.json`(패치 후) / `_pre_pig6.json`(패치 전); GemsPilot
`dorgems/p-gp-4-anchor-refresh`(`a05467b`)에서 mock 앵커 재생성 + `agent_qa_generated_TINN_v4.yaml`
(실기 15개) 추가. GemsPilot의 `Test-dat.lst` 실기 앵커 18개는 그 시스템이 있는 PC에서 재생성 필요.

## 머지 후 상태 (2026-09-03)

- InverseGems 로컬 `master`: `c4600b2`(P-IG-6 머지) → `50beefb`(기본 OPC를 CEM I 42.5 대표 조성
  CaO 63.3 / SiO2 20.2 / Al2O3 5.0 / Fe2O3 3.0 / MgO 1.8 / SO3 3.0 / Na2O 0.2 / K2O 0.8, Bogue C3S 57.7 %) →
  `5969d38`(테스트). **원격 미푸시.** 스위트 226 passed / 1 failed(기존).
- GemsPilot 로컬 `main`: `2924b96`(P-GP-4 머지) — mock 앵커 18개 + TINN_v4 실기 앵커 15개(`agent_qa_generated_TINN_v4.yaml`)
  를 머지 커널로 재생성. **원격 미푸시.** 48 passed / 1 failed(기존 하드코딩 경로).
- DoRGems 실기 앵커 `docs/real_anchors_TINN_v4.json` 재갱신: OPC100 w/b 0.5 28 d → pH 13.534, porosity 0.343,
  CH 23.5 g, 결합수 28.4 g, 화학수축 0.064 mL/g; OPC60/slag40 → pH 13.322, porosity 0.348, CH 6.7 g.

G2-3(등급 A 39건/11편) 커널별 비교:

| 커널 | 모델 CH 중앙 | median r | 판정 |
|---|---|---|---|
| base `e84d7a9` | 24.8 | +12.4 | tension |
| P-IG-6만 | 25.7 | +13.3 | tension |
| **머지(P-IG-6 + CEM I)** | **23.2** | **+11.5** | tension |

조성 교체로 CH가 2.5 g 내려왔지만 문헌 중앙(11.8)과의 차이는 여전히 크다. 논문별로는 3편이 ±2.4 g 안에
들고 8편은 +8~+20 g. 남은 후보: CSHQ Ca/Si 상한(≈ 3 g), PK 28 d C3S α 0.82의 과대 가능성, 그리고 무엇보다
등급 A로 남은 39건 자체의 basis 신뢰성(11편 중 3편만 모델과 정합). 사용자 결정에 따라 DB basis 재판정은
보류하고, 시나리오 B의 1차 검증량을 결합수·화학수축·QXRD 상 비율로 옮기는 방향을 검토한다.
