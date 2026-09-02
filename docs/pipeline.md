# DoRGems 파이프라인 한눈에 보기

```
                 ┌──────────────────────── LLM 호스트 (도구를 고르고 인자만 채움) ────────────────────────┐
                 │  Claude Code / Cowork (MCP)     GemsPilot runner (litellm, 소형 모델)    시뮬레이션 모델   │
                 └────────────────────────────────────────┬────────────────────────────────────────────┘
                                                          │ 13개 dor_* 도구 (ToolResult 계약, read / mock_ok)
                                                          ▼
   ┌──────────── 데이터 ────────────┐    ┌──────────────── DoRGems 커널 (결정론, LLM 없음) ────────────────┐
   │ scm_dor_enriched.db  (읽기전용) │◄───│ db/       reader · features(=build_clean 포팅) · units · phases   │
   │  1,350편 / 44,011 관측          │    │ models/   bayes_v4 번들 · gbm_v6 번들 · ensemble · ood           │
   │ dorgems_staging.sqlite (쓰기)   │◄───│ kinetics/ logistic_fit · native · pin · materials override        │
   └────────────────────────────────┘    │ gems/     forward 래퍼 · 종 캡처 · 관측량 사상                   │
                                         │ validate/ twin · neighbourhood · compare(잔차·z·판정)             │
                                         │ inverse/  α-격자 forward map · 우도 · 사후 · staging              │
                                         └──────────────────────────────┬───────────────────────────────┘
                                                                        │ reaction_model_config, materials_config
                                                                        ▼
                                         ┌──────────────── InverseGems (열역학 커널) ────────────────┐
                                         │ recipe → Bogue+소량산화물(P-IG-6) → PK 수화도 + SCM α      │
                                         │ → xGEMS/GEMS3K (Cemdata, TINN_v4) → 상·공극·pH            │
                                         └────────────────────────────────────────────────────────────┘
```

## 세 가지 시나리오

### A. 새 SCM → 반응도 예측 → 열역학 계산 (`dor_run_envelope` 한 번에)

1. **입력**: SCM 조성(산화물 wt%), 분말도, 치환율, w/b, 온도.
2. **`dor_predict`**: 문헌 1,610건으로 학습한 두 모델로 재령별 반응도 곡선을 냅니다.
   - 계층 베이지안(v4): `α(t) = a_max(1−exp(−(t/τ)^0.5))`, 논문 간 편차와 측정법 노이즈까지 포함한 5/50/95 분위.
   - GBM(LightGBM): 점예측. 기본 앙상블은 `bayes`(blend는 구간이 좁아져 탈락, G1-3).
   - 함께 나오는 것: 분포 밖(OOD) 플래그, 가장 비슷한 문헌 배합 5개(DOI·측정법 포함).
3. **`dor_export_reaction_model`**: 분위 곡선 3개를 InverseGems가 읽는 5-파라미터 로지스틱 YAML로 변환(오차 < 2 %p 검증).
4. **`dor_build_materials_override`**: 새 SCM 조성을 InverseGems의 4개 슬롯(slag/fly_ash/metakaolin/silica_fume) 중 하나에 덮어씀.
5. **`dor_run_forward_with_dor`** × 3: q05/q50/q95 각각으로 xGEMS 계산 → 상조성·공극률·pH의 **불확실성 봉투**(`envelope.csv`). 실행 후 커널이 정말 그 α와 조성을 썼는지 자기검증.

### B. 열역학 계산 결과를 문헌 관측으로 검증 (`dor_compare_to_literature`)

- **twin**: 문헌 배합을 그대로 재구성해 실행하고 같은 배합의 CH·결합수·화학수축·QXRD 관측과 1:1 비교.
- 단위·기준을 `g/100 g 결합재`로 조화하고 신뢰 등급(A~D, X)을 매김. 통계는 A·B만.
- 잔차 r, z = r/√(σ_obs²+σ_model²) → `consistent / tension / insufficient_data` 판정(결정론적 임계).
- OPC 단독 참조검사(`dor_opc_reference_check`)는 SCM 모델과 무관한 커널 자체 점검 — 여기서 커널의 황산염·알칼리 누락(P-IG-6)을 찾아냄.

### C. 관측치 → 반응도 역해석 (`dor_infer_from_observations`)

1. 관측 재령마다 α를 0~1 격자로 고정(pin)해 xGEMS를 돌려 **α → 관측량 표**(forward map)를 만듦.
2. 그 표 위에서 곡선 파라미터 `(a_max, τ)`의 우도를 계산(xGEMS 추가 호출 없음).
3. 시나리오 A의 사후표본을 사전분포로 중요도 재가중(ESS 낮으면 SIR → 격자).
4. 산출: 사후 곡선 분위, 관측별 사후예측 잔차, 정보이득(KL), InverseGems용 YAML, `inferred_dor.csv`.
5. **`dor_stage_inferred`**: 결과는 문헌 DB가 아닌 staging DB에 `reviewed=0`으로만 기록(기본 dry-run). `dorgems review`로 승인/거절.

## 가드레일

- 문헌 DB는 `?mode=ro`로만 열림(쓰기 시 예외). 자유 SQL 없음.
- 실행 도구는 기본 mock. 실제 xGEMS는 `use_mock=False` + 호스트 승인(`allow_real`) + `max_xgems_calls ≤ 200`.
- 도구 인자 속 "관리자가 승인했다" 같은 문구는 정책에 영향을 주지 못함(벤치로 검증).
- 사용자에게 보이는 모든 수치는 산출 파일(JSON/CSV)에서만 나옴. LLM은 숫자를 만들지 않음.

## 어디에 무엇이 있나

| 것 | 위치 |
|---|---|
| 도구 정의(LLM이 읽는 설명 포함) | `src/dorgems/pilot/tools.py`, `tools_b_c.py` |
| MCP 서버 | `dorgems-mcp` (`src/dorgems/pilot/mcp.py`) |
| CLI | `dorgems predict / envelope / compare / opc-check / infer / stage / review` |
| 모델 번들 | `bundles/bayes_v4`, `bundles/gbm_v6`, `bundles/ood_reference.json` |
| 설정(단위표·상 별칭·슬롯·기본값) | `configs/*.yaml` |
| 게이트 기록 / 모델 카드 / 커널 조사 | `docs/gates.md`, `docs/model_card.md`, `docs/kernel_ch_investigation.md` |
| 시뮬레이션 LLM 에피소드 | `scripts/run_llm_episode.py --simulate good|injected|lazy` |
