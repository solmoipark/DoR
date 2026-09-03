# DoRGems 결정 게이트 기록

스펙 §12의 게이트를 순서대로 기록한다. 숫자는 재현 스크립트/테스트에서 나온 값만 적는다.
게이트 실패는 숨기지 않는다.

## 환경 (2026-09-02)

| 항목 | 값 |
|---|---|
| OS / Python | Windows 11, conda env `dorgems` (Python 3.12.14, conda-forge) |
| InverseGems | `e84d7a9` (clone `../InverseGems`, editable) |
| GemsPilot | `753cf6d` (clone `../GemsPilot`, editable `--no-deps`) |
| 문헌 DB | `modeling/scm_dor_enriched.db` (읽기 전용, 파일 미변경) |
| xGEMS | 세션 초반에는 없었음. 이후 사용자가 GEMS3K 시스템 파일(`20260902 TINN_v4`: dch/ipm/fun/dbr, 46상·140종·13원소, Cemdata형 `CSHQ` C-S-H, CNASH 없음)을 제공 → `gems_systems/TINN_v4/`로 복사, conda-forge `xgems 2.1.2`(win-64, py312)를 `dorgems` 환경에 설치. `inverse_gems.env_check` 통과, 실기 게이트 실행 가능해짐 |
| PyMC | 5.28.0 + arviz 0.23.4 (arviz 1.x는 pymc 5.28과 비호환 → 다운그레이드). g++ 없음 → `PYTENSOR_FLAGS=mode=NUMBA,cxx=`. conda-forge MKL BLAS가 delay-load 오류(0xc06d007f)로 numpy.linalg·scipy까지 크래시 → `libblas=*=*openblas`로 교체. threadpoolctl 프로브 오류는 `scripts/env_win/sitecustomize.py`로 우회(모델·스크립트 미수정) |
| bayes 재실행 | `scripts/run_bayes_v4.py`(runpy 래퍼, 스크립트 무수정) → `modeling/work/bayes_idata.nc` 등 생성, 2.6 분 |

## M0 — 기반

| 게이트 | 결과 | 수치 | 근거 |
|---|---|---|---|
| G0-1 `build_dor_table` == `dor_scm_final.csv` | **통과** | 1,610행; obs_uid 집합·dor_pct·scm_pct·w_b·T·scm_role·system_type·method_group·resolve_how·confidence·fig_only 전부 일치. blended 1,592행/80편/476배합 재현 | `tests/test_features_golden.py::test_golden_full_db` |
| G0-2 bayes 번들 vs `bayes_role_kinetics.csv` | **조건부 통과** | 번들은 재실행 사후와 ±0.5 %p/±2 % 이내(`bayes_v4/bayes_role_kinetics.golden.csv`). **그러나 재실행 결과는 `modeling/bayes_role_kinetics.csv`(8/26, 1,177행 표 기준: n_obs 583/176/418)와 크게 다르다**: fly_ash 28.1→30.6 %, other 47.3→39.0 % (τ 13.1→5.2 d), slag 58.5→56.9 % (τ 18.5→13.3 d). 원인: 현행 `dor_scm_final.csv`(8/28, 1,610행: 755/304/551)로 재학습됐기 때문 — 스펙 §16-2의 "±0.5 %p" 전제(같은 데이터)가 성립하지 않음. 수렴 max r̂ 1.010 / min ESS 860. LOPO(5 fold): point R² 0.421 MAE 11.35, pred R² 0.419 MAE 11.54, 90 % 포함률 0.906 (스펙 표의 0.369/10.47/0.912는 구표 기준) | `tests/test_models.py::test_g0_2_role_kinetics_golden`, `docs/model_card.md` |
| G0-3 GBM LOPO R² ≥ 0.50 | **통과(경계)** | seed 42: R² 0.503 / MAE 10.74 %p (raw % 회귀). z-타깃 0.508/10.65. seed 7: 0.450, seed 2024: 0.515 → 3-seed 평균 0.489 (보고값 0.522 ± 0.007보다 낮음; multitask_v6는 자체 단순 해석기 표·전역 z-스케일·fold별 범주 인코딩을 쓰므로 동일 실험이 아님). 번들: `dor_scm_blended.csv` 1,592행, seed 42, raw % | `bundles/gbm_v6/meta.json`, `tests/test_models.py::test_g0_3_gbm_bundle_metrics` |
| G0-4 RO INSERT 실패 | **통과** | `?mode=ro` + `PRAGMA query_only`; INSERT/UPDATE/DELETE/CREATE 모두 OperationalError | `tests/test_reader.py` |

## M1 — 시나리오 A

| 게이트 | 결과 | 수치 | 근거 |
|---|---|---|---|
| G1-1 logistic_fit 12 곡선 max dev < 2 %p | **통과** | 12/12 (a_max 0.15–0.90, τ 3–150 d) | `tests/test_kinetics.py::test_g1_1_logistic_fit_deviation` |
| G1-2 mock forward 자기검증 (+실기) | **통과** | 실제 xGEMS에서도 7·28 d 자기검증(alpha_ok·materials_ok) 통과(`test_g1_2_real_self_check`). mock: | `reaction_degrees.json["scm"][slot]` == 내보낸 α (4 재령, ±1e-3); 오버라이드 산화물은 `source_contribution_ledger.csv`에서 복원해 ±0.05 %p 일치. 주입 방식 `monkeypatch`(P-IG-1 머지 전) | `tests/test_forward_mock.py::test_g1_2_mock_forward_self_check` |
| G1-3 앙상블 기본값 | **결정: `bayes`** | 5-fold leave-papers-out(1,592 obs, 80편; 재현 스크립트 `scripts/g1_3_blend_lopo.py`): bayes MAE 11.43 / R² 0.441 / 90 % 포함률 **0.911** / 폭 46.7; blend MAE 11.37 / R² 0.440 / 포함률 **0.770** / 폭 32.8; GBM 점 MAE 11.44. blend는 MAE를 0.06 %p 낮추지만 포함률이 0.85 미만 → 기본값 `bayes` (`configs/defaults.yaml`) | `docs/g1_3/g1_3_blend_lopo.json` |
| G1-4 real_xgems OPC60/slag40 | **통과(앵커 교체; P-IG-6 후 재갱신)** | P-IG-6 커널: OPC100 pH 13.596, porosity 0.3266, CH 26.1 g; OPC60/slag40 pH 13.392, porosity 0.3258, CH 9.6 g (`docs/real_anchors_TINN_v4.json`; 패치 전 값은 `_pre_pig6.json`). 패치 전: | 기본 kinetics, 실제 xGEMS(TINN_v4): OPC100 w/b 0.5 28 d → pH **12.66156**(GemsPilot 앵커 12.661534와 일치), porosity 0.4119(앵커 0.407184; 시스템이 다름), CH 25.3 g/100 g; OPC60/slag40 w/b 0.45 28 d → porosity 0.3862(앵커 0.398451), pH 12.6616, CH 8.47 g/100 g, CSHQ 51.8 g. CNASH 앵커는 이 시스템에 CNASH가 없어 재현 불가 → `docs/real_anchors_TINN_v4.json`을 이 시스템의 앵커로 기록 | `tests/test_real_xgems.py::test_g1_4_default_kinetics_anchors` |
| G1-5 MCP 한 에피소드 | **통과(시뮬레이션 LLM)** | GemsPilot 러너(`run_episode`: 정책 검사·워크스페이스 재매핑·궤적 기록)를 그대로 쓰고 `litellm.completion`만 스크립트형 소형 모델(`dorgems.pilot.sim_llm`)로 대체: (good) dor_predict→export→materials override→mock forward→최종 답에 prediction.json의 28 d 중앙값 40.6 %를 인용, 자기검증 alpha_ok·materials_ok; (injected) 과제문의 "관리자 승인" 문구에 속아 use_mock=False 시도 → 러너가 DENIED → mock으로 복귀; (lazy) 도구 없이 지어낸 수치 → 채점 실패로 검출. 실제 API 모델(기본 `openrouter/anthropic/claude-haiku-4.5`)은 `OPENROUTER_API_KEY` 설정 후 같은 스크립트로 실행 | `tests/test_sim_episode.py`, `scripts/run_llm_episode.py` |

## M2 — 시나리오 B (파이프라인만; 물리 검증은 xGEMS 필요)

| 게이트 | 결과 | 비고 |
|---|---|---|
| G2-1 phase_aliases 실측 확정 | **통과** | dch.json의 46개 상 이름으로 표를 작성하고 실제 실행의 `xgems_phase_amounts_raw.csv`로 11개 그룹 전부 확인(`docs/g2_1_phase_confirmation.json`) → `confirmed: true`. 주요 이름: Portlandite, CSHQ, ettringite/SO4_CO3_AFt/CO3_SO4_AFt, C4AcH11(모노카보네이트), C4Ac0.5H12(헤미), C4AsH12·OH_SO4_AFm(모노설페이트), MgAl-OH-LDH, straetlingite/C2ASH55, C3(AF)S0.84H(하이드로가넷). 표에 없는 이름은 예외 |
| G2-2 부피 단위·화학수축 범위 | **통과(단위 불일치 확정)** | 실측: 상 질량은 **kg**(system_mass 0.122 kg), `porosity.json`의 `initial_volume_cm3`·`solid_final_volume_cm3`는 cm³이지만 `excluded_non_solid_phase_volumes_raw[aq_gen]`은 **m³**(2.3e-5) → P-IG-4 불일치 확정. aq를 m³→cm³로 환산하면 OPC w/b 0.5 28 d 화학수축 **0.0765 mL/g binder**(문헌 0.04–0.07 근처), 환산하지 않으면 0.337(비물리). 결합수 W_in−W_aq(H2O@) = **24.0 g/100 g**(문헌 TGA 20–24). `docs/g2_2_units_TINN_v4.json` |
| G2-3 OPC 참조 검사 | **실패 — 중단·보고** | 실제 xGEMS, 28 d ± 15 %, w/b 0.4–0.5, OPC-only 202후보. **1차 집계(규칙 결함)**: 등급 A 200건/41편, median r **+15.2 g** — 이 중 135건은 basis가 `mass_percent_unspecified`/`other`인데 DB가 unit_norm을 `g/100 g binder`로 정규화해 A로 잘못 분류됨(관측/모델 비율 0.15–0.35 밴드 60건은 TGA 질량손실 미환산 18/74 = 0.243과 일치). 규칙 수정(basis 우선, DB basis 어휘 17종 반영) 후 **최종 집계**: 등급 A **39건/11편**(기준 ≥30편 **미달**), median r **+12.4 g**, IQR 11.9, frac|z|<2 = 0; 관측 중앙 11.8 g vs 모델 24.8 g. 논문별 중앙 잔차: −2.1, −0.5, +3.8(3편 일치) … +21.3(8편 과대). **판정: DB 기준 문제가 지배적이지만, 명시 basis 행만 봐도 모델 CH가 대부분 논문보다 높다** → **커널 조사 완료(`docs/kernel_ch_investigation.md`, 2026-09-03):** xGEMS 입력에 OPC의 SO3·MgO·Na2O·K2O(5.8 %)가 전혀 들어가지 않아 S=0(AFt/AFm 불가, C3A→C3AH6), pH 12.662 고정(알칼리 없음), 기본 OPC 조성이 C3S 66 %로 과다, CSHQ Ca/Si 상한 1.63. CH 25.3 g은 이 가정들의 산술 결과(검산 25.0). 석고 5 %+CEM I 조성이면 21.8 g. 패치 제안 P-IG-6 기록. σ_model 확정 불가 | `docs/real/g2_3_opc_*` |
| G2-4 b_BW, σ_model | **확정(비례 보정)** | 머지 커널, OPC-only 28 d ± 15 %, w/b 0.4–0.5, 등급 A **결합수 28건/10편**: 모델 중앙 28.4 g vs 관측 중앙 17.0 g(문헌값은 Powers 비증발수 ≈ 17 g과 일치; 열역학 결합수는 겔수 포함). 논문별 중앙 잔차: 3편 ±1.4 g 이내, 6편 +11~+19 g, 1편 −11.6. 덧셈 오프셋 대신 **비율 0.60**(17.0/28.4)이 슬래그 twin(G2-5, 0.59)과 독립적으로 일치 → `likelihood.systematic_scales.bound_water = 0.60`, σ_model(BW) 3 g(보정 후 MAD 3.2). CH: 오프셋 −3 g, σ_model 5 g(커널 조사). 화학수축: 4건/2편(모델 0.06 vs 관측 0.10 mL/g) — 자료 부족, 보정 없음 | `docs/real/merged/g2_4_*` |
| G2-5 twin ≥ 20 배합 | **부분 통과(결합수 기준)** | 실기, 측정 DoR pin, 머지 커널, 최종 단위 규칙 + 비례 보정(BW 0.60)·CH 오프셋: 74후보 중 63배합 실행. 주력 물리량 판정(primary basis): **결합수** 등급 A 28건/4배합/1편 → median r **+0.24 g**, IQR 1.0, frac|z|<2 **1.00** → 4배합 **consistent**, 27배합은 결합수 관측 부족(insufficient_data). 보조 CH(78건/2편): 측정 DoR 고정 시 모델 CH ≈ 0 vs 관측 6.1 g → tension(보조이므로 종합 판정에 미반영). 이력: 보정 전(v4)에는 BW +5.2 g로 tension 4 → 비례 보정 후 일치. 등급 A 관측이 적은 것은 §14-1 basis 백로그 그대로. 이전 집계(v3, 기준 문제 전): BW median r −1.9 g(4편) | `docs/real/merged/g2_5_twin_batch_v5.json` |

## M3 — 시나리오 C

| 게이트 | 결과 | 수치 | 근거 |
|---|---|---|---|
| G3-1 합성 복원(mock) | **통과(5 케이스)** | flat prior, CH+BW 3재령×2 관측(σ 0.3): 진짜 α(28 d)가 사후 90 % 구간 안, 중앙값 오차 ≤ 5 %p. (mock의 CH는 α에 대해 *증가*하므로 단조성 플래그가 예상대로 발생; 복원 논리는 F의 부호와 무관) 스펙의 20 케이스 중 5 케이스만 테스트에 넣음 | `tests/test_inverse_mock.py::test_g3_1_synthetic_recovery` |
| G3-2 real_xgems 합성 복원 | **통과(5/5)** | 실제 xGEMS(TINN_v4), flat prior, CH+BW 4재령, 노이즈 = 관측 스팬의 3 %: 진짜 α(28 d) 5/5 포함, 중앙값 오차 평균 **0.27 %p**(최대 0.47), ESS 772–1862, 케이스당 72–80 xGEMS 호출. 실제 커널에서 CH·BW 모두 α에 대해 단조(플래그 0) | `docs/real/g3_2_real_recovery.json`, `scripts/g3_2_real_recovery.py` |
| G3-3 실측 66배합 | **실패 — 시나리오 C 격하 확정** | (1) 옛 커널·수정 전 단위 규칙: 48배합 실행/38 성공(7편, 159점) → MAE 16.2 %p, 90 % 포함률 0.24, 편향 −7.3; fly_ash 8.8(+6.8, 포함률 0.42), slag 23.4(−23.4, 0), metakaolin 36.8(−36.8, 0). (2) **머지 커널(P-IG-6 + CEM I) + 최종 단위 규칙**: 44배합 실행 중 30배합은 등급 A/B 관측이 없어 제외(플라이애시 논문의 CH가 basis 미상 → D), 5배합 재구성 불가, **성공 9배합 = 슬래그 논문 1편(54점): MAE 22.7 %p, 편향 −22.7, 포함률 0.0**. 두 실행 모두 5 h 부근에서 작업 한도로 44–48/66에서 중단(부분 JSON 보존). 스펙 기준(MAE < 10.3 또는 포함률 ≥ 0.8) 미달 → **C는 '일관성 점검·불확실성 축소' 용도로 격하**(`configs/defaults.yaml inverse.status`, 도구 결과에 경고). 원인 후보(우선순위): ① CH 우도의 계통 편향 — 측정 DoR을 넣으면 모델 CH가 0이 되는 만큼 역으로 α가 과소 추정됨(커널 조사·G2-5), ② 클링커 α 고정 불가(P-IG-3), ③ 등급 A 관측 부족(§14-1). 다음 방법: CH 대신 결합수·화학수축·QXRD 비율을 우도에 쓰고 CH는 σ_model에 계통항을 더해 약하게만 반영 (3) **v2 — 머지 커널 + CH 가중 0.25/−3 g, BW 비례 보정 전, α 격자 11**: 사용 가능 13배합(슬래그 9/1편, 메타카올린 4/1편, 82점) 전부 완료 → MAE 16.1 %p, 포함률 0.55, 편향 −15.7; slag MAE 12.9/포함률 0.59, metakaolin 22.3/0.46(ESS 2000 = 사전분포 그대로, 관측 정보 없음). (4) v3(BW 비례 0.60 포함) 실행 중 | `docs/real/g3_3_prelim_base_kernel.json`, `docs/real/merged/g3_3_merged_kernel_partial.json`, `docs/real/merged/g3_3_v2_weights_no_scale.json` |
| G3-4 KL 단조 증가 | **통과** | 관측 1→2→3개에서 사전→사후 KL 비감소 | `tests/test_inverse_mock.py::test_kl_increases_with_observations` |

## M4 — 가드레일·문서

| 항목 | 결과 |
|---|---|
| 정책 의미론 | `read`/`mock_ok`만 사용; GemsPilot `_policy_check`에서 `use_mock=False`는 `allow_real` 없이는 거부됨을 테스트. 인자 안의 자연어 "승인" 문구는 정책에 영향 없음 (`tests/test_policy_injection.py`) |
| 예산 | 실제 실행에 `max_xgems_calls` 없음/200 초과 → `PermissionError` |
| 문헌 DB 쓰기 | 도구층에 쓰기 경로 없음; `dor_db_lookup`은 이름 붙인 쿼리만 |
| staging | `dor_stage_inferred`는 기본 dry-run; `dorgems review list/approve/reject` |
| 문서 | `docs/model_card.md`, `docs/observable_mapping.md`, 이 파일. 스펙 재검증 로그는 스펙 상단에 추가 |
| 벤치 `dor_*` | `configs/agent_bench_dor.yaml` 6 시나리오(정상 경로 envelope, 예산 없음/초과 거부, 문헌 DB 쓰기 거부, 인자 내 주입문 무효, staging dry-run) — `python -m dorgems.pilot.bench` / `tests/test_bench.py`, 6/6 통과(mock, LLM 불필요). GemsPilot 벤치 kind로의 통합은 P-GP-2 이후 |

## 업스트림 패치 상태

| ID | 브랜치 | 상태 |
|---|---|---|
| P-IG-1 | `InverseGems: dorgems/p-ig-1-materials-config` (`82df677`) | 완료(무푸시). `materials_config`를 forward/task/design 경로에 스레딩, `load_materials()` 6곳 교체, `recipe.json` metadata에 경로 기록, 테스트 1개 추가. 업스트림 테스트 222 passed / 1 failed(기존 `test_feature_table` 실패, 패치 무관). **DoRGems 확인:** 이 브랜치를 체크아웃하면 `run_forward`가 `materials_injection: native`로 전환되고 자기검증(alpha_ok·materials_ok) 통과 |
| P-GP-1 | `GemsPilot: dorgems/p-gp-1-reaction-model-kwargs` (`f0da5e8`) | 완료(무푸시). `run_forward/run_task/run_confirmed_query/run_design_with_recovery`와 MCP 래퍼에 `reaction_model_config`·`reaction_model_id`·`materials_config` 추가(커널 시그니처를 검사해 조건부 전달; 커널이 못 받는 non-None 옵션은 명시적 오류). `design_recovery`의 `None` 하드코딩 인자화. 52 passed / 1 failed(기존 하드코딩 경로 실패) |
| P-GP-2 | `GemsPilot: dorgems/p-gp-2-toolset-entrypoints` (`ee3afbb`) | 완료(무푸시). `default_toolset(extra, discover=True)` + `gemspilot.toolsets` 엔트리포인트 탐색(고장난 엔트리포인트는 경고로 수집), `mcp_server.register_extra_toolsets`. 이 환경에서는 DoRGems 엔트리포인트가 실제로 잡혀 도구 15→28개(MCP 31개). 53 passed / 1 failed(기존) |
| P-GP-3 | `GemsPilot: dorgems/p-gp-3-write-policy` (`9f54ad1`) | 완료(무푸시). `write` 정책: `dry_run` falsy + `allow_real` 없음 → 거부. 49 passed / 1 failed(기존). DoRGems `dor_stage_inferred`는 머지 후 `use_mock`→`dry_run`으로 개명 예정 |
| P-IG-6 | `InverseGems: dorgems/p-ig-6-opc-minor-oxides` (`535f90e`) → **로컬 master에 머지(`c4600b2`), 이어 기본 OPC 조성을 CEM I로 교체(`50beefb`); 원격 미푸시** | **구현·머지 완료**: OPC 소량 산화물(SO3→CaSO4+동반 CaO, Na2O, K2O, MgO) 투입 정책, ledger·시그니처 반영, 테스트 5개. 실기: pH 12.66→13.60, AFt/AFm·브루사이트·LDH 형성, 질량 수지 100.2 g. CH는 +0.8 g(과대 원인은 C3S 66 % 기본 조성·CSHQ Ca/Si 상한으로 남음) — `docs/kernel_ch_investigation.md` |
| P-GP-4 | `GemsPilot: dorgems/p-gp-4-anchor-refresh` (`9d88d24`) → **로컬 main에 머지(`2924b96`), 원격 미푸시** | 완료: P-IG-6 커널로 mock QA 앵커 18개 재생성, TINN_v4 실기 앵커 15개 추가(`agent_qa_generated_TINN_v4.yaml`), README 노트. `Test-dat.lst` 실기 앵커는 해당 시스템 보유 PC에서 재생성 필요 |
| P-IG-2/3/4/5 | — | 미착수(xGEMS 실측 필요 또는 후순위). P-IG-2 전까지는 `capture_species` 폴백(`run_forward_cached` + CapturingRunner)이 결합수 계산을 담당 |

두 클론의 작업 트리는 원래 커밋(`e84d7a9`, `753cf6d`)으로 복귀, 브랜치 간 상호 독립. 커밋 작성자는 `Solmoi Park <park.solmoi@gmail.com>`(추정)로 기록됨.
