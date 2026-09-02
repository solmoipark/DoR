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
| xGEMS | **없음** — `py313-xgems` 환경과 GEMS3K 시스템 파일이 이 PC에 없음. `real_xgems` 마커 테스트는 전부 skip. **M2의 물리 게이트(G2-1..G2-5), G1-4, G3-2는 실행 불가** |
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
| G1-2 mock forward 자기검증 | **통과** | `reaction_degrees.json["scm"][slot]` == 내보낸 α (4 재령, ±1e-3); 오버라이드 산화물은 `source_contribution_ledger.csv`에서 복원해 ±0.05 %p 일치. 주입 방식 `monkeypatch`(P-IG-1 머지 전) | `tests/test_forward_mock.py::test_g1_2_mock_forward_self_check` |
| G1-3 앙상블 기본값 | **결정: `bayes`** | 5-fold leave-papers-out(1,592 obs, 80편; 재현 스크립트 `scripts/g1_3_blend_lopo.py`): bayes MAE 11.43 / R² 0.441 / 90 % 포함률 **0.911** / 폭 46.7; blend MAE 11.37 / R² 0.440 / 포함률 **0.770** / 폭 32.8; GBM 점 MAE 11.44. blend는 MAE를 0.06 %p 낮추지만 포함률이 0.85 미만 → 기본값 `bayes` (`configs/defaults.yaml`) | `docs/g1_3/g1_3_blend_lopo.json` |
| G1-4 real_xgems OPC60/slag40 | **미실행** | xGEMS 없음. GemsPilot 회귀 앵커(porosity 0.398451, CNASH 0.045535) 미확인 | — |
| G1-5 MCP 한 에피소드 | **부분** | `dorgems-mcp` 서버 구현 + `dor_run_envelope`가 predict→export→override→forward×3를 한 호출로 완료(mock, CLI로 실증). LLM 호스트가 실제로 연결된 에피소드는 미실행(API 키 없음) | `tests/test_forward_mock.py`, CLI `dorgems envelope` |

## M2 — 시나리오 B (파이프라인만; 물리 검증은 xGEMS 필요)

| 게이트 | 결과 | 비고 |
|---|---|---|
| G2-1 phase_aliases 실측 확정 | **미실행** | `configs/phase_aliases.yaml`은 `confirmed: false`; `phases.confirm_from_raw_names()`로 실제 `xgems_phase_amounts_raw.csv`를 읽어 확정하는 절차만 준비 |
| G2-2 부피 단위·화학수축 범위 | **미실행** | `observables.chem_shrink_ml_per_g(volume_unit=…)` 인자로 cm³/m³ 전환 가능; 실측 후 확정 |
| G2-3 OPC 참조 검사 | **미실행(파이프라인 통과)** | 후보: 28 d ± 15 %, w/b 0.4–0.5, SCM 포함 binder JSON 제외 → 241행/229배합/51편, 등급 A 202행(190배합/41편; 스펙 상한 393은 w/b·재령 창 적용 전 값). `dorgems opc-check --max-mixes 8` mock 실행 → comparison.csv/json/summary.md 생성. T 결측 배합은 20 °C 가정+플래그(§8.2). mock 수치는 물리적 의미 없음 |
| G2-4 b_BW, σ_model | **미확정** | `defaults.yaml`의 초기값(σ_model CH 2.5, BW 3.0, CS 0.01) 유지 |
| G2-5 twin ≥ 20 배합 | **미실행(파이프라인 통과)** | `dorgems compare`(twin batch, mock): 후보 74배합(DoR ≥ 3 재령 & CH/BW/CS 공통 재령 ≥ 1) 중 63배합 실행(11배합은 w/b 또는 T 결측으로 제외), 전부 측정 DoR pin. 판정 분포(mock, 의미 없음): insufficient_data 49 / tension 14 — 대부분 관측이 등급 C·D(paste 기준·basis 미상)라 통계에서 제외됨 → §14-1의 basis 재판정 백로그가 실제 병목 |

## M3 — 시나리오 C

| 게이트 | 결과 | 수치 | 근거 |
|---|---|---|---|
| G3-1 합성 복원(mock) | **통과(5 케이스)** | flat prior, CH+BW 3재령×2 관측(σ 0.3): 진짜 α(28 d)가 사후 90 % 구간 안, 중앙값 오차 ≤ 5 %p. (mock의 CH는 α에 대해 *증가*하므로 단조성 플래그가 예상대로 발생; 복원 논리는 F의 부호와 무관) 스펙의 20 케이스 중 5 케이스만 테스트에 넣음 | `tests/test_inverse_mock.py::test_g3_1_synthetic_recovery` |
| G3-2 real_xgems 합성 복원 | **미실행** | — |
| G3-3 실측 66배합 | **미실행** | 실행 경로(`dor_infer_from_observations --mix-uid`)는 구현·mock 통과; 판정은 실제 xGEMS 필요 |
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
| P-IG-2/3/4/5 | — | 미착수(xGEMS 실측 필요 또는 후순위). P-IG-2 전까지는 `capture_species` 폴백(`run_forward_cached` + CapturingRunner)이 결합수 계산을 담당 |

두 클론의 작업 트리는 원래 커밋(`e84d7a9`, `753cf6d`)으로 복귀, 브랜치 간 상호 독립. 커밋 작성자는 `Solmoi Park <park.solmoi@gmail.com>`(추정)로 기록됨.
