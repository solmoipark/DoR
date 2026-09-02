# DoRGems — SCM 반응도(DoR) 에이전트 설계 스펙 v1

작성 2026-09-02 · 코딩 세션(Codex / Claude Code)용 작업 명세
대상 자산: `DoR of SCMs in blended cements/` (문헌 DB + 모델링), `solmoipark/InverseGems` (커널), `solmoipark/GemsPilot` (에이전트 레이어)

> **이름은 임시.** 이 문서에서 새 패키지를 `DoRGems`(파이썬 패키지 `dorgems`)로 부른다. 확정 시 전체 치환.
>
> **재검증 로그** — 2026-09-02: §1·§6·§15의 코드 인용을 InverseGems `e84d7a9`, GemsPilot `753cf6d` 클론과 로컬 파일에 대해 독립 대조하여 수정 반영(정책 의미론, `runner_factory` 범위, `config_path` 패치 지점, `signature_files` 키, `system_type` 유래, OPC 참조군 상한, HDI 수치). 구현 세션은 이 줄 아래에 자기 재검증 결과를 추가한다.
>
> **구현 세션 재검증** — 2026-09-02 (Claude Code, InverseGems `e84d7a9`, GemsPilot `753cf6d` 신규 클론, 문헌 DB 미변경):
> 1. §1.2·§1.3의 파일:라인 인용(`scm_reaction.py` 레지스트리·clip, `materials.py` SCM_NAMES, `utils.config_path`, `reaction_parameters.py` raw_config 보존, `reaction_model.py` signature_files, `chem_hash.py`, `cached_forward.py` runner_factory, 무인자 `load_materials()` 6곳, GemsPilot `_policy_check` 219-229·`ToolSpec`·15개 도구·`run_forward/run_task/run_confirmed_query` 인자 갭)은 **모두 일치**.
> 2. **정정(§1.2, §6.3, §12 G1-2):** `forward_query`/`run_forward_cached` 경로는 `input_reaction_degrees.json`·`input_materials_used.json`을 만들지 **않는다**. 이 두 파일은 `xgems_output_capture.save_run_outputs`(독립 hydration 스크립트 경로)에서만 나온다. 캐시 경로의 산출물은 `<db>/recipe_runs/<recipe_id>/{reaction_degrees.json(키 age_days·opc·scm·reaction_parameter_set), recipe.json, unreacted_masses.json, porosity.json, volumes.json, source_contribution_ledger.csv}`와 `<db>/chemistry_runs/<chem_hash>/xgems_raw/*`. 또한 `db` 인자는 sqlite 파일이 아니라 **디렉터리**(`InverseGemsDatabase(db_dir)`). DoRGems의 자기검증은 `reaction_degrees.json["scm"]`과 ledger의 산화물 몰(→ 조성 복원)로 수행한다.
> 3. **정정(§1.1, §16-2):** `modeling/bayes_role_kinetics.csv`·`bayes_oof.csv`(8/26)는 1,177행 구표 기준이고 현행 `dor_scm_final.csv`(8/28)는 1,610행이라, 재실행 결과가 ±0.5 %p 안에 들 수 없다(fly_ash 28.1→30.6, other 47.3→39.0, slag 58.5→56.9 %). 재실행 사후를 번들 골든으로 채택하고 `docs/gates.md` G0-2에 기록. `work/bayes_idata.nc`는 예상대로 없었고 재실행으로 생성.
> 4. **정정(§1.1 표):** GBM DoR-only LOPO는 seed 42에서 R² 0.503(스펙대로 `dor_scm_blended.csv`·raw % 회귀); 3-seed 평균 0.489. 보고값 0.522는 `multitask_v6.py`의 자체 표·전역 z-스케일 실험이라 동일 조건이 아님.
> 5. **환경:** 이 PC에는 `py313-xgems` 환경·GEMS3K 파일이 없어 real_xgems 게이트(G1-4, G2-*, G3-2/3)는 실행 불가. mock 파이프라인은 전 구간 통과.
> 6. **§5.4 결정(G1-3):** blend는 포함률 0.77 < 0.85 → 기본값 `bayes`.
> 7. **xGEMS 확보 후 추가 재검증(2026-09-02 저녁, TINN_v4 시스템, xgems 2.1.2):** (a) §1.2 "`porosity.yaml`은 cm3, `backfill.py`/`cli.py` 기본은 m3" 의혹은 실측으로 확정 — `porosity.json`의 `initial_volume_cm3`·`solid_final_volume_cm3`는 cm³, `excluded_non_solid_phase_volumes_raw`는 raw xGEMS **m³**(P-IG-4 대상). 상 질량은 kg. (b) §12 G1-4의 GemsPilot 앵커는 `Test-dat.lst`(CNASH 시스템) 기준이라 이 시스템에서는 pH만 재현(12.66156 vs 12.661534); porosity 0.4119/0.3862, CNASH 없음(`CSHQ`). (c) §4.3 표의 "`mass_percent_unspecified` → D"는 DB가 unit_norm을 `g/100 g binder`로 정규화해 두어 unit 기준 매칭으로는 A로 새는 함정이 있다 — basis 우선 매칭이 필수. DB의 basis 어휘는 17종(`per g binder/cement/clinker/paste/sample/MK`, `per_100g_ignited` 등)이며 chem_shrink는 unit_norm이 `%`로 무의미해 unit_reported로 판정해야 한다. (d) G2-3는 최종 규칙에서도 실패(모델 CH가 문헌보다 +12 g) — 스펙대로 중단·보고.

---

## 0. 이 문서를 쓰는 방법

이 문서는 구현 세션에 그대로 넘기는 작업 명세다. 구현 에이전트는 다음 규칙을 따른다.

1. **§1의 "검증된 사실"은 2026-09-02 시점 코드에서 직접 확인한 것이다.** 파일:라인이 적힌 항목은 구현 전에 다시 열어 확인하고, 달라졌으면 문서를 고친 뒤 진행한다. 추측으로 채우지 않는다.
2. **마일스톤은 순서대로, 결정 게이트(§12)를 통과해야 다음으로 간다.** 게이트는 pytest로 자동화되는 것을 원칙으로 한다.
3. **문헌 DB는 절대 쓰지 않는다.** `scm_dor_enriched.db`는 읽기 전용(`?mode=ro`)으로만 연다. 에이전트가 만들어내는 값(추론 DoR 등)은 별도 staging DB에만 쓴다(§11).
4. **LLM은 숫자를 만들지 않는다.** GemsPilot의 원칙("The language model never produces a scientific conclusion")을 그대로 계승한다. 모든 예측·비교·판정은 `dorgems` 커널 함수가 결정론적으로 계산하고, LLM은 도구를 선택·파라미터화만 한다.
5. **하지 말 것(§14 참조)**: 비정질 함량을 100−Σ결정상으로 계산, 존재하지 않는 Cemdata 상 이름 발명, GemsPilot 가드레일 기본값(`use_mock=True`) 완화, GEMS3K 시스템 파일 커밋, 문헌 표에 대한 어떤 UPDATE/INSERT.

---

## 1. 배경: 세 자산의 현재 상태 (검증된 사실)

### 1.1 문헌 DB와 DoR 모델링 (`DoR of SCMs in blended cements/`)

**현행 DB는 `modeling/scm_dor_enriched.db`** (루트의 `scm_dor.db`는 확장 전 구버전, 스키마 동일). SQLite 4개 표:

| 표 | 행 | 핵심 컬럼 |
|---|---|---|
| `papers` | 1,350 | `doi` PK, title, year, journal, oa_status, retrieval_route, extraction_notes |
| `materials` | 1,634 | `material_uid`(= doi::material_id) PK, paper_doi, material_id, **role**, name_in_paper, CaO SiO2 Al2O3 Fe2O3 MgO SO3 Na2O K2O TiO2 LOI, blaine_m2_kg, d50_um, bet_m2_g, amorphous_pct, source_locator |
| `mixes` | 6,870 | `mix_uid`(= doi::mix_id) PK, paper_doi, **binder_composition_json**(material_id→질량%), scm_total_pct, primary_scm_role, primary_scm_pct, w_b, sand_binder, curing_temp_C, **curing_type**, notes |
| `observations` | 44,011 | `obs_uid` PK, paper_doi, mix_uid, **age_d, quantity, phase_name, value_norm, unit_norm**, basis_reported, method, method_detail, uncertainty, fig_only, extraction_confidence, sanity_ok, reviewed |

`quantity` 값과 유효(`value_norm IS NOT NULL`) 건수: `cum_heat` 17,637 · `QXRD_phase` 5,956 · `CH_TGA` 4,203 · `bound_water` 3,468 · `DoR_SCM` 1,700 (`%` 1,676 + `fraction` 24) · `DoR_clinker` 1,407 · `chem_shrink` 532 · `CH_XRD` 391. `mixes.curing_type = 'model_system'`은 합성 모델계이며 blended-cement 모집단에서 제외한다(`multitask_v6.py:45-60`).

**단위 기준(basis)은 이질적이다.** CH_TGA는 `g/100 g binder` 3,180건이 주류지만 `g/100 g paste`(219), `g/100 g clinker`(201), `g/100 g cement`(175), `g/100 g sample`(112), `g/100 g PC/OPC`, `g/100 g MK`, `%`(44) 등이 섞여 있고, `basis_reported`는 `mass_percent_unspecified`가 유효 행 기준 1,790건(전체 행 3,074건)이다. bound_water도 같은 구조(유효 행 1,071건). QXRD_phase는 전부 `wt% (as reported)`. chem_shrink는 `%` 488 / `mL/g binder` 44. 시나리오 B·C의 단위 조화(§4.3)가 필수인 이유다.

**모델링 표 생성 로직은 `modeling/build_clean.py`에 있다** (스크립트 형태, import 불가). SCM 역할 해석기(`role_from_text`, `CEMENT_ROLES`, `FILLER_ROLES`, phase_match → json → unique_mat → name_match → common_role → keyword_role 순의 fallback), `method_group()` 매핑, 단위 `fraction`→×100, 0–100 범위 필터, 중복 제거 키 `(mix_uid, age_d, method_group, dor_pct)`. 출력 `dor_scm_final.csv` 1,610행(28컬럼: obs_uid, paper_doi, mix_uid, age_d, dor_pct, method_group, confidence, fig_only, scm_role, scm_pct, resolve_how, w_b, curing_temp_C, system_type, scm_CaO…scm_amorphous_pct). blended 기준선 `dor_scm_blended.csv`는 1,592행 / 80편 / 476배합. 결측률: scm_bet 91%, d50 90%, amorphous 76%, blaine 63%, 화학조성 16–20%, w_b 19%, T 20%.

**모델 성능(leave-one-paper-out, 7차 보고서 기준)**

| 모델 | 스크립트 | R² | MAE (%p) | 비고 |
|---|---|---|---|---|
| GBM tabular (LightGBM, TIGHT) | `multitask_v6.py` DoR-only | 0.522 ± 0.007 (blended 0.531) | 10.3–10.4 | 점예측 최고 |
| 계층 베이지안 v4 (PyMC) | `bayes_hier_v4.py` | 0.369 | 10.47 | **90% 예측구간 포함률 91.2%**, 평균 폭 44 %p |
| Kinetics v2 (앵커 재매개변수화) | `train_kinetics_v2.py` | 0.380 | 11.8 | 앵커 예측 자체가 약함 |
| 멀티태스크(보조 타깃을 task로 풀링) | `multitask_v6.py` | 0.437 | 11.66 | **유해함으로 확정** — 쓰지 않는다 |

베이지안 v4 구조(`bayes_hier_v4.py:32-38, 72-87`):

```
alpha(t) = a_max · (1 − exp(−(t/τ)^0.5))                      # β 고정 0.5
logit(a_max/100) = a0[role] + X·β_a + u_paper,  u_paper ~ N(0, sd_paper)
log τ            = t0[role] + X·β_t                            # τ에는 논문효과 없음
y ~ Normal(alpha(t), σ[method_group])
X = 표준화된 [scm_pct, w_b, curing_temp_C, CaO_SiO2]  (결측은 중앙값 대치 후 표준화)
role ∈ {slag, fly_ash, other}
```

추정치(평균 배합 조건): slag a_max 58.5% [50.8, 66.3], τ 18.5 d · fly_ash 28.1% [21.7, 35.5], τ 27.8 d · other 47.3% [34.2, 72.8], τ 13.1 d (`bayes_role_kinetics.csv` 기준; 3차 보고서의 구간은 이전 실행값이라 소폭 다름). 측정법별 관측 노이즈 σ: selective_dissolution 5.5, XRD 5.6, mass_balance 8.0, SEM_BSE 13.0, NMR 22.1 %p (위치 편차는 논문과 교락되어 식별 불가 — 3차 보고서 §3). 사후표본은 스크립트상 `work/bayes_idata.nc`에 저장되지만 **현재 폴더 목록(`modeling/work/`)에 이 파일이 없어 재실행이 필요할 가능성이 크고**, **특징 표준화 상수(평균·표준편차·중앙값)는 어디에도 저장되지 않는다** → 번들 내보내기 때 재계산 필요(§5.1).

GBM v6 특징(`multitask_v6.py:22-24, 94-99`): `log_age`(log10), `scm_pct, w_b, curing_temp_C`, 조성 파생 `CaO_SiO2, basicity=(C+M+A)/S, pozz_sum=S+A+F, Al_Si, amorph, fineness(blaine)`, 범주형 `scm_role, method_group`(결측 → `'missing'`). 하이퍼파라미터 `TIGHT = dict(objective='regression', num_leaves=7, max_depth=4, learning_rate=0.03, n_estimators=500, min_child_samples=40, subsample=0.7, subsample_freq=1, colsample_bytree=0.6, reg_lambda=20.0, reg_alpha=2.0, verbose=-1)`.

**검증 설계에 쓸 중첩 규모(model_system 제외, 같은 배합·같은 재령):** DoR_SCM과 함께 관측된 CH_TGA 79배합/19편, bound_water 75배합/8편, QXRD_phase 39배합/10편, DoR_clinker 109배합/16편. DoR_SCM + (CH 또는 결합수)가 **3개 이상 공통 재령**에 있는 배합 66개. OPC 단독 배합 중 28일 CH_TGA와 w/b가 있는 것 393배합/84편(`scm_total_pct`가 NULL 또는 0; 이 중 29배합은 binder JSON에 SCM 재료가 있어 제외하면 364배합/79편 — 393은 상한).

### 1.2 InverseGems (커널, `src/inverse_gems/`, ~26k LOC, BSD-3)

- **"inverse"는 inverse *design*이다.** 목표 상조성(surrogate 예측) → 배합 탐색. 관측치 → 반응도를 푸는 inverse *analysis*는 없다. 유일한 실험데이터 입력은 `kinetics_calibration.calibrate_scm_kinetics(data_csv[scm, age_d, dor], ...)` 곡선 피팅이다.
- **SCM 반응도는 쿼리 파라미터가 아니다.** `reaction_model_config` YAML → `load_reaction_parameters()` → `scm_alpha(age, params)`로만 들어간다. 기본 모델 `five_param_logistic`: `alpha = D + (A − D)/(1 + (t/C)^B)^G` (`scm_reaction.py:70-74`), 결과는 **분율**로 `np.clip(alpha, 0, 1)` (`scm_reaction.py:191-193`) — 퍼센트를 넣으면 1.0으로 잘린다. 기본값(`configs/scm_reaction.yaml`): slag `{A:0,B:0.75,C:20,D:0.55,G:1}`, fly_ash `{A:0,B:1.05,C:35,D:0.40,G:1}`, metakaolin `{A:0,B:0.95,C:5,D:0.55,G:1}`, silica_fume `{A:0,B:0.80,C:3,D:0.85,G:1}`.
- **상수 α 고정(pin)은 가능하다(검증됨):** `{A: α, B: 1, C: 1, D: α, G: 1}` + `availability_modifier: {enabled: false}` → 모든 재령에서 alpha = α.
- **커스텀 kinetics 등록:** `register_scm_kinetics(name, *, required, asymptote_key)` 데코레이터, 모듈 레벨 딕셔너리 `_KINETICS_REGISTRY` (`scm_reaction.py:40-63`). YAML의 `model:` 키로 선택(`scm_reaction.py:152-163`). **엔트리포인트 탐색이 없어 같은 프로세스에서 `load_reaction_parameters` 전에 등록해야 한다.** 가용성 보정기(`availability_modifier.py:94-99`)는 `asymptote_key`가 가리키는 값을 `.D`/`.with_D`로 다시 쓴다(분율, [0,1] 클램프).
- **산화물 조성은 `configs/materials.yaml`에 고정.** 경로 해석은 패키지 상대(`project_root()/configs`) → CWD 순(`utils.py:16-28`). 상위 API 어디에도 materials 경로 인자가 없고 `load_materials()`가 무인자로 6곳에서 호출된다(`cached_forward.py:229`, `chemistry_candidate_table.py:81`, `forward_query.py:88`, `xgems_preflight.py:252`, `backfill.py:335`, `cli.py:113`). 하위 함수(`parse_recipe`, `build_xgems_input`, `compute_porosity`, `build_source_ledger`, `C3S_C2SAvailabilityModifier`)만 `materials=`를 받는다.
- **SCM 이름은 고정 집합.** `SCM_NAMES = {slag, fly_ash, metakaolin, silica_fume}`, `BINDER_COMPONENTS = SCM_NAMES ∪ {OPC, limestone, gypsum}` (`materials.py:10-11`). `materials.yaml`에 없는 이름은 `materials.py:92`(`canonicalize_material_name`)에서, 있어도 `BINDER_COMPONENTS` 밖이면 `recipe.py:88-89`·`forward_query.py:92-93`에서 거부. `materials.yaml`의 `aliases:`로 기존 슬롯에 새 이름만 붙일 수 있다. 5번째 SCM을 추가하려면 두 집합 + 하드코딩 리스트 **14곳 이상**(`sampling.py:14, 471-479`, `chemistry_design_query_runner.py:127`, `xgems_quality_cases.py:15`, `task_query_preview.py:53`, `feature_diagnostics.py:13`, `design_run_report.py:13`, `design_query.py:18-26`, `inverse_forward_workflow.py:24-32`, `inverse_query.py:116-122`, `candidate_review.py:294-300, 348-354`, `candidate_selection.py:540-546`, `cli.py:1258-1266`, `forward_query.py:57-63, 166-174`, `configs/c3s_c2s_availability.yaml`)을 고쳐야 한다.
- **캐시 키 `chem_hash`**(`chem_hash.py:84-94`)는 `chem_hash_version`, 원소 몰 벡터(Ca Si Al Fe Mg S Na K C H O, 12유효숫자), water_mol, T, P, dat_lst sha256, species_map hash, run_mode, `thermodynamic_database_identifier`(항상 None)로 구성. 반응모델 시그니처·재령·배합은 포함되지 않고 `recipe_runs`에 별도 컬럼으로 남는다. → **산화물·DoR을 바꿔도 캐시는 안전**(원소 벡터가 바뀌면 다른 키, 같으면 물리적으로 같은 입력).
- **클링커 수화도는 Parrot-Killoh로만 결정**되고 공개 경로로 고정할 수 없다(`build_xgems_input`에 precomputed degrees 인자 없음, `xgems_input_builder.py:107-123, 166`).
- **원시 출력 디렉터리에 `phase_species_moles`가 저장되지 않는다**(`database.py:452-457`는 phase_masses, phase_volumes, aqueous_species, scalars, attribute_report만 기록). 수화물 내 결합수는 인프로세스에서 `runner.capture_raw_state()`를 가로채야 계산할 수 있다. 수용액 상 이름은 `aq_gen`. `porosity.json`에는 `initial_volume_cm3, solid_final_volume_cm3, xgems_solid_phase_volume_cm3, unreacted_binder_volumes_cm3, excluded_non_solid_phase_volumes_raw` 등이 있어 화학수축을 유도할 수 있다. **단위 주의:** `porosity.yaml`은 `cm3`, `backfill.py:329`/`cli.py:921` 기본은 `m3` — 실제 xGEMS 부피 단위를 M2에서 반드시 실측 확인.
- 시스템 기준은 **결합재 100 g + 물**, xGEMS 입력은 kg. `reacted_only` 모드에서는 반응한 분율만 계에 넣는다.
- GEMS3K 시스템 파일(`Test-dch/ipm/fun/dbr.json`)은 레포에 없다(.gitignore). xGEMS는 별도 conda 환경(`py313-xgems`). Mock runner(`MockXGEMSRunner`)로 파이프라인 전체를 xGEMS 없이 돌릴 수 있다.

### 1.3 GemsPilot (에이전트 레이어, `src/gemspilot/`, ~3.5k LOC, BSD-3)

- 이미 MCP 서버다: `gemspilot-mcp`(stdio), `mcp_server.py`에 18개 `@app.tool()`. litellm 기반 러너(`runner.py`), `ToolSpec(name, func, policy)` 정책 `read | mock_ok | real_gated`, 워크스페이스 격리, `INVERSE_GEMS_ARTIFACT_ROOTS` 아티팩트 허용목록, 세션 메모리(JSONL), 관찰-재계획 루프, GEMS-Agent-Bench.
- 도구 반환 계약: `{"contract": "inverse-gems-tool/1.0", "tool", "ok", "summary": {status, task_type, run_dir, answer_available, missing_outputs, result_summary}, "artifacts": {name: path}, "warnings", "error"}`.
- **갭:** `agent_tools.run_forward / run_task / run_confirmed_query`와 MCP 래퍼(`mcp_server.py:72-141`)는 커널이 받는 `reaction_model_config`, `reaction_model_id` 인자를 받지도 전달하지도 않는다(`agent_tools.py:258-268, 289-299, 320-329`; `design_recovery.py:54-56, 347-358`은 `None` 하드코딩). 예외적으로 **design** 작업은 task_query YAML 안의 `design_query.reaction_model.config`를 커널이 직접 읽으므로(`design_query.py:353-366`) `run_task`로 우회 전달이 가능하다 — 갭은 **forward 실행과 diagnose 경로**에서 완전하다. 유일한 DoR 경로는 `calibrate_scm_kinetics` → `reaction_parameters.<id>.yaml` 아티팩트인데 그것을 forward에 넘길 수단이 없다.
- **정책 검사의 실제 의미**(`runner.py:219-229 _policy_check`): `read`는 항상 허용; 그 외는 인자 `use_mock`(기본 True)만 본다. `use_mock=False`일 때 `mock_ok`는 `Episode.allow_real`이면 허용, **`real_gated`는 무조건 거부**. `use_mock` 인자가 없는 도구는 정책과 무관하게 항상 통과한다. 현재 `default_toolset()` 15개 도구 중 `real_gated`는 없다. → DoRGems 도구 정책 설계는 이 의미론을 따른다(§10, §11).
- 환경변수 `INVERSE_GEMS_ROOT`(커널 트리), `INVERSE_GEMS_ARTIFACT_ROOTS`, `GEMSPILOT_ROOT`; LLM은 `OPENROUTER_API_KEY`/`OPENAI_API_KEY`. 테스트에 하드코딩 경로 `C:\Users\solmo\InverseGems v2` 존재.

### 1.4 결론 — 무엇이 이어지고 무엇이 비어 있나

이어지는 것: 베이지안 v4의 `(a_max, τ)` 곡선은 InverseGems `reaction_model_config`로 **직접 변환**된다(§6). DB의 CH·결합수·QXRD·화학수축은 xGEMS forward 출력과 **같은 물리량**이라 대조 가능하다(§8).

비어 있는 것: (i) 새 SCM 조성 주입 경로, (ii) GemsPilot의 `reaction_model_config` 전달, (iii) 관측치→DoR 역해석 루프, (iv) 결합수 계산용 원시 출력, (v) 클링커 DoH 고정. (i)(ii)(iv)는 소규모 업스트림 패치(§6.4), (iii)은 DoRGems의 본체(§9), (v)는 후순위.

---

## 2. 설계 원칙

1. **커널이 계산하고 파일럿이 조율한다.** `dorgems`(커널)는 LLM·네트워크 의존이 없는 순수 파이썬 패키지다. `dorgems.pilot`는 GemsPilot `ToolSpec` 규약을 따르는 얇은 래퍼만 담는다. 커널은 CLI로도 완전히 동작해야 한다.
2. **결정론과 출처(provenance).** 같은 입력·같은 번들·같은 시드는 같은 출력을 낸다. 모든 산출물은 `manifest.json`에 DoRGems 버전, 모델 번들 해시, 학습 DB sha256, InverseGems/GemsPilot 커밋, 입력 해시를 기록한다.
3. **불확실성은 1급 시민.** 모든 DoR 예측은 점추정과 함께 5/50/95 분위(잠재 곡선 기준과 관측 기준 두 가지)를 낸다. InverseGems로 내보낼 때도 중앙값 하나가 아니라 `{lo90, median, hi90}` 세 config를 만든다.
4. **문헌 표는 불변, 추론값은 격리.** 문헌 DB는 읽기 전용. 에이전트가 만든 값은 `dorgems_staging.sqlite`의 `inferred_dor` 표에만 쓰고 `reviewed=0`으로 시작한다. 검토 승인 없이는 학습 데이터로 재유입되지 않는다.
5. **GemsPilot 가드레일 계승.** 예측·조회는 `read`; xGEMS를 부르는 도구는 `mock_ok`(기본 `use_mock=True`, 실제 실행은 호스트의 `allow_real`에서만); staging 쓰기 도구도 같은 게이트를 타도록 `use_mock`(dry-run 의미) 인자를 둔다(§1.3 정책 의미론, P-GP-3). `real_gated`는 GemsPilot에서 "실제 실행 무조건 거부"를 뜻하므로 쓰지 않는다. `XGEMSCallBudget` 예산은 모든 xGEMS 호출 경로에 부착한다.
6. **분포 밖(OOD)은 숨기지 않는다.** 학습 분포 밖의 SCM(예: 소성점토·메타카올린·실리카퓸은 각 1–8편)에 대해서는 넓은 구간과 명시적 `ood_flags`를 내고, 가장 가까운 문헌 유사사례를 근거로 함께 제시한다.
7. **업스트림 최소 침습.** InverseGems·GemsPilot 수정은 인자 전달·훅 추가 수준으로 제한하고 각 패치를 별도 PR로 분리한다(§6.4). DoRGems는 패치 없이도 mock 파이프라인이 돌아가야 한다(폴백 경로 유지).

---

## 3. 아키텍처

### 3.1 구성도

```
                ┌──────────────────────────── LLM 호스트 ────────────────────────────┐
                │  Claude Code / Cowork (MCP 클라이언트)   ·   GemsPilot runner (litellm) │
                └───────────────┬─────────────────────────────────┬──────────────────┘
                                │ MCP stdio                        │ python tool-calling
                ┌───────────────▼──────────────┐   ┌───────────────▼──────────────┐
                │ gemspilot-mcp (18 tools)      │   │ gemspilot.runner.default_toolset│
                │  + dorgems.pilot tools (P-GP-2)│   │  + dorgems.pilot TOOLSETS      │
                └───────────────┬──────────────┘   └───────────────┬──────────────┘
                                └──────────────┬───────────────────┘
                                 ┌─────────────▼─────────────┐
                                 │ dorgems  (결정론적 커널)    │
                                 │  db/  models/  kinetics/   │
                                 │  gems/ validate/ inverse/  │
                                 └──┬──────────┬──────────┬───┘
            read-only              │          │          │  reaction_model_config / materials_config
   ┌───────────────────────────┐   │          │   ┌──────▼─────────────────────────────┐
   │ scm_dor_enriched.db (RO)  │◄──┘          │   │ inverse_gems.api  (forward / design)│
   │ dorgems_staging.sqlite(RW)│◄─────────────┘   │   └─ XGEMSRunner ─ xGEMS/GEMS3K      │
   └───────────────────────────┘                  │      (또는 MockXGEMSRunner)          │
                                                  └────────────────────────────────────┘
```

### 3.2 레포 레이아웃 (`DoRGems/`)

```
DoRGems/
  pyproject.toml            # 패키지 dorgems; extras: [kernel]=inverse-gems, [pilot]=gemspilot,mcp, [export]=arviz,pymc
  README.md  LICENSE(BSD-3)  docs/  (이 스펙, model_card.md, observable_mapping.md)
  src/dorgems/
    __init__.py             # __version__, BUNDLE_SCHEMA_VERSION
    config.py               # 환경변수·경로 해석 (DORGEMS_DB, DORGEMS_STAGING_DB, DORGEMS_BUNDLE, INVERSE_GEMS_ROOT)
    db/
      reader.py             # 읽기 전용 커넥션 (sqlite3 URI ?mode=ro), 파라미터화 쿼리만
      features.py           # build_clean.py 해석기의 import 가능한 포팅 (§4.2)
      units.py              # basis/unit 조화 (§4.3)
      phases.py             # DB phase_name ↔ Cemdata18 원시 상 이름 별칭표 (§4.4)
      analogues.py          # 유사 배합 검색 (§8.2)
    models/
      bundle.py             # 번들 로드/검증 (§5.1)
      bayes.py              # numpy 사후예측 (§5.2)
      gbm.py                # LightGBM 추론 (§5.3)
      ensemble.py           # 중요도 재가중 앙상블 (§5.4)
      ood.py                # 분포 밖 판정 (§5.5)
      export_bundle.py      # 1회성: work/bayes_idata.nc + 재학습 GBM → bundles/ (arviz/pymc 필요, extras[export])
    kinetics/
      curves.py             # stretched_exp, five_param_logistic (분율 단위)
      fit.py                # logistic_fit: 곡선→(A,B,C,D,G), 최대편차 보고
      reaction_model.py     # reaction_model_config YAML 작성 (modes: logistic_fit / native / pin), 사이드카 provenance
      registry.py           # register_scm_kinetics("dorgems_stretched_exp", …) — pilot import 시 등록
      materials_override.py # materials_config 생성 (슬롯 매핑 + alias + 산화물 덮어쓰기)
    gems/
      forward.py            # inverse_gems.api 래퍼: reaction_model_config·materials_config 전달, 예산, runner_factory 캡처 훅
      capture.py            # phase_species_moles 캡처 → bound_water_g 계산 (§8.3)
      observables.py        # run_dir → {CH, bound_water, chem_shrink, QXRD phases, unreacted clinker} (§8.3)
    validate/
      twin.py               # DB 배합 → InverseGems 레시피 변환 + 실행 (§8.4)
      compare.py            # 잔차·z·판정 (§8.5)
    inverse/
      alpha_grid.py         # α 격자 forward map (§9.2)
      likelihood.py         # 관측 우도 (§9.3)
      posterior.py          # 사전표본 중요도 재가중 / 격자 사후 (§9.4)
      staging.py            # inferred_dor 쓰기 (§11)
    pilot/
      tools.py              # ToolSpec 목록 (§10)
      mcp.py                # 독립 실행 MCP 서버 `dorgems-mcp` (GemsPilot 패치 전 폴백)
      schemas.py            # pydantic 입력 스키마 (scm, mix, observation)
    cli.py                  # dorgems predict | export | materials | analogues | run-forward | compare | infer | stage | model-card
  configs/
    slots.yaml              # role → InverseGems 슬롯 매핑 규칙
    phase_aliases.yaml      # §4.4
    unit_basis.yaml         # §4.3
    analogue_tolerances.yaml# §8.2
    defaults.yaml           # α 격자, 예산 기본값, 앙상블 σ_g 등
  bundles/
    bayes_v4/  posterior.npz  scaler.json  manifest.json
    gbm_v6/    model.txt  meta.json  manifest.json
  tests/
    fixtures/  mini_scm_dor.sqlite (5편 서브셋)  golden/*.csv
    test_*.py  # markers: default(mock) / real_xgems (DORGEMS_REAL_XGEMS=1일 때만)
```

### 3.3 실행 환경 가정

- 사용자 PC는 Windows(`C:\Users\User\...`), conda 환경 `py313-xgems`에 xGEMS. 경로는 `pathlib`로만 다루고 심볼릭 링크를 쓰지 않는다.
- InverseGems는 editable 설치(`pip install -e`)를 가정한다 → `configs/`는 패키지 상대 경로로 해석된다(§1.2). CWD 기반 오버라이드에 의존하지 않는다.
- 실제 xGEMS는 `dat_lst`(`Test-dat.lst`)와 비공개 GEMS3K 파일이 있는 트리에서만 돈다. 테스트 기본은 mock.

---

## 4. 데이터 계층 (`dorgems.db`)

### 4.1 읽기 전용 리더 (`reader.py`)

- `open_ro(path) -> sqlite3.Connection`: `sqlite3.connect(f"file:{path}?mode=ro", uri=True)`, `row_factory = sqlite3.Row`. 쓰기 시도는 `sqlite3.OperationalError`로 실패해야 하며 이를 테스트한다.
- 자유 SQL은 노출하지 않는다. 이름 붙인 파라미터화 쿼리만 제공: `dor_observations(role=None, paper=None)`, `mix(mix_uid)`, `materials_for_paper(doi)`, `observations_for_mix(mix_uid, quantities, age_window)`, `opc_only_reference(quantity, age_d, tol)`, `paper(doi)`.
- 모든 결과에 `paper_doi`, `source_locator`, `extraction_confidence`, `fig_only`를 붙여 반환한다(근거 추적용).

### 4.2 특징 빌더 (`features.py`) — `build_clean.py`의 포팅

`build_clean.py`와 `multitask_v6.py`의 해석 로직을 **동일 동작**으로 함수화한다. 골든 테스트가 이를 강제한다(§13).

- `method_group(method: str|None) -> str` — 원본 `build_clean.py:14-25` 그대로.
- `role_from_text(*texts) -> str|None` — 원본 `:35-47` 그대로. `CEMENT_ROLES`, `FILLER_ROLES` 상수 유지.
- `resolve_scm(obs_row, mix_row, materials_of_paper) -> ResolvedSCM(material, role, scm_pct, how)` — 원본 `:76-124`의 fallback 사슬(phase_match → json → unique_mat → name_match → common_role → keyword_role → primary_scm_role) 그대로.
- `build_dor_table(con) -> DataFrame` — `dor_scm_final.csv`와 같은 28컬럼. 단, `build_clean.py`는 27컬럼만 만들고 `system_type`은 별도 단계에서 붙은 것(스크립트 미확인)이므로 `mixes.curing_type == 'model_system'`이면 `'model_system'`, 아니면 `'blended_cement'`로 유도한다(이 규칙으로 `system_type != 'model_system'` 필터가 1,592행/80편/476배합의 blended 기준선을 재현함을 확인했다). **골든: `scm_dor_enriched.db`에서 실행하면 `obs_uid` 집합·`dor_pct`·`scm_role`·`system_type`이 `modeling/dor_scm_final.csv`(1,610행)와 완전히 일치**해야 한다.
- `derive_composition_features(df) -> df` — `CaO_SiO2, basicity, pozz_sum, Al_Si, amorph, fineness, log_age` (`multitask_v6.py:94-99`).
- `build_aux_table(con, quantities) -> DataFrame` — 보조 관측(CH_TGA, CH_XRD, bound_water, QXRD_phase, chem_shrink, DoR_clinker)을 같은 해석기로 배합 특징에 결합. `model_system` 제외, 범위 필터 `RANGES`(`multitask_v6.py:28-29`) 적용, 여기에 `QXRD_phase (0,100)`, `chem_shrink (0, 0.2 mL/g)`, `CH_XRD (0,40)` 추가.
- `scm_input_to_features(scm: SCMSpec, mix: MixSpec) -> FeatureRow` — **새 SCM 입력을 같은 특징 공간으로 사상**하는 함수. 학습 표와 같은 결측 처리(중앙값 대치는 번들 `scaler.json`의 값 사용).

입력 스키마(`pilot/schemas.py`, pydantic, `extra="forbid"`):

```python
class SCMSpec(BaseModel):
    name: str                                  # 사용자 이름 (예: "GGBS-Pohang-2026")
    role: Literal["slag","fly_ash","metakaolin","calcined_clay","silica_fume","limestone",
                  "natural_pozzolan","glass_powder","steel_slag","other"]
    oxides: dict[str, float]                   # wt%: CaO SiO2 Al2O3 Fe2O3 MgO SO3 Na2O K2O TiO2 (LOI 선택)
    blaine_m2_kg: float | None = None
    d50_um: float | None = None
    bet_m2_g: float | None = None
    amorphous_pct: float | None = None         # 측정값만. 100−Σ결정상 계산값은 넣지 말 것
    density_kg_m3: float | None = None         # materials_config용 (없으면 슬롯 기본값)

class MixSpec(BaseModel):
    scm_pct: float                             # 결합재 중 SCM 질량% (0–100)
    w_b: float
    curing_temp_C: float = 20.0
    opc_oxides: dict[str, float] | None = None # 없으면 InverseGems OPC 기본 조성
    other_components: dict[str, float] = {}    # limestone/gypsum 등 (슬롯 이름)

class ObservationSpec(BaseModel):
    age_d: float
    quantity: Literal["CH_TGA","CH_XRD","bound_water","QXRD_phase","chem_shrink","DoR_SCM","DoR_clinker"]
    value: float
    unit: str                                  # unit_basis.yaml의 키
    phase_name: str | None = None              # QXRD_phase, DoR_clinker에서 필수
    method: str | None = None
    uncertainty: float | None = None
```

### 4.3 단위·기준 조화 (`units.py`, `configs/unit_basis.yaml`)

목표 정규 단위: 질량계 **g / 100 g 무수 결합재**, 화학수축 **mL / g 결합재**, DoR **분율**. 변환 규칙과 신뢰 등급:

| DB unit_norm / basis | 변환 (→ g/100 g binder) | 신뢰 |
|---|---|---|
| `g/100 g binder`, `g/100 g anhydrous binder`, `g_per_100g_binder` | ×1 | A |
| `g/100 g cement`, `g/100 g PC`, `g/100 g OPC`, `g/100 g clinker` | × (OPC 분율 f_OPC = 1 − scm_total_pct/100); clinker 기준은 f_OPC × (1 − 석고 분율)로 근사하고 등급 B | A / B |
| `g/100 g paste`, `g/100 g dry paste`, `g/100 g sample` | × (1 + w/b) — **건조 페이스트 기준으로 가정**(105 °C 건조 후 질량이면 ×(1 + w_bound/100) 등 반복 필요) → 등급 C, `basis_assumption` 기록 | C |
| `per_100g_ignited` | 결합수 W가 있어야 역산: ×(1 − W/100)^−1 … 순환 의존 → v1에서는 **비교 제외**, 플래그 | X |
| `mass_percent_unspecified`, `%`, `other`, NULL | 값 범위와 논문 내 다른 행의 basis로 추정 시도, 실패하면 등급 D(참고용, 잔차 통계에서 제외) | D |
| `g/100 g MK` 등 SCM 기준 | × (scm_pct/100) | B |
| QXRD `wt% (as reported)` | 그대로 두고 **상 사이의 비율만** 비교(§8.3) 또는 논문이 명시한 기준이 anhydrous면 A | C |
| chem_shrink `mL/g binder` | ×1 | A |
| chem_shrink `%` | 기준 불명(페이스트 부피? 결합재 질량?) → D | D |

`harmonize(obs, mix) -> HarmonizedObs(value_g_per_100g_binder, grade, assumptions[])`. 등급 C·D는 `compare`에서 가중 0 또는 별도 표기. **변환 규칙은 코드가 아니라 `unit_basis.yaml`에 두고 테스트로 고정한다.**

### 4.4 상 이름 별칭표 (`phases.py`, `configs/phase_aliases.yaml`)

DB `phase_name`(소문자 정규화) ↔ Cemdata18 원시 상 이름(InverseGems 출력의 `phase_mass__<RawName>`)과 InverseGems `configs/output_selection.yaml`의 그룹. **첫 실제 xGEMS 실행에서 `xgems_phase_amounts_raw.csv`의 실제 이름을 읽어 표를 확정**하고, 표에 없는 이름은 절대 추측하지 않는다(M2 게이트 항목).

| DB phase_name | Cemdata18 raw (확인 필요) | 비고 |
|---|---|---|
| portlandite, ch | `Portlandite` | CH_TGA·CH_XRD와 동일 물리량 |
| ettringite, aft | `ettringite` (+ `SO4_CO3_AFt` 합산) | |
| monocarbonate | `C4AcH11` | |
| hemicarbonate | `C4Ac0.5H12` | |
| monosulfate | `OH_SO4_AFm` 그룹(C4AsH12 계열) | 상 이름 실측 확인 |
| calcite | `Calcite` / `Cal` | 석회석 미반응분 포함 |
| hydrotalcite | `OH-hydrotalcite` | |
| strätlingite | `straetlingite` | |
| c-s-h, amorphous | `CNASH` (+ 기타 비정질) | QXRD amorphous는 미반응 유리질 SCM도 포함 → 직접 비교 금지, §8.3 |
| alite, c3s / belite, c2s / aluminate, c3a(-cubic, -ortho) / ferrite, c4af | 미반응 클링커: `unreacted_masses_g[OPC] × Bogue 분율 × (1−α_phase)` | `input_reaction_degrees.json`의 `opc` α 사용 |
| gypsum, bassanite, anhydrite | 초기 황산염 → 28 d 이상에서는 0 기대 | |
| periclase, quartz, hematite, perovskite, arcanite, aphthitalite | 비교 대상 아님(불활성/알칼리 황산염) | 무시 목록 |

---

## 5. DoR 모델 계층 (`dorgems.models`)

### 5.1 모델 번들 (`bundle.py`, `export_bundle.py`)

학습 스크립트는 그대로 두고, **추론에 필요한 것만** 번들로 동결한다. 추론 시 PyMC·ArviZ 의존을 없앤다(numpy·lightgbm만).

`bundles/bayes_v4/`
- `posterior.npz`: 사후표본을 체인 병합 후 균일 시닝하여 **S = 2,000** 개. 배열: `a0_role (S,3)`, `t0_role (S,3)`, `beta_a (S,4)`, `beta_t (S,4)`, `sd_paper_amax (S,)`, `sigma_method (S,n_m)`. 원본 `work/bayes_idata.nc`에서 추출.
- `scaler.json`: `feats = ["scm_pct","w_b","curing_temp_C","CaO_SiO2"]`, 각 특징의 `median`(대치용), `mean`, `std`(표준화용) — **학습 표(`dor_scm_final.csv`의 `age_d>0` 행)에서 재계산**(`bayes_hier_v4.py:50-60`과 동일 순서: 중앙값 대치 → 평균/표준편차). `roles = ["fly_ash","other","slag"]`(pd.Categorical 정렬 순서와 일치하는지 검증), `methods = [...]`(method_group 카테고리 순서), `beta_shape = 0.5`.
- `manifest.json`: 학습 DB 파일명·sha256, 행 수(1,610 또는 재학습 시 실제값), PyMC 버전, 수렴 진단(max r̂, min ESS), LOPO 지표(R², MAE, 90% 포함률), 생성 일시, dorgems 버전.

`bundles/gbm_v6/`
- `model.txt`: LightGBM booster(`Booster.save_model`). **재학습 규정:** `dor_scm_blended.csv`(model_system 제외) 전체 행, 특징 `BASE + COMP + CAT`, `TIGHT`, `random_state=42`, 타깃 `dor_pct`(z-표준화 없이 % 단위 직접 회귀 — 단일 태스크이므로). 범주형 카테고리 목록을 `meta.json`에 저장(추론 시 같은 코드로 인코딩).
- `meta.json`: 특징 순서, 카테고리 사전, LOPO 지표(R² 0.53, MAE 10.3 %p → `sigma_point_pct = 12.0` 기본), 학습 행 수.

`bundle.load(path) -> Bundle`은 스키마 버전·해시를 검증하고 실패 시 명시적 예외.

### 5.2 베이지안 사후예측 (`bayes.py`)

```python
def predict_curve(bundle, x: FeatureRow, role: str, ages: np.ndarray, *,
                  new_study: bool = True, method_group: str | None = None,
                  rng_seed: int = 0) -> CurveDraws
```
알고리즘(`bayes_hier_v4.py:146-163`의 OOF 계산과 동일하게):
1. `role`을 {slag, fly_ash} 외에는 `other`로 사상하고 `ood_flags`에 `role_pooled_as_other` 기록.
2. `x`를 `scaler.json`으로 대치·표준화 → `xs (4,)`.
3. 각 사후표본 s: `eta_s = a0[s,role] + xs·beta_a[s]`; `new_study`면 `u_s ~ N(0, sd_paper[s])` 1개(곡선 전체가 같은 오프셋을 공유); `a_max_s = 100·sigmoid(eta_s + u_s)`; `tau_s = exp(t0[s,role] + xs·beta_t[s])`; `alpha_s(t) = a_max_s·(1 − exp(−(t/tau_s)^0.5))`.
4. 반환 `CurveDraws(a_max (S,), tau (S,), alpha (S,len(ages)), sigma_obs (S,) | None)` — `method_group`이 주어지면 `sigma_method[s, m]`도 반환(관측 기준 구간용).
5. 요약: 재령별 `mean, q05, q50, q95`를 **잠재 곡선**(σ 없음)과 **관측**(σ 포함) 두 가지로. 값은 [0,100]으로 클립.

골든 테스트: 평균 조건(`xs = 0`), `new_study=False`에서 `a_max`의 사후평균이 `modeling/bayes_role_kinetics.csv`의 `a_max_pct`와 ±0.5 %p, `tau_days`와 ±2% 이내.

### 5.3 GBM 추론 (`gbm.py`)

```python
def predict_points(bundle, x: FeatureRow, role: str, ages: np.ndarray, *,
                   method_group: str = "missing") -> np.ndarray  # (len(ages),) %, clipped [0,100]
```
- GBM은 **측정법 조건부** 예측이다(`method_group` 특징). 새 SCM에는 기본 `"missing"`을 쓰고, 사용자가 특정 측정법 기준을 원하면 지정. 결과 메타에 `method_group_used` 기록.
- 세부 역할(metakaolin, calcined_clay, silica_fume …)은 GBM 카테고리에 그대로 들어간다(베이지안보다 세분화됨). 학습에 없는 역할은 `'missing'`으로 대치하고 플래그.

### 5.4 앙상블 (`ensemble.py`) — GBM 앵커 중요도 재가중

목적: 파라메트릭 곡선(베이지안)의 형태·불확실성은 유지하면서 점예측 정확도가 더 좋은 GBM 정보를 반영한다.

```
anchors  t_k ∈ {3, 7, 28, 90, 180} d          (defaults.yaml)
g_k      = GBM 예측 (%),  σ_g = sigma_point_pct (기본 12 %p; 재령별 상수)
w_s     ∝ Π_k N(g_k | alpha_s(t_k), σ_g²)       (사후표본별 중요도 가중)
ESS     = (Σw)² / Σw²
```
- `ESS ≥ 100`이면 가중 분위수로 요약하고 `(a_max, tau)` 표본도 가중 리샘플링(S개 유지)한다. `ESS < 100`이면 **재가중을 포기하고 베이지안 단독 결과**를 내며 `warnings += ["model_disagreement: GBM anchors far from Bayesian prior (ESS=…)"]`.
- 모드: `ensemble ∈ {"blend"(기본), "bayes", "gbm_anchor_only"}`. `gbm_anchor_only`는 GBM 앵커에 stretched-exp를 직접 피팅(불확실성 없음, 진단용).
- **한계를 문서화한다.** 두 모델은 같은 데이터로 학습되어 독립 정보가 아니다. blend는 휴리스틱이며, **M1 게이트 G1-3**(§12)에서 LOPO OOF(`bayes_oof.csv` + GBM OOF)로 blend가 MAE를 낮추고 90% 포함률을 0.85 이상 유지하는지 확인해 기본값을 확정한다. 실패하면 기본값을 `bayes`로 바꾼다.

### 5.5 분포 밖 판정과 유사사례 (`ood.py`, `db/analogues.py`)

- 역할별 학습 분포에서 `(scm_pct, w_b, curing_temp_C, CaO_SiO2, Al_Si, basicity)`의 1–99 백분위 범위와 표준화 마할라노비스 거리를 번들에 저장. 범위 밖 특징 목록 + 거리 백분위를 `ood_flags`, `ood_score`로 반환.
- 역할별 학습 표본 수가 적으면(`n_papers < 5`: metakaolin 8편 근처, calcined_clay 4, silica_fume 4, limestone 3) `sparse_role` 플래그.
- `find_analogues(scm, mix, k=5)`: 같은 역할 + 표준화 화학 거리 + |Δscm_pct| + |Δw_b| 가중합으로 가까운 **배합** k개와 그 DoR 관측(재령별), 논문 DOI, 측정법을 반환 → 예측 결과에 `evidence`로 첨부. 이것이 시나리오 A의 근거 제시이자 (보류한) 시나리오 D의 최소 구현이다.

### 5.6 예측 결과 스키마 (`prediction.json`)

```json
{
  "schema": "dorgems-prediction/1.0",
  "input": {"scm": {...}, "mix": {...}, "ages_d": [1,3,7,28,90,180,365]},
  "features": {"scm_pct": 40, "w_b": 0.45, "curing_temp_C": 20, "CaO_SiO2": 1.05, "imputed": ["curing_temp_C"]},
  "role_bayes": "slag", "role_gbm": "slag",
  "bayes":   {"a_max": {"q05":..,"q50":..,"q95":..}, "tau_d": {...},
              "alpha_pct": {"latent": {"q05":[..],"q50":[..],"q95":[..]}, "observed": {...}}},
  "gbm":     {"alpha_pct": [..], "method_group_used": "missing", "sigma_point_pct": 12.0},
  "ensemble":{"mode": "blend", "ess": 412.3, "a_max": {...}, "tau_d": {...}, "alpha_pct": {...}},
  "recommended": {"source": "ensemble", "alpha_pct_q50": [..], "alpha_pct_q05": [..], "alpha_pct_q95": [..]},
  "ood": {"flags": [], "score_pct": 37.2, "sparse_role": false},
  "evidence": {"analogues": [{"mix_uid": "...", "paper_doi": "...", "scm_pct": 40, "w_b": 0.5,
                              "dor": [{"age_d": 28, "value_pct": 33.1, "method_group": "selective_dissolution"}]}]},
  "warnings": [],
  "provenance": {"dorgems": "0.1.0", "bundle_bayes": "sha256:…", "bundle_gbm": "sha256:…", "seed": 0}
}
```

---

## 6. InverseGems 연동 계층 (`dorgems.kinetics`, `dorgems.gems`)

### 6.1 반응모델 내보내기 (`reaction_model.py`)

`export_reaction_model(prediction, out_dir, *, mode, slot, quantiles=(0.05,0.5,0.95), config_id) -> {q: path}`

세 모드. **모두 α는 분율**로 쓴다(퍼센트 금지 — §1.2).

**mode = `logistic_fit`(기본, 호환성 최대).** 분위 곡선(q05/q50/q95) 각각에 `five_param_logistic`을 피팅한다.
- 격자: `t ∈ logspace(log10 0.25, log10 730, 48)` d, 가중 균일(log t).
- 초기값·경계: `A∈[0,0.2]`, `B∈[0.05,5]`, `C∈[0.1,500]`, `D∈[0.05,1]`, `G∈[0.05,5]` (커널 `kinetics_calibration._DEFAULT_BOUNDS`와 동일). `A`는 0으로 고정하는 옵션 기본 on(stretched-exp는 t→0에서 0).
- 보고: `max_abs_dev_pct` on [1, 365] d. **> 2 %p면 warning, > 5 %p면 실패**(mode `native` 권고).
- YAML:
  ```yaml
  id: dorgems_<config_id>_q50
  scm_reaction:
    slag: {A: 0.0, B: 0.62, C: 24.1, D: 0.585, G: 1.31}
  availability_modifier: {enabled: false}     # 이유: 치환율(scm_pct)·CH 가용성 효과는 이미 DoR 모델에 포함 → 이중 계산 방지
  ```
  **provenance는 사이드카 `<id>.provenance.json`에만 쓴다.** `load_reaction_parameters`(`reaction_parameters.py:117-140`)는 미지 최상위 키를 무시하지만 파일 전체를 `raw_config`로 보존해 `signature_payload()`에 넣으므로(`:35-55, 132`), YAML에 provenance 텍스트를 넣으면 같은 파라미터라도 `reaction_model_signature`가 달라진다. YAML은 파라미터만 담는다.

**mode = `native`(인프로세스 전용).** `dorgems.kinetics.registry`가 import 시 `register_scm_kinetics("dorgems_stretched_exp", required=("a_max","tau","beta"), asymptote_key="a_max")`로 `alpha = a_max·(1 − exp(−(t/tau)^beta))`(a_max 분율)를 등록. YAML은 `slag: {model: dorgems_stretched_exp, a_max: 0.585, tau: 18.5, beta: 0.5}`. `inverse-gems` CLI 단독으로는 열 수 없으므로 YAML 헤더 주석에 명시하고, `dorgems` CLI/pilot 경로에서만 사용.

**mode = `pin`(시나리오 C용).** `{A: α, B: 1, C: 1, D: α, G: 1}` + `availability_modifier: {enabled: false}`. α 하나당 YAML 하나(`…_pin_a0.350.yaml`).

세 분위 config를 InverseGems에 각각 돌리면 상조성·공극률의 **불확실성 봉투**가 나온다 — 이것이 시나리오 A의 표준 산출물이다.

### 6.2 재료 오버라이드 (`materials_override.py`, `configs/slots.yaml`)

새 SCM은 InverseGems 슬롯 4개 중 하나로 **사상**된다(`SCM_NAMES` 고정 — §1.2):

| SCMSpec.role | 슬롯 | 비고 |
|---|---|---|
| slag, steel_slag | `slag` | steel_slag는 OOD 플래그 |
| fly_ash, natural_pozzolan, glass_powder | `fly_ash` | 저CaO 포졸란 계열 |
| metakaolin, calcined_clay | `metakaolin` | calcined_clay는 MK 함량 미지 → 경고 |
| silica_fume | `silica_fume` | |
| limestone | `limestone` | 반응도 모델 대상 아님(충전재 취급) |
| other | 화학 거리 최소 슬롯 | 반드시 경고 |

`build_materials_config(scm, slot, alias, cement=None, base=<InverseGems configs/materials.yaml>) -> path`: 기본 materials.yaml을 읽어 **해당 슬롯의 `oxide_mass_percent`(및 density)를 SCMSpec 값으로 덮어쓰고**, `aliases:`에 사용자 이름(`scm.name` 정규화)을 추가해 `materials.dorgems_<hash>.yaml`로 저장. 산화물 합이 100 ± 3이 아니면 정규화하고 경고(LOI 제외). OPC 조성이 주어지면 `OPC` 항목도 덮어쓴다.

주입 방법 두 가지:
- **정식(P-IG-1 이후):** `run_forward_request(..., materials_config=path)`.
- **폴백(패치 전):** `dorgems.gems.forward`가 같은 프로세스에서 **`inverse_gems.materials.config_path`**를 컨텍스트 매니저로 임시 교체(monkeypatch)해 `materials.yaml` 요청만 오버라이드 파일로 돌린다. `load_materials` 자체를 패치하면 안 된다 — 각 모듈이 `from .materials import load_materials`로 import 시점에 바인딩하므로 `cached_forward`·`forward_query` 등에 미치지 않는다(`config_path`는 `load_materials` 내부에서 호출 시점에 해석됨). 결과 `manifest.json`에 `materials_injection: "monkeypatch"`를 기록하고 warning을 낸다. 실측 검증: 실행 후 `input_materials_used.json`이 오버라이드 값을 담고 있는지 확인(assert).

캐시 안전성: 산화물이 바뀌면 원소 벡터가 바뀌어 `chem_hash`가 달라진다(§1.2). `reaction_model_signature`는 기본 `configs/materials.yaml`을 해시하는데(`reaction_model.py:16-26 DEFAULT_SIGNATURE_FILES`, materials.yaml은 22행) **reaction_model_config YAML의 `signature_files:` 키(`reaction_model.py:70`)나 `extra_signature_files` 인자로 오버라이드 파일을 시그니처에 추가할 수 있어 커널 패치가 필요 없다.** 내보내는 YAML마다 `signature_files: [<materials override path>]`를 넣는다(이 키는 파라미터가 아니라 시그니처 구성 지시이므로 §6.1의 provenance 제외 원칙과 충돌하지 않는다).

### 6.3 Forward 래퍼 (`gems/forward.py`)

```python
def run_forward(forward_query: dict, *, out, db, reaction_model_config, materials_config=None,
                use_mock=True, dat_lst=None, max_xgems_calls=None, capture_species=False) -> ForwardRunResult
```
두 경로를 가진다. 어느 쪽이든 `materials_config`는 §6.2 규칙으로 주입한다.

- **표준 경로(`capture_species=False`):** forward_query dict를 `out/forward_query.yaml`로 **먼저 파일로 쓴 뒤**(커널 `run_forward_request`의 `forward_query`는 dict가 아니라 YAML 파일 경로다 — `forward_query.py:407-408`; GemsPilot의 `_materialize_query`와 같은 처리) `inverse_gems.api.run_forward_request(forward_query=path, out=…, db=…, reaction_model_config=…, use_mock=…, dat_lst=…, max_xgems_calls=…, disable_plots=True)`를 호출한다.
- **캡처 경로(`capture_species=True`, 결합수 계산용 §8.3):** `runner_factory`는 `run_forward_cached`(`cached_forward.py:222`)만 받고 `run_forward_query`/`run_forward_request`로는 스레딩되지 않는다. 따라서 재령별로 `run_forward_cached(recipe_text="OPC 60, slag 40, w/b 0.45, age 28", db=…, dat_lst=…, use_mock=…, runner_factory=CapturingRunnerFactory, reaction_model_config=…, recipe_id=…)`를 직접 호출하고, 반환 dict의 `chemistry_dir`에 `dorgems_capture.json`(`phase_species_moles`, `species_in_phase`, `species_molar_masses`, `phase_masses`)을 추가로 쓴다. `CapturingRunner`는 `XGEMSRunner`(또는 mock)를 감싸 `capture_raw_state()` 결과를 가로챈다. `runner_factory` 호출 시그니처: `runner_factory(dat_lst_path=…, temperature_celsius=…)` (`cached_forward.py:35, 396`). 시계열은 DoRGems가 재령별 결과를 모아 `time_series.csv`와 같은 컬럼 규약으로 합성한다. P-IG-2가 머지되면 표준 경로만 쓴다.
- 실행 후 **자기검증:** `input_reaction_degrees.json["scm"][slot]`이 내보낸 곡선의 해당 재령 값과 1e-3 이내인지 확인(아니면 config가 적용되지 않은 것 → 실패). `materials_config`를 썼으면 `input_materials_used.json`의 슬롯 산화물이 오버라이드 값과 일치하는지 확인.
- 반환에 GemsPilot ToolResult 호환 dict + `time_series.csv` 경로 + 재령별 `phase_mass__*`, `scalar__porosity`, `scalar__pH`.

### 6.4 업스트림 패치 목록 (각각 별도 PR, 최소 침습)

| ID | 레포 | 내용 | 필요 시점 |
|---|---|---|---|
| **P-GP-1** | GemsPilot | `agent_tools.run_forward / run_task / run_confirmed_query / run_design_with_recovery` 및 MCP 래퍼에 `reaction_model_config: str\|None`, `reaction_model_id: str\|None`, `materials_config: str\|None` 키워드 추가 → 커널에 전달. `design_recovery.py`의 `reaction_model_config=None` 하드코딩을 인자화. 경로는 `_resolve_artifact_path`로 검사. | M1 |
| **P-GP-2** | GemsPilot | 툴셋 플러그인 훅: `runner.default_toolset(extra: Iterable[ToolSpec] = ())` + `importlib.metadata.entry_points(group="gemspilot.toolsets")` 탐색; `mcp_server`도 같은 훅으로 도구 등록. DoRGems는 `pyproject`에 `[project.entry-points."gemspilot.toolsets"] dorgems = "dorgems.pilot.tools:TOOLSET"` 선언. | M1 (없으면 `dorgems-mcp` 독립 서버로 폴백) |
| **P-GP-3** | GemsPilot | `_policy_check`에 `write` 정책 추가: `use_mock`이 아닌 `dry_run`(기본 True) 인자를 보고 `allow_real`로 게이트. 현재는 `use_mock` 인자가 없는 도구가 정책과 무관하게 통과한다(§1.3). | M3 (폴백: 쓰기 도구에 `use_mock` 인자를 dry-run 의미로 둠) |
| **P-IG-1** | InverseGems | `materials_config: str\|Path\|None` 인자를 `run_request / run_forward_request / run_forward_query / run_forward_cached / run_task_query / run_chemistry_design_query`에 스레딩; 6곳의 `load_materials()`를 `load_materials(materials_config)`로(`materials or load_materials()` 폴백 7곳은 상위에서 `materials=`로 내려주면 자동 해결). 시그니처 포함은 `signature_files:` 키로 이미 가능하므로 범위 밖. | M1 (폴백: `config_path` monkeypatch) |
| **P-IG-2** | InverseGems | 원시 캡처에 `xgems_phase_species_moles_raw.json`(상별 종 몰), `species_molar_masses` 저장; 선택적으로 `scalar__bound_water_g`(고체 내 H2O 질량) 계산·저장. 아울러 `runner_factory`를 `run_forward_query`/`run_forward_request`까지 스레딩. | M2 (폴백: `run_forward_cached` + `capture_species` 훅) |
| **P-IG-3** | InverseGems | `build_xgems_input(..., clinker_degrees: dict[str,float]\|None)`로 PK를 우회하는 명시적 클링커 α 입력. | M3 이후(선택) — DoH 역해석 확장용 |
| **P-IG-4** | InverseGems | 부피 단위(cm³/m³) 불일치 감사(`porosity.yaml` vs `backfill.py:329`, `cli.py:921`); 실측 xGEMS 단위로 통일. | M2 게이트에서 결과에 따라 |
| P-IG-5 | InverseGems | `custom_scm_1..2` 슬롯 추가(두 집합 + 하드코딩 리스트 14곳 이상, §1.2). | 후순위 |

---

## 7. 시나리오 A — 새 SCM → DoR 예측 → InverseGems 실행

### 7.1 흐름

```
SCMSpec + MixSpec + ages
   │  dor_predict                       (read)
   ▼
prediction.json  ── evidence(analogues), ood, q05/q50/q95 곡선
   │  dor_export_reaction_model         (read)   mode=logistic_fit, 3분위
   ▼
reaction_parameters.dorgems_<id>_q{05,50,95}.yaml (+ provenance.json)
   │  dor_build_materials_override      (read)
   ▼
materials.dorgems_<hash>.yaml            (슬롯 산화물 덮어쓰기 + alias)
   │  dor_run_forward_with_dor          (mock_ok; 실제 실행은 allow_real)  ×3 분위
   ▼
run_dir/{q05,q50,q95}/forward/time_series.csv  ── 상조성·공극률·pH 봉투
   │  (선택) GemsPilot run_task(design_query.reaction_model.config=q50)  → inverse design
```

### 7.2 forward_query 생성 규칙

- `recipe.binders`: `{OPC: 100 − scm_pct − Σother, <slot>: scm_pct, **other_components}` (합 100). `w_b`, `temperature_celsius = curing_temp_C`, `age_grid.values = ages`.
- `response_summary.phases`는 **실측 상 이름 목록**(§4.4 확정 후)만 사용. mock 실행에서는 `Mock …` 이름으로 자동 대체.
- `material_system`은 지정하지 않는다(라우팅 필터 회피). design_query로 갈 때만 `design_space.allowed_materials`에 슬롯 포함.

### 7.3 산출물

`out/A_<id>/`: `prediction.json`, `reaction_models/*.yaml`, `materials.*.yaml`, `runs/q05|q50|q95/…`, `envelope.csv`(재령 × {phase, porosity, pH} × {q05,q50,q95}), `summary.md`(숫자는 전부 CSV에서 렌더), `manifest.json`.

### 7.4 GemsPilot과의 관계

- 사용자가 DoR 측정 CSV를 이미 갖고 있으면 기존 `calibrate_scm_kinetics`가 맞는 도구다. DoRGems는 **측정이 없는 새 SCM** 또는 **문헌 사전(prior)이 필요한 경우**를 맡는다. 두 결과의 비교(`dor_compare_reaction_models`)를 도구로 제공한다: 같은 재령 격자에서 두 config의 α 차이.
- design_query에는 `reaction_model: {config: <q50 path>, mismatch_policy: …}`로 전달(커널 `DesignQuerySpec.reaction_model`, `design_query.py:136-150`). 이 경로는 task_query YAML에 담기면 **현재 GemsPilot `run_task`로도 이미 통과**한다(§1.3) — P-GP-1은 forward·diagnose 경로용이다.

---

## 8. 시나리오 B — xGEMS forward 결과를 문헌 관측으로 검증

### 8.1 두 가지 모드

- **`twin`(기본, 정확):** DB의 실제 배합을 골라 **그 배합을 InverseGems 레시피로 재구성해 실행**하고, 같은 배합·같은 재령의 관측과 1:1 비교한다. 재료 조성은 DB `materials`에서 오버라이드(§6.2), DoR은 (a) 그 배합의 측정 DoR_SCM이 있으면 pin, (b) 없으면 DoR 모델 q50. 클링커 α는 PK.
- **`neighbourhood`(느슨):** 사용자의 run_dir(임의 배합)을 유사 배합(§8.2)의 관측 분포와 비교. 배합이 정확히 같지 않으므로 등급을 한 단계 낮춘다.

### 8.2 유사 배합 선택 (`analogues.py`, `configs/analogue_tolerances.yaml`)

기본 허용: 같은 슬롯/역할, `|Δscm_pct| ≤ 10 %p`, `|Δw_b| ≤ 0.05`, `|ΔT| ≤ 5 °C`(결측 T는 20 °C로 간주 + 플래그), 재령 비 `t_obs/t_run ∈ [0.75, 1.33]`, `curing_type ≠ model_system`, 화학 거리(표준화 `CaO_SiO2, Al_Si, basicity`) ≤ 1.5. 가중 `w = exp(−d²/2)`. 결과에 `n_mixes, n_papers`를 반드시 표기(1편이면 `single_paper` 플래그).

### 8.3 관측량 사상 (`gems/observables.py`) — xGEMS 출력 → DB 물리량

시스템 기준은 결합재 100 g. xGEMS 질량(kg) × 1000 = g/100 g binder.

| DB quantity | 모델 측 계산 | 비교 가능성 |
|---|---|---|
| **CH_TGA, CH_XRD** | `phase_mass__Portlandite` × 1000 | A. 가장 직접적. TGA CH는 탄산화 보정 여부(method_detail)에 따라 하향 편향 가능 → 플래그 |
| **bound_water** | `W_bound = W_in − W_aq`, `W_aq` = `aq_gen` 상의 H2O 질량(`dorgems_capture.json`의 상별 종 몰 × 몰질량; 없으면 `aq_gen` 상 질량 − 용존 이온 질량 근사). 대안 정의 `W_solids` = 고체 상들의 H2O-당량 몰 × 18.015 (동일해야 함; 불일치 시 경고) | B. TGA 결합수(105–1000 °C 또는 105–550 °C)는 105 °C에서 일부 겔/층간수를 잃음 → **계통 오프셋** 예상. OPC 단독 참조군(393배합)으로 오프셋 `b_BW(t)`를 추정해 `compare`의 기준선으로 씀 |
| **QXRD_phase** (결정상) | 별칭표(§4.4)로 `phase_mass__<raw>` 합산 × 1000. 미반응 클링커: `unreacted_masses_g["OPC"] × f_Bogue,phase × (1 − α_phase)` (`input_reaction_degrees.json`, `input_recipe.json`, Bogue는 커널 `bogue_phases`) | C. `wt% (as reported)` 기준 불명 → 논문이 anhydrous 기준을 명시한 경우만 절대 비교(A), 아니면 **상 간 비율**(예: 에트린자이트/포틀란다이트, 모노카보네이트/헤미카보네이트)과 재령 추세 방향만 비교 |
| QXRD `amorphous` | 비교 금지. 미반응 유리질 SCM + C-S-H + 미검출상 혼합 | X |
| **chem_shrink** | `CS = V_initial − (V_solid_final + V_aq)` (cm³/100 g binder) ÷ 100 → mL/g binder. `porosity.json`의 `initial_volume_cm3, solid_final_volume_cm3`, `excluded_non_solid_phase_volumes_raw[aq_gen]` 사용. **부피 단위(cm³/m³) 실측 확인 필수(P-IG-4)** | B (mL/g binder 44건만) |
| DoR_clinker | `input_reaction_degrees.json["opc"]` (PK) vs 관측 | 참고용: PK 상수 적합성 진단 |
| cum_heat | **v1 범위 밖.** 후속: Bogue 상별 반응열 × α (Parrot-Killoh 경로) 또는 GEMS 엔탈피 | — |

### 8.4 twin 실행 (`validate/twin.py`)

`db_mix_to_recipe(mix_row, materials_of_paper) -> (forward_query, materials_config, warnings)`
- `binder_composition_json`의 material_id를 `materials.role`로 슬롯에 사상(cement/clinker → OPC, gypsum → gypsum, limestone → limestone, SCM 역할 → §6.2 표). 사상 불가 성분이 5% 초과면 배합 제외.
- OPC 산화물은 DB `materials(role='cement')` 값이 완전하면(CaO SiO2 Al2O3 Fe2O3 SO3 필수, 241건 가능) 오버라이드, 아니면 기본값 + 플래그. Blaine은 `fineness_m2_kg`로(있으면; 74건).
- w/b, curing_temp_C 필수(없으면 제외).
- DoR: `pin`(측정치 보간: 관측 재령 사이는 stretched-exp 피팅, 재령 3개 미만이면 모델 q50) 또는 `model`.
- 예산: 배합당 재령 수만큼 xGEMS 호출. `max_xgems_calls`로 상한. 캐시 적중은 예산 미소모(커널 규칙).

### 8.5 비교·판정 (`validate/compare.py`)

관측 하나마다 `r = model − obs_harmonized`(정규 단위), `σ_obs` = 관측 `uncertainty`가 있으면 그것, 없으면 물리량별 기본(CH 1.5 g/100 g, BW 1.5 g/100 g, CS 0.005 mL/g; `defaults.yaml`), `σ_model` = 물리량별 모델 불일치 항(초기값: OPC 참조군 잔차 표준편차에서 추정, M2에서 확정). `z = r / sqrt(σ_obs² + σ_model²)`.

집계(등급 A·B 관측만): `n, n_papers, median_r, IQR_r, frac_|z|<2, bias_sign_test_p`. 판정(결정론적 임계, `defaults.yaml`):
- `consistent`: `frac_|z|<2 ≥ 0.7` 그리고 `|median_r| ≤ σ_model`
- `tension`: 그 외, `n ≥ 5`
- `insufficient_data`: `n < 5` 또는 등급 A·B 관측 없음

출력 `comparison.json` + `comparison.csv`(관측 단위 행: obs_uid, quantity, age, obs, model, r, z, grade, assumptions) + `summary.md`. **판정 문구는 템플릿, 숫자는 CSV에서만.**

### 8.6 OPC 단독 참조 검사(모델 자체 점검)

SCM이 없는 364–393배합(79–84편; §1.1의 상한 주의 — binder JSON에 SCM이 있는 29배합은 제외)은 DoR 모델과 무관하게 **커널(Cemdata18 + PK + 물 정책)의 기준 성능**을 보여준다. M2에서 먼저 실행: 28 d, w/b 0.4–0.5, 20 ± 3 °C, CH_TGA 등급 A 관측에 대해 `phase_mass__Portlandite`의 잔차 분포를 보고 `σ_model[CH]`, `b_BW`를 확정한다. 이 결과는 InverseGems 자체의 검증 자료로도 가치가 있다(별도 노트로 남김).

---

## 9. 시나리오 C — 관측치 → DoR 역해석 루프

### 9.1 문제 정의

주어진 배합(MixSpec + SCMSpec 또는 DB mix_uid)과 재령별 관측 `{(t_i, q, y_iq)}`(CH, 결합수, 화학수축, QXRD 상)에서 SCM 반응도 곡선 `α(t) = a_max(1 − exp(−(t/τ)^0.5))`의 사후분포를 구한다. 사전분포는 시나리오 A의 DoR 모델(사후표본 `(a_max_s, τ_s)`), 또는 `flat`. 클링커 α는 v1에서 PK 고정(σ_model에 흡수; P-IG-3 후 확장).

### 9.2 α-격자 forward map (`inverse/alpha_grid.py`) — 핵심 아이디어

관측 재령이 정해지면 xGEMS가 필요한 것은 **각 재령에서 α의 함수로서의 관측량** `F_q(α; t_i)`뿐이다. 곡선 파라미터 `(a_max, τ)`는 그 다음 단계에서 xGEMS 호출 없이 평가된다.

```
for t_i in observed ages:
    for α_j in grid (기본 0.00, 0.05, …, 1.00 → 21점; 사전 90% 구간 안은 0.025 간격으로 세분)
        config = pin(α_j)                      # §6.1 mode=pin
        run_forward(recipe @ t_i, config)      # 캐시: 같은 (원소벡터, T) 재사용
        F_q[i, j] = observables(run_dir)[q]    # §8.3
    F_q(·; t_i) ← 단조(PCHIP) 보간
```
- 호출 수 ≤ 21 × n_ages(≤ 6) ≈ 126, 캐시로 감소. 예산 `max_xgems_calls` 필수.
- 물리적 단조성 점검: CH는 α에 대해 비증가여야 한다(포졸란 반응이 CH 소비). 위반 시 해당 격자점 플래그(솔버 문제 가능).
- 결과 `forward_map.npz` + `forward_map_report.md`.

### 9.3 우도 (`inverse/likelihood.py`)

관측 `y_iq`(정규 단위, 등급 A·B만), 모델 `μ_iq = F_q(α(t_i); t_i) + b_q`:
```
log L(a_max, τ) = Σ_iq  log N(y_iq | μ_iq, σ_iq²),   σ_iq² = σ_obs,iq² + σ_model,q²
```
- `b_q`: 물리량별 계통 오프셋. CH는 0 고정, bound_water는 §8.6의 `b_BW`(사전 N(b̂, sd)로 주변화 — 격자 5점), chem_shrink 0.
- DoR_SCM 관측이 입력에 포함되면 **역해석 대상이 아니라 검증용**으로만 쓰고 우도에서 제외(옵션 `use_direct_dor=False` 기본).

### 9.4 사후 (`inverse/posterior.py`)

1. 사전 = A의 사후표본 `(a_max_s, τ_s)`, S=2,000 (모드 `prior="model"`); `prior="flat"`이면 `a_max ~ U(0.02, 1.0)`, `log τ ~ U(log 0.5, log 3000)` 격자 40×40.
2. 중요도 가중 `w_s ∝ L(a_max_s, τ_s)`; `ESS = (Σw)²/Σw²`. `ESS < 50`이면 사전표본을 중심으로 가우시안 커널로 표본을 2배 확장(SIR 1회)한 뒤 재계산; 그래도 낮으면 격자 사후로 전환하고 `posterior_method` 기록.
3. 요약: `(a_max, τ)` 가중 분위, 재령 격자별 α 분위, 관측별 사후예측 잔차(`ppc.csv`), `prior_vs_posterior_kl`(정보 이득 지표), 우도 기여도 분해(어느 관측이 결론을 끌었는가).
4. 출력 `inference.json`(schema `dorgems-inference/1.0`), `reaction_parameters.dorgems_inferred_<id>_q{05,50,95}.yaml`(§6.1 logistic_fit), `inferred_dor.csv`(`scm, age_d, dor` — 커널 `calibrate_scm_kinetics` 호환 포맷).

### 9.5 피드백 경로

- **InverseGems로:** 위 YAML을 시나리오 A와 같은 방식으로 forward/design에 전달.
- **DB로:** `dorgems_staging.sqlite.inferred_dor`에 행 추가(§11). 문헌 표에는 절대 쓰지 않는다. 검토 승인(`reviewed=1`)된 행만 향후 재학습 스크립트가 `method_group='inverse_thermodynamic'`으로 읽어 갈 수 있다 — 재학습 자체는 이 스펙 범위 밖.
- **DoR 모델 갱신 판단 근거:** 역해석 결과와 모델 사전의 KL, 그리고 `ESS`를 누적 기록해 어느 SCM 역할·조성 영역에서 모델이 체계적으로 틀리는지 리포트(`staging` 뷰).

### 9.6 검증 데이터셋(내장)

- 합성: mock 및 실제 xGEMS로 알려진 `(a_max*, τ*)`에서 CH·BW·CS를 생성하고 노이즈를 더해 복원(20 케이스, 역할·치환율·w/b 조합).
- 실측: DoR_SCM과 CH/BW가 ≥3 공통 재령에 있는 **66배합**. DoR_SCM은 숨기고 CH/BW/(QXRD)로 역해석 → 숨긴 DoR와 비교(§12 G3).

---

## 10. 도구 계약 (`dorgems.pilot.tools`)

모든 도구는 GemsPilot `ToolResult` 계약을 그대로 반환한다(`contract: "inverse-gems-tool/1.0"`, `tool`, `ok`, `summary`, `artifacts`, `warnings`, `error`). 쿼리 인자는 GemsPilot 관례처럼 dict·YAML 문자열·파일 경로를 모두 받는다(`_load_query_payload` 재사용 또는 동일 구현). `out/db/session` 경로는 GemsPilot 러너의 워크스페이스 재매핑을 그대로 통과해야 하므로 **키워드 이름을 `out`, `db`, `session`으로 통일**한다.

| 도구 | 시그니처(요약) | 정책 | 내용 |
|---|---|---|---|
| `dor_predict` | `(scm, mix, ages, out, *, ensemble="blend", method_group=None, seed=0, session=None)` | read | §5 → `prediction.json` |
| `dor_export_reaction_model` | `(prediction, out, *, mode="logistic_fit", slot=None, quantiles=[0.05,0.5,0.95], config_id=None)` | read | §6.1 → YAML ×3 + provenance |
| `dor_build_materials_override` | `(scm, out, *, slot=None, alias=None, cement=None)` | read | §6.2 |
| `dor_find_analogues` | `(scm, mix, *, age_days=None, quantities=None, limit=20)` | read | §5.5 / §8.2, 근거 표 |
| `dor_run_forward_with_dor` | `(forward_query, reaction_model_config, out, db, *, materials_config=None, use_mock=True, dat_lst=None, max_xgems_calls=None, capture_species=False, session=None)` | mock_ok | §6.3 |
| `dor_run_envelope` | `(scm, mix, ages, out, db, *, use_mock=True, dat_lst=None, max_xgems_calls=None, session=None)` | mock_ok | 시나리오 A 전체(예측→내보내기→3분위 실행→envelope.csv) 한 번에 |
| `dor_compare_to_literature` | `(out, db, *, run_dir=None, mix_uid=None, scm=None, mix=None, mode="twin", quantities=None, use_mock=True, dat_lst=None, max_xgems_calls=None, session=None)` | mock_ok | §8 |
| `dor_opc_reference_check` | `(out, db, *, age_days=28, w_b_range=[0.4,0.5], use_mock=True, dat_lst=None, max_xgems_calls=None)` | mock_ok | §8.6 (mock에서는 파이프라인 점검만 의미 있음) |
| `dor_infer_from_observations` | `(mix, observations, out, db, *, scm=None, mix_uid=None, prior="model", alpha_grid=21, use_mock=True, dat_lst=None, max_xgems_calls=None, session=None)` | mock_ok | §9 |
| `dor_stage_inferred` | `(inference, staging_db, *, use_mock=True, note=None)` | mock_ok | §11 — `use_mock=True`면 미리보기 JSON만 생성(dry-run), `False`면 staging에 기록. GemsPilot 정책 검사가 `use_mock`만 보기 때문에 이 이름을 쓴다(P-GP-3 후 `dry_run`으로 개명) |
| `dor_compare_reaction_models` | `(config_a, config_b, *, ages=None)` | read | 두 YAML의 α(t) 차이표 |
| `dor_model_card` | `()` | read | 번들 메타·지표·학습 DB 해시 |
| `dor_db_lookup` | `(query_name, params, *, limit=50)` | read | §4.1의 이름 붙인 쿼리만 |

- `dorgems.pilot.tools.TOOLSET: list[ToolSpec]` — GemsPilot P-GP-2 엔트리포인트용. `ToolSpec(name, func, policy)` 형식은 `gemspilot.runner`에서 import(없으면 동일 구조의 로컬 dataclass로 폴백).
- `dorgems.pilot.mcp`: `dorgems-mcp` 콘솔 스크립트. FastMCP/MCPServer 버전 분기는 GemsPilot `mcp_server.py`의 `_ServerApp` 패턴을 그대로 따른다. 모든 실행 도구 기본 `use_mock=True`.
- **정책 표기 규칙:** `read` = 계산·조회만, 부작용 없음. `mock_ok` = 기본 mock, `use_mock=False`는 호스트 `allow_real`에서만 허용(GemsPilot `_policy_check` 의미론). `real_gated`는 사용하지 않는다(현재 구현에서 "실제 실행 불가"를 뜻함).
- **도구 설명문(docstring)은 LLM이 읽는다.** 각 도구에 언제 쓰고 언제 쓰지 말지, 실제 실행이 예산·승인을 요구함을 명시. 예: `dor_infer_from_observations`는 "측정 DoR이 이미 있으면 `calibrate_scm_kinetics`를 쓸 것".

---

## 11. 가드레일과 거버넌스

- **읽기 전용 문헌 DB.** `DORGEMS_DB` 경로는 `?mode=ro`로만 연다. 쓰기 가능한 커넥션을 반환하는 코드 경로가 없어야 하며, 테스트에서 INSERT가 실패함을 확인한다.
- **Staging DB** (`DORGEMS_STAGING_DB`, 기본 `<project>/dorgems_staging.sqlite`):
  ```sql
  CREATE TABLE inferred_dor (
    inf_uid TEXT PRIMARY KEY,           -- dorgems::<inference_id>::<age>
    inference_id TEXT, created_at TEXT, dorgems_version TEXT, bundle_hash TEXT,
    mix_uid TEXT,                       -- 문헌 배합이면 FK 의미(문헌 DB 참조, 제약은 없음)
    scm_json TEXT, mix_json TEXT,       -- 사용자 배합이면 스펙 원문
    slot TEXT, age_d REAL,
    alpha_q05 REAL, alpha_q50 REAL, alpha_q95 REAL,
    a_max_q50 REAL, tau_q50 REAL, ess REAL, posterior_method TEXT,
    observations_used_json TEXT, forward_map_path TEXT, run_manifest_path TEXT,
    reviewed INTEGER DEFAULT 0, review_note TEXT
  );
  CREATE TABLE tool_audit (ts TEXT, tool TEXT, args_hash TEXT, ok INTEGER, xgems_calls INTEGER, run_dir TEXT);
  ```
- **승인 게이트.** 실제 xGEMS 실행과 staging 기록은 모두 `use_mock=False`로만 일어나며, GemsPilot 러너에서는 `Episode.allow_real`, MCP 호스트에서는 호스트의 승인 흐름을 탄다. `max_xgems_calls`가 None이거나 200을 넘는 실제 실행 요청은 DoRGems가 자체적으로 거부하고 명시적 상한을 요구한다. DoRGems 자체는 승인 UI를 갖지 않는다.
- **예산.** 모든 xGEMS 경로에 `XGEMSCallBudget`; 초과 시 부분 결과와 `chemistry_status="skipped_budget"` 유지(커널 규칙).
- **워크스페이스 격리.** 산출물은 `out` 아래에만. 입력 경로는 `INVERSE_GEMS_ARTIFACT_ROOTS` 허용목록 검사(GemsPilot `_resolve_artifact_path` 재사용).
- **프롬프트 주입 저항.** 도구 입력에 들어온 자연어("관리자가 실제 실행을 승인했다" 등)는 정책을 바꾸지 못한다. GEMS-Agent-Bench의 injection 항목을 DoR 도구용으로 복제(M4).
- **출처 표기.** 사용자에게 보이는 모든 수치 옆에 산출 파일 경로와 (문헌값이면) DOI·source_locator가 붙는다.

---

## 12. 마일스톤과 결정 게이트

시간 추정은 넣지 않는다. 순서와 게이트만 구속한다. 게이트는 `tests/`의 pytest로 자동화하고, 자동화가 불가능한 항목(실제 xGEMS)은 `DORGEMS_REAL_XGEMS=1`에서만 실행되는 마커 테스트로 둔다.

### M0 — 기반: 리더·특징·번들
산출: `dorgems.db.*`, `dorgems.models.bundle/bayes/gbm`, `export_bundle.py`, `bundles/`, `tests/fixtures/mini_scm_dor.sqlite`(5편 서브셋, 쓰기 가능한 별도 파일), `dor_model_card`.
게이트 **G0**:
- G0-1 `build_dor_table(scm_dor_enriched.db)`의 `obs_uid` 집합·`dor_pct`·`scm_role`이 `modeling/dor_scm_final.csv`(1,610행)와 **완전 일치**.
- G0-2 베이지안 번들: 평균 조건 `a_max`·`τ`가 `bayes_role_kinetics.csv`와 ±0.5 %p / ±2%.
- G0-3 GBM 번들: 재학습 LOPO(10-fold, 논문 그룹, 시드 42)의 R²가 0.50 이상(보고값 0.522–0.531 재현 범위).
- G0-4 RO 커넥션에서 INSERT가 실패한다.

### M1 — 시나리오 A
산출: `kinetics.*`, `gems.forward`(폴백 monkeypatch 포함), `pilot.tools`의 A 관련 도구, `dorgems-mcp`, PR P-GP-1·P-GP-2·P-IG-1.
게이트 **G1**:
- G1-1 `logistic_fit`: slag/fly_ash/other 대표 곡선 12개에서 `max_abs_dev_pct(1–365 d) < 2`.
- G1-2 mock forward 실행 후 `input_reaction_degrees.json["scm"][slot]` == 내보낸 α(±1e-3), `input_materials_used.json`이 오버라이드 산화물을 담음.
- G1-3 **앙상블 기본값 결정:** `bayes_oof.csv`와 GBM OOF(같은 fold 규칙 재생성)로 blend vs bayes를 LOPO 비교. blend가 MAE를 낮추고 90% 포함률 ≥ 0.85면 `blend` 유지, 아니면 기본값 `bayes`로 변경하고 문서 갱신.
- G1-4 (real_xgems) 사용자 PC `py313-xgems`에서 OPC 60 / slag 40, w/b 0.45, 28 d 실제 실행 성공, q05/q50/q95 봉투 산출. GemsPilot 회귀 앵커(OPC60/slag40 w/b 0.45 28 d porosity 0.398451, `phase_mass__CNASH` 0.045535 — 기본 kinetics 기준)와 **기본 config로 돌린 결과가 일치**함을 먼저 확인한 뒤 DoRGems config로 차이를 기록.
- G1-5 GemsPilot MCP에서 `dor_predict → dor_export_reaction_model → run_forward(reaction_model_config=…)`가 한 에피소드로 완료(P-GP-1·2 머지 전이면 `dorgems-mcp` 경로).

### M2 — 시나리오 B
산출: `units`, `phases`(실측 확정), `gems.capture/observables`, `validate.*`, `dor_compare_to_literature`, `dor_opc_reference_check`, PR P-IG-2(·P-IG-4).
게이트 **G2**:
- G2-1 실제 xGEMS 1회 실행의 `xgems_phase_amounts_raw.csv` 상 이름으로 `phase_aliases.yaml` 확정; 표에 없는 이름 참조 시 예외.
- G2-2 부피 단위 실측(cm³/m³) 확정, `chem_shrink` 계산이 OPC 페이스트 28 d에서 문헌 범위(0.04–0.07 mL/g binder 부근)에 들어옴.
- G2-3 OPC 참조 검사(§8.6, 등급 A CH_TGA, ≥ 30편): `phase_mass__Portlandite` 잔차 `|median_r| ≤ 4 g/100 g`. 초과하면 **중단·보고**(DoR 문제가 아니라 커널/PK/DB 기준 문제) — 게이트 실패는 정보다.
- G2-4 결합수 오프셋 `b_BW` 추정치와 불확실성 확정, `σ_model` 표 확정.
- G2-5 twin 모드로 ≥ 20 SCM 배합(측정 DoR pin)에서 `comparison.json` 생성, 판정 분포 리포트.

### M3 — 시나리오 C
산출: `inverse.*`, `dor_infer_from_observations`, `dor_stage_inferred`, staging 스키마.
게이트 **G3**:
- G3-1 합성 복원(mock 20 케이스): 진짜 `α(28 d)`가 사후 90% 구간에 ≥ 85% 포함, 사후 중앙값 MAE ≤ 5 %p.
- G3-2 (real_xgems) 합성 복원 5 케이스 동일 기준.
- G3-3 실측 66배합(또는 단위 등급 A·B 관측이 있는 부분집합): 숨긴 DoR_SCM 대비 사후 중앙값 MAE와 90% 포함률 보고. **결정:** MAE가 GBM LOPO(10.3 %p)보다 낮거나 포함률 ≥ 0.8이면 C를 정식 기능으로 승격; 둘 다 아니면 C는 "일관성 점검·불확실성 축소" 용도로 격하하고 §9.1의 클링커 α 고정(P-IG-3)을 다음 원인 후보로 기록.
- G3-4 사전→사후 정보 이득(KL)이 관측 수와 함께 단조 증가하는지 확인(우도가 실제로 작동).

### M4 — 피드백·벤치·문서
산출: 검토 CLI(`dorgems review list|approve|reject`), staging 리포트 뷰, GemsPilot 벤치 시나리오 `dor_*`(정상 경로, 예산 초과, 문헌 DB 쓰기 요청 거부, 주입 저항), `docs/model_card.md`, `docs/observable_mapping.md`.
게이트 **G4**: 벤치 항목 전부 통과(mock), 문서에 §1의 사실 재검증 일자 기록.

---

## 13. 테스트 전략

- **마커:** 기본(mock, 네트워크·xGEMS 불필요) / `real_xgems`(`DORGEMS_REAL_XGEMS=1`이고 `dat_lst` 존재 시) / `slow`(GBM 재학습, LOPO).
- **골든:** `tests/golden/` — `dor_scm_final.csv` 서브셋(fixture DB 5편에 대응), `bayes_role_kinetics.csv`, 대표 곡선 12개의 logistic 피팅 결과, `unit_basis` 변환 표 케이스.
- **속성 테스트:** 곡선 단조성, α∈[0,1], 분위 순서(q05 ≤ q50 ≤ q95), 산화물 정규화 합 100, RO 위반 예외.
- **회귀 앵커:** GemsPilot `configs/agent_qa_generated.yaml`의 frozen 값(OPC100 w/b 0.5 28 d porosity 0.407184, pH 12.661534 등)을 real_xgems 마커에서 재현해 **커널 환경이 정상인지**를 먼저 확인.
- **합성 복원:** §9.6.
- **CI:** mock 경로만. real_xgems는 로컬 수동.

---

## 14. 알려진 리스크와 미결 사항

1. **단위·기준의 이질성**이 B·C의 유효 표본을 크게 깎는다(유효 CH_TGA 행 중 `mass_percent_unspecified` 1,790건, bound_water 1,071건). 등급 D 관측을 살리려면 논문별 basis 재판정 작업(DB 측)이 필요 — 이 스펙 범위 밖, 별도 백로그.
2. **결합수 정의 불일치**(TGA 105 °C 기준 vs 열역학 결합수). 오프셋 모델로 흡수하되 재령 의존성이 있을 수 있음.
3. **QXRD 기준 불명**으로 절대 비교는 소수 논문만 가능. 비율 비교로 후퇴.
4. **클링커 DoH를 고정할 수 없음**(PK만). C에서 SCM α와 클링커 α의 교락이 CH 우도에 들어간다. P-IG-3이 해법.
5. **앙상블의 이중 계산.** G1-3으로 실증 결정.
6. **OOD 역할**(MK·CC·SF·LS)은 베이지안에서 `other`로 풀링되어 구간이 매우 넓다. 정직한 결과지만 사용자 기대와 다를 수 있음 → 도구 설명과 출력에 명시.
7. **측정법 편차는 식별 불가**(위치 편차 없음, 노이즈 척도만). B의 판정에서 측정법 차이를 오프셋으로 넣지 않는다.
8. **부피 단위 불일치**(cm³/m³) — M2 첫 실행에서 실측.
9. **`register_scm_kinetics`는 프로세스 로컬** — native 모드는 CLI 상호운용성이 없음. 기본은 logistic_fit.
10. **InverseGems 슬롯 4개 한계** — 새 역할은 슬롯에 사상되며 이름은 alias로만 보존. custom 슬롯은 후순위(P-IG-5).
11. **cum_heat(17,637건)은 v1에서 미사용.** 가장 큰 보조 데이터인데 열역학 forward에서 직접 나오지 않는다. 후속: 클링커 상별 반응열 기반 추정 → 별도 스펙.
12. **재학습 파이프라인은 범위 밖.** 이 스펙은 기존 모델을 서빙·활용한다. 7차 보고서의 "보조값을 특징으로 쓰는 2단 모델"은 별도 트랙.

---

## 15. 부록

### 15.1 참조 파일

| 항목 | 경로 |
|---|---|
| 현행 DB | `DoR of SCMs in blended cements/modeling/scm_dor_enriched.db` |
| 모델링 표 | `modeling/dor_scm_final.csv`(1,610), `modeling/dor_scm_blended.csv`(1,592), `modeling/dor_scm_blended_oof.csv` |
| 표 빌더·해석기 | `modeling/build_clean.py` |
| 베이지안 v4 | `modeling/bayes_hier_v4.py`, `modeling/work/bayes_idata.nc`(폴더 목록에 없음 → 재실행 예상), `modeling/bayes_oof.csv`, `modeling/bayes_role_kinetics.csv`, `modeling/bayes_feature_effects.csv`, `modeling/method_noise.csv` |
| GBM v6 | `modeling/multitask_v6.py`(DoR-only 경로만 사용) |
| 보고서 | `modeling/DoR_모델링_3차보고서.md`(모델 구조·교락 논증), `DoR_6차보고서.md`, `DoR_7차보고서.md`(멀티태스크 결론), `SCM별_OPC함량_반응도.md`(OPC 효과), `BASELINE_PRE_AUX.md` |
| InverseGems | `src/inverse_gems/{api.py, forward_query.py, cached_forward.py, xgems_input_builder.py, scm_reaction.py, reaction_parameters.py, availability_modifier.py, materials.py, recipe.py, chem_hash.py, xgems_runner.py, database.py, porosity.py, kinetics_calibration.py}`, `configs/{materials.yaml, scm_reaction.yaml, reaction_parameters.example.yaml, porosity.yaml, output_selection.yaml, species_map.yaml, formula_map.yaml}` |
| GemsPilot | `src/gemspilot/{agent_tools.py, mcp_server.py, runner.py, design_recovery.py, agent_bench.py}`, `configs/{agent_qa_generated.yaml, models.yaml}`, `tests/test_agent_tools.py` |

### 15.2 InverseGems 핵심 시그니처(발췌)

```python
inverse_gems.api.run_forward_request(*, forward_query,  # ← YAML 파일 경로(dict 아님)
    out, db, dat_lst=None, use_mock=False,
    run_mode="reacted_only", gems_class_path="xgems:ChemicalEngineDicts", xgems_input_mode="formula",
    xgems_water_mode="initial", ..., max_xgems_calls=None,
    reaction_model_id=None, reaction_model_config=None, disable_plots=False, fail_fast=False) -> RequestResult
inverse_gems.cached_forward.run_forward_cached(*, recipe_text, db, dat_lst, use_mock, ..., runner_factory=None,
    recipe_id, recipe_metadata, reaction_model_id, reaction_model_config)
inverse_gems.scm_reaction.register_scm_kinetics(name, *, required, asymptote_key)
inverse_gems.scm_reaction.scm_alpha(age_days, params)               # 분율, clip [0,1]
inverse_gems.reaction_parameters.load_reaction_parameters(path=None, *, reaction_model_id=None)
inverse_gems.xgems_input_builder.build_xgems_input(recipe, *, materials=None, species_map=None,
    run_mode="reacted_only", opc_phase_mass_percent=None, temperature_celsius=20.0, relative_humidity=None,
    fineness_m2_kg=None, apply_availability_modifier=None, reaction_parameters=None, xgems_water_mode="initial",
    xgems_water_factor=1.0, xgems_water_g=None, xgems_water_w_b=None) -> XGEMSInput
inverse_gems.kinetics_calibration.calibrate_scm_kinetics(*, data_csv, out, model="five_param_logistic",
    config_id=None, scms=None, fixed_params=None, param_init=None, param_bounds=None, min_points=8, make_plot=True)
```

### 15.3 GemsPilot 핵심 시그니처(발췌)

```python
gemspilot.agent_tools.run_forward(forward_query, out, db, *, use_mock=False, dat_lst=None, disable_plots=True,
    retry_water_on_failure=False, retry_water_policy="diagnosis", max_xgems_calls=None, session=None) -> dict
gemspilot.agent_tools.calibrate_scm_kinetics(data_csv, out, *, model="five_param_logistic", config_id=None,
    scms=None, session=None) -> dict
gemspilot.runner.ToolSpec(name, func, policy: str, description="")   # policy ∈ {"read","mock_ok","real_gated"}; 검증 없음(runner.py:59-64)
gemspilot.runner.default_toolset() -> list[ToolSpec]
gemspilot.runner.Episode(model, workspace, allow_real=False, protocol="full", max_steps=12, temperature=0.0, ...)
```

### 15.4 forward_query 예시(mock)

```yaml
name: dorgems_A_example
task: forward_time_series
recipe: {binders: {OPC: 60, slag: 40}, w_b: 0.45}
age_grid: {values: [1, 3, 7, 28, 90, 365]}
temperature_celsius: 20.0
outputs: {phase_masses: all, phase_volumes: all, phase_volumes_reconstructed: all, aqueous_species: all, scalars: all}
plots: []
response_summary: {enabled: true, phases: [Portlandite, CNASH, ettringite, C4AcH11], scalars: [pH, porosity], narrative_language: ko}
```

---

## 16. 코딩 세션 시작 프롬프트(붙여넣기용)

```
당신은 DoRGems 구현 세션이다. 작업 명세는 `DoRGems_agent_spec_v1.md`이며 그 규칙(§0)을 따른다.

1. 먼저 §1의 파일:라인 사실을 현재 코드에서 재확인하고 달라진 점을 문서 상단 "재검증 로그"에 기록한다.
2. M0을 구현한다: 레포 스캐폴드(§3.2), 읽기 전용 리더, build_clean.py의 함수 포팅, 번들 내보내기, 골든 테스트 G0-1..G0-4. 
   - 문헌 DB(`modeling/scm_dor_enriched.db`)는 절대 수정하지 않는다.
   - `modeling/work/bayes_idata.nc`는 현재 없을 가능성이 크다. 없으면 `bayes_hier_v4.py`를 `modeling/` 디렉터리에서 그대로 재실행해 만든다(입력 `dor_scm_final.csv`, 시드 42, PyMC 필요). 스크립트는 수정하지 않는다. 재실행 결과의 `bayes_role_kinetics.csv`가 기존 파일과 ±0.5 %p 이내인지 먼저 확인한다.
3. G0가 전부 통과하면 M1로 간다. 업스트림 패치(P-GP-1, P-GP-2, P-IG-1)는 각각 별도 브랜치·PR로 만들고, 머지 전에는 폴백 경로로 동작시킨다.
4. 게이트마다 결과를 `docs/gates.md`에 표로 남긴다(통과/실패/수치). 게이트 실패는 숨기지 말고 원인 가설과 함께 보고한 뒤 멈춘다.
5. 모든 실행 도구의 기본은 use_mock=True다. 실제 xGEMS 실행은 사용자가 명시적으로 허용한 경우에만, DORGEMS_REAL_XGEMS=1로 수행한다.
```
