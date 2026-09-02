# DoRGems model card (bundles v1, 2026-09-02)

## Data

| item | value |
|---|---|
| source | `DoR of SCMs in blended cements/modeling/scm_dor_enriched.db` (read-only) |
| modelling table | `dor_scm_final.csv` — 1,610 rows / 80 papers (reproduced exactly by `dorgems.db.features.build_dor_table`, G0-1) |
| blended baseline | 1,592 rows / 80 papers / 476 mixes (`system_type != model_system`) |
| roles (rows) | fly_ash 755, slag 551, metakaolin 118, calcined_clay 57, limestone 46, silica_fume 37, other 25, steel_slag 10, other_scm 6, natural_pozzolan 5 |
| roles (papers) | slag 41, fly_ash 37, metakaolin 8, calcined_clay 4, silica_fume 4, limestone 3, other 2, natural_pozzolan 1, steel_slag 1, other_scm 1 |
| missingness | BET 91 %, d50 90 %, amorphous 76 %, Blaine 63 %, chemistry 16–20 %, w/b 19 %, T 20 % |

## bayes_v4 (hierarchical Bayesian, PyMC 5.28)

`alpha(t) = a_max (1 − exp(−(t/τ)^0.5))`, `logit(a_max/100) = a0[role] + X·β_a + u_paper`,
`log τ = t0[role] + X·β_t`, `y ~ N(alpha, σ[method])`, `X = z(scm_pct, w_b, curing_temp_C, CaO/SiO2)`
with median imputation. Roles: `slag`, `fly_ash`, everything else pooled as `other`.

Re-run 2026-09-02 on the 1,610-row table (the archived `bayes_role_kinetics.csv`
in `modeling/` came from the earlier 1,177-row table — see gates.md G0-2):

| role | n_obs | a_max % (mean, 94 % HDI) | τ d |
|---|---|---|---|
| fly_ash | 755 | 30.6 [25.7, 36.4] | 30.7 |
| other | 304 | 39.0 [32.9, 45.7] | 5.2 |
| slag | 551 | 56.9 [51.2, 62.9] | 13.3 |

Convergence: max r̂ 1.010, min ESS 860. Method noise σ (%p): selective dissolution 6.3,
XRD 6.8, unknown 7.6, PONKCS 9.5, mass balance 10.4, other 10.6, SEM-BSE 12.8,
calorimetry 14.3, NMR 14.3. Location bias by method is *not* identifiable (method ≈ paper).

Leave-papers-out (5 folds): point R² 0.421 / MAE 11.35 %p; predictive R² 0.419 /
MAE 11.54 %p; 90 % interval coverage 0.906 (mean width 46 %p).

Bundle: 2,000 posterior draws (uniform thinning of 6,000), `scaler.json`, sha256 manifest.

## gbm_v6 (LightGBM 4.7, DoR-only)

Features `log_age, scm_pct, w_b, curing_temp_C, CaO_SiO2, basicity, pozz_sum, Al_Si, amorph, fineness`
+ categorical `scm_role, method_group`; `TIGHT` hyper-parameters; raw-% target; seed 42;
trained on the 1,592-row blended table. LOPO (10 paper folds): R² 0.503, MAE 10.74 %p
(seed 7: 0.450, seed 2024: 0.515). `sigma_point_pct = 12`.

## Ensemble

Default **`bayes`**. The GBM-anchored importance re-weighting (`blend`) was tested
leave-papers-out (docs/g1_3): MAE 11.37 vs 11.43 %p but 90 % coverage 0.77 vs 0.91 —
it narrows intervals without earning it, because both models share the training data.
`blend` and `gbm_anchor_only` remain available as explicit options.

## Known limitations

1. Roles other than slag/fly_ash are pooled → wide intervals; `ood.sparse_role` flags roles
   with < 5 papers (metakaolin 8 is the only non-sparse one).
2. Measurement-method bias is absorbed into paper effects; only noise scales differ.
3. InverseGems has four fixed SCM slots; new SCMs are mapped by role/chemistry and the
   user's name survives only as an alias (`materials.dorgems_<hash>.yaml`).
4. Clinker hydration is Parrot–Killoh only; scenario C absorbs it into σ_model.
5. No real xGEMS run was possible on the build machine; every physical comparison
   (M2 gates, phase aliases, volume units) is still to be confirmed (see gates.md).
