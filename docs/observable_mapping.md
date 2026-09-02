# Observable mapping: literature DB ↔ xGEMS output

Basis: 100 g anhydrous binder. Real xGEMS phase masses are kg → ×1000 = g/100 g binder; the mock runner reports grams. `observables.mass_factor()` detects the unit from `scalar__system_mass` (> 10 ⇒ grams) and must be re-checked on the first real run (P-IG-4).
Implemented in `dorgems.gems.observables`, aliases in `configs/phase_aliases.yaml`
(**provisional** — `confirmed: false` until a real `xgems_phase_amounts_raw.csv` is read, G2-1),
units in `configs/unit_basis.yaml` (`dorgems.db.units`).

| DB quantity | model side | grade | notes |
|---|---|---|---|
| CH_TGA, CH_XRD | `phase_mass__Portlandite` × 1000 | A | TGA CH may be biased low without carbonation correction (method_detail) |
| bound_water | `W_in − W_aq`; `W_aq` from `aq_gen` species moles captured in-process (`dorgems_capture.json`), else aqueous phase mass (approx.) | B | TGA bound water (105 °C) loses gel/interlayer water → systematic offset `b_BW(t)`; estimated from the OPC-only reference set (G2-4, pending) |
| QXRD_phase (crystalline) | Σ `phase_mass__<raw>` of the alias group × 1000; unreacted clinker from `unreacted_masses.json` × Bogue × (1 − α) | C | `wt% (as reported)` basis unknown → ratios and trends only |
| QXRD amorphous | — | X | contains unreacted glassy SCM + C-S-H + undetected phases |
| chem_shrink | `(V_initial − (V_solid_final + V_aq)) / 100` from `porosity.json` | B | cm³ vs m³ unit to be confirmed on a real run (P-IG-4) |
| DoR_clinker | `reaction_degrees.json["opc"]` (PK) | ref | diagnostic for the Parrot–Killoh constants |
| cum_heat | — | — | out of scope in v1 |

## Unit / basis grades

A exact (`g/100 g binder`); A/B cement / clinker basis with `f_OPC = 1 − scm_total/100`
(clinker × (1 − 0.05)); C paste basis (`× (1 + w/b)`, dry-paste assumption);
D unspecified (`mass_percent_unspecified`, `%`) — reference only, excluded from statistics;
X ignited basis (circular) — excluded.

Counts in the current DB (valid rows): CH_TGA `g/100 g binder` 3,180 (A) vs
`mass_percent_unspecified` 1,790 (D); bound_water valid 1,071; chem_shrink `mL/g binder` 44 (A) vs `%` 488 (D).

## Kernel artefact names (re-verified 2026-09-02)

The `forward_query` / `run_forward_cached` path writes, per age,
`<db>/recipe_runs/<recipe_id>/{reaction_degrees.json, recipe.json, unreacted_masses.json, porosity.json, volumes.json, source_contribution_ledger.csv}`
and `<db>/chemistry_runs/<chem_hash>/xgems_raw/{xgems_phase_amounts_raw.csv, xgems_phase_volumes_raw.csv, xgems_scalars_raw.json, …}`.
`input_reaction_degrees.json` / `input_materials_used.json` exist only on the stand-alone
`xgems_output_capture.save_run_outputs` path. DoRGems verifies the applied reaction model from
`reaction_degrees.json["scm"]` and the applied materials from the source ledger's oxide moles.
