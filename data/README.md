# data/

`scm_dor_enriched.db` — the literature database DoRGems is trained on and validates
against (SQLite, 24 MB, snapshot of 2026-08-28; sha256 `de071b0c…`, the same file
recorded as `training_db_sha256` in `bundles/*/manifest.json`).

| table | rows | content |
|---|---|---|
| papers | 1,350 | doi, title, year, journal, retrieval route |
| materials | 1,634 | per-paper materials: role, oxides (wt%), Blaine, d50, BET, amorphous % |
| mixes | 6,870 | binder composition (JSON, material_id → mass %), w/b, curing |
| observations | 44,011 | age, quantity (DoR_SCM, DoR_clinker, CH_TGA, CH_XRD, bound_water, QXRD_phase, chem_shrink, cum_heat), value, unit, basis, method, source locator, extraction confidence |

Rules (spec §0, §11): DoRGems opens this file **read-only** (`?mode=ro`, no free SQL);
agent-produced values go to `dorgems_staging.sqlite`, never here. Discovery order:
`DORGEMS_DB` → this file → `../DoR of SCMs in blended cements/modeling/`.

Known limitation: `unit_norm` is normalised to `g/100 g binder` even when
`basis_reported` is `mass_percent_unspecified`; DoRGems grades such rows D (reference only).
See `configs/unit_basis.yaml` and `docs/gates.md` G2-3.
