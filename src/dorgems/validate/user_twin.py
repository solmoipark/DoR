"""User-supplied material + observations → forward run at the observation ages → 1:1
comparison (the ``dorgems compare --input`` path). Uses the DoR model q50 for the SCM
(scenario A) unless the input carries measured DoR_SCM at ≥ 3 ages, in which case the
measured curve is pinned like the literature twin."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from ..db.units import harmonize
from ..envelope import build_forward_query
from ..gems.forward import run_forward
from ..gems.observables import observables_from_run
from ..kinetics.materials_override import build_materials_config, slot_for_role
from ..kinetics.reaction_model import export_reaction_model
from ..predict import predict
from .compare import OBS_TO_MODEL, aggregate, compare_rows, write_comparison
from .twin import dor_pin_curve


def user_compare(
    scm: Any,
    mix: Any,
    observations: list[Any],
    *,
    out: str | Path,
    ig_db: str | Path,
    use_mock: bool = True,
    dat_lst: str | Path | None = None,
    max_xgems_calls: int | None = None,
    lit_db: str | Path | None = None,
) -> dict[str, Any]:
    from ..pilot.schemas import coerce_mix, coerce_observations, coerce_scm

    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    scm_m, mix_m = coerce_scm(scm), coerce_mix(mix)
    obs = coerce_observations(observations)
    obs_cmp = [o for o in obs if o.quantity in OBS_TO_MODEL]
    if not obs_cmp:
        return {"ok": False, "error": "no comparable observations (CH_TGA, CH_XRD, bound_water, chem_shrink)"}
    ages = sorted({float(o.age_d) for o in obs_cmp})
    warnings: list[str] = []
    slot, w, _ = slot_for_role(scm_m.role, scm_m.oxides)
    warnings += w
    mat = build_materials_config(scm_m, out, slot=slot, cement=mix_m.opc_oxides)
    warnings += mat["warnings"]
    # DoR source: measured curve if ≥3 ages, else model q50
    dor_rows = pd.DataFrame([{"age_d": o.age_d, "dor_pct": o.value if (o.unit or "").strip() != "fraction" else o.value * 100} for o in obs if o.quantity == "DoR_SCM"])
    curve, w2 = dor_pin_curve(dor_rows, np.asarray(ages, float)) if not dor_rows.empty else (None, [])
    warnings += w2
    if curve is not None:
        from scipy.optimize import least_squares

        from ..kinetics.curves import stretched_exp

        t = np.asarray(ages, float)
        r = least_squares(lambda p: stretched_exp(t, p[0], np.exp(p[1]), 0.5) - curve, x0=[float(curve.max()), np.log(20.0)], bounds=([0.01, np.log(0.1)], [1.0, np.log(5000.0)]))
        pred = {"id": "user_pin", "input": {"ages_d": ages}, "beta_shape": 0.5, "bayes": {"a_max": {"q50": float(r.x[0]) * 100}, "tau_d": {"q50": float(np.exp(r.x[1]))}}, "recommended": {"source": "bayes"}, "provenance": {"source": "measured DoR_SCM pinned"}}
        dor_source = "pin"
    else:
        pred = predict(scm_m, mix_m, ages, db_path=lit_db)
        dor_source = f"model_q50({pred['recommended']['source']})"
    rm = export_reaction_model(pred, out / "rm", mode="logistic_fit", slot=slot, quantiles=(0.5,), config_id=pred["id"], signature_files=[mat["path"]])
    fq = build_forward_query(mix_m, slot, ages, name="user_twin")
    res = run_forward(fq, out=out / "run", db=ig_db, reaction_model_config=rm["q50"]["path"], materials_config=mat["path"], slot=slot, use_mock=use_mock, dat_lst=dat_lst, max_xgems_calls=max_xgems_calls, capture_species=True)
    warnings += res.warnings
    if not res.ok:
        return {"ok": False, "error": res.error or "forward run failed", "warnings": warnings}
    model = observables_from_run(res.forward_dir)
    pairs = []
    for i, o in enumerate(obs_cmp):
        h = harmonize({"quantity": o.quantity, "value": o.value, "unit": o.unit, "unit_reported": o.unit, "basis_reported": None}, {"scm_total_pct": mix_m.scm_pct, "w_b": mix_m.w_b}, scm_pct=mix_m.scm_pct)
        row = model[model["age_d"] == float(o.age_d)]
        col = OBS_TO_MODEL[o.quantity]
        mv = float(row.iloc[0][col]) if (not row.empty and row.iloc[0][col] is not None and not pd.isna(row.iloc[0][col])) else None
        pairs.append({"obs_uid": f"user_{i}", "paper_doi": "user", "mix_uid": scm_m.name, "quantity": o.quantity, "phase_name": o.phase_name, "age_d": float(o.age_d), "method": o.method, "grade": h.grade, "assumptions": "; ".join(h.assumptions), "obs_value": h.value, "model_value": mv, "uncertainty": o.uncertainty})
    df = compare_rows(pairs)
    agg = aggregate(df)
    files = write_comparison(df, agg, out, header={"mode": "user_twin", "target": scm_m.name, "dor_source": dor_source, "use_mock": use_mock, "slot": slot, "materials_injection": res.materials_injection})
    files["run_dir"] = str(res.run_dir)
    (out / "manifest.json").write_text(json.dumps({"scm": scm_m.model_dump(), "mix": mix_m.model_dump(), "dor_source": dor_source, "slot": slot, "use_mock": use_mock, "warnings": warnings}, indent=2, default=str), encoding="utf-8")
    return {"ok": True, "slot": slot, "dor_source": dor_source, "n_obs": len(pairs), "aggregate": agg, "files": files, "warnings": warnings}
