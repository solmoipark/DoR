"""Scenario C orchestration (spec §9): mix + observations → forward map → posterior → artefacts."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

import numpy as np

from .. import __version__
from ..db.units import harmonize
from ..envelope import build_forward_query
from ..kinetics.materials_override import build_materials_config, slot_for_role
from ..kinetics.reaction_model import export_reaction_model
from ..models.bayes import predict_curve
from ..models.bundle import load_bundle
from ..predict import load_defaults
from ..db.features import scm_input_to_features
from .alpha_grid import build_forward_map, default_alpha_grid
from .likelihood import Likelihood, build_points
from .posterior import infer, summarise, write_inference


def _hid(payload: Any) -> str:
    return hashlib.sha256(json.dumps(payload, sort_keys=True, default=str).encode()).hexdigest()[:12]


def harmonise_user_observations(observations: list[Any], mix: dict[str, Any]) -> list[dict[str, Any]]:
    out = []
    for i, o in enumerate(observations):
        d = o.model_dump() if hasattr(o, "model_dump") else dict(o)
        h = harmonize({"quantity": d["quantity"], "value": d["value"], "unit": d.get("unit"), "basis_reported": d.get("basis")}, {"scm_total_pct": mix.get("scm_pct"), "w_b": mix.get("w_b")}, scm_pct=mix.get("scm_pct"))
        out.append({"obs_uid": d.get("obs_uid", f"user_{i}"), "age_d": float(d["age_d"]), "quantity": d["quantity"], "phase_name": d.get("phase_name"), "value": h.value, "grade": h.grade, "assumptions": h.assumptions, "uncertainty": d.get("uncertainty"), "method": d.get("method"), "raw_value": d["value"], "unit": d.get("unit")})
    return out


def infer_from_observations(
    mix: Any,
    observations: list[Any],
    *,
    out: str | Path,
    ig_db: str | Path,
    scm: Any | None = None,
    mix_uid: str | None = None,
    lit_db: str | Path | None = None,
    prior: str = "model",
    alpha_grid_n: int = 21,
    refine: bool = True,
    use_mock: bool = True,
    dat_lst: str | Path | None = None,
    max_xgems_calls: int | None = None,
    seed: int = 0,
    b_bw_prior: tuple[float, float] | None = None,
    quantities: tuple[str, ...] = ("CH_TGA", "CH_XRD", "bound_water", "chem_shrink"),
) -> dict[str, Any]:
    from ..pilot.schemas import coerce_mix, coerce_observations, coerce_scm

    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    defaults = load_defaults()
    warnings: list[str] = []
    # --- inputs: user spec or literature mix ---
    if mix_uid is not None:
        from ..db.reader import LiteratureDB
        from ..validate.twin import db_mix_to_recipe

        dbp = Path(lit_db) if lit_db else None
        if dbp is None:
            from ..config import literature_db_path

            dbp = literature_db_path()
        with LiteratureDB(dbp) as db:
            mrow = db.mix(mix_uid)
            if mrow is None:
                raise ValueError(f"mix {mix_uid} not found")
            mats = db.materials_for_paper(mrow["paper_doi"])
            obs_rows = db.observations_for_mix(mix_uid, list(quantities))
        obs_rows = [o for o in obs_rows if o["age_d"] and o["age_d"] > 0]
        ages = sorted({float(o["age_d"]) for o in obs_rows})
        rec = db_mix_to_recipe(mrow, mats, ages=ages, out_dir=out)
        if rec.get("excluded_reason"):
            raise ValueError(f"mix {mix_uid} cannot be rebuilt: {rec['excluded_reason']}")
        warnings += rec["warnings"]
        fq_base = rec["forward_query"]
        slot = rec["slot"]
        mat_path = rec["materials_config"]
        scm_pct = float(rec.get("scm_pct") or 0.0)
        obs_h = []
        for o in obs_rows:
            h = harmonize(o, mrow, scm_pct=scm_pct)
            obs_h.append({"obs_uid": o["obs_uid"], "age_d": float(o["age_d"]), "quantity": o["quantity"], "phase_name": o.get("phase_name"), "value": h.value, "grade": h.grade, "assumptions": h.assumptions, "uncertainty": o.get("uncertainty"), "method": o.get("method")})
        m = rec["scm_material"] or {}
        ox = {k: float(m[k]) for k in ("CaO", "SiO2", "Al2O3", "Fe2O3", "MgO", "SO3", "Na2O", "K2O", "TiO2") if m.get(k) is not None}
        scm_m = coerce_scm({"name": str(m.get("name_in_paper") or "scm"), "role": m.get("role") if m.get("role") in ("slag", "fly_ash", "metakaolin", "calcined_clay", "silica_fume", "limestone", "natural_pozzolan", "glass_powder", "steel_slag") else "other", "oxides": ox or {"SiO2": 40.0, "CaO": 40.0}})
        mix_m = coerce_mix({"scm_pct": scm_pct, "w_b": float(mrow["w_b"]), "curing_temp_C": float(mrow["curing_temp_C"] if mrow.get("curing_temp_C") is not None else 20.0)})
        if mrow.get("curing_temp_C") is None:
            warnings.append("curing_temp_C missing in DB; 20 °C assumed")
    else:
        if scm is None:
            raise ValueError("scm is required when mix_uid is not given")
        scm_m, mix_m = coerce_scm(scm), coerce_mix(mix)
        obs_list = coerce_observations(observations)
        obs_h = harmonise_user_observations(obs_list, mix_m.model_dump())
        ages = sorted({o["age_d"] for o in obs_h})
        slot, w, _ = slot_for_role(scm_m.role, scm_m.oxides)
        warnings += w
        mat = build_materials_config(scm_m, out, slot=slot, cement=mix_m.opc_oxides)
        mat_path = mat["path"]
        warnings += mat["warnings"]
        fq_base = build_forward_query(mix_m, slot, ages, name="dorgems_C")
    if slot is None:
        raise ValueError("the mix has no SCM slot to infer")
    # --- prior ---
    bundle = load_bundle(require_gbm=False)
    feats = scm_input_to_features(scm_m, mix_m)
    draws = predict_curve(bundle.bayes, {k: feats.get(k) for k in bundle.bayes.feats}, scm_m.role, np.array([28.0]), new_study=True, rng_seed=seed)
    prior_a, prior_t = draws.a_max / 100.0, draws.tau
    warnings += [f"prior:{f}" for f in draws.flags]
    a28 = draws.alpha[:, 0] / 100.0
    lo, hi = float(np.quantile(a28, 0.05)), float(np.quantile(a28, 0.95))
    grid_cfg = defaults.get("alpha_grid", {})
    alphas = default_alpha_grid(int(alpha_grid_n or grid_cfg.get("n_points", 21)), refine_interval=(lo, hi) if (prior == "model" and refine) else None, refine_step=float(grid_cfg.get("refine_step", 0.025)))
    # --- forward map ---
    fmap = build_forward_map(fq_base, slot=slot, ages=ages, alphas=alphas, out=out / "forward_map", ig_db=ig_db, materials_config=mat_path, use_mock=use_mock, dat_lst=dat_lst, max_xgems_calls=max_xgems_calls, quantities=tuple(q for q in quantities if q != "CH_XRD") + ("CH_XRD",))
    warnings += [str(w) for w in fmap.meta.get("warnings", [])]
    mono = fmap.monotonicity_report()
    for k, flags in mono.items():
        if flags:
            warnings.append(f"monotonicity:{k}:{flags}")
    # --- likelihood & posterior ---
    so = dict(defaults.get("sigma_obs_default", {}))
    sm = dict(defaults.get("sigma_model_initial", {}))
    lik_cfg = defaults.get("likelihood", {})
    q_weights = {k: float(v) for k, v in (lik_cfg.get("quantity_weights") or {}).items()}
    q_offsets = {k: float(v) for k, v in (lik_cfg.get("systematic_offsets") or {}).items()}
    q_scales = {k: float(v) for k, v in (lik_cfg.get("systematic_scales") or {}).items()}
    pts, skipped = build_points(obs_h, ages, sigma_obs_default=so, sigma_model=sm, weights=q_weights)
    warnings += skipped
    if not pts:
        raise ValueError("no usable (grade A/B) observations for the likelihood")
    if q_weights or q_offsets or q_scales:
        warnings.append(f"likelihood policy: weights={q_weights}, systematic_scales={q_scales}, systematic_offsets={q_offsets}")
    lik = Likelihood(fmap, pts, offsets=dict(q_offsets), beta=bundle.bayes.beta_shape, scales=dict(q_scales))
    inv_cfg = defaults.get("inverse", {})
    post = infer(lik, prior_a_max=prior_a, prior_tau=prior_t, prior=prior, ess_min=float(inv_cfg.get("ess_min", 50)), sir_rounds=int(inv_cfg.get("sir_rounds", 1)), grid_n=int(inv_cfg.get("flat_prior_grid", 40)), rng_seed=seed, b_bw_prior=b_bw_prior, b_bw_points=int(inv_cfg.get("b_bw_marginal_points", 5)))
    ages_out = np.asarray(defaults.get("ages_d_default", [1, 3, 7, 28, 90, 180, 365]), float)
    summ = summarise(post, lik, ages_out, prior_a_max=prior_a, prior_tau=prior_t)
    iid = _hid({"scm": scm_m.model_dump(), "mix": mix_m.model_dump(), "obs": obs_h, "mix_uid": mix_uid, "seed": seed})
    # direct DoR (validation only)
    direct = [o for o in obs_h if o["quantity"] == "DoR_SCM" and o.get("value") is not None]
    validation = None
    if direct:
        validation = []
        for o in direct:
            j = int(np.argmin(np.abs(ages_out - o["age_d"])))
            validation.append({"age_d": o["age_d"], "dor_measured": o["value"], "posterior_q50_at_nearest_age": summ["alpha"]["q50"][j], "nearest_age": float(ages_out[j])})
    inv_cfg_status = defaults.get("inverse", {})
    header = {"status": inv_cfg_status.get("status", "validated"), "status_by_slot": (inv_cfg_status.get("status_by_slot") or {}).get(slot), "status_note": inv_cfg_status.get("status_note"), "mix_uid": mix_uid, "scm": scm_m.model_dump(), "mix": mix_m.model_dump(), "slot": slot, "materials_config": mat_path, "observations": obs_h, "prior": prior, "alpha_grid": alphas.tolist(), "use_mock": use_mock, "direct_dor_validation": validation, "warnings": warnings, "provenance": {"dorgems": __version__, **bundle.provenance, "seed": seed}}
    files = write_inference(summ, out, inference_id=iid, header=header)
    files.update(fmap.save(out / "forward_map"))
    # exports
    pred_like = {"id": f"inferred_{iid}", "input": {"ages_d": ages_out.tolist()}, "beta_shape": bundle.bayes.beta_shape, "bayes": {"a_max": {k: v * 100 for k, v in summ["a_max"].items()}, "tau_d": summ["tau_d"]}, "recommended": {"source": "bayes"}, "provenance": header["provenance"]}
    try:
        rm = export_reaction_model(pred_like, out / "reaction_models", mode="logistic_fit", slot=slot, config_id=f"inferred_{iid}", signature_files=[mat_path] if mat_path else None)
        files.update({f"reaction_model_{k}": v["path"] for k, v in rm.items()})
        warnings += [w for r in rm.values() for w in r["warnings"]]
    except Exception as exc:  # noqa: BLE001
        warnings.append(f"reaction model export failed: {exc}")
    csv = out / "inferred_dor.csv"
    with open(csv, "w", encoding="utf-8") as f:
        f.write("scm,age_d,dor\n")
        for age, v in zip(summ["alpha"]["ages_d"], summ["alpha"]["q50"]):
            f.write(f"{slot},{age:g},{v:.5f}\n")
    files["inferred_dor_csv"] = str(csv)
    manifest = {"dorgems": __version__, "inference_id": iid, "files": files, "xgems_calls": fmap.meta.get("xgems_calls"), "use_mock": use_mock, "warnings": warnings}
    (out / "manifest.json").write_text(json.dumps(manifest, indent=2, default=str), encoding="utf-8")
    files["manifest"] = str(out / "manifest.json")
    inference = json.loads(Path(files["inference"]).read_text(encoding="utf-8"))
    inference["files"] = files
    Path(files["inference"]).write_text(json.dumps(inference, indent=2, ensure_ascii=False, default=str), encoding="utf-8")
    return {"ok": True, "inference_id": iid, "slot": slot, "summary": {"a_max": summ["a_max"], "tau_d": summ["tau_d"], "alpha_q50": summ["alpha"]["q50"], "alpha_q05": summ["alpha"]["q05"], "alpha_q95": summ["alpha"]["q95"], "ages_d": summ["alpha"]["ages_d"], "ess": summ["ess"], "posterior_method": summ["posterior_method"], "kl": summ["prior_vs_posterior_kl"], "n_observations_used": summ["n_observations_used"], "xgems_calls": fmap.meta.get("xgems_calls")}, "files": files, "warnings": warnings, "validation": validation}
