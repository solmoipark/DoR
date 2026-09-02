"""Twin mode (spec §8.1, §8.4): rebuild a literature mix as an InverseGems recipe,
run it with the mix's own DoR (pinned) or the model q50, and compare 1:1 with
the observations of the same mix. Plus the OPC-only reference check (§8.6).
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from ..db.features import CEMENT_ROLES, FILLER_ROLES, blended_only, build_dor_table, parse_binder_composition, resolve_scm
from ..db.reader import LiteratureDB
from ..db.units import harmonize
from ..envelope import build_forward_query
from ..gems.forward import run_forward
from ..gems.observables import observables_from_run
from ..kinetics.curves import stretched_exp
from ..kinetics.materials_override import SLOTS, build_materials_config, slot_for_role
from ..kinetics.reaction_model import export_reaction_model, pin_reaction_model
from ..pilot.schemas import SCMSpec
from .compare import OBS_TO_MODEL, aggregate, compare_rows, write_comparison

OPC_REQUIRED = ("CaO", "SiO2", "Al2O3", "Fe2O3", "SO3")
COMPARE_QUANTITIES = ("CH_TGA", "CH_XRD", "bound_water", "chem_shrink")


def _role_to_slot(role: str | None, oxides: dict[str, float] | None) -> tuple[str | None, list[str]]:
    if role in ("cement", "clinker"):
        return "OPC", []
    if role in ("gypsum", "sulfate_source"):
        return "gypsum", []
    if role == "limestone":
        return "limestone", []
    if role is None or role in CEMENT_ROLES or role in FILLER_ROLES:
        return None, [f"role {role!r} has no InverseGems slot"]
    try:
        slot, w, _ = slot_for_role(role, oxides)
        return slot, w
    except Exception as exc:  # noqa: BLE001
        return None, [str(exc)]


def db_mix_to_recipe(mix: dict[str, Any], materials: list[dict[str, Any]], *, ages: list[float], out_dir: Path) -> dict[str, Any]:
    """→ {forward_query, materials_config, slot, scm_material, warnings, excluded_reason}."""
    warnings: list[str] = []
    bc = parse_binder_composition(mix.get("binder_composition_json"), mix.get("notes"))
    by_id = {m["material_id"]: m for m in materials}
    binders: dict[str, float] = {}
    unmapped = 0.0
    scm_slot, scm_mat = None, None
    for mid, pct in bc.items():
        m = by_id.get(mid)
        if m is None:
            unmapped += float(pct or 0)
            warnings.append(f"binder component {mid!r} not in materials table")
            continue
        ox = {k: m[k] for k in ("CaO", "SiO2", "Al2O3", "Fe2O3", "MgO", "SO3", "Na2O", "K2O", "TiO2") if m.get(k) is not None}
        slot, w = _role_to_slot(m.get("role"), ox if ox.get("SiO2") else None)
        warnings += w
        if slot is None:
            unmapped += float(pct or 0)
            continue
        binders[slot] = binders.get(slot, 0.0) + float(pct or 0)
        if slot in SLOTS and (scm_mat is None or float(pct or 0) > float(scm_mat[1] or 0)):
            scm_slot, scm_mat = slot, (m, pct)
    total = sum(binders.values())
    if total <= 0:
        return {"excluded_reason": "no mappable binder components", "warnings": warnings}
    if unmapped > 5.0:
        return {"excluded_reason": f"unmappable components {unmapped:.1f} % > 5 %", "warnings": warnings}
    if mix.get("w_b") is None or mix.get("curing_temp_C") is None:
        return {"excluded_reason": "w_b or curing_temp_C missing", "warnings": warnings}
    binders = {k: v * 100.0 / total for k, v in binders.items()}
    # materials override: SCM oxides from the DB, OPC if complete
    mat_cfg = None
    cement = None
    opc_rows = [m for m in materials if m.get("role") in ("cement", "clinker")]
    if opc_rows and all(opc_rows[0].get(k) is not None for k in OPC_REQUIRED):
        cement = {k: float(opc_rows[0][k]) for k in ("CaO", "SiO2", "Al2O3", "Fe2O3", "MgO", "SO3", "Na2O", "K2O") if opc_rows[0].get(k) is not None}
    else:
        warnings.append("OPC composition incomplete in DB; kernel default OPC used")
    if scm_mat is not None:
        m = scm_mat[0]
        ox = {k: float(m[k]) for k in ("CaO", "SiO2", "Al2O3", "Fe2O3", "MgO", "SO3", "Na2O", "K2O", "TiO2") if m.get(k) is not None}
        if ox.get("SiO2") is not None and sum(ox.values()) > 50:
            spec = SCMSpec(name=str(m.get("name_in_paper") or m["material_id"]), role=_scm_role_literal(m.get("role")), oxides=ox, blaine_m2_kg=m.get("blaine_m2_kg"), amorphous_pct=m.get("amorphous_pct"))
            res = build_materials_config(spec, out_dir, slot=scm_slot, cement=cement)
            mat_cfg = res["path"]
            warnings += res["warnings"]
        else:
            warnings.append("SCM oxides incomplete in DB; slot default composition used")
    elif cement is not None:
        spec = SCMSpec(name="opc_only", role="slag", oxides={"SiO2": 36.49, "Al2O3": 12.26, "CaO": 41.79, "MgO": 7.48, "SO3": 1.98})
        res = build_materials_config(spec, out_dir, slot="slag", cement=cement)
        mat_cfg = res["path"]
    fq = {
        "name": f"twin_{mix['mix_uid']}".replace("/", "_").replace(":", "_"),
        "task": "forward_time_series",
        "recipe": {"binders": {k: round(v, 6) for k, v in binders.items() if v > 0}, "w_b": float(mix["w_b"])},
        "age_grid": {"values": [float(a) for a in ages]},
        "temperature_celsius": float(mix["curing_temp_C"]),
        "outputs": {"phase_masses": "all", "phase_volumes": "all", "phase_volumes_reconstructed": "all", "aqueous_species": "all", "scalars": "all"},
        "plots": [],
        "response_summary": {"enabled": False},
    }
    return {"forward_query": fq, "materials_config": mat_cfg, "slot": scm_slot, "scm_material": scm_mat[0] if scm_mat else None, "scm_pct": binders.get(scm_slot) if scm_slot else 0.0, "warnings": warnings, "excluded_reason": None}


def _scm_role_literal(role: str | None) -> str:
    from ..pilot.schemas import ROLES

    return role if role in ROLES else "other"


def dor_pin_curve(dor_obs: pd.DataFrame, ages: np.ndarray) -> tuple[np.ndarray | None, list[str]]:
    """Measured DoR of the mix → alpha(ages): stretched-exp fit if ≥3 ages, else None."""
    if dor_obs.empty:
        return None, ["no measured DoR for this mix"]
    g = dor_obs.groupby("age_d")["dor_pct"].mean()
    if len(g) < 3:
        return None, [f"only {len(g)} DoR ages; need ≥3 for pin interpolation"]
    from scipy.optimize import least_squares

    t, y = g.index.values.astype(float), g.values.astype(float) / 100.0

    def resid(p: np.ndarray) -> np.ndarray:
        return stretched_exp(t, p[0], np.exp(p[1]), 0.5) - y

    r = least_squares(resid, x0=[min(max(y.max() * 1.2, 0.05), 1.0), np.log(20.0)], bounds=([0.01, np.log(0.1)], [1.0, np.log(5000.0)]))
    a, tau = r.x[0], np.exp(r.x[1])
    return stretched_exp(ages, a, tau, 0.5), [f"DoR pinned from {len(g)} measured ages (a_max={a:.3f}, tau={tau:.1f} d)"]


def twin_compare_mix(
    db: LiteratureDB,
    mix_uid: str,
    *,
    out: Path,
    ig_db: str | Path,
    use_mock: bool = True,
    dat_lst: str | Path | None = None,
    max_xgems_calls: int | None = None,
    quantities: tuple[str, ...] = COMPARE_QUANTITIES,
    dor_source: str = "pin",
    bundle: Any = None,
) -> dict[str, Any]:
    out.mkdir(parents=True, exist_ok=True)
    mix = db.mix(mix_uid)
    if mix is None:
        return {"ok": False, "error": f"mix {mix_uid} not found"}
    mats = db.materials_for_paper(mix["paper_doi"])
    obs = db.observations_for_mix(mix_uid, list(quantities))
    obs = [o for o in obs if o["age_d"] is not None and o["age_d"] > 0]
    if not obs:
        return {"ok": False, "error": "no comparable observations for this mix", "mix_uid": mix_uid}
    ages = sorted({float(o["age_d"]) for o in obs})
    rec = db_mix_to_recipe(mix, mats, ages=ages, out_dir=out)
    if rec.get("excluded_reason"):
        return {"ok": False, "error": rec["excluded_reason"], "warnings": rec["warnings"], "mix_uid": mix_uid}
    warnings = list(rec["warnings"])
    slot = rec["slot"]
    # --- DoR for the SCM slot ---
    rm_path = None
    pin_ages = np.asarray(ages, float)
    if slot is not None:
        table = build_dor_table(db.con)
        dor_mix = table[table["mix_uid"] == mix_uid]
        curve, w = dor_pin_curve(dor_mix, pin_ages)
        warnings += w
        if curve is not None and dor_source == "pin":
            # one config per age is exact for pins; use the fitted curve via native-free logistic fit instead
            pred = {"id": f"twin_{abs(hash(mix_uid)) % 10**8}", "input": {"ages_d": ages}, "beta_shape": 0.5, "bayes": {"a_max": {"q50": float(curve.max() / max(1e-9, (1 - np.exp(-((pin_ages.max() / 20.0) ** 0.5)))) * 100)}, "tau_d": {"q50": 20.0}}, "recommended": {"source": "bayes"}, "provenance": {}}
            # simpler and exact: re-fit params for export
            from scipy.optimize import least_squares

            def resid(p: np.ndarray) -> np.ndarray:
                return stretched_exp(pin_ages, p[0], np.exp(p[1]), 0.5) - curve

            r = least_squares(resid, x0=[float(curve.max()), np.log(20.0)], bounds=([0.01, np.log(0.1)], [1.0, np.log(5000.0)]))
            pred["bayes"] = {"a_max": {"q50": float(r.x[0]) * 100}, "tau_d": {"q50": float(np.exp(r.x[1]))}}
            rm = export_reaction_model(pred, out / "rm", mode="logistic_fit", slot=slot, quantiles=(0.5,), config_id=pred["id"], signature_files=[rec["materials_config"]] if rec["materials_config"] else None)
            rm_path = rm["q50"]["path"]
            dor_used = "pin"
        else:
            from ..predict import predict

            m = rec["scm_material"]
            ox = {k: float(m[k]) for k in ("CaO", "SiO2", "Al2O3", "Fe2O3", "MgO", "SO3", "Na2O", "K2O", "TiO2") if m and m.get(k) is not None}
            if not ox.get("SiO2"):
                ox = {"SiO2": 40.0, "CaO": 40.0, "Al2O3": 12.0}
                warnings.append("SCM oxides unknown; generic composition used for the DoR prior")
            spec = SCMSpec(name=str((m or {}).get("name_in_paper") or "scm"), role=_scm_role_literal((m or {}).get("role")), oxides=ox)
            p = predict(spec, {"scm_pct": float(rec["scm_pct"] or 0), "w_b": float(mix["w_b"]), "curing_temp_C": float(mix["curing_temp_C"])}, ages, bundle=bundle, db_path=db.path)
            rm = export_reaction_model(p, out / "rm", mode="logistic_fit", slot=slot, quantiles=(0.5,), config_id=p["id"], signature_files=[rec["materials_config"]] if rec["materials_config"] else None)
            rm_path = rm["q50"]["path"]
            dor_used = "model_q50"
    else:
        dor_used = "none(opc_only)"
    res = run_forward(rec["forward_query"], out=out / "run", db=ig_db, reaction_model_config=rm_path, materials_config=rec["materials_config"], slot=slot, use_mock=use_mock, dat_lst=dat_lst, max_xgems_calls=max_xgems_calls, capture_species=True)
    warnings += res.warnings
    if not res.ok:
        return {"ok": False, "error": res.error or "forward run failed", "warnings": warnings, "mix_uid": mix_uid, "run_dir": str(res.run_dir)}
    model = observables_from_run(res.forward_dir)
    pairs = []
    for o in obs:
        h = harmonize(o, mix, scm_pct=rec.get("scm_pct"))
        col = OBS_TO_MODEL.get(o["quantity"])
        row = model[model["age_d"] == float(o["age_d"])]
        mv = float(row.iloc[0][col]) if (col and not row.empty and row.iloc[0][col] is not None and not pd.isna(row.iloc[0][col])) else None
        pairs.append({"obs_uid": o["obs_uid"], "paper_doi": o["paper_doi"], "mix_uid": mix_uid, "quantity": o["quantity"], "phase_name": o.get("phase_name"), "age_d": float(o["age_d"]), "method": o.get("method"), "grade": h.grade, "assumptions": "; ".join(h.assumptions), "obs_value": h.value, "model_value": mv, "uncertainty": o.get("uncertainty"), "source_locator": o.get("source_locator"), "fig_only": o.get("fig_only"), "extraction_confidence": o.get("extraction_confidence")})
    df = compare_rows(pairs)
    agg = aggregate(df)
    files = write_comparison(df, agg, out, header={"mode": "twin", "target": mix_uid, "dor_source": dor_used, "use_mock": use_mock, "slot": slot, "materials_injection": res.materials_injection})
    return {"ok": True, "mix_uid": mix_uid, "paper_doi": mix["paper_doi"], "slot": slot, "dor_source": dor_used, "n_obs": len(pairs), "aggregate": agg, "files": files, "run_dir": str(res.run_dir), "warnings": warnings}


def opc_reference_candidates(db: LiteratureDB, *, age_days: float = 28, w_b_range: tuple[float, float] = (0.4, 0.5), temp_C: float = 20.0, temp_tol: float = 3.0, quantity: str = "CH_TGA") -> pd.DataFrame:
    rows = db.opc_only_reference(quantity, age_days, 0.15, w_b_range=w_b_range)
    mats_all = db.materials_by_paper()
    keep = []
    for r in rows:
        bc = parse_binder_composition(r.get("binder_composition_json"), None)
        pm = mats_all.get(r["paper_doi"], {})
        if any(mid in pm and pm[mid]["role"] not in CEMENT_ROLES and pm[mid]["role"] not in FILLER_ROLES for mid in bc):
            continue  # SCM present in the binder JSON (spec §1.1: 393 is an upper bound)
        T = r.get("curing_temp_C")
        if T is not None and abs(float(T) - temp_C) > temp_tol:
            continue
        h = harmonize(r, r)
        keep.append({**r, "grade": h.grade, "obs_harmonised": h.value, "assumptions": "; ".join(h.assumptions)})
    return pd.DataFrame(keep)


def opc_reference_check(db: LiteratureDB, *, out: Path, ig_db: str | Path, age_days: float = 28, w_b_range: tuple[float, float] = (0.4, 0.5), use_mock: bool = True, dat_lst: str | Path | None = None, max_xgems_calls: int | None = None, max_mixes: int | None = None, grades: tuple[str, ...] = ("A",)) -> dict[str, Any]:
    out.mkdir(parents=True, exist_ok=True)
    cands = opc_reference_candidates(db, age_days=age_days, w_b_range=w_b_range)
    cands = cands[cands["grade"].isin(grades)] if not cands.empty else cands
    n_all = int(len(cands))
    if max_mixes:
        cands = cands.drop_duplicates("mix_uid").head(int(max_mixes))
    pairs = []
    warnings: list[str] = []
    for _, r in cands.iterrows():
        mix = db.mix(r["mix_uid"])
        mats = db.materials_for_paper(r["paper_doi"])
        rec = db_mix_to_recipe(mix, mats, ages=[float(r["age_d"])], out_dir=out / "materials")
        if rec.get("excluded_reason"):
            warnings.append(f"{r['mix_uid']}: {rec['excluded_reason']}")
            continue
        run_out = out / "runs" / str(r["mix_uid"]).replace("/", "_").replace(":", "_")
        res = run_forward(rec["forward_query"], out=run_out, db=ig_db, reaction_model_config=None, materials_config=rec["materials_config"], slot=None, use_mock=use_mock, dat_lst=dat_lst, max_xgems_calls=max_xgems_calls, capture_species=True)
        if not res.ok:
            warnings.append(f"{r['mix_uid']}: run failed {res.error}")
            continue
        model = observables_from_run(res.forward_dir)
        mv = model.iloc[0]["CH_g"] if not model.empty else None
        pairs.append({"obs_uid": r["obs_uid"], "paper_doi": r["paper_doi"], "mix_uid": r["mix_uid"], "quantity": "CH_TGA", "age_d": float(r["age_d"]), "method": r.get("method"), "grade": r["grade"], "assumptions": r["assumptions"], "obs_value": r["obs_harmonised"], "model_value": None if mv is None or pd.isna(mv) else float(mv), "uncertainty": r.get("uncertainty"), "source_locator": r.get("source_locator"), "fig_only": r.get("fig_only"), "extraction_confidence": r.get("extraction_confidence")})
    df = compare_rows(pairs)
    agg = aggregate(df)
    ch = agg["by_quantity"].get("CH_TGA", {})
    gate = {"n_papers": ch.get("n_papers", 0), "median_r": ch.get("median_r"), "pass": (ch.get("n_papers", 0) >= 30 and ch.get("median_r") is not None and abs(ch["median_r"]) <= 4.0)}
    if use_mock:
        gate["note"] = "mock runner: numbers are not physical; pipeline check only"
    sigma_model_est = {"CH_TGA": float(df.loc[df["usable"], "r"].std()) if df["usable"].sum() >= 3 else None}
    files = write_comparison(df, agg, out, header={"mode": "opc_reference", "target": f"OPC-only, {age_days} d, w/b {w_b_range}", "use_mock": use_mock, "n_candidates": n_all, "gate_G2_3": gate, "sigma_model_estimate": sigma_model_est})
    return {"ok": True, "n_candidates": n_all, "n_run": len(pairs), "aggregate": agg, "gate_G2_3": gate, "sigma_model_estimate": sigma_model_est, "files": files, "warnings": warnings}


__all__ = ["db_mix_to_recipe", "twin_compare_mix", "opc_reference_check", "opc_reference_candidates", "blended_only", "resolve_scm"]
