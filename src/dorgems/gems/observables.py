"""xGEMS run output → DB physical quantities (spec §8.3).

Basis: 100 g binder. Kernel masses are kg → ×1000 = g/100 g binder.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from ..db import phases as P

KG_TO_G_PER_100G = 1000.0
M_H2O = 18.015


def _is_mock(row: dict[str, Any]) -> bool:
    return any(c.startswith("phase_mass__Mock") for c in row)


def mass_factor(row: dict[str, Any]) -> float:
    """Factor to g/100 g binder. Real xGEMS reports kg (spec §8.3) → ×1000; the mock
    runner reports grams. Detect from ``scalar__system_mass`` (≈ 100 g binder + water:
    > 10 means grams) and fall back to the mock/real distinction. To be re-checked on
    the first real run (P-IG-4 unit audit)."""
    sm = row.get("scalar__system_mass")
    try:
        if sm is not None and float(sm) > 0:
            return 1.0 if float(sm) > 10.0 else KG_TO_G_PER_100G
    except (TypeError, ValueError):
        pass
    return 1.0 if _is_mock(row) else KG_TO_G_PER_100G


def portlandite_g(row: dict[str, Any]) -> float | None:
    v = P.phase_mass_of_group(row, "portlandite", mock=_is_mock(row), strict=False)
    return None if v is None else v * mass_factor(row)


def hydrate_group_g(row: dict[str, Any], group: str) -> float | None:
    if P.compare_mode(group) == "forbidden":
        return None
    v = P.phase_mass_of_group(row, group, mock=_is_mock(row), strict=False)
    return None if v is None else v * mass_factor(row)


def bound_water_g(capture: dict[str, Any] | None, row: dict[str, Any], *, water_in_g: float | None) -> dict[str, Any]:
    """W_bound = W_in − W_aq, W_aq from the aqueous phase species (capture) or, without
    capture, the aqueous phase mass minus dissolved ions (approximation).
    Returns {value_g, method, warnings}; value_g is None when nothing is available."""
    warnings: list[str] = []
    if water_in_g is None:
        return {"value_g": None, "method": None, "warnings": ["initial water unknown"]}
    w_aq = None
    method = None
    w_solids = None
    if capture:
        psm = capture.get("phase_species_moles") or {}
        aq = psm.get("aq_gen")
        if isinstance(aq, dict):
            h2o_mol = sum(float(mol) for sp, mol in aq.items() if sp.startswith("H2O"))
            if h2o_mol > 0:
                w_aq = h2o_mol * M_H2O  # mol × g/mol = g (100 g binder basis)
                method = "capture:aq_gen_H2O_moles"
        # cross-check: H in solid phases as H2O equivalents (phase_elements_amounts, mol)
        pea = capture.get("phase_elements_amounts") or {}
        if pea:
            h_solid = sum(float((el or {}).get("H", 0.0)) for ph, el in pea.items() if ph not in ("aq_gen", "gas_gen"))
            w_solids = h_solid / 2.0 * M_H2O
        if w_aq is None:
            warnings.append("capture present but aq_gen species moles missing")
    if w_aq is None:
        aq_mass = None
        for c, v in row.items():
            if c.startswith("phase_mass__") and c[len("phase_mass__"):] in ("aq_gen", "aq", "Aqueous"):
                aq_mass = float(v) * mass_factor(row)
        if aq_mass is not None:
            w_aq = aq_mass
            method = "aqueous_phase_mass_approx"
            warnings.append("W_aq approximated by the aqueous phase mass (dissolved ions included)")
    if w_aq is None:
        return {"value_g": None, "method": None, "warnings": warnings + ["no aqueous water information in the run output"], "w_solids_g": w_solids}
    w_bound = float(water_in_g) - float(w_aq)
    if w_solids is not None and abs(w_solids - w_bound) > max(1.0, 0.1 * abs(w_bound)):
        warnings.append(f"bound-water definitions disagree: W_in-W_aq={w_bound:.2f} g vs H-in-solids={w_solids:.2f} g")
    return {"value_g": w_bound, "method": method, "warnings": warnings, "w_solids_g": w_solids}


def chem_shrink_ml_per_g(porosity: dict[str, Any] | None, *, volume_unit: str = "auto") -> dict[str, Any]:
    """CS = V_initial − (V_solid_final + V_aq) per 100 g binder ÷ 100 → mL/g binder.

    Units (verified on the TINN_v4 real run, 2026-09-02, P-IG-4): ``initial_volume_cm3`` and
    ``solid_final_volume_cm3`` are cm³, but ``excluded_non_solid_phase_volumes_raw`` holds the
    raw xGEMS phase volumes in **m³** (e.g. aq_gen 2.3e-5). ``auto`` converts the raw aqueous
    volume to cm³ when it is < 1e-2 (a 100 g-binder paste can never have < 0.01 cm³ of water)."""
    if not porosity:
        return {"value": None, "warnings": ["porosity.json missing"]}
    v0 = porosity.get("initial_volume_cm3")
    vs = porosity.get("solid_final_volume_cm3")
    excl = porosity.get("excluded_non_solid_phase_volumes_raw") or {}
    v_aq = None
    for k, v in excl.items():
        if str(k).startswith("aq"):
            v_aq = float(v)
    if v0 is None or vs is None or v_aq is None:
        return {"value": None, "warnings": ["porosity.json lacks initial/solid/aqueous volumes"]}
    warnings: list[str] = []
    if volume_unit == "auto":
        if v_aq < 1e-2:
            v_aq *= 1.0e6
            warnings.append("raw aqueous volume interpreted as m3 and converted to cm3 (P-IG-4 unit mismatch in porosity.json)")
    elif volume_unit == "m3":
        v_aq *= 1.0e6
    cs = (float(v0) - (float(vs) + v_aq)) / 100.0
    return {"value": cs, "warnings": warnings, "components_cm3": {"initial": float(v0), "solid_final": float(vs), "aqueous": v_aq}}


def unreacted_clinker_phase_g(run_dir: Path, bogue_phase: str) -> float | None:
    """unreacted_masses_g[OPC] × Bogue fraction × (1 − α_phase) from input_reaction_degrees.json."""
    rd_path = run_dir / "input_reaction_degrees.json"
    if not rd_path.is_file():
        return None
    rd = json.loads(rd_path.read_text(encoding="utf-8"))
    opc_alpha = (rd.get("opc") or {}).get(bogue_phase)
    pct = (rd.get("opc_phase_mass_percent") or {}).get(bogue_phase)
    recipe = run_dir / "input_recipe.json"
    opc_mass = None
    if recipe.is_file():
        r = json.loads(recipe.read_text(encoding="utf-8"))
        opc_mass = (r.get("binders") or r.get("masses_g") or {}).get("OPC")
    if opc_alpha is None or pct is None or opc_mass is None:
        return None
    return float(opc_mass) * float(pct) / 100.0 * (1.0 - float(opc_alpha))


def observables_for_row(row: dict[str, Any], *, capture: dict[str, Any] | None = None, porosity: dict[str, Any] | None = None, water_in_g: float | None = None) -> dict[str, Any]:
    out: dict[str, Any] = {"age_d": float(row.get("age_days", np.nan))}
    out["CH_g"] = portlandite_g(row)
    for g in ("ettringite", "monocarbonate", "hemicarbonate", "monosulfate", "calcite", "hydrotalcite", "straetlingite"):
        out[f"{g}_g"] = hydrate_group_g(row, g)
    bw = bound_water_g(capture, row, water_in_g=water_in_g)
    out["bound_water_g"] = bw["value_g"]
    out["bound_water_method"] = bw["method"]
    out["bound_water_solids_g"] = bw.get("w_solids_g")
    cs = chem_shrink_ml_per_g(porosity)
    out["chem_shrink_ml_g"] = cs["value"]
    out["porosity"] = row.get("scalar__porosity", row.get("porosity"))
    out["pH"] = row.get("scalar__pH")
    out["warnings"] = bw["warnings"] + cs["warnings"]
    return out


def observables_from_run(forward_dir: Path) -> pd.DataFrame:
    ts = pd.read_csv(forward_dir / "time_series.csv")
    caps = {}
    cp = forward_dir / "dorgems_captures.json"
    if cp.is_file():
        for c in json.loads(cp.read_text(encoding="utf-8")):
            caps[float(c["age_d"])] = c
    rows = []
    for _, r in ts.iterrows():
        row = r.to_dict()
        age = float(row.get("age_days"))
        cap = caps.get(age)
        porosity = None
        water_in = None
        if cap:
            # porosity.json lives in recipe_runs/<recipe_id>/ (cached_forward.py:642)
            for key in ("recipe_dir", "chemistry_dir"):
                if cap.get(key):
                    pj = sorted(Path(cap[key]).rglob("porosity.json"))
                    if pj:
                        porosity = json.loads(pj[-1].read_text(encoding="utf-8"))
                        break
        if row.get("xgems_water_g") is not None:
            water_in = float(row["xgems_water_g"])
        rows.append(observables_for_row(row, capture=(cap or {}).get("capture"), porosity=porosity, water_in_g=water_in))
    return pd.DataFrame(rows)
