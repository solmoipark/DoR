"""One-file input for a new material (templates/new_material_template.yaml):
``scm``, ``mix``, optional ``ages_d``, ``observations``, ``notes``, ``source``.

``load_input_file`` validates against the pydantic schemas and returns the parts;
``validate_report`` additionally grades every observation with the unit/basis rules
so the user sees, before any computation, which observations will count.
"""

from __future__ import annotations

import math
from pathlib import Path
from typing import Any

import yaml

from ..db.units import harmonize
from .schemas import MixSpec, ObservationSpec, SCMSpec, coerce_mix, coerce_observations, coerce_scm

PLACEHOLDER = "?"


def _strip_placeholders(obj: Any, path: str = "", missing: list[str] | None = None) -> Any:
    """Replace '?' placeholders by None and record their paths."""
    if missing is None:
        missing = []
    if isinstance(obj, dict):
        return {k: _strip_placeholders(v, f"{path}.{k}" if path else k, missing) for k, v in obj.items()}
    if isinstance(obj, list):
        return [_strip_placeholders(v, f"{path}[{i}]", missing) for i, v in enumerate(obj)]
    if obj == PLACEHOLDER:
        missing.append(path)
        return None
    return obj


def load_input_file(path: str | Path) -> dict[str, Any]:
    raw = yaml.safe_load(Path(path).read_text(encoding="utf-8")) or {}
    if not isinstance(raw, dict) or "scm" not in raw or "mix" not in raw:
        raise ValueError("input file must contain 'scm' and 'mix' mappings")
    missing: list[str] = []
    data = _strip_placeholders(raw, "", missing)
    scm_raw = dict(data["scm"] or {})
    scm_raw["oxides"] = {k: v for k, v in (scm_raw.get("oxides") or {}).items() if v is not None}
    scm_raw = {k: v for k, v in scm_raw.items() if v is not None}
    mix_raw = {k: v for k, v in dict(data["mix"] or {}).items() if v is not None}
    if isinstance(mix_raw.get("opc_oxides"), dict):
        mix_raw["opc_oxides"] = {k: v for k, v in mix_raw["opc_oxides"].items() if v is not None}
    obs_raw = [o for o in (data.get("observations") or []) if isinstance(o, dict) and o.get("value") is not None]
    scm = coerce_scm(scm_raw)
    mix = coerce_mix(mix_raw)
    observations = coerce_observations([{k: v for k, v in o.items() if v is not None} for o in obs_raw])
    return {"scm": scm, "mix": mix, "ages_d": data.get("ages_d"), "observations": observations, "notes": data.get("notes"), "source": data.get("source"), "placeholders_left": missing, "observations_dropped": len(data.get("observations") or []) - len(obs_raw)}


def validate_report(path: str | Path) -> dict[str, Any]:
    """Never raises for schema problems: returns {ok, errors, warnings, grades, summary}."""
    errors: list[str] = []
    warnings: list[str] = []
    try:
        parts = load_input_file(path)
    except Exception as exc:  # noqa: BLE001
        return {"ok": False, "errors": [f"{type(exc).__name__}: {exc}"], "warnings": [], "grades": [], "summary": {}}
    scm: SCMSpec = parts["scm"]
    mix: MixSpec = parts["mix"]
    if parts["placeholders_left"]:
        warnings.append(f"'?' left unfilled at: {parts['placeholders_left']} (treated as missing)")
    if parts["observations_dropped"]:
        warnings.append(f"{parts['observations_dropped']} observation(s) without a value were dropped")
    ox_sum = sum(v for k, v in scm.oxides.items() if k != "LOI")
    if abs(ox_sum - 100) > 3:
        warnings.append(f"oxide sum {ox_sum:.1f} % (LOI excluded) will be renormalised to 100")
    if scm.amorphous_pct is not None:
        warnings.append("amorphous_pct given: make sure it is a measured value, not 100 − Σ crystalline")
    from ..kinetics.materials_override import slot_for_role

    try:
        slot, w, reactive = slot_for_role(scm.role, scm.oxides)
        warnings += w
    except Exception as exc:  # noqa: BLE001
        errors.append(f"slot mapping failed: {exc}")
        slot, reactive = None, None
    grades = []
    usable = 0
    for o in parts["observations"]:
        o: ObservationSpec
        h = harmonize({"quantity": o.quantity, "value": o.value, "unit": o.unit, "unit_reported": o.unit, "basis_reported": None}, {"scm_total_pct": mix.scm_pct, "w_b": mix.w_b}, scm_pct=mix.scm_pct)
        role = "validation_only" if o.quantity in ("DoR_SCM", "DoR_clinker") else ("primary" if o.quantity in ("bound_water", "chem_shrink", "QXRD_phase") else "secondary")
        grades.append({"age_d": o.age_d, "quantity": o.quantity, "value": o.value, "unit": o.unit, "grade": h.grade, "harmonised": None if h.value is None or (isinstance(h.value, float) and math.isnan(h.value)) else round(h.value, 4), "target_unit": h.unit, "role": role, "assumptions": h.assumptions})
        if h.usable and role != "validation_only":
            usable += 1
    ages_obs = sorted({o.age_d for o in parts["observations"]})
    summary = {"name": scm.name, "role": scm.role, "slot": slot, "reactive_slot": reactive, "scm_pct": mix.scm_pct, "w_b": mix.w_b, "curing_temp_C": mix.curing_temp_C, "n_observations": len(parts["observations"]), "n_usable_for_likelihood": usable, "observation_ages_d": ages_obs, "can_run": {"A_envelope": not errors, "B_compare": usable > 0, "C_infer": usable > 0 and len(ages_obs) >= 3}}
    if parts["observations"] and usable == 0:
        warnings.append("no observation has a usable (grade A/B) unit — write the basis into 'unit', e.g. 'g/100 g binder'")
    return {"ok": not errors, "errors": errors, "warnings": warnings, "grades": grades, "summary": summary}
