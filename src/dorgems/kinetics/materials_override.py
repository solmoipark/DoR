"""Map a new SCM onto an InverseGems binder slot and write a materials override YAML (spec §6.2).

InverseGems fixes ``SCM_NAMES = {slag, fly_ash, metakaolin, silica_fume}``; a new
SCM can only *replace the oxide composition of one slot* and be reachable under
its own name through ``aliases:``.
"""

from __future__ import annotations

import hashlib
import json
import re
from pathlib import Path
from typing import Any

import numpy as np
import yaml

from ..config import configs_dir, inverse_gems_root

SLOTS = ("slag", "fly_ash", "metakaolin", "silica_fume")
OXIDE_KEYS = ("SiO2", "Al2O3", "Fe2O3", "CaO", "MgO", "SO3", "Na2O", "K2O")


def load_slot_rules() -> dict[str, Any]:
    p = configs_dir() / "slots.yaml"
    return yaml.safe_load(p.read_text(encoding="utf-8")) if p.is_file() else {"slots": {}}


def default_materials_path() -> Path:
    root = inverse_gems_root(required=True)
    return root / "configs" / "materials.yaml"


def load_default_materials(path: str | Path | None = None) -> dict[str, Any]:
    p = Path(path) if path else default_materials_path()
    return yaml.safe_load(p.read_text(encoding="utf-8"))


def _norm_oxides(oxides: dict[str, float]) -> tuple[dict[str, float], list[str]]:
    warnings: list[str] = []
    ox = {k: float(oxides.get(k, 0.0) or 0.0) for k in OXIDE_KEYS}
    extra = {k: v for k, v in oxides.items() if k not in OXIDE_KEYS and k != "LOI"}
    if extra:
        warnings.append(f"oxides not representable in InverseGems dropped: {sorted(extra)}")
    total = sum(ox.values())
    if total <= 0:
        raise ValueError("oxide sum is zero")
    if abs(total - 100.0) > 3.0:
        warnings.append(f"oxide sum {total:.1f} (LOI excluded) renormalised to 100")
        ox = {k: v * 100.0 / total for k, v in ox.items()}
    else:
        # keep values, but make the sum exactly 100 for a clean element vector
        ox = {k: v * 100.0 / total for k, v in ox.items()}
    return {k: round(v, 4) for k, v in ox.items()}, warnings


def chemical_distance_to_slots(oxides: dict[str, float], materials: dict[str, Any]) -> dict[str, float]:
    x = np.array([float(oxides.get(k, 0.0) or 0.0) for k in OXIDE_KEYS])
    out = {}
    for s in SLOTS:
        ref = materials[s]["oxide_mass_percent"]
        r = np.array([float(ref.get(k, 0.0)) for k in OXIDE_KEYS])
        out[s] = float(np.sqrt(np.sum((x - r) ** 2)))
    return out


def slot_for_role(role: str, oxides: dict[str, float] | None = None, materials: dict[str, Any] | None = None) -> tuple[str, list[str], bool]:
    """(slot, warnings, reactive). 'other' picks the chemically nearest slot."""
    rules = load_slot_rules().get("slots", {})
    rule = rules.get(role)
    warnings: list[str] = []
    if rule is None:
        warnings.append(f"role {role!r} has no slot rule; treated as 'other'")
        rule = rules.get("other", {"slot": "nearest"})
    slot = rule["slot"]
    reactive = bool(rule.get("reactive", True))
    if rule.get("warn"):
        warnings.append(str(rule["warn"]))
    if slot == "nearest":
        if not oxides:
            raise ValueError("role 'other' needs oxides to choose the nearest slot")
        mats = materials or load_default_materials()
        d = chemical_distance_to_slots(oxides, mats)
        slot = min(d, key=d.get)
        warnings.append(f"nearest slot by oxide distance: {slot} (distances {json.dumps({k: round(v, 1) for k, v in d.items()})})")
    return slot, warnings, reactive


def normalise_alias(name: str) -> str:
    s = re.sub(r"[^A-Za-z0-9]+", "_", name.strip()).strip("_")
    return s or "user_scm"


def build_materials_config(
    scm: Any,
    out_dir: str | Path,
    *,
    slot: str | None = None,
    alias: str | None = None,
    cement: dict[str, float] | None = None,
    base: str | Path | None = None,
) -> dict[str, Any]:
    """Write ``materials.dorgems_<hash>.yaml`` with the slot's oxides (and density)
    replaced by the SCM's, plus an alias for the user's name. Returns
    {path, slot, alias, oxides, warnings, base_path, hash}."""
    from ..db.features import _get

    role = _get(scm, "role")
    ox_in = dict(_get(scm, "oxides", {}) or {})
    mats = load_default_materials(base)
    if slot is None:
        slot, warnings, _ = slot_for_role(role, ox_in, mats)
    else:
        warnings = []
        if slot not in SLOTS:
            raise ValueError(f"slot must be one of {SLOTS}")
    ox, w2 = _norm_oxides(ox_in)
    warnings += w2
    entry = dict(mats[slot])
    entry["oxide_mass_percent"] = ox
    dens = _get(scm, "density_kg_m3")
    if dens:
        entry["density_g_cm3"] = round(float(dens) / 1000.0, 4)
    alias = alias or normalise_alias(str(_get(scm, "name", slot)))
    aliases = list(entry.get("aliases", []))
    if alias not in aliases and alias != slot:
        aliases.append(alias)
    entry["aliases"] = aliases
    mats[slot] = entry
    if cement:
        cox, w3 = _norm_oxides(cement)
        warnings += [f"OPC: {w}" for w in w3]
        opc = dict(mats["OPC"])
        opc["oxide_mass_percent"] = cox
        mats["OPC"] = opc
    payload = yaml.safe_dump(mats, sort_keys=False)
    h = hashlib.sha256(payload.encode("utf-8")).hexdigest()[:12]
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    path = out / f"materials.dorgems_{h}.yaml"
    path.write_text(payload, encoding="utf-8")
    return {"path": str(path), "slot": slot, "alias": alias, "oxides": ox, "warnings": warnings, "base_path": str(Path(base) if base else default_materials_path()), "hash": h, "opc_overridden": bool(cement)}
