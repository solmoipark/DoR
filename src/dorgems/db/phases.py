"""DB phase_name ↔ Cemdata18 raw phase alias table (spec §4.4). Unknown names raise."""

from __future__ import annotations

from typing import Any

import yaml

from ..config import configs_dir


class UnknownPhaseError(KeyError):
    pass


_TABLE: dict[str, Any] | None = None


def load_table() -> dict[str, Any]:
    global _TABLE
    if _TABLE is None:
        _TABLE = yaml.safe_load((configs_dir() / "phase_aliases.yaml").read_text(encoding="utf-8"))
    return _TABLE


def _norm(s: str) -> str:
    return " ".join(s.strip().lower().replace("_", " ").split())


def group_of_db_phase(phase_name: str) -> str | None:
    """Hydrate group key for a DB phase_name, 'clinker:<C3S|C2S|C3A|C4AF>' for clinker
    phases, 'sulfate_initial', 'ignore', or None when unknown."""
    t = load_table()
    n = _norm(phase_name)
    for key, g in t["groups"].items():
        if n in [_norm(x) for x in g["db"]]:
            return key
    for key, g in t["clinker"].items():
        if n in [_norm(x) for x in g["db"]]:
            return f"clinker:{g['bogue']}"
    if n in [_norm(x) for x in t["sulfates_initial"]]:
        return "sulfate_initial"
    if n in [_norm(x) for x in t["ignore"]]:
        return "ignore"
    return None


def raw_names(group: str, *, mock: bool = False) -> list[str]:
    t = load_table()
    if mock:
        # the mock runner only produces a few phases; groups it lacks are simply absent
        return list(t.get("mock", {}).get(group) or [])
    if group not in t["groups"]:
        raise UnknownPhaseError(f"unknown phase group {group!r}")
    return list(t["groups"][group]["raw"])


def compare_mode(group: str) -> str:
    return str(load_table()["groups"][group].get("compare", "ratio"))


def is_confirmed() -> bool:
    return bool(load_table().get("confirmed", False))


def phase_mass_of_group(row: dict[str, Any], group: str, *, mock: bool = False, strict: bool = True) -> float | None:
    """Sum of ``phase_mass__<raw>`` columns of a time-series row for a group (kernel units, kg)."""
    names = raw_names(group, mock=mock)
    if not names:
        return None
    cols = [f"phase_mass__{n}" for n in names]
    present = [c for c in cols if c in row]
    if not present:
        if strict and not mock and is_confirmed():
            raise UnknownPhaseError(f"none of {names} present in the run output for group {group!r}")
        return None
    return float(sum(float(row[c] or 0.0) for c in present))


def confirm_from_raw_names(raw_phase_names: list[str], source: str) -> dict[str, Any]:
    """G2-1 helper: report which table entries exist in a real run's raw phase list.
    Does not edit the YAML (a human confirms and flips ``confirmed: true``)."""
    t = load_table()
    have = set(raw_phase_names)
    report = {"source": source, "groups": {}, "unmatched_raw": sorted(have - {r for g in t["groups"].values() for r in g["raw"]})}
    for key, g in t["groups"].items():
        found = [r for r in g["raw"] if r in have]
        report["groups"][key] = {"expected": g["raw"], "found": found, "ok": bool(found)}
    return report
