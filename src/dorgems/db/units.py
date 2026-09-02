"""Unit / basis harmonisation (spec §4.3). Rules live in configs/unit_basis.yaml."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

import yaml

from ..config import configs_dir

_RULES: dict[str, Any] | None = None


def load_rules() -> dict[str, Any]:
    global _RULES
    if _RULES is None:
        p = configs_dir() / "unit_basis.yaml"
        _RULES = yaml.safe_load(p.read_text(encoding="utf-8"))
    return _RULES


@dataclass
class HarmonizedObs:
    value: float | None  # in the target unit
    unit: str
    grade: str  # A B C D X
    factor: float | None
    assumptions: list[str] = field(default_factory=list)
    source_unit: str | None = None
    source_basis: str | None = None

    @property
    def usable(self) -> bool:
        return self.grade in ("A", "B") and self.value is not None


def _norm(s: str | None) -> str:
    return (s or "").strip().lower().replace("  ", " ")


def _match(rules: list[dict[str, Any]], unit: str | None, basis: str | None) -> dict[str, Any] | None:
    u, b = _norm(unit), _norm(basis)
    for rule in rules:
        keys = [_norm(k) if k is not None else "" for k in rule["match"]]
        if u in keys or (b and b in keys):
            return rule
    return None


def _eval_factor(expr: str | None, ctx: dict[str, float]) -> float | None:
    if expr is None:
        return None
    return float(eval(expr, {"__builtins__": {}}, ctx))  # noqa: S307 — expressions come from our own YAML


def harmonize(obs: dict[str, Any], mix: dict[str, Any] | None = None, *, scm_pct: float | None = None) -> HarmonizedObs:
    """``obs`` needs quantity, value_norm (or value), unit_norm (or unit), basis_reported;
    ``mix`` provides scm_total_pct and w_b for the basis factors."""
    rules = load_rules()
    q = obs.get("quantity")
    value = obs.get("value_norm", obs.get("value"))
    unit = obs.get("unit_norm", obs.get("unit"))
    basis = obs.get("basis_reported")
    mix = mix or {}
    scm_total = mix.get("scm_total_pct")
    w_b = mix.get("w_b")
    ctx = {
        "f_opc": 1.0 - float(scm_total) / 100.0 if scm_total is not None else 1.0,
        "w_b": float(w_b) if w_b is not None else 0.0,
        "scm_pct": (float(scm_pct) / 100.0) if scm_pct is not None else (float(scm_total) / 100.0 if scm_total else 1.0),
        "gypsum_frac": float(rules.get("gypsum_frac_default", 0.05)),
    }
    assumptions: list[str] = []
    if value is None:
        return HarmonizedObs(None, "", "X", None, ["no value"], unit, basis)
    if q in rules["mass_quantities"]:
        target = "g/100 g binder"
        rule = _match(rules["rules"], unit, basis)
        if rule is None:
            rule = _match(rules["rules"], None, "mass_percent_unspecified")
            assumptions.append(f"unit {unit!r}/basis {basis!r} not in unit_basis.yaml → grade D")
        if rule["grade"] == "A" and rule["factor"] == "f_opc" and scm_total is None:
            assumptions.append("scm_total_pct missing; f_OPC assumed 1")
        if rule["factor"] == "1 + w_b" and w_b is None:
            return HarmonizedObs(None, target, "X", None, ["paste basis but w/b unknown"], unit, basis)
    elif q == "QXRD_phase":
        target = "wt% (as reported)"
        rule = _match(rules["qxrd"], unit, basis) or rules["qxrd"][0]
    elif q == "chem_shrink":
        target = "mL/g binder"
        rule = _match(rules["chem_shrink"], unit, basis) or rules["chem_shrink"][-1]
    elif q in ("DoR_SCM", "DoR_clinker"):
        target = "fraction"
        rule = _match(rules["dor"], unit, basis) or rules["dor"][0]
    else:
        return HarmonizedObs(None, "", "X", None, [f"quantity {q!r} not harmonised in v1"], unit, basis)
    factor = _eval_factor(rule.get("factor"), ctx)
    if rule.get("assumption"):
        assumptions.append(str(rule["assumption"]))
    val = None if factor is None else float(value) * factor
    return HarmonizedObs(val, target, str(rule["grade"]), factor, assumptions, unit, basis)
