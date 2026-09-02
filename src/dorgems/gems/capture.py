"""Species-level capture for bound-water accounting (spec §6.3 capture path, §8.3).

The kernel's raw output directory does not store ``phase_species_moles``
(database.py:452-457). ``runner_factory`` is only accepted by
``run_forward_cached`` (cached_forward.py:222), so this path calls it per age
with a ``CapturingRunner`` that wraps the real/mock runner and copies what
``capture_raw_state()`` returns into ``<chemistry_dir>/dorgems_capture.json``.
The per-age results are then assembled into a ``time_series.csv`` with the
same column conventions as the standard path. P-IG-2 makes this obsolete.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from .. import __version__

_LAST_CAPTURES: list[dict[str, Any]] = []


class CapturingRunner:
    """Proxy that forwards every call to the wrapped xGEMS (or mock) runner."""

    def __init__(self, inner: Any):
        self._inner = inner
        self.captured: dict[str, Any] | None = None

    def __getattr__(self, name: str) -> Any:
        return getattr(self._inner, name)

    def capture_raw_state(self) -> dict[str, Any]:
        raw = self._inner.capture_raw_state()
        keep = {}
        for k in ("phase_species_moles", "species_in_phase", "species_molar_masses", "phase_masses", "phase_volumes", "aqueous_species", "scalars"):
            if k in raw:
                keep[k] = raw[k]
        # element-level amounts per phase (real xgems exposes these as properties; the
        # kernel's capture does not record them) — used for the bound-water cross-check
        gems = getattr(self._inner, "gems", None)
        for attr in ("phase_elements_amounts", "aq_elements_amounts", "element_molar_masses"):
            try:
                v = getattr(gems, attr, None) if gems is not None else None
                if callable(v):
                    v = v()
                if isinstance(v, dict):
                    keep[attr] = {str(a): ({str(b): float(c) for b, c in x.items()} if isinstance(x, dict) else float(x)) for a, x in v.items()}
            except Exception:  # noqa: BLE001
                pass
        # generic fallbacks for real xGEMS objects exposing these names differently
        for k in list(raw.keys()):
            if "species" in k.lower() and "mol" in k.lower() and k not in keep:
                keep[k] = raw[k]
        self.captured = keep
        _LAST_CAPTURES.append(keep)
        return raw


def make_runner_factory(use_mock: bool):
    from inverse_gems.xgems_runner import MockXGEMSRunner, XGEMSRunner

    def factory(dat_lst_path: Any = None, temperature_celsius: float = 20.0, **kw: Any) -> CapturingRunner:
        if use_mock:
            inner = MockXGEMSRunner(dat_lst_path=dat_lst_path, temperature_celsius=temperature_celsius)
        else:
            # mirror cached_forward.py:400-405 (the kernel's own construction)
            inner = XGEMSRunner(dat_lst_path, temperature_celsius=temperature_celsius, gems_class_path=kw.get("gems_class_path", "xgems:ChemicalEngineDicts"), input_mode=kw.get("input_mode", "formula"))
        return CapturingRunner(inner)

    return factory


def recipe_text(binders: dict[str, float], w_b: float, age_d: float) -> str:
    parts = [f"{name} {float(v):g}" for name, v in binders.items() if float(v) > 0]
    return ", ".join(parts) + f", w/b {float(w_b):g}, age {float(age_d):g}"


def run_forward_capturing(
    forward_query: dict[str, Any],
    *,
    out: Path,
    db: str | Path,
    reaction_model_config: str | Path | None,
    materials_config: str | Path | None,
    slot: str | None,
    use_mock: bool,
    dat_lst: str | Path | None,
    max_xgems_calls: int | None,
    reaction_model_id: str | None = None,
):
    from inverse_gems.cached_forward import run_forward_cached
    from inverse_gems.call_budget import XGEMSCallBudget

    from .forward import ForwardRunResult, materials_config_override, self_check

    recipe = forward_query["recipe"]
    binders = dict(recipe.get("binders") or {})
    w_b = float(recipe["w_b"])
    ages = list(map(float, forward_query["age_grid"]["values"]))
    T = float(forward_query.get("temperature_celsius", 20.0))
    budget = XGEMSCallBudget(max_xgems_calls) if max_xgems_calls else None
    forward_dir = out / "forward"
    forward_dir.mkdir(parents=True, exist_ok=True)
    rows: list[dict[str, Any]] = []
    captures: list[dict[str, Any]] = []
    warnings: list[str] = []
    injection = None
    # recipe ids must be unique within the kernel DB (recipe_runs.recipe_id is UNIQUE):
    # hash the run location + configs and add a uuid tail like the kernel does.
    import hashlib
    import uuid

    run_tag = hashlib.sha256(f"{out.resolve()}|{reaction_model_config}|{materials_config}".encode("utf-8")).hexdigest()[:10] + "_" + uuid.uuid4().hex[:6]
    try:
        with materials_config_override(materials_config) as how:
            injection = how
            for i, age in enumerate(ages, 1):
                _LAST_CAPTURES.clear()
                res = run_forward_cached(
                    recipe_text=recipe_text(binders, w_b, age),
                    db=db,
                    dat_lst=dat_lst,
                    use_mock=use_mock,
                    temperature_celsius=T,
                    xgems_call_budget=budget,
                    # a cache hit would skip the runner (no species capture): force the solve,
                    # it costs ~5 ms per equilibrate on the real kernel
                    force_rerun_xgems=True,
                    runner_factory=make_runner_factory(use_mock),
                    recipe_id=f"dorgems_capture_{run_tag}_age_{i:04d}",
                    recipe_metadata={"dorgems_capture": True, "template_name": forward_query.get("name", "dorgems")},
                    reaction_model_id=reaction_model_id,
                    reaction_model_config=str(reaction_model_config) if reaction_model_config else None,
                )
                chem_dir = Path(res["chemistry_dir"])
                cap = _LAST_CAPTURES[-1] if _LAST_CAPTURES else None
                if cap is None:
                    warnings.append(f"age {age}: cache hit or runner not invoked — no species capture (chemistry_status={res.get('chemistry_status')})")
                    existing = chem_dir / "dorgems_capture.json"
                    cap = json.loads(existing.read_text(encoding="utf-8")) if existing.is_file() else {}
                else:
                    (chem_dir / "dorgems_capture.json").write_text(json.dumps(cap, indent=1, default=_jsonable), encoding="utf-8")
                captures.append({"age_d": age, "chemistry_dir": str(chem_dir), "recipe_dir": res.get("recipe_dir"), "capture": cap})
                rows.append(_row(age, res, chem_dir))
    except Exception as exc:  # noqa: BLE001
        return ForwardRunResult(ok=False, run_dir=out, forward_dir=forward_dir, time_series=None, status="error", materials_injection=injection, self_check={}, warnings=warnings, error=f"{type(exc).__name__}: {exc}")
    ts = pd.DataFrame(rows)
    ts.to_csv(forward_dir / "time_series.csv", index=False)
    (forward_dir / "dorgems_captures.json").write_text(json.dumps(captures, indent=1, default=_jsonable), encoding="utf-8")
    db_root = Path(db)  # InverseGemsDatabase(db) is a *directory* (chemistry_runs/, recipe_runs/, …)
    check = self_check(forward_dir, reaction_model_config=reaction_model_config, slot=slot, ages=ages, materials_config=materials_config, db_dir=db_root)
    ok = all(r.get("chemistry_status") in ("complete", "ok", "success", "mock_success", "cached") or r.get("solver_status") in ("mock_success", "success", "ok") for r in rows) and bool(rows)
    if reaction_model_config and check.get("alpha_ok") is False:
        ok = False
        warnings.append("self-check failed: kernel reaction degrees do not match the exported reaction model")
    if materials_config and check.get("materials_ok") is False:
        ok = False
        warnings.append("self-check failed: input_materials_used.json does not carry the override oxides")
    if injection == "monkeypatch":
        warnings.append("materials_config injected by monkeypatching inverse_gems.materials.config_path")
    (out / "dorgems_forward_manifest.json").write_text(json.dumps({"dorgems": __version__, "path": "capture_species", "materials_injection": injection, "self_check": check, "use_mock": use_mock, "max_xgems_calls": max_xgems_calls}, indent=2, default=str), encoding="utf-8")
    return ForwardRunResult(ok=ok, run_dir=out, forward_dir=forward_dir, time_series=ts, status="complete" if ok else "incomplete", materials_injection=injection, self_check=check, warnings=warnings, result_files={"captures": str(forward_dir / "dorgems_captures.json")}, summary={"n_ages": len(rows)})


def _row(age: float, res: dict[str, Any], chem_dir: Path) -> dict[str, Any]:
    from inverse_gems.database import read_name_value_csv  # type: ignore

    row: dict[str, Any] = {"age_days": age, "recipe_id": res.get("recipe_id"), "chem_hash": res.get("chem_hash"), "chemistry_status": res.get("chemistry_status"), "solver_status": res.get("solver_status"), "reused_cache": res.get("reused_cache"), "porosity": res.get("porosity"), "scalar__porosity": res.get("porosity"), "xgems_water_g": res.get("xgems_water_g"), "xgems_w_b": res.get("xgems_w_b"), "xgems_water_mode": res.get("xgems_water_mode")}
    raw_dirs = sorted(chem_dir.rglob("xgems_phase_amounts_raw.csv"))
    if raw_dirs:
        raw = raw_dirs[-1].parent
        for k, v in read_name_value_csv(raw / "xgems_phase_amounts_raw.csv").items():
            row[f"phase_mass__{k}"] = v
        pv = raw / "xgems_phase_volumes_raw.csv"
        if pv.is_file():
            for k, v in read_name_value_csv(pv).items():
                row[f"phase_volume__{k}"] = v
        sc = raw / "xgems_scalars_raw.json"
        if sc.is_file():
            for k, v in (json.loads(sc.read_text(encoding="utf-8")) or {}).items():
                if isinstance(v, (int, float)) or v is None:
                    row[f"scalar__{k}"] = v
    return row


def _jsonable(o: Any) -> Any:
    if isinstance(o, (np.floating, np.integer)):
        return o.item()
    if isinstance(o, np.ndarray):
        return o.tolist()
    return str(o)
