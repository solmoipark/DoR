"""Forward wrapper around ``inverse_gems.api`` (spec §6.3).

* writes the forward_query dict to ``out/forward_query.yaml`` first (the kernel
  takes a path, not a dict — forward_query.py:407-408);
* injects ``materials_config`` either natively (P-IG-1, if the kernel accepts
  the keyword) or by temporarily re-pointing ``inverse_gems.materials.config_path``
  for ``materials.yaml`` (fallback; recorded in the manifest as ``monkeypatch``);
* self-checks after the run: ``input_reaction_degrees.json["scm"][slot]`` must
  equal the exported alpha(age) within 1e-3 and ``input_materials_used.json``
  must carry the override oxides;
* the ``capture_species`` path drives ``run_forward_cached`` per age with a
  capturing runner factory (bound-water calculation, §8.3).
"""

from __future__ import annotations

import contextlib
import inspect
import json
from collections.abc import Iterator
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import yaml

from .. import TOOL_CONTRACT, __version__
from ..kinetics.reaction_model import alpha_from_config
from ..kinetics.registry import register as _register_kinetics

HARD_CAP_XGEMS_CALLS = 200
SCM_SLOTS = ("slag", "fly_ash", "metakaolin", "silica_fume")


@dataclass
class ForwardRunResult:
    ok: bool
    run_dir: Path
    forward_dir: Path
    time_series: pd.DataFrame | None
    status: str
    materials_injection: str | None
    self_check: dict[str, Any]
    warnings: list[str] = field(default_factory=list)
    error: str | None = None
    result_files: dict[str, str] = field(default_factory=dict)
    summary: dict[str, Any] = field(default_factory=dict)

    def to_tool_result(self, tool: str = "dor_run_forward_with_dor") -> dict[str, Any]:
        return {
            "contract": TOOL_CONTRACT,
            "tool": tool,
            "ok": self.ok,
            "summary": {
                "status": self.status,
                "task_type": "forward_time_series",
                "run_dir": str(self.run_dir),
                "answer_available": self.time_series is not None,
                "missing_outputs": {},
                "result_summary": self.summary,
                "materials_injection": self.materials_injection,
                "self_check": self.self_check,
            },
            "artifacts": {"time_series": str(self.forward_dir / "time_series.csv"), **self.result_files},
            "warnings": list(self.warnings),
            "error": self.error,
        }


# ---------------------------------------------------------------------------
# materials override injection
# ---------------------------------------------------------------------------


def kernel_accepts_materials_config() -> bool:
    try:
        from inverse_gems.api import run_forward_request
    except ImportError:
        return False
    return "materials_config" in inspect.signature(run_forward_request).parameters


@contextlib.contextmanager
def materials_config_override(path: str | Path | None) -> Iterator[str | None]:
    """Fallback injection: re-point ``inverse_gems.materials.config_path`` for
    ``materials.yaml`` only. ``load_materials`` itself is *not* patched (each module
    binds it at import time); ``config_path`` is resolved inside ``load_materials``
    at call time, so the patch reaches every caller."""
    if path is None:
        yield None
        return
    import inverse_gems.materials as im

    target = Path(path).resolve()
    if not target.is_file():
        raise FileNotFoundError(target)
    original = im.config_path

    def patched(filename: str) -> Path:
        if filename == "materials.yaml":
            return target
        return original(filename)

    im.config_path = patched  # type: ignore[assignment]
    try:
        yield "monkeypatch"
    finally:
        im.config_path = original  # type: ignore[assignment]


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------


def write_forward_query(forward_query: dict[str, Any], out: Path) -> Path:
    out.mkdir(parents=True, exist_ok=True)
    p = out / "forward_query.yaml"
    p.write_text(yaml.safe_dump(forward_query, sort_keys=False, allow_unicode=True), encoding="utf-8")
    return p


OXIDE_MOLAR_MASS = {"SiO2": 60.084, "Al2O3": 101.961, "Fe2O3": 159.688, "CaO": 56.077, "MgO": 40.304, "SO3": 80.063, "Na2O": 61.979, "K2O": 94.196}


def _load_reaction_degrees(forward_dir: Path, db_dir: Path | None) -> list[dict[str, Any]]:
    """Every run's reaction degrees. Re-verified 2026-09-02: the forward_query /
    cached_forward path writes ``recipe_runs/<recipe_id>/reaction_degrees.json``
    (cached_forward.py:632) — ``input_reaction_degrees.json`` and
    ``input_materials_used.json`` (xgems_output_capture.py) are only produced by the
    stand-alone hydration script path. Materials actually used are recovered from
    ``source_contribution_ledger.csv`` (per-material oxide moles)."""
    out = []
    roots = [forward_dir] + ([db_dir] if db_dir else [])
    seen: set[Path] = set()
    for root in roots:
        if not root or not root.exists():
            continue
        for name in ("reaction_degrees.json", "input_reaction_degrees.json"):
            for p in root.rglob(name):
                if p in seen:
                    continue
                seen.add(p)
                try:
                    rd = json.loads(p.read_text(encoding="utf-8"))
                except Exception:  # noqa: BLE001
                    continue
                mu = p.parent / "input_materials_used.json"
                materials_used = json.loads(mu.read_text(encoding="utf-8")) if mu.is_file() else _materials_from_ledger(p.parent / "source_contribution_ledger.csv")
                out.append({"path": p, "age_days": rd.get("age_days"), "scm": rd.get("scm", {}), "opc": rd.get("opc"), "materials_used": materials_used, "reaction_parameter_set": rd.get("reaction_parameter_set")})
    return out


def _materials_from_ledger(ledger: Path) -> dict[str, Any] | None:
    """{material: {"oxide_mass_percent": {...}}} recovered from the source ledger:
    pct_k = oxide_equivalent_mol_k × M_k / reacted_mass_g × 100 (per material)."""
    if not ledger.is_file():
        return None
    try:
        df = pd.read_csv(ledger)
    except Exception:  # noqa: BLE001
        return None
    out: dict[str, Any] = {}
    for mat, g in df.groupby("source_material"):
        g2 = g[g["oxide_equivalent"].notna() & (g["oxide_equivalent"].astype(str) != "")]
        if g2.empty:
            continue
        reacted = float(g2["reacted_mass_g"].sum()) if "source_phase_or_oxide" in g2 and g2["source_phase_or_oxide"].nunique() > 1 else None
        # per (phase/oxide) rows share reacted_mass_g; take the sum over distinct sources
        per_src = g2.drop_duplicates("source_phase_or_oxide")[["source_phase_or_oxide", "reacted_mass_g"]]
        reacted = float(per_src["reacted_mass_g"].sum())
        if reacted <= 0:
            continue
        pct: dict[str, float] = {}
        for ox, gg in g2.groupby("oxide_equivalent"):
            M = OXIDE_MOLAR_MASS.get(str(ox))
            if M is None:
                continue
            pct[str(ox)] = float(gg["oxide_equivalent_mol"].sum()) * M / reacted * 100.0
        out[str(mat)] = {"oxide_mass_percent": pct, "reacted_mass_g": reacted}
    return out


def self_check(
    forward_dir: Path,
    *,
    reaction_model_config: str | Path | None,
    slot: str | None,
    ages: list[float],
    materials_config: str | Path | None,
    db_dir: Path | None = None,
    tol: float = 1e-3,
) -> dict[str, Any]:
    """Verify that the config actually reached the kernel (spec §6.3)."""
    result: dict[str, Any] = {"alpha_ok": None, "materials_ok": None, "details": []}
    runs = _load_reaction_degrees(forward_dir, db_dir)
    cfg_id = None
    if reaction_model_config:
        cfg = yaml.safe_load(Path(reaction_model_config).read_text(encoding="utf-8")) or {}
        cfg_id = cfg.get("id")
    mine = [r for r in runs if cfg_id is None or (r.get("reaction_parameter_set") or {}).get("id") == cfg_id]
    if reaction_model_config and slot:
        expected = alpha_from_config(reaction_model_config, slot, np.asarray(ages, float))
        exp_by_age = {float(a): float(v) for a, v in zip(ages, expected)}
        ok_all = bool(mine)
        for r in mine:
            age = float(r["age_days"]) if r.get("age_days") is not None else None
            got = (r.get("scm") or {}).get(slot)
            if age is None or age not in exp_by_age or got is None:
                continue
            d = abs(float(got) - exp_by_age[age])
            ok = d <= tol
            ok_all &= ok
            result["details"].append({"age_days": age, "alpha_expected": exp_by_age[age], "alpha_in_kernel": float(got), "ok": ok})
        result["alpha_ok"] = ok_all and any(d["ok"] for d in result["details"]) if result["details"] else False
    if materials_config and slot:
        want = yaml.safe_load(Path(materials_config).read_text(encoding="utf-8"))[slot]["oxide_mass_percent"]
        ok_all = bool(mine)
        found = False
        for r in mine:
            mu = r.get("materials_used") or {}
            used = (mu.get(slot) or {}).get("oxide_mass_percent")
            if used is None:
                continue
            found = True
            # ledger-recovered compositions carry float noise; 0.05 %p tolerance on each oxide
            diffs = {k: abs(float(used.get(k, 0.0)) - float(v)) for k, v in want.items()}
            ok = all(d < 0.05 for d in diffs.values())
            ok_all &= ok
            result["details"].append({"materials_check": {"age_days": r.get("age_days"), "max_abs_diff_pct": max(diffs.values()), "ok": ok}})
        result["materials_ok"] = found and ok_all
    return result


def _budget_guard(use_mock: bool, max_xgems_calls: int | None) -> None:
    if not use_mock and (max_xgems_calls is None or max_xgems_calls > HARD_CAP_XGEMS_CALLS):
        raise PermissionError(f"real xGEMS runs need an explicit max_xgems_calls <= {HARD_CAP_XGEMS_CALLS} (got {max_xgems_calls!r})")


# ---------------------------------------------------------------------------
# main entry
# ---------------------------------------------------------------------------


def run_forward(
    forward_query: dict[str, Any],
    *,
    out: str | Path,
    db: str | Path,
    reaction_model_config: str | Path | None,
    materials_config: str | Path | None = None,
    slot: str | None = None,
    use_mock: bool = True,
    dat_lst: str | Path | None = None,
    max_xgems_calls: int | None = None,
    capture_species: bool = False,
    reaction_model_id: str | None = None,
) -> ForwardRunResult:
    _budget_guard(use_mock, max_xgems_calls)
    _register_kinetics()
    out = Path(out)
    warnings: list[str] = []
    query_path = write_forward_query(forward_query, out)
    ages = list(map(float, (forward_query.get("age_grid") or {}).get("values") or []))
    binders = (forward_query.get("recipe") or {}).get("binders") or {}
    if slot is None:
        slot = next((s for s in SCM_SLOTS if float(binders.get(s, 0) or 0) > 0), None)
    if capture_species:
        from .capture import run_forward_capturing

        return run_forward_capturing(forward_query, out=out, db=db, reaction_model_config=reaction_model_config, materials_config=materials_config, slot=slot, use_mock=use_mock, dat_lst=dat_lst, max_xgems_calls=max_xgems_calls, reaction_model_id=reaction_model_id)

    from inverse_gems.api import run_forward_request

    native = kernel_accepts_materials_config()
    injection: str | None = None
    kwargs: dict[str, Any] = dict(
        forward_query=query_path,
        out=out,
        db=db,
        dat_lst=dat_lst,
        use_mock=use_mock,
        max_xgems_calls=max_xgems_calls,
        reaction_model_id=reaction_model_id,
        reaction_model_config=str(reaction_model_config) if reaction_model_config else None,
        disable_plots=True,
    )
    try:
        if materials_config and native:
            kwargs["materials_config"] = str(materials_config)
            injection = "native"
            result = run_forward_request(**kwargs)
        else:
            with materials_config_override(materials_config) as how:
                injection = how
                result = run_forward_request(**kwargs)
        if injection == "monkeypatch":
            warnings.append("materials_config injected by monkeypatching inverse_gems.materials.config_path (kernel patch P-IG-1 not merged)")
    except Exception as exc:  # noqa: BLE001
        return ForwardRunResult(ok=False, run_dir=out, forward_dir=out / "forward", time_series=None, status="error", materials_injection=injection, self_check={}, warnings=warnings, error=f"{type(exc).__name__}: {exc}")

    data = result.to_dict()
    forward_dir = out / "forward"
    ts_path = forward_dir / "time_series.csv"
    ts = pd.read_csv(ts_path) if ts_path.is_file() else None
    check = self_check(forward_dir, reaction_model_config=reaction_model_config, slot=slot, ages=ages, materials_config=materials_config, db_dir=Path(db))
    ok = data.get("status") == "complete"
    if reaction_model_config and check.get("alpha_ok") is False:
        ok = False
        warnings.append("self-check failed: kernel reaction degrees do not match the exported reaction model")
    if materials_config and check.get("materials_ok") is False:
        ok = False
        warnings.append("self-check failed: input_materials_used.json does not carry the override oxides")
    summary = dict(data.get("summary") or {})
    manifest = {
        "dorgems": __version__,
        "materials_injection": injection,
        "reaction_model_config": str(reaction_model_config) if reaction_model_config else None,
        "materials_config": str(materials_config) if materials_config else None,
        "use_mock": use_mock,
        "max_xgems_calls": max_xgems_calls,
        "self_check": check,
        "kernel_status": data.get("status"),
    }
    (out / "dorgems_forward_manifest.json").write_text(json.dumps(manifest, indent=2, default=str), encoding="utf-8")
    return ForwardRunResult(ok=ok, run_dir=out, forward_dir=forward_dir, time_series=ts, status=str(data.get("status")), materials_injection=injection, self_check=check, warnings=warnings + [str(w) for w in data.get("warnings") or []], error=data.get("error"), result_files={k: str(v) for k, v in (data.get("files") or {}).items()}, summary=summary)


def series_columns(ts: pd.DataFrame) -> dict[str, list[str]]:
    return {
        "phase_mass": [c for c in ts.columns if c.startswith("phase_mass__")],
        "phase_volume": [c for c in ts.columns if c.startswith("phase_volume__")],
        "scalar": [c for c in ts.columns if c.startswith("scalar__")],
    }
