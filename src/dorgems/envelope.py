"""Scenario A end-to-end (spec §7): predict → export 3 quantile configs → materials
override → forward ×3 → envelope.csv + summary.md + manifest.json."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from . import __version__
from .gems.forward import run_forward
from .kinetics.materials_override import build_materials_config, slot_for_role
from .kinetics.reaction_model import export_reaction_model
from .models.bundle import sha256_of
from .predict import predict, write_prediction


def build_forward_query(mix: Any, slot: str, ages: list[float], *, name: str = "dorgems_A", curing_temp_C: float | None = None) -> dict[str, Any]:
    """spec §7.2: binders sum to 100; no material_system; response_summary minimal."""
    from .db.features import _get

    scm_pct = float(_get(mix, "scm_pct"))
    others = dict(_get(mix, "other_components", {}) or {})
    opc = 100.0 - scm_pct - sum(float(v) for v in others.values())
    if opc < 0:
        raise ValueError("binder components exceed 100 %")
    binders: dict[str, float] = {"OPC": round(opc, 6)}
    if scm_pct > 0:
        binders[slot] = round(scm_pct, 6)
    for k, v in others.items():
        binders[k] = binders.get(k, 0.0) + float(v)
    T = curing_temp_C if curing_temp_C is not None else float(_get(mix, "curing_temp_C", 20.0))
    return {
        "name": name,
        "task": "forward_time_series",
        "recipe": {"binders": binders, "w_b": float(_get(mix, "w_b"))},
        "age_grid": {"values": [float(a) for a in ages]},
        "temperature_celsius": T,
        "outputs": {"phase_masses": "all", "phase_volumes": "all", "phase_volumes_reconstructed": "all", "aqueous_species": "all", "scalars": "all"},
        "plots": [],
        "response_summary": {"enabled": True, "phases": [], "scalars": ["pH", "porosity"], "narrative_enabled": False, "narrative_language": "ko"},
    }


def envelope_from_runs(runs: dict[str, pd.DataFrame]) -> pd.DataFrame:
    """Long table: age_days × variable × {q05,q50,q95}."""
    frames = []
    for q, ts in runs.items():
        if ts is None:
            continue
        ts = ts.copy()
        if "porosity" in ts.columns and "scalar__porosity" not in ts.columns:
            ts["scalar__porosity"] = ts["porosity"]
        cols = [c for c in ts.columns if c.startswith(("phase_mass__", "phase_volume__", "scalar__")) and pd.api.types.is_numeric_dtype(pd.to_numeric(ts[c], errors="coerce"))]
        for c in cols:
            ts[c] = pd.to_numeric(ts[c], errors="coerce")
        cols = [c for c in cols if ts[c].notna().any()]
        m = ts.melt(id_vars=["age_days"], value_vars=cols, var_name="variable", value_name="value")
        m["quantile"] = q
        frames.append(m)
    if not frames:
        return pd.DataFrame(columns=["age_days", "variable", "q05", "q50", "q95"])
    long = pd.concat(frames, ignore_index=True)
    wide = long.pivot_table(index=["age_days", "variable"], columns="quantile", values="value", aggfunc="first").reset_index()
    wide.columns.name = None
    return wide


def run_envelope(
    scm: Any,
    mix: Any,
    ages: list[float] | None,
    *,
    out: str | Path,
    db: str | Path,
    use_mock: bool = True,
    dat_lst: str | Path | None = None,
    max_xgems_calls: int | None = None,
    ensemble: str | None = None,
    seed: int = 0,
    export_mode: str = "logistic_fit",
    lit_db: str | Path | None = None,
    quantiles: tuple[float, ...] = (0.05, 0.5, 0.95),
) -> dict[str, Any]:
    from .pilot.schemas import coerce_mix, coerce_scm

    scm_m, mix_m = coerce_scm(scm), coerce_mix(mix)
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    warnings: list[str] = []
    pred = predict(scm_m, mix_m, ages, ensemble=ensemble, seed=seed, db_path=lit_db)
    pred_path = write_prediction(pred, out)
    ages_used = pred["input"]["ages_d"]
    slot, w_slot, reactive = slot_for_role(scm_m.role, scm_m.oxides)
    warnings += w_slot
    if not reactive:
        warnings.append(f"role {scm_m.role} is not reactive in InverseGems; the DoR model is still applied to the slot")
    mat = build_materials_config(scm_m, out, slot=slot, cement=mix_m.opc_oxides)
    warnings += mat["warnings"]
    rm = export_reaction_model(pred, out / "reaction_models", mode=export_mode, slot=slot, quantiles=quantiles, config_id=pred["id"], signature_files=[mat["path"]])
    for r in rm.values():
        warnings += r["warnings"]
    fq = build_forward_query(mix_m, slot, ages_used, name=f"dorgems_A_{pred['id']}")
    runs: dict[str, pd.DataFrame] = {}
    run_results: dict[str, Any] = {}
    for key, r in rm.items():
        res = run_forward(fq, out=out / "runs" / key, db=db, reaction_model_config=r["path"], materials_config=mat["path"], slot=slot, use_mock=use_mock, dat_lst=dat_lst, max_xgems_calls=max_xgems_calls)
        run_results[key] = res.to_tool_result()
        warnings += [f"{key}: {w}" for w in res.warnings]
        if not res.ok:
            warnings.append(f"{key}: forward run failed: {res.error}")
        runs[key] = res.time_series
    env = envelope_from_runs(runs)
    env_path = out / "envelope.csv"
    env.to_csv(env_path, index=False)
    summary_path = out / "summary.md"
    summary_path.write_text(_summary_md(pred, env, slot, mat, rm, run_results, use_mock), encoding="utf-8")
    manifest = {
        "dorgems": __version__,
        "prediction": str(pred_path),
        "prediction_provenance": pred["provenance"],
        "slot": slot,
        "materials_config": mat["path"],
        "materials_sha256": sha256_of(Path(mat["path"])),
        "reaction_models": {k: {"path": v["path"], "sha256": sha256_of(Path(v["path"]))} for k, v in rm.items()},
        "forward_query": fq,
        "runs": {k: {"ok": v["ok"], "run_dir": v["summary"]["run_dir"], "materials_injection": v["summary"].get("materials_injection"), "self_check": v["summary"].get("self_check")} for k, v in run_results.items()},
        "use_mock": use_mock,
        "max_xgems_calls": max_xgems_calls,
        "kernel_versions": _kernel_versions(),
        "warnings": warnings,
    }
    (out / "manifest.json").write_text(json.dumps(manifest, indent=2, ensure_ascii=False, default=str), encoding="utf-8")
    ok = all(v["ok"] for v in run_results.values()) and bool(run_results)
    return {"ok": ok, "out": str(out), "prediction": str(pred_path), "envelope": str(env_path), "summary": str(summary_path), "manifest": str(out / "manifest.json"), "slot": slot, "materials_config": mat["path"], "reaction_models": {k: v["path"] for k, v in rm.items()}, "runs": run_results, "warnings": warnings, "n_rows_envelope": int(len(env))}


def _kernel_versions() -> dict[str, Any]:
    out: dict[str, Any] = {}
    for name in ("inverse_gems", "gemspilot"):
        try:
            mod = __import__(name)
            out[name] = {"path": str(Path(mod.__file__).resolve().parents[2]), "git": _git_head(Path(mod.__file__).resolve().parents[2])}
        except Exception:  # noqa: BLE001
            out[name] = None
    return out


def _git_head(root: Path) -> str | None:
    head = root / ".git" / "HEAD"
    if not head.is_file():
        return None
    ref = head.read_text().strip()
    if ref.startswith("ref:"):
        p = root / ".git" / ref.split(" ", 1)[1]
        return p.read_text().strip()[:12] if p.is_file() else None
    return ref[:12]


def _summary_md(pred: dict[str, Any], env: pd.DataFrame, slot: str, mat: dict[str, Any], rm: dict[str, Any], runs: dict[str, Any], use_mock: bool) -> str:
    """Template text only; every number is rendered from the JSON/CSV objects."""
    rec = pred["recommended"]
    ages = pred["input"]["ages_d"]
    lines = [
        f"# DoRGems 시나리오 A 요약 — {pred['input']['scm']['name']}",
        "",
        f"- 역할 `{pred['input']['scm']['role']}` → InverseGems 슬롯 `{slot}` (materials override: `{Path(mat['path']).name}`)",
        f"- DoR 예측 출처: `{rec['source']}`; 재령(일): {ages}",
        f"- 실행 모드: {'mock' if use_mock else 'real xGEMS'}; 분위 실행 {', '.join(f'{k}: {'ok' if v['ok'] else 'FAIL'}' for k, v in runs.items())}",
        "",
        "| 재령 d | α q05 % | α q50 % | α q95 % |",
        "|---|---|---|---|",
    ]
    for a, lo, md, hi in zip(ages, rec["alpha_pct_q05"], rec["alpha_pct_q50"], rec["alpha_pct_q95"]):
        lines.append(f"| {a:g} | {lo:.1f} | {md:.1f} | {hi:.1f} |")
    lines += ["", f"- OOD: flags={pred['ood']['flags']}, score_pct={pred['ood'].get('score_pct')}, sparse_role={pred['ood']['sparse_role']}", f"- 근거 유사배합: {pred['evidence'].get('n_mixes', 0)}배합 / {pred['evidence'].get('n_papers', 0)}편 (prediction.json → evidence)", ""]
    if not env.empty:
        keys = [v for v in env["variable"].unique() if v.startswith("scalar__porosity") or v.startswith("scalar__pH") or "Portlandite" in v or "CNASH" in v or "C-S-H" in v]
        sub = env[env["variable"].isin(keys)]
        lines += ["| 재령 d | 변수 | q05 | q50 | q95 |", "|---|---|---|---|---|"]
        for _, r in sub.iterrows():
            lines.append(f"| {r['age_days']:g} | {r['variable']} | {_fmt(r.get('q05'))} | {_fmt(r.get('q50'))} | {_fmt(r.get('q95'))} |")
    lines += ["", "모든 수치의 원본: `prediction.json`, `envelope.csv`, `runs/*/forward/time_series.csv`. 경고는 `manifest.json`."]
    return "\n".join(lines) + "\n"


def _fmt(v: Any) -> str:
    try:
        return f"{float(v):.4g}"
    except (TypeError, ValueError):
        return "—"


__all__ = ["run_envelope", "build_forward_query", "envelope_from_runs"]
