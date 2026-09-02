"""Residuals, z-scores and the deterministic verdict (spec §8.5).

    r = model − obs_harmonised,  z = r / sqrt(σ_obs² + σ_model²)
    consistent : frac(|z|<2) ≥ 0.7 and |median r| ≤ σ_model     (grade A/B only)
    tension    : otherwise, n ≥ 5
    insufficient_data : n < 5
Verdict text is a template; numbers come only from comparison.csv/json.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from scipy.stats import binomtest

from ..predict import load_defaults

OBS_TO_MODEL = {"CH_TGA": "CH_g", "CH_XRD": "CH_g", "bound_water": "bound_water_g", "chem_shrink": "chem_shrink_ml_g"}


def sigma_tables() -> tuple[dict[str, float], dict[str, float]]:
    d = load_defaults()
    return dict(d.get("sigma_obs_default", {})), dict(d.get("sigma_model_initial", {}))


def compare_rows(pairs: list[dict[str, Any]], *, sigma_model: dict[str, float] | None = None, sigma_obs_default: dict[str, float] | None = None, offsets: dict[str, float] | None = None) -> pd.DataFrame:
    """``pairs``: dicts with obs_uid, quantity, age_d, obs_value (harmonised), grade,
    assumptions, model_value, uncertainty (optional), paper_doi, mix_uid, phase_name."""
    so_def, sm_def = sigma_tables()
    so_def.update(sigma_obs_default or {})
    sm_def.update(sigma_model or {})
    offsets = offsets or {}
    rows = []
    for p in pairs:
        q = p["quantity"]
        obs, mod = p.get("obs_value"), p.get("model_value")
        b = float(offsets.get(q, 0.0))
        so = float(p["uncertainty"]) if p.get("uncertainty") not in (None, 0) else float(so_def.get(q, np.nan))
        sm = float(sm_def.get(q, np.nan))
        r = z = None
        if obs is not None and mod is not None and np.isfinite(obs) and np.isfinite(mod):
            r = float(mod + b - obs)
            denom = np.sqrt(so**2 + sm**2)
            z = float(r / denom) if np.isfinite(denom) and denom > 0 else None
        rows.append({**{k: p.get(k) for k in ("obs_uid", "paper_doi", "mix_uid", "quantity", "phase_name", "age_d", "method", "grade", "assumptions", "source_locator", "fig_only", "extraction_confidence")}, "obs": obs, "model": mod, "offset_b": b, "r": r, "sigma_obs": so, "sigma_model": sm, "z": z, "usable": p.get("grade") in ("A", "B") and r is not None})
    return pd.DataFrame(rows)


def aggregate(df: pd.DataFrame, *, min_n: int | None = None, frac_thr: float | None = None) -> dict[str, Any]:
    d = load_defaults().get("compare", {})
    min_n = int(min_n if min_n is not None else d.get("min_n", 5))
    frac_thr = float(frac_thr if frac_thr is not None else d.get("consistent_frac_z_lt2", 0.7))
    out: dict[str, Any] = {"by_quantity": {}, "verdict": {}}
    if df.empty:
        return {**out, "overall": "insufficient_data", "n_usable": 0}
    for q, g in df.groupby("quantity"):
        u = g[g["usable"]]
        n = int(len(u))
        entry: dict[str, Any] = {"n": n, "n_all_grades": int(len(g)), "n_papers": int(u["paper_doi"].nunique()) if n else 0, "grades": g["grade"].value_counts().to_dict()}
        if n:
            r = u["r"].astype(float)
            z = u["z"].astype(float)
            entry.update({"median_r": float(r.median()), "iqr_r": float(r.quantile(0.75) - r.quantile(0.25)), "mean_r": float(r.mean()), "frac_abs_z_lt2": float((z.abs() < 2).mean()), "sigma_model": float(u["sigma_model"].iloc[0]), "rmse": float(np.sqrt((r**2).mean()))})
            pos = int((r > 0).sum())
            nz = int((r != 0).sum())
            entry["bias_sign_test_p"] = float(binomtest(pos, nz, 0.5).pvalue) if nz else None
            if n >= min_n:
                verdict = "consistent" if (entry["frac_abs_z_lt2"] >= frac_thr and abs(entry["median_r"]) <= entry["sigma_model"]) else "tension"
            else:
                verdict = "insufficient_data"
        else:
            verdict = "insufficient_data"
        entry["verdict"] = verdict
        out["by_quantity"][q] = entry
        out["verdict"][q] = verdict
    vs = list(out["verdict"].values())
    out["overall"] = "tension" if "tension" in vs else ("consistent" if "consistent" in vs else "insufficient_data")
    out["n_usable"] = int(df["usable"].sum())
    return out


def write_comparison(df: pd.DataFrame, agg: dict[str, Any], out: Path, *, header: dict[str, Any] | None = None) -> dict[str, str]:
    out.mkdir(parents=True, exist_ok=True)
    csv = out / "comparison.csv"
    df.to_csv(csv, index=False)
    js = out / "comparison.json"
    payload = {"schema": "dorgems-comparison/1.0", **(header or {}), "aggregate": agg}
    js.write_text(json.dumps(payload, indent=2, ensure_ascii=False, default=_jd), encoding="utf-8")
    md = out / "summary.md"
    md.write_text(summary_md(agg, header or {}), encoding="utf-8")
    return {"comparison_csv": str(csv), "comparison_json": str(js), "summary": str(md)}


def summary_md(agg: dict[str, Any], header: dict[str, Any]) -> str:
    lines = [f"# 문헌 대조 요약 — mode `{header.get('mode', '?')}`", "", f"- 대상: {header.get('target', '')}", f"- 종합 판정: **{agg.get('overall')}** (등급 A·B 관측 {agg.get('n_usable', 0)}건)", "", "| 물리량 | n(A·B) | 편수 | median r | IQR r | frac |z|<2 | σ_model | 판정 |", "|---|---|---|---|---|---|---|---|"]
    for q, e in agg.get("by_quantity", {}).items():
        lines.append(f"| {q} | {e['n']} | {e.get('n_papers', 0)} | {_f(e.get('median_r'))} | {_f(e.get('iqr_r'))} | {_f(e.get('frac_abs_z_lt2'))} | {_f(e.get('sigma_model'))} | {e['verdict']} |")
    lines += ["", "판정 문구는 템플릿이며 수치는 `comparison.csv`/`comparison.json`에서만 온다. 등급 C·D 관측은 참고용으로 CSV에 남지만 통계에서 제외된다."]
    return "\n".join(lines) + "\n"


def _f(v: Any) -> str:
    try:
        return "—" if v is None else f"{float(v):.3g}"
    except (TypeError, ValueError):
        return "—"


def _jd(o: Any) -> Any:
    if isinstance(o, (np.floating, np.integer)):
        return o.item()
    if isinstance(o, (np.bool_,)):
        return bool(o)
    return str(o)
