"""Nearest literature mixes for a new SCM/mix (spec §5.5, §8.2).

Tolerances live in ``configs/analogue_tolerances.yaml``; the distance is
    d² = Σ_c ((Δc) / scale_c)²  over standardized chemistry + |Δscm_pct| + |Δw_b| (+ |ΔT|)
and the weight is ``w = exp(−d²/2)``.
"""

from __future__ import annotations

import sqlite3
from typing import Any

import numpy as np
import pandas as pd
import yaml

from ..config import configs_dir
from .features import blended_only, build_dor_table, derive_composition_features
from .reader import MODEL_SYSTEM

DEFAULT_TOL = {
    "scm_pct": 10.0,
    "w_b": 0.05,
    "curing_temp_C": 5.0,
    "age_ratio": [0.75, 1.33],
    "chem_distance_max": 1.5,
    "chem_features": ["CaO_SiO2", "Al_Si", "basicity"],
    "default_curing_temp_C": 20.0,
}


def load_tolerances() -> dict[str, Any]:
    p = configs_dir() / "analogue_tolerances.yaml"
    tol = dict(DEFAULT_TOL)
    if p.is_file():
        tol.update(yaml.safe_load(p.read_text(encoding="utf-8")) or {})
    return tol


_TABLE_CACHE: dict[str, pd.DataFrame] = {}


def dor_table_with_features(con: sqlite3.Connection, key: str | None = None) -> pd.DataFrame:
    if key and key in _TABLE_CACHE:
        return _TABLE_CACHE[key]
    df = derive_composition_features(blended_only(build_dor_table(con)))
    if key:
        _TABLE_CACHE[key] = df
    return df


def _role_scales(df: pd.DataFrame, chem: list[str]) -> dict[str, float]:
    scales = {}
    for c in chem:
        s = float(pd.to_numeric(df[c], errors="coerce").std())
        scales[c] = s if s > 0 and np.isfinite(s) else 1.0
    return scales


def find_analogues(
    con: sqlite3.Connection,
    feats: dict[str, Any],
    role: str | None,
    *,
    k: int = 5,
    age_days: float | None = None,
    same_role_only: bool = True,
    tol: dict[str, Any] | None = None,
    cache_key: str | None = None,
) -> dict[str, Any]:
    """Return {mixes: [...], n_mixes, n_papers, flags}. ``feats`` needs scm_pct, w_b,
    curing_temp_C (optional) and the chemistry features from derive_composition_features."""
    tol = {**load_tolerances(), **(tol or {})}
    df = dor_table_with_features(con, cache_key)
    flags: list[str] = []
    if same_role_only and role is not None:
        sub = df[df["scm_role"] == role]
        if sub.empty:
            flags.append(f"no_training_rows_for_role:{role}")
            sub = df
    else:
        sub = df
    chem = [c for c in tol["chem_features"] if c in sub.columns]
    scales = _role_scales(sub, chem)
    T_default = float(tol["default_curing_temp_C"])

    # one row per mix (features are mix-level)
    mix_cols = ["mix_uid", "paper_doi", "scm_role", "scm_pct", "w_b", "curing_temp_C"] + chem
    mixes = sub.groupby("mix_uid", as_index=False)[mix_cols[1:]].first()
    d2 = np.zeros(len(mixes))
    used_terms = 0
    x_scm = feats.get("scm_pct")
    x_wb = feats.get("w_b")
    x_T = feats.get("curing_temp_C")
    if x_scm is not None:
        d2 += ((pd.to_numeric(mixes["scm_pct"], errors="coerce").fillna(x_scm) - x_scm) / tol["scm_pct"]) ** 2
        used_terms += 1
    if x_wb is not None:
        d2 += ((pd.to_numeric(mixes["w_b"], errors="coerce").fillna(x_wb) - x_wb) / tol["w_b"]) ** 2
        used_terms += 1
    T_ref = pd.to_numeric(mixes["curing_temp_C"], errors="coerce")
    t_missing = T_ref.isna()
    if x_T is None:
        x_T = T_default
    d2 += ((T_ref.fillna(T_default) - x_T) / tol["curing_temp_C"]) ** 2
    chem_d2 = np.zeros(len(mixes))
    chem_terms = 0
    for c in chem:
        v = feats.get(c)
        if v is None or (isinstance(v, float) and np.isnan(v)):
            continue
        ref = pd.to_numeric(mixes[c], errors="coerce")
        chem_d2 += ((ref.fillna(v) - v) / scales[c]) ** 2
        chem_terms += 1
    chem_dist = np.sqrt(chem_d2 / max(chem_terms, 1))
    d2 += chem_d2
    mixes = mixes.assign(d2=d2, chem_distance=chem_dist, weight=np.exp(-d2 / 2.0), curing_temp_missing=t_missing.values)
    within = mixes[(mixes["chem_distance"] <= tol["chem_distance_max"]) | (chem_terms == 0)]
    if within.empty:
        flags.append("no_mix_within_chemistry_tolerance")
        within = mixes
    top = within.sort_values("d2").head(int(k))

    out_mixes: list[dict[str, Any]] = []
    for _, m in top.iterrows():
        obs = sub[sub["mix_uid"] == m["mix_uid"]]
        if age_days is not None:
            lo, hi = tol["age_ratio"]
            obs_age = obs[(obs["age_d"] / float(age_days) >= lo) & (obs["age_d"] / float(age_days) <= hi)]
            if not obs_age.empty:
                obs = obs_age
        out_mixes.append(
            {
                "mix_uid": m["mix_uid"],
                "paper_doi": m["paper_doi"],
                "scm_role": m["scm_role"],
                "scm_pct": _f(m["scm_pct"]),
                "w_b": _f(m["w_b"]),
                "curing_temp_C": _f(m["curing_temp_C"]),
                "chem_distance": float(m["chem_distance"]),
                "distance": float(np.sqrt(m["d2"])),
                "weight": float(m["weight"]),
                "dor": [
                    {"age_d": float(o["age_d"]), "value_pct": float(o["dor_pct"]), "method_group": o["method_group"], "obs_uid": o["obs_uid"], "fig_only": _i(o["fig_only"]), "confidence": o["confidence"]}
                    for _, o in obs.sort_values("age_d").iterrows()
                ],
            }
        )
    n_papers = int(top["paper_doi"].nunique())
    if n_papers == 1:
        flags.append("single_paper")
    if bool(top["curing_temp_missing"].any()):
        flags.append("curing_temp_assumed_20C_for_some_analogues")
    return {"mixes": out_mixes, "n_mixes": int(len(top)), "n_papers": n_papers, "flags": flags, "tolerances": {k_: v for k_, v in tol.items()}}


def _f(v: Any) -> float | None:
    try:
        f = float(v)
    except (TypeError, ValueError):
        return None
    return None if np.isnan(f) else f


def _i(v: Any) -> int | None:
    try:
        return None if v is None or (isinstance(v, float) and np.isnan(v)) else int(v)
    except (TypeError, ValueError):
        return None


__all__ = ["find_analogues", "load_tolerances", "dor_table_with_features", "MODEL_SYSTEM"]
