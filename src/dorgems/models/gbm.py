"""LightGBM point predictions from the gbm_v6 bundle (spec §5.3).

Feature order and categorical codes come from ``meta.json`` so inference encodes
exactly like training (``export_bundle.encode_frame``).
"""

from __future__ import annotations

from typing import Any

import numpy as np
import pandas as pd

from ..db.features import CAT_FEATURES, GBM_FEATURES, derive_composition_features
from .bundle import GBMBundle

MISSING = "missing"


def encode_frame(df: pd.DataFrame, features: list[str], categories: dict[str, list[str]]) -> pd.DataFrame:
    """Categorical → integer codes (unknown/NaN → the 'missing' code). Same function at train and predict."""
    X = pd.DataFrame(index=df.index)
    for f in features:
        if f in categories:
            cats = categories[f]
            miss = cats.index(MISSING) if MISSING in cats else -1
            codes = pd.Categorical(df[f].astype(object).where(df[f].notna(), MISSING), categories=cats).codes.astype(int)
            codes = np.where(codes < 0, miss, codes)
            X[f] = codes
        else:
            X[f] = pd.to_numeric(df[f], errors="coerce").astype(float)
    return X


def predict_points(
    bundle: GBMBundle,
    x: dict[str, Any],
    role: str | None,
    ages: np.ndarray | list[float],
    *,
    method_group: str = MISSING,
) -> tuple[np.ndarray, dict[str, Any]]:
    """(len(ages),) DoR in %, clipped to [0,100], plus meta (method_group_used, role_used, flags)."""
    ages = np.asarray(ages, dtype=float)
    flags: list[str] = []
    role_used = role if role in bundle.categories.get("scm_role", []) else MISSING
    if role_used == MISSING and role is not None:
        flags.append(f"gbm_role_unseen:{role}")
    mg = method_group if method_group in bundle.categories.get("method_group", []) else MISSING
    if mg != method_group:
        flags.append(f"gbm_method_group_unseen:{method_group}")
    rows = []
    for a in ages:
        r = dict(x)
        r["age_d"] = float(a)
        r["scm_role"] = role_used
        r["method_group"] = mg
        rows.append(r)
    df = pd.DataFrame(rows)
    for c in ["scm_CaO", "scm_SiO2", "scm_Al2O3", "scm_Fe2O3", "scm_MgO", "scm_amorphous_pct", "scm_blaine_m2_kg"]:
        if c not in df.columns:
            df[c] = np.nan
    df = derive_composition_features(df)
    X = encode_frame(df, bundle.features, bundle.categories)
    pred = bundle.booster.predict(X.values, num_threads=1)
    pred = np.clip(np.asarray(pred, dtype=float), 0.0, 100.0)
    return pred, {"method_group_used": mg, "role_used": role_used, "flags": flags}


__all__ = ["predict_points", "encode_frame", "GBM_FEATURES", "CAT_FEATURES", "MISSING"]
