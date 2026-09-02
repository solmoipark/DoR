"""One-off bundle export (spec §5.1). Needs extras[export] (arviz, lightgbm, scikit-learn).

    python -m dorgems.models.export_bundle bayes  --idata modeling/work/bayes_idata.nc --table modeling/dor_scm_final.csv
    python -m dorgems.models.export_bundle gbm    --table modeling/dor_scm_blended.csv [--lopo]
    python -m dorgems.models.export_bundle ood    --table modeling/dor_scm_blended.csv

The training scripts (bayes_hier_v4.py, multitask_v6.py) are not modified; this
module only freezes what inference needs.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from .. import __version__
from ..config import bundles_dir, modeling_dir
from ..db.features import BAYES_FEATURES, CAT_FEATURES, GBM_FEATURES, derive_composition_features
from .bundle import sha256_of, write_manifest
from .gbm import MISSING, encode_frame

TIGHT = dict(
    objective="regression",
    num_leaves=7,
    max_depth=4,
    learning_rate=0.03,
    n_estimators=500,
    min_child_samples=40,
    subsample=0.7,
    subsample_freq=1,
    colsample_bytree=0.6,
    reg_lambda=20.0,
    reg_alpha=2.0,
    verbose=-1,
)
N_DRAWS = 2000


# ---------------------------------------------------------------------------
# Bayes
# ---------------------------------------------------------------------------


def bayes_training_frame(table: Path) -> pd.DataFrame:
    """bayes_hier_v4.py:50-53 — rows with age_d > 0 and CaO/SiO2 ratio."""
    df = pd.read_csv(table)
    df = df[df["age_d"] > 0].copy().reset_index(drop=True)
    C, S = [pd.to_numeric(df["scm_" + k], errors="coerce") for k in ["CaO", "SiO2"]]
    df["CaO_SiO2"] = C / S
    return df


def compute_scaler(df: pd.DataFrame) -> dict[str, Any]:
    """bayes_hier_v4.py:55-65 — median-impute → mean/std; role/method category order."""
    med, mean, std = {}, {}, {}
    for c in BAYES_FEATURES:
        col = df[c].copy()
        m = float(col.median())
        col = col.fillna(m)
        med[c] = m
        mean[c] = float(col.mean())
        std[c] = float(col.std())  # pandas ddof=1, as in the script
    methods = list(pd.Categorical(df["method_group"]).categories)
    _r = df["scm_role"].fillna("unknown")
    roles = list(pd.Categorical(_r.where(_r.isin(["slag", "fly_ash"]), "other")).categories)
    return {
        "feats": BAYES_FEATURES,
        "median": med,
        "mean": mean,
        "std": std,
        "roles": roles,
        "methods": methods,
        "beta_shape": 0.5,
        "n_rows": int(len(df)),
        "n_papers": int(df["paper_doi"].nunique()),
    }


def export_bayes(idata_path: Path, table: Path, out: Path, *, n_draws: int = N_DRAWS, seed: int = 0, db_path: Path | None = None) -> Path:
    import arviz as az

    idata = az.from_netcdf(idata_path)
    po = idata.posterior
    df = bayes_training_frame(table)
    scaler = compute_scaler(df)
    n_r, n_m, n_f = len(scaler["roles"]), len(scaler["methods"]), len(BAYES_FEATURES)
    if po["a0_role"].shape[-1] != n_r:
        raise ValueError(f"idata a0_role has {po['a0_role'].shape[-1]} roles, table implies {n_r}")
    if po["sigma_method"].shape[-1] != n_m:
        raise ValueError(f"idata sigma_method has {po['sigma_method'].shape[-1]} methods, table implies {n_m}")

    def flat(name: str) -> np.ndarray:
        v = po[name].values  # (chain, draw, ...)
        return v.reshape((-1,) + v.shape[2:])

    total = flat("a0_role").shape[0]
    rng = np.random.default_rng(seed)
    idx = np.sort(rng.choice(total, size=min(n_draws, total), replace=False)) if total > n_draws else np.arange(total)
    arrays = {
        "a0_role": flat("a0_role")[idx],
        "t0_role": flat("t0_role")[idx],
        "beta_a": flat("beta_a")[idx],
        "beta_t": flat("beta_t")[idx],
        "sd_paper_amax": flat("sd_paper_amax")[idx],
        "sigma_method": flat("sigma_method")[idx],
        "draw_index": idx,
    }
    out.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(out / "posterior.npz", **arrays)
    (out / "scaler.json").write_text(json.dumps(scaler, indent=2), encoding="utf-8")

    summ = az.summary(idata, var_names=["a0_role", "t0_role", "sd_paper_amax", "beta_a", "beta_t", "sigma_method"])
    extra = {
        "model": "bayes_hier_v4",
        "training_table": table.name,
        "training_table_sha256": sha256_of(table),
        "training_db": db_path.name if db_path else None,
        "training_db_sha256": sha256_of(db_path) if db_path else None,
        "n_rows": scaler["n_rows"],
        "n_papers": scaler["n_papers"],
        "idata_source": str(idata_path),
        "idata_sha256": sha256_of(idata_path),
        "posterior_total_draws": int(total),
        "posterior_kept_draws": int(len(idx)),
        "thinning_seed": seed,
        "convergence": {"max_r_hat": float(summ["r_hat"].max()), "min_ess_bulk": float(summ["ess_bulk"].min())},
        "pymc_version": str(idata.posterior.attrs.get("created_at", "")) and str(idata.attrs.get("pymc_version", idata.posterior.attrs.get("inference_library_version", ""))),
        "created_at": pd.Timestamp.now(tz="UTC").isoformat(),
    }
    write_manifest(out, extra, ["posterior.npz", "scaler.json"])
    return out


# ---------------------------------------------------------------------------
# GBM
# ---------------------------------------------------------------------------


def gbm_training_frame(table: Path) -> pd.DataFrame:
    df = pd.read_csv(table)
    if "system_type" in df.columns:
        df = df[df["system_type"] != "model_system"]
    df = df[df["age_d"] > 0].copy().reset_index(drop=True)
    df = derive_composition_features(df)
    for c in CAT_FEATURES:
        df[c] = df[c].fillna(MISSING)
    return df


def gbm_categories(df: pd.DataFrame) -> dict[str, list[str]]:
    cats = {}
    for c in CAT_FEATURES:
        vals = sorted(set(df[c].astype(str)))
        if MISSING not in vals:
            vals.append(MISSING)
        cats[c] = vals
    return cats


def _fit(X: pd.DataFrame, y: np.ndarray, seed: int, cat_idx: list[int]):
    import lightgbm as lgb

    m = lgb.LGBMRegressor(**{**TIGHT, "random_state": seed})
    m.fit(X.values, y, categorical_feature=cat_idx)
    return m


def gbm_lopo(df: pd.DataFrame, cats: dict[str, list[str]], *, seed: int = 42, n_folds: int = 10, target_scaling: str = "raw") -> tuple[float, float, pd.Series]:
    """multitask_v6.py:124-151 DoR-only path: paper-grouped 10 folds, seed-shuffled."""
    from sklearn.metrics import mean_absolute_error, r2_score

    papers = np.array(sorted(df["paper_doi"].unique()))
    rs = np.random.RandomState(seed)
    rs.shuffle(papers)
    folds = np.array_split(papers, n_folds)
    X = encode_frame(df, GBM_FEATURES, cats)
    cat_idx = [GBM_FEATURES.index(c) for c in CAT_FEATURES]
    y = df["dor_pct"].values.astype(float)
    oof = pd.Series(index=df.index, dtype=float)
    for f in folds:
        te = df["paper_doi"].isin(set(f)).values
        tr = ~te
        if te.sum() == 0:
            continue
        ytr = y[tr]
        if target_scaling == "z":
            mu, sd = ytr.mean(), ytr.std()
            m = _fit(X[tr], (ytr - mu) / sd, seed, cat_idx)
            oof[te] = m.predict(X[te].values) * sd + mu
        else:
            m = _fit(X[tr], ytr, seed, cat_idx)
            oof[te] = m.predict(X[te].values)
    oof = oof.clip(0, 100)
    ok = oof.notna()
    return float(r2_score(y[ok], oof[ok])), float(mean_absolute_error(y[ok], oof[ok])), oof


def gbm_oof_bayes_folds(df: pd.DataFrame, cats: dict[str, list[str]], table_for_folds: pd.DataFrame, *, seed: int = 42, target_scaling: str = "raw") -> pd.DataFrame:
    """GBM OOF on the *Bayesian* 5-fold split (bayes_hier_v4.py:135-138) for G1-3."""
    papers = np.array(sorted(table_for_folds["paper_doi"].unique()))
    rng = np.random.default_rng(seed)
    rng.shuffle(papers)
    X = encode_frame(df, GBM_FEATURES, cats)
    cat_idx = [GBM_FEATURES.index(c) for c in CAT_FEATURES]
    y = df["dor_pct"].values.astype(float)
    oof = pd.Series(index=df.index, dtype=float)
    for f in np.array_split(papers, 5):
        te = df["paper_doi"].isin(set(f)).values
        tr = ~te
        if te.sum() == 0:
            continue
        ytr = y[tr]
        if target_scaling == "z":
            mu, sd = ytr.mean(), ytr.std()
            m = _fit(X[tr], (ytr - mu) / sd, seed, cat_idx)
            oof[te] = m.predict(X[te].values) * sd + mu
        else:
            m = _fit(X[tr], ytr, seed, cat_idx)
            oof[te] = m.predict(X[te].values)
    return pd.DataFrame({"obs_uid": df["obs_uid"], "gbm_oof": oof.clip(0, 100)})


def export_gbm(table: Path, out: Path, *, seed: int = 42, lopo: bool = True, target_scaling: str = "raw", db_path: Path | None = None, oof_out: Path | None = None) -> dict[str, Any]:
    df = gbm_training_frame(table)
    cats = gbm_categories(df)
    X = encode_frame(df, GBM_FEATURES, cats)
    cat_idx = [GBM_FEATURES.index(c) for c in CAT_FEATURES]
    y = df["dor_pct"].values.astype(float)
    metrics: dict[str, Any] = {}
    if lopo:
        r2, mae, _ = gbm_lopo(df, cats, seed=seed, target_scaling=target_scaling)
        metrics = {"lopo_r2": r2, "lopo_mae_pct": mae, "lopo_folds": 10, "lopo_seed": seed}
        print(f"LOPO(10-fold, seed {seed}, target={target_scaling}): R2={r2:.3f} MAE={mae:.2f} %p")
    scale = {"mode": target_scaling}
    if target_scaling == "z":
        scale.update({"mean": float(y.mean()), "std": float(y.std())})
        m = _fit(X, (y - y.mean()) / y.std(), seed, cat_idx)
    else:
        m = _fit(X, y, seed, cat_idx)
    out.mkdir(parents=True, exist_ok=True)
    m.booster_.save_model(str(out / "model.txt"))
    meta = {
        "features": GBM_FEATURES,
        "categorical": CAT_FEATURES,
        "categories": cats,
        "target": "dor_pct",
        "target_scaling": scale,
        "hyperparameters": {**TIGHT, "random_state": seed},
        "n_rows": int(len(df)),
        "n_papers": int(df["paper_doi"].nunique()),
        "sigma_point_pct": 12.0,
        **metrics,
    }
    (out / "meta.json").write_text(json.dumps(meta, indent=2), encoding="utf-8")
    if oof_out is not None:
        bayes_tbl = bayes_training_frame(table.parent / "dor_scm_final.csv") if (table.parent / "dor_scm_final.csv").is_file() else df
        gbm_oof_bayes_folds(df, cats, bayes_tbl, seed=seed, target_scaling=target_scaling).to_csv(oof_out, index=False)
    import lightgbm as lgb

    write_manifest(
        out,
        {
            "model": "gbm_v6_dor_only",
            "training_table": table.name,
            "training_table_sha256": sha256_of(table),
            "training_db": db_path.name if db_path else None,
            "training_db_sha256": sha256_of(db_path) if db_path else None,
            "lightgbm_version": lgb.__version__,
            "created_at": pd.Timestamp.now(tz="UTC").isoformat(),
        },
        ["model.txt", "meta.json"],
    )
    return meta


# ---------------------------------------------------------------------------
# OOD reference (spec §5.5)
# ---------------------------------------------------------------------------

OOD_FEATURES = ["scm_pct", "w_b", "curing_temp_C", "CaO_SiO2", "Al_Si", "basicity"]


def export_ood_reference(table: Path, out_file: Path) -> dict[str, Any]:
    df = gbm_training_frame(table)
    ref: dict[str, Any] = {"features": OOD_FEATURES, "roles": {}}
    for role, g in df.groupby("scm_role"):
        sub = g[OOD_FEATURES].apply(pd.to_numeric, errors="coerce")
        med = sub.median()
        filled = sub.fillna(med)
        entry: dict[str, Any] = {
            "n_rows": int(len(g)),
            "n_papers": int(g["paper_doi"].nunique()),
            "n_mixes": int(g["mix_uid"].nunique()),
            "p01": {c: (None if np.isnan(v) else float(v)) for c, v in sub.quantile(0.01).items()},
            "p99": {c: (None if np.isnan(v) else float(v)) for c, v in sub.quantile(0.99).items()},
            "median": {c: (None if np.isnan(v) else float(v)) for c, v in med.items()},
        }
        ok_cols = [c for c in OOD_FEATURES if filled[c].notna().all() and filled[c].std() > 0]
        if len(g) >= 10 and len(ok_cols) >= 2:
            Z = filled[ok_cols].values.astype(float)
            mu = Z.mean(axis=0)
            cov = np.cov(Z, rowvar=False) + np.eye(len(ok_cols)) * 1e-6
            inv = np.linalg.pinv(cov)
            d2 = np.einsum("ij,jk,ik->i", Z - mu, inv, Z - mu)
            entry["mahalanobis"] = {
                "cols": ok_cols,
                "mean": mu.tolist(),
                "cov_inv": inv.tolist(),
                "d2_quantiles": {str(q): float(np.quantile(d2, q)) for q in (0.5, 0.75, 0.9, 0.95, 0.99)},
                "d2_sorted": np.sort(d2).tolist(),
            }
        ref["roles"][str(role)] = entry
    out_file.write_text(json.dumps(ref, indent=1), encoding="utf-8")
    return ref


# ---------------------------------------------------------------------------


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("what", choices=["bayes", "gbm", "ood", "all"])
    ap.add_argument("--idata")
    ap.add_argument("--table")
    ap.add_argument("--db")
    ap.add_argument("--out")
    ap.add_argument("--no-lopo", action="store_true")
    ap.add_argument("--target-scaling", default="raw", choices=["raw", "z"])
    ap.add_argument("--oof-out")
    a = ap.parse_args(argv)
    md = modeling_dir(required=False)
    root = Path(a.out) if a.out else bundles_dir()
    db = Path(a.db) if a.db else (md / "scm_dor_enriched.db" if md else None)
    if a.what in ("bayes", "all"):
        idata = Path(a.idata) if a.idata else md / "work" / "bayes_idata.nc"
        table = Path(a.table) if a.table else md / "dor_scm_final.csv"
        export_bayes(idata, table, root / "bayes_v4", db_path=db)
        print("bayes bundle ->", root / "bayes_v4")
    if a.what in ("gbm", "all"):
        table = Path(a.table) if a.table else md / "dor_scm_blended.csv"
        export_gbm(table, root / "gbm_v6", lopo=not a.no_lopo, target_scaling=a.target_scaling, db_path=db, oof_out=Path(a.oof_out) if a.oof_out else None)
        print("gbm bundle ->", root / "gbm_v6")
    if a.what in ("ood", "all"):
        table = Path(a.table) if a.table else md / "dor_scm_blended.csv"
        export_ood_reference(table, root / "ood_reference.json")
        print("ood reference ->", root / "ood_reference.json")
    return 0


if __name__ == "__main__":
    sys.exit(main())
