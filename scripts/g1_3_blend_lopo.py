"""G1-3: decide the ensemble default by leave-papers-out comparison of
``bayes`` vs ``blend`` (GBM-anchored importance re-weighting).

Reproduces the 5-fold paper split of bayes_hier_v4.py (sorted papers, seeded
shuffle) and, per fold:
  * re-samples the v4 model on the training papers (same PyMC model, fewer
    draws than the reference run — only for the fold comparison), producing
    posterior draws exactly like the production bundle;
  * trains the GBM (spec settings) on the same training papers;
  * for every held-out observation: Bayesian predictive draws (new-study paper
    effect + method noise) → point/interval; GBM predictions at the anchor ages
    → importance weights (dorgems.models.ensemble.reweight) → blended
    point/interval.
Reports MAE, R² and 90 % coverage for both, writes work/g1_3_blend_lopo.csv/json.

Usage (needs the dorgems conda env, NUMBA pytensor, sitecustomize on PYTHONPATH):
    python scripts/g1_3_blend_lopo.py <modeling_dir> [--draws 600 --tune 1000]
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

SEED, BETA = 42, 0.5
ANCHORS = np.array([3.0, 7.0, 28.0, 90.0, 180.0])


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("modeling_dir")
    ap.add_argument("--draws", type=int, default=600)
    ap.add_argument("--tune", type=int, default=1000)
    ap.add_argument("--sigma-g", type=float, default=12.0)
    ap.add_argument("--ess-min", type=float, default=100.0)
    ap.add_argument("--out", default=None)
    a = ap.parse_args()

    import pymc as pm
    import pytensor.tensor as pt

    from dorgems.db.features import BAYES_FEATURES, CAT_FEATURES, GBM_FEATURES
    from dorgems.models.bayes import CurveDraws, stretched_exp_pct
    from dorgems.models.ensemble import reweight
    from dorgems.models.export_bundle import _fit, bayes_training_frame, gbm_categories, gbm_training_frame
    from dorgems.models.gbm import encode_frame

    md = Path(a.modeling_dir)
    out_dir = Path(a.out) if a.out else md / "work"
    df = bayes_training_frame(md / "dor_scm_final.csv")
    # --- identical preprocessing to bayes_hier_v4.py:55-69 ---
    X = df[BAYES_FEATURES].copy()
    for c in BAYES_FEATURES:
        X[c] = X[c].fillna(X[c].median())
        if X[c].std() > 0:
            X[c] = (X[c] - X[c].mean()) / X[c].std()
    Xv = X.fillna(0.0).values.astype("float64")
    methods = pd.Categorical(df["method_group"])
    _r = df["scm_role"].fillna("unknown")
    roles = pd.Categorical(_r.where(_r.isin(["slag", "fly_ash"]), "other"))
    role_idx = roles.codes
    t = df["age_d"].values.astype("float64")
    y = df["dor_pct"].values.astype("float64")
    n_f, n_m, n_r = Xv.shape[1], len(methods.categories), len(roles.categories)

    def build(Xv, t, y, p_idx, m_idx, r_idx, n_p):
        with pm.Model() as mdl:
            a0 = pm.Normal("a0_role", 0.0, 1.0, shape=n_r)
            beta_a = pm.Normal("beta_a", 0.0, 0.4, shape=n_f)
            sd_pa = pm.HalfNormal("sd_paper_amax", 0.8)
            z_pa = pm.Normal("z_paper_amax", 0.0, 1.0, shape=n_p)
            a_max = 100.0 * pm.math.sigmoid(a0[r_idx] + pt.dot(Xv, beta_a) + sd_pa * z_pa[p_idx])
            t0 = pm.Normal("t0_role", np.log(30.0), 0.8, shape=n_r)
            beta_t = pm.Normal("beta_t", 0.0, 0.4, shape=n_f)
            tau = pm.math.exp(t0[r_idx] + pt.dot(Xv, beta_t))
            alpha = a_max * (1.0 - pm.math.exp(-((t / tau) ** BETA)))
            sigma = pm.HalfNormal("sigma_method", 8.0, shape=n_m)
            pm.Normal("obs", mu=alpha, sigma=sigma[m_idx], observed=y)
        return mdl

    # --- GBM frame (blended, same rows minus model_system) ---
    gdf = gbm_training_frame(md / "dor_scm_blended.csv")
    cats = gbm_categories(gdf)
    Xg = encode_frame(gdf, GBM_FEATURES, cats)
    cat_idx = [GBM_FEATURES.index(c) for c in CAT_FEATURES]
    g_by_uid = {u: i for i, u in enumerate(gdf["obs_uid"])}

    rng = np.random.default_rng(SEED)
    pap = np.array(sorted(df["paper_doi"].unique()))
    rng.shuffle(pap)
    rows = []
    t_start = time.time()
    for k, f in enumerate(np.array_split(pap, 5)):
        te = df["paper_doi"].isin(set(f)).values
        tr = ~te
        p_tr = pd.Categorical(df.loc[tr, "paper_doi"])
        m = build(Xv[tr], t[tr], y[tr], p_tr.codes, methods.codes[tr], role_idx[tr], len(p_tr.categories))
        with m:
            it = pm.sample(a.draws, tune=a.tune, chains=2, cores=2, target_accept=0.93, random_seed=SEED + k, progressbar=False, idata_kwargs={"log_likelihood": False})
        q = it.posterior
        fa = q["beta_a"].values.reshape(-1, n_f)
        ft = q["beta_t"].values.reshape(-1, n_f)
        ra = q["a0_role"].values.reshape(-1, n_r)
        rt = q["t0_role"].values.reshape(-1, n_r)
        sdp = q["sd_paper_amax"].values.ravel()
        sig = q["sigma_method"].values.reshape(-1, n_m)
        S = len(fa)
        # GBM on the same training papers
        gtr = ~gdf["paper_doi"].isin(set(f)).values
        gm = _fit(Xg[gtr], gdf.loc[gtr, "dor_pct"].values.astype(float), SEED, cat_idx)
        # anchor-age frames for held-out rows
        te_idx = np.where(te)[0]
        for i in te_idx:
            uid = df.loc[i, "obs_uid"]
            gi = g_by_uid.get(uid)
            if gi is None:
                continue  # model_system rows are not in the GBM table
            xs = Xv[i]
            ri, mi = role_idx[i], methods.codes[i]
            u = rng.normal(0.0, sdp)
            eta = ra[:, ri] + fa @ xs + u
            am = 100.0 / (1.0 + np.exp(-eta))
            tau = np.exp(rt[:, ri] + ft @ xs)
            ages = np.array([t[i]])
            alpha = np.clip(stretched_exp_pct(ages, am, tau, BETA), 0, 100)
            noise = rng.normal(0.0, 1.0, size=(S, 1)) * sig[:, mi][:, None]
            draws_obs = np.clip(alpha + noise, 0, 100)[:, 0]
            b_mean, b_lo, b_hi = draws_obs.mean(), np.percentile(draws_obs, 5), np.percentile(draws_obs, 95)
            # GBM anchors for this mix
            base = Xg.iloc[[gi]].copy()
            frames = []
            for an in ANCHORS:
                r = base.copy()
                r["log_age"] = np.log10(an)
                frames.append(r)
            g_anch = np.clip(gm.predict(pd.concat(frames).values), 0, 100)
            g_pt = float(np.clip(gm.predict(base.values), 0, 100)[0])
            cd = CurveDraws(ages=ages, a_max=am, tau=tau, alpha=alpha, sigma_obs=sig[:, mi], role_used=str(roles.categories[ri]), imputed=[])
            new, info = reweight(cd, ANCHORS, g_anch, sigma_g=a.sigma_g, beta=BETA, ess_min=a.ess_min, rng_seed=SEED + i)
            if info["applied"]:
                noise2 = rng.normal(0.0, 1.0, size=(S, 1)) * new.sigma_obs[:, None]
                d2 = np.clip(new.alpha + noise2, 0, 100)[:, 0]
                e_mean, e_lo, e_hi = d2.mean(), np.percentile(d2, 5), np.percentile(d2, 95)
            else:
                e_mean, e_lo, e_hi = b_mean, b_lo, b_hi
            rows.append(dict(obs_uid=uid, fold=k, paper_doi=df.loc[i, "paper_doi"], age_d=t[i], role=str(roles.categories[ri]), method_group=str(methods.categories[mi]), observed=y[i], bayes_mean=b_mean, bayes_lo=b_lo, bayes_hi=b_hi, gbm_point=g_pt, blend_mean=e_mean, blend_lo=e_lo, blend_hi=e_hi, ess=info["ess"], applied=info["applied"]))
        print(f"fold {k+1}/5: {int(te.sum())} rows, {(time.time()-t_start)/60:.1f} min", flush=True)

    res = pd.DataFrame(rows)
    res.to_csv(out_dir / "g1_3_blend_lopo.csv", index=False)
    from sklearn.metrics import mean_absolute_error, r2_score

    def metrics(prefix: str) -> dict:
        m_ = res[f"{prefix}_mean"]
        return {"mae": float(mean_absolute_error(res["observed"], m_)), "r2": float(r2_score(res["observed"], m_)), "coverage90": float(np.mean((res["observed"] >= res[f"{prefix}_lo"]) & (res["observed"] <= res[f"{prefix}_hi"]))), "width90": float(np.mean(res[f"{prefix}_hi"] - res[f"{prefix}_lo"]))}

    summary = {"n": int(len(res)), "n_papers": int(res["paper_doi"].nunique()), "bayes": metrics("bayes"), "blend": metrics("blend"), "gbm_point_mae": float(mean_absolute_error(res["observed"], res["gbm_point"])), "gbm_point_r2": float(r2_score(res["observed"], res["gbm_point"])), "blend_applied_frac": float(res["applied"].mean()), "ess_median": float(res["ess"].median()), "sigma_g": a.sigma_g, "ess_min": a.ess_min, "draws": a.draws, "tune": a.tune}
    for role, g in res.groupby("role"):
        summary[f"role_{role}"] = {"n": int(len(g)), "bayes_mae": float(np.mean(np.abs(g["observed"] - g["bayes_mean"]))), "blend_mae": float(np.mean(np.abs(g["observed"] - g["blend_mean"]))), "bayes_cov": float(np.mean((g["observed"] >= g["bayes_lo"]) & (g["observed"] <= g["bayes_hi"]))), "blend_cov": float(np.mean((g["observed"] >= g["blend_lo"]) & (g["observed"] <= g["blend_hi"])))}
    decision = "blend" if (summary["blend"]["mae"] < summary["bayes"]["mae"] and summary["blend"]["coverage90"] >= 0.85) else "bayes"
    summary["decision"] = decision
    (out_dir / "g1_3_blend_lopo.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))
    print("DECISION:", decision)


if __name__ == "__main__":
    main()
