"""G0-2 (Bayes bundle golden), G0-3 (GBM bundle metrics), bundle integrity, ensemble properties."""

from __future__ import annotations

import json

import numpy as np
import pandas as pd
import pytest

from dorgems.models.bayes import predict_curve, role_kinetics_at_mean, weighted_quantile
from dorgems.models.bundle import BundleError, load_bayes_bundle, load_bundle
from dorgems.models.ensemble import fit_stretched_exp, reweight
from dorgems.models.gbm import predict_points
from dorgems.models.ood import assess


@pytest.fixture(scope="module")
def bundle(bundles_dir):
    return load_bundle(bundles_dir, require_gbm=False)


def test_bundle_hash_check_detects_edit(bundles_dir, tmp_path):
    import shutil

    dst = tmp_path / "bayes_v4"
    shutil.copytree(bundles_dir / "bayes_v4", dst)
    scaler = json.loads((dst / "scaler.json").read_text())
    scaler["beta_shape"] = 0.7
    (dst / "scaler.json").write_text(json.dumps(scaler))
    with pytest.raises(BundleError):
        load_bayes_bundle(dst)


def test_g0_2_role_kinetics_golden(bundle, modeling_dir):
    """Average-condition a_max/tau of the bundle vs the training run's bayes_role_kinetics.csv."""
    golden_path = bundle.bayes.path / "bayes_role_kinetics.golden.csv"
    if not golden_path.is_file():
        golden_path = modeling_dir / "work" / "bayes_role_kinetics.csv"
    gold = pd.read_csv(golden_path).set_index("role")
    got = role_kinetics_at_mean(bundle.bayes)
    for role, row in gold.iterrows():
        assert abs(got[role]["a_max_pct"] - row["a_max_pct"]) <= 0.5, (role, got[role], row.to_dict())
        assert abs(got[role]["tau_days"] / row["tau_days"] - 1.0) <= 0.02, (role, got[role], row.to_dict())


def test_predict_curve_properties(bundle):
    ages = np.array([1, 3, 7, 28, 90, 365], float)
    x = {"scm_pct": 40.0, "w_b": 0.45, "curing_temp_C": 20.0, "CaO_SiO2": 1.1}
    d = predict_curve(bundle.bayes, x, "slag", ages, new_study=True, method_group="selective_dissolution", rng_seed=1)
    assert d.alpha.shape == (bundle.bayes.n_draws, len(ages))
    assert np.all(d.alpha >= 0) and np.all(d.alpha <= 100)
    assert np.all(np.diff(d.alpha, axis=1) >= -1e-9), "curves must be monotone in age"
    s = d.summary()
    for q in ("q05", "q50", "q95"):
        assert len(s["alpha_pct"]["latent"][q]) == len(ages)
    lat = s["alpha_pct"]["latent"]
    assert all(a <= b <= c for a, b, c in zip(lat["q05"], lat["q50"], lat["q95"]))
    obs = s["alpha_pct"]["observed"]
    assert obs["q95"][3] - obs["q05"][3] >= lat["q95"][3] - lat["q05"][3] - 1e-9
    # pooled role flag
    d2 = predict_curve(bundle.bayes, x, "metakaolin", ages)
    assert any(f.startswith("role_pooled_as_other") for f in d2.flags)
    # determinism
    d3 = predict_curve(bundle.bayes, x, "slag", ages, new_study=True, method_group="selective_dissolution", rng_seed=1)
    assert np.allclose(d.alpha, d3.alpha)
    # imputation flagged
    d4 = predict_curve(bundle.bayes, {"scm_pct": 40.0, "w_b": 0.45}, "slag", ages)
    assert set(d4.imputed) == {"curing_temp_C", "CaO_SiO2"}


def test_weighted_quantile_matches_unweighted():
    v = np.random.default_rng(0).normal(size=5000)
    w = np.ones_like(v)
    assert np.allclose(weighted_quantile(v, [0.05, 0.5, 0.95], w), np.quantile(v, [0.05, 0.5, 0.95]), atol=0.02)


def test_g0_3_gbm_bundle_metrics(bundle):
    if bundle.gbm is None:
        pytest.skip("gbm bundle not exported")
    meta = bundle.gbm.meta
    assert meta["lopo_r2"] >= 0.50, meta
    ages = np.array([1, 3, 7, 28, 90, 365], float)
    x = {"scm_pct": 40.0, "w_b": 0.45, "curing_temp_C": 20.0, "scm_CaO": 41.8, "scm_SiO2": 36.5, "scm_Al2O3": 12.3, "scm_Fe2O3": 0.0, "scm_MgO": 7.5}
    p, info = predict_points(bundle.gbm, x, "slag", ages)
    assert p.shape == (len(ages),) and np.all(p >= 0) and np.all(p <= 100)
    assert info["method_group_used"] == "missing"
    p2, info2 = predict_points(bundle.gbm, x, "martian_dust", ages)
    assert info2["role_used"] == "missing" and info2["flags"]


def test_ensemble_reweight_and_fallback(bundle):
    ages = np.array([3, 7, 28, 90, 180], float)
    x = {"scm_pct": 40.0, "w_b": 0.45, "curing_temp_C": 20.0, "CaO_SiO2": 1.1}
    d = predict_curve(bundle.bayes, x, "slag", ages, rng_seed=0)
    # anchors at the posterior median → high ESS, applied
    med = np.median(d.alpha, axis=0)
    new, info = reweight(d, ages, med, sigma_g=12.0)
    assert info["applied"] and info["ess"] >= 100
    assert new.alpha.shape == d.alpha.shape
    # absurd anchors → fallback
    new2, info2 = reweight(d, ages, np.full(len(ages), 100.0), sigma_g=1.0)
    assert not info2["applied"] and info2["warnings"]
    assert new2 is d
    fit = fit_stretched_exp(ages, med)
    assert 0 < fit["a_max"] <= 100 and fit["tau_d"] > 0


def test_ood_assess(bundle):
    if not bundle.ood:
        pytest.skip("ood reference not exported")
    r = assess(bundle.ood, "slag", {"scm_pct": 40, "w_b": 0.45, "curing_temp_C": 20, "CaO_SiO2": 1.1, "Al_Si": 0.3, "basicity": 1.7})
    assert r["role_in_training"] and r["score_pct"] is not None
    r2 = assess(bundle.ood, "slag", {"scm_pct": 99, "w_b": 2.0, "curing_temp_C": 95, "CaO_SiO2": 1.1, "Al_Si": 0.3, "basicity": 1.7})
    assert any("above_p99" in f for f in r2["flags"])
    r3 = assess(bundle.ood, "unobtainium", {})
    assert r3["sparse_role"] and not r3["role_in_training"]
