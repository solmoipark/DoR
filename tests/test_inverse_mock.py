"""G3-1 (mock): synthetic recovery of alpha(28 d) from CH/bound-water observations,
staging dry-run vs write, and the scenario-C tool contract.

The mock runner's CH is proportional to reacted solids (it *increases* with alpha,
the reverse of a real pozzolanic system) — the monotonicity flag must fire, but
the recovery logic is agnostic to the sign of F(alpha).
"""

from __future__ import annotations

import json
import sqlite3
from pathlib import Path

import numpy as np
import pytest

pytest.importorskip("inverse_gems")

from dorgems.envelope import build_forward_query  # noqa: E402
from dorgems.inverse.alpha_grid import build_forward_map, default_alpha_grid  # noqa: E402
from dorgems.inverse.likelihood import Likelihood, build_points  # noqa: E402
from dorgems.inverse.posterior import infer, summarise  # noqa: E402
from dorgems.inverse.staging import open_staging, review_approve, review_list, stage_inference  # noqa: E402
from dorgems.kinetics.curves import stretched_exp  # noqa: E402

SCM = {"name": "FA-test", "role": "fly_ash", "oxides": {"CaO": 5.0, "SiO2": 55.0, "Al2O3": 25.0, "Fe2O3": 6.0, "MgO": 1.5}}
CASES = [(0.30, 25.0), (0.45, 12.0), (0.20, 60.0), (0.60, 8.0), (0.35, 30.0)]


@pytest.fixture(scope="module")
def fmap(tmp_path_factory, bundles_dir):
    tmp = tmp_path_factory.mktemp("fmap")
    fq = build_forward_query({"scm_pct": 30, "w_b": 0.5}, "fly_ash", [7, 28, 90])
    return build_forward_map(fq, slot="fly_ash", ages=[7.0, 28.0, 90.0], alphas=default_alpha_grid(11), out=tmp, ig_db=tmp / "igdb", materials_config=None, use_mock=True, quantities=("CH_TGA", "bound_water"))


def test_forward_map_shape_and_flags(fmap):
    assert fmap.table["CH_g"].shape == (3, 11)
    assert np.all(np.isfinite(fmap.table["CH_g"]))
    rep = fmap.monotonicity_report()
    assert rep["CH_g"], "mock CH increases with alpha → flag expected"
    v = fmap.value("CH_TGA", 1, np.array([0.0, 0.5, 1.0]))
    assert v.shape == (3,)


@pytest.mark.parametrize("a_true,tau_true", CASES)
def test_g3_1_synthetic_recovery(fmap, a_true, tau_true, bundles_dir):
    rng = np.random.default_rng(int(a_true * 1000))
    ages = [7.0, 28.0, 90.0]
    obs = []
    for i, t in enumerate(ages):
        a = float(stretched_exp(np.array([t]), a_true, tau_true, 0.5)[0])
        for q, sig in (("CH_TGA", 0.3), ("bound_water", 0.3)):
            y = float(fmap.value(q, i, a)) + rng.normal(0, sig)
            obs.append({"obs_uid": f"{q}@{t}", "age_d": t, "quantity": q, "value": y, "grade": "A", "uncertainty": sig})
    pts, skipped = build_points(obs, ages, sigma_obs_default={}, sigma_model={"CH_TGA": 0.2, "bound_water": 0.2})
    assert len(pts) == 6 and not skipped
    lik = Likelihood(fmap, pts, beta=0.5)
    post = infer(lik, prior_a_max=None, prior_tau=None, prior="flat", grid_n=40, rng_seed=1)
    summ = summarise(post, lik, np.array([7.0, 28.0, 90.0]))
    a28_true = float(stretched_exp(np.array([28.0]), a_true, tau_true, 0.5)[0])
    q05, q50, q95 = summ["alpha"]["q05"][1], summ["alpha"]["q50"][1], summ["alpha"]["q95"][1]
    assert q05 - 0.02 <= a28_true <= q95 + 0.02, (a28_true, q05, q50, q95, post["method"])
    assert abs(q50 - a28_true) <= 0.05, (a28_true, q50)
    assert summ["n_observations_used"] == 6 and len(summ["ppc"]) == 6


def test_kl_increases_with_observations(fmap, bundles_dir):
    """G3-4: information gain grows with the number of observations (model prior)."""
    from dorgems.models.bundle import load_bundle
    from dorgems.models.bayes import predict_curve

    b = load_bundle(bundles_dir, require_gbm=False)
    d = predict_curve(b.bayes, {"scm_pct": 30, "w_b": 0.5, "curing_temp_C": 20, "CaO_SiO2": 0.09}, "fly_ash", np.array([28.0]), rng_seed=0)
    ages = [7.0, 28.0, 90.0]
    a_true, tau_true = 0.30, 25.0
    obs = []
    for i, t in enumerate(ages):
        a = float(stretched_exp(np.array([t]), a_true, tau_true, 0.5)[0])
        obs.append({"obs_uid": f"CH@{t}", "age_d": t, "quantity": "CH_TGA", "value": float(fmap.value("CH_TGA", i, a)), "grade": "A", "uncertainty": 0.3})
    kls = []
    for n in (1, 2, 3):
        pts, _ = build_points(obs[:n], ages, sigma_obs_default={}, sigma_model={"CH_TGA": 0.2})
        post = infer(Likelihood(fmap, pts), prior_a_max=d.a_max / 100.0, prior_tau=d.tau, prior="model", ess_min=1.0, rng_seed=0)
        kls.append(post["prior_vs_posterior_kl"])
    assert kls[0] <= kls[1] + 1e-9 <= kls[2] + 2e-9, kls


def test_staging_dry_run_and_write(tmp_path):
    inf = {"id": "abc123", "alpha": {"ages_d": [7, 28], "q05": [0.1, 0.2], "q50": [0.2, 0.3], "q95": [0.3, 0.4]}, "a_max": {"q50": 0.5}, "tau_d": {"q50": 20.0}, "ess": 120.0, "posterior_method": "importance", "ppc": [{"label": "x"}], "slot": "fly_ash", "provenance": {"bundle_bayes": "sha256:1"}}
    db = tmp_path / "staging.sqlite"
    dry = stage_inference(inf, path=db, dry_run=True)
    assert dry["dry_run"] and dry["n_rows"] == 2
    con = open_staging(db)
    assert con.execute("SELECT COUNT(*) FROM inferred_dor").fetchone()[0] == 0
    con.close()
    w = stage_inference(inf, path=db, dry_run=False, note="test")
    assert not w["dry_run"]
    rows = review_list(db)
    assert rows[0]["inference_id"] == "abc123" and rows[0]["reviewed"] == 0
    review_approve(db, "abc123", note="ok")
    assert review_list(db, reviewed=1)[0]["inference_id"] == "abc123"
    con = sqlite3.connect(db)
    assert con.execute("SELECT DISTINCT reviewed FROM inferred_dor").fetchall() == [(1,)]


def test_infer_tool_end_to_end_mock(tmp_path, bundles_dir):
    from dorgems.pilot.tools_b_c import dor_infer_from_observations, dor_stage_inferred

    obs = [{"age_d": 7, "quantity": "CH_TGA", "value": 4.0, "unit": "g/100 g binder"}, {"age_d": 28, "quantity": "CH_TGA", "value": 5.0, "unit": "g/100 g binder"}, {"age_d": 28, "quantity": "bound_water", "value": 12.0, "unit": "g/100 g paste"}, {"age_d": 28, "quantity": "DoR_SCM", "value": 30.0, "unit": "%"}]
    r = dor_infer_from_observations({"scm_pct": 30, "w_b": 0.5}, obs, str(tmp_path / "inf"), str(tmp_path / "igdb"), scm=SCM, alpha_grid=6, use_mock=True)
    assert r["ok"], r
    s = r["summary"]
    assert s["n_observations_used"] == 2  # bound water in paste basis → grade C excluded; DoR kept for validation only
    assert s["validation"] and s["validation"][0]["dor_measured"] == pytest.approx(0.30)
    inf_path = r["artifacts"]["inference"]
    payload = json.loads(Path(inf_path).read_text(encoding="utf-8"))
    assert payload["schema"] == "dorgems-inference/1.0"
    assert Path(r["artifacts"]["inferred_dor_csv"]).read_text().startswith("scm,age_d,dor")
    st = dor_stage_inferred(inf_path, str(tmp_path / "staging.sqlite"), use_mock=True)
    assert st["ok"] and st["summary"]["dry_run"]
    st2 = dor_stage_inferred(inf_path, str(tmp_path / "staging.sqlite"), use_mock=False)
    assert st2["ok"] and not st2["summary"]["dry_run"]
