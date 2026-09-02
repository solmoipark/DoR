"""G1-1 (logistic_fit deviation < 2 %p on 12 representative curves) and kinetics properties."""

from __future__ import annotations

import numpy as np
import pytest
import yaml

from dorgems.kinetics.curves import five_param_logistic, pin_params, stretched_exp
from dorgems.kinetics.fit import default_grid, logistic_fit
from dorgems.kinetics.reaction_model import alpha_from_config, compare_reaction_models, export_reaction_model, pin_reaction_model

# (role, a_max fraction, tau d): spans the posterior ranges of §1.1 plus OOD-ish extremes
REPRESENTATIVE = [
    ("slag", 0.585, 18.5), ("slag", 0.50, 8.0), ("slag", 0.66, 40.0), ("slag", 0.40, 60.0),
    ("fly_ash", 0.281, 27.8), ("fly_ash", 0.20, 10.0), ("fly_ash", 0.36, 80.0), ("fly_ash", 0.15, 150.0),
    ("other", 0.473, 13.1), ("other", 0.34, 3.0), ("other", 0.73, 30.0), ("other", 0.90, 5.0),
]


@pytest.mark.parametrize("role,a_max,tau", REPRESENTATIVE)
def test_g1_1_logistic_fit_deviation(role, a_max, tau):
    t = default_grid()
    y = stretched_exp(t, a_max, tau, 0.5)
    fit = logistic_fit(t, y)
    assert fit.max_abs_dev_pct < 2.0, (role, a_max, tau, fit.to_dict())
    assert fit.status == "ok"
    p = fit.params
    assert p["A"] == 0.0 and 0.05 <= p["D"] <= 1.0 and p["C"] > 0


def test_curve_properties():
    t = np.logspace(-1, 3, 50)
    a = stretched_exp(t, 0.6, 20.0)
    assert np.all(np.diff(a) >= 0) and a.min() >= 0 and a.max() <= 1
    pp = pin_params(0.35)
    assert np.allclose(five_param_logistic(t, **pp), 0.35)
    assert np.all(five_param_logistic(t, 0, 0.75, 20, 0.55, 1) <= 1.0)


def _fake_prediction():
    return {
        "id": "t1",
        "input": {"ages_d": [1, 3, 7, 28, 90, 365]},
        "role_bayes": "slag",
        "role_gbm": "slag",
        "beta_shape": 0.5,
        "bayes": {"a_max": {"q05": 40.0, "q50": 58.5, "q95": 70.0}, "tau_d": {"q05": 10.0, "q50": 18.5, "q95": 35.0}},
        "recommended": {"source": "bayes"},
        "provenance": {"dorgems": "test"},
    }


def test_export_modes_and_alpha_roundtrip(tmp_path):
    pred = _fake_prediction()
    res = export_reaction_model(pred, tmp_path, mode="logistic_fit", slot="slag", config_id="t1")
    assert set(res) == {"q05", "q50", "q95"}
    y = yaml.safe_load(open(res["q50"]["path"], encoding="utf-8"))
    assert y["id"] == "dorgems_t1_q50" and y["availability_modifier"] == {"enabled": False}
    assert set(y["scm_reaction"]["slag"]) == {"A", "B", "C", "D", "G"}
    assert "provenance" not in y and (tmp_path / "dorgems_t1_q50.provenance.json").is_file()
    ages = np.array([1, 3, 7, 28, 90, 365], float)
    a = alpha_from_config(res["q50"]["path"], "slag", ages)
    ref = stretched_exp(ages, 0.585, 18.5, 0.5)
    assert np.max(np.abs(a - ref)) * 100 < 2.0
    # quantile ordering of exported curves
    a05 = alpha_from_config(res["q05"]["path"], "slag", ages)
    a95 = alpha_from_config(res["q95"]["path"], "slag", ages)
    assert np.all(a05 <= a + 0.02) and np.all(a95 >= a - 0.02)
    # native
    resn = export_reaction_model(pred, tmp_path / "native", mode="native", slot="slag", config_id="t1")
    yn = yaml.safe_load(open(resn["q50"]["path"], encoding="utf-8"))
    assert yn["scm_reaction"]["slag"]["model"] == "dorgems_stretched_exp"
    an = alpha_from_config(resn["q50"]["path"], "slag", ages)
    assert np.allclose(an, ref, atol=1e-9)
    # pin
    p = pin_reaction_model(0.35, "fly_ash", tmp_path / "pin")
    assert np.allclose(alpha_from_config(p, "fly_ash", ages), 0.35)
    cmp = compare_reaction_models(res["q50"]["path"], resn["q50"]["path"], "slag", ages)
    assert cmp["max_abs_diff_pct"] < 2.0
    with pytest.raises(ValueError):
        export_reaction_model(pred, tmp_path, slot="limestone")
