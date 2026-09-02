"""G1-2: mock forward run with an exported reaction model and a materials override.

input_reaction_degrees.json["scm"][slot] == exported alpha(age) (±1e-3) and
input_materials_used.json carries the override oxides.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest
import yaml

pytest.importorskip("inverse_gems")

from dorgems.envelope import build_forward_query  # noqa: E402
from dorgems.gems.forward import kernel_accepts_materials_config, materials_config_override, run_forward  # noqa: E402
from dorgems.kinetics.materials_override import build_materials_config, slot_for_role  # noqa: E402
from dorgems.kinetics.reaction_model import alpha_from_config, export_reaction_model, pin_reaction_model  # noqa: E402

SCM = {"name": "GGBS Pohang 2026", "role": "slag", "oxides": {"CaO": 42.5, "SiO2": 34.0, "Al2O3": 13.5, "MgO": 6.5, "Fe2O3": 0.5, "SO3": 2.0, "TiO2": 0.7}, "blaine_m2_kg": 450}
MIX = {"scm_pct": 40, "w_b": 0.45, "curing_temp_C": 20}


def _pred():
    return {"id": "mock1", "input": {"ages_d": [1, 7, 28, 90], "scm": SCM, "mix": MIX}, "beta_shape": 0.5, "bayes": {"a_max": {"q05": 45.0, "q50": 58.5, "q95": 68.0}, "tau_d": {"q05": 12.0, "q50": 18.5, "q95": 30.0}}, "recommended": {"source": "bayes"}, "provenance": {}}


def test_materials_override_and_slot(tmp_path):
    slot, warns, reactive = slot_for_role("slag", SCM["oxides"])
    assert slot == "slag" and reactive
    slot2, warns2, _ = slot_for_role("other", {"SiO2": 95, "Al2O3": 1, "CaO": 1})
    assert slot2 == "silica_fume" and warns2
    res = build_materials_config(SCM, tmp_path)
    y = yaml.safe_load(Path(res["path"]).read_text(encoding="utf-8"))
    ox = y["slag"]["oxide_mass_percent"]
    assert abs(sum(ox.values()) - 100) < 1e-3  # values are rounded to 4 decimals
    assert "GGBS_Pohang_2026" in y["slag"]["aliases"]
    assert set(y) >= {"OPC", "slag", "fly_ash", "metakaolin", "silica_fume", "limestone", "gypsum", "water"}
    with materials_config_override(res["path"]):
        from inverse_gems.materials import load_materials

        mats = load_materials()
        assert abs(mats["slag"].oxide_mass_percent["CaO"] - ox["CaO"]) < 1e-9
        assert "GGBS_Pohang_2026" in mats["slag"].aliases
    from inverse_gems.materials import load_materials as lm2

    assert abs(lm2()["slag"].oxide_mass_percent["CaO"] - 41.79) < 1e-6, "patch must be undone"


def test_g1_2_mock_forward_self_check(tmp_path):
    pred = _pred()
    mat = build_materials_config(SCM, tmp_path, slot="slag")
    rm = export_reaction_model(pred, tmp_path / "rm", slot="slag", config_id="mock1", signature_files=[mat["path"]])
    fq = build_forward_query(MIX, "slag", [1, 7, 28, 90], name="g1_2")
    assert fq["recipe"]["binders"] == {"OPC": 60.0, "slag": 40.0}
    res = run_forward(fq, out=tmp_path / "run", db=tmp_path / "ig.sqlite", reaction_model_config=rm["q50"]["path"], materials_config=mat["path"], slot="slag", use_mock=True)
    assert res.ok, (res.error, res.warnings, res.self_check)
    assert res.self_check["alpha_ok"] is True and res.self_check["materials_ok"] is True, res.self_check
    assert res.materials_injection in ("monkeypatch", "native")
    if not kernel_accepts_materials_config():
        assert res.materials_injection == "monkeypatch"
    ts = res.time_series
    assert ts is not None and len(ts) == 4
    assert any(c.startswith("phase_mass__") for c in ts.columns)
    # the kernel's alpha equals the exported curve at every age
    ages = np.array([1, 7, 28, 90], float)
    expected = alpha_from_config(rm["q50"]["path"], "slag", ages)
    alpha_details = [d for d in res.self_check["details"] if "alpha_in_kernel" in d]
    assert len(alpha_details) == 4
    for d in alpha_details:
        assert abs(d["alpha_in_kernel"] - expected[list(ages).index(d["age_days"])]) <= 1e-3
    assert any("materials_check" in d and d["materials_check"]["ok"] for d in res.self_check["details"])
    tr = res.to_tool_result()
    assert tr["contract"] == "inverse-gems-tool/1.0" and tr["ok"]
    manifest = json.loads((tmp_path / "run" / "dorgems_forward_manifest.json").read_text())
    assert manifest["materials_injection"] == res.materials_injection


def test_pin_config_gives_constant_alpha_in_kernel(tmp_path):
    p = pin_reaction_model(0.35, "fly_ash", tmp_path)
    fq = build_forward_query({"scm_pct": 30, "w_b": 0.5}, "fly_ash", [3, 28, 365])
    res = run_forward(fq, out=tmp_path / "run", db=tmp_path / "ig.sqlite", reaction_model_config=p, slot="fly_ash", use_mock=True)
    assert res.ok and res.self_check["alpha_ok"], res.self_check
    assert all(abs(d["alpha_in_kernel"] - 0.35) < 1e-9 for d in res.self_check["details"])


def test_real_run_requires_budget(tmp_path):
    fq = build_forward_query(MIX, "slag", [28])
    with pytest.raises(PermissionError):
        run_forward(fq, out=tmp_path / "r", db=tmp_path / "ig.sqlite", reaction_model_config=None, use_mock=False, max_xgems_calls=None)
    with pytest.raises(PermissionError):
        run_forward(fq, out=tmp_path / "r", db=tmp_path / "ig.sqlite", reaction_model_config=None, use_mock=False, max_xgems_calls=10_000)


def test_capture_path_mock(tmp_path):
    pred = _pred()
    rm = export_reaction_model(pred, tmp_path / "rm", slot="slag", config_id="cap")
    fq = build_forward_query(MIX, "slag", [7, 28])
    res = run_forward(fq, out=tmp_path / "run", db=tmp_path / "ig.sqlite", reaction_model_config=rm["q50"]["path"], slot="slag", use_mock=True, capture_species=True)
    assert res.ok, (res.error, res.warnings)
    assert res.time_series is not None and len(res.time_series) == 2
    caps = json.loads(Path(res.result_files["captures"]).read_text())
    assert len(caps) == 2 and "phase_masses" in caps[0]["capture"]
