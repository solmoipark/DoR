"""Real-xGEMS gates (DORGEMS_REAL_XGEMS=1, xgems importable, a GEMS3K dat.lst present).

G1-4  default-kinetics OPC60/slag40 and OPC100 runs; values recorded as anchors for THIS
      system (the GemsPilot anchors were produced on a different system with CNASH).
G2-1  phase names of xgems_phase_amounts_raw.csv confirm configs/phase_aliases.yaml.
G2-2  mass/volume units and chemical shrinkage of OPC paste at 28 d in the literature range.
G1-2r self-check of an exported DoRGems reaction model on the real kernel.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from dorgems import config
from dorgems.db import phases as P
from dorgems.envelope import build_forward_query
from dorgems.gems.forward import run_forward
from dorgems.gems.observables import chem_shrink_ml_per_g, mass_factor, observables_from_run
from dorgems.kinetics.materials_override import build_materials_config
from dorgems.kinetics.reaction_model import export_reaction_model

pytestmark = pytest.mark.real_xgems

SCM = {"name": "GGBS-real-test", "role": "slag", "oxides": {"CaO": 42.5, "SiO2": 34.0, "Al2O3": 13.5, "MgO": 6.5, "Fe2O3": 0.5, "SO3": 2.0}}


@pytest.fixture(scope="module")
def dat_lst() -> Path:
    if not config.xgems_available():
        pytest.skip("xgems not importable")
    p = config.dat_lst_path(required=False)
    if p is None or not p.is_file():
        pytest.skip("no dat.lst (set DORGEMS_DAT_LST)")
    return p


@pytest.fixture(scope="module")
def anchors_dir(tmp_path_factory) -> Path:
    return tmp_path_factory.mktemp("real")


def _run(fq, out, db, dat_lst, **kw):
    res = run_forward(fq, out=out, db=db, reaction_model_config=kw.pop("reaction_model_config", None), use_mock=False, dat_lst=dat_lst, max_xgems_calls=kw.pop("max_xgems_calls", 10), **kw)
    assert res.ok, (res.error, res.warnings, res.self_check)
    return res


def test_g1_4_default_kinetics_anchors(dat_lst, anchors_dir):
    db = anchors_dir / "igdb"
    fq1 = build_forward_query({"scm_pct": 0, "w_b": 0.5}, "slag", [28.0], name="opc100")
    r1 = _run(fq1, anchors_dir / "opc100", db, dat_lst)
    fq2 = build_forward_query({"scm_pct": 40, "w_b": 0.45}, "slag", [28.0], name="opc60_slag40")
    r2 = _run(fq2, anchors_dir / "opc60_slag40", db, dat_lst)
    rows = {"opc100_wb0.5_28d": r1.time_series.iloc[0].to_dict(), "opc60_slag40_wb0.45_28d": r2.time_series.iloc[0].to_dict()}
    anchors = {}
    for k, row in rows.items():
        anchors[k] = {c: row[c] for c in row if c.startswith(("porosity", "scalar__porosity", "scalar__pH", "scalar__system_mass", "scalar__system_volume", "phase_mass__Portlandite", "phase_mass__CSHQ", "phase_mass__CNASH"))}
        por = row.get("scalar__porosity", row.get("porosity"))
        assert por is not None and 0.0 < float(por) < 1.0, anchors[k]
    out = config.repo_root() / "docs" / "real_anchors_TINN_v4.json"
    out.write_text(json.dumps({"dat_lst": str(dat_lst), "anchors": anchors}, indent=2, default=str), encoding="utf-8")
    assert (anchors_dir / "opc100" / "forward" / "time_series.csv").is_file()


def test_g2_1_phase_names_confirm_alias_table(dat_lst, anchors_dir):
    db = anchors_dir / "igdb"
    fq = build_forward_query({"scm_pct": 30, "w_b": 0.5, "other_components": {"limestone": 5}}, "slag", [28.0], name="phases")
    res = _run(fq, anchors_dir / "phases", db, dat_lst)
    raw = sorted(Path(db).rglob("xgems_phase_amounts_raw.csv"))
    assert raw, "no raw phase file written"
    names = list(pd.read_csv(raw[-1]).iloc[:, 0].astype(str))
    report = P.confirm_from_raw_names(names, str(raw[-1]))
    (config.repo_root() / "docs" / "g2_1_phase_confirmation.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
    missing = [g for g, e in report["groups"].items() if not e["ok"] and g not in ("hydrogarnet", "brucite")]
    assert not missing, (missing, report["unmatched_raw"])
    assert "Portlandite" in names


def test_g2_2_units_and_chem_shrinkage(dat_lst, anchors_dir):
    db = anchors_dir / "igdb"
    fq = build_forward_query({"scm_pct": 0, "w_b": 0.5}, "slag", [28.0], name="opc_cs")
    res = _run(fq, anchors_dir / "opc_cs", db, dat_lst, capture_species=True)
    row = res.time_series.iloc[0].to_dict()
    f = mass_factor(row)
    ch = float(row["phase_mass__Portlandite"]) * f
    # 28-d OPC paste w/b 0.5: portlandite ~ 15–28 g/100 g binder in the literature
    assert 10.0 < ch < 35.0, (ch, f, row.get("scalar__system_mass"))
    obs = observables_from_run(res.forward_dir)
    caps = json.loads((res.forward_dir / "dorgems_captures.json").read_text(encoding="utf-8"))
    pj = sorted(Path(caps[0]["recipe_dir"] or caps[0]["chemistry_dir"]).rglob("porosity.json"))
    porosity = json.loads(pj[-1].read_text(encoding="utf-8")) if pj else None
    cs_cm3 = chem_shrink_ml_per_g(porosity, volume_unit="cm3")["value"]
    cs_m3 = chem_shrink_ml_per_g(porosity, volume_unit="m3")["value"]
    result = {"CH_g_per_100g": ch, "mass_factor": f, "chem_shrink_if_cm3": cs_cm3, "chem_shrink_if_m3": cs_m3, "bound_water_g": obs.iloc[0]["bound_water_g"], "porosity_keys": list(porosity) if porosity else None}
    (config.repo_root() / "docs" / "g2_2_units_TINN_v4.json").write_text(json.dumps(result, indent=2, default=str), encoding="utf-8")
    ok_unit = [u for u, v in (("cm3", cs_cm3), ("m3", cs_m3)) if v is not None and 0.02 <= v <= 0.10]
    assert ok_unit, result  # literature 0.04–0.07 mL/g binder at 28 d


def test_g1_2_real_self_check(dat_lst, anchors_dir):
    db = anchors_dir / "igdb"
    pred = {"id": "real1", "input": {"ages_d": [7, 28]}, "beta_shape": 0.5, "bayes": {"a_max": {"q50": 56.9}, "tau_d": {"q50": 13.3}}, "recommended": {"source": "bayes"}, "provenance": {}}
    mat = build_materials_config(SCM, anchors_dir, slot="slag")
    rm = export_reaction_model(pred, anchors_dir / "rm", slot="slag", quantiles=(0.5,), config_id="real1", signature_files=[mat["path"]])
    fq = build_forward_query({"scm_pct": 40, "w_b": 0.45}, "slag", [7.0, 28.0], name="real_selfcheck")
    res = _run(fq, anchors_dir / "real_selfcheck", db, dat_lst, reaction_model_config=rm["q50"]["path"], materials_config=mat["path"], slot="slag")
    assert res.self_check["alpha_ok"] is True and res.self_check["materials_ok"] is True, res.self_check
    assert np.all(res.time_series["age_days"].values == [7.0, 28.0])
