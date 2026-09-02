"""Scenario B on the fixture DB with the mock kernel: recipe reconstruction, twin
comparison artefacts, OPC reference pipeline, and the read-only guarantee."""

from __future__ import annotations

import sqlite3
from pathlib import Path

import pandas as pd
import pytest

pytest.importorskip("inverse_gems")

from dorgems.db.reader import LiteratureDB  # noqa: E402
from dorgems.validate.twin import db_mix_to_recipe, opc_reference_candidates, twin_compare_mix  # noqa: E402
from dorgems.validate.twin_batch import candidate_mixes  # noqa: E402


def _a_mix_with_obs(db: LiteratureDB) -> str | None:
    rows = db.con.execute("SELECT mix_uid, COUNT(*) n FROM observations WHERE quantity IN ('CH_TGA','bound_water') AND value_norm IS NOT NULL AND age_d > 0 GROUP BY mix_uid ORDER BY n DESC").fetchall()
    for r in rows:
        m = db.mix(r[0])
        if m and m.get("w_b") is not None and m.get("curing_temp_C") is not None:
            return r[0]
    return None


def test_db_mix_to_recipe(fixture_db_path, tmp_path):
    with LiteratureDB(fixture_db_path) as db:
        for m in db.mixes_for_paper("10.1016/j.cemconcomp.2012.02.004")[:5]:
            rec = db_mix_to_recipe(m, db.materials_for_paper(m["paper_doi"]), ages=[28.0], out_dir=tmp_path)
            if rec.get("excluded_reason"):
                continue
            b = rec["forward_query"]["recipe"]["binders"]
            assert abs(sum(b.values()) - 100) < 1e-6
            assert set(b) <= {"OPC", "slag", "fly_ash", "metakaolin", "silica_fume", "limestone", "gypsum"}
            break
        else:
            pytest.skip("no reconstructible mix in fixture paper")


def test_twin_compare_mock(fixture_db_path, tmp_path, bundles_dir):
    with LiteratureDB(fixture_db_path) as db:
        mix_uid = _a_mix_with_obs(db)
        if mix_uid is None:
            pytest.skip("fixture has no mix with CH/BW observations and w/b + T")
        res = twin_compare_mix(db, mix_uid, out=tmp_path / "twin", ig_db=tmp_path / "igdb", use_mock=True)
    if not res["ok"]:
        pytest.skip(f"mix not reconstructible: {res.get('error')}")
    assert Path(res["files"]["comparison_csv"]).is_file()
    df = pd.read_csv(res["files"]["comparison_csv"])
    assert {"obs_uid", "quantity", "age_d", "obs", "model", "r", "z", "grade", "assumptions"} <= set(df.columns)
    assert res["aggregate"]["overall"] in ("consistent", "tension", "insufficient_data")
    # literature DB untouched
    con = sqlite3.connect(fixture_db_path)
    assert con.execute("SELECT COUNT(*) FROM observations").fetchone()[0] == 283


def test_candidates_and_opc_reference(fixture_db_path, full_db_path):
    with LiteratureDB(full_db_path) as db:
        c = candidate_mixes(db, min_dor_ages=3, min_common_ages=3)
        assert len(c) >= 50, len(c)  # spec §1.1: 66 mixes with ≥3 common ages (CH or BW)
        opc = opc_reference_candidates(db, age_days=28, w_b_range=(0.4, 0.5))
        assert 0 < opc["mix_uid"].nunique() <= 393
        assert set(opc["grade"]) <= {"A", "B", "C", "D", "X"}
