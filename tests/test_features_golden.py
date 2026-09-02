"""G0-1: build_dor_table reproduces modeling/dor_scm_final.csv.

The fixture test runs always; the full-DB test needs the literature project.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from dorgems.db.features import (
    DOR_TABLE_COLUMNS,
    blended_only,
    build_aux_table,
    build_dor_table,
    derive_composition_features,
    method_group,
    role_from_text,
    scm_input_to_features,
)
from dorgems.db.reader import open_ro

from pathlib import Path

GOLDEN = Path(__file__).resolve().parent / "golden"

COMPARE_NUM = ["dor_pct", "scm_pct", "w_b", "curing_temp_C", "age_d"]
COMPARE_STR = ["scm_role", "system_type", "method_group", "resolve_how", "confidence", "mix_uid", "paper_doi"]


def _assert_matches_golden(df: pd.DataFrame, gold: pd.DataFrame) -> None:
    assert list(df.columns) == DOR_TABLE_COLUMNS == list(gold.columns)
    assert set(df["obs_uid"]) == set(gold["obs_uid"]), "obs_uid set differs"
    assert len(df) == len(gold)
    m = gold.merge(df, on="obs_uid", suffixes=("_g", "_d"))
    for c in COMPARE_NUM:
        a, b = m[c + "_g"], m[c + "_d"]
        assert (a.isna() == b.isna()).all(), c
        assert np.nanmax(np.abs(a - b).fillna(0)) < 1e-9, c
    for c in COMPARE_STR:
        a, b = m[c + "_g"].fillna("<NA>").astype(str), m[c + "_d"].fillna("<NA>").astype(str)
        assert (a == b).all(), c
    assert (m["fig_only_g"].fillna(-1) == m["fig_only_d"].fillna(-1)).all()


def test_golden_subset_on_fixture(fixture_db_path):
    con = open_ro(fixture_db_path)
    df = build_dor_table(con)
    gold = pd.read_csv(GOLDEN / "dor_scm_final_subset.csv")
    _assert_matches_golden(df, gold)
    # every resolver branch is exercised by the fixture
    assert {"json", "phase_match", "name_match", "common_role", "keyword_role", "unique_mat"} <= set(df["resolve_how"])
    assert "model_system" in set(df["system_type"])
    assert len(blended_only(df)) < len(df)


@pytest.mark.golden_full
def test_golden_full_db(full_db_path, modeling_dir):
    con = open_ro(full_db_path)
    df = build_dor_table(con)
    gold = pd.read_csv(modeling_dir / "dor_scm_final.csv")
    assert len(gold) == 1610
    _assert_matches_golden(df, gold)
    b = blended_only(df)
    assert (len(b), b["paper_doi"].nunique(), b["mix_uid"].nunique()) == (1592, 80, 476)


def test_method_group_and_role_text():
    assert method_group(None) == "unknown"
    assert method_group("Selective dissolution (EDTA)") == "selective_dissolution"
    assert method_group("SEM-BSE image analysis") == "SEM_BSE"
    assert method_group("XRD PONKCS") == "XRD_PONKCS"
    assert method_group("29Si NMR") == "NMR"
    assert method_group("Rietveld") == "XRD_other"
    assert method_group("mass balance") == "mass_balance"
    assert method_group("TGA") == "TGA"
    assert method_group("???") == "other_unknown"
    assert role_from_text("GGBS 40%") == "slag"
    assert role_from_text("class F fly ash") == "fly_ash"
    assert role_from_text("MK") == "metakaolin"
    assert role_from_text(None, "") is None


def test_aux_table_and_features(fixture_db_path):
    con = open_ro(fixture_db_path)
    aux = build_aux_table(con, ["CH_TGA", "bound_water", "QXRD_phase", "DoR_clinker"])
    assert set(aux["quantity"]) <= {"CH_TGA", "bound_water", "QXRD_phase", "DoR_clinker"}
    assert "source_locator" in aux.columns
    df = derive_composition_features(build_dor_table(con))
    for c in ["log_age", "CaO_SiO2", "basicity", "pozz_sum", "Al_Si", "amorph", "fineness"]:
        assert c in df.columns
    row = scm_input_to_features(
        {"role": "slag", "oxides": {"CaO": 40.0, "SiO2": 35.0, "Al2O3": 12.0, "MgO": 8.0}},
        {"scm_pct": 40.0, "w_b": 0.45, "curing_temp_C": 20.0},
        age_d=28,
    )
    assert abs(row["CaO_SiO2"] - 40 / 35) < 1e-12
    assert abs(row["basicity"] - (40 + 8 + 12) / 35) < 1e-12
    assert np.isnan(row["fineness"])
    assert abs(row["log_age"] - np.log10(28)) < 1e-12
