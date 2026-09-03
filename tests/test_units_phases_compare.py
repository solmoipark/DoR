"""Unit/basis harmonisation table cases, phase alias behaviour and the deterministic verdict."""

from __future__ import annotations

import numpy as np
import pytest

from dorgems.db.phases import UnknownPhaseError, group_of_db_phase, phase_mass_of_group, raw_names
from dorgems.db.units import harmonize
from dorgems.validate.compare import aggregate, compare_rows


@pytest.mark.parametrize(
    "unit,basis,mix,expected,grade",
    [
        ("g/100 g binder", None, {"scm_total_pct": 30, "w_b": 0.5}, 10.0, "A"),
        ("g/100 g cement", None, {"scm_total_pct": 30, "w_b": 0.5}, 7.0, "A"),
        ("g/100 g clinker", None, {"scm_total_pct": 30, "w_b": 0.5}, 7.0 * 0.95, "B"),
        ("g/100 g paste", None, {"scm_total_pct": 30, "w_b": 0.5}, 15.0, "C"),
        ("per_100g_ignited", None, {"scm_total_pct": 30, "w_b": 0.5}, None, "X"),
        ("g/100 g MK", None, {"scm_total_pct": 30, "w_b": 0.5}, 3.0, "B"),
        (None, "mass_percent_unspecified", {"scm_total_pct": 30, "w_b": 0.5}, 10.0, "D"),
        ("%", None, {}, 10.0, "D"),
    ],
)
def test_mass_basis_rules(unit, basis, mix, expected, grade):
    h = harmonize({"quantity": "CH_TGA", "value_norm": 10.0, "unit_norm": unit, "basis_reported": basis}, mix)
    assert h.grade == grade
    if expected is None:
        assert h.value is None
    else:
        assert abs(h.value - expected) < 1e-9
    assert h.usable == (grade in ("A", "B"))


def test_basis_precedence_over_normalised_unit():
    """G2-3 finding: the DB normalises unit_norm to 'g/100 g binder' even when the basis was
    unspecified; the explicit basis (or, when NULL, the reported unit text) must decide."""
    mix = {"scm_total_pct": 0, "w_b": 0.5}
    h = harmonize({"quantity": "CH_TGA", "value_norm": 9.4, "unit_norm": "g/100 g binder", "unit_reported": "wt.%", "basis_reported": "mass_percent_unspecified"}, mix)
    assert h.grade == "D"
    h2 = harmonize({"quantity": "CH_TGA", "value_norm": 17.2, "unit_norm": "g/100 g binder", "unit_reported": "%", "basis_reported": None}, mix)
    assert h2.grade == "D"
    h3 = harmonize({"quantity": "CH_TGA", "value_norm": 15.3, "unit_norm": "g/100 g binder", "unit_reported": "g/100 g binder", "basis_reported": None}, mix)
    assert h3.grade == "A"
    h4 = harmonize({"quantity": "CH_TGA", "value_norm": 15.9, "unit_norm": "g/100 g binder", "unit_reported": "g/g of cement", "basis_reported": "per g cement"}, {"scm_total_pct": 20, "w_b": 0.5})
    assert h4.grade == "A" and abs(h4.value - 15.9 * 0.8) < 1e-9


def test_other_quantities():
    assert harmonize({"quantity": "chem_shrink", "value_norm": 0.05, "unit_norm": "mL/g binder"}, {}).grade == "A"
    assert harmonize({"quantity": "chem_shrink", "value_norm": 5.0, "unit_norm": "%"}, {}).grade == "D"
    h = harmonize({"quantity": "DoR_SCM", "value_norm": 35.0, "unit_norm": "%"}, {})
    assert abs(h.value - 0.35) < 1e-12 and h.unit == "fraction"
    assert harmonize({"quantity": "QXRD_phase", "value_norm": 8.0, "unit_norm": "wt% (as reported)"}, {}).grade == "C"
    assert harmonize({"quantity": "cum_heat", "value_norm": 300.0, "unit_norm": "J/g"}, {}).grade == "X"
    # paste basis without w/b → excluded, never silently guessed
    assert harmonize({"quantity": "bound_water", "value_norm": 20.0, "unit_norm": "g/100 g paste"}, {}).grade == "X"


def test_phase_aliases():
    assert group_of_db_phase("Portlandite") == "portlandite"
    assert group_of_db_phase("CH") == "portlandite"
    assert group_of_db_phase("Alite") == "clinker:C3S"
    assert group_of_db_phase("quartz") == "ignore"
    assert group_of_db_phase("gypsum") == "sulfate_initial"
    assert group_of_db_phase("unobtainite") is None
    assert "Portlandite" in raw_names("portlandite")
    with pytest.raises(UnknownPhaseError):
        raw_names("unobtainite")
    row = {"phase_mass__Mock Portlandite": 0.01, "phase_mass__Mock C-S-H raw phase": 0.02}
    assert abs(phase_mass_of_group(row, "portlandite", mock=True) - 0.01) < 1e-12
    assert phase_mass_of_group({"phase_mass__X": 1.0}, "portlandite", mock=False, strict=False) is None


def test_verdict_rules():
    rng = np.random.default_rng(0)
    # bound_water is a primary quantity with no systematic offset
    pairs = [{"obs_uid": f"o{i}", "paper_doi": f"p{i % 7}", "mix_uid": f"m{i}", "quantity": "bound_water", "age_d": 28, "grade": "A", "obs_value": 10.0, "model_value": 10.0 + rng.normal(0, 1.0)} for i in range(20)]
    df = compare_rows(pairs)
    agg = aggregate(df)
    assert agg["verdict"]["bound_water"] == "consistent" and agg["overall"] == "consistent" and agg["overall_basis"] == "primary"
    pairs2 = [dict(p, model_value=p["obs_value"] + 12.0) for p in pairs]
    agg2 = aggregate(compare_rows(pairs2))
    assert agg2["verdict"]["bound_water"] == "tension"
    pairs3 = [dict(p, grade="D") for p in pairs]
    agg3 = aggregate(compare_rows(pairs3))
    assert agg3["verdict"]["bound_water"] == "insufficient_data" and agg3["n_usable"] == 0
    assert aggregate(compare_rows(pairs[:3]))["verdict"]["bound_water"] == "insufficient_data"
    # CH is secondary: its verdict is reported but does not decide the overall when a primary exists
    ch = [dict(p, quantity="CH_TGA", model_value=p["obs_value"] + 15.0) for p in pairs]
    agg4 = aggregate(compare_rows(pairs + ch))
    assert agg4["verdict"]["CH_TGA"] == "tension" and agg4["overall"] == "consistent" and agg4["secondary_verdicts"] == {"CH_TGA": "tension"}
    # the kernel offset (−3 g) is applied to CH model values by default
    ch_only = compare_rows([dict(ch[0])])
    assert abs(ch_only.iloc[0]["offset_b"] + 3.0) < 1e-9 and abs(ch_only.iloc[0]["r"] - 12.0) < 1e-9
