from __future__ import annotations

import os
import shutil
from pathlib import Path

import pytest

TESTS = Path(__file__).resolve().parent
FIXTURE_DB = TESTS / "fixtures" / "mini_scm_dor.sqlite"
GOLDEN = TESTS / "golden"


@pytest.fixture(scope="session")
def fixture_db_path() -> Path:
    assert FIXTURE_DB.is_file(), "run scripts/make_fixture.py first"
    return FIXTURE_DB


@pytest.fixture()
def fixture_db_copy(tmp_path: Path, fixture_db_path: Path) -> Path:
    """A writable copy — used only to prove that open_ro refuses writes even on a writable file."""
    dst = tmp_path / "mini_copy.sqlite"
    shutil.copy(fixture_db_path, dst)
    return dst


@pytest.fixture(scope="session")
def full_db_path() -> Path:
    from dorgems import config

    p = config.literature_db_path(required=False)
    if p is None or not p.is_file():
        pytest.skip("full literature DB not available (set DORGEMS_DB)")
    return p


@pytest.fixture(scope="session")
def modeling_dir() -> Path:
    from dorgems import config

    md = config.modeling_dir(required=False)
    if md is None:
        pytest.skip("modeling dir not available (set DORGEMS_MODELING_DIR)")
    return md


@pytest.fixture(scope="session")
def bundles_dir() -> Path:
    from dorgems import config

    b = config.bundles_dir()
    if not (b / "bayes_v4" / "posterior.npz").is_file():
        pytest.skip("bayes_v4 bundle not exported yet")
    return b


def pytest_collection_modifyitems(config, items):  # noqa: ANN001
    real = os.environ.get("DORGEMS_REAL_XGEMS", "") == "1"
    for item in items:
        if "real_xgems" in item.keywords and not real:
            item.add_marker(pytest.mark.skip(reason="DORGEMS_REAL_XGEMS != 1"))
