"""Environment variables and path resolution.

Environment
-----------
DORGEMS_DB            literature DB (read-only). Default: discovered (see :func:`literature_db_path`).
DORGEMS_STAGING_DB    staging DB for agent-produced values (read-write).
DORGEMS_BUNDLE        directory holding ``bayes_v4/`` and ``gbm_v6/`` bundles. Default: ``<repo>/bundles``.
DORGEMS_MODELING_DIR  the ``modeling/`` directory of the literature project (for export/golden tests).
INVERSE_GEMS_ROOT     InverseGems checkout (its ``configs/`` is used for materials/kinetics defaults).
DORGEMS_REAL_XGEMS    "1" enables the real-xGEMS test marker.
"""

from __future__ import annotations

import os
from pathlib import Path

PROJECT_DIRNAME = "DoR of SCMs in blended cements"
LITERATURE_DB_NAME = "scm_dor_enriched.db"


def repo_root() -> Path:
    return Path(__file__).resolve().parents[2]


def configs_dir() -> Path:
    return repo_root() / "configs"


def bundles_dir() -> Path:
    env = os.environ.get("DORGEMS_BUNDLE")
    return Path(env) if env else repo_root() / "bundles"


def _modeling_candidates() -> list[Path]:
    cands: list[Path] = []
    env = os.environ.get("DORGEMS_MODELING_DIR")
    if env:
        cands.append(Path(env))
    root = repo_root()
    for base in (root.parent, root.parent.parent, Path.home()):
        cands.append(base / PROJECT_DIRNAME / "modeling")
    return cands


def modeling_dir(required: bool = True) -> Path | None:
    for cand in _modeling_candidates():
        if cand.is_dir():
            return cand
    if required:
        raise FileNotFoundError(
            "Literature modeling directory not found. Set DORGEMS_MODELING_DIR to "
            f"'<...>/{PROJECT_DIRNAME}/modeling'."
        )
    return None


def data_dir() -> Path:
    return repo_root() / "data"


def literature_db_path(required: bool = True) -> Path | None:
    """DORGEMS_DB > <repo>/data/scm_dor_enriched.db (shipped with the repository) >
    the literature project's modeling/ directory."""
    env = os.environ.get("DORGEMS_DB")
    if env:
        return Path(env)
    shipped = data_dir() / LITERATURE_DB_NAME
    if shipped.is_file():
        return shipped
    md = modeling_dir(required=False)
    if md is not None and (md / LITERATURE_DB_NAME).is_file():
        return md / LITERATURE_DB_NAME
    if required:
        raise FileNotFoundError("Literature DB not found. Set DORGEMS_DB to the scm_dor_enriched.db path.")
    return None


def staging_db_path() -> Path:
    env = os.environ.get("DORGEMS_STAGING_DB")
    return Path(env) if env else repo_root() / "dorgems_staging.sqlite"


def inverse_gems_root(required: bool = False) -> Path | None:
    env = os.environ.get("INVERSE_GEMS_ROOT")
    if env:
        return Path(env)
    try:
        import inverse_gems  # type: ignore

        return Path(inverse_gems.__file__).resolve().parents[2]
    except Exception:  # noqa: BLE001
        pass
    cand = repo_root().parent / "InverseGems"
    if cand.is_dir():
        return cand
    if required:
        raise FileNotFoundError("InverseGems checkout not found. Set INVERSE_GEMS_ROOT.")
    return None


def real_xgems_enabled() -> bool:
    return os.environ.get("DORGEMS_REAL_XGEMS", "") == "1"


GEMS_SYSTEMS_DIRNAME = "gems_systems"


def dat_lst_path(required: bool = False) -> Path | None:
    """GEMS3K ``*-dat.lst`` for real runs: ``DORGEMS_DAT_LST`` or the first ``*-dat.lst``
    under ``<parent of repo>/gems_systems/*/`` (system files are never committed)."""
    env = os.environ.get("DORGEMS_DAT_LST")
    if env:
        return Path(env)
    base = repo_root().parent / GEMS_SYSTEMS_DIRNAME
    if base.is_dir():
        for cand in sorted(base.glob("*/*-dat.lst")):
            return cand
    if required:
        raise FileNotFoundError("no GEMS3K dat.lst found; set DORGEMS_DAT_LST")
    return None


def xgems_available() -> bool:
    try:
        import xgems  # type: ignore  # noqa: F401

        return True
    except Exception:  # noqa: BLE001
        return False
