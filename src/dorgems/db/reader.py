"""Read-only access to the literature DB.

Rules (spec §4.1, §11):

* the DB is opened with ``?mode=ro`` only. There is no code path that returns a
  writable connection to the literature DB;
* no free SQL is exposed. Only named, parameterised queries;
* every row carries provenance columns (``paper_doi``, ``source_locator``,
  ``extraction_confidence``, ``fig_only``) where the table has them.
"""

from __future__ import annotations

import sqlite3
from collections.abc import Iterable, Sequence
from pathlib import Path
from typing import Any

MODEL_SYSTEM = "model_system"

QUANTITIES = (
    "DoR_SCM",
    "DoR_clinker",
    "CH_TGA",
    "CH_XRD",
    "bound_water",
    "QXRD_phase",
    "chem_shrink",
    "cum_heat",
)


def _ro_uri(path: str | Path) -> str:
    p = Path(path).resolve()
    if not p.is_file():
        raise FileNotFoundError(f"Literature DB not found: {p}")
    return p.as_uri() + "?mode=ro"


def open_ro(path: str | Path) -> sqlite3.Connection:
    """Open the literature DB read-only. Any write raises ``sqlite3.OperationalError``."""
    con = sqlite3.connect(_ro_uri(path), uri=True)
    con.row_factory = sqlite3.Row
    # Belt and braces: even if the file were writable, refuse writes at the
    # connection level too.
    con.execute("PRAGMA query_only = 1")
    return con


def rows_to_dicts(rows: Iterable[sqlite3.Row]) -> list[dict[str, Any]]:
    return [dict(r) for r in rows]


class LiteratureDB:
    """Named-query facade over a read-only connection."""

    def __init__(self, path: str | Path):
        self.path = Path(path).resolve()
        self.con = open_ro(self.path)

    # -- context management -------------------------------------------------
    def close(self) -> None:
        self.con.close()

    def __enter__(self) -> "LiteratureDB":
        return self

    def __exit__(self, *exc: Any) -> None:
        self.close()

    # -- helpers ------------------------------------------------------------
    def _q(self, sql: str, params: Sequence[Any] | dict[str, Any] = ()) -> list[dict[str, Any]]:
        return rows_to_dicts(self.con.execute(sql, params))

    def counts(self) -> dict[str, int]:
        out = {}
        for t in ("papers", "materials", "mixes", "observations"):
            out[t] = int(self.con.execute(f"SELECT COUNT(*) FROM {t}").fetchone()[0])
        return out

    # -- named queries (spec §4.1) -----------------------------------------
    def paper(self, doi: str) -> dict[str, Any] | None:
        rows = self._q("SELECT * FROM papers WHERE doi = ?", (doi,))
        return rows[0] if rows else None

    def mix(self, mix_uid: str) -> dict[str, Any] | None:
        rows = self._q("SELECT * FROM mixes WHERE mix_uid = ?", (mix_uid,))
        return rows[0] if rows else None

    def mixes_for_paper(self, doi: str) -> list[dict[str, Any]]:
        return self._q("SELECT * FROM mixes WHERE paper_doi = ? ORDER BY mix_uid", (doi,))

    def materials_for_paper(self, doi: str) -> list[dict[str, Any]]:
        return self._q("SELECT * FROM materials WHERE paper_doi = ? ORDER BY material_uid", (doi,))

    def materials_by_paper(self) -> dict[str, dict[str, dict[str, Any]]]:
        """{paper_doi: {material_id: material_row}} — same shape as build_clean.mats_by_paper."""
        out: dict[str, dict[str, dict[str, Any]]] = {}
        for r in self.con.execute("SELECT * FROM materials"):
            out.setdefault(r["paper_doi"], {})[r["material_id"]] = dict(r)
        return out

    def dor_observations(
        self,
        role: str | None = None,
        paper: str | None = None,
        *,
        include_model_system: bool = True,
    ) -> list[dict[str, Any]]:
        """DoR_SCM observations joined with their mix (raw; role filtering is on
        ``mixes.primary_scm_role`` because the resolved role needs the feature
        builder — use :mod:`dorgems.db.features` for the modelling table)."""
        sql = (
            "SELECT o.obs_uid, o.paper_doi, o.mix_uid, o.age_d, o.quantity, o.phase_name, "
            "o.value_norm, o.unit_norm, o.basis_reported, o.method, o.method_detail, o.uncertainty, "
            "o.source_locator, o.fig_only, o.extraction_confidence, o.sanity_ok, o.reviewed, "
            "m.curing_type, m.primary_scm_role, m.scm_total_pct, m.w_b, m.curing_temp_C "
            "FROM observations o JOIN mixes m ON o.mix_uid = m.mix_uid "
            "WHERE o.quantity = 'DoR_SCM' AND o.value_norm IS NOT NULL"
        )
        params: list[Any] = []
        if role is not None:
            sql += " AND m.primary_scm_role = ?"
            params.append(role)
        if paper is not None:
            sql += " AND o.paper_doi = ?"
            params.append(paper)
        if not include_model_system:
            sql += " AND (m.curing_type IS NULL OR m.curing_type != ?)"
            params.append(MODEL_SYSTEM)
        sql += " ORDER BY o.obs_uid"
        return self._q(sql, params)

    def observations_for_mix(
        self,
        mix_uid: str,
        quantities: Sequence[str] | None = None,
        age_window: tuple[float, float] | None = None,
    ) -> list[dict[str, Any]]:
        sql = (
            "SELECT obs_uid, paper_doi, mix_uid, age_d, quantity, phase_name, value_reported, unit_reported, "
            "basis_reported, value_norm, unit_norm, norm_note, method, method_detail, uncertainty, "
            "source_locator, fig_only, extraction_confidence, sanity_ok, reviewed "
            "FROM observations WHERE mix_uid = ? AND value_norm IS NOT NULL"
        )
        params: list[Any] = [mix_uid]
        if quantities:
            qs = [q for q in quantities if q in QUANTITIES]
            if not qs:
                return []
            placeholders = ",".join(["?"] * len(qs))
            sql += f" AND quantity IN ({placeholders})"
            params.extend(qs)
        if age_window is not None:
            sql += " AND age_d >= ? AND age_d <= ?"
            params.extend([float(age_window[0]), float(age_window[1])])
        sql += " ORDER BY quantity, age_d, obs_uid"
        return self._q(sql, params)

    def opc_only_reference(
        self,
        quantity: str,
        age_d: float,
        tol: float = 0.15,
        *,
        w_b_range: tuple[float, float] | None = None,
    ) -> list[dict[str, Any]]:
        """OPC-only mixes (``scm_total_pct`` NULL or 0, not model_system) with an
        observation of ``quantity`` within ``age_d * (1 ± tol)``.

        Mixes whose binder JSON still names an SCM material are returned with
        ``binder_composition_json`` so the caller can exclude them (spec §1.1: the
        393-mix figure is an upper bound)."""
        if quantity not in QUANTITIES:
            raise ValueError(f"unknown quantity {quantity!r}")
        lo, hi = float(age_d) * (1.0 - tol), float(age_d) * (1.0 + tol)
        sql = (
            "SELECT o.obs_uid, o.paper_doi, o.mix_uid, o.age_d, o.quantity, o.phase_name, o.value_norm, "
            "o.unit_norm, o.basis_reported, o.method, o.method_detail, o.uncertainty, o.source_locator, "
            "o.fig_only, o.extraction_confidence, m.binder_composition_json, m.scm_total_pct, m.w_b, "
            "m.curing_temp_C, m.curing_type "
            "FROM observations o JOIN mixes m ON o.mix_uid = m.mix_uid "
            "WHERE o.quantity = ? AND o.value_norm IS NOT NULL AND o.age_d BETWEEN ? AND ? "
            "AND (m.scm_total_pct IS NULL OR m.scm_total_pct = 0) "
            "AND (m.curing_type IS NULL OR m.curing_type != ?) AND m.w_b IS NOT NULL"
        )
        params: list[Any] = [quantity, lo, hi, MODEL_SYSTEM]
        if w_b_range is not None:
            sql += " AND m.w_b BETWEEN ? AND ?"
            params.extend([float(w_b_range[0]), float(w_b_range[1])])
        sql += " ORDER BY o.paper_doi, o.mix_uid, o.age_d"
        return self._q(sql, params)

    def model_system_mixes(self) -> set[str]:
        return {r[0] for r in self.con.execute("SELECT mix_uid FROM mixes WHERE curing_type = ?", (MODEL_SYSTEM,))}


NAMED_QUERIES = {
    "paper": ("doi",),
    "mix": ("mix_uid",),
    "mixes_for_paper": ("doi",),
    "materials_for_paper": ("doi",),
    "dor_observations": ("role", "paper", "include_model_system"),
    "observations_for_mix": ("mix_uid", "quantities", "age_window"),
    "opc_only_reference": ("quantity", "age_d", "tol", "w_b_range"),
}


def run_named_query(db: LiteratureDB, name: str, params: dict[str, Any] | None = None, *, limit: int = 50) -> Any:
    """Dispatch for the ``dor_db_lookup`` tool: only the names above are callable."""
    if name not in NAMED_QUERIES:
        raise KeyError(f"unknown query {name!r}; allowed: {sorted(NAMED_QUERIES)}")
    params = dict(params or {})
    unknown = set(params) - set(NAMED_QUERIES[name])
    if unknown:
        raise ValueError(f"query {name!r} does not accept {sorted(unknown)}")
    if params.get("age_window") is not None:
        params["age_window"] = tuple(params["age_window"])
    if params.get("w_b_range") is not None:
        params["w_b_range"] = tuple(params["w_b_range"])
    result = getattr(db, name)(**params)
    if isinstance(result, list):
        return result[: int(limit)]
    return result
