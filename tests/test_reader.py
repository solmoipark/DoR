"""G0-4: the literature DB is opened read-only; every write path fails."""

from __future__ import annotations

import sqlite3

import pytest

from dorgems.db.reader import LiteratureDB, open_ro, run_named_query


def test_open_ro_refuses_insert(fixture_db_copy):
    con = open_ro(fixture_db_copy)
    with pytest.raises(sqlite3.OperationalError):
        con.execute("INSERT INTO papers (doi) VALUES ('x')")
    with pytest.raises(sqlite3.OperationalError):
        con.execute("UPDATE observations SET value_norm = 0")
    with pytest.raises(sqlite3.OperationalError):
        con.execute("DELETE FROM mixes")
    with pytest.raises(sqlite3.OperationalError):
        con.execute("CREATE TABLE t (x)")
    # the file itself is untouched: a fresh writable connection still sees the same counts
    n_before = sqlite3.connect(fixture_db_copy).execute("SELECT COUNT(*) FROM papers").fetchone()[0]
    assert n_before == 5


def test_named_queries_and_provenance(fixture_db_path):
    with LiteratureDB(fixture_db_path) as db:
        counts = db.counts()
        assert counts["papers"] == 5
        rows = db.dor_observations()
        assert rows and {"paper_doi", "source_locator", "extraction_confidence", "fig_only"} <= set(rows[0])
        mix = db.mix(rows[0]["mix_uid"])
        assert mix is not None and mix["paper_doi"] == rows[0]["paper_doi"]
        obs = db.observations_for_mix(rows[0]["mix_uid"], ["DoR_SCM", "CH_TGA"])
        assert all(o["quantity"] in ("DoR_SCM", "CH_TGA") for o in obs)
        assert db.observations_for_mix(rows[0]["mix_uid"], ["not_a_quantity"]) == []
        mats = db.materials_for_paper(rows[0]["paper_doi"])
        assert mats
        ms = db.model_system_mixes()
        assert isinstance(ms, set)
        with pytest.raises(ValueError):
            db.opc_only_reference("nope", 28)


def test_run_named_query_is_allowlisted(fixture_db_path):
    with LiteratureDB(fixture_db_path) as db:
        with pytest.raises(KeyError):
            run_named_query(db, "DROP TABLE papers")
        with pytest.raises(ValueError):
            run_named_query(db, "paper", {"doi": "x", "sql": "1"})
        out = run_named_query(db, "dor_observations", {}, limit=3)
        assert len(out) == 3
