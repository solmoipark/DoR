"""Staging DB for agent-produced values (spec §11). The literature DB is never touched.

Tables: inferred_dor (reviewed=0 by default), tool_audit.
"""

from __future__ import annotations

import hashlib
import json
import sqlite3
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from .. import __version__
from ..config import staging_db_path

SCHEMA = """
CREATE TABLE IF NOT EXISTS inferred_dor (
  inf_uid TEXT PRIMARY KEY,
  inference_id TEXT, created_at TEXT, dorgems_version TEXT, bundle_hash TEXT,
  mix_uid TEXT,
  scm_json TEXT, mix_json TEXT,
  slot TEXT, age_d REAL,
  alpha_q05 REAL, alpha_q50 REAL, alpha_q95 REAL,
  a_max_q50 REAL, tau_q50 REAL, ess REAL, posterior_method TEXT,
  observations_used_json TEXT, forward_map_path TEXT, run_manifest_path TEXT,
  reviewed INTEGER DEFAULT 0, review_note TEXT
);
CREATE TABLE IF NOT EXISTS tool_audit (ts TEXT, tool TEXT, args_hash TEXT, ok INTEGER, xgems_calls INTEGER, run_dir TEXT);
"""


def open_staging(path: str | Path | None = None) -> sqlite3.Connection:
    p = Path(path) if path else staging_db_path()
    p.parent.mkdir(parents=True, exist_ok=True)
    con = sqlite3.connect(p)
    con.row_factory = sqlite3.Row
    con.executescript(SCHEMA)
    return con


def _now() -> str:
    return datetime.now(timezone.utc).isoformat()


def audit(tool: str, args: dict[str, Any], ok: bool, *, xgems_calls: int | None = None, run_dir: str | None = None, path: str | Path | None = None) -> None:
    h = hashlib.sha256(json.dumps(args, sort_keys=True, default=str).encode()).hexdigest()[:16]
    con = open_staging(path)
    with con:
        con.execute("INSERT INTO tool_audit VALUES (?,?,?,?,?,?)", (_now(), tool, h, int(bool(ok)), xgems_calls, run_dir))
    con.close()


def rows_from_inference(inference: dict[str, Any]) -> list[dict[str, Any]]:
    iid = inference["id"]
    alpha = inference["alpha"]
    rows = []
    for j, age in enumerate(alpha["ages_d"]):
        rows.append(
            {
                "inf_uid": f"dorgems::{iid}::{float(age):g}",
                "inference_id": iid,
                "created_at": _now(),
                "dorgems_version": __version__,
                "bundle_hash": (inference.get("provenance") or {}).get("bundle_bayes"),
                "mix_uid": inference.get("mix_uid"),
                "scm_json": json.dumps(inference.get("scm"), default=str) if inference.get("scm") else None,
                "mix_json": json.dumps(inference.get("mix"), default=str) if inference.get("mix") else None,
                "slot": inference.get("slot"),
                "age_d": float(age),
                "alpha_q05": float(alpha["q05"][j]),
                "alpha_q50": float(alpha["q50"][j]),
                "alpha_q95": float(alpha["q95"][j]),
                "a_max_q50": float(inference["a_max"]["q50"]),
                "tau_q50": float(inference["tau_d"]["q50"]),
                "ess": float(inference.get("ess", 0.0)),
                "posterior_method": inference.get("posterior_method"),
                "observations_used_json": json.dumps([p["label"] for p in inference.get("ppc", [])]),
                "forward_map_path": (inference.get("files") or {}).get("forward_map"),
                "run_manifest_path": (inference.get("files") or {}).get("manifest"),
                "reviewed": 0,
                "review_note": None,
            }
        )
    return rows


def stage_inference(inference: dict[str, Any], *, path: str | Path | None = None, dry_run: bool = True, note: str | None = None) -> dict[str, Any]:
    rows = rows_from_inference(inference)
    if note:
        for r in rows:
            r["review_note"] = note
    if dry_run:
        return {"dry_run": True, "n_rows": len(rows), "preview": rows, "staging_db": str(Path(path) if path else staging_db_path())}
    con = open_staging(path)
    cols = list(rows[0].keys())
    with con:
        con.executemany(f"INSERT OR REPLACE INTO inferred_dor ({','.join(cols)}) VALUES ({','.join('?' * len(cols))})", [tuple(r[c] for c in cols) for r in rows])
    con.close()
    return {"dry_run": False, "n_rows": len(rows), "staging_db": str(Path(path) if path else staging_db_path()), "inference_id": inference["id"]}


def review_list(path: str | Path | None = None, *, reviewed: int | None = None) -> list[dict[str, Any]]:
    con = open_staging(path)
    sql = "SELECT inference_id, mix_uid, slot, COUNT(*) AS n_ages, MIN(age_d) AS age_min, MAX(age_d) AS age_max, AVG(ess) AS ess, MAX(reviewed) AS reviewed, MAX(created_at) AS created_at, MAX(review_note) AS review_note FROM inferred_dor"
    params: list[Any] = []
    if reviewed is not None:
        sql += " WHERE reviewed = ?"
        params.append(int(reviewed))
    sql += " GROUP BY inference_id ORDER BY created_at DESC"
    rows = [dict(r) for r in con.execute(sql, params)]
    con.close()
    return rows


def _set_review(path: str | Path | None, inference_id: str, value: int, note: str | None) -> dict[str, Any]:
    con = open_staging(path)
    with con:
        cur = con.execute("UPDATE inferred_dor SET reviewed = ?, review_note = COALESCE(?, review_note) WHERE inference_id = ?", (value, note, inference_id))
    n = cur.rowcount
    con.close()
    return {"inference_id": inference_id, "reviewed": value, "rows_updated": n}


def review_approve(path: str | Path | None, inference_id: str, *, note: str | None = None) -> dict[str, Any]:
    return _set_review(path, inference_id, 1, note)


def review_reject(path: str | Path | None, inference_id: str, *, note: str | None = None) -> dict[str, Any]:
    return _set_review(path, inference_id, -1, note)


def staging_report(path: str | Path | None = None) -> dict[str, Any]:
    """Where does the model disagree with observations? KL and ESS by slot (spec §9.5)."""
    con = open_staging(path)
    rows = [dict(r) for r in con.execute("SELECT slot, COUNT(DISTINCT inference_id) AS n_inferences, AVG(ess) AS ess_mean, MIN(ess) AS ess_min, SUM(reviewed=1) AS approved, SUM(reviewed=-1) AS rejected FROM inferred_dor GROUP BY slot")]
    audit_rows = [dict(r) for r in con.execute("SELECT tool, COUNT(*) AS n, SUM(ok) AS n_ok, SUM(COALESCE(xgems_calls,0)) AS xgems_calls FROM tool_audit GROUP BY tool")]
    con.close()
    return {"by_slot": rows, "tool_audit": audit_rows}
