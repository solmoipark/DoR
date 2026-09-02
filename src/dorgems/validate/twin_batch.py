"""G2-5: twin comparisons over many literature SCM mixes with a pinned (measured) DoR."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pandas as pd

from ..db.features import build_dor_table
from ..db.reader import LiteratureDB
from .twin import COMPARE_QUANTITIES, twin_compare_mix


def candidate_mixes(db: LiteratureDB, quantities: tuple[str, ...] = COMPARE_QUANTITIES, *, min_dor_ages: int = 3, min_common_ages: int = 1) -> pd.DataFrame:
    """Mixes with ≥ min_dor_ages measured DoR ages and observations of ``quantities`` at
    ≥ min_common_ages of those ages (spec §1.1: 66 mixes with ≥3 common ages)."""
    table = build_dor_table(db.con)
    table = table[table["system_type"] != "model_system"]
    dor_ages = table.groupby("mix_uid")["age_d"].agg(lambda s: sorted(set(s)))
    rows = []
    ph = ",".join(["?"] * len(quantities))
    for mix_uid, ages in dor_ages.items():
        if len(ages) < min_dor_ages:
            continue
        obs = db.con.execute(f"SELECT DISTINCT age_d FROM observations WHERE mix_uid = ? AND value_norm IS NOT NULL AND quantity IN ({ph})", (mix_uid, *quantities)).fetchall()
        common = sorted(set(float(a[0]) for a in obs if a[0] is not None) & set(map(float, ages)))
        if len(common) >= min_common_ages:
            rows.append({"mix_uid": mix_uid, "paper_doi": table.loc[table["mix_uid"] == mix_uid, "paper_doi"].iloc[0], "n_dor_ages": len(ages), "n_common_ages": len(common)})
    return pd.DataFrame(rows)


def twin_batch(lit_db: Path, *, out: Path, ig_db: str | Path, use_mock: bool = True, dat_lst: str | Path | None = None, max_xgems_calls: int | None = None, quantities: tuple[str, ...] = COMPARE_QUANTITIES, max_mixes: int | None = None) -> dict[str, Any]:
    out.mkdir(parents=True, exist_ok=True)
    results = []
    warnings: list[str] = []
    with LiteratureDB(lit_db) as db:
        cands = candidate_mixes(db, quantities)
        if max_mixes:
            cands = cands.head(int(max_mixes))
        for _, c in cands.iterrows():
            sub = out / str(c["mix_uid"]).replace("/", "_").replace(":", "_")
            try:
                r = twin_compare_mix(db, c["mix_uid"], out=sub, ig_db=ig_db, use_mock=use_mock, dat_lst=dat_lst, max_xgems_calls=max_xgems_calls, quantities=quantities)
            except Exception as exc:  # noqa: BLE001
                r = {"ok": False, "mix_uid": c["mix_uid"], "error": f"{type(exc).__name__}: {exc}"}
            results.append({"mix_uid": c["mix_uid"], "paper_doi": c["paper_doi"], "ok": r.get("ok"), "error": r.get("error"), "dor_source": r.get("dor_source"), "verdict": (r.get("aggregate") or {}).get("overall"), "by_quantity": {q: e.get("verdict") for q, e in ((r.get("aggregate") or {}).get("by_quantity") or {}).items()}, "n_obs": r.get("n_obs")})
    df = pd.DataFrame(results)
    df.to_csv(out / "twin_batch.csv", index=False)
    dist = df["verdict"].value_counts(dropna=False).to_dict() if not df.empty else {}
    summary = {"n_candidates": int(len(cands)), "n_ok": int(df["ok"].fillna(False).astype(bool).sum()) if not df.empty else 0, "verdict_distribution": {str(k): int(v) for k, v in dist.items()}, "use_mock": use_mock}
    (out / "twin_batch.json").write_text(json.dumps({**summary, "results": results}, indent=2, default=str), encoding="utf-8")
    md = ["# Twin batch (G2-5)", "", f"- candidates: {summary['n_candidates']}, ran ok: {summary['n_ok']}, mock: {use_mock}", "", "| verdict | n |", "|---|---|"] + [f"| {k} | {v} |" for k, v in summary["verdict_distribution"].items()]
    (out / "summary.md").write_text("\n".join(md) + "\n", encoding="utf-8")
    return {"ok": True, **summary, "files": {"twin_batch_csv": str(out / "twin_batch.csv"), "twin_batch_json": str(out / "twin_batch.json"), "summary": str(out / "summary.md")}, "warnings": warnings}
