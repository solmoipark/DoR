"""Neighbourhood mode (spec §8.1): compare an existing run_dir (any mix) with the
observation distribution of analogous literature mixes. One grade lower than twin."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from ..db.analogues import find_analogues, load_tolerances
from ..db.features import build_aux_table, scm_input_to_features
from ..db.reader import open_ro
from ..db.units import harmonize
from ..gems.observables import observables_from_run
from .compare import OBS_TO_MODEL, aggregate, compare_rows, write_comparison

GRADE_DOWN = {"A": "B", "B": "C", "C": "D", "D": "D", "X": "X"}


def neighbourhood_compare(lit_db: Path, *, run_dir: str | None, scm: Any, mix: Any, out: Path, quantities: tuple[str, ...], k: int = 10) -> dict[str, Any]:
    from ..pilot.schemas import coerce_mix, coerce_scm

    if not run_dir:
        raise ValueError("neighbourhood mode needs run_dir")
    fdir = Path(run_dir)
    fdir = fdir / "forward" if (fdir / "forward").is_dir() else fdir
    model = observables_from_run(fdir)
    scm_m, mix_m = coerce_scm(scm), coerce_mix(mix)
    feats = scm_input_to_features(scm_m, mix_m)
    con = open_ro(lit_db)
    tol = load_tolerances()
    try:
        ana = find_analogues(con, feats, scm_m.role, k=k, cache_key=str(lit_db))
        aux = build_aux_table(con, [q for q in quantities])
    finally:
        con.close()
    mix_uids = {m["mix_uid"] for m in ana["mixes"]}
    aux = aux[aux["mix_uid"].isin(mix_uids)]
    pairs = []
    lo, hi = tol["age_ratio"]
    for _, o in aux.iterrows():
        col = OBS_TO_MODEL.get(o["quantity"])
        if not col:
            continue
        ratio = model["age_d"] / float(o["age_d"])
        cand = model[(ratio >= lo) & (ratio <= hi)]
        if cand.empty:
            continue
        row = cand.iloc[int(np.argmin(np.abs(np.log(cand["age_d"].values / float(o["age_d"])))))]
        mv = row[col]
        h = harmonize({"quantity": o["quantity"], "value_norm": o["value"], "unit_norm": o["unit_norm"], "basis_reported": o["basis_reported"]}, {"scm_total_pct": o["scm_total_pct"], "w_b": o["w_b"]}, scm_pct=o["scm_pct"])
        pairs.append({"obs_uid": o["obs_uid"], "paper_doi": o["paper_doi"], "mix_uid": o["mix_uid"], "quantity": o["quantity"], "phase_name": o.get("phase_name"), "age_d": float(o["age_d"]), "method": o.get("method"), "grade": GRADE_DOWN.get(h.grade, "D"), "assumptions": "; ".join(h.assumptions + ["neighbourhood mode: grade lowered one step"]), "obs_value": h.value, "model_value": None if mv is None or pd.isna(mv) else float(mv), "uncertainty": o.get("uncertainty"), "source_locator": o.get("source_locator"), "fig_only": o.get("fig_only"), "extraction_confidence": o.get("extraction_confidence")})
    df = compare_rows(pairs)
    agg = aggregate(df)
    files = write_comparison(df, agg, out, header={"mode": "neighbourhood", "target": run_dir, "n_analogue_mixes": ana["n_mixes"], "n_analogue_papers": ana["n_papers"], "analogue_flags": ana["flags"]})
    return {"ok": True, "n_obs": len(pairs), "n_analogue_mixes": ana["n_mixes"], "n_analogue_papers": ana["n_papers"], "aggregate": agg, "files": files, "warnings": ana["flags"]}
