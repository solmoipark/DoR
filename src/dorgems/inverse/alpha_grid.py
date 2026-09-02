"""alpha-grid forward map F_q(alpha; t_i) (spec §9.2).

For every observed age the kernel is run with alpha pinned on a grid; each
observable is then a monotone (PCHIP) interpolant in alpha, so curve parameters
(a_max, tau) are evaluated afterwards without further xGEMS calls.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from scipy.interpolate import PchipInterpolator

from ..gems.forward import run_forward
from ..gems.observables import observables_from_run
from ..kinetics.reaction_model import pin_reaction_model

OBS_KEYS = {"CH_TGA": "CH_g", "CH_XRD": "CH_g", "bound_water": "bound_water_g", "chem_shrink": "chem_shrink_ml_g"}


def default_alpha_grid(n: int = 21, *, refine_interval: tuple[float, float] | None = None, refine_step: float = 0.025) -> np.ndarray:
    g = np.linspace(0.0, 1.0, int(n))
    if refine_interval is not None:
        lo, hi = max(0.0, refine_interval[0]), min(1.0, refine_interval[1])
        g = np.union1d(g, np.arange(lo, hi + 1e-9, refine_step))
    return np.round(np.unique(np.clip(g, 0, 1)), 4)


class ForwardMap:
    def __init__(self, ages: np.ndarray, alphas: np.ndarray, table: dict[str, np.ndarray], meta: dict[str, Any]):
        self.ages = np.asarray(ages, float)
        self.alphas = np.asarray(alphas, float)
        self.table = {k: np.asarray(v, float) for k, v in table.items()}  # (n_ages, n_alphas)
        self.meta = meta
        self._interp: dict[tuple[str, int], PchipInterpolator] = {}

    def value(self, quantity: str, age_index: int, alpha: np.ndarray | float) -> np.ndarray:
        key = OBS_KEYS.get(quantity, quantity)
        row = self.table[key][age_index]
        ok = np.isfinite(row)
        if ok.sum() < 2:
            return np.full(np.shape(alpha), np.nan)
        k = (key, age_index)
        if k not in self._interp:
            self._interp[k] = PchipInterpolator(self.alphas[ok], row[ok], extrapolate=True)
        return self._interp[k](np.clip(np.asarray(alpha, float), 0.0, 1.0))

    def monotonicity_report(self) -> dict[str, Any]:
        rep: dict[str, Any] = {}
        for key, tab in self.table.items():
            flags = []
            for i, row in enumerate(tab):
                d = np.diff(row[np.isfinite(row)])
                if key == "CH_g" and np.any(d > 1e-9):
                    flags.append({"age_d": float(self.ages[i]), "issue": "CH increases with alpha (pozzolanic reaction should consume CH)"})
                if key == "bound_water_g" and np.any(d < -1e-9):
                    flags.append({"age_d": float(self.ages[i]), "issue": "bound water decreases with alpha"})
            rep[key] = flags
        return rep

    def save(self, out: Path) -> dict[str, str]:
        out.mkdir(parents=True, exist_ok=True)
        np.savez_compressed(out / "forward_map.npz", ages=self.ages, alphas=self.alphas, **self.table)
        (out / "forward_map_meta.json").write_text(json.dumps(self.meta, indent=2, default=str), encoding="utf-8")
        rep = self.monotonicity_report()
        lines = ["# alpha-grid forward map", "", f"ages (d): {self.ages.tolist()}", f"alpha grid: {self.alphas.tolist()}", f"xGEMS calls (non-cached): {self.meta.get('xgems_calls')}", ""]
        for key, tab in self.table.items():
            lines += [f"## {key}", "", "| age d | " + " | ".join(f"α={a:g}" for a in self.alphas) + " |", "|---|" + "---|" * len(self.alphas)]
            for i, row in enumerate(tab):
                lines.append(f"| {self.ages[i]:g} | " + " | ".join("—" if not np.isfinite(v) else f"{v:.3g}" for v in row) + " |")
            if rep.get(key):
                lines += ["", f"monotonicity flags: {rep[key]}"]
            lines.append("")
        (out / "forward_map_report.md").write_text("\n".join(lines), encoding="utf-8")
        return {"forward_map": str(out / "forward_map.npz"), "forward_map_report": str(out / "forward_map_report.md")}

    @classmethod
    def load(cls, path: Path) -> "ForwardMap":
        z = np.load(path)
        table = {k: z[k] for k in z.files if k not in ("ages", "alphas")}
        meta_p = path.parent / "forward_map_meta.json"
        meta = json.loads(meta_p.read_text(encoding="utf-8")) if meta_p.is_file() else {}
        return cls(z["ages"], z["alphas"], table, meta)


def build_forward_map(
    forward_query_base: dict[str, Any],
    *,
    slot: str,
    ages: list[float],
    alphas: np.ndarray,
    out: Path,
    ig_db: str | Path,
    materials_config: str | Path | None,
    use_mock: bool = True,
    dat_lst: str | Path | None = None,
    max_xgems_calls: int | None = None,
    quantities: tuple[str, ...] = ("CH_TGA", "bound_water", "chem_shrink"),
) -> ForwardMap:
    """Runs every (alpha_j) config over all ages at once (one forward query per alpha);
    the kernel cache dedups identical element vectors."""
    keys = sorted({OBS_KEYS[q] for q in quantities if q in OBS_KEYS})
    table = {k: np.full((len(ages), len(alphas)), np.nan) for k in keys}
    calls = 0
    warnings: list[str] = []
    budget_left = max_xgems_calls
    for j, a in enumerate(alphas):
        cfg = pin_reaction_model(float(a), slot, out / "pins", config_id="grid", signature_files=[str(materials_config)] if materials_config else None)
        fq = dict(forward_query_base)
        fq["age_grid"] = {"values": [float(t) for t in ages]}
        fq["name"] = f"alpha_grid_{j:02d}"
        res = run_forward(fq, out=out / "runs" / f"a{a:.3f}", db=ig_db, reaction_model_config=cfg, materials_config=materials_config, slot=slot, use_mock=use_mock, dat_lst=dat_lst, max_xgems_calls=budget_left, capture_species=True)
        if not res.ok:
            warnings.append(f"alpha={a:.3f}: {res.error or res.warnings}")
            continue
        obs = observables_from_run(res.forward_dir)
        n_new = int((res.time_series["reused_cache"] == False).sum()) if res.time_series is not None and "reused_cache" in res.time_series else len(ages)  # noqa: E712
        calls += n_new
        if budget_left is not None:
            budget_left = max(0, budget_left - n_new)
        for i, t in enumerate(ages):
            row = obs[obs["age_d"] == float(t)]
            if row.empty:
                continue
            for k in keys:
                v = row.iloc[0].get(k)
                table[k][i, j] = np.nan if v is None or pd.isna(v) else float(v)
    fm = ForwardMap(np.asarray(ages, float), np.asarray(alphas, float), table, {"slot": slot, "use_mock": use_mock, "xgems_calls": calls, "warnings": warnings, "materials_config": str(materials_config) if materials_config else None})
    fm.save(out)
    return fm
