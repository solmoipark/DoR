"""Fit the InverseGems five-parameter logistic to an arbitrary alpha(t) curve (spec §6.1, mode logistic_fit)."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np
from scipy.optimize import least_squares

from .curves import five_param_logistic

# identical to inverse_gems.kinetics_calibration._DEFAULT_BOUNDS
BOUNDS = {"A": (0.0, 0.2), "B": (0.05, 5.0), "C": (0.1, 500.0), "D": (0.05, 1.0), "G": (0.05, 5.0)}
INIT = {"A": 0.0, "B": 0.8, "C": 20.0, "D": 0.6, "G": 1.0}


@dataclass
class LogisticFit:
    params: dict[str, float]
    max_abs_dev_pct: float  # on report_range
    rmse_pct: float
    report_range_d: tuple[float, float]
    status: str  # ok | warn | fail
    message: str

    def to_dict(self) -> dict[str, Any]:
        return {
            "params": dict(self.params),
            "max_abs_dev_pct": self.max_abs_dev_pct,
            "rmse_pct": self.rmse_pct,
            "report_range_d": list(self.report_range_d),
            "status": self.status,
            "message": self.message,
        }


def default_grid(t_min: float = 0.25, t_max: float = 730.0, n: int = 48) -> np.ndarray:
    return np.logspace(np.log10(t_min), np.log10(t_max), n)


def logistic_fit(
    t: np.ndarray,
    alpha_fraction: np.ndarray,
    *,
    fix_A_zero: bool = True,
    report_range_d: tuple[float, float] = (1.0, 365.0),
    warn_dev_pct: float = 2.0,
    fail_dev_pct: float = 5.0,
    multi_start: bool = True,
) -> LogisticFit:
    t = np.asarray(t, float)
    y = np.clip(np.asarray(alpha_fraction, float), 0.0, 1.0)
    if t.shape != y.shape or t.ndim != 1:
        raise ValueError("t and alpha must be 1-D arrays of equal length")
    if np.any(t <= 0):
        raise ValueError("t must be positive")
    names = ["B", "C", "D", "G"] if fix_A_zero else ["A", "B", "C", "D", "G"]
    lo = np.array([BOUNDS[n][0] for n in names])
    hi = np.array([BOUNDS[n][1] for n in names])

    def unpack(p: np.ndarray) -> dict[str, float]:
        d = dict(zip(names, map(float, p)))
        if fix_A_zero:
            d["A"] = 0.0
        return d

    def resid(p: np.ndarray) -> np.ndarray:
        d = unpack(p)
        return five_param_logistic(t, d["A"], d["B"], d["C"], d["D"], d["G"]) - y

    ymax = float(y.max())
    starts = [[INIT[n] for n in names]]
    if multi_start:
        for B0 in (0.5, 1.0, 1.5):
            for C0 in (5.0, 20.0, 60.0):
                d0 = float(np.clip(ymax * 1.05, BOUNDS["D"][0], BOUNDS["D"][1]))
                s = {"A": 0.0, "B": B0, "C": C0, "D": d0, "G": 1.0}
                starts.append([s[n] for n in names])
    best = None
    for s0 in starts:
        x0 = np.clip(np.array(s0, float), lo, hi)
        try:
            r = least_squares(resid, x0=x0, bounds=(lo, hi), method="trf", max_nfev=2000)
        except Exception:  # noqa: BLE001
            continue
        if best is None or r.cost < best.cost:
            best = r
    if best is None:
        raise RuntimeError("logistic fit failed for every start")
    params = unpack(best.x)
    mask = (t >= report_range_d[0]) & (t <= report_range_d[1])
    dev = np.abs(resid(best.x)) * 100.0
    mad = float(dev[mask].max()) if mask.any() else float(dev.max())
    rmse = float(np.sqrt(np.mean(dev[mask] ** 2))) if mask.any() else float(np.sqrt(np.mean(dev**2)))
    if mad > fail_dev_pct:
        status, msg = "fail", f"max |dev| {mad:.2f} %p > {fail_dev_pct} %p on {report_range_d} d; use mode=native"
    elif mad > warn_dev_pct:
        status, msg = "warn", f"max |dev| {mad:.2f} %p > {warn_dev_pct} %p on {report_range_d} d"
    else:
        status, msg = "ok", f"max |dev| {mad:.2f} %p"
    return LogisticFit(params={k: params[k] for k in ["A", "B", "C", "D", "G"]}, max_abs_dev_pct=mad, rmse_pct=rmse, report_range_d=report_range_d, status=status, message=msg)
