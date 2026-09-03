"""Observation likelihood on the alpha-grid forward map (spec §9.3).

    log L(a_max, τ) = Σ_iq log N(y_iq | F_q(α(t_i); t_i) + b_q, σ_obs,iq² + σ_model,q²)
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

import numpy as np

from ..kinetics.curves import stretched_exp
from .alpha_grid import OBS_KEYS, ForwardMap


@dataclass
class ObsPoint:
    age_index: int
    age_d: float
    quantity: str
    y: float  # harmonised (target units)
    sigma_obs: float
    sigma_model: float
    grade: str
    label: str = ""
    weight: float = 1.0

    @property
    def var(self) -> float:
        return self.sigma_obs**2 + self.sigma_model**2


@dataclass
class Likelihood:
    fmap: ForwardMap
    points: list[ObsPoint]
    offsets: dict[str, float] = field(default_factory=dict)  # b_q
    beta: float = 0.5

    def mu(self, a_max: np.ndarray, tau: np.ndarray, p: ObsPoint) -> np.ndarray:
        a = np.clip(a_max * (1.0 - np.exp(-((p.age_d / tau) ** self.beta))), 0.0, 1.0)
        return self.fmap.value(p.quantity, p.age_index, a) + float(self.offsets.get(p.quantity, 0.0))

    def loglik(self, a_max: np.ndarray, tau: np.ndarray, *, per_point: bool = False) -> np.ndarray | tuple[np.ndarray, np.ndarray]:
        a_max = np.asarray(a_max, float)
        tau = np.asarray(tau, float)
        total = np.zeros_like(a_max)
        contrib = np.zeros((len(self.points), a_max.shape[0]))
        for k, p in enumerate(self.points):
            m = self.mu(a_max, tau, p)
            ll = -0.5 * ((p.y - m) ** 2 / p.var + np.log(2 * np.pi * p.var))
            ll = np.where(np.isfinite(ll), ll, -1e6)
            contrib[k] = p.weight * ll
            total += contrib[k]
        return (total, contrib) if per_point else total

    def predicted(self, a_max: np.ndarray, tau: np.ndarray) -> np.ndarray:
        """(n_points, S) model means — for posterior predictive checks."""
        return np.stack([self.mu(a_max, tau, p) for p in self.points])


def build_points(observations: list[dict[str, Any]], ages: list[float], *, sigma_obs_default: dict[str, float], sigma_model: dict[str, float], use_direct_dor: bool = False, weights: dict[str, float] | None = None) -> tuple[list[ObsPoint], list[str]]:
    """``observations``: harmonised dicts with age_d, quantity, value, grade, uncertainty, obs_uid.
    ``weights`` (per quantity) scale each observation's log-likelihood contribution; 0 excludes."""
    pts: list[ObsPoint] = []
    skipped: list[str] = []
    weights = dict(weights or {})
    for o in observations:
        q = o["quantity"]
        if q in weights and float(weights[q]) <= 0.0:
            skipped.append(f"{o.get('obs_uid', q)}: quantity {q} has weight 0")
            continue
        if q in ("DoR_SCM", "DoR_clinker") and not use_direct_dor:
            skipped.append(f"{o.get('obs_uid', q)}: direct DoR kept for validation only")
            continue
        if q not in OBS_KEYS:
            skipped.append(f"{o.get('obs_uid', q)}: quantity {q} not in the forward map")
            continue
        if o.get("grade") not in ("A", "B") or o.get("value") is None:
            skipped.append(f"{o.get('obs_uid', q)}: grade {o.get('grade')} excluded")
            continue
        try:
            i = ages.index(float(o["age_d"]))
        except ValueError:
            skipped.append(f"{o.get('obs_uid', q)}: age {o['age_d']} not in the forward map")
            continue
        so = float(o["uncertainty"]) if o.get("uncertainty") else float(sigma_obs_default.get(q, 1.5))
        pts.append(ObsPoint(i, float(o["age_d"]), q, float(o["value"]), so, float(sigma_model.get(q, 2.5)), str(o.get("grade")), label=str(o.get("obs_uid", f"{q}@{o['age_d']}")), weight=float(weights.get(q, 1.0))))
    return pts, skipped
