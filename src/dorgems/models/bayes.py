"""Posterior-predictive DoR curves from the hierarchical Bayesian v4 bundle (spec §5.2).

alpha(t) = a_max · (1 − exp(−(t/τ)^β)),  β fixed (0.5)
logit(a_max/100) = a0[role] + xs·beta_a + u_paper,  u_paper ~ N(0, sd_paper)   (new study)
log τ = t0[role] + xs·beta_t
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

import numpy as np

from .bundle import BayesBundle


def stretched_exp_pct(t: np.ndarray, a_max_pct: np.ndarray, tau: np.ndarray, beta: float) -> np.ndarray:
    """(S,) params × (T,) ages → (S,T) alpha in %."""
    t = np.asarray(t, dtype=float)[None, :]
    return a_max_pct[:, None] * (1.0 - np.exp(-((t / tau[:, None]) ** beta)))


@dataclass
class CurveDraws:
    ages: np.ndarray
    a_max: np.ndarray  # (S,) %
    tau: np.ndarray  # (S,) days
    alpha: np.ndarray  # (S,T) %
    sigma_obs: np.ndarray | None  # (S,) %p, method-conditional
    role_used: str
    imputed: list[str]
    flags: list[str] = field(default_factory=list)
    weights: np.ndarray | None = None  # importance weights (ensemble), None = uniform

    def summary(self, quantiles: tuple[float, ...] = (0.05, 0.5, 0.95), *, rng_seed: int = 0) -> dict[str, Any]:
        w = self.weights
        q = list(quantiles)
        out: dict[str, Any] = {
            "a_max": _wq(self.a_max, q, w),
            "tau_d": _wq(self.tau, q, w),
            "alpha_pct": {"latent": _wq_cols(self.alpha, q, w)},
        }
        out["a_max"]["mean"] = float(np.average(self.a_max, weights=w))
        out["tau_d"]["mean"] = float(np.average(self.tau, weights=w))
        out["alpha_pct"]["latent"]["mean"] = [float(v) for v in np.average(self.alpha, axis=0, weights=w)]
        if self.sigma_obs is not None:
            rng = np.random.default_rng(rng_seed)
            noise = rng.normal(0.0, 1.0, size=self.alpha.shape) * self.sigma_obs[:, None]
            obs = np.clip(self.alpha + noise, 0.0, 100.0)
            out["alpha_pct"]["observed"] = _wq_cols(obs, q, w)
            out["alpha_pct"]["observed"]["mean"] = [float(v) for v in np.average(obs, axis=0, weights=w)]
        return out


def weighted_quantile(values: np.ndarray, qs: list[float], weights: np.ndarray | None = None) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    if weights is None:
        return np.quantile(values, qs)
    order = np.argsort(values)
    v = values[order]
    w = np.asarray(weights, dtype=float)[order]
    cw = np.cumsum(w)
    cw = cw / cw[-1]
    return np.interp(qs, cw - 0.5 * w / w.sum(), v, left=v[0], right=v[-1])


def _wq(v: np.ndarray, qs: list[float], w: np.ndarray | None) -> dict[str, float]:
    vals = weighted_quantile(v, qs, w)
    return {f"q{int(round(q * 100)):02d}": float(x) for q, x in zip(qs, vals)}


def _wq_cols(m: np.ndarray, qs: list[float], w: np.ndarray | None) -> dict[str, list[float]]:
    cols = [weighted_quantile(m[:, j], qs, w) for j in range(m.shape[1])]
    arr = np.asarray(cols)  # (T, nq)
    return {f"q{int(round(q * 100)):02d}": [float(x) for x in np.clip(arr[:, i], 0.0, 100.0)] for i, q in enumerate(qs)}


def predict_curve(
    bundle: BayesBundle,
    x: dict[str, Any],
    role: str | None,
    ages: np.ndarray | list[float],
    *,
    new_study: bool = True,
    method_group: str | None = None,
    rng_seed: int = 0,
) -> CurveDraws:
    """Posterior draws of (a_max, tau, alpha(ages)) — the OOF algorithm of bayes_hier_v4.py:146-163."""
    ages = np.asarray(ages, dtype=float)
    if np.any(ages < 0):
        raise ValueError("ages must be non-negative")
    ridx, pooled = bundle.role_index(role)
    role_used = bundle.roles[ridx]
    flags: list[str] = []
    if pooled:
        flags.append(f"role_pooled_as_other:{role}")
    xs, imputed = bundle.standardize(x)
    rng = np.random.default_rng(rng_seed)
    eta = bundle.a0_role[:, ridx] + bundle.beta_a @ xs
    if new_study:
        eta = eta + rng.normal(0.0, 1.0, size=eta.shape) * bundle.sd_paper_amax
    a_max = 100.0 / (1.0 + np.exp(-eta))
    tau = np.exp(bundle.t0_role[:, ridx] + bundle.beta_t @ xs)
    alpha = np.clip(stretched_exp_pct(ages, a_max, tau, bundle.beta_shape), 0.0, 100.0)
    sigma_obs = None
    if method_group is not None:
        if method_group in bundle.methods:
            sigma_obs = bundle.sigma_method[:, bundle.methods.index(method_group)]
        else:
            # unknown method label → the 'unknown' noise scale if trained, else pooled mean
            if "unknown" in bundle.methods:
                sigma_obs = bundle.sigma_method[:, bundle.methods.index("unknown")]
            else:
                sigma_obs = bundle.sigma_method.mean(axis=1)
            flags.append(f"method_group_unseen:{method_group}")
    return CurveDraws(ages=ages, a_max=a_max, tau=tau, alpha=alpha, sigma_obs=sigma_obs, role_used=role_used, imputed=imputed, flags=flags)


def role_kinetics_at_mean(bundle: BayesBundle) -> dict[str, dict[str, float]]:
    """Average-condition (xs = 0, no paper effect) posterior means — the golden of bayes_role_kinetics.csv."""
    out = {}
    for i, r in enumerate(bundle.roles):
        a0 = bundle.a0_role[:, i]
        t0 = bundle.t0_role[:, i]
        out[r] = {
            "a_max_pct": float(100.0 / (1.0 + np.exp(-a0.mean()))),
            "a_max_hdi3": float(100.0 / (1.0 + np.exp(-np.percentile(a0, 3)))),
            "a_max_hdi97": float(100.0 / (1.0 + np.exp(-np.percentile(a0, 97)))),
            "tau_days": float(np.exp(t0.mean())),
        }
    return out
