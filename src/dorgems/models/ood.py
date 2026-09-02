"""Out-of-distribution assessment against the training feature distribution (spec §5.5)."""

from __future__ import annotations

from typing import Any

import numpy as np

SPARSE_ROLE_MIN_PAPERS = 5


def assess(ood_ref: dict[str, Any], role: str | None, feats: dict[str, Any]) -> dict[str, Any]:
    """Return {flags, score_pct, sparse_role, n_papers_role, role_in_training}."""
    flags: list[str] = []
    roles = (ood_ref or {}).get("roles", {})
    entry = roles.get(str(role)) if role is not None else None
    result: dict[str, Any] = {"flags": flags, "score_pct": None, "sparse_role": False, "n_papers_role": None, "role_in_training": entry is not None}
    if entry is None:
        flags.append(f"role_not_in_training:{role}")
        result["sparse_role"] = True
        return result
    n_papers = int(entry.get("n_papers", 0))
    result["n_papers_role"] = n_papers
    if n_papers < SPARSE_ROLE_MIN_PAPERS:
        result["sparse_role"] = True
        flags.append(f"sparse_role:{role}:{n_papers}_papers")
    for c in ood_ref.get("features", []):
        v = feats.get(c)
        if v is None or (isinstance(v, float) and np.isnan(v)):
            continue
        lo, hi = entry["p01"].get(c), entry["p99"].get(c)
        if lo is not None and v < lo:
            flags.append(f"{c}_below_p01:{v:.3g}<{lo:.3g}")
        elif hi is not None and v > hi:
            flags.append(f"{c}_above_p99:{v:.3g}>{hi:.3g}")
    mh = entry.get("mahalanobis")
    if mh:
        cols = mh["cols"]
        z = np.array([feats.get(c) if feats.get(c) is not None and not (isinstance(feats.get(c), float) and np.isnan(feats.get(c))) else entry["median"][c] for c in cols], float)
        mu = np.asarray(mh["mean"], float)
        inv = np.asarray(mh["cov_inv"], float)
        d2 = float((z - mu) @ inv @ (z - mu))
        ref = np.asarray(mh["d2_sorted"], float)
        pct = float(100.0 * np.searchsorted(ref, d2, side="right") / len(ref))
        result["score_pct"] = pct
        result["mahalanobis_d2"] = d2
        if pct >= 99.0:
            flags.append(f"mahalanobis_beyond_p99:{d2:.2f}")
    return result
