"""Model bundles (spec §5.1): frozen inference artefacts, no PyMC/ArviZ at run time.

bundles/bayes_v4/
    posterior.npz   a0_role (S,3)  t0_role (S,3)  beta_a (S,4)  beta_t (S,4)
                    sd_paper_amax (S,)  sigma_method (S,n_m)
    scaler.json     feats, median, mean, std, roles, methods, beta_shape
    manifest.json   provenance + convergence + LOPO metrics
bundles/gbm_v6/
    model.txt       LightGBM booster
    meta.json       feature order, category codes, sigma_point_pct, metrics
    manifest.json
"""

from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np

from .. import BUNDLE_SCHEMA_VERSION
from ..config import bundles_dir


class BundleError(RuntimeError):
    pass


def sha256_of(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _read_json(path: Path) -> dict[str, Any]:
    if not path.is_file():
        raise BundleError(f"missing bundle file: {path}")
    return json.loads(path.read_text(encoding="utf-8"))


@dataclass
class BayesBundle:
    path: Path
    a0_role: np.ndarray
    t0_role: np.ndarray
    beta_a: np.ndarray
    beta_t: np.ndarray
    sd_paper_amax: np.ndarray
    sigma_method: np.ndarray
    feats: list[str]
    median: dict[str, float]
    mean: dict[str, float]
    std: dict[str, float]
    roles: list[str]
    methods: list[str]
    beta_shape: float
    manifest: dict[str, Any]
    hash: str

    @property
    def n_draws(self) -> int:
        return int(self.a0_role.shape[0])

    def role_index(self, role: str | None) -> tuple[int, bool]:
        """Return (index, pooled) — roles outside {slag, fly_ash} are pooled as 'other'."""
        r = role if role in ("slag", "fly_ash") else "other"
        return self.roles.index(r), r != role

    def standardize(self, x: dict[str, Any]) -> tuple[np.ndarray, list[str]]:
        """Median-impute then z-score in the training order (bayes_hier_v4.py:55-60)."""
        xs = np.zeros(len(self.feats))
        imputed: list[str] = []
        for i, f in enumerate(self.feats):
            v = x.get(f)
            if v is None or (isinstance(v, float) and np.isnan(v)):
                v = self.median[f]
                imputed.append(f)
            sd = self.std[f]
            xs[i] = (float(v) - self.mean[f]) / sd if sd > 0 else 0.0
        return xs, imputed


@dataclass
class GBMBundle:
    path: Path
    booster: Any
    features: list[str]
    categories: dict[str, list[str]]
    sigma_point_pct: float
    meta: dict[str, Any]
    manifest: dict[str, Any]
    hash: str


@dataclass
class Bundle:
    bayes: BayesBundle | None = None
    gbm: GBMBundle | None = None
    ood: dict[str, Any] = field(default_factory=dict)

    @property
    def provenance(self) -> dict[str, Any]:
        return {
            "bundle_bayes": f"sha256:{self.bayes.hash}" if self.bayes else None,
            "bundle_gbm": f"sha256:{self.gbm.hash}" if self.gbm else None,
            "training_db_sha256": (self.bayes.manifest.get("training_db_sha256") if self.bayes else None),
        }


def _check_schema(manifest: dict[str, Any], where: Path) -> None:
    v = manifest.get("bundle_schema_version")
    if v != BUNDLE_SCHEMA_VERSION:
        raise BundleError(f"{where}: bundle_schema_version {v!r} != {BUNDLE_SCHEMA_VERSION}")


def _check_hashes(manifest: dict[str, Any], folder: Path) -> str:
    files = manifest.get("files") or {}
    if not files:
        raise BundleError(f"{folder}: manifest lists no files")
    combined = hashlib.sha256()
    for name, expected in sorted(files.items()):
        p = folder / name
        if not p.is_file():
            raise BundleError(f"{folder}: missing {name}")
        actual = sha256_of(p)
        if actual != expected:
            raise BundleError(f"{folder}/{name}: sha256 mismatch (bundle corrupted or edited)")
        combined.update(actual.encode())
    return combined.hexdigest()


def load_bayes_bundle(folder: Path) -> BayesBundle:
    manifest = _read_json(folder / "manifest.json")
    _check_schema(manifest, folder)
    h = _check_hashes(manifest, folder)
    post = np.load(folder / "posterior.npz")
    scaler = _read_json(folder / "scaler.json")
    required = ["a0_role", "t0_role", "beta_a", "beta_t", "sd_paper_amax", "sigma_method"]
    for k in required:
        if k not in post.files:
            raise BundleError(f"posterior.npz lacks {k}")
    S = post["a0_role"].shape[0]
    if post["a0_role"].shape != (S, len(scaler["roles"])):
        raise BundleError("a0_role shape does not match scaler.roles")
    if post["beta_a"].shape != (S, len(scaler["feats"])):
        raise BundleError("beta_a shape does not match scaler.feats")
    if post["sigma_method"].shape != (S, len(scaler["methods"])):
        raise BundleError("sigma_method shape does not match scaler.methods")
    return BayesBundle(
        path=folder,
        a0_role=post["a0_role"],
        t0_role=post["t0_role"],
        beta_a=post["beta_a"],
        beta_t=post["beta_t"],
        sd_paper_amax=post["sd_paper_amax"],
        sigma_method=post["sigma_method"],
        feats=list(scaler["feats"]),
        median={k: float(v) for k, v in scaler["median"].items()},
        mean={k: float(v) for k, v in scaler["mean"].items()},
        std={k: float(v) for k, v in scaler["std"].items()},
        roles=list(scaler["roles"]),
        methods=list(scaler["methods"]),
        beta_shape=float(scaler.get("beta_shape", 0.5)),
        manifest=manifest,
        hash=h,
    )


def load_gbm_bundle(folder: Path) -> GBMBundle:
    manifest = _read_json(folder / "manifest.json")
    _check_schema(manifest, folder)
    h = _check_hashes(manifest, folder)
    meta = _read_json(folder / "meta.json")
    import lightgbm as lgb

    booster = lgb.Booster(model_file=str(folder / "model.txt"))
    feats = list(meta["features"])
    if booster.num_feature() != len(feats):
        raise BundleError("model.txt feature count != meta.features")
    return GBMBundle(
        path=folder,
        booster=booster,
        features=feats,
        categories={k: list(v) for k, v in meta["categories"].items()},
        sigma_point_pct=float(meta.get("sigma_point_pct", 12.0)),
        meta=meta,
        manifest=manifest,
        hash=h,
    )


def load_bundle(path: str | Path | None = None, *, require_gbm: bool = True) -> Bundle:
    root = Path(path) if path else bundles_dir()
    bayes = load_bayes_bundle(root / "bayes_v4")
    gbm = None
    if (root / "gbm_v6" / "manifest.json").is_file():
        gbm = load_gbm_bundle(root / "gbm_v6")
    elif require_gbm:
        raise BundleError(f"gbm_v6 bundle missing under {root}")
    ood_path = root / "ood_reference.json"
    ood = _read_json(ood_path) if ood_path.is_file() else {}
    return Bundle(bayes=bayes, gbm=gbm, ood=ood)


def write_manifest(folder: Path, extra: dict[str, Any], files: list[str]) -> dict[str, Any]:
    from .. import __version__

    manifest = {
        "bundle_schema_version": BUNDLE_SCHEMA_VERSION,
        "dorgems_version": __version__,
        "files": {name: sha256_of(folder / name) for name in files},
    }
    manifest.update(extra)
    (folder / "manifest.json").write_text(json.dumps(manifest, indent=2, ensure_ascii=False), encoding="utf-8")
    return manifest
