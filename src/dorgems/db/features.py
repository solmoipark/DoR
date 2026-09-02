"""Importable port of ``modeling/build_clean.py`` and the feature logic of
``modeling/multitask_v6.py`` (spec §4.2).

The resolver chain, the method grouping and the deduplication key are kept
byte-for-byte equivalent to the audited scripts; golden tests (G0-1) enforce
that :func:`build_dor_table` reproduces ``modeling/dor_scm_final.csv``.
"""

from __future__ import annotations

import json
import re
import sqlite3
from dataclasses import dataclass
from typing import Any

import numpy as np
import pandas as pd

from .reader import MODEL_SYSTEM

# ---------------------------------------------------------------------------
# build_clean.py:14-25
# ---------------------------------------------------------------------------


def method_group(m: str | None) -> str:
    if m is None:
        return "unknown"
    s = m.lower().replace(" ", "_")
    if "selective" in s or "dissolution" in s or s.startswith("hcl"):
        return "selective_dissolution"
    if "sem" in s or "bse" in s or "image" in s:
        return "SEM_BSE"
    if "ponkcs" in s:
        return "XRD_PONKCS"
    if "nmr" in s:
        return "NMR"
    if "calorimetry" in s or "heat" in s:
        return "calorimetry"
    if "xrd" in s or "rietveld" in s:
        return "XRD_other"
    if "mass_balance" in s or "balance" in s:
        return "mass_balance"
    if "tga" in s:
        return "TGA"
    return "other_unknown"


# build_clean.py:31-33
CEMENT_ROLES = {
    "cement",
    "clinker",
    "gypsum",
    "sulfate_source",
    "hydrated_lime",
    "lime",
    "activator",
    "accelerator_blend",
    "expansive_agent",
    "expansive_additive",
}
FILLER_ROLES = {"quartz_filler", "quartz"}


# build_clean.py:35-47
def role_from_text(*texts: str | None) -> str | None:
    t = " ".join(x.lower() for x in texts if x)
    if not t:
        return None
    if "fly ash" in t or re.search(r"\bfa\d*\b", t) or "pfa" in t:
        return "fly_ash"
    if "ggbs" in t or "ggbfs" in t or "slag" in t or re.search(r"\bbfs\b", t):
        return "slag"
    if "metakaolin" in t or re.search(r"\bmk\d*\b", t):
        return "metakaolin"
    if "silica fume" in t or re.search(r"\bsf\d*\b", t) or "microsilica" in t:
        return "silica_fume"
    if "calcined clay" in t or "lc3" in t or re.search(r"\bcc\d*\b", t):
        return "calcined_clay"
    if "limestone" in t or re.search(r"\bls\d*\b", t):
        return "limestone"
    if "pozzolan" in t or "zeolite" in t or "pumice" in t:
        return "natural_pozzolan"
    if "rice husk" in t or "rha" in t:
        return "rice_husk_ash"
    if "glass" in t:
        return "glass_powder"
    return None


# build_clean.py:49-50
_PAT_PROPS = re.compile(r"([A-Z]+\d*)\s+(\d+(?:\.\d+)?)%")
_PAT_REPL = re.compile(r"(\d+(?:\.\d+)?)\s*%\s*(?:cement\s+)?replacement\s+by\s+(\w+)", re.I)

SCM_MATERIAL_COLUMNS = [
    "CaO",
    "SiO2",
    "Al2O3",
    "Fe2O3",
    "MgO",
    "SO3",
    "Na2O",
    "K2O",
    "TiO2",
    "LOI",
    "blaine_m2_kg",
    "d50_um",
    "bet_m2_g",
    "amorphous_pct",
]

DOR_TABLE_COLUMNS = [
    "obs_uid",
    "paper_doi",
    "mix_uid",
    "age_d",
    "dor_pct",
    "method_group",
    "confidence",
    "fig_only",
    "scm_role",
    "scm_pct",
    "resolve_how",
    "w_b",
    "curing_temp_C",
    "system_type",
] + [f"scm_{k}" for k in SCM_MATERIAL_COLUMNS]


def parse_binder_composition(binder_json: str | None, notes: str | None) -> dict[str, float]:
    """build_clean.py:66-74 — binder JSON with the notes-regex fallback."""
    try:
        bc = json.loads(binder_json or "{}")
    except Exception:  # noqa: BLE001
        bc = {}
    if not isinstance(bc, dict):
        bc = {}
    bc = {k: v for k, v in bc.items() if v is not None}
    if not bc and notes:
        props = _PAT_PROPS.findall(notes)
        if len(props) >= 2:
            bc = {k: float(v) for k, v in props}
        else:
            mrep = _PAT_REPL.search(notes)
            if mrep:
                bc = {mrep.group(2): float(mrep.group(1))}
    return bc


@dataclass
class ResolvedSCM:
    material: dict[str, Any] | None
    role: str | None
    scm_pct: float | None
    how: str


def resolve_scm(obs: Any, mix: Any, materials_of_paper: dict[str, dict[str, Any]]) -> ResolvedSCM:
    """build_clean.py:76-124 — the audited fallback chain.

    ``obs`` needs ``phase_name``; ``mix`` needs ``binder_composition_json``,
    ``notes``, ``name_in_paper``, ``primary_scm_role``, ``primary_scm_pct``,
    ``scm_total_pct``. Both accept ``sqlite3.Row`` or ``dict``.
    """
    bc = parse_binder_composition(mix["binder_composition_json"], mix["notes"])
    pm = materials_of_paper
    scm_mat, scm_pct_json, how, common_role = None, None, None, None
    cands_in_mix = [
        (pm[mid], pct)
        for mid, pct in bc.items()
        if mid in pm and pm[mid]["role"] not in CEMENT_ROLES and pm[mid]["role"] not in FILLER_ROLES
    ]
    if cands_in_mix:
        named = None
        if obs["phase_name"]:
            pn = obs["phase_name"].strip().lower()
            for m, pct in cands_in_mix:
                nm = (m["name_in_paper"] or "").strip().lower()
                if (nm and (nm == pn or nm in pn or pn in nm)) or m["role"] == role_from_text(obs["phase_name"]):
                    named = (m, pct)
                    break
        if named is not None:
            scm_mat, scm_pct_json = named
            how = "phase_match"
        else:
            cands_in_mix.sort(key=lambda x: -(x[1] or 0))
            scm_mat, scm_pct_json = cands_in_mix[0]
            how = "json"
    if scm_mat is None:
        cands = [m for m in pm.values() if m["role"] not in CEMENT_ROLES and m["role"] not in FILLER_ROLES]
        if len(cands) == 1:
            scm_mat = cands[0]
            how = "unique_mat"
        elif len(cands) > 1:
            ctx = " ".join(x for x in [mix["name_in_paper"], obs["phase_name"]] if x).lower()
            hits = [
                m
                for m in cands
                if m["name_in_paper"]
                and re.search(r"(?<![a-z0-9])" + re.escape(m["name_in_paper"].lower()) + r"(?![a-z0-9])", ctx)
            ]
            if not hits:
                hits = [
                    m
                    for m in cands
                    if m["material_id"]
                    and re.search(r"(?<![a-z0-9])" + re.escape(m["material_id"].lower()) + r"(?![a-z0-9])", ctx)
                ]
            if len(hits) == 1:
                scm_mat = hits[0]
                how = "name_match"
            elif len({m["role"] for m in cands}) == 1:
                common_role = cands[0]["role"]
                how = "common_role"

    scm_role = scm_mat["role"] if scm_mat else common_role
    if scm_role in (None, "other", "other_scm"):
        kw = role_from_text(obs["phase_name"], scm_mat["name_in_paper"] if scm_mat else None, mix["name_in_paper"])
        if kw:
            scm_role = kw
            how = how or "keyword_role"
    scm_role = scm_role or mix["primary_scm_role"]
    if how is None:
        how = "fail" if scm_role is None else "keyword_role"
    scm_pct = scm_pct_json if scm_pct_json is not None else (mix["primary_scm_pct"] or mix["scm_total_pct"])
    return ResolvedSCM(material=scm_mat, role=scm_role, scm_pct=scm_pct, how=how)


# Exactly the SQL of build_clean.py:53-60 (row order matters for the dedup step).
_DOR_SQL = """
SELECT o.obs_uid, o.paper_doi, o.mix_uid, o.age_d, o.value_norm, o.unit_norm, o.method,
       o.extraction_confidence, o.fig_only, o.phase_name,
       m.binder_composition_json, m.scm_total_pct, m.primary_scm_role, m.primary_scm_pct,
       m.w_b, m.curing_temp_C, m.notes, m.name_in_paper, m.curing_type
FROM observations o JOIN mixes m ON o.mix_uid=m.mix_uid
WHERE o.quantity='DoR_SCM' AND o.value_norm IS NOT NULL AND o.age_d IS NOT NULL
"""


def _materials_by_paper(con: sqlite3.Connection) -> dict[str, dict[str, dict[str, Any]]]:
    out: dict[str, dict[str, dict[str, Any]]] = {}
    for r in con.execute("SELECT * FROM materials"):
        out.setdefault(r["paper_doi"], {})[r["material_id"]] = dict(r)
    return out


def system_type_of(curing_type: str | None) -> str:
    return MODEL_SYSTEM if curing_type == MODEL_SYSTEM else "blended_cement"


def build_dor_table(con: sqlite3.Connection, *, with_stats: bool = False) -> pd.DataFrame | tuple[pd.DataFrame, dict[str, int]]:
    """Port of build_clean.py main loop → the 28-column ``dor_scm_final`` table."""
    if con.row_factory is not sqlite3.Row:
        con.row_factory = sqlite3.Row
    mats = _materials_by_paper(con)
    rows_out: list[dict[str, Any]] = []
    stats: dict[str, int] = {}
    for r in con.execute(_DOR_SQL):
        val, unit = r["value_norm"], (r["unit_norm"] or "").strip()
        if unit == "fraction":
            val *= 100.0
        if val < 0 or val > 100:
            continue
        res = resolve_scm(r, r, mats.get(r["paper_doi"], {}))
        stats[res.how] = stats.get(res.how, 0) + 1
        out = dict(
            obs_uid=r["obs_uid"],
            paper_doi=r["paper_doi"],
            mix_uid=r["mix_uid"],
            age_d=r["age_d"],
            dor_pct=val,
            method_group=method_group(r["method"]),
            confidence=r["extraction_confidence"] or "unknown",
            fig_only=r["fig_only"],
            scm_role=res.role,
            scm_pct=res.scm_pct,
            resolve_how=res.how,
            w_b=r["w_b"],
            curing_temp_C=r["curing_temp_C"],
            system_type=system_type_of(r["curing_type"]),
        )
        for k in SCM_MATERIAL_COLUMNS:
            out[f"scm_{k}"] = res.material[k] if res.material else None
        rows_out.append(out)
    df = pd.DataFrame(rows_out, columns=DOR_TABLE_COLUMNS)
    df = df.drop_duplicates(subset=["mix_uid", "age_d", "method_group", "dor_pct"]).reset_index(drop=True)
    if with_stats:
        return df, stats
    return df


def blended_only(df: pd.DataFrame) -> pd.DataFrame:
    return df[df["system_type"] != MODEL_SYSTEM].reset_index(drop=True)


# ---------------------------------------------------------------------------
# multitask_v6.py:22-29, 94-99 — feature derivation
# ---------------------------------------------------------------------------

CAT_FEATURES = ["scm_role", "method_group"]
BASE_FEATURES = ["log_age", "scm_pct", "w_b", "curing_temp_C"]
COMP_FEATURES = ["CaO_SiO2", "basicity", "pozz_sum", "Al_Si", "amorph", "fineness"]
GBM_FEATURES = BASE_FEATURES + COMP_FEATURES + CAT_FEATURES
BAYES_FEATURES = ["scm_pct", "w_b", "curing_temp_C", "CaO_SiO2"]

RANGES = {
    "DoR_SCM": (0, 100),
    "DoR_clinker": (0, 100),
    "CH_TGA": (0, 40),
    "bound_water": (0, 40),
    "cum_heat": (0, 600),
    # spec §4.2 additions
    "QXRD_phase": (0, 100),
    "chem_shrink": (0, 0.2),
    "CH_XRD": (0, 40),
}


def derive_composition_features(df: pd.DataFrame) -> pd.DataFrame:
    """Add ``log_age, CaO_SiO2, basicity, pozz_sum, Al_Si, amorph, fineness`` (multitask_v6.py:94-99)."""
    df = df.copy()
    if "age_d" in df.columns:
        df["log_age"] = np.log10(pd.to_numeric(df["age_d"], errors="coerce").clip(lower=0.01))
    C, S, A, F, M = [pd.to_numeric(df["scm_" + k], errors="coerce") for k in ["CaO", "SiO2", "Al2O3", "Fe2O3", "MgO"]]
    df["CaO_SiO2"] = C / S
    df["basicity"] = (C + M + A) / S
    df["pozz_sum"] = S + A + F
    df["Al_Si"] = A / S
    df["amorph"] = pd.to_numeric(df["scm_amorphous_pct"], errors="coerce")
    df["fineness"] = pd.to_numeric(df["scm_blaine_m2_kg"], errors="coerce")
    return df


_AUX_SQL = """
SELECT o.obs_uid, o.paper_doi, o.mix_uid, o.age_d, o.value_norm, o.unit_norm, o.unit_reported, o.quantity,
       o.method, o.method_detail, o.phase_name, o.basis_reported, o.uncertainty,
       o.source_locator, o.fig_only, o.extraction_confidence,
       m.binder_composition_json, m.scm_total_pct, m.primary_scm_role, m.primary_scm_pct,
       m.w_b, m.curing_temp_C, m.notes, m.name_in_paper, m.curing_type
FROM observations o JOIN mixes m ON o.mix_uid=m.mix_uid
WHERE o.value_norm IS NOT NULL AND o.age_d IS NOT NULL AND o.quantity IN ({placeholders})
"""


def build_aux_table(con: sqlite3.Connection, quantities: list[str] | tuple[str, ...]) -> pd.DataFrame:
    """Auxiliary observations (CH_TGA, CH_XRD, bound_water, QXRD_phase, chem_shrink,
    DoR_clinker, …) joined to mix features with the *same* resolver as the DoR table.
    ``model_system`` mixes are excluded and ``RANGES`` are applied (multitask_v6.py:59-65).
    """
    if con.row_factory is not sqlite3.Row:
        con.row_factory = sqlite3.Row
    qs = [q for q in quantities if q in RANGES]
    if not qs:
        return pd.DataFrame()
    mats = _materials_by_paper(con)
    placeholders = ",".join(["?"] * len(qs))
    rows: list[dict[str, Any]] = []
    for r in con.execute(_AUX_SQL.format(placeholders=placeholders), qs):
        if r["curing_type"] == MODEL_SYSTEM:
            continue
        val = r["value_norm"]
        if (r["unit_norm"] or "").strip() == "fraction" and r["quantity"].startswith("DoR"):
            val *= 100.0
        lo, hi = RANGES[r["quantity"]]
        if r["quantity"] == "chem_shrink" and (r["unit_norm"] or "").strip() != "mL/g binder":
            pass  # range check only meaningful in mL/g; '%' rows are grade D anyway (units.py)
        elif not (lo <= val <= hi):
            continue
        res = resolve_scm(r, r, mats.get(r["paper_doi"], {}))
        d = dict(
            obs_uid=r["obs_uid"],
            paper_doi=r["paper_doi"],
            mix_uid=r["mix_uid"],
            quantity=r["quantity"],
            phase_name=r["phase_name"],
            age_d=r["age_d"],
            value=val,
            unit_norm=r["unit_norm"],
            unit_reported=r["unit_reported"],
            basis_reported=r["basis_reported"],
            method=r["method"],
            method_detail=r["method_detail"],
            method_group=method_group(r["method"]),
            uncertainty=r["uncertainty"],
            source_locator=r["source_locator"],
            fig_only=r["fig_only"],
            extraction_confidence=r["extraction_confidence"],
            scm_role=res.role,
            scm_pct=res.scm_pct,
            resolve_how=res.how,
            scm_total_pct=r["scm_total_pct"],
            w_b=r["w_b"],
            curing_temp_C=r["curing_temp_C"],
            binder_composition_json=r["binder_composition_json"],
        )
        for k in SCM_MATERIAL_COLUMNS:
            d["scm_" + k] = res.material[k] if res.material else None
        rows.append(d)
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# New-SCM input → the same feature space (spec §4.2 ``scm_input_to_features``)
# ---------------------------------------------------------------------------


def _get(obj: Any, key: str, default: Any = None) -> Any:
    if obj is None:
        return default
    if isinstance(obj, dict):
        return obj.get(key, default)
    return getattr(obj, key, default)


def scm_input_to_features(scm: Any, mix: Any, age_d: float | None = None) -> dict[str, Any]:
    """Map an ``SCMSpec``/``MixSpec`` (pydantic models or plain dicts) to one
    feature row in the training feature space. Missing values stay ``NaN``;
    imputation is the bundle's job (``scaler.json``)."""
    ox = dict(_get(scm, "oxides", {}) or {})
    row: dict[str, Any] = {
        "scm_role": _get(scm, "role"),
        "scm_pct": float(_get(mix, "scm_pct")),
        "w_b": float(_get(mix, "w_b")),
        "curing_temp_C": _get(mix, "curing_temp_C", 20.0),
        "method_group": None,
    }
    for k in SCM_MATERIAL_COLUMNS:
        row["scm_" + k] = None
    for k in ["CaO", "SiO2", "Al2O3", "Fe2O3", "MgO", "SO3", "Na2O", "K2O", "TiO2", "LOI"]:
        if k in ox and ox[k] is not None:
            row["scm_" + k] = float(ox[k])
    for k in ["blaine_m2_kg", "d50_um", "bet_m2_g", "amorphous_pct"]:
        v = _get(scm, k)
        row["scm_" + k] = None if v is None else float(v)
    if age_d is not None:
        row["age_d"] = float(age_d)
    df = derive_composition_features(pd.DataFrame([row]))
    out = df.iloc[0].to_dict()
    # keep NaN (not None) for numeric; role/method as-is
    return out
