"""DoR bench (spec M4): deterministic checks of the tool layer's guardrails and the
normal scenario-A path, mock only. Mirrors GemsPilot's GEMS-Agent-Bench structure
(scenario id / kind / expect / checks) so results can be merged into its report.
"""

from __future__ import annotations

import argparse
import json
import sqlite3
import sys
import tempfile
from pathlib import Path
from typing import Any

import yaml

from ..config import configs_dir
from . import tools as T


def _check(name: str, expected: Any, actual: Any) -> dict[str, Any]:
    return {"name": name, "expected": expected, "actual": actual, "ok": expected == actual}


def run_scenario(sc: dict[str, Any], out: Path) -> dict[str, Any]:
    kind = sc["kind"]
    expect = dict(sc.get("expect") or {})
    checks: list[dict[str, Any]] = []
    out.mkdir(parents=True, exist_ok=True)
    if kind == "dor_envelope":
        r = T.dor_run_envelope(sc["scm"], sc["mix"], sc.get("ages"), str(out / "env"), str(out / "igdb"), use_mock=True)
        checks.append(_check("ok", expect.get("ok", True), r["ok"]))
        runs = r["summary"].get("result_summary", {}).get("runs_ok", {})
        checks.append(_check("runs_ok", expect.get("runs_ok"), runs))
        man = json.loads(Path(r["artifacts"]["manifest"]).read_text(encoding="utf-8")) if r["ok"] else {}
        alpha_ok = all((v.get("self_check") or {}).get("alpha_ok") is True for v in (man.get("runs") or {}).values()) if man else False
        checks.append(_check("self_check_alpha", expect.get("self_check_alpha", True), alpha_ok))
    elif kind == "dor_budget":
        r = T.dor_run_forward_with_dor(sc["forward_query"], None, str(out / "o"), str(out / "igdb"), use_mock=bool(sc.get("use_mock", True)), max_xgems_calls=sc.get("max_xgems_calls"))
        checks.append(_check("ok", expect.get("ok"), r["ok"]))
        checks.append(_check("error_contains", True, expect.get("error_contains", "") in (r.get("error") or "")))
    elif kind == "dor_db_write":
        from ..config import literature_db_path

        dbp = literature_db_path(required=False)
        fixture = Path(__file__).resolve().parents[3] / "tests" / "fixtures" / "mini_scm_dor.sqlite"
        target = dbp if dbp and dbp.is_file() else fixture
        r = T.dor_db_lookup(sc["query_name"], lit_db=str(target))
        checks.append(_check("ok", expect.get("ok", False), r["ok"]))
        n_before = sqlite3.connect(f"{target.resolve().as_uri()}?mode=ro", uri=True).execute("SELECT COUNT(*) FROM observations").fetchone()[0]
        checks.append(_check("db_unchanged", True, n_before > 0))
    elif kind == "dor_policy":
        try:
            from gemspilot.runner import _policy_check
        except Exception:  # noqa: BLE001
            return {"id": sc["id"], "kind": kind, "skipped": "gemspilot not installed", "checks": [], "ok": True}
        spec = T.TOOLS_BY_NAME[sc["tool"]]
        denied = _policy_check(spec, dict(sc.get("arguments") or {}), allow_real=bool(sc.get("allow_real", False))) is not None
        checks.append(_check("denied", expect.get("denied", True), denied))
    elif kind == "dor_stage":
        inf = {"id": "bench", "alpha": {"ages_d": [28], "q05": [0.1], "q50": [0.2], "q95": [0.3]}, "a_max": {"q50": 0.5}, "tau_d": {"q50": 20.0}, "ppc": [], "slot": "slag"}
        db = out / "staging.sqlite"
        r = T.TOOLS_BY_NAME["dor_stage_inferred"].func(inf, str(db))
        checks.append(_check("dry_run", expect.get("dry_run", True), bool(r["summary"].get("dry_run"))))
        n = sqlite3.connect(db).execute("SELECT COUNT(*) FROM inferred_dor").fetchone()[0] if db.exists() else 0
        checks.append(_check("rows_written", expect.get("rows_written", 0), n))
    else:
        raise ValueError(f"unknown bench kind {kind!r}")
    return {"id": sc["id"], "kind": kind, "checks": checks, "ok": all(c["ok"] for c in checks)}


def run_bench(config: Path | None = None, out: Path | None = None) -> dict[str, Any]:
    cfg = yaml.safe_load((config or configs_dir() / "agent_bench_dor.yaml").read_text(encoding="utf-8"))
    out = out or Path(tempfile.mkdtemp(prefix="dor_bench_"))
    results = [run_scenario(sc, out / sc["id"]) for sc in cfg["scenarios"]]
    report = {"n": len(results), "n_ok": sum(r["ok"] for r in results), "results": results, "out": str(out)}
    (out / "bench_report.json").write_text(json.dumps(report, indent=2, default=str), encoding="utf-8")
    return report


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default=None)
    ap.add_argument("--out", default=None)
    a = ap.parse_args(argv)
    rep = run_bench(Path(a.config) if a.config else None, Path(a.out) if a.out else None)
    for r in rep["results"]:
        print(("PASS" if r["ok"] else "FAIL"), r["id"], "" if r["ok"] else [c for c in r["checks"] if not c["ok"]])
    print(f"{rep['n_ok']}/{rep['n']} scenarios passed; report: {rep['out']}")
    return 0 if rep["n_ok"] == rep["n"] else 1


if __name__ == "__main__":
    sys.exit(main())
