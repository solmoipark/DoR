"""Build tests/fixtures/mini_scm_dor.sqlite (5-paper subset, writable, separate file)
and the golden CSV subset tests/golden/dor_scm_final_subset.csv.

Usage: python scripts/make_fixture.py [full_db] [dor_scm_final.csv]
The five papers are chosen to exercise every resolver fallback
(json, phase_match, name_match, common_role, keyword_role, unique_mat) and
the model_system exclusion.
"""

from __future__ import annotations

import sqlite3
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from dorgems import config  # noqa: E402

FIXTURE_PAPERS = [
    "10.1016/j.cemconres.2014.06.011",  # keyword_role + phase_match, metakaolin/slag, model_system mixes
    "10.1016/S0008-8846(03)00213-8",  # common_role, fly_ash
    "10.1016/j.cemconres.2019.04.015",  # name_match, fly_ash
    "10.1016/S0008-8846(01)00581-6",  # phase_match + unique_mat, fly_ash
    "10.1016/j.cemconcomp.2012.02.004",  # json; fly_ash, slag, steel_slag
]


def main() -> None:
    md = config.modeling_dir()
    full_db = Path(sys.argv[1]) if len(sys.argv) > 1 else md / "scm_dor_enriched.db"
    final_csv = Path(sys.argv[2]) if len(sys.argv) > 2 else md / "dor_scm_final.csv"
    root = Path(__file__).resolve().parents[1]
    out_db = root / "tests" / "fixtures" / "mini_scm_dor.sqlite"
    out_csv = root / "tests" / "golden" / "dor_scm_final_subset.csv"
    out_db.parent.mkdir(parents=True, exist_ok=True)
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    if out_db.exists():
        out_db.unlink()

    src = sqlite3.connect(f"{full_db.resolve().as_uri()}?mode=ro", uri=True)
    dst = sqlite3.connect(out_db)
    for (sql,) in src.execute("SELECT sql FROM sqlite_master WHERE type='table' AND sql IS NOT NULL"):
        dst.execute(sql)
    ph = ",".join(["?"] * len(FIXTURE_PAPERS))
    for table, col in [("papers", "doi"), ("materials", "paper_doi"), ("mixes", "paper_doi"), ("observations", "paper_doi")]:
        rows = src.execute(f"SELECT * FROM {table} WHERE {col} IN ({ph})", FIXTURE_PAPERS).fetchall()
        ncol = len(src.execute(f"SELECT * FROM {table} LIMIT 1").description)
        dst.executemany(f"INSERT INTO {table} VALUES ({','.join(['?'] * ncol)})", rows)
        print(f"{table}: {len(rows)} rows")
    dst.commit()
    dst.close()

    df = pd.read_csv(final_csv)
    sub = df[df["paper_doi"].isin(FIXTURE_PAPERS)]
    sub.to_csv(out_csv, index=False)
    print(f"golden subset: {len(sub)} rows, hows={sorted(sub['resolve_how'].unique())}, "
          f"system_type={sub['system_type'].value_counts().to_dict()}")


if __name__ == "__main__":
    main()
