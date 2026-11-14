#!/usr/bin/env python3
"""Export the raw and processed model outputs as CSV/TSV.

The canonical archive of the AlphaGenome 2 output is Parquet, which preserves
dtypes exactly and is ~7x smaller. This script additionally writes CSV copies
for tools that cannot read Parquet. CSVs are gzip-compressed because the main
table has 162,144 rows.

Run:  python scripts/export_csv.py
"""

from __future__ import annotations

import pathlib
import sys

import pandas as pd

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "src"))

from ghsr_ag2 import config  # noqa: E402

EXPORTS = {
    "phase3/phase3_raw_scores_main.parquet": (
        "raw_scores_main_1Mb.csv.gz",
        "Raw per-track AlphaGenome 2 variant scores, 1 Mb context "
        "(the headline analysis).",
    ),
    "phase3/phase3_raw_scores_SEQUENCE_LENGTH_16KB.parquet": (
        "raw_scores_16kb.csv.gz",
        "Raw per-track scores at 16 kb sequence context.",
    ),
    "phase3/phase3_raw_scores_SEQUENCE_LENGTH_100KB.parquet": (
        "raw_scores_100kb.csv.gz",
        "Raw per-track scores at 131 kb sequence context.",
    ),
    "phase3/phase3_raw_scores_SEQUENCE_LENGTH_500KB.parquet": (
        "raw_scores_500kb.csv.gz",
        "Raw per-track scores at 524 kb sequence context.",
    ),
    "phase4/phase4_ism_reference_distribution.parquet": (
        "ism_reference_distribution.csv.gz",
        "Raw per-track scores for all 1,197 single-base substitutions across "
        "the sequenced amplicon, in four conditions.",
    ),
}


def main() -> None:
    out_dir = config.DATA_PROCESSED / "csv_exports"
    out_dir.mkdir(parents=True, exist_ok=True)
    manifest = []
    for rel, (name, desc) in EXPORTS.items():
        src = config.DATA_RAW / rel
        if not src.exists():
            print(f"[skip] {src} not present")
            continue
        df = pd.read_parquet(src)
        dst = out_dir / name
        df.to_csv(dst, index=False, compression="gzip")
        manifest.append(
            {
                "file": name,
                "rows": len(df),
                "columns": df.shape[1],
                "size_mb": round(dst.stat().st_size / 1e6, 2),
                "source_parquet": rel,
                "description": desc,
            }
        )
        print(f"[csv] {name}: {len(df):,} rows, {dst.stat().st_size / 1e6:.1f} MB")

    pd.DataFrame(manifest).to_csv(out_dir / "MANIFEST.csv", index=False)
    print(f"\nWrote {len(manifest)} CSV exports to {out_dir}")
    print(
        "Note: results/tables/*.csv already holds every derived/summary table "
        "in plain CSV."
    )


if __name__ == "__main__":
    main()
