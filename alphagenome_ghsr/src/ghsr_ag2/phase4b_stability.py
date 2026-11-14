"""Phase IV-b - cross-sequence-context stability.

The headline analysis uses a 1 Mb sequence context. Phase IV repeats every
prediction at 16 kb, 100 kb and 500 kb. This module asks the question that
matters for interpretation: **does the sign of the predicted effect survive a
change of sequence context?**

That is the sharpest available internal test of whether an effect of this
magnitude is meaningful. A genuine cis-regulatory effect driven by the local
sequence should not reverse direction when flanking sequence far from the
variant is added or removed. An effect that does reverse is operating inside
the model's own context-dependent numerical variability and should not be
interpreted as a direction of effect at all.

Outputs:
  table_S9_context_stability_ghsr_expression.csv
  table_S10_context_stability_observed_ctcf.csv
  data/processed/phase4b_stability_summary.json

Run:  python -m ghsr_ag2.phase4b_stability
"""

from __future__ import annotations

import json

import numpy as np
import pandas as pd

from . import config, provenance

CONTEXTS = (
    "SEQUENCE_LENGTH_16KB",
    "SEQUENCE_LENGTH_100KB",
    "SEQUENCE_LENGTH_500KB",
    "SEQUENCE_LENGTH_1MB",
)
HEADLINE = "SEQUENCE_LENGTH_1MB"


def load_all_contexts() -> pd.DataFrame:
    raw = config.DATA_RAW / "phase3"
    paths = {
        "SEQUENCE_LENGTH_1MB": raw / "phase3_raw_scores_main.parquet",
        "SEQUENCE_LENGTH_16KB": raw
        / "phase3_raw_scores_SEQUENCE_LENGTH_16KB.parquet",
        "SEQUENCE_LENGTH_100KB": raw
        / "phase3_raw_scores_SEQUENCE_LENGTH_100KB.parquet",
        "SEQUENCE_LENGTH_500KB": raw
        / "phase3_raw_scores_SEQUENCE_LENGTH_500KB.parquet",
    }
    missing = [k for k, p in paths.items() if not p.exists()]
    if missing:
        raise RuntimeError(
            "missing sequence-context runs: "
            f"{missing}. Run `python -m ghsr_ag2.phase4_qc --seqlen` first."
        )
    return pd.concat(
        [pd.read_parquet(p) for p in paths.values()], ignore_index=True
    )


def _stability(pivot: pd.DataFrame) -> pd.DataFrame:
    cols = [c for c in CONTEXTS if c in pivot.columns]
    out = pivot[cols].copy()
    vals = out.to_numpy(dtype=float)
    signs = np.sign(vals)
    out["n_contexts"] = len(cols)
    out["sign_consistent_across_contexts"] = (
        np.abs(np.nansum(signs, axis=1)) == len(cols)
    )
    out["all_negative"] = (vals < 0).all(axis=1)
    out["all_positive"] = (vals > 0).all(axis=1)
    out["min_across_contexts"] = np.nanmin(vals, axis=1)
    out["max_across_contexts"] = np.nanmax(vals, axis=1)
    out["range_across_contexts"] = (
        out["max_across_contexts"] - out["min_across_contexts"]
    )
    out["max_abs_across_contexts"] = np.nanmax(np.abs(vals), axis=1)
    out["headline_1mb_value"] = (
        out[HEADLINE] if HEADLINE in out.columns else np.nan
    )
    return out.reset_index()


def ghsr_expression_stability(df: pd.DataFrame) -> pd.DataFrame:
    sel = df[
        (df.scorer_key == "RNA_SEQ")
        & (df.gene_name == config.GHSR_SYMBOL)
        & df.is_primary_configuration
    ]
    pivot = sel.pivot_table(
        index=["rsid", "condition_label", "panel"],
        columns="sequence_length_key",
        values="raw_score",
        aggfunc="mean",
    )
    return _stability(pivot)


def observed_ctcf_stability(df: pd.DataFrame) -> pd.DataFrame:
    sel = df[
        (df.scorer_key == "CHIP_TF")
        & (df.imputation_mode == "observed")
        & (df.transcription_factor == "CTCF")
        & df.is_primary_configuration
    ]
    pivot = sel.pivot_table(
        index=["rsid", "condition_label"],
        columns="sequence_length_key",
        values="raw_score",
        aggfunc="mean",
    )
    return _stability(pivot)


def summarise(name: str, tbl: pd.DataFrame) -> dict[str, object]:
    by_rsid = {}
    for rsid, g in tbl.groupby("rsid"):
        by_rsid[rsid] = {
            "n_cells": int(len(g)),
            "n_sign_consistent": int(g.sign_consistent_across_contexts.sum()),
            "fraction_sign_consistent": round(
                float(g.sign_consistent_across_contexts.mean()), 3
            ),
            "n_negative_in_all_contexts": int(g.all_negative.sum()),
            "n_positive_in_all_contexts": int(g.all_positive.sum()),
            "median_range_across_contexts": round(
                float(g.range_across_contexts.median()), 6
            ),
            "max_abs_any_context": round(
                float(g.max_abs_across_contexts.max()), 6
            ),
        }
    return {
        "readout": name,
        "n_cells_total": int(len(tbl)),
        "n_sign_consistent_total": int(
            tbl.sign_consistent_across_contexts.sum()
        ),
        "fraction_sign_consistent_total": round(
            float(tbl.sign_consistent_across_contexts.mean()), 3
        ),
        "by_rsid": by_rsid,
    }


def main() -> None:
    print("Phase IV-b: cross-sequence-context stability")
    df = load_all_contexts()

    ghsr = ghsr_expression_stability(df)
    ctcf = observed_ctcf_stability(df)

    p1 = config.TABLES / "table_S9_context_stability_ghsr_expression.csv"
    p2 = config.TABLES / "table_S10_context_stability_observed_ctcf.csv"
    ghsr.to_csv(p1, index=False)
    ctcf.to_csv(p2, index=False)

    s1 = summarise("GHSR RNA-seq log-fold-change", ghsr)
    s2 = summarise("observed CTCF ChIP log2 fold change", ctcf)
    (config.DATA_PROCESSED / "phase4b_stability_summary.json").write_text(
        json.dumps({"ghsr_expression": s1, "observed_ctcf": s2}, indent=2)
    )

    for s in (s1, s2):
        print(f"\n--- {s['readout']} ---")
        print(
            f"sign-consistent across all 4 contexts: "
            f"{s['n_sign_consistent_total']}/{s['n_cells_total']} "
            f"({100 * s['fraction_sign_consistent_total']:.0f}%)"
        )
        for rsid, v in s["by_rsid"].items():
            print(
                f"  {rsid:<13} {v['n_sign_consistent']}/{v['n_cells']} "
                f"consistent  (neg in all: {v['n_negative_in_all_contexts']}, "
                f"pos in all: {v['n_positive_in_all_contexts']}, "
                f"max|effect| {v['max_abs_any_context']:.4f})"
            )
    print(f"\nWrote {p1}\nWrote {p2}")

    provenance.write(
        "phase4b",
        {
            "contexts_compared": list(CONTEXTS),
            "headline_context": HEADLINE,
            "ghsr_expression_summary": s1,
            "observed_ctcf_summary": s2,
            "rationale": (
                "sign reversal under a change of sequence context indicates "
                "that an effect of this magnitude lies within the model's own "
                "context-dependent variability and carries no interpretable "
                "direction"
            ),
        },
    )


if __name__ == "__main__":
    main()
