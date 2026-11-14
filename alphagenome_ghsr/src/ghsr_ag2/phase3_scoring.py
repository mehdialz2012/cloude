"""Phase III - AlphaGenome 2 variant-effect prediction.

For every variant configuration x every selected condition, requests the full
set of officially recommended AlphaGenome 2 variant scorers and stores the
complete per-track output.

Design decisions
----------------
* All four variant configurations (three rsIDs; rs34136054 contributes two
  alternate alleles) are scored identically. No variant is omitted.
* The complete officially recommended scorer set is requested for every
  variant/condition pair, rather than a hand-picked subset. This guarantees
  that no modality can be selectively omitted or selectively reported, and the
  pre-registered split between primary (mechanistically plausible) and
  secondary (completeness) modalities lives in `config.py`, fixed before the
  results were seen.
* Sequence context is 1 Mb (the largest supported), which is required for the
  gene-mask scorers to see GHSR and all neighbouring genes.
* Scores are stored exactly as returned by the API, at full track resolution,
  in `data/raw/phase3/`. Every downstream table is derived from that archive.
* `*_ACTIVE` scorers return absolute activity (max of REF and ALT) and are
  never treated as effect sizes; they are retained for QC only and flagged in
  the output via the `is_activity_scorer` column.

Run:  python -m ghsr_ag2.phase3_scoring
"""

from __future__ import annotations

import json
import time

import numpy as np
import pandas as pd
from alphagenome.data import genome
from alphagenome.models.v2 import dna_client, variant_scorers

from . import config, provenance

RAW = config.DATA_RAW / "phase3"
RAW.mkdir(parents=True, exist_ok=True)


def build_variant(spec: config.VariantSpec) -> genome.Variant:
    return genome.Variant(
        chromosome=config.CHROM,
        position=spec.position,
        reference_bases=spec.ref,
        alternate_bases=spec.alt,
        name=f"{spec.rsid}_{spec.ref}>{spec.alt}",
    )


def recommended_scorers() -> list:
    """Officially recommended human scorers, excluding cell-conditioning-only."""
    return variant_scorers.get_recommended_scorers(
        organism=dna_client.Organism.HOMO_SAPIENS.to_proto(),
        cell_conditioning=False,
    )


def score_one(
    client: dna_client.DnaClient,
    spec: config.VariantSpec,
    condition: config.ConditionSpec,
    sequence_length_key: str,
    scorers: list,
) -> pd.DataFrame:
    variant = build_variant(spec)
    interval = variant.reference_interval.resize(
        dna_client.SUPPORTED_SEQUENCE_LENGTHS[sequence_length_key]
    )
    outs = client.score_variant(
        interval=interval,
        variant=variant,
        condition=condition.curie,
        variant_scorers=scorers,
    )
    tidy = variant_scorers.tidy_scores(outs, match_gene_strand=True)
    if tidy is None or len(tidy) == 0:
        return pd.DataFrame()
    tidy = tidy.copy()
    tidy["variant_key"] = spec.key
    tidy["rsid"] = spec.rsid
    tidy["manuscript_label"] = spec.manuscript_label
    tidy["is_primary_configuration"] = spec.is_primary
    tidy["ref_plus_strand"] = spec.ref
    tidy["alt_plus_strand"] = spec.alt
    tidy["condition_curie"] = condition.curie
    tidy["condition_label"] = condition.label
    tidy["panel"] = condition.panel
    tidy["target_priority"] = condition.target_priority
    tidy["surrogate_for"] = condition.surrogate_for or ""
    tidy["sequence_length_key"] = sequence_length_key
    tidy["sequence_length_bp"] = dna_client.SUPPORTED_SEQUENCE_LENGTHS[
        sequence_length_key
    ]
    return tidy


def make_parquet_safe(df: pd.DataFrame) -> pd.DataFrame:
    """Stringify non-scalar object columns so the frame can be serialised.

    `tidy_scores` returns rich objects in some columns (e.g. `variant_id`
    holds a `genome.Variant`, `scored_interval` a `genome.Interval`). Their
    string forms are the canonical, lossless textual representations used
    throughout the AlphaGenome documentation, so converting them preserves all
    information needed to reconstruct the request.
    """
    df = df.copy()
    for col in df.columns:
        if df[col].dtype != object:
            continue
        sample = df[col].dropna()
        if sample.empty:
            continue
        first = sample.iloc[0]
        if isinstance(first, (str, bool, int, float, np.integer, np.floating)):
            continue
        df[col] = df[col].map(lambda v: None if v is None else str(v))
    return df


def annotate_scorer_classes(df: pd.DataFrame) -> pd.DataFrame:
    """Label each scorer as an effect-size or an activity (QC) readout."""
    rec = variant_scorers.RECOMMENDED_VARIANT_SCORERS
    name_to_key = {v.name: k for k, v in rec.items()}
    df = df.copy()
    df["scorer_key"] = df["variant_scorer"].map(name_to_key)
    df["is_activity_scorer"] = df["scorer_key"].isin(
        config.ACTIVITY_SCORER_KEYS
    )
    df["is_signed_score"] = ~df["is_activity_scorer"]
    df["modality_class"] = "unclassified"
    df.loc[
        df["scorer_key"].isin(config.PRIMARY_SCORER_KEYS), "modality_class"
    ] = "primary"
    df.loc[
        df["scorer_key"].isin(config.SECONDARY_SCORER_KEYS), "modality_class"
    ] = "secondary_completeness"
    df.loc[df["is_activity_scorer"], "modality_class"] = "activity_qc"
    return df


def run(
    sequence_length_key: str = config.SEQUENCE_LENGTH_KEY,
    tag: str = "main",
) -> pd.DataFrame:
    client = dna_client.create(config.api_key(), timeout=120)
    scorers = recommended_scorers()
    print(
        f"[phase3:{tag}] {len(scorers)} recommended scorers x "
        f"{len(config.VARIANTS)} variant configurations x "
        f"{len(config.CONDITIONS)} conditions @ {sequence_length_key}"
    )
    frames: list[pd.DataFrame] = []
    failures: list[dict] = []
    t0 = time.time()
    for spec in config.VARIANTS:
        for cond in config.CONDITIONS:
            try:
                frames.append(
                    score_one(client, spec, cond, sequence_length_key, scorers)
                )
                print(
                    f"[phase3:{tag}] ok  {spec.key:<18} {cond.label:<26} "
                    f"({time.time() - t0:6.1f}s)"
                )
            except Exception as exc:  # noqa: BLE001 - recorded, not hidden
                failures.append(
                    {
                        "variant_key": spec.key,
                        "condition_curie": cond.curie,
                        "condition_label": cond.label,
                        "sequence_length_key": sequence_length_key,
                        "error_type": type(exc).__name__,
                        "error": str(exc)[:600],
                    }
                )
                print(f"[phase3:{tag}] FAIL {spec.key} {cond.label}: {exc}")

    if not frames:
        raise RuntimeError("no scores were returned; refusing to continue")
    df = annotate_scorer_classes(
        make_parquet_safe(pd.concat(frames, ignore_index=True))
    )

    out = RAW / f"phase3_raw_scores_{tag}.parquet"
    df.to_parquet(out, index=False)
    print(f"[phase3:{tag}] wrote {len(df):,} raw track-level scores -> {out}")

    fail_path = RAW / f"phase3_failures_{tag}.json"
    fail_path.write_text(json.dumps(failures, indent=2))
    if failures:
        print(f"[phase3:{tag}] {len(failures)} failed request(s); see {fail_path}")

    provenance.write(
        f"phase3_{tag}",
        {
            "sequence_length_key": sequence_length_key,
            "n_variant_configurations": len(config.VARIANTS),
            "n_conditions": len(config.CONDITIONS),
            "scorers_requested": [s.name for s in scorers],
            "n_scorers": len(scorers),
            "n_rows": int(len(df)),
            "n_failed_requests": len(failures),
            "raw_output": str(out),
            "elapsed_seconds": round(time.time() - t0, 1),
            "data_sent_to_api": (
                "public GRCh38 variant coordinates and public ontology CURIEs "
                "only; no participant-level or clinical data"
            ),
        },
    )
    return df


def main() -> None:
    print("Phase III: AlphaGenome 2 variant-effect prediction")
    df = run()
    # Quick orientation printout: GHSR expression effect in each condition.
    ghsr = df[
        (df.gene_name == config.GHSR_SYMBOL)
        & (df.scorer_key == "RNA_SEQ")
        & df.is_primary_configuration
    ]
    if len(ghsr):
        summary = (
            ghsr.groupby(["rsid", "condition_label"])["raw_score"]
            .agg(["mean", "min", "max", "count"])
            .reset_index()
        )
        print("\nGHSR RNA-seq log-fold-change (ALT - REF), per condition:")
        print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
