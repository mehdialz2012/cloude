"""Phase IV - robustness, quality control and interpretive context.

This module does not re-interpret Phase III; it stress-tests it.

QC1  REF/ALT integrity. Re-verifies that every scored REF base matches the
     GRCh38 reference sequence, that the variant identity sent to the API is
     the one intended, and that all variant configurations were scored in
     every condition with no silent omissions.

QC2  Local null / sensitivity context. The central interpretive problem is
     that an AlphaGenome score of, say, -0.0002 has no intrinsic meaning. In
     silico mutagenesis (ISM) over the sequenced amplicon generates the
     distribution of effects achievable by *any* single-base substitution in
     that same window, in the same condition, under the same scorer. The three
     observed SNPs are then placed as percentiles of that empirical
     distribution. This is an empirical reference distribution, not a
     statistical test, and no P-value is derived from it.

     The same ISM run doubles as a positive control for model sensitivity: if
     some positions in the amplicon produce large predicted effects while the
     studied SNPs do not, then a near-zero score is an informative negative
     rather than a silent pipeline failure.

QC3  Sequence-context sensitivity. The headline analysis uses a 1 Mb context.
     Scores are recomputed at 16 kb, 100 kb and 500 kb to confirm that
     conclusions are not an artefact of context length.

QC4  Dependence on imputed tracks. Every summary is recomputed using only
     `observed` tracks, so that no conclusion rests on imputed-exploratory
     data without that being visible.

QC5  Aggregated-score / raw-track consistency. Confirms that the aggregated
     centre-mask scores from Phase III agree in sign and order of magnitude
     with the base-resolution REF/ALT tracks from Phase III-b.

Run:  python -m ghsr_ag2.phase4_qc  [--ism]
"""

from __future__ import annotations

import argparse
import json

import numpy as np
import pandas as pd
from alphagenome.data import genome
from alphagenome.models.v2 import dna_client, variant_scorers

from . import config, phase1_annotation, phase3_scoring, provenance

RAW = config.DATA_RAW / "phase4"
RAW.mkdir(parents=True, exist_ok=True)

# Conditions and scorers used for the ISM reference distribution. Chosen
# before inspecting effect sizes: the two highest-GHSR brain/endocrine
# conditions, the best-supported bone condition and the fully observed
# reference tissue.
ISM_CONDITIONS = (
    "UBERON:0000007",             # pituitary gland
    "UBERON:0001898",             # hypothalamus
    "UBERON:0002481_CL:0000062",  # osteoblast
    "UBERON:0000945",             # stomach (all output types observed)
)
ISM_SCORER_KEYS = ("RNA_SEQ", "CAGE", "DNASE", "ATAC")


# --------------------------------------------------------------------------
# QC1 - REF/ALT integrity
# --------------------------------------------------------------------------
def qc1_ref_alt_integrity(scores: pd.DataFrame) -> dict[str, object]:
    seq_start, seq_end = 172_449_100, 172_449_700
    seq = phase1_annotation.fetch_sequence(
        seq_start, seq_end, f"{seq_start}_{seq_end}"
    )
    per_variant = []
    for spec in config.VARIANTS:
        observed_ref = seq[spec.position - seq_start]
        variant_obj = phase3_scoring.build_variant(spec)
        subset = scores[scores.variant_key == spec.key]
        per_variant.append(
            {
                "variant_key": spec.key,
                "position": spec.position,
                "ref_declared": spec.ref,
                "ref_in_grch38": observed_ref,
                "ref_matches_grch38": observed_ref == spec.ref,
                "variant_string_sent_to_api": str(variant_obj),
                "variant_id_returned_by_api": sorted(
                    set(subset.get("variant_id", pd.Series(dtype=str)))
                ),
                "n_conditions_scored": int(subset.condition_curie.nunique()),
                "n_conditions_expected": len(config.CONDITIONS),
                "n_scorers_returned": int(subset.variant_scorer.nunique()),
                "n_rows": int(len(subset)),
            }
        )
    df = pd.DataFrame(per_variant)
    return {
        "per_variant": df.to_dict("records"),
        "all_ref_match_grch38": bool(df.ref_matches_grch38.all()),
        "all_variants_scored_in_all_conditions": bool(
            (df.n_conditions_scored == df.n_conditions_expected).all()
        ),
        "identical_scorer_count_across_variants": bool(
            df.n_scorers_returned.nunique() == 1
        ),
    }


# --------------------------------------------------------------------------
# QC2 - ISM reference distribution
# --------------------------------------------------------------------------
def qc2_ism_reference_distribution(
    ism_start: int = config.AMPLICON_START,
    ism_end: int = config.AMPLICON_END,
) -> pd.DataFrame:
    """ISM across the sequenced amplicon; returns one row per ISM variant."""
    client = dna_client.create(config.api_key(), timeout=300)
    scorers = [
        variant_scorers.RECOMMENDED_VARIANT_SCORERS[k] for k in ISM_SCORER_KEYS
    ]
    # Interval must be a supported sequence length centred on the amplicon.
    centre = (ism_start + ism_end) // 2
    interval = genome.Interval(
        config.CHROM, centre - 1, centre
    ).resize(dna_client.SUPPORTED_SEQUENCE_LENGTHS["SEQUENCE_LENGTH_100KB"])
    # genome.Interval is 0-based half-open; convert the 1-based inclusive
    # amplicon coordinates accordingly.
    ism_interval = genome.Interval(config.CHROM, ism_start - 1, ism_end)
    print(
        f"[phase4] ISM over {ism_interval} "
        f"({ism_interval.width} bp) x {len(ISM_CONDITIONS)} conditions"
    )

    # `score_ism_variants` does not stamp `uns['variant_scorer']` on its
    # outputs (unlike `score_variant`), so it is attached here by request
    # order. The assertion below independently verifies the assignment
    # against the track names returned by the server, so a change in ordering
    # would fail loudly rather than silently mislabel modalities.
    expected_assay_substring = {
        "RNA_SEQ": "RNA-seq",
        "CAGE": "CAGE",
        "DNASE": "DNase-seq",
        "ATAC": "ATAC-seq",
    }

    def _attach_scorers(adatas: list) -> list:
        if len(adatas) != len(scorers):
            raise RuntimeError(
                f"expected {len(scorers)} ISM score objects per variant, "
                f"got {len(adatas)}"
            )
        for adata, key, scorer in zip(adatas, ISM_SCORER_KEYS, scorers):
            names = " ".join(adata.var["name"].astype(str))
            needle = expected_assay_substring[key]
            if needle.lower() not in names.lower():
                raise RuntimeError(
                    "ISM scorer/output mismatch: expected tracks containing "
                    f"{needle!r} for scorer {key}, got {names[:200]!r}. "
                    "Refusing to proceed with possibly mislabelled modalities."
                )
            adata.uns["variant_scorer"] = scorer
        return adatas

    frames = []
    for curie in ISM_CONDITIONS:
        label = config.CONDITION_BY_CURIE[curie].label
        per_variant_scores = client.score_ism_variants(
            interval=interval,
            ism_interval=ism_interval,
            condition=curie,
            variant_scorers=scorers,
            progress_bar=False,
        )
        for adatas in per_variant_scores:
            tidy = variant_scorers.tidy_scores(
                _attach_scorers(list(adatas)), match_gene_strand=True
            )
            if tidy is None or not len(tidy):
                continue
            tidy = tidy.copy()
            tidy["condition_curie"] = curie
            tidy["condition_label"] = label
            frames.append(tidy)
        print(f"[phase4] ISM done: {label} ({len(per_variant_scores)} variants)")

    df = phase3_scoring.annotate_scorer_classes(
        phase3_scoring.make_parquet_safe(pd.concat(frames, ignore_index=True))
    )
    out = RAW / "phase4_ism_reference_distribution.parquet"
    df.to_parquet(out, index=False)
    print(f"[phase4] wrote {len(df):,} ISM score rows -> {out}")
    return df


def qc2_percentiles(ism: pd.DataFrame) -> pd.DataFrame:
    """Place each observed SNP inside its own local ISM distribution.

    Because the ISM interval spans the whole sequenced amplicon, and ISM
    enumerates all three possible substitutions at every position in that
    interval, the three studied SNPs are themselves members of the ISM set.
    Their scores are therefore read out of the same run, under an identical
    sequence context, scorer and condition as the reference distribution they
    are compared against - no cross-run normalisation is involved.
    """
    rows = []
    endpoints = [
        ("RNA_SEQ", config.GHSR_SYMBOL),
        ("CAGE", None),
        ("DNASE", None),
        ("ATAC", None),
    ]
    # Canonical variant_id strings for the studied configurations.
    wanted = {
        f"{config.CHROM}:{s.position}:{s.ref}>{s.alt}": s for s in config.VARIANTS
    }
    for curie in ISM_CONDITIONS:
        label = config.CONDITION_BY_CURIE[curie].label
        for scorer_key, gene in endpoints:
            m = (ism.scorer_key == scorer_key) & (ism.condition_curie == curie)
            if gene is not None:
                m &= ism.gene_name == gene
            sel = ism[m]
            if not len(sel):
                continue
            # Collapse tracks to one value per ISM variant (mean across
            # tracks), applied identically to every variant in the set.
            by_variant = sel.groupby("variant_id")["raw_score"].mean()
            vals = by_variant.to_numpy()
            abs_vals = np.abs(vals)
            for vid, spec in wanted.items():
                if vid not in by_variant.index:
                    continue
                value = float(by_variant.loc[vid])
                pct_abs = float((abs_vals <= abs(value)).mean() * 100)
                rows.append(
                    {
                        "condition_curie": curie,
                        "condition_label": label,
                        "scorer_key": scorer_key,
                        "gene": gene or "(variant-centred window)",
                        "variant_key": spec.key,
                        "rsid": spec.rsid,
                        "variant_id": vid,
                        "observed_score": value,
                        "ism_n_variants": int(vals.size),
                        "ism_median_abs": float(np.median(abs_vals)),
                        "ism_p95_abs": float(np.percentile(abs_vals, 95)),
                        "ism_max_abs": float(abs_vals.max()),
                        "percentile_of_abs_effect_in_amplicon": pct_abs,
                        "ratio_to_ism_max_abs": (
                            abs(value) / float(abs_vals.max())
                            if abs_vals.max() > 0
                            else np.nan
                        ),
                    }
                )
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------
# QC3 - sequence-context sensitivity
# --------------------------------------------------------------------------
def qc3_sequence_length_sensitivity() -> pd.DataFrame:
    frames = []
    for key in config.SENSITIVITY_SEQUENCE_LENGTH_KEYS:
        if key == config.SEQUENCE_LENGTH_KEY:
            path = config.DATA_RAW / "phase3" / "phase3_raw_scores_main.parquet"
            frames.append(pd.read_parquet(path))
            continue
        frames.append(phase3_scoring.run(sequence_length_key=key, tag=key))
    df = pd.concat(frames, ignore_index=True)
    out = RAW / "phase4_sequence_length_sensitivity.parquet"
    df.to_parquet(out, index=False)
    return df


# --------------------------------------------------------------------------
# QC4 - dependence on imputed tracks
# --------------------------------------------------------------------------
def qc4_observed_vs_imputed(scores: pd.DataFrame) -> pd.DataFrame:
    sel = scores[
        (~scores.is_activity_scorer)
        & (scores.modality_class == "primary")
    ].copy()
    sel["is_ghsr"] = sel.gene_name == config.GHSR_SYMBOL
    rows = []
    for (vk, rsid, cond, scorer), g in sel.groupby(
        ["variant_key", "rsid", "condition_label", "scorer_key"],
        observed=True,
    ):
        for subset_name, sub in (
            ("all_tracks", g),
            ("observed_only", g[g.imputation_mode == "observed"]),
            (
                "imputed_exploratory_only",
                g[g.imputation_mode == "imputed_exploratory"],
            ),
        ):
            if not len(sub):
                continue
            rows.append(
                {
                    "variant_key": vk,
                    "rsid": rsid,
                    "condition_label": cond,
                    "scorer_key": scorer,
                    "track_subset": subset_name,
                    "n_tracks": int(len(sub)),
                    "mean_raw_score": float(sub.raw_score.mean()),
                    "max_abs_raw_score": float(sub.raw_score.abs().max()),
                }
            )
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------
# QC5 - aggregated score vs raw track consistency
# --------------------------------------------------------------------------
def qc5_track_consistency(scores: pd.DataFrame) -> pd.DataFrame:
    """Compare Phase III centre-mask DNASE/ATAC scores with Phase III-b tracks."""
    idx_path = config.DATA_PROCESSED / "phase3b_track_index.csv"
    if not idx_path.exists():
        return pd.DataFrame()
    raw3b = config.DATA_RAW / "phase3b"
    rows = []
    for spec in config.VARIANTS:
        for curie in ("UBERON:0000007", "UBERON:0000945"):
            for ot, scorer_key in (("DNASE", "DNASE"), ("ATAC", "ATAC")):
                npz = (
                    raw3b
                    / f"tracks_{spec.key}_{curie.replace(':', '-')}_{ot}.npz"
                )
                if not npz.exists():
                    continue
                d = np.load(npz)
                ref, alt = d["ref"], d["alt"]
                res = int(d["resolution"])
                start = int(d["interval_start"])
                # Centre-mask scorers for DNASE/ATAC use a 501 bp window.
                half = 501 // 2
                lo = (spec.position - 1 - half - start) // res
                hi = (spec.position - 1 + half - start) // res + 1
                lo, hi = max(lo, 0), min(hi, ref.shape[0])
                ref_sum = ref[lo:hi].sum(axis=0)
                alt_sum = alt[lo:hi].sum(axis=0)
                with np.errstate(divide="ignore", invalid="ignore"):
                    lfc = np.log2(alt_sum) - np.log2(ref_sum)
                track_derived = float(np.nanmean(lfc))
                agg = scores[
                    (scores.variant_key == spec.key)
                    & (scores.condition_curie == curie)
                    & (scores.scorer_key == scorer_key)
                ]
                if not len(agg):
                    continue
                api_score = float(agg.raw_score.mean())
                rows.append(
                    {
                        "variant_key": spec.key,
                        "condition_curie": curie,
                        "condition_label": config.CONDITION_BY_CURIE[
                            curie
                        ].label,
                        "output_type": ot,
                        "api_aggregated_mean_score": api_score,
                        "recomputed_from_tracks_mean_log2fc": track_derived,
                        "absolute_difference": abs(api_score - track_derived),
                        "same_sign_or_both_negligible": bool(
                            np.sign(api_score) == np.sign(track_derived)
                            or max(abs(api_score), abs(track_derived)) < 1e-3
                        ),
                    }
                )
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------
def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--ism",
        action="store_true",
        help="run the ISM reference distribution (slowest QC step)",
    )
    ap.add_argument(
        "--seqlen",
        action="store_true",
        help="run the sequence-context sensitivity analysis",
    )
    args = ap.parse_args()

    print("Phase IV: robustness and quality control")
    scores = pd.read_parquet(
        config.DATA_RAW / "phase3" / "phase3_raw_scores_main.parquet"
    )

    qc1 = qc1_ref_alt_integrity(scores)
    (config.DATA_PROCESSED / "phase4_qc1_ref_alt_integrity.json").write_text(
        json.dumps(qc1, indent=2, default=str)
    )
    print("QC1 REF/ALT integrity:")
    print(
        json.dumps(
            {k: v for k, v in qc1.items() if k != "per_variant"}, indent=2
        )
    )

    qc4 = qc4_observed_vs_imputed(scores)
    qc4.to_csv(
        config.TABLES / "table_S6_observed_vs_imputed_tracks.csv", index=False
    )
    print(f"QC4 observed-vs-imputed: {len(qc4)} summary rows")

    qc5 = qc5_track_consistency(scores)
    if len(qc5):
        qc5.to_csv(
            config.TABLES / "table_S7_track_consistency.csv", index=False
        )
        print(
            "QC5 aggregated-vs-track consistency: "
            f"{int(qc5.same_sign_or_both_negligible.sum())}/{len(qc5)} consistent"
        )

    ism_pct = pd.DataFrame()
    if args.ism:
        ism = qc2_ism_reference_distribution()
        ism_pct = qc2_percentiles(ism)
        ism_pct.to_csv(
            config.TABLES / "table_S5_ism_reference_distribution.csv",
            index=False,
        )
        print("\nQC2 ISM reference distribution (percentile of |effect|):")
        print(
            ism_pct[
                [
                    "condition_label",
                    "scorer_key",
                    "rsid",
                    "observed_score",
                    "ism_max_abs",
                    "percentile_of_abs_effect_in_amplicon",
                ]
            ].to_string(index=False)
        )

    if args.seqlen:
        seq_df = qc3_sequence_length_sensitivity()
        ghsr = seq_df[
            (seq_df.gene_name == config.GHSR_SYMBOL)
            & (seq_df.scorer_key == "RNA_SEQ")
        ]
        pivot = ghsr.pivot_table(
            index=["rsid", "condition_label"],
            columns="sequence_length_key",
            values="raw_score",
            aggfunc="mean",
        )
        pivot.to_csv(
            config.TABLES / "table_S8_sequence_length_sensitivity.csv"
        )
        print("\nQC3 GHSR RNA-seq LFC by sequence context length:")
        print(pivot.to_string())

    provenance.write(
        "phase4",
        {
            "qc1_summary": {
                k: v for k, v in qc1.items() if k != "per_variant"
            },
            "ism_run": bool(args.ism),
            "ism_conditions": list(ISM_CONDITIONS) if args.ism else [],
            "ism_scorers": list(ISM_SCORER_KEYS) if args.ism else [],
            "seqlen_run": bool(args.seqlen),
            "n_qc4_rows": int(len(qc4)),
            "n_qc5_rows": int(len(qc5)),
        },
    )


if __name__ == "__main__":
    main()
