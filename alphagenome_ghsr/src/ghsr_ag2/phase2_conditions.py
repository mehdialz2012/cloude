"""Phase II - condition (cell type / tissue) selection and data-quality audit.

AlphaGenome 2 is a condition-conditioned model: predictions are requested for
one reference condition at a time, identified by a composite ontology CURIE
(CLO / UBERON / CL / EFO / NTR). `DnaClient.supported_conditions()` returns the
full set of reference conditions together with, for each output type, whether
the corresponding tracks are:

  observed            - real experimental data for this condition exists in the
                        model's training corpus
  imputed_supported   - no data for this condition/output pair, but the
                        condition has related observed data supporting
                        imputation
  imputed_exploratory - imputed with no closely related observed data; the
                        weakest evidence class

This module:
  1. Downloads and archives the complete supported-condition table.
  2. Confirms every pre-specified condition in `config.CONDITIONS` exists.
  3. Records the imputation mode of each output type for each selected
     condition, so that no result can be reported without its evidence class.
  4. Documents pre-specified biological priorities that the model does not
     support at all (e.g. osteoclasts).
  5. Measures predicted GHSR activity per condition (`RNA_SEQ_ACTIVE`, an
     absolute-activity readout, not an effect size) to verify that GHSR
     regulation is interpretable in each selected condition.

Run:  python -m ghsr_ag2.phase2_conditions
"""

from __future__ import annotations

import json

import numpy as np
import pandas as pd
from alphagenome.data import genome
from alphagenome.models.v2 import dna_client, variant_scorers

from . import config, provenance

OUTPUT_TYPE_COLUMNS = (
    "rna_seq",
    "cage",
    "dnase",
    "atac",
    "chip_tf",
    "chip_histone",
    "sc_atac",
    "sc_gene_expression",
    "methylation",
    "splice_sites",
    "splice_site_usage",
    "splice_junctions",
    "contact_maps",
)

# Pre-specified biological priorities, with the regex used to search the
# supported-condition catalogue for each.
PRIORITY_SEARCHES = {
    "hippocampal neurons": r"hippocamp",
    "hypothalamic neurons": r"hypothalam",
    "other brain neuronal populations": r"neuron",
    "pituitary cells": r"pituitar",
    "osteoblasts": r"osteoblast",
    "osteocytes": r"osteocyte",
    "osteoclast-lineage cells": r"osteoclast",
}


def get_client() -> dna_client.DnaClient:
    return dna_client.create(config.api_key(), timeout=120)


def fetch_supported_conditions(
    client: dna_client.DnaClient,
) -> pd.DataFrame:
    df = client.supported_conditions()
    out = config.DATA_RAW / "alphagenome2_supported_conditions_human.csv"
    df.to_csv(out, index=False)
    print(f"[phase2] archived {len(df)} supported conditions -> {out}")
    return df


def audit_priorities(supported: pd.DataFrame) -> pd.DataFrame:
    """Record, for each pre-specified priority, what the model actually has."""
    rows = []
    for priority, pattern in PRIORITY_SEARCHES.items():
        hits = supported[
            supported.biosample_name.str.contains(
                pattern, case=False, na=False, regex=True
            )
        ]
        selected = [
            c.label
            for c in config.CONDITIONS
            if c.curie in set(hits.condition)
        ]
        rows.append(
            {
                "biological_priority": priority,
                "n_conditions_in_model": len(hits),
                "available_in_model": len(hits) > 0,
                "n_with_observed_rna_seq": int(
                    (hits.rna_seq == "observed").sum()
                ),
                "conditions_selected_for_this_study": "; ".join(selected)
                or "none",
                "resolution": (
                    "no condition of this type exists in the AlphaGenome 2 "
                    "human condition set; a surrogate was used and is "
                    "flagged in the condition table"
                    if len(hits) == 0
                    else "available and represented in the analysis panel"
                ),
            }
        )
    return pd.DataFrame(rows)


def build_condition_table(supported: pd.DataFrame) -> pd.DataFrame:
    """Join the pre-specified panel against the model's own metadata."""
    idx = supported.set_index("condition")
    missing = [c.curie for c in config.CONDITIONS if c.curie not in idx.index]
    if missing:
        raise RuntimeError(
            "Pre-specified conditions absent from the model's supported "
            f"condition set: {missing}. Do not invent condition availability; "
            "update config.CONDITIONS instead."
        )

    rows = []
    for spec in config.CONDITIONS:
        rec = idx.loc[spec.curie]
        modes = {col: rec.get(col) for col in OUTPUT_TYPE_COLUMNS}
        n_obs = sum(1 for v in modes.values() if v == "observed")
        n_sup = sum(1 for v in modes.values() if v == "imputed_supported")
        n_exp = sum(1 for v in modes.values() if v == "imputed_exploratory")
        rows.append(
            {
                "condition_curie": spec.curie,
                "label": spec.label,
                "panel": spec.panel,
                "model_biosample_name": rec["biosample_name"],
                "model_biosample_type": rec["biosample_type"],
                "target_priority": spec.target_priority,
                "surrogate_for": spec.surrogate_for or "",
                **{f"mode_{k}": v for k, v in modes.items()},
                "n_output_types_observed": n_obs,
                "n_output_types_imputed_supported": n_sup,
                "n_output_types_imputed_exploratory": n_exp,
                "evidence_tier": (
                    "A: observed RNA-seq"
                    if modes["rna_seq"] == "observed"
                    else "B: imputed-supported RNA-seq"
                    if modes["rna_seq"] == "imputed_supported"
                    else "C: imputed-exploratory RNA-seq"
                ),
                "rationale": spec.rationale,
            }
        )
    df = pd.DataFrame(rows)
    df["panel"] = pd.Categorical(
        df["panel"], categories=list(config.PANEL_ORDER), ordered=True
    )
    return df.sort_values(["panel", "evidence_tier", "label"]).reset_index(
        drop=True
    )


def measure_ghsr_activity(
    client: dna_client.DnaClient,
) -> pd.DataFrame:
    """Predicted absolute GHSR activity per condition (QC, not an effect size).

    Uses the `RNA_SEQ_ACTIVE` / `CAGE_ACTIVE` gene-mask active-allele scorers,
    which return max(ALT, REF) aggregated signal over the gene's exons. Because
    the variants studied here have negligible predicted effects, this value is
    effectively the model's predicted reference expression level and is used
    only to confirm that GHSR is expressed (and therefore that a regulatory
    effect on it would be detectable) in each condition.
    """
    spec = config.PRIMARY_VARIANTS[0]
    variant = genome.Variant(
        chromosome=config.CHROM,
        position=spec.position,
        reference_bases=spec.ref,
        alternate_bases=spec.alt,
        name=spec.rsid,
    )
    interval = variant.reference_interval.resize(
        dna_client.SUPPORTED_SEQUENCE_LENGTHS[config.SEQUENCE_LENGTH_KEY]
    )
    # RNA_SEQ_ACTIVE is a gene-mask scorer: rows are genes, so GHSR activity
    # can be read out directly. CAGE_ACTIVE / DNASE_ACTIVE are centre-mask
    # scorers with no gene axis; for those the value describes the 501 bp
    # window centred on the variant itself, not the GHSR promoter.
    gene_scorer_keys = ("RNA_SEQ_ACTIVE",)
    local_scorer_keys = ("CAGE_ACTIVE", "DNASE_ACTIVE")
    scorers = [
        variant_scorers.RECOMMENDED_VARIANT_SCORERS[k]
        for k in (*gene_scorer_keys, *local_scorer_keys)
    ]

    def _summarise(adata, row_index: int, prefix: str, rec: dict) -> None:
        x = np.asarray(adata.X)[row_index, :].astype(float)
        obs_mask = (adata.var["imputation_mode"] == "observed").to_numpy()
        rec[f"{prefix}_max"] = float(np.nanmax(x)) if x.size else np.nan
        rec[f"{prefix}_mean"] = float(np.nanmean(x)) if x.size else np.nan
        rec[f"{prefix}_n_tracks"] = int(x.size)
        rec[f"{prefix}_max_observed_tracks"] = (
            float(np.nanmax(x[obs_mask])) if obs_mask.any() else np.nan
        )
        rec[f"{prefix}_n_observed_tracks"] = int(obs_mask.sum())

    rows = []
    for cond in config.CONDITIONS:
        outs = client.score_variant(
            interval=interval,
            variant=variant,
            condition=cond.curie,
            variant_scorers=scorers,
        )
        rec: dict[str, object] = {
            "condition_curie": cond.curie,
            "label": cond.label,
            "panel": cond.panel,
        }
        for adata, key in zip(outs, (*gene_scorer_keys, *local_scorer_keys)):
            prefix_key = key.lower()
            if "gene_name" in adata.obs.columns:
                mask = (adata.obs["gene_name"] == config.GHSR_SYMBOL).to_numpy()
                if not mask.any():
                    rec[f"ghsr_{prefix_key}_max"] = np.nan
                    rec[f"ghsr_{prefix_key}_n_tracks"] = 0
                    continue
                _summarise(
                    adata, int(np.flatnonzero(mask)[0]), f"ghsr_{prefix_key}", rec
                )
                # Within-condition reference distribution: the same readout
                # for every other gene in the same 1 Mb window. This
                # self-normalises the arbitrary model units so that GHSR
                # activity can be judged without inventing an absolute cut-off.
                x_all = np.nanmax(np.asarray(adata.X).astype(float), axis=1)
                ghsr_val = x_all[int(np.flatnonzero(mask)[0])]
                others = x_all[~mask]
                pc = (adata.obs["gene_type"] == "protein_coding").to_numpy()
                rec[f"ghsr_{prefix_key}_n_genes_in_window"] = int(x_all.size)
                rec[f"ghsr_{prefix_key}_rank_in_window"] = int(
                    (x_all > ghsr_val).sum() + 1
                )
                rec[f"ghsr_{prefix_key}_rank_among_protein_coding"] = int(
                    (x_all[pc] > ghsr_val).sum() + 1
                )
                rec[f"ghsr_{prefix_key}_n_protein_coding_in_window"] = int(
                    pc.sum()
                )
                rec[f"ghsr_{prefix_key}_window_median"] = float(
                    np.nanmedian(others)
                )
                rec[f"ghsr_{prefix_key}_window_max"] = float(np.nanmax(others))
                rec[f"ghsr_{prefix_key}_over_window_median"] = (
                    float(ghsr_val / np.nanmedian(others))
                    if np.nanmedian(others) > 0
                    else np.nan
                )
                rec[f"ghsr_{prefix_key}_over_window_max"] = (
                    float(ghsr_val / np.nanmax(others))
                    if np.nanmax(others) > 0
                    else np.nan
                )
            else:
                _summarise(adata, 0, f"variantlocus_{prefix_key}", rec)
        rows.append(rec)
        print(
            f"[phase2] {cond.label:<28} "
            f"GHSR RNA-seq activity = "
            f"{rec.get('ghsr_rna_seq_active_max', float('nan')):.3f}"
        )
    return pd.DataFrame(rows)


def main() -> None:
    print("Phase II: condition selection and data-quality audit")
    client = get_client()
    supported = fetch_supported_conditions(client)

    priorities = audit_priorities(supported)
    conditions = build_condition_table(supported)
    activity = measure_ghsr_activity(client)

    merged = conditions.merge(
        activity.drop(columns=["label", "panel"]),
        on="condition_curie",
        how="left",
    )
    # GHSR interpretability banding.
    #
    # The gene-mask active RNA-seq score is in the model's own arbitrary
    # units, so there is no meaningful absolute cut-off and none is invented
    # here. Instead GHSR activity is expressed relative to the highest value
    # observed across this condition panel, and banded descriptively. This is
    # a transparent descriptive banding of predicted expression, NOT a
    # statistical test and NOT a claim about measured expression.
    col = "ghsr_rna_seq_active_max"
    panel_max = float(merged[col].max())
    merged["ghsr_activity_relative_to_panel_max"] = merged[col] / panel_max

    def _band(rel: float) -> str:
        if rel >= 0.25:
            return "high (>=25% of panel maximum)"
        if rel >= 0.05:
            return "intermediate (5-25% of panel maximum)"
        return "low / near-zero (<5% of panel maximum)"

    merged["ghsr_predicted_expression_band"] = merged[
        "ghsr_activity_relative_to_panel_max"
    ].map(_band)
    merged["ghsr_regulation_interpretable"] = (
        merged["ghsr_activity_relative_to_panel_max"] >= 0.05
    )
    merged["ghsr_interpretability_basis"] = (
        "predicted gene-mask active RNA-seq signal expressed relative to the "
        f"panel maximum ({panel_max:.4f}, pituitary gland); descriptive "
        "banding only, not a statistical test"
    )

    p_cond = config.TABLES / "table_S3_condition_selection_and_qc.csv"
    merged.to_csv(p_cond, index=False)
    p_prio = config.TABLES / "table_S4_priority_availability_audit.csv"
    priorities.to_csv(p_prio, index=False)
    supported_sel = supported[
        supported.condition.isin([c.curie for c in config.CONDITIONS])
    ]
    supported_sel.to_csv(
        config.DATA_PROCESSED / "phase2_selected_conditions_raw_modes.csv",
        index=False,
    )

    print()
    print(
        merged[
            [
                "label",
                "panel",
                "evidence_tier",
                "mode_rna_seq",
                "mode_dnase",
                "mode_chip_tf",
                "ghsr_rna_seq_active_max",
                "ghsr_rna_seq_active_rank_among_protein_coding",
                "ghsr_activity_relative_to_panel_max",
                "ghsr_predicted_expression_band",
            ]
        ].to_string(index=False)
    )
    print()
    print(priorities.to_string(index=False))
    print(f"\nWrote {p_cond}\nWrote {p_prio}")

    provenance.write(
        "phase2",
        {
            "n_supported_conditions_human": int(len(supported)),
            "n_conditions_selected": int(len(config.CONDITIONS)),
            "unavailable_priorities": priorities.loc[
                ~priorities.available_in_model, "biological_priority"
            ].tolist(),
            "ghsr_interpretability_basis": "relative to panel maximum",
            "scorers_used": ["RNA_SEQ_ACTIVE", "CAGE_ACTIVE", "DNASE_ACTIVE"],
            "sequence_length": config.SEQUENCE_LENGTH_KEY,
        },
    )
    (config.DATA_PROCESSED / "phase2_summary.json").write_text(
        json.dumps(
            {
                "n_supported_conditions": int(len(supported)),
                "selected": [c.curie for c in config.CONDITIONS],
                "unavailable_priorities": priorities.loc[
                    ~priorities.available_in_model, "biological_priority"
                ].tolist(),
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
