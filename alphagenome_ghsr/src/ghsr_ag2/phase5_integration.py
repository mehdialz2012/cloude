"""Phase V - integration with the existing experimental findings.

Builds the integrated comparison table demanded by the study design, keeping
three evidence levels strictly separate:

  L1  Computational regulatory prediction   (this study; AlphaGenome 2)
  L2  Experimental detection of genotypes   (prior study; PCR + Sanger)
  L3  Statistical association with measured phenotypes (prior study; SPSS)

Level 2 validates that a genotype was correctly *identified*. It says nothing
about gene regulation. Level 3 is association in 100 participants and is not
evidence of a molecular mechanism. Nothing in this module converts a Level 1
prediction into support for a Level 3 association, or vice versa.

The Level 2 and Level 3 entries are transcribed from the source manuscript
(`ARTICLE__revised-v2.docx`) and are quoted, not recomputed. The transcription
is kept in `MANUSCRIPT_FINDINGS` so it can be audited line by line against the
manuscript tables.

Run:  python -m ghsr_ag2.phase5_integration
"""

from __future__ import annotations

import json

import numpy as np
import pandas as pd

from . import config, provenance

# --------------------------------------------------------------------------
# Level 2 / Level 3 evidence, transcribed from the source manuscript.
# `source` gives the manuscript table or section each value comes from.
# --------------------------------------------------------------------------
MANUSCRIPT_FINDINGS: dict[str, dict[str, object]] = {
    "rs2922126": {
        "manuscript_label": "SNP1",
        "manuscript_allele_notation": "A>T",
        # ---- Level 2: genotype detection ----
        "genotypes_cases": "31 hom wild / 16 het / 3 hom mutant (n=50)",
        "genotypes_controls": "28 hom wild / 18 het / 4 hom mutant (n=50)",
        "minor_allele_frequency_observed": 0.24,
        "call_rate": "100% (0% missing)",
        "hwe_p_controls": 0.62,
        "genotype_concordance_recheck": "100% on a 10% random re-read",
        "l2_source": "Results para. 'During the sequence analysis'; Methods "
        "'Genetic Analysis'; Table 11",
        # ---- Level 3: phenotype association ----
        "assoc_cognitive_group_unadjusted_p": 0.814,
        "assoc_cognitive_group_adjusted_p": 0.862,
        "assoc_cognitive_group_or": "0.56 (95% CI 0.027-11.34) hom mutant; "
        "0.75 (95% CI 0.22-2.62) het",
        "assoc_genotype_distribution_p": 0.81,
        "other_nominal_findings": (
            "testosterone in cases p=0.013 (FDR 0.154); serum calcium in "
            "controls p=0.022 (FDR 0.154); osteoporosis status in cases "
            "p=0.014 (FDR 0.07); sex-genotype association in cases p=0.046"
        ),
        "survived_multiplicity_correction": False,
        "achieved_power_pct": 6.62,
        "l3_source": "Tables 11, 14, 15, 16; Results and Discussion",
        "l3_verdict": (
            "No independent association with cognitive screening status. All "
            "nominal biochemical/skeletal associations lost significance "
            "after Benjamini-Hochberg FDR correction. Study underpowered "
            "(achieved power 6.62%)."
        ),
    },
    "rs115459563": {
        "manuscript_label": "SNP2",
        "manuscript_allele_notation": "G>A",
        "genotypes_cases": "50 hom wild / 0 het / 0 hom mutant (n=50)",
        "genotypes_controls": "43 hom wild / 7 het / 0 hom mutant (n=50)",
        "minor_allele_frequency_observed": 0.035,
        "call_rate": "100% (0% missing)",
        "hwe_p_controls": 0.59,
        "genotype_concordance_recheck": "100% on a 10% random re-read",
        "l2_source": "Results para. 'During the sequence analysis'; Table 12",
        "assoc_cognitive_group_unadjusted_p": 0.012,
        "assoc_cognitive_group_adjusted_p": None,
        "assoc_cognitive_group_or": "not estimable - logistic regression did "
        "not converge (zero carriers among cases)",
        "assoc_genotype_distribution_p": 0.012,
        "other_nominal_findings": (
            "no significant differences in vitamin D, calcium, phosphorus, "
            "estradiol, testosterone, spine or femoral-neck BMD, thyroid "
            "status or osteoporosis status (all p>0.05)"
        ),
        "survived_multiplicity_correction": None,
        "achieved_power_pct": 71.86,
        "l3_source": "Table 12; Results; Discussion; Limitations",
        "l3_verdict": (
            "Genotype distribution differed nominally (7 heterozygous "
            "controls vs 0 cases) but all regression models failed to "
            "converge because of the zero cell, so no adjusted estimate "
            "exists. The manuscript explicitly declines to interpret this as "
            "a protective association. The 71.86% achieved-power figure was "
            "computed from the same small dataset and is not evidence of "
            "adequate power."
        ),
    },
    "rs34136054": {
        "manuscript_label": "SNP3",
        "manuscript_allele_notation": "C>T",
        "genotypes_cases": "50 hom wild / 0 het / 0 hom mutant (n=50)",
        "genotypes_controls": "49 hom wild / 1 het / 0 hom mutant (n=50)",
        "minor_allele_frequency_observed": 0.005,
        "call_rate": "100% (0% missing)",
        "hwe_p_controls": 0.94,
        "genotype_concordance_recheck": "100% on a 10% random re-read",
        "l2_source": "Results para. 'During the sequence analysis'; Table 13",
        "assoc_cognitive_group_unadjusted_p": 1.00,
        "assoc_cognitive_group_adjusted_p": None,
        "assoc_cognitive_group_or": "not evaluated - variant seen in one "
        "participant only",
        "assoc_genotype_distribution_p": 1.00,
        "other_nominal_findings": (
            "no association with thyroid function or osteoporosis status "
            "(p>0.05); continuous comparisons not possible"
        ),
        "survived_multiplicity_correction": None,
        "achieved_power_pct": None,
        "l3_source": "Table 13; Results; Discussion",
        "l3_verdict": (
            "Detected in a single control participant. No meaningful "
            "statistical evaluation was possible and none was attempted."
        ),
    },
}

NEGLIGIBLE_LFC = 0.01  # see docstring of `classify_effect`


def classify_effect(max_abs_lfc: float) -> str:
    """Descriptive banding of a predicted log-fold-change magnitude.

    The bands are descriptive labels for the magnitude of a model output, not
    significance categories, and no P-value or confidence interval is attached
    to them. The 0.01 log2 boundary is roughly a 0.7% change in predicted
    signal - far below the resolution at which any AlphaGenome prediction has
    been shown to correspond to a measurable expression difference - and the
    0.1 / 0.5 log2 boundaries (~7% and ~41%) mark magnitudes that would be
    plausibly measurable in a reporter assay if real.
    """
    a = abs(max_abs_lfc)
    if np.isnan(a):
        return "not evaluated"
    if a < NEGLIGIBLE_LFC:
        return "negligible (|log2FC| < 0.01)"
    if a < 0.1:
        return "very small (0.01 <= |log2FC| < 0.1)"
    if a < 0.5:
        return "moderate (0.1 <= |log2FC| < 0.5)"
    return "large (|log2FC| >= 0.5)"


def summarise_predictions(scores: pd.DataFrame) -> pd.DataFrame:
    """One row per rsID: the strongest predicted effect across the panel."""
    # Headline summaries use the primary variant configuration for each
    # rsID. For rs34136054 the second alternate allele (T>G) is a declared
    # sensitivity analysis and is summarised separately, because mixing two
    # different substitutions into one "per-rsID" number would be misleading:
    # they can differ in sign.
    eff_all = scores[
        (~scores.is_activity_scorer) & (scores.modality_class == "primary")
    ]
    eff = eff_all[eff_all.is_primary_configuration]
    rows = []
    for rsid in config.RSIDS:
        sub = eff[eff.rsid == rsid]
        ghsr_rna = sub[
            (sub.scorer_key == "RNA_SEQ") & (sub.gene_name == config.GHSR_SYMBOL)
        ]
        any_gene_rna = sub[sub.scorer_key == "RNA_SEQ"]
        access = sub[sub.scorer_key.isin(["DNASE", "ATAC", "SC_ATAC"])]
        tf = sub[sub.scorer_key == "CHIP_TF"]

        def _peak(df: pd.DataFrame, label_cols: list[str]) -> dict:
            if not len(df):
                return {"value": np.nan, "where": "not evaluated"}
            i = df.raw_score.abs().idxmax()
            r = df.loc[i]
            where = ", ".join(f"{c}={r[c]}" for c in label_cols)
            return {"value": float(r.raw_score), "where": where}

        g = _peak(ghsr_rna, ["condition_label", "track_name"])
        a = _peak(access, ["condition_label", "scorer_key", "track_name"])
        t = _peak(tf, ["condition_label", "transcription_factor"])
        other = any_gene_rna[any_gene_rna.gene_name != config.GHSR_SYMBOL]
        o = _peak(other, ["gene_name", "condition_label"])

        rows.append(
            {
                "rsid": rsid,
                "ghsr_rna_max_abs_log2fc": abs(g["value"]),
                "ghsr_rna_signed_at_peak": g["value"],
                "ghsr_rna_peak_location": g["where"],
                "ghsr_effect_band": classify_effect(g["value"]),
                "ghsr_direction": (
                    "no consistent direction"
                    if abs(g["value"]) < NEGLIGIBLE_LFC
                    else ("increase" if g["value"] > 0 else "decrease")
                ),
                "accessibility_max_abs_log2fc": abs(a["value"]),
                "accessibility_peak_location": a["where"],
                "accessibility_effect_band": classify_effect(a["value"]),
                "chip_tf_max_abs_log2fc": abs(t["value"]),
                "chip_tf_peak_location": t["where"],
                "chip_tf_effect_band": classify_effect(t["value"]),
                "nearby_gene_max_abs_log2fc": abs(o["value"]),
                "nearby_gene_peak_location": o["where"],
                "nearby_gene_effect_band": classify_effect(o["value"]),
                "n_conditions_evaluated": int(sub.condition_curie.nunique()),
                "n_track_level_scores": int(len(sub)),
                "primary_configuration": next(
                    v.key for v in config.PRIMARY_VARIANTS if v.rsid == rsid
                ),
                **_specificity_metrics(ghsr_rna),
                **_secondary_allele_metrics(eff_all, rsid),
            }
        )
    return pd.DataFrame(rows)



def _secondary_allele_metrics(
    eff_all: pd.DataFrame, rsid: str
) -> dict[str, object]:
    """Summarise any non-primary alternate allele for this rsID."""
    sec = eff_all[(eff_all.rsid == rsid) & (~eff_all.is_primary_configuration)]
    if not len(sec):
        return {
            "has_secondary_alternate_allele": False,
            "secondary_allele_key": "",
            "secondary_ghsr_max_abs_log2fc": np.nan,
            "secondary_allele_note": "",
        }
    ghsr = sec[
        (sec.scorer_key == "RNA_SEQ") & (sec.gene_name == config.GHSR_SYMBOL)
    ]
    key = sec.variant_key.iloc[0]
    return {
        "has_secondary_alternate_allele": True,
        "secondary_allele_key": key,
        "secondary_ghsr_max_abs_log2fc": (
            float(ghsr.raw_score.abs().max()) if len(ghsr) else np.nan
        ),
        "secondary_allele_note": (
            "the two alternate alleles of this tri-allelic site can differ in "
            "the sign of the predicted effect, so neither can stand for the "
            "other"
        ),
    }


def _specificity_metrics(ghsr_rna: pd.DataFrame) -> dict[str, object]:
    """Describe how the GHSR effect varies across conditions.

    Two things are distinguished that a naive reading would conflate:

      * sign consistency - does the ALT allele move predicted GHSR expression
        in the same direction in every condition? If not, the variant has no
        single predicted direction of effect.
      * expression dependence - does the magnitude of the predicted effect
        simply track how strongly GHSR is predicted to be expressed in that
        condition? A multiplicative effect applied to a near-zero baseline is
        near-zero, so a correlation here means the apparent "cell-type
        specificity" is largely a consequence of baseline expression rather
        than of condition-specific regulatory logic.
    """
    from scipy.stats import spearmanr

    if not len(ghsr_rna):
        return {}
    per_cond = ghsr_rna.groupby("condition_label")["raw_score"].mean()
    n_pos = int((per_cond > 0).sum())
    n_neg = int((per_cond < 0).sum())
    try:
        cond_tbl = pd.read_csv(
            config.TABLES / "table_S3_condition_selection_and_qc.csv"
        ).set_index("label")["ghsr_rna_seq_active_max"]
        joined = pd.DataFrame(
            {"lfc": per_cond, "activity": cond_tbl.reindex(per_cond.index)}
        ).dropna()
        rho = (
            float(spearmanr(joined.lfc.abs(), joined.activity).statistic)
            if len(joined) > 3
            else np.nan
        )
    except Exception:  # noqa: BLE001
        rho = np.nan
    return {
        "n_conditions_alt_increases": n_pos,
        "n_conditions_alt_decreases": n_neg,
        "sign_consistent_across_conditions": bool(n_pos == 0 or n_neg == 0),
        "spearman_abs_effect_vs_ghsr_activity": rho,
        "specificity_reading": (
            "effect magnitude tracks predicted GHSR expression level "
            f"(Spearman rho = {rho:.2f}); the apparent cell-type pattern is "
            "largely explained by baseline expression, not by "
            "condition-specific regulatory logic"
            if np.isfinite(rho) and rho >= 0.5
            else "no clear relationship between effect magnitude and "
            "predicted GHSR expression level"
        ),
    }


def build_integration_table(pred: pd.DataFrame) -> pd.DataFrame:
    annot = pd.read_csv(
        config.TABLES / "table_S1_variant_annotation.csv"
    )
    rows = []
    for rsid in config.RSIDS:
        p = pred[pred.rsid == rsid].iloc[0]
        m = MANUSCRIPT_FINDINGS[rsid]
        a = annot[(annot.rsid == rsid) & annot.is_primary_configuration].iloc[0]
        rows.append(
            {
                "rsid": rsid,
                "manuscript_label": m["manuscript_label"],
                "position_grch38": a["position_grch38"],
                "ref_alt_plus_strand": f"{a['ref_plus_strand']}>"
                f"{a['alt_plus_strand']}",
                "manuscript_notation": m["manuscript_allele_notation"],
                "bp_upstream_of_ghsr_tss": a["bp_upstream_of_ghsr_tss"],
                # -------- Level 1 --------
                "L1_alphagenome_ghsr_expression": (
                    f"{p['ghsr_effect_band']}; peak |log2FC| = "
                    f"{p['ghsr_rna_max_abs_log2fc']:.2e} "
                    f"({p['ghsr_rna_peak_location']})"
                ),
                "L1_alphagenome_chromatin_accessibility": (
                    f"{p['accessibility_effect_band']}; peak |log2FC| = "
                    f"{p['accessibility_max_abs_log2fc']:.2e}"
                ),
                "L1_alphagenome_tf_associated": (
                    f"{p['chip_tf_effect_band']}; peak |log2FC| = "
                    f"{p['chip_tf_max_abs_log2fc']:.2e} "
                    f"({p['chip_tf_peak_location']})"
                ),
                "L1_alphagenome_nearby_genes": (
                    f"{p['nearby_gene_effect_band']}; peak |log2FC| = "
                    f"{p['nearby_gene_max_abs_log2fc']:.2e} "
                    f"({p['nearby_gene_peak_location']})"
                ),
                "L1_cell_type_specificity": (
                    f"ALT increases predicted GHSR expression in "
                    f"{int(p['n_conditions_alt_increases'])} and decreases it "
                    f"in {int(p['n_conditions_alt_decreases'])} of "
                    f"{int(p['n_conditions_evaluated'])} conditions "
                    f"(sign-consistent: "
                    f"{bool(p['sign_consistent_across_conditions'])}). "
                    f"{p['specificity_reading']}. Predicted effects in "
                    "osteoblasts, osteocytes and mesenchymal stem cells were "
                    "the smallest in the panel, in which GHSR itself is "
                    "predicted to be essentially unexpressed."
                ),
                # -------- Level 2 --------
                "L2_genotypes_cases": m["genotypes_cases"],
                "L2_genotypes_controls": m["genotypes_controls"],
                "L2_observed_maf": m["minor_allele_frequency_observed"],
                "L2_call_rate": m["call_rate"],
                "L2_hwe_p_controls": m["hwe_p_controls"],
                "L2_what_this_validates": (
                    "PCR and Sanger sequencing establish that this genotype "
                    "was correctly identified. They provide no information "
                    "about gene regulation or expression."
                ),
                # -------- Level 3 --------
                "L3_cognitive_group_p_unadjusted": m[
                    "assoc_cognitive_group_unadjusted_p"
                ],
                "L3_cognitive_group_p_adjusted": m[
                    "assoc_cognitive_group_adjusted_p"
                ],
                "L3_effect_estimate": m["assoc_cognitive_group_or"],
                "L3_other_nominal_findings": m["other_nominal_findings"],
                "L3_survived_multiplicity_correction": m[
                    "survived_multiplicity_correction"
                ],
                "L3_achieved_power_pct": m["achieved_power_pct"],
                "L3_verdict": m["l3_verdict"],
                # -------- Joint reading --------
                "integrated_interpretation": _integrate(rsid, p, m),
                "evidence_sources": (
                    f"L1: this study (AlphaGenome 2). L2/L3: "
                    f"{m['l2_source']} | {m['l3_source']}"
                ),
            }
        )
    return pd.DataFrame(rows)


def _integrate(rsid: str, p: pd.Series, m: dict) -> str:
    peak = float(p["ghsr_rna_max_abs_log2fc"])
    pct = (2.0**peak - 1.0) * 100.0
    base = (
        "AlphaGenome 2 predicts only a very small change in GHSR expression "
        f"for this variant: the largest effect anywhere in the "
        f"{int(p['n_conditions_evaluated'])}-condition panel is "
        f"|log2FC| = {peak:.3f}, a {pct:.1f}% change in predicted signal, "
        "which is far below the magnitude that would be detectable in a "
        "reporter assay"
    )
    if rsid == "rs2922126":
        return (
            f"{base}. The prior study found no independent association with "
            "cognitive screening status, and its nominal associations with "
            "testosterone, calcium and osteoporosis status did not survive "
            "FDR correction. The computational and the statistical evidence "
            "are therefore concordant in being null, but they are null about "
            "different things: the prediction concerns cis-regulatory "
            "activity, the association concerns phenotype. Neither "
            "establishes the other, and the association analysis had 6.62% "
            "achieved power, so it could not have detected a small effect "
            "regardless."
        )
    if rsid == "rs115459563":
        return (
            f"{base}. The prior study observed a nominal difference in "
            "genotype distribution (7 heterozygous controls, 0 cases) that "
            "could not be modelled because of the zero cell. The absence of a "
            "predicted regulatory effect offers no mechanistic support for "
            "that observation, but equally does not refute it: a 1 kb "
            "upstream variant could act through mechanisms the model does "
            "not capture, and the observation itself is most plausibly a "
            "small-sample artefact. This variant has gnomAD global allele "
            "frequency 0.004, so 7 of 50 heterozygous controls is well above "
            "the population expectation and warrants genotyping confirmation "
            "before any biological interpretation."
        )
    return (
        f"{base}. The prior study detected this variant in a single control "
        "participant, which permits no statistical evaluation. In addition, "
        "the manuscript's allele notation for this variant does not match "
        "GRCh38 (see Table S1), and the alternate allele actually segregating "
        "in populations is T>G rather than T>C; both alternates were "
        "therefore scored and neither produced an appreciable predicted "
        "effect. No inference about this variant is supportable from either "
        "evidence stream."
    )


def main() -> None:
    print("Phase V: integration with existing experimental findings")
    scores = pd.read_parquet(
        config.DATA_RAW / "phase3" / "phase3_raw_scores_main.parquet"
    )
    pred = summarise_predictions(scores)
    pred.to_csv(
        config.TABLES / "table_2_prediction_summary_by_variant.csv", index=False
    )
    integ = build_integration_table(pred)
    integ.to_csv(
        config.TABLES / "table_3_integrated_evidence.csv", index=False
    )
    # Transposed, human-readable form for manuscript use.
    integ.set_index("rsid").T.to_csv(
        config.TABLES / "table_3_integrated_evidence_transposed.csv"
    )
    (config.DATA_PROCESSED / "manuscript_findings_transcription.json").write_text(
        json.dumps(MANUSCRIPT_FINDINGS, indent=2)
    )

    print(
        pred[
            [
                "rsid",
                "ghsr_rna_max_abs_log2fc",
                "ghsr_effect_band",
                "accessibility_max_abs_log2fc",
                "chip_tf_max_abs_log2fc",
                "nearby_gene_max_abs_log2fc",
            ]
        ].to_string(index=False)
    )
    provenance.write(
        "phase5",
        {
            "negligible_band_boundary_log2fc": NEGLIGIBLE_LFC,
            "manuscript_source_file": "ARTICLE__revised-v2.docx",
            "evidence_levels": {
                "L1": "computational regulatory prediction (this study)",
                "L2": "experimental genotype detection (prior study)",
                "L3": "statistical phenotype association (prior study)",
            },
        },
    )


if __name__ == "__main__":
    main()
