#!/usr/bin/env python3
"""Generate the executable analysis notebook.

The notebook is generated from this script so that it stays in step with the
pipeline modules rather than drifting from them. Run:

    python scripts/build_notebook.py          # write the notebook
    python scripts/build_notebook.py --execute  # write and execute it
"""

from __future__ import annotations

import argparse
import pathlib
import sys

import nbformat as nbf

ROOT = pathlib.Path(__file__).resolve().parents[1]
OUT = ROOT / "notebooks" / "AlphaGenome2_GHSR_analysis.ipynb"

MD = nbf.v4.new_markdown_cell
CODE = nbf.v4.new_code_cell

CELLS = [
    MD(
        """# AlphaGenome 2 functional analysis of three *GHSR* upstream variants

A **retrospective, complementary in silico** investigation of the three *GHSR*
polymorphisms genotyped in Saeidi et al., *"Ghrelin Receptor Gene
Polymorphisms, Cognitive Screening Status, and Osteoporosis-Related
Parameters: An Exploratory Case-Control Study"*.

The computational analysis was run **after** the genotyping was complete. It
did not guide the experimental design, the variant selection or the statistics
of that study.

**Hypothesis under test:** that variants in regulatory regions of *GHSR*
influence its transcription in a cell-type-specific manner, which would
supply a molecular mechanism relevant to cognition and bone metabolism.

**What this notebook can and cannot show.** AlphaGenome returns model
predictions, not measurements. No P-value, confidence interval or significance
threshold is attached to any score, because none exists. A near-zero predicted
effect is evidence about the model's expectation, not proof of biological
inertness; a non-zero predicted effect is a hypothesis, not a finding.

---

### Requirements

```bash
pip install -r ../requirements.txt
export ALPHA_GENOME_API_KEY="<your key>"
```

The key is read from the environment and is never printed or stored."""
    ),
    CODE(
        """import os
import sys

sys.path.insert(0, os.path.abspath("../src"))

import numpy as np
import pandas as pd

from ghsr_ag2 import config, provenance

pd.set_option("display.width", 200)
pd.set_option("display.max_columns", 60)

print("API key present:", bool(os.environ.get(config.API_KEY_ENV)))
for k, v in provenance.environment_record().items():
    if k != "packages":
        print(f"{k}: {v}")
print()
for k, v in provenance.package_versions().items():
    print(f"  {k}: {v}")"""
    ),
    MD(
        """## Phase I - what are these variants, exactly?

Nothing here is taken on trust from the manuscript. Coordinates and alleles
come from NCBI dbSNP and Ensembl independently, the REF base is checked
against the GRCh38 reference sequence, and the manuscript's own PCR primers
are mapped onto the genome to confirm which region was actually sequenced."""
    ),
    CODE(
        """from ghsr_ag2 import phase1_annotation

annot, extra = phase1_annotation.build_variant_annotation()
annot[
    [
        "rsid",
        "position_grch38",
        "ref_plus_strand",
        "alt_plus_strand",
        "grch38_ref_base_observed",
        "ref_matches_reference_genome",
        "manuscript_notation",
        "manuscript_orientation_resolved",
    ]
]"""
    ),
    MD(
        """### Strand: the manuscript's allele notation is not internally consistent

*GHSR* is on the **minus** strand, so "upstream" means **higher** coordinates.
Reconciling each reported notation against GRCh38:

* **rs2922126 "A>T"** - the same substitution written in *GHSR* (minus-strand)
  orientation. Consistent. The manuscript's "forward" primer anneals to the
  plus strand, i.e. it reads in the gene's direction, which explains this.
* **rs115459563 "G>A"** - written in plus-strand orientation. Also correct,
  but the opposite convention to rs2922126.
* **rs34136054 "C>T"** - matches neither strand as written. GRCh38 carries
  **T** at this position, so the manuscript's "wild type C" is an alternate
  allele and REF/ALT are inverted. The genotype *identity* is unaffected (a
  C/T heterozygote is a T/C heterozygote), only the labelling of which allele
  is wild type.

This is a reporting inconsistency, not a genotyping error, and it is reported
as such."""
    ),
    CODE(
        """amplicon = phase1_annotation.build_amplicon_check()
for k, v in amplicon.items():
    print(f"{k}: {v}")"""
    ),
    MD(
        """The published primer pair maps to **chr3:172,449,261-172,449,659**, a
product of exactly **399 bp** - the size the manuscript reports - and it
contains all three SNPs. This is an independent confirmation of the source
study's experimental design.

Relative to the canonical *GHSR* TSS (chr3:172,448,456, minus strand), the
three variants sit **916-1,015 bp upstream**. Ensembl VEP classifies all three
as `upstream_gene_variant` for every *GHSR* transcript. None of them overlaps
an Ensembl Regulatory Build feature; the nearest is an annotated **enhancer**
(ENSR3_B3983C, chr3:172,449,561-172,450,143), 90 bp away from rs2922126, and a
**CTCF binding site** (ENSR3_53HK6KS) 349 bp away.

So "promoter/regulatory region" as described in the manuscript is
approximately right but imprecise: this is the *GHSR* **proximal upstream
region**, adjacent to - not inside - an annotated enhancer."""
    ),
    CODE(
        """genes = phase1_annotation.build_gene_window_table(config.VARIANTS[0].position)
genes[genes.biotype == "protein_coding"].head(8)[
    [
        "gene_symbol",
        "biotype",
        "strand",
        "tss",
        "signed_distance_variant_to_tss_bp",
        "variant_inside_gene_body",
    ]
]"""
    ),
    MD(
        """*GHSR* is by far the nearest gene. The next protein-coding TSS,
*TNFSF10*, is 74 kb away. Any predicted effect on another gene has to be read
with that distance in mind."""
    ),
    MD(
        """## Phase II - which cell types does the model actually support?

AlphaGenome 2 is **condition-conditioned**: a prediction is requested for one
reference condition at a time, named by a composite ontology CURIE. For every
condition and output type the API also reports whether the underlying tracks
are `observed`, `imputed_supported`, or `imputed_exploratory`. That evidence
class is carried through every table in this analysis - no score is reported
without it."""
    ),
    CODE(
        """from alphagenome.models.v2 import dna_client

client = dna_client.create(config.api_key(), timeout=120)
supported = client.supported_conditions()
print(f"{len(supported)} human reference conditions available")
supported.biosample_type.value_counts()"""
    ),
    CODE(
        """from ghsr_ag2 import phase2_conditions

priorities = phase2_conditions.audit_priorities(supported)
priorities[
    [
        "biological_priority",
        "n_conditions_in_model",
        "available_in_model",
        "n_with_observed_rna_seq",
        "conditions_selected_for_this_study",
    ]
]"""
    ),
    MD(
        """**Osteoclast-lineage cells do not exist in the AlphaGenome 2 human
condition set.** Blood monocyte is used as a declared surrogate (osteoclasts
differentiate from circulating monocytes) and is flagged as a surrogate
everywhere it appears. Hippocampal *neurons* are also unavailable: the two
hippocampal conditions are a tissue layer and an astrocyte, neither with
observed RNA-seq. These are limitations of the model's coverage, not choices,
and they are not papered over."""
    ),
    CODE(
        """conditions = pd.read_csv(
    config.TABLES / "table_S3_condition_selection_and_qc.csv"
)
conditions[
    [
        "label",
        "panel",
        "evidence_tier",
        "mode_rna_seq",
        "mode_dnase",
        "mode_chip_tf",
        "ghsr_rna_seq_active_max",
        "ghsr_predicted_expression_band",
    ]
]"""
    ),
    MD(
        """### A useful sanity check falls out of this

The model's predicted *GHSR* expression ranks **pituitary gland >
hypothalamus > testis / temporal lobe / hippocampus > other brain >>
bone and peripheral tissues**. That recapitulates known *GHSR* biology without
being told to, which is reassuring about the predictions at this locus.

It also carries a finding in its own right: *GHSR* is predicted to be
**essentially unexpressed in osteoblasts, osteocytes and mesenchymal stem
cells**. A cell-autonomous *GHSR*-regulatory mechanism in bone has very little
room to operate."""
    ),
    MD(
        """## Phase III - variant-effect prediction

All four variant configurations (three rsIDs; rs34136054 contributes two
alternate alleles, T>C and T>G) are scored identically, in all 16 conditions,
with the **complete** officially recommended scorer set. Requesting everything
rather than a hand-picked subset is deliberate: it makes selective omission or
selective reporting impossible after the fact.

Sequence context is 1 Mb, the largest supported, which the gene-mask scorers
need in order to see *GHSR* and its neighbours."""
    ),
    CODE(
        """# The full run takes ~4 minutes. Load the archived output if present.
raw_path = config.DATA_RAW / "phase3" / "phase3_raw_scores_main.parquet"
if raw_path.exists():
    scores = pd.read_parquet(raw_path)
    print(f"loaded {len(scores):,} archived per-track scores")
else:
    from ghsr_ag2 import phase3_scoring

    scores = phase3_scoring.run()

scores.groupby(["scorer_key", "modality_class"]).size()"""
    ),
    MD(
        """### Four scorers returned nothing - and that is the correct answer

`SPLICE_SITES`, `SPLICE_SITE_USAGE`, `SPLICE_JUNCTIONS` and
`POLYADENYLATION` returned well-formed **empty** results (shape
`(0, n_tracks)`): the request succeeded, and no splicing or polyadenylation
feature met the scorer's criteria. That is expected for a variant ~1 kb
upstream of the TSS and far from any splice junction. It is recorded as an
empty result, not as a failure, and not quietly dropped."""
    ),
    CODE(
        """ghsr = scores[
    (scores.scorer_key == "RNA_SEQ")
    & (scores.gene_name == config.GHSR_SYMBOL)
    & scores.is_primary_configuration
]
table = ghsr.pivot_table(
    index=["panel", "condition_label"],
    columns="rsid",
    values="raw_score",
    aggfunc="mean",
).round(5)
act = conditions.set_index("label")["ghsr_rna_seq_active_max"]
table["GHSR_predicted_activity"] = [act.get(i[1], np.nan) for i in table.index]
table.sort_values("GHSR_predicted_activity", ascending=False)"""
    ),
    MD(
        """### Reading this table

The units are log2 fold change of predicted signal, ALT minus REF. The largest
value anywhere is about **0.027**, i.e. a **~1.9%** change in predicted *GHSR*
expression. For comparison, a variant that genuinely disrupted a promoter
would be expected to move this by a large fraction of a log2 unit.

Three things stand out:

1. **Magnitude is uniformly very small** - nothing approaches a size that a
   reporter assay could resolve.
2. **Effect size tracks baseline expression.** The largest effects appear in
   exactly the conditions where *GHSR* is predicted to be most expressed.
   A multiplicative effect on a near-zero baseline is near-zero, so this is
   largely arithmetic, not condition-specific regulatory logic. Quantified
   below.
3. **Bone conditions are flat** - |log2FC| < 0.0004 for all three variants in
   osteoblasts, osteocytes and mesenchymal stem cells."""
    ),
    CODE(
        """from scipy.stats import spearmanr

# Computed directly from the per-track scores, restricted to the primary
# variant configuration of each rsID. rs34136054 is tri-allelic and its
# second alternate allele (T>G) is summarised separately below: mixing two
# different substitutions into one per-rsID number would be misleading,
# because they can differ in sign.
for rsid in config.RSIDS:
    per_cond = (
        ghsr[ghsr.rsid == rsid].groupby("condition_label")["raw_score"].mean()
    )
    joined = pd.DataFrame(
        {"lfc": per_cond, "activity": act.reindex(per_cond.index)}
    ).dropna()
    rho = spearmanr(joined.lfc.abs(), joined.activity).statistic
    print(
        f"{rsid:<13} rho(|log2FC|, GHSR activity) = {rho:+.2f}   "
        f"ALT raises in {int((per_cond > 0).sum())}, "
        f"lowers in {int((per_cond < 0).sum())} "
        f"of {len(per_cond)} conditions"
    )"""
    ),
    MD(
        """**No variant has a consistent direction of effect across conditions.**
That is the opposite of what a real, shared cis-regulatory mechanism would
look like, and it argues directly against the cell-type-specific regulatory
hypothesis as stated."""
    ),
    MD(
        """## The evidence-class problem at this locus

Before reading anything into the chromatin or transcription-factor
predictions, it is worth seeing how much of the TF signal is imputed."""
    ),
    CODE(
        """tf = scores[scores.scorer_key == "CHIP_TF"]
print(tf.imputation_mode.value_counts().to_string())
print()
observed_tf = tf[tf.imputation_mode == "observed"]
print("Observed TF ChIP tracks available at this locus:")
print(
    observed_tf.groupby(["transcription_factor", "condition_label"])
    .size()
    .to_string()
)"""
    ),
    MD(
        """Only **60 of the 43,728** TF ChIP values for the three variants rest on observed data, and they cover
just a handful of factors - dominated by **CTCF**. Any statement about
transcription-factor occupancy at this locus that relies on the other 43,668
rows is resting on imputed-exploratory tracks, the model's weakest evidence
class. The two signals below are reported *because* they survive restriction
to observed tracks."""
    ),
    CODE(
        """ctcf = observed_tf[observed_tf.transcription_factor == "CTCF"]
for rsid in config.RSIDS:
    s = ctcf[ctcf.rsid == rsid]
    print(
        f"{rsid:<13} n={len(s):2d} observed CTCF tracks, "
        f"{int((s.raw_score < 0).sum()):2d} negative, "
        f"mean={s.raw_score.mean():+.5f}, min={s.raw_score.min():+.5f}"
    )

print()
hist = scores[
    (scores.scorer_key == "CHIP_HISTONE")
    & (scores.imputation_mode == "observed")
    & scores.is_primary_configuration
]
top = hist.reindex(hist.raw_score.abs().sort_values(ascending=False).index)
top.head(6)[["rsid", "condition_label", "histone_mark", "raw_score"]]"""
    ),
    MD(
        """### Two hypotheses worth testing at the bench - and nothing more

* **rs2922126 lowers predicted CTCF ChIP signal in 12 of 12 observed CTCF
  tracks**, across 8 different conditions (mean log2FC -0.010, strongest
  -0.025 in osteoblast). This is the most directionally consistent result in
  the whole study, and it is on observed tracks. Note carefully: the 501 bp
  scoring window centred on rs2922126 does **not** reach the annotated CTCF
  element 349 bp away, so this concerns predicted local CTCF occupancy near
  the variant, not that element. **No motif analysis supports a disrupted CTCF
  binding site, and we do not claim one.** The effect is also *shared* across
  cell types rather than specific to any - again pointing away from the
  cell-type-specificity hypothesis.
* **rs34136054 raises predicted active-promoter marks in brain** - H3K4me3
  (+0.032), H3K9ac (+0.027) and H3K27ac (+0.027) in temporal lobe, on observed
  tracks, with the same direction in brain and hippocampus.

Both are small. Both are predictions. Either could be tested directly - see
the Discussion for the experiments that would do it."""
    ),
    MD(
        """## Phase IV - is "very small" actually small?

A log2FC of 0.027 means nothing in isolation. The reference frame used here is
**in silico mutagenesis**: every possible single-base substitution across the
whole 399 bp sequenced amplicon (1,197 variants) is scored in the same
condition with the same scorer. Because ISM enumerates all three substitutions
at every position in the interval, **the three studied SNPs are themselves
members of that set** - so their scores are read from the same run, under an
identical sequence context, with no cross-run normalisation.

This also serves as a positive control for model sensitivity. If some
positions in this amplicon produce large predicted effects while the studied
SNPs do not, then a near-zero score is an informative negative rather than a
silent pipeline failure."""
    ),
    CODE(
        """ism_path = (
    config.DATA_RAW / "phase4" / "phase4_ism_reference_distribution.parquet"
)
if ism_path.exists():
    ism = pd.read_parquet(ism_path)
    pct = pd.read_csv(
        config.TABLES / "table_S5_ism_reference_distribution.csv"
    )
    display(
        pct[
            [
                "condition_label",
                "scorer_key",
                "rsid",
                "observed_score",
                "ism_median_abs",
                "ism_max_abs",
                "percentile_of_abs_effect_in_amplicon",
            ]
        ]
    )
else:
    print(
        "Run `python -m ghsr_ag2.phase4_qc --ism` to generate the ISM "
        "reference distribution (~25 min)."
    )"""
    ),
    MD(
        """## Phase V - integrating three kinds of evidence, without mixing them

Three separate things are known about these variants, and they are not
interchangeable:

| | What it is | What it establishes | What it does **not** establish |
|---|---|---|---|
| **Level 1** | AlphaGenome 2 prediction (this analysis) | predicted cis-regulatory consequence | measured expression, phenotype, diagnosis |
| **Level 2** | PCR + Sanger sequencing (prior study, n=100) | the genotype was correctly read | any effect on gene regulation |
| **Level 3** | Case-control association (prior study, n=100) | exploratory, underpowered association | causation, mechanism, or an AD diagnosis |

Sanger sequencing validates **genotype identification**. It is not functional
validation of a predicted regulatory effect, and nothing in this analysis
should be described that way."""
    ),
    CODE(
        """integ_path = config.TABLES / "table_3_integrated_evidence.csv"
if integ_path.exists():
    integ = pd.read_csv(integ_path)
    display(integ.set_index("rsid").T)
else:
    print("Run `python -m ghsr_ag2.phase5_integration` first.")"""
    ),
    MD(
        """## Conclusion

Within the limits of what a sequence model can say:

1. **No variant shows a predicted regulatory effect of plausibly functional
   magnitude.** The maximum predicted change in *GHSR* expression across 16
   conditions, 4 variant configurations and 162,144 track-level scores is
   ~1.9%.
2. **The hypothesis of cell-type-specific regulation is not supported.** Effect
   magnitude tracks baseline *GHSR* expression (rho 0.63-0.83), no variant has
   a consistent direction across conditions, and the single most consistent
   signal (rs2922126 -> CTCF) is shared across cell types rather than specific
   to any.
3. **The bone arm of the hypothesis is the weakest.** *GHSR* is predicted to be
   essentially unexpressed in osteoblasts, osteocytes and mesenchymal stem
   cells, and predicted effects there are the smallest in the panel.
4. ***GHSR* remains the only plausible target gene.** The nearest other
   protein-coding TSS is 74 kb away and predicted effects on neighbours are no
   larger than on *GHSR*.
5. **This is concordant with, but logically independent of, the prior study's
   null findings.** Both are null; they are null about different things. The
   prior association analysis had 6.62% achieved power for rs2922126, so it
   could not have detected a small effect either way.

A negative result is a result. What it does not license is the inverse claim:
these predictions do not demonstrate that the variants are biologically
inert, only that this model, at this locus, in these conditions, expects
little."""
    ),
]


def build() -> nbf.NotebookNode:
    nb = nbf.v4.new_notebook(cells=CELLS)
    nb.metadata = {
        "kernelspec": {
            "display_name": "Python 3",
            "language": "python",
            "name": "python3",
        },
        "language_info": {"name": "python", "version": "3.13"},
    }
    return nb


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--execute", action="store_true")
    args = ap.parse_args()

    nb = build()
    OUT.parent.mkdir(parents=True, exist_ok=True)
    nbf.write(nb, OUT)
    print(f"wrote {OUT} ({len(nb.cells)} cells)")

    if args.execute:
        from nbclient import NotebookClient

        client = NotebookClient(
            nb,
            timeout=3600,
            kernel_name="python3",
            resources={"metadata": {"path": str(OUT.parent)}},
        )
        client.execute()
        nbf.write(nb, OUT)
        print(f"executed and re-wrote {OUT}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
