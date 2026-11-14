# AlphaGenome 2 functional analysis of three GHSR upstream variants

A retrospective, complementary **in silico** investigation of the three *GHSR*
polymorphisms genotyped in *"Ghrelin Receptor Gene Polymorphisms, Cognitive
Screening Status, and Osteoporosis-Related Parameters: An Exploratory
Case-Control Study"* (Saeidi et al., `ARTICLE__revised-v2.docx`).

The computational analysis was performed **after** the genotyping was complete.
It did not guide the experimental design, the variant selection, or the
statistical analysis of the source study, and nothing in this repository should
be described as having done so.

| Variant | GRCh38 (chr3) | REF>ALT (+ strand) | Position relative to *GHSR* |
|---|---|---|---|
| rs2922126 | 172,449,471 | T>A | 1,015 bp upstream of the TSS |
| rs115459563 | 172,449,444 | G>A | 988 bp upstream of the TSS |
| rs34136054 | 172,449,372 | T>C *and* T>G | 916 bp upstream of the TSS |

*GHSR* (ENSG00000121853) is on the **minus** strand, so "upstream" means higher
coordinates. All three variants are `upstream_gene_variant` for every *GHSR*
transcript and lie inside the 399 bp amplicon the source study sequenced.

## Headline result

AlphaGenome 2 predicts only **very small** allelic effects for all three
variants, in every one of the 16 conditions examined: the largest predicted
change in *GHSR* expression anywhere in the panel is |log2FC| = 0.027
(a ~1.9% change in predicted signal). No variant produced a sign-consistent
direction of effect across all conditions. The effect magnitudes track how
strongly *GHSR* is predicted to be expressed in each condition
(Spearman rho 0.63-0.83), so the apparent cell-type pattern is largely a
consequence of baseline expression rather than condition-specific regulatory
logic.

Two directionally consistent signals on **observed** (non-imputed) tracks are
reported as hypotheses for wet-lab testing, not as established effects:
rs2922126 lowered predicted CTCF ChIP signal in 12 of 12 observed CTCF tracks
across 8 conditions, and rs34136054 raised predicted active-promoter histone
marks (H3K4me3, H3K9ac, H3K27ac) in brain tissue. Both are of small magnitude.

Predicted effects in osteoblasts, osteocytes and mesenchymal stem cells were
the smallest in the panel. *GHSR* itself is predicted to be essentially
unexpressed in those cell types, which removes the straightforward
cell-autonomous route by which these variants could affect bone.

Full reasoning, caveats and limitations: [`docs/REPORT_EN.md`](docs/REPORT_EN.md).

## Installation

Requires Python >= 3.11 (3.13.18 was used) and an AlphaGenome API key.

```bash
git clone <this repository>
cd alphagenome_ghsr
python -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
```

Provide the API key through the environment. **Never** put it in a file:

```bash
export ALPHA_GENOME_API_KEY="<your key>"
```

Register for a key at <https://deepmind.google.com/science/alphagenome>.

### Network access

The pipeline needs outbound HTTPS to:

| Host | Used for |
|---|---|
| `gdmscience.googleapis.com:443` | AlphaGenome 2 gRPC API |
| `rest.ensembl.org` | Ensembl coordinates, VEP, sequence, regulatory features |
| `api.ncbi.nlm.nih.gov` | NCBI dbSNP RefSNP records |

Only public variant coordinates and public ontology identifiers are ever sent
to any external service. No participant-level, clinical, biochemical or
densitometric data leaves the machine.

## Reproducing the analysis

```bash
export PYTHONPATH=src

python -m ghsr_ag2.phase1_annotation       # ~1 min   coordinates & annotation
python -m ghsr_ag2.phase2_conditions       # ~2 min   condition selection & QC
python -m ghsr_ag2.phase3_scoring          # ~4 min   variant scoring (main)
python -m ghsr_ag2.phase3b_tracks          # ~2 min   REF/ALT track profiles
python -m ghsr_ag2.phase4_qc --ism --seqlen  # ~30 min  robustness & QC
python -m ghsr_ag2.phase4b_stability       # <1 min   context stability
python -m ghsr_ag2.phase5_integration      # <1 min   integrated evidence
python -m ghsr_ag2.phase6_figures          # ~1 min   figures
python scripts/export_csv.py               # <1 min   CSV exports
```

or in one step:

```bash
bash scripts/run_all.sh
```

Phase I caches every raw API response under `data/raw/phase1/`, so re-runs are
offline and deterministic. Delete that directory to force a fresh download.

The executable notebook `notebooks/AlphaGenome2_GHSR_analysis.ipynb` walks
through the same analysis with narrative commentary.

## Layout

```
config/                     (reserved for user overrides)
src/ghsr_ag2/
  config.py                 variants, conditions, scorer sets, paths
  provenance.py             version / environment / database recording
  phase1_annotation.py      Phase I   coordinate & annotation verification
  phase2_conditions.py      Phase II  condition selection and data quality
  phase3_scoring.py         Phase III variant-effect scoring
  phase3b_tracks.py         Phase III-b REF vs ALT base-resolution tracks
  phase4_qc.py              Phase IV  robustness, ISM reference distribution
  phase5_integration.py     Phase V   three-level evidence integration
  phase6_figures.py         Phase VI  figures
  style.py                  validated figure palette and styling
data/raw/                   unmodified API responses and model outputs
data/processed/             intermediate derived artefacts
results/tables/             publication tables (CSV)
results/figures/            figures (PDF, SVG, 600 dpi PNG)
docs/                       report, manuscript sections, legends, citations
notebooks/                  executable notebook
MANUSCRIPT_INTEGRATION_RECOMMENDATIONS.md
```

## Outputs

### Tables

| File | Contents |
|---|---|
| `table_S1_variant_annotation.csv` | Phase I annotation, strand reconciliation, allele frequencies |
| `table_S2_genes_in_1mb_window.csv` | All genes within the 1 Mb scoring window |
| `table_S3_condition_selection_and_qc.csv` | Conditions, per-output-type evidence class, predicted *GHSR* activity |
| `table_S4_priority_availability_audit.csv` | Which pre-specified cell types the model does and does not support |
| `table_S5_ism_reference_distribution.csv` | Each SNP's percentile within the local ISM distribution |
| `table_S6_observed_vs_imputed_tracks.csv` | Every summary recomputed on observed tracks only |
| `table_S7_track_consistency.csv` | Aggregated scores vs base-resolution tracks |
| `table_S8_sequence_length_sensitivity.csv` | Scores at 16 kb / 131 kb / 524 kb / 1 Mb context |
| `table_S9_context_stability_ghsr_expression.csv` | Sign stability of the *GHSR* expression effect across contexts |
| `table_S10_context_stability_observed_ctcf.csv` | Sign stability of the observed CTCF ChIP effect |
| `table_2_prediction_summary_by_variant.csv` | Per-variant predicted-effect summary |
| `table_3_integrated_evidence.csv` | Three-level integrated evidence comparison |

`data/raw/phase3/phase3_raw_scores_main.parquet` holds the **unaggregated**
per-track model output (162,144 rows). Everything in `results/tables/` is
derived from it. `data/processed/csv_exports/` carries the same raw outputs as
gzipped CSV for tools that cannot read Parquet (`python scripts/export_csv.py`).

Base-resolution REF/ALT prediction arrays (~590 MB) are **not** committed; they
are regenerated by `python -m ghsr_ag2.phase3b_tracks`. The per-track
aggregated scores and the track index are committed.

### Figures

| Figure | Contents |
|---|---|
| `figure1_ghsr_locus` | *GHSR* transcripts, TSS, regulatory features, amplicon, the three SNPs |
| `figure2_predicted_effects` | Effect magnitudes per modality, against the ISM reference distribution |
| `figure3_celltype_specificity` | Signed effects across conditions and modalities, two scales |
| `figure4_integrated_evidence` | The three evidence levels and the inferences they do not license |

## Interpretation guardrails

These are the rules the analysis was written to, and they apply to anyone
reusing these outputs:

* An AlphaGenome score is a **model prediction**, not a measurement. No
  P-value, confidence interval or significance threshold is attached to any
  score in this repository, because none exists.
* The effect bands ("negligible", "very small", "moderate", "large") are
  descriptive labels for predicted magnitude, defined in
  `phase5_integration.classify_effect`. They are not significance categories.
* PCR and Sanger sequencing in the source study validate **genotype
  identification**. They do not validate any predicted regulatory effect.
* The source study measured a cognitive **screening** score (6-CIT), not a
  clinical diagnosis of Alzheimer's disease. Nothing here is evidence about
  diagnosed Alzheimer's disease.
* Tracks are labelled `observed`, `imputed_supported` or
  `imputed_exploratory`. Conclusions that depend on imputed-exploratory tracks
  are flagged as such; at this locus, 43,668 of the 43,728 TF ChIP values for
  the three variants are imputed-exploratory.
* Predictions are not clinical or diagnostic evidence.

## Citation

See [`docs/CITATIONS.md`](docs/CITATIONS.md) for the model, software and
database citations, with access dates.

## License

Analysis code: Apache-2.0 (matching the AlphaGenome SDK). The source
manuscript and its data are the property of their authors.
