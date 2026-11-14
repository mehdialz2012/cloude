# Methods subsection (drafted for *Scientific Reports*)

> Intended as a new subsection of the existing Methods, placed after
> "Genetic Analysis" and before "Statistical Analysis". See
> `MANUSCRIPT_INTEGRATION_RECOMMENDATIONS.md` for placement rationale.

---

## In silico variant-effect prediction

### Study design and timing

A retrospective, complementary in silico analysis of the three genotyped
*GHSR* variants was performed after the genotyping, biochemical and
densitometric analyses described above had been completed. The computational
analysis did not inform the selection of variants, the design of the
experimental work, or the statistical analysis plan, and is reported here as a
separate, hypothesis-generating layer of evidence. No participant-level
genetic, clinical, biochemical or densitometric data were used in, or
transmitted during, the computational analysis; only public variant
coordinates and public ontology identifiers were submitted to any external
service.

### Variant definition and coordinate verification

Genomic coordinates and alleles for rs2922126, rs115459563 and rs34136054
were retrieved independently from NCBI dbSNP (build 157, accessed
14 November 2026) and from the Ensembl REST API (GRCh38, accessed
14 November 2026), and were required to agree. Reference alleles were verified
base-by-base against the GRCh38 reference sequence retrieved from the Ensembl
`sequence/region` endpoint. Variant consequences, transcript annotation and
distances to the *GHSR* transcription start site were obtained from Ensembl
VEP; overlap with annotated regulatory elements was assessed against the
Ensembl Regulatory Build; and population allele frequencies were taken from
gnomAD genomes and 1000 Genomes phase 3 as served by Ensembl.

All coordinates are reported on the GRCh38 plus strand. Because *GHSR*
(ENSG00000121853) is transcribed from the minus strand, upstream positions
correspond to higher coordinates. Allele notations in the present manuscript
were reconciled against the reference genome individually, and residual
inconsistencies are reported in Supplementary Table S1 rather than silently
harmonised.

To confirm which genomic interval was experimentally sequenced, the published
primer pair (Table 1) was mapped onto the GRCh38 reference sequence. The
resulting product spans chr3:172,449,261-172,449,659 (399 bp), matching the
reported product size, and contains all three variants.

The alternate allele of rs34136054 could not be determined unambiguously from
the sequencing records available for this analysis, and dbSNP lists this site
as tri-allelic (T>C and T>G). Both alternate alleles were therefore scored
independently, and both are reported.

### Model and sequence context

Predictions were generated with AlphaGenome 2 (Google DeepMind) through the
official Python SDK `alphagenome` version 2.0.0 (`alphagenome.models.v2`),
accessed at `gdmscience.googleapis.com:443` on 14 November 2026 with the
server's current default model checkpoint. AlphaGenome 2 is conditioned on a
biological condition supplied as a composite ontology CURIE, and returns, for
each requested output type, predictions for the reference and alternate
alleles together with metadata declaring whether the underlying tracks are
observed experimental data, imputed with support from related observed data
(`imputed_supported`), or imputed without closely related observed data
(`imputed_exploratory`).

Each variant was centred in a 1,048,576 bp (1 Mb) sequence context, the
largest supported, which is required for the gene-mask scorers to evaluate
*GHSR* together with all genes in the surrounding window. As a sensitivity
analysis, all predictions were repeated at 16,384 bp, 131,072 bp and
524,288 bp contexts.

### Condition selection

The complete set of 2,799 supported human reference conditions was retrieved
from the API and archived. Conditions were pre-specified to cover the
biological priorities implied by the phenotypes measured in the present study:
hippocampal neurons, hypothalamic neurons, other brain neuronal populations,
pituitary cells, osteoblasts, osteocytes and osteoclast-lineage cells. Where a
priority was not represented in the model's condition set, this was recorded
rather than substituted silently.

Sixteen conditions were selected (Supplementary Table S3): seven
brain/neuronal (hypothalamus, layer of hippocampus, hippocampal astrocyte,
brain, temporal lobe, brain neuron, glutamatergic neuron), one pituitary
(pituitary gland), four bone-related (bone-tissue osteoblast, osteocyte,
mesenchymal stem cell, blood monocyte) and four fully observed context tissues
providing a specificity reference and corresponding to other phenotypes
measured in the present study (stomach, testis, thyroid gland, ovary).

Two pre-specified priorities were not fully available. No osteoclast or
osteoclast-precursor condition exists in the AlphaGenome 2 human condition
set; blood monocyte was used as an explicitly declared surrogate, on the basis
that osteoclasts differentiate from circulating monocytes, and is flagged as a
surrogate wherever it appears. No hippocampal *neuron* condition exists; the
two available hippocampal conditions (a tissue layer and an astrocyte) were
both included, neither of which has observed RNA-seq data. These limitations
are reported in Supplementary Table S4.

To confirm that *GHSR* regulation is interpretable in each selected condition,
predicted absolute *GHSR* activity was measured using the gene-mask
active-allele RNA-seq scorer and expressed relative to the highest value in
the panel. Because the scorer returns values in the model's own arbitrary
units, no absolute expression threshold was applied; conditions were banded
descriptively relative to the panel maximum.

### Variant scoring

For every variant configuration and every condition, the complete set of 22
officially recommended human variant scorers was requested in a single call
via `variant_scorers.get_recommended_scorers`. The full recommended set was
used, rather than a selected subset, so that no modality could be omitted or
chosen after the results were seen. Scorers were classified before the
analysis was run into: primary modalities mechanistically plausible for a
non-coding upstream single-nucleotide variant (RNA-seq, CAGE, DNase, ATAC,
single-cell ATAC, transcription-factor ChIP, histone ChIP); secondary
modalities computed for completeness but not expected a priori (splice sites,
splice-site usage, splice junctions, polyadenylation, contact maps,
methylation, single-cell gene expression); and `*_ACTIVE` scorers, which
return the absolute activity of the more active allele rather than an allelic
difference and were used only for quality control, never as effect sizes.

Gene-level expression effects were scored with the gene-mask log-fold-change
scorer (difference of means, alternate minus reference); chromatin
accessibility, CAGE, transcription-factor ChIP and histone ChIP with
centre-mask scorers (501 bp window for DNase, ATAC, CAGE and TF ChIP;
1,201 bp for single-cell ATAC; 2,001 bp for histone ChIP) aggregating as
log2 of the summed signal. Scores are therefore log2 fold changes of predicted
signal and are signed, with positive values indicating higher predicted signal
for the alternate allele.

All 162,144 per-track scores were retained unaggregated and are provided as
the primary computational data; all reported summaries are derived from that
archive. Base-resolution reference and alternate prediction tracks were
additionally retrieved with `predict_variant` at a 131,072 bp context for six
representative conditions.

### Reference distribution for effect magnitude

AlphaGenome scores have no intrinsic scale and no associated P-value,
confidence interval or significance threshold; none was computed, and none
should be inferred. To provide an interpretable frame of reference, in silico
mutagenesis was performed across the entire 399 bp sequenced amplicon,
scoring all three possible single-base substitutions at every position
(1,197 variants) in four conditions (pituitary gland, hypothalamus,
bone-tissue osteoblast and stomach) with the RNA-seq, CAGE, DNase and ATAC
scorers. Because in silico mutagenesis enumerates every substitution in the
interval, the three studied variants are themselves members of this set, so
their scores were read from the same run under an identical sequence context,
scorer and condition, without cross-run normalisation. Each variant's
predicted effect is reported as a percentile of the absolute-effect
distribution of this local reference set. This is an empirical reference
distribution, not a statistical test. The same analysis provides a positive
control for model sensitivity at this locus.

### Quality control

The following checks were pre-specified and are reported in full. Reference
alleles were re-verified against GRCh38 after scoring, and the variant strings
transmitted to the API were compared against those returned. Completeness was
verified by confirming that every variant configuration was scored in every
condition with an identical scorer count. Every summary was recomputed using
only `observed` tracks, and separately using only `imputed_exploratory`
tracks, so that no conclusion could rest on imputed data without that being
visible. Aggregated centre-mask DNase and ATAC scores were independently
recomputed from the base-resolution reference and alternate tracks and checked
for agreement in sign and order of magnitude. Scorers that returned no rows
were investigated individually rather than dropped.

Effect magnitudes are described using fixed descriptive bands defined before
the results were inspected (negligible, |log2 fold change| < 0.01; very small,
0.01-0.1; moderate, 0.1-0.5; large, >= 0.5). These are labels for predicted
magnitude and are explicitly not significance categories.

### Reproducibility

Software versions, the model family and API endpoint, the genome build, the
verified variant coordinates, condition identifiers, scorer definitions,
requested output modalities, parameter settings, database sources with access
dates, execution timestamps and every technical failure or empty result are
recorded automatically for each analysis step. The analysis code, all raw
model outputs at per-track resolution, derived tables, figure source data and
an executable notebook are available at <repository URL>. The AlphaGenome API
key was supplied through an environment variable and is not contained in the
code, outputs or logs.

### Statistical treatment of the computational results

No inferential statistics were applied to the model predictions, and no
P-values, confidence intervals or functional-significance thresholds were
derived from them. Two descriptive, clearly labelled summaries are reported:
the number of conditions in which the alternate allele raised versus lowered
predicted *GHSR* expression (a test of directional consistency), and the
Spearman rank correlation between the absolute predicted effect and the
predicted *GHSR* expression level across conditions (used to assess whether
apparent cell-type differences simply follow baseline expression). These
correlations describe the structure of the model's own output across 16
conditions and are not association tests involving participant data.
