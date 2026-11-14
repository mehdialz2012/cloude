# Figure legends

Each figure is provided as vector PDF and SVG and as a 600 dpi PNG in
`results/figures/`. The exact values plotted are in
`results/tables/figure<N>_source_data.csv`.

---

## Figure 1. Genomic organisation of the *GHSR* locus and the three genotyped variants.

**(A)** The *GHSR* locus on chromosome 3 (GRCh38). The three annotated *GHSR*
transcripts (ENSG00000121853) are drawn with exons as filled boxes and introns
as connecting lines; arrows give the direction of transcription. *GHSR* is
transcribed from the **minus** strand, so upstream positions lie at **higher**
genomic coordinates. The vertical line marks the transcription start site of
the canonical transcript GHSR-201 (ENST00000241256) at chr3:172,448,456. The
399 bp amplicon sequenced in the source study is marked, and lies entirely
upstream of the gene.

**(B)** Enlargement of the upstream region. The amplicon was mapped onto GRCh38
from the primer pair published in Table 1 and spans
chr3:172,449,261-172,449,659; the mapped product size reproduces the reported
399 bp exactly and contains all three variants. Coloured triangles mark the
three variants, labelled with their GRCh38 plus-strand reference and alternate
alleles, their position, and their distance upstream of the transcription
start site. Elements of the Ensembl Regulatory Build in the window are shown
above: an annotated enhancer (ENSR3_B3983C, chr3:172,449,561-172,450,143) and
a CTCF binding site (ENSR3_53HK6KS, chr3:172,449,820-172,449,841). All three
variants are classified by Ensembl VEP as `upstream_gene_variant` for every
*GHSR* transcript, and none overlaps an annotated regulatory element, the
closest approach being 90 bp between rs2922126 and the enhancer.

---

## Figure 2. Predicted allelic effects are small in absolute terms and unremarkable within the local distribution of mutational effects.

**(A)** Largest absolute predicted effect for each variant and each primary
output modality, taken across all 16 conditions and all tracks within each
modality; values are printed beside each bar and the axis is logarithmic.
Effects on *GHSR* expression are gene-mask log-fold-change scores;
accessibility, CAGE and ChIP effects are centre-mask log2-fold-change scores.
The dashed line marks |log2 fold change| = 0.1 (a ~7% change in predicted
signal), included as a reference magnitude; no variant reaches it for any
modality. Because these are maxima over many tracks and conditions, they
represent the most extreme value observed anywhere rather than a typical
effect.

**(B)** The same variants placed inside their own local reference
distribution. In silico mutagenesis was used to score all three possible
single-base substitutions at every position of the 399 bp sequenced amplicon
(1,197 variants). Violins show the distribution of absolute predicted effects
for that complete substitution set on a log10 scale; filled circles mark the
three studied variants, which are themselves members of that set, so violins
and circles come from the same run and share sequence context, scorer and
condition (pituitary gland is shown; the other three conditions analysed are
in Supplementary Table S5). Circles are offset vertically within each row so
that all three variants remain visible.

Across the amplicon the model predicts effects 4.7- to 95-fold larger than the
largest effect produced by any of the three studied variants, confirming that
it is responsive at this locus and that the near-null results are informative
negatives rather than a failure of the analysis. Across 48 condition x
endpoint combinations none of the three variants exceeded the 90th percentile
of its local distribution; median percentiles were 43.7 (rs2922126), 42.1
(rs34136054) and 18.9 (rs115459563).

No P-values, confidence intervals or significance thresholds are shown,
because AlphaGenome scores are model predictions for which none is defined.

---

## Figure 3. Predicted allelic effects across conditions and modalities.

Signed predicted effects (log2 fold change, alternate minus reference) for
each variant (column groups) and each primary output modality (columns) in
each of the 16 selected conditions (rows). Rows are grouped by biological
panel, separated by horizontal rules and labelled at the left, and ordered
within each panel by predicted *GHSR* expression. Each cell shows the signed
score of the track with the largest absolute effect within that modality; red
indicates higher predicted signal for the alternate allele, blue lower, and
the neutral grey midpoint indicates no change. A grey dot marks cells whose
peak track is `imputed_exploratory`, the model's weakest evidence class; at
this locus 43,668 of the 43,728 transcription-factor ChIP values for the
three variants fall into that class, so that column is almost entirely imputed.

The colour scale is **fixed** at +/- 0.1 log2 fold change (a ~7% change in
predicted signal) rather than scaled to the data, so that the near-neutral
colour of nearly every cell reflects the true magnitude of the effects. The
largest absolute effect anywhere in the matrix is 0.051, stated on the figure.

The right-hand column gives the predicted *GHSR* expression level of each
condition as a percentage of the panel maximum (reached in the pituitary
gland). It is included because it is the main explanation for what structure
the matrix does contain: the magnitude of the predicted effect tracks baseline
*GHSR* expression (Spearman rho 0.63-0.83), and the bone conditions, in which
*GHSR* is predicted to be essentially unexpressed, show the smallest predicted
effects in the panel.

---

## Figure 4. The three levels of evidence and the inferences they do not license.

Schematic separating the three kinds of evidence that bear on these variants.
**Level 1** is the computational regulatory prediction reported here
(AlphaGenome 2; 16 conditions, 22 scorers, 4 variant configurations, GRCh38);
it establishes a predicted cis-regulatory consequence and does not establish
measured expression, phenotype or diagnosis. **Level 2** is the experimental
detection of genotypes in the source study (PCR and Sanger sequencing,
n = 100, 399 bp amplicon, 100% call rate, Hardy-Weinberg equilibrium met in
controls); it establishes that the genotype at each position was read
correctly and establishes nothing about gene regulation. **Level 3** is the
statistical association analysis of the source study (case-control design,
n = 100, groups defined by 6-CIT screening score); no association survived
false-discovery-rate correction, achieved power for rs2922126 was 6.62%, the
rs115459563 regression models did not converge, and rs34136054 was observed in
a single participant, so it establishes nothing beyond exploratory,
underpowered association, and in particular establishes neither causation,
mechanism, nor a clinical diagnosis of Alzheimer's disease.

Barred arrows mark two inferences the study design does not support: a null
computational prediction does not refute a reported association, and
experimental genotyping does not functionally validate a predicted regulatory
effect. The computational analysis was performed after the genotyping and did
not guide the experimental design.
