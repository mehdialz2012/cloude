# Results subsection (drafted for *Scientific Reports*)

> Intended as the final Results subsection, after the logistic-regression
> results. Placing it last reinforces that the computational analysis followed
> the experimental work. See `MANUSCRIPT_INTEGRATION_RECOMMENDATIONS.md` §3.2.

---

## In silico prediction of regulatory effects

### Variant localisation and annotation

To assess whether the genotyped variants could plausibly affect *GHSR*
transcription, we performed a retrospective in silico analysis using
AlphaGenome 2. Coordinates and alleles retrieved independently from NCBI dbSNP
(build 157) and Ensembl were in complete agreement, and the reference allele of
each variant matched the GRCh38 reference sequence (Supplementary Table S1).
Mapping the primer pair used for amplification (Table 1) onto GRCh38 yielded a
product spanning chr3:172,449,261-172,449,659, 399 bp in length, reproducing
the observed product size and containing all three variants.

*GHSR* (ENSG00000121853) is transcribed from the minus strand, so the three
variants, which lie at higher coordinates than the gene, are upstream of it:
916 bp (rs34136054), 988 bp (rs115459563) and 1,015 bp (rs2922126) from the
transcription start site of the canonical transcript. Ensembl VEP classified
all three as `upstream_gene_variant` for every *GHSR* transcript, and none lies
within a *GHSR* exon or intron. None overlapped an element of the Ensembl
Regulatory Build; the closest annotated elements were an enhancer
(ENSR3_B3983C, 90 bp from rs2922126) and a CTCF binding site (ENSR3_53HK6KS,
349 bp from rs2922126) (Figure 1). *GHSR* was by a wide margin the nearest
gene, the next protein-coding transcription start site (*TNFSF10*) lying 74 kb
away.

Two discrepancies with the allele notation used elsewhere in this manuscript
were identified and are reported for transparency. The notation for rs2922126
(A>T) describes the substitution in *GHSR* (minus-strand) orientation, whereas
that for rs115459563 (G>A) is given on the plus strand. For rs34136054, GRCh38
carries T at chr3:172,449,372, so the reported "C>T" inverts reference and
alternate alleles; dbSNP additionally lists this site as tri-allelic (T>C and
T>G), and gnomAD indicates that the alternate allele segregating in populations
is G. Because genotype identity is unaffected by allele labelling, no
experimental result is altered, but both alternate alleles of rs34136054 were
scored and both are reported below.

Observed allele frequencies were broadly consistent with gnomAD estimates for
the Middle Eastern population for rs2922126 (0.24 observed vs 0.276) and
rs34136054 (0.005 vs 0.017). The carrier frequency of rs115459563 (7 of 50
controls; allele frequency 0.035) was approximately five-fold higher than the
gnomAD Middle Eastern estimate of 0.0068.

### Condition selection and predicted *GHSR* expression

Of the 2,799 human reference conditions supported by AlphaGenome 2, 16 were
selected to span the biological priorities implied by the measured phenotypes
(Supplementary Table S3). Two pre-specified priorities were not fully
available: no osteoclast or osteoclast-precursor condition exists in the model's
condition set, so blood monocyte was used as a declared surrogate, and no
hippocampal neuron condition exists, so the two available hippocampal
conditions were used instead (Supplementary Table S4).

Predicted *GHSR* expression across the panel was highest in the pituitary
gland, followed by the hypothalamus (41% of the pituitary value), testis
(24%), temporal lobe (23%) and hippocampus (22%), and was lowest in the three
bone-related conditions: bone-tissue osteoblast, mesenchymal stem cell and
osteocyte each reached 1% or less of the panel maximum. This ordering
reproduces the established tissue distribution of *GHSR* without having been
specified in advance.

### Predicted allelic effects

For each of the four variant configurations in each of the 16 conditions, the
complete set of 22 officially recommended variant scorers was requested,
yielding 162,144 track-level scores. All requests completed successfully. Four
scorers (splice sites, splice-site usage, splice junctions, polyadenylation)
returned valid empty result sets, indicating that no splicing or
polyadenylation feature was affected, as expected for variants approximately
1 kb upstream of the transcription start site and distant from any splice
junction.

Predicted allelic effects were very small throughout (Figure 2A). The largest
predicted change in *GHSR* expression anywhere in the panel was a log2 fold
change of 0.027 for rs2922126 (a 1.9% change in predicted signal), with
corresponding maxima of 0.016 for rs115459563 and 0.025 for rs34136054. The
largest effects on any modality were 0.051 for chromatin accessibility, 0.050
for transcription-factor ChIP and 0.032 for histone ChIP; no value in any
modality approached a log2 fold change of 0.1 (~7% change). Predicted effects
on neighbouring genes were no larger than those on *GHSR* (maximum 0.015), and
*GHSR* remained the only plausible target.

To establish what magnitude of effect is achievable at this locus, in silico
mutagenesis was used to score all 1,197 possible single-base substitutions
across the sequenced amplicon in four conditions. The amplicon contained
substitutions producing predicted effects 4.7- to 95-fold larger than the
largest effect of any studied variant, confirming that the model is responsive
at this locus. Across 48 condition-by-endpoint combinations, none of the three
variants exceeded the 90th percentile of its own local distribution, and
median percentiles were 43.7 (rs2922126), 42.1 (rs34136054) and 18.9
(rs115459563) (Figure 2B; Supplementary Table S5). The three variants are
therefore unremarkable members of the local distribution of mutational effects
rather than positions of unusual regulatory importance.

### Cell-type specificity

No variant produced a consistent direction of effect across conditions
(Figure 3). The alternate allele raised predicted *GHSR* expression in 15 of 16
conditions for rs2922126, in 15 of 16 for rs34136054 and in 5 of 16 for
rs115459563, with the remaining conditions showing the opposite sign. The
magnitude of the predicted effect correlated with the predicted *GHSR*
expression level of the condition (Spearman rho 0.63 for rs2922126, 0.78 for
rs115459563 and 0.83 for rs34136054), which is the expected consequence of a
multiplicative effect acting on a variable baseline. Predicted effects in the
three bone-related conditions were the smallest in the panel, with absolute
log2 fold changes below 0.0004 for all three variants.

Repeating all predictions at 16 kb, 131 kb, 524 kb and 1 Mb sequence contexts
showed that the sign of the predicted *GHSR* expression effect was preserved
across all four contexts in only 13 of 48 condition-by-variant combinations
(27%; 3 of 16 for rs2922126, 7 of 16 for rs115459563, 3 of 16 for rs34136054)
(Supplementary Tables S8 and S9). Effects of this magnitude therefore lie
within the model's own context-dependent variability and do not support
inference of a direction of effect. For rs34136054, the two alternate alleles
could differ in sign: in the pituitary gland the predicted DNase effect was
+0.024 for T>C and -0.033 for T>G.

### Analyses restricted to observed data

Because imputed tracks predominate at this locus, every summary was recomputed
using only tracks backed by observed experimental data (Supplementary
Table S6). This restriction is most consequential for transcription-factor
ChIP, where only 60 of the 43,728 track-level values for the three variants rest on
observed data.

Two observations survived this restriction and were more stable across
sequence contexts than the expression predictions, and are reported as
hypotheses for experimental testing rather than as established effects. First,
rs2922126 was predicted to reduce CTCF ChIP signal in all 12 observed CTCF
tracks available across eight conditions at the 1 Mb context (mean log2 fold
change -0.010; largest -0.025 in osteoblast), and this direction was preserved
across all four sequence contexts in 7 of 8 conditions (Supplementary
Table S10). The 501 bp scoring window centred on rs2922126 does not reach the
annotated CTCF element 349 bp away, no sequence-motif analysis was performed,
and no disruption of a CTCF binding site is claimed. Notably, this effect was
shared across cell types rather than restricted to any, which does not support
cell-type-specific regulation. Second, rs34136054 was predicted to increase
active-promoter histone marks in brain tissue on observed tracks (temporal
lobe H3K4me3 +0.032, H3K9ac +0.027, H3K27ac +0.027), with the same direction
in brain and hippocampus.

### Quality control

All reference alleles matched GRCh38 after scoring; every variant
configuration was scored in every condition with an identical scorer count;
and aggregated centre-mask accessibility scores agreed in sign and magnitude
with values recomputed independently from base-resolution reference and
alternate tracks in all 16 comparisons tested. No P-values, confidence
intervals or significance thresholds were derived from any model output.
