# Discussion insertions (drafted for *Scientific Reports*)

> Four insertions, each cautious and each keyed to an existing location in the
> manuscript. See `MANUSCRIPT_INTEGRATION_RECOMMENDATIONS.md` §3.3 for
> placement.

---

## §1. For "GHSR Polymorphisms and Cognitive-Group Classification"

> *Insert after the existing paragraph beginning "However, the present study
> evaluated genetic variants at the DNA level and did not measure GHSR
> expression, receptor isoform abundance, circulating ghrelin concentrations,
> or downstream signalling activity."*

To address this gap at least computationally, we retrospectively evaluated the
three genotyped variants with AlphaGenome 2, a sequence-based model that
predicts regulatory output directly from DNA. Across 16 cell types and
tissues, none of the three variants produced a predicted change in *GHSR*
expression of a magnitude that would plausibly be functional: the largest
effect anywhere in the panel was a 1.9% change in predicted signal
(|log2 fold change| = 0.027), and the corresponding figures for chromatin
accessibility and histone modification were similarly small. The
interpretation of such small values depends critically on knowing what a large
value would look like at this locus, so we generated a local reference
distribution by scoring every possible single-base substitution across the
same 399-bp sequenced amplicon. That analysis showed that the amplicon does
contain positions at which substitutions produce predicted effects 4.7- to
95-fold larger than anything the three studied variants produced, and that
across 48 condition-by-endpoint combinations none of the three variants
exceeded the 90th percentile of its local distribution. The model is therefore
demonstrably responsive at this locus, and the near-null predictions are
informative negatives rather than an absence of signal in the analysis itself.

These predictions also bear on the specific hypothesis that *GHSR* regulatory
variation might act in a cell-type-specific manner. We found no evidence for
such specificity. No variant moved predicted *GHSR* expression in a consistent
direction across conditions, and the magnitude of the predicted effect
correlated with how strongly *GHSR* itself was predicted to be expressed in
each condition (Spearman rho 0.63-0.83), which is the expected consequence of
applying a multiplicative effect to a variable baseline rather than evidence
of condition-specific regulatory logic. A further test reinforces this
reading: when the predictions were repeated across four sequence-context
lengths, the sign of the predicted expression effect was preserved in only 13
of 48 condition-by-variant combinations, indicating that effects of this size
lie within the model's own context-dependent variability and do not support
inference of a direction of effect.

The one observation that was both directionally consistent and stable ran
counter to the cell-type-specificity hypothesis rather than supporting it.
rs2922126 was predicted to reduce CTCF ChIP signal in all 12 observed CTCF
tracks available across eight conditions, and this direction was preserved
across all four sequence contexts in seven of those eight conditions. Because
the effect was shared across cell types rather than restricted to any, it is
evidence of a cell-type-*invariant* prediction, not a specific one. We
emphasise that the effect is small, that the 501 bp scoring window does not
reach the annotated CTCF element 349 bp away, that no sequence-motif analysis
was performed, and that we therefore make no claim that a CTCF binding site is
disrupted; it is reported solely as a hypothesis for experimental testing.

This computational result is concordant with, but logically independent of,
the absence of association reported above. The two analyses are null about
different things: the predictions concern cis-regulatory activity, whereas the
association analysis concerns phenotype in 100 participants, and neither
establishes the other. In particular, the low achieved power of the
association analysis for rs2922126 (6.62%) means that it could not have
detected a small genotype-phenotype effect regardless of what the model
predicts.

---

## §2. For "GHSR Polymorphisms, Bone Mineral Density, and Thyroid Function"

> *Insert after the existing discussion of ghrelin signalling in osteoblasts
> (refs 55-57) and the Amini/Lee contrast (refs 58, 59).*

The computational analysis adds a specific consideration to this literature.
Although ghrelin and its receptor have been reported in osteoblasts, and
ghrelin has been shown to influence osteoblast proliferation and
mineralisation (55-57), AlphaGenome 2 predicted *GHSR* to be essentially
unexpressed in bone-tissue osteoblasts, osteocytes and mesenchymal stem cells:
predicted *GHSR* activity in those cell types was below 1% of the panel
maximum, which was reached in the pituitary gland, with the hypothalamus
second. Reassuringly, that ranking recapitulates the established tissue
distribution of *GHSR* without having been specified in advance, which
supports the reliability of the model's output at this locus. Consistent with
the low predicted baseline, the predicted allelic effects in the three bone
conditions were the smallest in the entire panel (|log2 fold change|
< 0.0004 for all three variants).

If *GHSR* is expressed at very low levels in osteoblast-lineage cells, then a
cell-autonomous mechanism by which upstream *GHSR* regulatory variants could
alter bone cell behaviour has little room to operate, and any genuine effect
of the ghrelin-GHSR axis on bone would be more likely to act indirectly, for
example through growth hormone secretion from the pituitary, through
hypothalamic regulation of energy balance, or through GHSR-independent
pathways of the kind proposed by Delhanty and colleagues (56). This offers a
possible interpretation of why the present study, like the meta-analysis of
Lee et al. (59), found no consistent relationship between the ghrelin-GHSR
system and bone mineral density.

Two caveats apply. First, these are predicted and not measured expression
levels, and *GHSR* expression in human osteoblasts has been reported
experimentally; a model trained predominantly on bulk assays may underestimate
expression in cell types that are sparsely represented in its training data.
Second, no osteoclast or osteoclast-precursor condition exists in the model's
condition set, so the osteoclast arm of the hypothesis could not be assessed
directly, and blood monocyte was used only as a declared surrogate.

---

## §3. For "Strengths and Limitations"

> *Insert as an additional numbered limitation.*

Finally, the complementary in silico analysis has limitations of its own.
AlphaGenome 2 returns model predictions rather than measurements, and no
P-value, confidence interval or significance threshold can be attached to
them; the descriptive magnitude bands used here are labels for predicted
effect size and not statistical categories. The analysis was retrospective and
did not inform the experimental design. A sequence-based model cannot exclude
mechanisms it does not represent, including effects mediated by distal
elements outside the sequence window, by three-dimensional chromatin
interactions not captured by the contact-map predictions, by haplotype
combinations rather than single variants, by cell types or developmental
stages absent from the model's training data, or by processes operating at the
RNA or protein level. Coverage of the specific cell types of interest was also
incomplete: no osteoclast-lineage condition was available, and neither of the
two hippocampal conditions has observed RNA-seq data. At this locus the
transcription-factor ChIP predictions rest almost entirely on imputed tracks
(43,668 of the 43,728 track-level values for the three variants were in the
model's weakest evidence class), so the only transcription-factor observations we report are those that
persist when the analysis is restricted to tracks backed by observed data.
Absence of a predicted effect is therefore not evidence that these variants
are biologically inert, and the predictions are not clinical or diagnostic
evidence.

---

## §4. For "Future Research Directions"

> *Extend the existing paragraph, or add as a second paragraph.*

The computational analysis also indicates which functional experiments would
be informative, and which would not. Because the predicted effects on *GHSR*
transcription were uniformly very small, a conventional luciferase reporter
assay comparing the 399-bp upstream fragment carrying reference and alternate
alleles would most likely be underpowered to resolve differences of the
predicted size, and a negative result would be uninformative. Should such an
assay be attempted, it would need to be run in a cell type with appreciable
*GHSR* expression, such as a pituitary-derived or hypothalamic neuronal line
rather than an osteoblastic line, with enough replicates to resolve effects of
a few per cent, and with the three variants tested both individually and in
combination to capture possible haplotype effects.

Two more specific and better-powered experiments follow from the predictions
that survived restriction to observed data. First, the prediction that
rs2922126 reduces local CTCF occupancy could be tested directly by
allele-specific CTCF chromatin immunoprecipitation followed by quantitative
PCR, or by CUT&RUN, in cells heterozygous at this position; this would be a
direct test of a specific molecular prediction, and a motif analysis of the
surrounding sequence should precede it. Second, the prediction that
rs34136054 increases active-promoter histone marks (H3K4me3, H3K9ac, H3K27ac)
in brain tissue could be tested by allele-specific chromatin
immunoprecipitation for those marks, although the rarity of the variant makes
assembling a heterozygous sample set the limiting step.

More informative than either, for a variant of this kind, would be
allele-specific expression analysis in human tissue: quantifying the ratio of
*GHSR* transcripts derived from each allele in pituitary, hypothalamic or
cortical tissue from heterozygous donors internally controls for trans-acting
and environmental variation and can resolve small cis effects that a reporter
assay cannot. Public expression quantitative trait locus resources for brain
and pituitary tissue should be interrogated first, since the required data may
already exist. Finally, any attempt to establish a genotype-phenotype
relationship for these variants will require a substantially larger and
independent cohort with clinically adjudicated outcomes; the rarity of
rs115459563 and rs34136054 means that targeted sequencing of large cohorts, or
meta-analysis across cohorts, will be necessary before their effects can be
estimated at all.
