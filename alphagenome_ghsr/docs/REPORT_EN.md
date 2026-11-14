# AlphaGenome 2 functional analysis of three *GHSR* upstream variants

**A retrospective, complementary in silico investigation**

Execution date: 2026-11-14 (UTC) · Genome build: GRCh38/hg38 · Model:
AlphaGenome 2 via `alphagenome` 2.0.0 (`alphagenome.models.v2`)

---

## Summary

Three *GHSR* polymorphisms (rs2922126, rs115459563, rs34136054) previously
genotyped in 100 older adults were evaluated with AlphaGenome 2 to test
whether they are predicted to alter *GHSR* transcriptional regulation in a
cell-type-specific manner.

They are not. Across 16 cell types and tissues, 4 variant configurations and
162,144 track-level scores, the largest predicted change in *GHSR* expression
anywhere was a log2 fold change of 0.027 — a 1.9% change in predicted signal.
A local reference distribution built by scoring every possible single-base
substitution in the same 399 bp amplicon shows that the model predicts effects
4.7- to 95-fold larger for other positions in that same window, so the model
is demonstrably responsive here and the near-null result is an informative
negative rather than a failed analysis. None of the three variants reached the
90th percentile of its own local distribution in any of 48 condition-endpoint
combinations.

The central hypothesis of cell-type-specific regulation finds no support.
Effect magnitude simply tracks how strongly *GHSR* is predicted to be
expressed in each condition (Spearman rho 0.63–0.83); no variant has a
consistent direction of effect across conditions; and when sequence context is
varied, the sign of the predicted expression effect survives in only 13 of 48
condition-variant combinations. The single most stable finding — rs2922126
lowering predicted CTCF ChIP signal — is *shared across* cell types, which is
the opposite of cell-type specificity.

The bone arm of the hypothesis is the weakest. *GHSR* is predicted to be
essentially unexpressed in osteoblasts, osteocytes and mesenchymal stem cells
(≤1% of the panel maximum), and predicted allelic effects there were the
smallest observed.

This is concordant with, but logically independent of, the source study's null
association findings. Both are null; they are null about different things.

---

## 1. What was done

### 1.1 Design and timing

The genotyping, biochemical and densitometric work reported in the source
manuscript was completed before this analysis began. **The computational
analysis did not guide variant selection, primer design, or the statistical
plan.** It is a separate, retrospective layer of evidence, and is described as
such throughout.

No participant-level data were used. Only public GRCh38 variant coordinates
and public ontology identifiers were transmitted to any external service.

### 1.2 Pipeline

| Phase | What it does | Key outputs |
|---|---|---|
| I | Independent coordinate/allele verification; amplicon mapping; annotation | `table_S1`, `table_S2` |
| II | Condition selection against the model's own catalogue; evidence-class audit | `table_S3`, `table_S4` |
| III | Variant scoring: 4 configurations × 16 conditions × 22 scorers | 162,144 raw scores |
| III-b | Base-resolution REF/ALT tracks | `data/raw/phase3b/` |
| IV | REF/ALT integrity, ISM reference distribution, context sensitivity, observed-vs-imputed, track consistency | `table_S5`–`table_S8` |
| IV-b | Cross-context sign stability | `table_S9`, `table_S10` |
| V | Three-level evidence integration | `table_2`, `table_3` |
| VI | Figures | Figures 1–4 |

---

## 2. Phase I — where these variants actually are

Coordinates and alleles were retrieved independently from **NCBI dbSNP build
157** and **Ensembl (GRCh38)** and required to agree. They did. Every reference
allele matched the GRCh38 reference sequence.

| Variant | chr3 (GRCh38) | REF>ALT (+) | Upstream of *GHSR* TSS | gnomAD global AF | gnomAD Middle Eastern AF |
|---|---|---|---|---|---|
| rs2922126 | 172,449,471 | T>A | 1,015 bp | 0.280 | 0.276 |
| rs115459563 | 172,449,444 | G>A | 988 bp | 0.004 | 0.007 |
| rs34136054 | 172,449,372 | T>C / T>G | 916 bp | 0.005 (G) | 0.017 (G) |

**Independent confirmation of the experimental design.** The primer pair
published in Table 1 of the source manuscript maps onto GRCh38 at
chr3:172,449,261–172,449,659 — a product of exactly **399 bp**, reproducing the
reported size, and containing all three variants.

***GHSR* is on the minus strand** (chr3:172,442,551–172,448,456), so upstream
means *higher* coordinates. Ensembl VEP classifies all three variants as
`upstream_gene_variant` for all three *GHSR* transcripts. None lies in a *GHSR*
exon or intron, and none overlaps an Ensembl Regulatory Build element; the
nearest are an annotated **enhancer** (ENSR3_B3983C, 90 bp from rs2922126) and
a **CTCF binding site** (ENSR3_53HK6KS, 349 bp away).

The manuscript's description of this region as "promoter/regulatory" is
therefore approximately right but imprecise. It is the *GHSR* **proximal
upstream region**, adjacent to — not within — an annotated enhancer.

### 2.1 Discrepancies with the manuscript (documented, not silently fixed)

1. **rs34136054 "C>T" inverts REF and ALT.** GRCh38 carries **T** at this
   position. dbSNP lists the site as tri-allelic (T>C, T>G), and gnomAD shows
   that the allele segregating in populations is **G**, not C. Genotype
   *identity* is unaffected, so no experimental result changes — but the
   notation is wrong, and the manuscript does not identify which alternate
   allele was seen. **Both were therefore scored.**
2. **Strand convention is inconsistent.** rs2922126 ("A>T") is written in
   *GHSR* minus-strand orientation — consistent with the sequencing, since the
   manuscript's "forward" primer anneals to the plus strand and so reads in
   the gene's direction. rs115459563 ("G>A") is written on the plus strand.
   Two conventions, three variants.
3. **rs115459563 is ~5× more frequent in the sample than expected.** Observed
   carrier frequency 7/50 controls (allele frequency 0.035) versus a gnomAD
   Middle Eastern estimate of 0.0068. With n = 50 this may be sampling
   variation or population structure, but since this is the manuscript's only
   nominally significant genotype finding (p = 0.012), **re-reading those seven
   chromatograms before publication is the single most consequential check we
   can recommend.** This is a wet-lab question; no computation can settle it.

---

## 3. Phase II — what the model can and cannot represent

AlphaGenome 2 is condition-conditioned: predictions are requested one
condition at a time, and each output type is labelled `observed`,
`imputed_supported`, or `imputed_exploratory`. That label is carried through
every table here.

Sixteen of the 2,799 available human conditions were selected to span the
phenotypes the source study measured. Two pre-specified priorities could not
be met, and this is recorded rather than worked around:

* **No osteoclast or osteoclast-precursor condition exists.** Blood monocyte
  was used as a declared surrogate (osteoclasts differentiate from circulating
  monocytes) and is flagged as a surrogate wherever it appears.
* **No hippocampal *neuron* condition exists.** The two available hippocampal
  conditions are a tissue layer and an astrocyte; neither has observed
  RNA-seq.

### 3.1 A useful sanity check, and a finding

Predicted *GHSR* expression across the panel ranks:

**pituitary gland (100%) > hypothalamus (41%) > testis (24%) ≈ temporal lobe
(23%) ≈ hippocampus (22%) > brain neuron (12%) ≈ brain (10%) ≫ glutamatergic
neuron, bone and peripheral tissues (≤1%)**

This reproduces the known tissue distribution of *GHSR* without being told to,
which supports the reliability of the model's output at this locus.

It is also a finding in its own right: ***GHSR* is predicted to be essentially
unexpressed in osteoblasts, osteocytes and mesenchymal stem cells.**

---

## 4. Phase III — the predictions

All four variant configurations were scored in all 16 conditions with the
complete recommended scorer set — not a subset — so that no modality could be
omitted or selected after the results were seen. All requests succeeded.

**Four scorers returned nothing, and that is the correct answer.** Splice
sites, splice-site usage, splice junctions and polyadenylation returned
well-formed *empty* results (shape `(0, n_tracks)`): the requests succeeded and
no such feature was affected. That is expected for a variant ~1 kb upstream of
a TSS and far from any splice junction. It is recorded as an empty result, not
as a failure.

### 4.1 Magnitudes

| Endpoint | rs2922126 | rs115459563 | rs34136054 |
|---|---|---|---|
| *GHSR* expression (max \|log2FC\|) | **0.027** | 0.016 | 0.025 |
| Chromatin accessibility | 0.041 | 0.043 | **0.051** |
| TF ChIP (aggregate) | 0.039 | 0.029 | 0.050 |
| Histone ChIP | 0.031 | 0.013 | 0.032 |
| Nearby genes | 0.013 | 0.015 | 0.007 |

These are **maxima over 16 conditions and every track**, i.e. the most extreme
value found anywhere — not typical effects. The largest, 0.051, is a 3.6%
change in predicted signal. Nothing approaches 0.1 log2FC (~7%).

---

## 5. Phase IV — is "very small" actually small?

A score of 0.027 means nothing in isolation. Two controls were applied.

### 5.1 Local reference distribution (ISM)

Every possible single-base substitution across the whole 399 bp amplicon
(1,197 variants) was scored in four conditions. Because ISM enumerates all
three substitutions at every position, **the three studied SNPs are themselves
members of this set**, so their scores come from the same run under identical
context, scorer and condition.

| | rs2922126 | rs115459563 | rs34136054 |
|---|---|---|---|
| Median percentile of \|effect\| within the amplicon | 43.7 | 18.9 | 42.1 |
| Cells above the 90th percentile (of 16 each) | 0 | 0 | 0 |
| Median ratio to the local maximum effect | 0.045 | 0.021 | 0.044 |

**Positive control:** the amplicon contains substitutions whose predicted
effects are **4.7- to 95-fold larger** than anything the three studied variants
produce. The model is responsive at this locus. The three variants are
unremarkable members of the local distribution of mutational effects.

### 5.2 Sequence-context stability

Repeating every prediction at 16 kb, 131 kb, 524 kb and 1 Mb:

| Readout | Sign preserved across all four contexts |
|---|---|
| *GHSR* expression | **13 / 48 (27%)** — rs2922126 3/16, rs115459563 7/16, rs34136054 3/16 |
| Observed CTCF ChIP | **19 / 24 (79%)** — rs2922126 7/8, rs115459563 7/8, rs34136054 5/8 |

This is the sharpest internal test available. A genuine local cis-regulatory
effect should not reverse when distant flanking sequence is added or removed.
**For expression, most do.** Effects of this magnitude lie within the model's
own context-dependent variability and carry no interpretable direction.

For rs34136054 the two alternate alleles can differ in sign (pituitary DNase:
+0.024 for T>C, −0.033 for T>G), which is a further reason the manuscript's
allele ambiguity matters.

### 5.3 Dependence on imputed tracks

At this locus, **43,668 of the 43,728 TF ChIP values for the three variants are
`imputed_exploratory`** —
the model's weakest evidence class. Only 60 rest on observed data. Every
summary was therefore recomputed on observed tracks only. Any TF statement
resting on the other 43,668 values would be resting on imputation.

### 5.4 Other QC

* All reference alleles re-verified against GRCh38 after scoring: pass.
* Every variant scored in every condition, identical scorer count: pass.
* Aggregated centre-mask scores vs values recomputed independently from
  base-resolution tracks: consistent in sign and magnitude, 16/16.

---

## 6. Phase VII — answering the scientific questions

### 6.1 Which variants, if any, show a potentially meaningful predicted effect?

**None, at the level of *GHSR* expression.** All three fall in the "very small"
band (0.01 ≤ |log2FC| < 0.1), below any magnitude a reporter assay could
resolve, unremarkable within their local ISM distribution, and mostly
sign-unstable across sequence contexts.

Two observations are worth carrying forward as **hypotheses**, because they
survive restriction to observed tracks *and* are context-stable:

1. **rs2922126 → reduced predicted CTCF occupancy.** Negative in all 12
   observed CTCF tracks across 8 conditions at 1 Mb (mean −0.010, largest
   −0.025 in osteoblast), with the sign preserved across all four contexts in
   7 of 8 conditions. The most directionally consistent result in the study.
   **Caveats:** small magnitude; the 501 bp scoring window does not reach the
   annotated CTCF element 349 bp away; no motif analysis was performed. **We do
   not claim a CTCF binding site is disrupted.**
2. **rs34136054 → increased active-promoter marks in brain.** Temporal lobe
   H3K4me3 +0.032, H3K9ac +0.027, H3K27ac +0.027 on observed tracks, same
   direction in brain and hippocampus.

### 6.2 In which cell types?

For expression, the largest values appear in the conditions where *GHSR* is
most expressed (pituitary, hypothalamus, temporal lobe, glutamatergic neuron),
but §6.3 explains why that is not meaningful specificity. The CTCF
observation is *not* cell-type restricted — it appears across brain, bone,
gonadal and gastrointestinal conditions alike.

### 6.3 Is the cell-type-specificity hypothesis supported?

**No**, on three independent grounds:

1. **No consistent direction.** rs2922126 raises predicted expression in 15/16
   conditions and lowers it in 1; rs115459563 raises in 5 and lowers in 11;
   rs34136054 raises in 15 and lowers in 1. None is uniform.
2. **Magnitude tracks baseline expression** (Spearman rho 0.63–0.83). A
   multiplicative effect on a near-zero baseline is near-zero; this is
   arithmetic, not regulatory logic.
3. **The one stable signal is shared, not specific.** The CTCF effect appears
   across unrelated cell types.

### 6.4 Is *GHSR* the likely affected gene?

**Yes, if anything is.** *GHSR* is overwhelmingly the nearest gene; the next
protein-coding TSS (*TNFSF10*) is 74 kb away. Predicted effects on
neighbouring genes were no larger than on *GHSR* (max 0.015, and the top hit
was *TBPL1P1*, a processed pseudogene 72 kb away — not a credible target).
There is no evidence that a neighbouring gene is the real target; there is
simply little predicted effect on anything.

### 6.5 Are the results compatible with the experimental observations?

**Yes, but the concordance is weaker than it looks.** Both evidence streams
are null, but about different things: Level 1 concerns cis-regulatory
activity, Level 3 concerns phenotype in 100 participants. Neither establishes
the other. The association analysis had **6.62% achieved power** for
rs2922126, so it could not have detected a small effect regardless of what the
model predicts.

The *GHSR*-not-expressed-in-bone prediction does offer a mechanistic reading
of why the bone associations were null — but it is a prediction, and *GHSR*
has been reported experimentally in human osteoblasts, so it should be stated
as a possibility, not a conclusion.

### 6.6 Which findings are weak, inconclusive or contradictory?

* **Weak:** all *GHSR* expression effects — small, mostly sign-unstable.
* **Weak:** the entire TF ChIP layer except CTCF and a handful of others, being
  almost entirely imputed-exploratory.
* **Inconclusive:** rs34136054 generally — the manuscript does not identify
  which alternate allele was observed, and the two alternates can differ in
  sign.
* **Inconclusive:** rs115459563 — the nominal genotype distribution cannot be
  assessed computationally, and its elevated frequency in the sample needs
  laboratory confirmation.
* **Contradictory with the hypothesis:** the only stable signal (CTCF) is
  cell-type-invariant, which argues against cell-type-specific regulation.
* **Not contradictory, but important:** *GHSR* predicted unexpressed in
  osteoblasts conflicts with published experimental reports of *GHSR* in human
  osteoblasts. The model may be wrong here; bulk-assay training data
  under-represent such cell types.

### 6.7 What experiments would actually test this?

**Not a standard luciferase assay first.** Predicted effects are far below
what such an assay resolves, so a negative result would be uninformative. If
attempted, it must use a *GHSR*-expressing line (pituitary-derived or
hypothalamic neuronal, **not** osteoblastic), enough replicates to resolve a
few per cent, and test the variants both individually and as haplotypes.

In priority order:

1. **Allele-specific expression in human tissue** (pituitary, hypothalamus,
   cortex from heterozygous donors). Internally controls for trans and
   environmental variation; resolves small *cis* effects a reporter cannot.
   **Check public brain/pituitary eQTL resources first — the data may already
   exist.**
2. **Allele-specific CTCF ChIP-qPCR or CUT&RUN** in cells heterozygous at
   rs2922126 — a direct test of the one stable prediction. Run a sequence-motif
   analysis first.
3. **Allele-specific ChIP for H3K4me3/H3K9ac/H3K27ac** in brain tissue for
   rs34136054; the rarity of the variant makes sample assembly the limiting
   step.
4. **Re-read the seven rs115459563 chromatograms** (§2.1.3). Cheapest, and
   potentially the most consequential.
5. **Independent, adequately powered cohorts** with clinically adjudicated
   outcomes for any genotype–phenotype claim.

---

## 7. Limitations

**Of the computational analysis**

* Predictions are not measurements. No P-value, confidence interval or
  significance threshold exists for any score here, and none was invented. The
  magnitude bands are descriptive labels.
* A sequence model cannot exclude mechanisms it does not represent: distal
  elements outside the window, 3D chromatin effects beyond the contact-map
  predictions, haplotype effects, absent cell types or developmental stages,
  and anything acting at the RNA or protein level.
* Condition coverage is incomplete (no osteoclast lineage; no hippocampal
  neuron; neither hippocampal condition has observed RNA-seq).
* The TF ChIP layer is almost entirely imputed at this locus.
* The served model checkpoint cannot be pinned from the client; the SDK
  version, API host, parameters and execution date are the reproducibility
  anchor.
* The canonical literature citation for AlphaGenome 2 could not be verified
  from this environment (the documentation host was not reachable) and has
  deliberately not been reproduced from memory — see `docs/CITATIONS.md`.

**Inherited from the source study**

* n = 100, underpowered for its genetic questions (6.62% achieved power for
  rs2922126); rare genotype categories; non-convergent models for rs115459563;
  rs34136054 seen once.
* The outcome is a cognitive **screening** score (6-CIT), not a clinical
  Alzheimer's disease diagnosis. **Nothing here is evidence about diagnosed
  Alzheimer's disease, and none of it is clinical or diagnostic evidence.**

---

## 8. Conclusion

Within what a sequence model can say: these three *GHSR* upstream variants are
not predicted to have cell-type-specific regulatory consequences of a
magnitude that would plausibly matter. The result is a well-controlled
negative — the model is demonstrably responsive at this locus, other
substitutions in the same amplicon produce far larger predicted effects, and
the studied variants sit near the middle of that distribution.

The one stable, observed-data signal (rs2922126 lowering predicted CTCF
occupancy) is cell-type-*invariant*, which argues against the specific
hypothesis rather than for it, and is offered only as a wet-lab hypothesis.

A negative result is a result. What it does not license is the inverse claim:
this analysis does **not** show that these variants are biologically inert —
only that this model, at this locus, in these conditions, predicts little.
