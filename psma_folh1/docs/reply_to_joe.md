# Response to review — PSMA/FOLH1 manuscript

Dear Joe,

Thank you — and please thank Uta too. You were right on all three counts, and on two of them the
data turned out to be worse than your review implied. We took the manuscript apart, had the whole
analysis independently replicated from the raw data by someone who had not written the original
code, and then rebuilt it along the five-step structure you proposed. The short version: the
computational work was sound and fully reproducible, but the inferences were not supported, and
the paper is now a different paper.

Below: what we checked, what we found, and what the manuscript now says. The restructured draft is
attached as `manuscript_v2`, and everything is reproducible from the code and tables we have
included.

---

## 0. First, we checked whether the original results were even reproducible

Before changing anything we re-downloaded every cohort from its primary source and re-implemented
every analysis from the Methods section, without using any file from the original package. Of 78
numerical claims, **75 reproduce exactly** — usually to three decimal places, including the
external R²=0.662, the 1,002-gene differential-expression set, the Beltran AUC of 0.82, the
survival nulls, and the DepMap ranks 1/2/4/15 of 17,293. All 31 reference DOIs resolve correctly.

Two corrections came out of it, which we have made regardless of everything else:

- The Hallmark GSEA p-values were not reproducible as p<1e-3. With 10,000 permutations the correct
  values are p≈0.01 (FDR≈0.02). The enrichment scores themselves replicate exactly.
- The TCGA "primary tumours, n=550" cohort included 52 normal-prostate samples and one metastasis.
  The tumour-only FOLH1–AR correlation is 0.262, not 0.332.

So this is not a reproducibility problem. It is an inference problem, which is what you said.

---

## 1. "RNA-derived TF activities to predict RNA expression… does not establish independent regulatory mechanisms or clinically useful predictive performance"

Correct, and worse than that: the model was partly reading FOLH1 off itself.

**FOLH1 is a member of the CollecTRI regulons of ELF5, ELF3, AR and VDR.** ELF5's regulon contains
11 genes, one of which is FOLH1, so FOLH1 contributes about 9% of ELF5's inferred activity
directly — and ELF5 was the model's dominant positive weight (+0.655, next-largest coefficient
−0.118) and the strongest single TF–FOLH1 association in mCRPC (ρ=+0.67).

Removing FOLH1 from every regulon and rerunning the identical pipeline:

| | published | FOLH1 removed from regulons |
|---|---|---|
| external R² | 0.662 | **0.561** |
| ELF5 model weight | +0.655 | **+0.040** |
| ELF5–FOLH1 ρ (TCGA) | +0.267 | **−0.043** (sign reverses) |

We then asked what the remaining 0.561 means. Twelve sets of **size-matched random gene sets** used
as pseudo-regulons, run through the identical pipeline, reach a mean external R² of 0.466 (max
0.574). The expression of **766 randomly chosen genes** reaches 0.509. In other words, once the
leak is removed, the model performs about as well as knowing nothing about transcriptional
regulation at all.

On clinical utility you were also right, and we now say so explicitly: the model predicts FOLH1
RNA from transcription-factor activities inferred *from the same RNA sample*. Anyone who can
compute the predictors has already measured the outcome. Its 95% prediction interval for an
individual tumour spans about 72-fold in FPKM, and absolute calibration does not transport across
platforms (raw-space R² = −12).

**What we did.** The transcription-factor model is removed from the manuscript. If a reviewer asks
for it, we will supply the de-leaked version as a supplementary, descriptive observation about
regulatory co-variation — not as a predictor, and not as evidence of regulatory architecture.

---

## 2. "Tumor composition, protein expression, and lesion-level PSMA targetability are not adequately addressed. The external validation largely confirms known lineage differences"

Correct on every clause. Three things we found:

**(a) The external validation was weaker than a standard marker panel.** In Beltran, a routine
12-gene neuroendocrine score separates NEPC from CRPC-adenocarcinoma with AUC 0.929; our 100-gene
"externally validated PSMA-low signature" managed 0.820, and FOLH1 alone managed 0.782. The
signature was a worse NE classifier than the markers it was implicitly competing with.

**(b) The GSE126078 validation was close to uninformative as posed.** A signature derived from the
quartiles of *any* gene anti-correlates with that gene in another cohort. Repeating our whole
pipeline for 15 random anchor genes gives |ρ| between 0.38 and 0.85 (median 0.65) in GSE126078; our
published 0.824 sits inside that null. (The Beltran AUC did exceed the random-anchor null, so the
signature does carry real lineage information — it simply adds nothing to markers that already
exist.)

**(c) Tumour composition was doing a great deal of the work.** We rebuilt the analysis with FOLH1
as a continuous variable and explicit adjustment for biopsy site, immune/stromal/hepatocyte/
osteoblast content (PanglaoDB compartments) and neuroendocrine score. Of 6,438 genes associated
with FOLH1 in crude analysis, **2,140 are explained by composition or site and 3,232 more by
neuroendocrine differentiation; 1,557 remain.** Of those, 251 reproduce in GSE126078 at sample
level and 191 in all three cohorts.

The programme-level consequence is specific and, we think, genuinely interesting:
epithelial-mesenchymal transition **reverses sign** on adjustment (NES +1.7 → −1.2), as do
inflammatory response (+1.4 → −1.3) and hypoxia (+1.5 → −1.0). The EMT-high, inflamed PSMA-low
phenotype we described was a property of the biopsy, not of the tumour cells. What survives is
narrow: proliferation up, androgen response and oxidative phosphorylation down.

**On your clinical-construct point.** Using the published molecular-phenotype calls in GSE126078,
PSMA-low is overwhelmingly AR-low or AR-negative disease: **only 6 of 60 AR-positive,
non-neuroendocrine metastases fall in the cohort's bottom FOLH1 quartile**, and within that
subgroup only 15 genes are associated with FOLH1. We no longer claim a distinct AR-intact PSMA-low
biology — the data do not support one at this sample size, and saying so is more useful than the
alternative.

**On the quartile cut-off.** We took your advice and moved to a continuous model throughout. For
what it is worth, this was the one methodological criticism that did not bite: quartile, tercile,
median and decile contrasts all agree with the continuous model at ρ 0.91–0.94. We report the
continuous model anyway.

---

## 3. "The proposed alternative targets and TXNRD1 dependencies are based on computational associations without sufficient protein-level or functional validation"

Correct, and here the data point the other way from our conclusions.

**On the targets.** We rebuilt the surfaceome analysis around therapeutic suitability rather than
differential expression, with four pre-specified criteria: abundance in PSMA-low tumours,
tumour-to-normal contrast benchmarked against PSMA's own, tumour-cell specificity from HPA
single-cell data, and protein-level immunohistochemistry. The result:

- **uPAR (PLAUR) does not survive at all.** It loses significance once composition is adjusted for
  (FDR 0.38), draws 40% of its single-cell expression from immune and stromal cells (top cell
  types: neutrophils and monocytes), sits 8 log2 units below PSMA on tumour-to-normal contrast,
  and is undetectable by immunohistochemistry in 9 of 9 prostate cancers. It was a myeloid
  infiltration signal.
- **LAT1 and GPC2 do not reproduce in any independent cohort**; GPC2 is above threshold in fewer
  than a third of PSMA-low tumours. **MCT4** reproduces in one cohort but is 3.6 log2 units below
  PSMA on contrast with 21% non-tumour-cell expression.
- **The luminal antigens are reduced, not absent.** TROP2 (median 5.45, 87% of PSMA-low tumours
  positive), B7-H3 (4.46, 93%) and STEAP1 (4.14, 78%) are each more abundant in PSMA-low tumours
  than any of our nominated candidates — and more abundant than FOLH1 itself there (1.96, 49%).
  "Significantly down-regulated" was not the same as "unavailable", and we drew the wrong
  conclusion from a rank-biserial correlation.
- **B7-H3 (CD276) is the only antigen that meets all four criteria**: retained in 92.5% of PSMA-low
  tumours, tumour-to-normal contrast within 1.2 log2 units of PSMA's own, 14% non-tumour-cell
  expression, positive by IHC in all 32 prostate cancers tested, every one at medium-to-high intensity with an
  "Enhanced"-reliability antibody — and with clinical agents already in trials.

So the target-switching conclusion is inverted: the paper now argues that alternatives should be
sought among antigens that are **retained**, not among antigens that are **enriched**, because
enrichment in a compartment defined by loss of a luminal marker selects for whatever is high in
non-luminal cells — which in a bulk biopsy means stroma and infiltrate.

**On TXNRD1.** This one does not survive either, and we are dropping it. Three independent checks:

1. Lineage alone explains 30–38% of the signature score's variance across DepMap; prostate is the
   lowest-scoring of 23 lineages. That a haematopoietic-restricted gene (NCKAP1L) appeared among
   the top hits should have warned us.
2. Demeaning within lineage, the "axis" collapses to TXNRD1 alone (SEPSECS falls from rank 1 to 9,
   PSTK from 15 to 265, GPX4 from 13 to 602), and GAPDH — a housekeeping benchmark — outranks three
   of its four members.
3. Rescoring DepMap with the **cross-cohort-validated** programme instead of the discovery-only
   signature, the axis disappears entirely: TXNRD1 falls to rank 1,709 of 17,293, SEPSECS to 6,267.
   What appears instead is nucleotide biosynthesis and p53 status — proliferation, not a selective
   vulnerability.

And the decisive test: correlated **directly against FOLH1 expression**, TXNRD1 dependency gives
ρ = −0.09, i.e. the *wrong sign* — FOLH1-high lines are marginally the more dependent ones. In
prostate lines (n=10), ρ=+0.16, p=0.66.

The therapeutic dependency section is removed.

---

## 4. The AlphaGenome section

You did not raise this, but we checked it while we were at it. Everything in Figure 5 reproduces
exactly — the strongest predicted prostate DNase element (chr11:49,177,981–49,178,322), the
variant effects, the 86-position tiling scan (−0.23 at the core), the absence of enhancer–promoter
contact. Two things changed:

- **The "CpG disruption" result is not about CpGs.** At the interrogated cytosine, C>T gives
  +0.288 but C>G gives +0.086 and C>A +0.006 — all three abolish the CpG. AlphaGenome takes DNA
  sequence as input and has no methylation channel, so the "directional concordance with the
  negative methylation–expression correlation" was a coincidence. The claim is restated without
  the methylation framing.
- **In our favour, for once:** the nominated enhancer survives a control the original analysis
  lacked. The same 4-bp scan over a matched peak-free gene-body region gives SD 0.007 and maximum
  effect 0.029, versus 0.053 and 0.233 for the element (7.1-fold dispersion, Levene p=2.3e-8).
  That control is now in the paper.

The section is reduced to a short, clearly exploratory subsection.

---

## 5. What the paper now says

Following your five-step outline:

1. **Discovery (SU2C, n=266)** — FOLH1 continuous, adjusted for site, composition and NE. Three
   quarters of the crude programme is confounding; what remains is proliferative and AR-axis-low.
2. **Validation (GSE126078 patient-level; Beltran for the NE question)** — 251 genes reproduce,
   191 in all three cohorts; programme-level agreement ρ=0.61. PSMA-low is, in the main, AR-low
   disease.
3. **Surfaceome prioritisation** — the centrepiece, as you suggested. Retained-versus-lost rather
   than enriched; four pre-specified suitability criteria; established and nominated targets scored
   side by side. B7-H3 is the only antigen that passes.
4. **DepMap** — exploratory, and now negative: no selective vulnerability of the PSMA-low state.
5. **TCGA** — secondary context: the programme is already present in treatment-naive primary
   tumours (98% of programme genes keep their direction), so it is a property of low-PSMA prostate
   cancer rather than of treatment resistance.

The title is now *"What distinguishes PSMA-low metastatic prostate cancer, and which surface
targets remain available in it"*.

We think the paper is smaller and considerably more useful than it was. Two of its three headline
claims are now negative results, which is the right outcome given what the data show, and the
positive claim — retention rather than enrichment as the basis for alternative-target nomination,
with B7-H3 as the leading candidate — is one we can actually defend.

**On your suggestion B:** yes, please. A fresh reader from Uta's class would be very welcome, and
we would specifically ask them to attack the surfaceome prioritisation, since that is now where
the weight of the paper sits. The criteria and thresholds are pre-specified and the code is
available, so it should be an easy thing to break if it is breakable.

Two things we would value your view on:

- Whether the B7-H3 argument is strong enough to carry a paper, given that it rests on retention
  and contrast rather than on enrichment, and that the decisive evidence would be lesion-level
  IHC in PSMA-PET-negative metastases, which we do not have.
- Whether the AR-low framing is the right one, or whether we should go further and simply say that
  "PSMA-low" in bulk transcriptomes is largely a restatement of AR-indifferent disease, and that
  the clinically distinct entity — AR-intact, PSMA-PET-negative — cannot be studied without
  imaging-anchored sampling. We lean towards saying that plainly.

With thanks for a review that was more useful than a favourable one would have been,

Mehdi
