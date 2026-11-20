# Phase 1 — Independent replication and validation

**Manuscript under review:** *An externally validated transcription-factor model of PSMA (FOLH1)
expression defines the PSMA-low state and nominates target-switching and druggable-dependency
strategies in metastatic prostate cancer* (deliverable package v1.1).

**What this is.** A from-scratch replication. I did not run the authors' notebook and did not use
any file from the deliverable package as input. Every cohort was re-downloaded from its primary
source, and every analysis was re-implemented from the Methods section. The authors' notebook was
read only afterwards, to explain discrepancies and to audit choices that the Methods do not state.

**Two separate questions are answered below, and they have different answers:**

1. *Do the numbers reproduce?* — **Yes, essentially perfectly.** 69 numerical claims were
   catalogued. One is not testable here (AlphaGenome, no API key). Of the other 68, **65 reproduce
   exactly**, most to three decimal places; one reproduces in direction with a different p-value
   (Gleason); one is technically correct but uninformative as stated (ChIP-Atlas); and **one does
   not reproduce** (GSEA permutation p-values).
2. *Do the numbers support the conclusions?* — **For the three headline conclusions, no.** A set of
   validation tests the manuscript does not contain changes the interpretation of the TF model, the
   external validation and the therapeutic claims. These are new analyses, not failures to
   reproduce.

This distinction matters for the reply to Joe: the work is **computationally sound and
reproducible**; the problem is **inferential**, which is exactly what he says.

---

## 1. Replication set-up

| Dataset | Source used here | Matches the manuscript? |
|---|---|---|
| TCGA-PRAD | UCSC Xena `HiSeqV2` + `PRAD_clinicalMatrix` + `HumanMethylation450` + 450k probemap | yes (20,530 × 550) |
| SU2C/PCF mCRPC | cBioPortal datahub `prad_su2c_2019`, polyA FPKM | yes (19,270 × 266) |
| Beltran CRPC/NEPC | cBioPortal datahub `nepc_wcm_2016`, FPKM | yes (33,194 × 49; 34 adeno / 15 NEPC) |
| GSE126078 | GEO per-sample processed FPKM (153 files) | yes (98 patient / 43 LuCaP / 12 cell line) |
| DepMap 24Q4 | figshare `CRISPRGeneEffect`, `OmicsExpressionProteinCodingGenesTPMLogp1`, `Model` | yes (1,103 models with both) |
| HPA | `rna_tissue_consensus.tsv` (51 tissues) + `rna_single_cell_type.tsv` | yes (51 tissues) |
| Surfaceome | SURFY + CSPA high-confidence + UniProt cell-membrane w/ TM or GPI | 3,702 vs 3,701 genes |
| CollecTRI | decoupleR / OmniPath | yes (42,990 edges, 1,185 TFs) |

Code: `psma_folh1/code/` (`s01`–`s08`). Machine-readable outputs: `psma_folh1/results/`.
Claim-by-claim table: `psma_folh1/results/claim_by_claim_verdicts.csv`.

---

## 2. Replication verdict — the numbers hold

Representative examples (full table in the CSV; "manuscript → replication"):

- FOLH1–AR score: TCGA 0.332 → **0.332**; SU2C 0.546 → **0.546**; Beltran 0.752 → **0.752**
- TF–FOLH1 cross-cohort concordance 0.502 → **0.502**; 367/766, 594/766, 39 opposite → **all exact**
- AR regulon decoupling: −0.052 → +0.259, Fisher z = −4.22 → **identical, including z and p**
- Model A nested-CV R² 0.717 → **0.719**; TF-only 0.712 → **0.714**
- Model B external R² 0.662 (CI 0.559–0.736), ρ = 0.802 → **0.662 (CI 0.563–0.735), ρ = 0.802**
- AR-only baseline 0.382 → **0.382**; ELF5 weight +0.655 → **+0.655**
- PSMA-low DE: 1,002 genes (168 up / 834 down) → **1,002 (168 / 834)**; concordance 0.955 → **0.955**
- Beltran AUC 0.82 → **0.820**; FOLH1-free 0.814 → **0.814**
- GSE126078 ρ = −0.757 / −0.824 / −0.753 → **−0.757 / −0.824 / −0.753**
- Survival: log-rank 0.89 → **0.893**; HR 1.03 (0.77–1.38) → **identical**; HR 1.06 (0.76–1.47) → **identical**
- Surfaceome: 78 surface genes up → **78**; STEAP1 −0.65 → **−0.663**; DLL3 null → **null**; SSTR2 FDR 0.055 → **0.055**
- DepMap: TXNRD1 −0.227 → **−0.227**; SEPSECS −0.234 → **−0.234**; ranks 1, 2, 4, 15 of 17,293 → **1, 2, 4, 15 of 17,293**
- FOLH1 GRCh38 coordinates/strand and the three CpG positions → **verified against Ensembl**

This is an unusually clean replication. The pipeline, the audit corrections applied in v1.1
(coordinate liftover, bootstrap pairing, common-essential exclusion, FOLH1-free re-validation),
and the supplementary tables are all internally consistent with the data.

### 2.1 The one claim that does not reproduce

**Hallmark GSEA p-values.** The NES values reproduce exactly (E2F 3.60, G2M 3.12, EMT 2.99,
androgen response −2.17). The p-values do not. The Methods state 1,000 permutations with "p = 0
reported as p < 1e-3". In replication the up-regulated sets give **p ≈ 0.08** at 1,000 permutations
and **p ≈ 0.010 (BH-FDR 0.021)** at 10,000. The sets remain significant, but not at p < 1e-3: for
E2F/G2M/EMT the correct p-values are two orders of magnitude larger than reported, though the
enrichment itself (NES) is unchanged. Fix: report the 10,000-permutation p-values and FDR.

### 2.2 Reference list

All 31 DOIs in the reference list resolve through Crossref to papers whose titles, journals and
years match the citations as written. Given the disclosed use of generative AI in drafting, this is
worth stating in the reply — there are no fabricated or mismatched references.

### 2.3 Claims not testable here

**AlphaGenome (Figure 5, Results §6).** Requires a Google AlphaGenome API key. The coordinate work
*is* verifiable and is correct: FOLH1 GRCh38 chr11:49,145,092–49,208,641 minus strand confirmed
against Ensembl; all three scored CpG-disruption positions are genuine CpG cytosines on the plus
strand; the nominated element lies inside the gene body 30.3 kb from the TSS.

> **Security issue — act on this today.** A live `ALPHAGENOME_API_KEY` is hard-coded in
> `execution_trace/psma_folh1_analysis.ipynb`, cell 41, and that notebook has been sent to
> external collaborators. **Revoke and rotate the key.** Replace it with
> `os.environ["ALPHAGENOME_API_KEY"]` before the package is circulated again.

---

## 3. Independent validation — the tests the manuscript does not contain

### 3.1 The TF model reads FOLH1 off itself (Joe's point 1, in a stronger form)

`FOLH1` is a **member of the CollecTRI regulons** of ELF3, **ELF5**, AR and VDR. ELF5's regulon has
only 11 genes — `FOLH1, KLK3, LYZ, STAT5A, CCND2, KRT77, MUC4, ITGB3, FOXA1, EGFR, LARGE1` — so
FOLH1 contributes ~9% of ELF5's inferred activity directly. ELF5 is the model's dominant positive
weight (+0.655, next-largest coefficient −0.118) and the strongest single TF–FOLH1 association in
SU2C (ρ = +0.666).

Removing FOLH1 from every regulon and re-running the identical pipeline:

| | with FOLH1 in the regulons | FOLH1 removed |
|---|---|---|
| Model B external R² | 0.662 | **0.561** |
| Model B external ρ | 0.802 | **0.727** |
| ELF5 model weight | +0.655 | **+0.040** |
| ELF5–FOLH1 ρ (SU2C) | +0.666 | **+0.350** |
| ELF5–FOLH1 ρ (TCGA) | +0.267 | **−0.043** (sign reverses) |

The Discussion's claim that "the dominant positive weight of ELF5, a luminal-lineage factor,
further anchors PSMA within the luminal program" is therefore largely an artefact of the outcome
variable sitting inside the predictor's gene set. This is cheap to fix — rerun with FOLH1 excluded
from all regulons and report that version — but it must be fixed.

**Is the remaining 0.561 evidence of transcription-factor architecture?** I built 12 sets of
pseudo-regulons: for each of the 770 real CollecTRI regulons, a random gene set of the same size
drawn from the transcriptome (FOLH1 excluded), then ran the identical ULM → elastic-net →
TCGA-to-SU2C pipeline.

| | external R² |
|---|---|
| real CollecTRI TF model, as published | 0.662 |
| real CollecTRI TF model, FOLH1 removed from regulons | **0.561** |
| random size-matched gene sets (12 runs) | mean **0.466**, range 0.315–0.574 |

The de-leaked model is matched or exceeded by 1 of 12 random-regulon models. In other words, once
the leak is removed, knowing which genes a transcription factor actually regulates buys little over
knowing nothing at all: what is being measured is that FOLH1 is predictable from the transcriptome,
not that it is predictable from TF *activity*.

**Even plainer baseline.** Replacing the 766 TF activities with the raw expression of 766
*randomly chosen genes*, same elastic net, same TCGA→SU2C transfer, gives external R² = 0.509 on
average over 8 repeats (range 0.388–0.599). Predicting a gene's expression from several hundred
other transcriptome-derived features is not a hard problem, and the manuscript contains no null
model that establishes otherwise.

**How much is lineage?** Simple lineage scores transfer as follows (same z-space, TCGA→SU2C):
AR score alone R² = 0.382, NE score alone R² = −0.083, AR + NE R² = 0.290 (worse than AR alone —
the two-feature model does not transfer). The TF model's prediction correlates with observed FOLH1
at r = 0.817, and **r = 0.683 after partialling out the AR and NE scores**. So, in fairness to the
authors, the TF model does carry substantial signal beyond the two lineage scores — it is not
merely an AR read-out. What it is not is evidence of *regulatory* architecture, which is what §3.1's
first two tests address.

**Clinical performance, stated in clinical terms.** The manuscript reports R² and ρ and never a
metric a clinician could act on. Converting:

- AUC for identifying the bottom-FOLH1 quartile in SU2C: **0.929**
- Calling the model's own bottom quartile "PSMA-low": PPV 0.79, sensitivity 0.79
- Per-tumour uncertainty: residual SD 0.582 SD units → a 95% prediction interval of **±8.5-fold**
  (a 72-fold span) in FPKM, with no absolute calibration across platforms (raw-space transfer
  R² = **−11.98**, reproducing the authors' own "naive transfer R² = −12")

The ordering is good; the individual prediction is not usable. And the deeper problem is the use
case: the model predicts FOLH1 RNA from TF activities **inferred from the same RNA sample**. Anyone
who can compute the predictors has already measured the outcome. There is no clinical setting in
which this model would be run — which is exactly Joe's "does not establish ... clinically useful
predictive performance". If the model is kept, it should be framed as a descriptive statement about
regulatory co-variation, not as a predictor.

### 3.2 The external validation confirms lineage, not the proposed construct (Joe's point 2)

**Head-to-head against the standard markers** (Beltran, NEPC vs CRPC-adeno):

| classifier | AUC |
|---|---|
| 12-gene neuroendocrine score (standard panel) | **0.929** |
| 100-gene PSMA-low signature (this paper) | 0.820 |
| 20-gene AR activity score, inverted | 0.804 |
| FOLH1 expression alone, inverted | 0.782 |

The signature is a *worse* NEPC classifier than a routine NE panel. This is robust to all four ways
of handling the Beltran matrix's negative expression floor (NE 0.924–0.939 vs signature
0.814–0.833).

**The GSE126078 "validation" is close to uninformative.** A signature derived from the quartiles of
*any* gene will anti-correlate with that gene in another cohort. Repeating the whole pipeline for 15
random anchor genes gives |ρ| = 0.38 to 0.85 in GSE126078 patient samples (median 0.65); 1 of the
15 exceeds the published FOLH1 value (0.824), which sits inside that null. The Beltran AUC *does* exceed the
random-anchor null (median 0.51, max 0.77; 0 of 15 reach 0.820), so the signature carries real lineage information — but
that information is already captured, better, by existing markers.

**Tumour composition is not controlled** (this is Joe's phrasing, and it is testable):

- FOLH1 tracks **epithelial content**: ρ = +0.287 (p = 1.9e-6). PSMA-low tumours have lower tumour
  content (epithelial score −0.47 vs +0.33, p = 1.6e-6).
- The PSMA-low *score* tracks **immune** (ρ = +0.346) and **stromal** (ρ = +0.250) content.
- The hepatocyte adjustment does not remove the site confound: the signature still identifies liver
  biopsies after adjustment (AUC 0.724), and FOLH1 still correlates with the hepatocyte score at
  ρ = −0.494 **within non-liver samples** — i.e. the covariate used for adjustment is itself tied to
  the outcome, so residualizing all 12,138 genes on it removes real signal as well as contamination.

**But the state is not merely liver + NEPC.** Excluding liver biopsies and histologically NE
tumours (n = 215) the contrast gets *stronger*, not weaker: 2,469 DE genes, rbc concordance 0.896
with the full analysis, 884 of the original 1,002 genes retained. A real transcriptional contrast
exists. What is not established is that it is a *clinical* phenotype distinct from AR-low /
NE-like disease.

**Cohort-structure issues:**

- GSE126078's 98 patient samples come from **55 patients** (up to 3 metastases each). The reported
  p = 1.1e-29 / 1.8e-25 is inflated by pseudo-replication; the correlation itself is unaffected.
- The 12 "cell lines" pooled into the GSE126078 correlation are **siRNA/knockout perturbation
  experiments** (C42B/PC3/PacMetUT1 with siREST or siNCT). REST knockdown is a deliberate
  neuroendocrine-differentiation manipulation. They should not be in a validation cohort at all.
- SU2C pseudo-replication is *not* an issue: 266 samples from 262 patients; all key effect sizes
  are unchanged at the patient level (e.g. STEAP1 rbc −0.821, PLAUR +0.499).

**Is there a PSMA-low subpopulation?** FOLH1 in SU2C: bimodality coefficient 0.632 (> 0.555
suggests bimodality), Shapiro–Wilk p = 2.4e-11, Q1-to-Q3 span 11.7-fold. There is weak evidence of
a low-expressing tail, but the quartile cut-off remains arbitrary — Joe's suggestion to treat FOLH1
continuously is well founded, and the data support reporting both.

**The signature is not prognostic either.** The manuscript reports the FOLH1 null (HR 1.03/1.06).
The signature score itself is equally null: HR 1.03 per SD (0.74–1.46), p = 0.84. Worth stating
explicitly, because it closes the prognostic claim rather than leaving it open.

### 3.3 The therapeutic conclusions are contradicted by the data, not merely unsupported (Joe's point 3)

**(a) Tumour-to-normal contrast — the quantity that actually matters for a ligand.** The manuscript
ranks candidates by rank-biserial correlation between two tumour groups, which says nothing about
whether a ligand can be delivered. Computing median tumour expression minus the highest vital-organ
expression (HPA consensus, log2 units):

| target | contrast in PSMA-low tumours |
|---|---|
| **FOLH1 in PSMA-high tumours (reference)** | **+2.9** |
| CD276 (B7-H3) | **+0.3** |
| GPC2 | −0.4 |
| STEAP1 | −0.7 |
| STEAP2 | −1.3 |
| SLC12A4 | −2.2 |
| SLC16A3 (MCT4) | −2.4 |
| ATP2B4 | −2.5 |
| SLC7A5 (LAT1), ADAM19, MRC2 | −2.8 |
| TACSTD2 (TROP2) | −3.7 |
| NCAM1 (CD56) | −6.4 |
| **PLAUR (uPAR)** | **−6.8** |
| MUC1 | −7.0 |

Every de-novo nominee is out-expressed by at least one vital organ. uPAR — the flagship candidate,
proposed with a clinical imaging modality — is 6.8 log2 units (≈110-fold) below its worst normal
organ. GPC2 is the only nominee that is close to neutral, and its absolute tumour level is 0.26% of
FOLH1's level in PSMA-high disease.

**(b) "Luminal antigens fail" is the opposite of what the data show.** In PSMA-low tumours the
*highest-expressed* surface targets are exactly the ones the manuscript rules out:

| target | median in PSMA-low tumours (log2 FPKM+1) | fraction of PSMA-low samples ≥ 2 |
|---|---|---|
| TACSTD2 (TROP2) | 5.45 | 87% |
| CD276 (B7-H3) | 4.46 | 93% |
| STEAP1 | 4.14 | 78% |
| SLC7A5 (LAT1) | 3.91 | 90% |
| SLC16A3 (MCT4) | 2.94 | 78% |
| PLAUR (uPAR) | 2.57 | 73% |
| **FOLH1** | **1.96** | **49%** |
| GPC2 | 1.32 | 30% |

Statistically down-regulated is not the same as absent. STEAP1 falls from 7.75 to 4.14 — it is
reduced, but it remains more abundant than every newly nominated target, and B7-H3 is the only
surface protein in the panel with positive tumour-to-normal contrast in PSMA-low disease. The
manuscript's inference ("co-suppressed with PSMA, therefore prevalence-based nomination
overestimates their coverage") does not follow from a rank-biserial correlation.

**(c) The nominated targets are substantially non-tumour-cell genes.** HPA single-cell consensus,
fraction of a gene's total expression in immune cell types:

PLAUR **38.6%** (top cell types: neutrophils, monocytes) · ADAM19 22.4% · SLC16A3 19.2% ·
SLC7A5 19.0% — versus FOLH1 **0.2%** (prostatic glandular cells), STEAP1 0.3%, STEAP2 0.2%,
TACSTD2 0.3%. MRC2's top cell types are fibro-adipogenic progenitors and bone-marrow MSCs.

Residualizing on tumour (epithelial) content instead of hepatocyte content shrinks the effect
sizes accordingly: VCAN −58% (and loses significance), MRC2 −47%, MXRA5 −45%, ADAM19 −33%,
PLAUR −25%. MCT4 (−9%), LAT1 (slightly *stronger*, +4%) and GPC2 (−19%) are robust, and STEAP1/2 become stronger still. So
the tumour-composition problem is real and selective — it does not sink MCT4/LAT1/GPC2, but it
does sink uPAR, MRC2, VCAN and MXRA5 as tumour-cell targets.

**(d) The TXNRD1 dependency does not survive a direct test.** Three independent checks:

1. **Lineage.** 29.8% of the signature score's variance across DepMap is explained by lineage
   alone (Kruskal–Wallis H = 388 on 22 df, p ≈ 1e-68). Prostate is the lowest-scoring of 23 lineages (mean −1.03; next
   lowest, breast, −0.28); the highest are myeloid (+0.64) and lymphoid (+0.23). That NCKAP1L —
   a haematopoietic-restricted gene — appears among the top hits is the signature of this artefact.
   The four lowest-scoring models in all of DepMap (1,673 models) are exactly the four
   FOLH1-expressing prostate adenocarcinoma lines — MDA-PCa-2b −3.97, LNCaP −3.67, VCaP −3.55,
   22Rv1 −1.96, with a gap of 0.37 to the fifth-lowest line of any lineage. Spearman ranking limits
   how much four outliers can drive the result, but it does not remove the lineage structure.
   **Surrogate fidelity:** across DepMap the signature score correlates with FOLH1 expression itself
   at only ρ = −0.179 — in cell lines the surrogate is only weakly related to the construct it
   stands for.
2. **Within-lineage re-analysis.** Demeaning both score and gene effect within lineage
   (22 lineages, 1,070 models):

   | gene | pan-lineage ρ (rank) | within-lineage ρ (rank) |
   |---|---|---|
   | TXNRD1 | −0.227 (2) | **−0.186 (1)** |
   | SEPSECS | −0.234 (1) | −0.147 (9) |
   | EEFSEC | −0.215 (4) | −0.125 (38) |
   | PSTK | −0.186 (15) | −0.095 (265) |
   | GPX4 | −0.188 (13) | −0.081 (602) |
   | **GAPDH (housekeeping benchmark)** | −0.195 (10) | **−0.133 (24)** |

   The "selenocysteine/thioredoxin axis" collapses to TXNRD1 alone, and a housekeeping gene
   outranks three of its four members.
3. **The construct itself.** Correlating **FOLH1 mRNA** directly with dependency — the thing the
   claim is actually about — gives TXNRD1 ρ = **−0.090**, i.e. the *wrong sign*: FOLH1-**high**
   lines are the more TXNRD1-dependent ones. SEPSECS +0.015, EEFSEC −0.007, PSTK +0.002 are null.
   Restricting to the 150 FOLH1-expressing lines: TXNRD1 ρ = −0.099, still the wrong sign.
   Prostate lines only (n = 10): FOLH1 vs TXNRD1 ρ = +0.16, p = 0.66.

**(e) Effect size.** ρ = −0.227 is 5.1% of variance. For calibration, the correlation between a
gene's own expression and its own dependency across the same 1,103 lines: CCND1 −0.63, MDM2 −0.46,
ERBB2 −0.35, EGFR −0.31, FOXA1 −0.21. TXNRD1's own expression versus its own dependency is **+0.15**
— TXNRD1 dependency is not even predicted by TXNRD1 expression. In prostate lines GPX4 shows much
deeper dependency (Chronos −0.24 to −1.54) than TXNRD1 (0.00 to −0.70).

The manuscript already hedges this ("the top of a broader fitness/lineage structure rather than a
unique mechanism"). The hedge is not sufficient: the direct test points the other way.

---

## 4. Other issues found

| # | Issue | Severity |
|---|---|---|
| 1 | **Live AlphaGenome API key hard-coded** in the delivered notebook (cell 41) | **Act now** — revoke |
| 2 | TCGA "primary, n=550" includes **52 normal-prostate samples and 1 metastasis**. All TCGA statements rest on this cohort; FOLH1–AR ρ falls 0.332 → 0.262 in tumours only (conclusions unchanged, the number is not) | Must fix |
| 3 | GSEA permutation p-values not reproducible as p < 1e-3 (correct: p ≈ 0.01, FDR ≈ 0.02 at 10,000 permutations) | Must fix |
| 4 | "ChIP-Atlas confirms direct binding" — technically true, but those target lists contain 10,961 (AR) to 17,578 (EZH2) genes; the EZH2 list is ~88% of all genes. FOLH1's average binding score is 17.2 for AR vs 51.1 FOXA1 and 58.2 HOXB13. As written this is presented as mitigation for inferred TF activity; it is close to uninformative for AR and EZH2 | Must soften |
| 5 | Beltran matrix: 59% of values are below the pipeline's negative floor. Clipping at 0 makes 6 of the 100 signature genes constant (silently dropped from the score) and zeroes **every** SYP and MYCN value, degrading the NE panel in that cohort. Headline AUC is robust (0.814–0.833) so this is not consequential, but the NE comparisons in Beltran are affected | Fix, low impact |
| 6 | Two different hepatocyte covariates are used: gene-wise z-scored mean for the main DE, unscaled mean for the known-target panel. Both replicate; the inconsistency should go | Cosmetic |
| 7 | 12 siRNA/knockout perturbation samples pooled into the GSE126078 validation; 98 patient samples come from 55 patients | Fix |
| 8 | The Gleason claim ("non-monotonic, drop at GS10") rests on **4 patients** with GS10 | Must soften |
| 9 | Results §1 "ARSI exposure did not alter tumor FOLH1 (p=0.91)" is a cross-sectional comparison of different patients, not an exposure effect | Wording |
| 10 | Abstract reports the external R² without a null model and without a clinically interpretable metric; the model's inputs are derived from the same RNA sample as its output | Framing |

---

## 5. Bottom line

**What survives, unchanged:**

- FOLH1 expression is strongly and reproducibly coupled to AR activity, increasingly so with
  progression, and inversely to neuroendocrine programs. Every number holds.
- A coherent, reproducible transcriptional contrast distinguishes FOLH1-low from FOLH1-high
  metastatic tumours, and it is **not** reducible to liver-biopsy contamination or to
  histologically NE tumours (it strengthens when both are removed).
- Locus methylation adds essentially nothing to TF activity in primary tumours (R² 0.714 → 0.719;
  methylation alone 0.108). This is a genuine, well-supported negative result.
- FOLH1 mRNA is not prognostic for OS in SU2C, and neither is the signature. A clean informative
  null.
- The surfaceome differential expression is correct as computed, and the known-target audit
  (DLL3 null, NCAM1 subset-driven, CEACAM5 liver-driven, SSTR2 borderline) is honest and holds.

**What does not survive:**

- "Quantitatively predictable from TF activity" overstates the result: 15% of the external R² comes
  from FOLH1 sitting inside its own predictors (0.662 → 0.561), and what remains is matched by
  size-matched random gene sets (mean 0.466, max 0.574) and by 766 random genes' raw expression
  (mean 0.509). The model does carry signal beyond AR/NE lineage scores (partial r = 0.683), but
  nothing in the analysis shows that signal to be transcription-factor *architecture*, and the
  model has no clinical use case because its inputs are derived from the same RNA as its output.
- "Externally validated PSMA-low state" overstates the external evidence: the GSE126078 correlation
  is inside a random-anchor null, and the Beltran AUC is beaten by a standard 12-gene NE panel.
- "Luminal alternatives fail; MCT4/LAT1/uPAR/GPC2 emerge" is contradicted on the metric that
  matters: in PSMA-low tumours TROP2, B7-H3 and STEAP1 remain the best-expressed surface targets,
  B7-H3 is the only one with positive tumour-to-normal contrast, and uPAR is a myeloid gene that is
  110-fold out-expressed by normal tissue.
- "PSMA-low-like cells are selectively dependent on the selenocysteine/thioredoxin axis" does not
  survive lineage adjustment (axis → TXNRD1 alone, barely above a GAPDH benchmark) and reverses
  sign when tested against FOLH1 expression itself.

**Assessment of Joe's review.** All three of his concerns are correct, and two of them are more
serious than he states. He infers the TF circularity conceptually; it is also literal (FOLH1 is a
regulon member). He says the therapeutic conclusions "exceed the evidence"; on tumour-to-normal
contrast and cell-type origin the evidence points the other way. His proposed reframing — *what
distinguishes PSMA-low metastatic prostate cancer biologically, and which alternative targets are
consistently retained or enriched* — is the question this dataset can actually answer, and §3.3(b)
is a preview of what that analysis yields.

**Recommended posture for the reply:** concede the three points without defensiveness, report the
replication as evidence that the computational work is sound, and present the reframed target
analysis (retention + tumour-to-normal contrast + cell-type origin) as the new contribution. Do not
defend the TF model as a predictive tool; keep it, de-leaked, as descriptive regulatory context.

Phase 2 (restructuring and re-analysis along these lines) is scoped separately and has not been
started.
