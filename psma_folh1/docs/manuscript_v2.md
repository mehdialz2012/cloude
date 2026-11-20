# What distinguishes PSMA-low metastatic prostate cancer, and which surface targets remain available in it

**Draft v2 — restructured after review.** This replaces the v1.1 manuscript *"An externally
validated transcription-factor model of PSMA (FOLH1) expression…"*. The transcription-factor
model, the survival analysis and the AlphaGenome section are removed from the main narrative for
the reasons set out in `docs/reply_to_joe.md`; every number below comes from the re-analysis in
`code/p2_*.py`, which treats FOLH1 as a continuous variable and adjusts for biopsy site, tumour
composition and neuroendocrine differentiation throughout.

---

## Abstract

**Background.** PSMA-targeted radioligand therapy is standard of care in metastatic
castration-resistant prostate cancer (mCRPC), but 10–20% of patients are PSMA-PET-negative at
baseline and most harbour at least one PSMA-negative metastasis. What distinguishes this
compartment biologically, and which cell-surface targets remain available in it, is unresolved.

**Methods.** FOLH1 was analysed as a continuous variable in SU2C/PCF mCRPC (n=266) with explicit
adjustment for biopsy site, tumour composition (PanglaoDB cell-type compartments) and
neuroendocrine differentiation. Associations were validated in GSE126078 (98 metastases from 55
patients, with published molecular-phenotype calls) at both sample and patient level, and in the
Beltran CRPC/NEPC cohort. A curated plasma-membrane surfaceome (3,702 genes) was screened against
four pre-specified suitability criteria — abundance in PSMA-low tumours, tumour-to-normal contrast
benchmarked against PSMA's own, tumour-cell specificity from single-cell data, and protein-level
evidence from immunohistochemistry. DepMap was used only for exploratory hypothesis generation,
and TCGA-PRAD (n=497 tumours) only to ask whether the programme is specific to advanced disease.

**Results.** Of 6,438 genes associated with FOLH1 in crude analysis, 2,140 are explained by tumour
composition or biopsy site and a further 3,232 by neuroendocrine differentiation; 1,557 remain.
Of these, 251 reproduce in an independent patient cohort, and 191 reproduce in all three. The
validated programme is narrow: proliferation up (E2F NES +2.9, G2M +2.6), androgen response and
oxidative phosphorylation down. Epithelial-mesenchymal transition, inflammation and hypoxia —
prominent in crude analysis — reverse sign on adjustment and are properties of the biopsy, not of
the tumour cells. In GSE126078, PSMA-low is overwhelmingly AR-low or AR-negative disease: only 6
of 60 AR-positive, non-neuroendocrine metastases fall in the cohort's bottom FOLH1 quartile, and
within that subgroup only 15 genes are associated with FOLH1. On the surfaceome, most antigens are
retained rather than lost (782 of 1,184 unchanged genes remain expressed in ≥50% of PSMA-low
tumours), but only two of the 36 evaluated targets meet all suitability criteria, and only one —
B7-H3 (CD276) — is a credible therapeutic candidate: reduced but retained in 92.5% of PSMA-low
tumours, within 1.2 log2 units of PSMA's own tumour-to-normal contrast, 14% non-tumour-cell
expression, and positive by immunohistochemistry in all 32 prostate cancers tested. uPAR, MCT4, LAT1
and GPC2 each fail at least two criteria. Exploratory dependency mapping with the validated
programme nominates nucleotide biosynthesis and p53 status, not a selective antioxidant
vulnerability. The programme is already present in treatment-naive primary tumours (98% of genes
with the same direction), so it is a property of low-PSMA prostate cancer rather than of
treatment resistance.

**Conclusions.** The PSMA-low metastatic compartment is, in the main, AR-low/neuroendocrine
disease with a proliferative, AR-axis-low transcriptional programme. It is not covered by the
surface targets most often proposed for it; B7-H3 is the only antigen in this analysis that both
survives independent validation and satisfies basic therapeutic-suitability criteria, and it is
reduced rather than enriched in PSMA-low tumours. Alternative-target nomination for this
compartment should be based on retention and tumour-to-normal contrast, not on differential
expression.

---

## 1. Introduction

PSMA-PET governs patient selection for radioligand therapy, and ¹⁷⁷Lu-PSMA-617 prolongs survival
in mCRPC. The therapeutic question that follows is not what regulates PSMA, but what to do for the
patients whose tumours do not express enough of it. That compartment is growing: 10–20% of
patients are PSMA-PET-negative at baseline, most patients harbour at least one PSMA-negative
metastasis, and lineage plasticity under androgen-receptor signalling inhibition produces
AR-indifferent and neuroendocrine states in which PSMA is characteristically suppressed.

Three things are needed before an alternative target can be proposed for these patients: a
description of what the PSMA-low state actually is once the artefacts of bulk metastatic biopsy
are removed; evidence that the description reproduces in independent patient cohorts; and an
assessment of candidate targets on the criteria that determine whether a ligand can be delivered,
rather than on differential expression alone. This paper addresses those three in turn, and then
asks, exploratorily, whether the state carries any functional vulnerability, and whether it is
specific to advanced disease.

---

## 2. Methods (summary)

**Cohorts.** SU2C/PCF `prad_su2c_2019` (polyA FPKM, n=266 from 262 patients) is the discovery
cohort. GSE126078 contributes 98 patient metastases from 55 patients with published
molecular-phenotype calls (ARpos/ARlow/ARneg × NEpos/NEneg) and tumour-site annotation; the 43
LuCaP PDX and 12 siRNA-perturbation cell-line samples in that series are excluded from validation.
Beltran `nepc_wcm_2016` (n=49; 34 CRPC-adenocarcinoma, 15 NEPC) is used for the neuroendocrine
question only. TCGA-PRAD contributes 497 primary tumours (the 52 normal-prostate and one
metastatic sample in the Xena matrix are excluded).

**Statistical model.** FOLH1 is used as a continuous, within-cohort z-scored variable. For every
expressed gene three nested ordinary-least-squares models are fitted:

| model | covariates |
|---|---|
| M1 | FOLH1 |
| M2 | M1 + immune, stromal, hepatocyte and osteoblast compartment scores + biopsy site |
| M3 | M2 + neuroendocrine score |

Compartment scores are gene-wise z-scored means over PanglaoDB human cell-type markers. Each gene
is labelled by the model in which its association disappears, so that confounding is reported
rather than removed silently. Quartile, tercile, median and decile contrasts are reported
alongside as sensitivity analyses. Programme-level inference uses Hallmark GSEA with 10,000
permutations.

**Validation.** A gene is counted as reproduced when it is significant in SU2C M3 (FDR<0.05) and
has the same sign at nominal p<0.05 in the validation cohort. GSE126078 is analysed both at sample
level and after aggregating to one observation per patient.

**Surfaceome and suitability.** The surfaceome is SURFY ∪ CSPA high-confidence ∪ UniProt-reviewed
cell-membrane proteins with ≥1 transmembrane helix or GPI anchor (3,702 genes). Plasma-membrane
confidence requires UniProt curation or agreement between SURFY and CSPA. Four pre-specified
criteria are applied:

1. **Abundance** — expressed at log2(FPKM+1) ≥ 2 in ≥50% of PSMA-low tumours.
2. **Selectivity** — tumour-to-normal contrast (median tumour expression minus the highest
   vital-organ nTPM in HPA) within 2 log2 units of PSMA's own contrast in PSMA-high tumours.
   Platform offsets cancel in this comparison; absolute contrast values do not transport.
3. **Tumour-cell specificity** — <25% of a gene's single-cell expression in immune or stromal cell
   types (HPA single-cell consensus).
4. **Protein evidence** — detected by immunohistochemistry in ≥50% of prostate cancers (HPA),
   reported with the antibody-reliability flag because this evidence is antibody-dependent.

**Exploratory analyses.** DepMap 24Q4 dependencies are correlated with the *validated* programme
score, reported pan-lineage and after demeaning within lineage. TCGA is used only to compare
effect directions between primary and metastatic disease.

---

## 3. Results

### 3.1 Three quarters of what looks like a PSMA-low programme is confounding (Figure 1)

FOLH1 is continuous in mCRPC, with a long low tail: SD 2.71 log2 units and a 12-fold span between
the first and third quartiles. All four cut-offs tested agree closely with the continuous model
(ρ 0.91–0.94) and with each other (ρ 0.85–0.98), so the choice of threshold is not what drives
the result — but the continuous model is reported throughout.

In a bulk metastatic biopsy, FOLH1 tracks lineage first and composition second: AR activity
ρ=+0.55, neuroendocrine score ρ=−0.34, osteoblast content ρ=−0.20, immune content ρ=−0.18,
hepatocyte score ρ=−0.13, stromal content ρ=−0.11, and biopsy site differs strongly across FOLH1
levels (Kruskal p=1.2e-7).

Adjusting for these in turn removes most of the apparent programme. Of 6,438 genes associated with
FOLH1 in crude analysis, 2,140 do not survive adjustment for composition and biopsy site, and a
further 3,232 do not survive adjustment for neuroendocrine differentiation. **1,557 genes (24%)
remain.** Median effect sizes shrink by 23.5% and a further 35.4% across the two adjustments.

The programme-level consequences are specific. Proliferation (E2F NES +2.9, G2M +2.6) and loss of
androgen response (−2.3), protein secretion (−1.8) and oxidative phosphorylation (−1.7) are
unaffected by adjustment. Epithelial-mesenchymal transition **reverses sign** (+1.7 → −1.2), as do
inflammatory response (+1.4 → −1.3) and hypoxia (+1.5 → −1.0). The EMT-high, inflamed phenotype
that crude analysis attributes to PSMA-low tumours is a property of what else is in the biopsy.

### 3.2 Independent validation narrows the programme to 251 genes, and places it in AR-low disease (Figure 2)

In GSE126078, where metastases carry published molecular-phenotype calls, the PSMA-low compartment
is almost entirely AR-low or AR-negative disease:

| phenotype | n | median FOLH1 | in the cohort's bottom FOLH1 quartile |
|---|---|---|---|
| ARpos_NEpos | 11 | 7.86 | 0 / 11 |
| ARpos_NEneg | 60 | 7.14 | 6 / 60 |
| ARlow_NEneg | 9 | 5.03 | 4 / 9 |
| ARneg_NEpos | 10 | 1.60 | 9 / 10 |
| ARneg_NEneg | 8 | 1.23 | 6 / 8 |

Within the 60 AR-positive, non-neuroendocrine metastases — the compartment where a PSMA-low
phenotype would be most clinically surprising and most therapeutically interesting — only 15 genes
are associated with FOLH1 at FDR<0.05, and their effect sizes agree only weakly with the discovery
cohort (ρ=0.32). We therefore do not claim a distinct AR-intact PSMA-low biology; the data do not
support one at this sample size.

Applying the identical nested models in GSE126078 reproduces the attrition: 3,749 genes crude, 821
after composition and site, 144 after neuroendocrine adjustment (36 at patient level). Of the 715
SU2C-robust genes measured in all three cohorts, **251 (35%) reproduce in GSE126078, 202 at
patient level, and 191 in all three cohorts.** Effect sizes correlate at ρ=0.48 between SU2C and
GSE126078 (0.40 patient-level) and ρ=0.36 with Beltran.

At programme level the agreement is better than at gene level (NES ρ=0.61 between SU2C and
GSE126078): proliferation up, androgen response down, in both. Beltran does not agree (NES
ρ=0.12), which is expected — in a 49-sample cohort that is 31% NEPC, FOLH1 and the neuroendocrine
score are collinear (r=−0.46) and conditioning on NE acts as a suppressor (median |t| for the
FOLH1 term rises from 1.38 to 2.68). Beltran is therefore reported for the neuroendocrine question
only: FOLH1 median 0.00 in NEPC versus 4.89 in CRPC-adenocarcinoma (p=1.3e-3), with 11 of 15 NEPC
and 8 of 34 adenocarcinomas in the cohort's bottom quartile.

### 3.3 Most surface antigens are retained, and only one is a credible alternative target (Figure 3)

Of 1,425 surfaceome genes expressed in mCRPC, **55 are enriched in PSMA-low tumours, 186 are lost
and 1,184 are unchanged** after full adjustment. The therapeutically relevant observation is not
the 55: it is that 782 of the 1,184 unchanged genes, and 146 of the 186 "lost" genes, are still
expressed in at least half of PSMA-low tumours. Loss of PSMA is not accompanied by loss of the
surface proteome.

Scoring the 36 evaluated targets against the four pre-specified criteria, with PSMA's own
tumour-to-normal contrast (+1.21 log2) as the benchmark:

| target | state in PSMA-low | reproduced | abundance (median, % positive) | contrast vs PSMA | non-tumour fraction | PCa IHC | tier |
|---|---|---|---|---|---|---|---|
| **CD276 (B7-H3)** | reduced | 1 cohort | 4.46, 92.5% | **−1.24** | 0.14 | 32/32, all medium/high (Enhanced) | **A** |
| STEAP1 | reduced | 1 cohort | 4.14, 77.6% | −2.13 | 0.06 | 1/9 (antibody Uncertain) | C |
| STEAP2 | reduced | 1 cohort | 3.05, 62.7% | −2.54 | 0.04 | 12/12 (Uncertain) | C |
| TACSTD2 (TROP2) | unchanged | 0 | 5.45, 86.6% | −6.04 | 0.01 | 18/32 | C |
| SLC7A5 (LAT1) | enriched | **0** | 3.91, 89.6% | −4.02 | 0.23 | 12/12 | C |
| SLC16A3 (MCT4) | enriched | 1 cohort | 2.94, 77.6% | −3.64 | 0.21 | 10/21 | C |
| GPC2 | enriched | **0** | 1.32, **29.9%** | −1.57 | 0.09 | no data | D |
| PLAUR (uPAR) | **not significant** | 0 | 2.57, 73.1% | **−7.99** | **0.40** | **0/9** | E |
| PSCA | unchanged | 0 | 1.84, 49.3% | −10.92 | 0.00 | 0/12 | E |
| FOLH1 (reference) | lost | 2 cohorts | 1.96, 49.3% | −7.23 | 0.01 | 18/18 | E |

Three conclusions follow.

**First, the alternatives proposed on the basis of differential expression do not hold up.** uPAR
loses significance entirely once composition is adjusted for (FDR 0.38), draws 40% of its
single-cell expression from immune and stromal cells, is 8 log2 units below PSMA on
tumour-to-normal contrast, and is undetectable by immunohistochemistry in 9 of 9 prostate cancers.
LAT1 and GPC2 do not reproduce in any independent cohort; GPC2 is expressed above threshold in
fewer than a third of PSMA-low tumours. MCT4 reproduces in one cohort but is 3.6 log2 units below
PSMA on contrast with 21% non-tumour-cell expression.

**Second, the luminal antigens are reduced, not absent, and remain the best-expressed options.**
TROP2 (median 5.45, 87% positive), B7-H3 (4.46, 93%) and STEAP1 (4.14, 78%) are all more abundant
in PSMA-low tumours than any newly nominated candidate, and all three are more abundant there than
FOLH1 itself (1.96, 49%).

**Third, only B7-H3 satisfies every criterion.** It is reduced in PSMA-low disease — this is a
retention argument, not an enrichment argument — but it is retained in 92.5% of PSMA-low tumours,
its tumour-to-normal contrast is within 1.2 log2 units of PSMA's own, only 14% of its single-cell
expression is non-tumour, and it is positive by immunohistochemistry in every prostate cancer
tested at medium or high intensity with an "Enhanced"-reliability antibody. Clinical agents exist
(ifinatamab deruxtecan, CAR-T, radiolabelled omburtamab). The second gene meeting all criteria,
ABCB6, is a porphyrin transporter with no surface-targeting rationale and is reported for
completeness only.

A caution that applies throughout: HPA immunohistochemistry is antibody-dependent. STEAP1 and
PSCA, both validated clinical targets, score 1/9 and 0/12 positive here with antibodies flagged
"Uncertain" and "Approved" respectively. Protein evidence is therefore carried as supporting
rather than decisive, and a negative IHC result with an unreliable antibody is not treated as
evidence of absence.

### 3.4 Exploratory: the validated programme nominates proliferation, not an antioxidant vulnerability (Figure 4)

Scoring DepMap 24Q4 with the 251-gene validated programme rather than a discovery-only signature
changes the dependency landscape completely. The selenocysteine/thioredoxin axis, which ranks 1,
2, 4 and 15 of 17,293 non-essential genes when the unvalidated signature is used, falls to ranks
6,267 (SEPSECS), 1,709 (TXNRD1), 3,201 (EEFSEC) and 11,337 (PSTK); GPX4 falls from 13 to 9,881.
The top dependencies of programme-high cells are instead TP53BP1, TP53 and purine/one-carbon
biosynthesis (ADSL, MTHFD1, PPAT, GART, PFAS, CTPS1) — the expected correlates of a proliferative,
p53-variable state rather than a selective vulnerability.

Two further observations argue against a therapeutic reading. Lineage alone explains 38.4% of the
programme score across DepMap. And tested directly against FOLH1 mRNA — the construct the claim is
about — the correlations are null or of the wrong sign (TXNRD1 ρ=−0.09, i.e. FOLH1-high lines are
marginally the more dependent). We therefore report no therapeutic dependency for the PSMA-low
state.

### 3.5 Exploratory: the programme is not specific to advanced disease (Figure 4c)

FOLH1 heterogeneity nearly doubles with progression (SD 1.47 in TCGA primary tumours versus 2.71
in mCRPC). But the programme itself is already present before treatment: effect sizes correlate at
ρ=0.55 between primary and metastatic disease, and 242 of the 247 validated-programme genes
measurable in TCGA (98%) have the same direction there, 229 of them nominally significant. Loss of
androgen response is consistent across all four cohorts (NES −2.3, −2.7, −1.8, −2.0); proliferation
is weaker in primary disease (+1.4 versus +2.9).

PSMA-low prostate cancer is therefore better described as a state that exists across disease
stages and becomes more prevalent and more extreme with progression, than as a state created by
treatment resistance.

---

## 4. Discussion

The useful finding here is a negative one, and it is actionable. The surface targets most often
proposed as PSMA alternatives were nominated by differential expression in bulk tumours. When
FOLH1 is modelled continuously and the analysis is adjusted for what else is in a metastatic
biopsy, three quarters of the differential expression disappears, two thirds of the remainder
fails to reproduce in an independent patient cohort, and the specific candidates that survive both
filters fail on the criteria that determine whether a ligand can be delivered: abundance,
tumour-to-normal contrast, and whether the antigen is on the tumour cell at all.

The reframing that follows is that **alternative targets for PSMA-low disease should be sought
among antigens that are retained, not among antigens that are enriched.** Enrichment in a
compartment defined by loss of a luminal marker selects for genes that are high in non-luminal
cells — which in a bulk biopsy includes the stroma and the infiltrate. Retention plus contrast is
the property that matters, and on that basis B7-H3 is the leading candidate in this dataset, with
STEAP1 and STEAP2 behind it. That B7-H3 is *reduced* in PSMA-low tumours is not an argument
against it; it is retained in 93% of them, at a contrast close to PSMA's own.

Two clinical implications follow. First, because the PSMA-low compartment is overwhelmingly AR-low
or neuroendocrine disease, a targeting strategy for it is largely a targeting strategy for
AR-indifferent disease, and should be evaluated in that population. Second, because only 6 of 60
AR-positive, non-neuroendocrine metastases were PSMA-low, the clinically discordant patient —
AR-intact but PSMA-PET-negative — is rare in transcriptomic cohorts and will need lesion-level,
imaging-anchored sampling to characterise. No public cohort pairs PSMA-PET with transcriptomes;
that is the single most useful dataset the field is missing.

### Limitations

Bulk RNA is a proxy for an antigen-density question that is properly answered at the protein and
lesion level. Tumour composition is adjusted for statistically, not resolved by deconvolution or
single-cell sequencing. The tumour-to-normal contrast metric compares FPKM against HPA nTPM and is
interpretable only relative to PSMA's own contrast, which is why it is reported that way. HPA
immunohistochemistry is antibody-dependent and is used as supporting evidence. The
validation cohort contributes 98 metastases from 55 patients, and patient-level analysis reduces
the validated programme further (36 genes at FDR<0.05), so the 251-gene programme should be read
as a sample-level result with patient-level support rather than as a fixed signature. DepMap and
TCGA analyses are exploratory by design and are reported as such. Finally, the central claim is a
prioritisation, not a demonstration: B7-H3 is nominated on retention, contrast, specificity and
protein evidence, and requires lesion-level confirmation in PSMA-PET-negative disease before it
can be called a solution for this compartment.

---

## Figures

- **Figure 1** `p2_fig1_discovery` — FOLH1 as a continuous variable; what it tracks in a bulk
  biopsy; attrition of the association under adjustment; programme reversals.
- **Figure 2** `p2_fig2_validation` — FOLH1 by molecular phenotype; cross-cohort effect-size
  agreement; the 6,438 → 251 funnel; programme concordance across four cohorts.
- **Figure 3** `p2_fig3_targets` — retained/lost/enriched surface genes; abundance versus
  selectivity; the pre-specified criteria applied to established and newly nominated targets.
- **Figure 4** `p2_fig4_depmap_tcga` — the selenocysteine axis before and after validation; what
  the validated programme nominates instead; primary versus metastatic effect sizes.

## Data and code

All cohorts are public. The analysis is in `psma_folh1/code/p2_*.py`; machine-readable outputs,
including the full surfaceome table and the target-prioritisation table, are in
`psma_folh1/results/`.
