# Manuscript integration recommendations

**Scope.** This file proposes where and how the AlphaGenome 2 analysis could
be incorporated into *"Ghrelin Receptor Gene Polymorphisms, Cognitive
Screening Status, and Osteoporosis-Related Parameters: An Exploratory
Case-Control Study"* (`ARTICLE__revised-v2.docx`). **The submitted manuscript
has not been modified.** Every suggestion below is a recommendation for the
authors to accept, modify or reject.

**One framing point governs everything else.** The computational analysis was
performed *after* the genotyping was complete. It must never be presented as
having guided variant selection, primer design, or the statistical plan. The
honest and defensible framing is: *a retrospective, complementary in silico
analysis that asks whether a plausible molecular mechanism exists for the
variants that were already genotyped*. Any wording that implies the
computational work came first would misrepresent the study.

---

## 1. Decision first: should this go in at all?

There are three viable options. Our recommendation is **Option B**.

| Option | What it means | When to choose it |
|---|---|---|
| **A. Omit** | Publish the manuscript unchanged | If the journal or editor objects to adding a new analysis layer at revision stage, or if the added length is unacceptable |
| **B. Compact integration (recommended)** | One Methods subsection, one short Results subsection, two Discussion paragraphs, one figure, two supplementary tables | Best balance: materially strengthens the mechanistic discussion without changing the study's claims |
| **C. Companion paper** | Publish the computational analysis separately | If the authors want to report the CTCF and histone-mark hypotheses, the ISM reference distribution and the full condition panel in depth |

**Why Option B strengthens the paper.** The manuscript's central weakness is
that it reports a set of null and non-surviving associations without being
able to say *why* a signal might be absent. It currently argues the mechanism
from the ghrelin-GHSR literature (refs 38, 40, 47, 53, 55-57), which
establishes that the pathway matters but says nothing about whether *these
three variants* could plausibly act through it. The computational analysis
closes that gap in a specific, checkable way:

1. It confirms the variants are where the manuscript says they are, and
   independently reproduces the 399 bp amplicon from the published primers.
2. It shows that the model predicts no appreciable regulatory effect for any
   of them, with a local reference distribution demonstrating the model *is*
   responsive at this locus.
3. It identifies a concrete reason the bone arm of the hypothesis is weak:
   *GHSR* is predicted to be essentially unexpressed in osteoblasts and
   osteocytes.

That converts "we found nothing, and we are underpowered" into "we found
nothing, we are underpowered, **and** an orthogonal sequence-based analysis
also predicts little, for a reason we can name." It is a stronger and more
honest paper.

**Important caution.** This also commits the authors to reporting a second
null result. The computational analysis does *not* rescue any of the nominal
associations, and it must not be presented as if it did.

---

## 2. Required corrections, independent of whether the analysis is included

These are errors or imprecisions in the current manuscript that Phase I
identified. They should be fixed **regardless** of the decision above.

### 2.1 rs34136054 allele notation is inconsistent with GRCh38 (must fix)

The manuscript reports rs34136054 as **"C>T"**, with C as wild type. GRCh38
carries **T** at chr3:172,449,372. dbSNP build 157 lists the site as
tri-allelic with alternates C and G, and gnomAD shows that the allele actually
segregating in populations is **G** (global AF 0.005), with C essentially
absent.

The reported notation therefore inverts reference and alternate, and does not
identify which alternate allele was observed. The genotype *identity* is
unaffected - a C/T heterozygote is the same individual as a T/C heterozygote -
so **no result changes**. But as written the notation is incorrect.

*Suggested correction (Methods, "Genetic Analysis"):*

> The three target variants located on chromosome 3 (GRCh38/hg38),
> rs2922126 (chr3:172,449,471 T>A), rs115459563 (chr3:172,449,444 G>A) and
> rs34136054 (chr3:172,449,372 T>C/T>G), were captured within the amplified
> 399-bp upstream regulatory fragment of *GHSR*. Alleles are given on the
> GRCh38 plus strand. Because *GHSR* is transcribed from the minus strand,
> chromatograms read in the direction of transcription show the reverse
> complement of these bases.

If the original chromatogram for the single rs34136054 heterozygote is still
available, re-reading it to establish whether the alternate allele was C or G
would settle the point definitively; both were scored computationally and
neither produced an appreciable predicted effect, so the computational
conclusion does not depend on the answer.

### 2.2 Strand convention is not stated and is applied inconsistently (should fix)

Reconciled against GRCh38:

* **rs2922126 "A>T"** - correct, written in *GHSR* (minus-strand) orientation.
  This is internally consistent with the sequencing: the manuscript's
  "forward" primer anneals to the plus strand, i.e. it primes synthesis in the
  gene's direction.
* **rs115459563 "G>A"** - correct, but written in **plus-strand** orientation,
  the opposite convention.
* **rs34136054 "C>T"** - matches neither convention (see 2.1).

Three variants, three conventions. Adding one sentence stating that alleles
are reported on the plus strand, and that *GHSR* is on the minus strand,
removes the ambiguity. This matters for any reader attempting to replicate the
genotyping.

### 2.3 "Promoter" is imprecise (minor, optional)

The manuscript describes the amplicon as a "promoter/regulatory fragment"
(Methods, "Genetic Analysis"). The three variants lie **916-1,015 bp
upstream** of the canonical *GHSR* transcription start site (chr3:172,448,456)
and are classified by Ensembl VEP as `upstream_gene_variant` for every *GHSR*
transcript. None overlaps an annotated Ensembl Regulatory Build element; the
nearest is an annotated enhancer (ENSR3_B3983C, 90 bp from rs2922126) and a
CTCF binding site (349 bp from rs2922126).

"Upstream regulatory region" or "proximal upstream region" is accurate;
"promoter" overstates the annotation. This is worth correcting because the
manuscript's mechanistic argument rests partly on the region being
regulatory, and "adjacent to an annotated enhancer" is a better-supported and
still favourable description.

### 2.4 Supportive findings worth adding (optional but recommended)

Two Phase I results strengthen the existing Methods at no cost:

* The published primer pair maps to **chr3:172,449,261-172,449,659**, a product
  of exactly **399 bp**, matching the reported size and containing all three
  variants. This is an independent in silico confirmation of the amplification
  strategy and is worth one sentence.
* Observed allele frequencies are consistent with the Middle Eastern
  population in gnomAD: rs2922126 alternate allele 0.24 observed vs 0.276
  (gnomAD `mid`), rs34136054 0.005 vs 0.017. This supports genotyping
  accuracy.

**One frequency does not fit.** rs115459563 was observed in 7 of 50 controls
(carrier frequency 14%, allele frequency 0.035), against a gnomAD `mid` allele
frequency of **0.0068** and a global frequency of 0.004 - roughly five- to
nine-fold higher than expected. With n = 50 this could be sampling variation
or genuine population structure, but it could also indicate a genotyping
artefact at this position. Since the rs115459563 genotype distribution is the
manuscript's only nominally significant genotype finding (p = 0.012), we
recommend that the authors **re-examine those seven chromatograms** before
publication. If the calls are confirmed, the frequency excess itself is worth
a sentence as a finding; if they are not, the headline genotype result
changes. The computational analysis cannot settle this - it is a wet-lab
question - but it is the single most consequential checkable item we
identified.

---

## 3. Proposed placement, Option B

### 3.1 Methods

Insert a new subsection, **"In silico variant-effect prediction"**, after
"Genetic Analysis" and before "Statistical Analysis". Draft text:
[`docs/METHODS.md`](docs/METHODS.md).

Add to the existing "Statistical Analysis" subsection a short statement that
no inferential statistics were applied to the model predictions and that no
P-values or significance thresholds were derived from them. This pre-empts the
most likely reviewer objection.

### 3.2 Results

Insert a new subsection, **"In silico prediction of regulatory effects"**, as
the **last** Results subsection, after the logistic-regression results. Draft
text: [`docs/RESULTS.md`](docs/RESULTS.md).

Placing it last reinforces the chronology: the computational work followed the
experimental work. Opening it with "To assess whether the genotyped variants
could plausibly affect *GHSR* transcription, we performed a retrospective in
silico analysis" makes the retrospective framing explicit at first mention.

Suggested numbering if the existing scheme is kept:

* **Figure 6** - the integrated figure (see 3.4).
* **Table 17** - integrated evidence summary, or move to supplementary.
* **Supplementary Table S1** - variant annotation (`table_S1_variant_annotation.csv`).
* **Supplementary Table S2** - condition selection and data quality
  (`table_S3_condition_selection_and_qc.csv`).
* **Supplementary Table S3** - ISM reference distribution
  (`table_S5_ism_reference_distribution.csv`).

### 3.3 Discussion

Two insertions:

1. **A new paragraph** in "GHSR Polymorphisms and Cognitive-Group
   Classification", placed after the existing paragraph beginning *"However,
   the present study evaluated genetic variants at the DNA level and did not
   measure GHSR expression..."*. That paragraph currently concedes the
   mechanistic gap; the computational result speaks directly to it. Draft
   text: [`docs/DISCUSSION.md`](docs/DISCUSSION.md), §1.

2. **A paragraph** in "GHSR Polymorphisms, Bone Mineral Density, and Thyroid
   Function", after the existing discussion of ghrelin signalling in
   osteoblasts (refs 55-57). The prediction that *GHSR* is essentially
   unexpressed in osteoblasts and osteocytes is directly relevant to that
   argument and should be stated there, with its limitations. Draft text:
   [`docs/DISCUSSION.md`](docs/DISCUSSION.md), §2.

Also add to **Limitations** a sentence acknowledging that the computational
analysis is a model prediction, not a measurement, and that a sequence model
cannot exclude mechanisms it does not represent. Draft text:
[`docs/DISCUSSION.md`](docs/DISCUSSION.md), §3.

Extend **Future Research Directions** with the specific functional experiments
the computational analysis motivates: [`docs/DISCUSSION.md`](docs/DISCUSSION.md), §4.

### 3.4 Figures

Of the four figures generated, **one** belongs in the main text. We recommend
**Figure 2** (`figure2_predicted_effects`), because the ISM reference
distribution in panel B is the element a reviewer will actually need: it shows
simultaneously that the predicted effects are small and that the model is
responsive at this locus, which is what makes the negative result
interpretable rather than vacuous.

Suggested disposition:

| Figure | Recommendation |
|---|---|
| `figure2_predicted_effects` | **Main text** (new Figure 6) |
| `figure1_ghsr_locus` | Supplementary, or main text if space allows - it is the clearest way to show the amplicon/TSS relationship and would also support the §2.3 correction |
| `figure3_celltype_specificity` | Supplementary |
| `figure4_integrated_evidence` | Supplementary, or omit - useful for a reviewer response, arguably too didactic for a published paper |

Legends: [`docs/FIGURE_LEGENDS.md`](docs/FIGURE_LEGENDS.md).

### 3.5 Abstract

Only if the journal's word limit allows. One sentence, after the genotype
results and before the conclusion:

> A retrospective in silico analysis using AlphaGenome 2 predicted no
> appreciable allelic effect on *GHSR* regulation for any of the three
> variants across 16 cell types and tissues, and predicted *GHSR* to be
> essentially unexpressed in osteoblasts and osteocytes.

Do **not** add it at the cost of the existing power and limitation statements,
which are more important to the paper's integrity.

### 3.6 Title, keywords, authorship

* **Title** - no change needed. The computational analysis is supporting
  evidence, not the study's subject. If the authors want to signal it,
  *"...: An Exploratory Case-Control Study with Complementary In Silico
  Analysis"* is accurate, but it slightly oversells a null secondary analysis.
* **Keywords** - consider adding *in silico variant effect prediction* and
  *gene regulation*.
* **Authorship and acknowledgements** - whoever performed and takes
  responsibility for the computational analysis must meet the journal's
  authorship criteria for it, or be acknowledged. The existing statement that
  language editing was assisted by AI tools should be extended, in whatever
  form the journal requires, to disclose how the computational analysis was
  produced. Decide this explicitly rather than by omission.
* **Data availability** - extend the existing statement to cover the
  computational outputs; suggested wording in
  [`docs/CITATIONS.md`](docs/CITATIONS.md) §5.

---

## 4. Wording to avoid

These phrasings would misrepresent the work. Each has a correct alternative.

| Do not write | Write instead |
|---|---|
| "We used AlphaGenome to select candidate variants" | "We retrospectively evaluated the previously genotyped variants using AlphaGenome 2" |
| "PCR and Sanger sequencing validated the predicted regulatory effects" | "PCR and Sanger sequencing confirmed the genotypes; the predicted regulatory effects remain experimentally untested" |
| "AlphaGenome confirmed that these variants are not functional" | "AlphaGenome 2 predicted no appreciable regulatory effect; this does not exclude mechanisms the model does not capture" |
| "The variants disrupt a CTCF binding site" | "rs2922126 was predicted to reduce CTCF ChIP signal in all 12 observed CTCF tracks examined; no motif analysis was performed and no binding-site disruption is claimed" |
| "significantly altered expression" / "P < 0.05" for any prediction | "predicted |log2 fold change| of 0.027 (~1.9% change in predicted signal)"; no P-value exists |
| "patients with Alzheimer's disease" | "participants in the high 6-CIT score group" (the manuscript is already careful about this; keep it) |
| "cell-type-specific regulatory effects were observed" | "effect magnitude tracked predicted *GHSR* expression level across conditions; no variant showed a consistent direction of effect" |
| "functional validation" (for anything done so far) | "functional validation remains to be performed" |

---

## 5. Anticipated reviewer questions

| Likely question | Where it is answered |
|---|---|
| "Why should I believe a near-zero prediction isn't just a failed analysis?" | ISM reference distribution: the same amplicon contains substitutions with 4.7-95x larger predicted effects (Figure 2B; `table_S5`) |
| "Are the predictions robust to sequence context?" | Repeated at 16 kb, 100 kb, 500 kb and 1 Mb (`table_S8`) |
| "How much of this rests on imputed data?" | Every summary recomputed on observed tracks only (`table_S6`); 43,668 of the 43,728 TF ChIP values for the three variants are imputed-exploratory, and this is stated wherever TF results appear |
| "Did you cherry-pick modalities?" | The complete recommended scorer set (22 scorers) was requested for every variant and condition; primary/secondary classification was fixed in `config.py` before results were seen |
| "Why is osteoclast missing?" | No osteoclast condition exists in the model's condition set; blood monocyte is a declared surrogate (`table_S4`) |
| "Does this change the study's conclusions?" | No. It supplies a mechanistic reason consistent with them, and this must be stated as such |
| "Was the in silico analysis used to choose the variants?" | No - stated explicitly in Methods, Results and Figure 4 |

---

## 6. What this analysis cannot do for the manuscript

Stated plainly, so it is not overclaimed:

* It cannot rescue the nominal associations with testosterone, calcium or
  osteoporosis status. Those did not survive FDR correction and the
  computational analysis provides no mechanistic support for them.
* It cannot resolve the rs115459563 genotype distribution. That is a wet-lab
  question (see §2.4).
* It cannot substitute for power. The study remains underpowered for its
  genetic questions, and no computational analysis changes that.
* It cannot establish that these variants are biologically inert. It shows
  that one sequence model, at this locus, in these conditions, predicts little.
* It is not clinical or diagnostic evidence, and must not be described as
  bearing on Alzheimer's disease diagnosis.
