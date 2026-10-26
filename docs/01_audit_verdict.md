# Deliverable 1 — §4 Adversarial Audit Verdict

**Project:** Conditionally-activated lytic peptide actuator for prostate cancer
**PI:** Mehdi Alizadeh
**Audit date:** 26 October 2026
**Status:** Delivered before any production compute, per §4.

---

## 0. Headline verdict

**Proceed — but with three substantive changes, and with the pH gate promoted to a
first, cheap, go/no-go experiment that decides the architecture.**

The project is not killed. The single most important risk the brief flagged — that the
2001/2009 PSMA-activated lytic peptide line died for an undisclosed scientific reason —
**resolves in the project's favour**. The three changes I am recommending are:

1. **The novelty claim as written in §5 is no longer accurate and must be narrowed.**
   Two papers (one from 2026, one from 2025) occupy more of the claimed ground than the
   brief's search found. The project is still novel, but not for the reason stated. See §3.
2. **C9's design goal is internally inconsistent with C6 and should be inverted.** The
   mask should block *assembly* while *permitting* pH-gated membrane binding. As written
   ("block both membrane binding and inter-monomer assembly") the mask destroys the
   spatial-focusing effect that C6 identifies as the most interesting idea in the brief.
   This is the most consequential technical correction in this audit. See §6.
3. **The AND-vs-single-gate question (C3) cannot be answered from the literature, but it
   collapses to a single measurable quantity**: the clean selectivity window of the pH
   gate. If the pH gate delivers ≥10×, AND logic wins decisively. If it delivers ≤3×, one
   excellent enzyme gate beats two mediocre ones and the architecture should change. This
   is cheap to determine and should be Work Package 1. See §5.

I am **not** recommending a change of target or of warhead. Both survive the audit, though
PSMA survives for a different reason than C7 gives, and the brief's stated reason is partly
wrong. See §4.

---

## 1. The forensic question: why did the 2001/2009 line stop?

The brief correctly identifies this as materially risk-changing and asks for it to be
settled early. **It is settled, and the answer is institutional, not scientific.**

Both papers verify exactly as described (PubMed):

- **Warren P, Li L, Song W, Holle E, Wei Y, Wagner T, Yu X.** *Cancer Res*
  2001;61(18):6783–7. PMID **11559551**. Amoebapore H3 inactivated by two γ-linked Glu on
  the ε-amino group of the C-terminal Lys; little activity on PSMA-negative PC-3; strong
  concentration-dependent lysis of PSMA-positive LNCaP; blocked by the carboxypeptidase
  inhibitor 4,4′-phosphonicobis(butane-1,3-dicarboxylic acid); active on PSMA-transfected
  PC-3. Affiliation: Oncology Research Institute of the Greenville Hospital System.
- **Holle L, Song W, Holle E, Nilsson J, Wei Y, Li J, Wagner TE, Yu X.** *Mol Med Rep*
  2009;2(3):399–403. PMID **21475841**. DOI
  [10.3892/mmr_00000112](https://doi.org/10.3892/mmr_00000112). LNCaP xenografts:
  significantly decreased tumour size and PSA, **5/12 treated mice tumour-free**, tumour
  peptide levels steady ~6 h, **no toxicity at 30 mg/kg IP or IV** (single-dose).

**Finding: the trail ends because the principal investigator left and changed field, in the
same year as the last paper.** Thomas E. Wagner directed the Oncology Research Institute at
Greenville Hospital System until 2009 and founded Orbis Health Solutions in 2009, pivoting
to an autologous tumour-lysate particle (TLPO) cancer vaccine platform, which is what he has
worked on since. Greenville Health System itself was subsequently absorbed into Prisma
Health. A targeted author/subject search returns no further work on this chemistry by any of
Warren, Li, Song, Holle (L or E), Wei, Wagner or Yu; every other apparent hit is a different
person of the same surname (Yu Xu at Huazhong; Yu Xin at Southwest Medical; Yu Xiaolin at
Augusta). There is **no retraction, no correction, no published failure to replicate, and no
contradicting report.**

**Interpretation.** A programme that had just produced 5/12 complete responses with no
single-dose toxicity does not get dropped for scientific reasons. It was orphaned. The
secondary signal the brief noticed — the venue falling from *Cancer Research* (2001) to
*Molecular Medicine Reports* (2009) — is consistent with a lab winding down and clearing its
drawer rather than with a result that had become hard to defend; note the 2009 paper reports
*positive* in vivo data, which is not what a struggling programme publishes.

**Consequence for risk:** the go/no-go gate in §6 ("the forensic review finds no scientific
reason the prior line was abandoned") is **PASSED**. This is the single best piece of news in
this audit: the concept has positive, unrefuted, never-followed-up in vivo precedent.

**Two cautions that come with it, both of which belong in the risk register:**

- **The 2009 safety result under-models human risk.** Mouse GCPII is *not* expressed in
  mouse prostate, and its tissue distribution differs materially from human (highest in mouse
  kidney, brain, salivary gland) — Knedlík et al. 2017, *FEBS Open Bio*; Rovenská et al.
  2008, *Prostate*. A xenograft model supplies human PSMA only inside the implanted tumour.
  "No toxicity at 30 mg/kg in nude mice" therefore says very little about what a
  PSMA-unmasked lytic peptide does to human salivary gland, kidney and jejunum. Do not carry
  that number forward as a safety expectation.
- **Single-dose only.** The 2009 toxicity study is explicitly single-dose. Nothing is known
  about repeat dosing, which is where a surface-enzyme-activated lytic agent would be expected
  to accumulate injury in PSMA-positive normal tissue.

---

## 2. Claim-by-claim verdicts (C1–C9)

| Claim | Verdict | Basis |
|---|---|---|
| **C1** Conditional activation beats intrinsic selectivity | **Verified, strongly** | Confirmed almost verbatim in the current literature — see below. |
| **C2** A single activating condition leaks | **Verified, and now quantified** | ~4× margin only; see below. |
| **C3** Two gates multiply selectivity | **Unresolved — reduces to one measurable quantity** | See §5. |
| **C4** Gate oligomerization, not net charge | **Verified, strongly** | Melittin pore is cooperative; monomer is not stably lytic. See below. |
| **C5** Extracellular acidity is usable | **Plausible but the weakest empirical link** | See §5. |
| **C6** Surface enzyme beats secreted | **Mechanistically supported, but in tension with C9** | See §6. |
| **C7** PSMA is the best surface enzyme | **Conclusion survives; stated reasoning is partly wrong** | See §4. |
| **C8** The two gates are independent | **Verified with one real exception** | Inflamed bowel. See below. |
| **C9** Poly-γ-Glu is the right mask | **Chemistry right, design goal wrong** | See §6. |

### C1 — verified, essentially verbatim

The brief asked me to verify that the field's stated unmet need is exposure control rather
than potency. It is, and the wording is closer than the brief claims. Xiong et al., *Int J
Nanomedicine* 2026;21:599827, DOI [10.2147/IJN.S599827](https://doi.org/10.2147/IJN.S599827)
(PMC13222608):

> "the central question in melittin nanomedicine is no longer whether the peptide is
> bioactive, but how that activity can be spatially, temporally, and pharmacologically
> constrained in vivo"

and, naming the quantity directly:

> "the most decisive issues are not efficacy alone, but the extent to which toxicity
> attenuation, **exposure control**, PK/PD interpretability, and formulation tractability can
> be achieved within the same platform."

**An important corollary the brief did not draw.** Every current answer to this problem in the
literature is a *formulation* answer — liposomes, polymeric and inorganic nanocarriers,
Pluronic F-127 (Król et al., *Front Oncol* 2026;16:1823622, DOI
[10.3389/fonc.2026.1823622](https://doi.org/10.3389/fonc.2026.1823622), which reduces
haemolysis by co-formulation; and the nanocarrier review above, which reports haemolysis
falling from ~40–60% to <10% by delivery engineering). The 2026 reviews conclude melittin is
"most plausibly advanced in delivery-controlled settings, particularly localized,
route-constrained" applications. **Nobody is solving exposure control with an intrinsic
molecular switch on the peptide itself.** That is a real and defensible differentiator, and it
independently justifies the §2 constraint that the output be SPPS-synthesizable with no
nanoparticle formulation. This is a stronger novelty position than the one the brief wrote,
and it survives the prior-art problems in §3.

### C2 — verified, and the leak is now a number

The brief's cited source is real: Zhao H. et al., "The development of activatable lytic
peptides for targeting triple negative breast cancer," *Cell Death Discovery* 2017;3:17037 —
construct `Ac-w-βA-e8-XPLG*LAG-klUklUkklUklUk-NH2` (e8 = octa-glutamate polyanionic mask,
MMP-cleavable linker, D-amino-acid lytic core, Aib residues). The design is exactly as
described: single-protease-gated, charge-masked.

**A better quantification exists and the brief should use it instead.** Aguilera et al.,
*Integr Biol* 2009;1:371–81, measured the thing C3 actually needs: systemically administered
nonaarginine CPP "causes acute toxicity in mice at a dose 4-fold lower than the MMP cleavable
ACPP, a complication **not observed with an uncleavable ACPP**." That is the leak, isolated and
measured by the correct control:

- A single protease gate buys roughly **4×** of acute-toxicity margin.
- The *uncleavable* control shows **no** toxicity — proving the residual toxicity of the
  cleavable construct is caused by premature/nonspecific cleavage, not by the masked species
  being intrinsically toxic.

This is the cleanest available evidence that single-gate masking is mechanistically sound but
quantitatively insufficient, and it is the number C3 should be argued against. **4× is not a
therapeutic window for a warhead whose dose-limiting toxicity is haemolysis.**

### C4 — verified, strongly; this is the best-supported claim in the brief

C4 rests on melittin killing as an oligomer, with the membrane-bound monomer comparatively
harmless. The biophysics supports this more cleanly than the brief expected, and supplies the
cooperativity number the design needs.

- **Functional pores are supported by ~4–8 peptides.** Pino-Ángeles & Lazaridis, *Biophys J*
  2018 — "transmembrane toroidal pores supported by four to eight peptides." Park et al.,
  *BBRC* 2006 — oligomeric state varies between tetramer and octamer during toroidal pore
  formation (7–8 nm rings around 3.5–4.5 nm cavities). Single-channel recordings: melittin
  averages **3–9 monomers** per pore (Fennouri et al., *BBA Biomembranes* 2017).
- **The decisive result for C4** is Miyazaki & Shinoda, *BBA Biomembranes* 2019, free-energy
  analysis: **six** melittin peptides significantly stabilise a pore; **with five or fewer the
  formation barrier is much higher**; and critically, *"the free energy cost required to insert
  a single melittin into the membrane is too high to explain the one-by-one insertion
  mechanism for pore formation, which also supports the collective melittin mechanism."* A
  single melittin leaves the pore.

**Verdict: C4 holds. Monomeric melittin is not the lytic unit; poration is cooperative with a
threshold near ~6.** Assembly is therefore a genuinely deep lever, and a steep one — a
threshold-6 cooperative process gives a sharp switch, which is exactly what a gate wants. This
also retroactively supports C6: a process with ~6th-order dependence on local surface
concentration is extremely sensitive to whether unmasking is spatially focused or diffuse.

**But C4's support comes with the brief's own §7 Q1 as the price.** The same literature shows
melittin's pore is *dynamic*: monomers "frequently transition between transmembrane and
surface states" (Leveritt et al., *Biophys J* 2015), pore molecularity varies as monomers
diffuse freely (Pino-Ángeles 2018), and melittin permeabilises via large *transient* pores
(Ulmschneider et al., *Nat Commun* 2024). So the gate cannot merely prevent a pre-formed
oligomer from existing — it must suppress *dynamic recruitment* of surface-bound monomers
into a transient assembly. "Membrane-bound but non-lytic until unmasked" is achievable in
principle but the mask must defeat a continuously re-sampling equilibrium, not a static state.
**This is the hardest requirement in the project and I agree with the brief that it is Q1.**

One further mechanistic hazard the brief should know about, because it bears directly on the
§6 library design: Zhao Z-H et al. (bioRxiv 2025) identify **K7** as the residue that decides
whether a stable pore forms — when the K7 side chain ion-pairs with a lipid phosphate, melittin
is anchored flat to the surface and pore formation is *prevented*; K7A and K7Q mutants lose
that anchoring and porate from more starting configurations. **This means a K7→His substitution
is not a neutral member of the His-scan.** At pH 7.4 a neutral His7 mimics K7Q/K7A — i.e. it may
*remove* the anchoring that normally keeps melittin surface-bound and non-porating, which is
the opposite of what gate 1 wants. K7 should be treated as a special case in the scan, not
pooled with the C-terminal cationic cluster.

### C8 — verified, with exactly one real exception, and the brief guessed it correctly

The brief asked whether any normal or pathological tissue is both acidic and PSMA-active, and
named inflamed gut as the obvious suspect. **It is, and it is confirmed.** Vornov et al.,
*Neurochem Res* 2019: in inflammatory bowel disease, *"expression of the FOLH1 gene that codes
for GCPII is strongly upregulated, as is the enzymatic activity in diseased patient
biopsies."* Inflamed tissue is also acidic. Inflamed bowel is therefore simultaneously
PSMA-active and low-pH, and **gate independence fails there.**

Mitigating considerations: it is a disease state rather than normal tissue; it is anatomically
confined; and it is screenable (a patient with active IBD is excludable). **Recommendation:
active IBD becomes a stated exclusion criterion and a named entry in the risk register, not a
project-killer.** This is the one place where C8's "the two gates are genuinely independent"
should be written as "independent except in actively inflamed bowel."

Erythrocytes — the claim that actually matters — are clean: pH 7.4 and no PSMA, so **both**
gates are shut. The structural-gating argument for haemolysis survives intact, and it is the
strongest single argument in the brief.

---

## 3. Prior art: the §5 novelty claim is no longer accurate

Per §2.5 and the instruction not to claim any of this as novel, two papers must be added to
the §5 table and the novelty claim rewritten. **This is the most important correction in this
audit after §6.**

### 3.1 An AND-gated membranolytic peptide already exists (2026)

**Pramanik A. et al., "EGFR-Targeted and MMP-Activated Membranolytic Peptides Kill Cancer
Cells Specifically In Vitro and Reduce Tumor Growth In Vivo," *Mol Pharm* 2026.**

The membranolytic peptide MP1 (from *Polybia paulista*) was modified with an EGFR-binding
sequence **with or without** an MMP-2 cleavage sequence; toxicity correlated with EGFR
expression *and* MMP-2 activity; an MMP-2 inhibitor reduced killing by the cleavable
constructs; efficacy and acceptable systemic tolerability were shown in MDA-MB-468 xenografts.

**The brief's statement "We found no AND-gated lytic peptide" is therefore false as of 2026
and must be struck.** A two-input conditional membranolytic peptide has been built, published,
and taken in vivo.

**What this does and does not take away.** It does not take away as much as it first appears,
because the two inputs are of different kinds:

- Pramanik's first input is a **receptor-homing domain** (EGFR binding) — a *localisation*
  element. Its second is a protease — an *activation* element. The logic is
  "go here, then switch on."
- The brief's two inputs are both **microenvironmental state variables** (pHe and surface
  enzyme activity), and both act on the same step.

That is a real distinction and it is defensible, but it is a narrower claim than "no AND-gated
lytic peptide exists," and it must be argued rather than asserted. I would also note that
Pramanik's design is the closest competitive prior art to this project and the single most
important paper to read in full before writing any grant or paper on this concept.

### 3.2 A computational pipeline for cleavable peptide masks already exists (2025)

**Montaner C. et al., "De Novo Design of Peptide Masks Enables Rapid Generation of
Conditionally-Active Miniprotein Binders," *JACS* 2025.** A general workflow: extend the
C-terminus with a protease-cleavable linker plus a masking helix; ~half of 20 designs achieved
>100× affinity reduction, the best >1000×; cleavage restored binding in 19/20; a
photocleavable variant gave light control.

This bears on the brief's hoped-for methods contribution ("no published pipeline jointly
optimizes post-unmasking aggregation propensity together with lytic potency and a pH-switch
target"). **That specific conjunction does appear to remain open** — Montaner masks a
*binding interface* on a *miniprotein* by steric occlusion, with no aggregation objective and
no pH switch. But the generic idea "computationally design a cleavable mask that reversibly
inactivates a peptide, and validate reversibility" is now published, well-executed, and will be
the reviewer's first comparison. The methods claim must be stated as the *aggregation-gating +
pH-switch* conjunction specifically, and must cite Montaner as the methodological precedent.

### 3.3 A constraint-based objection to part of the §5 table

Several entries in the brief's prior-art table — the (KLAKLAK)₂ / KLA constructs — are
**mechanistically disqualified from being direct prior art for this project**, and the brief
does not notice this. KLA acts on the **mitochondrial** membrane and requires cellular
internalisation to work (see Kotalík et al., *Eur J Med Chem* 2025, review: KLA "disrupts
bacterial membrane as well as eukaryotic mitochondrial membrane," with the field's problem
being "low cell internalization ability"). The brief's own hard constraint #1 is that the
actuator works at the **plasma membrane from the outside**. So the masked-KLA literature is
prior art for *masking chemistry* but not for *this actuator class*, and the project's
differentiation from it is structural rather than incremental. This is worth stating
explicitly because it is a point in the project's favour that the brief left on the table — and
it is also a point in melittin's favour as warhead (§4.2).

### 3.4 Revised novelty claim

I propose replacing the §5 claim with:

> We place the **pore-forming oligomerisation step** of a plasma-membrane-active lytic peptide
> under AND-logic control of two **independent microenvironmental state variables** — extracellular
> pH and cell-surface enzyme activity — such that neither input alone restores lytic assembly;
> we show that the enzyme being *membrane-anchored rather than secreted* spatially focuses the
> unmasked actuator at the surface where it must assemble; and we show that this dual-gating
> structurally abolishes rather than attenuates the erythrocyte liability, because erythrocytes
> fail both gates simultaneously.

Defensible against §5 plus the two additions above, because: the gated step is assembly rather
than binding or uptake (unoccupied); both inputs are microenvironmental rather than one being a
receptor-homing domain (distinguishes Pramanik); the spatial-focusing claim for an anchored vs
secreted enzyme appears genuinely unarticulated (C6 — and this is the part I would build the
paper around); and the actuator acts at the plasma membrane without internalisation or
formulation (distinguishes the KLA and nanocarrier literatures).

**Not novel and must not be claimed:** γ-glutamate masking of a lytic peptide reactivated by
PSMA (Warren 2001 — this is the project's *direct* precedent and must be cited as such in every
document); pH-responsive His-substituted melittin (the 2024 *Asian J Pharm Sci* work);
protease-restored lytic aggregation (US7456146); dual-input conditional membranolytic peptides
in general (Pramanik 2026); computational design of cleavable inactivating masks (Montaner 2025).

---

## 4. Challenges to target and warhead

### 4.1 Is PSMA the right surface enzyme? — Yes, but not for C7's reasons

The brief asked me to challenge C7 hardest. Doing so changes the *reasoning* substantially
while leaving the *conclusion* standing.

**Where C7 is wrong.** C7 hopes that in kidney "the literature reportedly notes no known
specific enzymatic substrate, so expression may not equal activity." **This hope does not
survive.** PSMA/GCPII is enzymatically active in normal human tissue, not merely present:
it hydrolyses folyl-poly-γ-glutamate in jejunum and NAAG in brain, with activity tracking
expression (Šácha et al., *Neuroscience* 2007; Barinka et al., *Protein Sci* 2004 — consensus
normal sites: prostate, kidney, nervous system, small intestine). The clinical consequence is
visible: salivary gland and kidney uptake is the dose-limiting problem for ¹⁷⁷Lu-PSMA
radioligand therapy, serious enough that inhibitor-prodrug strategies are being developed
specifically to protect those organs (Vornov 2019). Whether that uptake is PSMA itself or
partly the homolog GCPIII is actively disputed (Lucaroni et al., *EJNMMI* 2022, proposing
GCPIII; Lee et al., *EJNMMI* 2023, rebutting on 720× selectivity grounds and on GCPIII's
testis/muscle distribution) — but for our purposes the dispute is irrelevant: **either way,
PSMA-directed γ-glutamate chemistry will meet active enzyme in salivary gland, kidney and
small intestine.**

**Why the conclusion nonetheless stands.** Three reasons, of which the brief only half-states
the decisive one:

1. **Those tissues are all at pH 7.4.** The PSMA-active normal compartment and the acidic
   compartment are disjoint (sole exception: inflamed bowel, §2/C8). So C7's liabilities are
   not arguments against PSMA — they are *arguments for the second gate*. The honest framing
   is: a single PSMA gate has a therapeutic-index ceiling set by salivary/renal/jejunal
   activity, which is precisely why the 2001/2009 design could not have been pushed to a
   clean window, and precisely what gate 1 exists to rescue. This inverts C7's logic in a way
   that strengthens the overall architecture.
2. **Patient selectability, which no competitor target has.** PSMA-PET lets you confirm the
   target is present, and at what intensity, *in the individual patient before dosing*;
   ¹⁷⁷Lu-PSMA-617 has already validated that target engagement translates to clinical benefit.
   For a conditionally-activated agent whose entire thesis is "the gate is open only in the
   tumour," the ability to *image the gate* is worth more than a marginally cleaner
   normal-tissue profile. I consider this the strongest argument for PSMA and the brief
   under-weights it.
3. **The androgen-deprivation synergy is real and is a genuine clinical-design advantage.**
   PSMA is inversely regulated by androgen, so patients on standard ADT present *more* target.
   This is well documented and I did not find contradicting evidence, though the effect is
   heterogeneous in CRPC and should not be quantified from in vitro numbers.

**The alternative I was asked to evaluate properly: testisin (PRSS21).** I take it seriously
and it is the best "one excellent gate" candidate available. Founding data verify the brief's
impression (Shigemasa et al., *J Soc Gynecol Investig* 2000;7(6):358–62, PMID **11111072**):
testisin transcript abundant in ovarian carcinoma, **not detected in normal ovary**, and among
all normal fetal and adult tissues examined **only adult testis** showed abundant transcript;
expression higher in advanced stage; GPI-anchored (hence genuinely cell-surface, satisfying
constraint #1 and C6).

Its normal-tissue profile is **cleaner than PSMA's** — and testis sits behind the blood–testis
barrier, so the accessible somatic expression is arguably near-zero, which is a better starting
point than salivary gland plus kidney plus jejunum. That is a real argument and I do not want
to dismiss it the way the brief half-expected me to.

**But it loses on four counts that together are decisive:**

- **No companion imaging, so no patient selection.** There is no testisin-PET. You cannot
  confirm the gate exists in the patient you are about to dose. For this architecture that is
  a severe loss.
- **Wrong enzyme class for C9's elegance.** Testisin is a serine *endo*peptidase. The mask
  becomes a conventional cleavable-linker + polyanion (i.e. the Zhao 2017 / Pramanik 2026
  design space, which is crowded), and you lose the feature that makes C9 attractive: that the
  mask *is* the enzyme's natural substrate, processed residue-by-residue by an exopeptidase.
- **Much thinner enzymology.** PSMA has extracellular-domain crystal structures, a mapped
  glutamate-specific S1′ pocket, characterised kinetics on poly-γ-glutamated folates and NAAG,
  and known inhibitors for use as negative controls. Testisin has little of this, and the §6
  enzyme-accessibility work would have nothing to dock into.
- **It costs the indication, and with it the PI's institutional advantage.** Ovarian rather
  than prostate: no urology collaborators, no existing PSMA/PET familiarity, no sample access.
  §1 sets a high bar for switching indications and testisin does not clear it.

**Verdict: keep PSMA. Document testisin as the designated fallback** — specifically, as the
target to pivot to *if* the pH gate fails (§5), because testisin is exactly the "one excellent
gate" that a collapsed-to-single-gate architecture would want. That makes it a hedge rather
than a rejected option.

### 4.2 Is melittin the right warhead? — Yes, and the audit strengthened the case

Melittin survives, partly for reasons the brief did not give.

**In favour.** It acts directly on the plasma membrane from the outside, satisfying constraint
#1 without internalisation — unlike KLA/(KLAKLAK)₂, which needs cell entry to reach
mitochondria (§3.3). Its poration is steeply cooperative with a threshold near six (§2/C4),
which is the sharp nonlinearity an assembly gate needs. Its reference data are abundant
(structures, pore stoichiometry, mechanism, a large analogue literature with paired potency and
haemolysis values). Potency is not the limiting factor — reported IC₅₀ ~0.8–2.0 µM in
aggressive breast lines (Xiong 2026 review) — which matters for §5's efficiency argument,
because it means the design can afford to lose some potency to gating.

**Against, and this is the real cost.** Its intrinsic membrane avidity is extreme and its pore
is a dynamic equilibrium (§2/C4), making §7 Q1 — bind without lysing — the hardest requirement
in the project. I agree with the brief's own assessment here and I am not able to make this
risk go away by choice of scaffold.

**One comparator I recommend adding, cheaply.** **MelP5** (the five-substitution melittin
variant) is *more* cooperative, not less: single-channel recordings give **10–12 monomers** per
pore versus 3–9 for melittin (Fennouri 2017), and simulations put MelP5 at pentameric vs
melittin tetrameric with higher water permeability (Li et al. 2023). A higher-order assembly
requirement means a **steeper switch** and therefore a *better gatable* warhead in principle —
the mask has more cooperativity to defeat and the on/off ratio should be sharper. The
counterweight is that MelP5 is active at ~10× lower concentration, so the consequences of a
leak are correspondingly worse. This is a genuine design trade-off that costs almost nothing to
evaluate in silico alongside melittin, and I would add it as a single comparator arm rather
than switching to it.

**Not recommended:** MP1 (now occupied by Pramanik 2026 in an adjacent design, and the PI has
no expertise), and any KLA-family scaffold (wrong compartment, crowded).

---

## 5. Is AND logic worth its cost? — The question reduces to one measurable number

The brief asks for C3 quantified and says either answer is acceptable. Here is the honest
analysis, and it does not resolve to a clean "yes."

**What each gate can plausibly deliver:**

- **Enzyme gate alone.** The best available measurement of a single-protease gate's benefit is
  ~**4×** acute-toxicity margin, with the residual toxicity proven to come from premature
  cleavage (Aguilera 2009, §2/C2). For PSMA specifically the ceiling is set by active enzyme in
  salivary gland, kidney and jejunum (§4.1) — generously, perhaps one to 1.5 orders of
  magnitude, degraded by nonspecific and systemic cleavage.
- **pH gate alone.** This is the number nobody has, and it is the weakest link. The physical
  constraint is unforgiving: tumour pHe 6.4–6.9 vs 7.4 systemic is **at most ~1.0 pH unit, and
  realistically ~0.6–0.9**. A single histidine with an apparent pKa sited optimally in the
  window gives a protonation-state ratio across 0.9 units of only a few-fold; you need either
  multiple histidines acting cooperatively, or a steeply cooperative downstream response (which
  C4's threshold-6 assembly conveniently supplies), to convert that into a useful binding
  window. Tumour acidity is also heterogeneous and confined to poorly perfused regions —
  which is a *coverage* problem, not just a contrast problem.

**Why AND logic is nonetheless attractive here — the argument the brief should make.** The
value is not primarily the multiplied number. It is that **the three dominant toxicity
compartments each defeat a *different* single gate**, so only the conjunction closes all three:

| Compartment | pH | PSMA | Defeats a pH-only gate? | Defeats a PSMA-only gate? | Defeats AND? |
|---|---|---|---|---|---|
| Erythrocytes / blood | 7.4 | absent | no | no | **no** |
| Acidic normal tissue (ischaemic muscle, gastric mucosa, inflamed soft tissue) | ~6.5 | absent | **yes** | no | **no** |
| PSMA⁺ normal tissue (salivary, kidney, jejunum) | 7.4 | **active** | no | **yes** | **no** |
| Actively inflamed bowel | ~6.5 | **active** | yes | yes | **yes — the one hole (§2/C8)** |

A single PSMA gate — i.e. the 2001/2009 design — is defeated by salivary gland, kidney and
jejunum. A single pH gate is defeated by every acidic normal niche. **Neither single gate
closes more than one of the three compartments; the conjunction closes all three but one
disease state.** That, rather than a multiplication of selectivity factors, is the defensible
case for AND logic, and I think it is a strong one.

**The cost, stated honestly.** Serial conditional steps multiply yield losses. If the pH gate
leaves only a fraction *f_pH* of actuator competent to bind at tumour pH, and PSMA processes
only a fraction *f_enz* to completion within the residence time (§7 Q2 — processivity on a
bulky cargo is genuinely unknown), the active fraction is the product. A plausible 0.5 × 0.5
costs 4× potency. Melittin can afford roughly that (sub-µM IC₅₀); it could not afford 20×.

**So the decision rule, which I recommend adopting explicitly:**

- If the pH gate delivers a **clean ≥10× window** (transition midpoint 6.6–6.8, <10%
  association at 7.4), **AND logic wins decisively** — ~10× × ~10–30× at a ~4× potency cost is
  a large net gain over any single gate, and it closes all three compartments.
- If the pH gate delivers **≤3×**, the arithmetic inverts: ~3× × ~10× at a 4× potency cost is
  roughly a 7× effective gain, which is **worse than one excellent enzyme gate at full
  potency** — and in that case the right move is a single-gate design against the cleanest
  available enzyme, which is where testisin (§4.1) comes back as the fallback.

**This is why I am recommending the pH gate be promoted to Work Package 1.** It is the cheapest
of the §6 work packages, it is the one with a published validation target (reproduce the 2024
His-melittin analogues and a pHLIP WT/Var3 benchmark to ≤0.5 pKa unit, per §6), and its answer
determines the architecture. Determining it before spending quota on the assembly free-energy
work — which is the expensive part — is the correct ordering. The brief's §6 has it second;
I am moving it first and making it a hard gate.

---

## 6. The architectural correction: C9 as written defeats C6

This is the most consequential technical finding in the audit, and it is an internal
inconsistency in the brief rather than a literature problem.

**C6** says a *membrane-anchored* enzyme is better than a secreted one because unmasking
happens at the membrane, so the liberated actuator is "locally concentrated exactly where it
must assemble." C4 supplies the reason this matters enormously: assembly is cooperative with a
threshold near six, so the pore rate depends very steeply (high-order) on *local* surface
concentration. Spatial focusing is therefore not a marginal effect — it is the difference
between crossing the threshold and not. **I think C6 is the best idea in the brief and the
right thing to build the paper around.**

**C9** says the poly-γ-Glu tail should be "anionic and bulky enough to block **both** membrane
binding **and** inter-monomer assembly."

**These cannot both be satisfied.** If the mask blocks membrane binding, the masked actuator
cannot dock at the surface, and therefore cannot present its tail to a membrane-anchored PSMA
on that same surface. Unmasking is then forced into the bulk fluid — which is exactly the
secreted-enzyme scenario C6 argues against, with the product free to diffuse away before it
reaches assembly threshold. **Blocking binding converts the surface enzyme into a functionally
secreted one and throws away the project's most interesting mechanism.**

This also resolves §7 Q7 ("does unmasking need to happen before or after membrane binding, and
does the answer force a different architecture?"). It does force one:

> **Design target (revised): the mask must block *assembly* while *permitting* pH-gated
> monomer binding.**
>
> Productive sequence: at tumour pH the His-protonated cationic face allows the *masked*
> monomer to dock to the membrane (bound but non-lytic, because assembly is blocked) → the
> bound monomer presents its γ-Glu tail to surface PSMA → local, processive cleavage →
> unmasked monomers accumulate *in the plane of the membrane* at high local density → threshold-6
> cooperative assembly → pore.

Every element of the brief's logic is preserved and the inconsistency disappears: gate 1
controls binding (pH), gate 2 controls assembly (PSMA), C6's spatial focusing is retained
because cleavage happens on a membrane-bound substrate, and C4's cooperative threshold is what
converts local focusing into a sharp switch.

**Consequences for the plan, which are substantial:**

1. **The mask design variable changes.** §6 asks whether tail *charge* or tail *length*
   governs the block, and proposes decomposing electrostatic repulsion vs steric exclusion.
   That decomposition becomes the central design question rather than a sub-analysis, because
   the goal is now *differential*: enough inter-monomer repulsion/steric bulk to prevent
   oligomerisation, while not so much net anionic charge that it cancels the pH-gated cationic
   face and abolishes binding. Short-and-bulky may beat long-and-highly-charged. The
   decomposition in §6 is well-posed to answer this and should be kept, with its purpose
   restated.
2. **"State 2" in the truth table becomes the direct test of the revised design goal** —
   acidic, masked: *binding but no oligomer, no pore*. The brief already identifies State 2 as
   decisive. Under the revised architecture it is decisive for a sharper reason: it is the
   state that demonstrates binding and assembly have been successfully *decoupled*, which is
   the entire mechanism. I endorse giving it the most replicates.
3. **It raises the enzyme-accessibility risk rather than lowering it.** §6 already flags that
   PSMA's active site sits at the bottom of a deep funnel in a homodimeric ectodomain and that
   a 26-residue amphipathic helix may not fit. The revised architecture makes this harder, not
   easier: the substrate is now *membrane-bound* and must present its tail into that funnel
   while the core stays docked in the bilayer, which is a severe geometric constraint on the
   linker. This becomes a priority structural question, not an afterthought — and note the 2001
   work succeeded with a *smaller* lytic fragment (amoebapore H3 domain), so melittin's bulk is
   genuinely unsettled territory. If the geometry proves impossible, the fallback is
   cleavage-then-binding, which works but costs C6.
4. **It changes what "n" is optimising for.** The §6 mask-length ladder (n = 2, 4, 6, 8, with
   n=2 reproducing the 2001 prior art as the charge-masking-only control) stays, and the n=2
   control becomes more valuable: it is the direct experimental contrast between
   charge-masking (prior art) and assembly-masking (this project's claim).

---

## 7. Changes I am making to the brief, and open items

### Changes (deliverable 9 — running record of disagreements)

| # | Brief says | I am doing instead | Why |
|---|---|---|---|
| 1 | §5 novelty claim as written | Narrowed claim (§3.4); "no AND-gated lytic peptide" struck | Pramanik 2026 and Montaner 2025 occupy claimed ground |
| 2 | C9: mask blocks binding **and** assembly | Mask blocks **assembly**, permits pH-gated binding | As written it defeats C6 and forces bulk unmasking (§6) |
| 3 | §6 order: pH calibration second | **pH gate is WP1 and a hard architecture gate** | Its value decides AND vs single-gate (§5) |
| 4 | C7: kidney expression may not equal activity | Rejected; PSMA is enzymatically active in normal tissue | Šácha 2007, Vornov 2019, Pluvicto toxicity (§4.1) |
| 5 | C7 rationale: cleanest normal-tissue profile | PSMA kept for **patient selectability** instead | Testisin's normal-tissue profile is cleaner (§4.1) |
| 6 | His-scan over K7, K21, R22, K23, R24 uniformly | **K7 treated as a special case** | K7 anchoring suppresses poration; H7 may mimic K7Q (§2/C4) |
| 7 | Melittin alone as warhead | Melittin primary + **MelP5 as one comparator arm** | Higher-order assembly (10–12) ⇒ steeper, better-gatable switch (§4.2) |
| 8 | C8: gates independent | "Independent **except in actively inflamed bowel**"; IBD = exclusion criterion | Vornov 2019 confirms FOLH1 expression *and activity* up in IBD (§2/C8) |
| 9 | 2009 mouse safety as reassurance | Explicitly **discounted** | Mouse GCPII absent from mouse prostate, different distribution (§1) |

Sequence numbering in `GIGAVLKVLTTGLPALISWIKRKRQQ` was checked and the brief is correct:
G1 I2 G3 A4 V5 L6 **K7** V8 L9 T10 T11 G12 L13 P14 A15 L16 I17 S18 W19 I20 **K21 R22 K23 R24**
Q25 Q26 — 26 residues, C-terminal `KRKRQQ` at 21–26.

### Open items not settled by this audit

- **FTO / patent scan is not done.** US7456146 claims (including dependent claims and family
  members) must be read in full, plus a scan on PSMA-activated prodrugs, γ-glutamate masking of
  cytotoxics, and masked lytic peptides — assignees, claim scope, priority dates, status. The
  2001/2009 work's IP position matters especially: if US7456146 or a sibling is assigned to
  Wagner/Greenville and has lapsed, that is a materially different situation than if it is live
  and assigned. Flagged as the next discrete task, before any publication or grant.
- **§7 Q2 (PSMA processivity on a bulky cargo) is unanswered and unanswerable from the
  literature.** Whether PSMA runs a long γ-Glu tail to completion or stalls partway, leaving an
  ambiguous partially-masked species, is a wet-lab question. It is also a design risk, because a
  stalled intermediate is neither safely masked nor fully active.
- **§7 Q3 (duodenal folate deglutamation vs a systemically dosed poly-γ-Glu drug).** The brief's
  argument that this activity faces the lumen rather than the blood is *plausible and
  consistent* with jejunal folate absorption being the enzyme's physiological role (Šácha 2007),
  but I did not find a source that directly addresses accessibility of intestinal brush-border
  PSMA to a systemically administered substrate. Marked open for the wet lab, per §2.5, rather
  than estimated.
- **§7 Q6 (does the masked monomer stay soluble and inert in plasma?)** A poly-anionic
  amphipathic peptide is not automatically well-behaved; this is a formulation/biophysics
  question the simulations can partially address (masked-state aggregation propensity) but not
  settle.
- **§7 Q8 (death modality / immunogenic cell death).** Agree with the brief's own framing: report
  as mechanistic observation if seen, do not frame as novelty — the 2025–2026 self-assembling
  peptide → pyroptosis literature has claimed that narrative.

---

## 8. Go/no-go status against §6's gates

| Gate | Status |
|---|---|
| Forensic review finds no scientific reason for abandonment | **PASS** (§1) — institutional, PI left in 2009, no retraction or failed replication |
| Published pKa values reproduced to ≤0.5 unit | **Pending — now WP1** (§5) |
| Candidate with midpoint 6.6–6.8 and <10% association at pH 7.4 | Pending — WP1; **this is the architecture decision point** |
| Mask length that both abolishes and restores pore formation | Pending |
| State 2 shows no pore above pre-registered threshold | Pending — the decisive simulation (§6) |
| Enzyme sterically accommodates the construct | Pending — **risk raised** by the §6 correction |

**Bottom line: the concept is alive, its direct precedent worked and was abandoned for
non-scientific reasons, and its novelty is real but narrower than claimed. One internal
inconsistency (C9 vs C6) needed fixing and has been. The architecture now hinges on a single
cheap measurement — the pH gate's window — which should be done first.**

---

## References

Primary sources opened and read for this audit. PubMed was used for retrieval and metadata.

1. Warren P, Li L, Song W, Holle E, Wei Y, Wagner T, Yu X. *Cancer Res* 2001;61(18):6783–7. PMID 11559551.
2. Holle L, Song W, Holle E, Nilsson J, Wei Y, Li J, Wagner TE, Yu X. *Mol Med Rep* 2009;2(3):399–403. DOI 10.3892/mmr_00000112. PMID 21475841.
3. Xiong Y, Cai H, Jia S, Gong B, Hou K, Wu H, Chen H. *Int J Nanomedicine* 2026;21:599827. DOI 10.2147/IJN.S599827.
4. Król G, Mańkowska A, Paprocka P, et al. *Front Oncol* 2026;16:1823622. DOI 10.3389/fonc.2026.1823622.
5. Shigemasa K, Underwood LJ, Beard J, Tanimoto H, Ohama K, Parmley TH, O'Brien TJ. *J Soc Gynecol Investig* 2000;7(6):358–62. PMID 11111072.
6. Miyazaki Y, Shinoda W. *BBA Biomembranes* 2019 — free-energy analysis of melittin pore formation.
7. Pino-Ángeles A, Lazaridis T. *Biophys J* 2018 — peptide charge, orientation, concentration in melittin pores.
8. Leveritt JM, Pino-Ángeles A, Lazaridis T. *Biophys J* 2015 — structure of a melittin-stabilized pore.
9. Fennouri A, et al. *BBA Biomembranes* 2017 — single-channel recordings, melittin vs MelP5.
10. Park S-C, et al. *BBRC* 2006 — toroidal pore and oligomerization by TEM.
11. Ulmschneider JP, et al. *Nat Commun* 2024 — melittin permeabilizes via large transient pores.
12. Zhao Z-H, et al. bioRxiv 2025 — unbiased all-atom MD; K7 as the pore-determining residue.
13. Zhao H, et al. *Cell Death Discovery* 2017;3:17037 — activatable masked lytic peptides, TNBC.
14. Aguilera TA, et al. *Integr Biol* 2009;1:371–81 — ACPP vs CPP systemic distribution and acute toxicity.
15. Pramanik A, et al. *Mol Pharm* 2026 — EGFR-targeted, MMP-activated membranolytic MP1 peptides.
16. Montaner C, et al. *JACS* 2025 — de novo design of peptide masks for conditionally-active miniproteins.
17. Kotalík K, et al. *Eur J Med Chem* 2025 — KLAK peptide review (mitochondrial mechanism).
18. Vornov JJ, et al. *Neurochem Res* 2019 — GCPII inhibitors outside the brain; FOLH1 in IBD.
19. Šácha P, et al. *Neuroscience* 2007;144(4) — GCPII expression and activity in human brain; jejunal folate hydrolase.
20. Barinka C, et al. *Protein Sci* 2004 — N-glycosylation sites required for GCPII proteolytic activity.
21. Knedlík T, et al. *FEBS Open Bio* 2017 — mouse vs human GCPII activity and tissue distribution.
22. Rovenská M, et al. *Prostate* 2008 — tissue expression and enzymology of human, rat, pig PSMA.
23. Lucaroni L, et al. *EJNMMI* 2022 — GCPIII cross-reactivity as cause of salivary/kidney uptake.
24. Lee Z, et al. *EJNMMI* 2023 — rebuttal: GCPIII is not the off-target for urea-based PSMA ligands.
25. Jeong S-H, et al. *Investig Clin Urol* 2021 — PSMA theragnostics; 100–1000× tumour vs normal prostate.
26. Ristau BT, et al. *Urol Oncol* 2013 — 20 years of PSMA research; consensus normal expression sites.
