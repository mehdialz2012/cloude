"""Central configuration for the AlphaGenome 2 GHSR variant-effect study.

All coordinates are GRCh38/hg38, 1-based, plus-strand (reference-genome)
orientation. REF alleles in `VARIANTS` are the GRCh38 reference bases and are
verified against the reference genome in Phase I (`phase1_annotation.py`).

Nothing in this module contains patient-level data: only public variant
coordinates and public ontology identifiers are ever sent to the AlphaGenome
API.
"""

from __future__ import annotations

import dataclasses
import os
import pathlib

# --------------------------------------------------------------------------
# Paths
# --------------------------------------------------------------------------
PKG_DIR = pathlib.Path(__file__).resolve().parent
PROJECT_DIR = PKG_DIR.parents[1]

DATA_RAW = PROJECT_DIR / "data" / "raw"
DATA_PROCESSED = PROJECT_DIR / "data" / "processed"
TABLES = PROJECT_DIR / "results" / "tables"
FIGURES = PROJECT_DIR / "results" / "figures"
CONFIG_DIR = PROJECT_DIR / "config"
DOCS = PROJECT_DIR / "docs"

for _p in (DATA_RAW, DATA_PROCESSED, TABLES, FIGURES, CONFIG_DIR, DOCS):
    _p.mkdir(parents=True, exist_ok=True)

API_KEY_ENV = "ALPHA_GENOME_API_KEY"

# --------------------------------------------------------------------------
# Genome build / locus
# --------------------------------------------------------------------------
GENOME_BUILD = "GRCh38"
ENSEMBL_ASSEMBLY = "GRCh38"
CHROM = "chr3"
ENSEMBL_CHROM = "3"

# GHSR, Ensembl release-independent identifiers (verified in Phase I).
GHSR_GENE_ID = "ENSG00000121853"
GHSR_SYMBOL = "GHSR"
GHSR_CANONICAL_TX = "ENST00000241256"

# Experimentally sequenced amplicon from the source manuscript (Table 1
# primers), mapped to GRCh38 in Phase I.
AMPLICON_START = 172_449_261
AMPLICON_END = 172_449_659

PRIMER_FORWARD = "ACGCGTCTTCTATTTATTGAGTG"
PRIMER_REVERSE = "AGTAAACCCTGTCCAGAATCC"
REPORTED_AMPLICON_BP = 399


# --------------------------------------------------------------------------
# Variants
# --------------------------------------------------------------------------
@dataclasses.dataclass(frozen=True)
class VariantSpec:
    """One REF->ALT configuration to be scored.

    Attributes:
      key: Short unique key used in output tables and file names.
      rsid: dbSNP reference SNP identifier.
      position: 1-based GRCh38 plus-strand position.
      ref: GRCh38 plus-strand reference base.
      alt: Alternate base, plus-strand orientation.
      manuscript_notation: Allele notation exactly as written in the source
        manuscript (not necessarily plus-strand oriented).
      manuscript_label: SNP label used in the source manuscript.
      is_primary: True for the configuration used as the headline result for
        this rsID; False for additional alternate alleles analysed as a
        sensitivity analysis.
      note: Free-text provenance note.
    """

    key: str
    rsid: str
    position: int
    ref: str
    alt: str
    manuscript_notation: str
    manuscript_label: str
    is_primary: bool
    note: str


VARIANTS: tuple[VariantSpec, ...] = (
    VariantSpec(
        key="rs2922126_T_A",
        rsid="rs2922126",
        position=172_449_471,
        ref="T",
        alt="A",
        manuscript_notation="A>T",
        manuscript_label="SNP1",
        is_primary=True,
        note=(
            "Only ALT allele in dbSNP b157. Manuscript 'A>T' is the same "
            "substitution written in GHSR (minus-strand) orientation."
        ),
    ),
    VariantSpec(
        key="rs115459563_G_A",
        rsid="rs115459563",
        position=172_449_444,
        ref="G",
        alt="A",
        manuscript_notation="G>A",
        manuscript_label="SNP2",
        is_primary=True,
        note=(
            "Only ALT allele in dbSNP b157. Manuscript 'G>A' matches the "
            "plus-strand (reference-genome) orientation."
        ),
    ),
    VariantSpec(
        key="rs34136054_T_C",
        rsid="rs34136054",
        position=172_449_372,
        ref="T",
        alt="C",
        manuscript_notation="C>T",
        manuscript_label="SNP3",
        is_primary=True,
        note=(
            "dbSNP b157 ALT #1 of a tri-allelic site (T>C, T>G). The "
            "manuscript's 'C>T' assigns C as wild type, but GRCh38 carries T; "
            "REF/ALT are therefore inverted in the manuscript notation."
        ),
    ),
    VariantSpec(
        key="rs34136054_T_G",
        rsid="rs34136054",
        position=172_449_372,
        ref="T",
        alt="G",
        manuscript_notation="C>T",
        manuscript_label="SNP3",
        is_primary=False,
        note=(
            "dbSNP b157 ALT #2 and the only rs34136054 alternate allele "
            "segregating appreciably in gnomAD v4 (global AF 0.005). Analysed "
            "as a sensitivity analysis because the manuscript's notation does "
            "not identify which alternate allele was observed."
        ),
    ),
)

PRIMARY_VARIANTS = tuple(v for v in VARIANTS if v.is_primary)
RSIDS = ("rs2922126", "rs115459563", "rs34136054")


# --------------------------------------------------------------------------
# Conditions (AlphaGenome 2 reference conditions)
# --------------------------------------------------------------------------
@dataclasses.dataclass(frozen=True)
class ConditionSpec:
    """One AlphaGenome 2 reference condition selected for analysis.

    Attributes:
      curie: Composite ontology CURIE accepted by the v2 API `condition=` arg.
      label: Short display label.
      panel: Biological panel grouping for figures.
      rationale: Why this condition is in the panel, tied to the manuscript.
      target_priority: True if this condition addresses one of the
        pre-specified biological priorities (brain / pituitary / bone);
        False for context / specificity-reference tissues.
      surrogate_for: Non-empty when this condition is used as a defensible
        surrogate for a priority cell type that the model does not support.
    """

    curie: str
    label: str
    panel: str
    rationale: str
    target_priority: bool
    surrogate_for: str = ""


CONDITIONS: tuple[ConditionSpec, ...] = (
    # ---- Brain / neuronal -------------------------------------------------
    ConditionSpec(
        "UBERON:0001898",
        "hypothalamus",
        "Brain / neuronal",
        "Primary site of GHSR-mediated ghrelin signalling and the "
        "pre-specified hypothalamic-neuron priority.",
        True,
    ),
    ConditionSpec(
        "UBERON:0002305",
        "layer of hippocampus",
        "Brain / neuronal",
        "Closest available match to the pre-specified hippocampal-neuron "
        "priority; memory-relevant region.",
        True,
    ),
    ConditionSpec(
        "UBERON:0002421_CL:0002604",
        "hippocampal astrocyte",
        "Brain / neuronal",
        "Second hippocampal condition; non-neuronal hippocampal cell type.",
        True,
    ),
    ConditionSpec(
        "UBERON:0000955",
        "brain",
        "Brain / neuronal",
        "Best-supported brain condition (RNA-seq, DNase, CAGE and TF ChIP all "
        "observed); anchors the brain panel.",
        True,
    ),
    ConditionSpec(
        "UBERON:0001871",
        "temporal lobe",
        "Brain / neuronal",
        "Region in which reduced GHSR1a mRNA was reported in Alzheimer's "
        "disease (Gahete et al. 2010), cited by the manuscript.",
        True,
    ),
    ConditionSpec(
        "UBERON:0000955_CL:0000540",
        "brain neuron",
        "Brain / neuronal",
        "Generic brain neuronal population.",
        True,
    ),
    ConditionSpec(
        "CL:0000679",
        "glutamatergic neuron",
        "Brain / neuronal",
        "Best-supported neuronal condition (RNA-seq, DNase and TF ChIP "
        "observed); excitatory neurons of hippocampus/cortex.",
        True,
    ),
    # ---- Pituitary --------------------------------------------------------
    ConditionSpec(
        "UBERON:0000007",
        "pituitary gland",
        "Pituitary",
        "Canonical high-GHSR tissue (growth-hormone secretagogue receptor); "
        "pre-specified pituitary priority and the natural positive control "
        "for GHSR expression.",
        True,
    ),
    # ---- Bone -------------------------------------------------------------
    ConditionSpec(
        "UBERON:0002481_CL:0000062",
        "osteoblast (bone tissue)",
        "Bone",
        "Pre-specified osteoblast priority; best-supported bone condition "
        "(RNA-seq, CAGE, TF ChIP and histone ChIP observed).",
        True,
    ),
    ConditionSpec(
        "CL:0000137",
        "osteocyte",
        "Bone",
        "Pre-specified osteocyte priority (RNA-seq, DNase and TF ChIP "
        "observed).",
        True,
    ),
    ConditionSpec(
        "CL:0000134",
        "mesenchymal stem cell",
        "Bone",
        "Osteoblast precursor; RNA-seq, DNase and histone ChIP observed.",
        True,
    ),
    ConditionSpec(
        "UBERON:0000178_CL:0000576",
        "blood monocyte",
        "Bone",
        "Osteoclast-lineage surrogate: no osteoclast or pre-osteoclast "
        "condition exists in the AlphaGenome 2 condition set, and osteoclasts "
        "differentiate from circulating monocytes.",
        True,
        surrogate_for="osteoclast / osteoclast precursor",
    ),
    # ---- Context / specificity-reference tissues --------------------------
    ConditionSpec(
        "UBERON:0000945",
        "stomach",
        "Context tissue",
        "Ghrelin-producing tissue; all seven primary output types observed, "
        "so it is also the strongest specificity reference.",
        False,
    ),
    ConditionSpec(
        "UBERON:0000473",
        "testis",
        "Context tissue",
        "Relevant to the manuscript's nominal rs2922126-testosterone "
        "observation (GHSR reported in Leydig cells); fully observed.",
        False,
    ),
    ConditionSpec(
        "UBERON:0002046",
        "thyroid gland",
        "Context tissue",
        "Relevant to the manuscript's thyroid-function analyses; fully "
        "observed.",
        False,
    ),
    ConditionSpec(
        "UBERON:0000992",
        "ovary",
        "Context tissue",
        "Relevant to the manuscript's estradiol / menopause variables; fully "
        "observed.",
        False,
    ),
)

CONDITION_BY_CURIE = {c.curie: c for c in CONDITIONS}

PANEL_ORDER = ("Brain / neuronal", "Pituitary", "Bone", "Context tissue")

# --------------------------------------------------------------------------
# Scoring configuration
# --------------------------------------------------------------------------
# Sequence context used for the headline analysis. 1 Mb is the largest
# supported context and is required for the gene-mask scorers to see GHSR and
# its neighbours.
SEQUENCE_LENGTH_KEY = "SEQUENCE_LENGTH_1MB"

# Sensitivity analysis over sequence context length.
SENSITIVITY_SEQUENCE_LENGTH_KEYS = (
    "SEQUENCE_LENGTH_16KB",
    "SEQUENCE_LENGTH_100KB",
    "SEQUENCE_LENGTH_500KB",
    "SEQUENCE_LENGTH_1MB",
)

# Modalities that are mechanistically plausible a priori for a non-coding
# upstream SNV. Reported as the primary endpoints.
PRIMARY_SCORER_KEYS = (
    "RNA_SEQ",
    "CAGE",
    "DNASE",
    "ATAC",
    "SC_ATAC",
    "CHIP_TF",
    "CHIP_HISTONE",
)

# Modalities computed for completeness / to avoid selective omission, but not
# expected a priori for an upstream intergenic SNV.
SECONDARY_SCORER_KEYS = (
    "SPLICE_SITES",
    "SPLICE_SITE_USAGE",
    "SPLICE_JUNCTIONS",
    "POLYADENYLATION",
    "CONTACT_MAPS",
    "METHYLATION",
    "SC_GENE_EXPRESSION",
)

# `_ACTIVE` scorers return absolute activity (max of REF and ALT), not a
# REF/ALT difference. They are used for QC (is the gene/track active at all?),
# never as effect sizes.
ACTIVITY_SCORER_KEYS = (
    "RNA_SEQ_ACTIVE",
    "CAGE_ACTIVE",
    "DNASE_ACTIVE",
    "ATAC_ACTIVE",
    "SC_ATAC_ACTIVE",
    "CHIP_TF_ACTIVE",
    "CHIP_HISTONE_ACTIVE",
    "SC_GENE_EXPRESSION_ACTIVE",
)

# Genes of interest for the "nearby genes" objective: every protein-coding
# gene whose TSS lies within the 1 Mb scoring window (populated in Phase I).
FOCUS_GENES = ("GHSR", "TNFSF10", "FNDC3B", "NCEH1", "ECT2")


def api_key() -> str:
    """Return the AlphaGenome API key from the environment.

    The key is never written to disk, logged, or embedded in any output.
    """
    key = os.environ.get(API_KEY_ENV)
    if not key:
        raise RuntimeError(
            f"Environment variable {API_KEY_ENV} is not set. Export the "
            "AlphaGenome API key before running the pipeline; never hard-code "
            "it into source files."
        )
    return key
