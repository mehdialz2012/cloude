# PSMA / FOLH1 manuscript — independent replication and validation

Independent re-analysis of the manuscript *"An externally validated transcription-factor model of
PSMA (FOLH1) expression defines the PSMA-low state and nominates target-switching and
druggable-dependency strategies in metastatic prostate cancer"* (deliverable package v1.1),
prompted by the reviewer comments from Joe.

> This directory is a separate project from the lytic-peptide-actuator work in the repository
> root; it shares only the PI.

## Phase 1 (this directory) — done

Replication from primary data with independent code, then validation tests the manuscript does
not contain.

**Read first:** [`docs/phase1_replication_report.md`](docs/phase1_replication_report.md)

Headline: of 69 catalogued numerical claims, 65 reproduce exactly and 1 does not (GSEA permutation
p-values); the computational work is sound. Six new tests nonetheless overturn the three headline
conclusions, in the direction Joe's review indicates.

![summary](figures/fig_validation_summary.png)

## Phase 2 — not started

Restructuring and re-analysis along the reframed question, plus the reply to Joe. Scoped in §5 of
the report; awaiting go-ahead.

## Layout

```
code/      the replication pipeline (s00-s08), runnable top to bottom
docs/      the Phase 1 report
results/   machine-readable outputs; claim_by_claim_verdicts.csv is the audit trail
figures/   summary figure (PNG + SVG)
```

## Reproducing

Data are re-downloaded from primary sources (TCGA/Xena, cBioPortal datahub, GEO, DepMap figshare,
HPA, UniProt, OmniPath/CollecTRI, Ensembl, ChIP-Atlas); `code/common.py` documents every URL and
file. Nothing from the authors' deliverable package is used as input.

```bash
python -m venv venv && ./venv/bin/pip install numpy scipy pandas scikit-learn statsmodels \
    lifelines decoupler pyreadr matplotlib
# place the downloaded cohorts in ./data, then:
PYTHONPATH=code python code/s00_meth_extract.py
PYTHONPATH=code python code/s01_landscape.py
...                                    # s02, s02b, s03 ... s08 in order
```

`s02_tf_model.py` and `s02b_leakage_nulls.py` are the slow steps (nested cross-validation and the
random-regulon nulls); budget roughly an hour each on 4 cores.

## Note for the authors

A live `ALPHAGENOME_API_KEY` is hard-coded in the delivered notebook
(`execution_trace/psma_folh1_analysis.ipynb`, cell 41), which has already been circulated to
external collaborators. **Revoke and rotate it.**
