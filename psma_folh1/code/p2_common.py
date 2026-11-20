"""Phase 2 — shared machinery for the restructured analysis.

Design follows the reviewer's five-step outline:
  1 discovery in SU2C with FOLH1 as a continuous variable and explicit adjustment for
    biopsy site, tumour composition and neuroendocrine differentiation
  2 independent validation in GSE126078 (patient-level) and Beltran
  3 surfaceome prioritisation on therapeutic suitability, not statistical enrichment
  4 DepMap as exploratory hypothesis generation
  5 TCGA as a secondary primary-versus-metastatic comparison
"""
import os
import numpy as np
import pandas as pd
from scipy import stats as st
from statsmodels.stats.multitest import multipletests

from common import DATA, OUT, AR20, NE12, LIVER10, zrow, sig_score, load_su2c, load_beltran

# ---------------------------------------------------------------- composition model
# PanglaoDB cell-type markers (OmniPath), grouped into the compartments that confound
# a bulk metastatic biopsy.
COMPARTMENTS = {
    "immune": ["T cells", "Macrophages", "B cells", "NK cells", "Monocytes",
               "Dendritic cells", "Mast cells", "Neutrophils", "Plasma cells"],
    "stromal": ["Fibroblasts", "Endothelial cells", "Smooth muscle cells", "Adipocytes",
                "Pericytes"],
    "hepatocyte": ["Hepatocytes"],
    "osteoblast": ["Osteoblasts"],
    "epithelial": ["Epithelial cells"],
}


def panglao_net(min_specificity=0.0):
    p = pd.read_csv(f"{DATA}/panglaodb.csv")
    p = p[p["human"] == True]
    if min_specificity:
        p = p[p["human_specificity"] <= min_specificity]
    rows = []
    for comp, cts in COMPARTMENTS.items():
        g = sorted(set(p.loc[p["cell_type"].isin(cts), "genesymbol"]))
        rows.append(pd.DataFrame({"source": comp, "target": g, "weight": 1.0}))
    return pd.concat(rows, ignore_index=True)


def composition_scores(logmat, net=None):
    """Per-sample compartment scores (gene-wise z, then mean over markers present)."""
    net = panglao_net() if net is None else net
    out = {}
    z = zrow(logmat)
    for comp, sub in net.groupby("source"):
        gs = [g for g in sub["target"] if g in logmat.index]
        out[comp] = z.loc[gs].mean(axis=0)
    df = pd.DataFrame(out)
    df["AR"] = sig_score(logmat, AR20)[0]
    df["NE"] = sig_score(logmat, NE12)[0]
    df["hepatocyte_10g"] = zrow(logmat.loc[[g for g in LIVER10 if g in logmat.index]]).mean(axis=0)
    return df


# ---------------------------------------------------------------- fast multi-gene OLS
def ols_beta(Y, X, coef=1):
    """OLS of every row of Y (genes x samples) on design X (samples x p).

    Returns beta, se, t, p for column `coef` of X.
    """
    X = np.asarray(X, float)
    Y = np.asarray(Y, float)
    n, p = X.shape
    XtX_inv = np.linalg.pinv(X.T @ X)
    B = XtX_inv @ X.T @ Y.T                      # p x genes
    resid = Y.T - X @ B                          # n x genes
    dof = n - np.linalg.matrix_rank(X)
    sigma2 = (resid ** 2).sum(axis=0) / dof
    se = np.sqrt(sigma2 * XtX_inv[coef, coef])
    beta = B[coef]
    t = beta / np.where(se == 0, np.nan, se)
    pv = 2 * st.t.sf(np.abs(t), dof)
    return beta, se, t, pv


def design(df, cols, add_const=True, dummies=None):
    """Build a design matrix with the FOLH1 term first (after the intercept)."""
    parts = []
    if add_const:
        parts.append(pd.Series(1.0, index=df.index, name="const"))
    for c in cols:
        parts.append(df[c].astype(float).rename(c))
    if dummies is not None:
        for name, s in dummies.items():
            d = pd.get_dummies(s.reindex(df.index).astype(str), prefix=name, drop_first=True)
            for c in d.columns:
                parts.append(d[c].astype(float))
    X = pd.concat(parts, axis=1)
    return X


def fdr(p):
    p = np.asarray(p, float)
    out = np.full_like(p, np.nan)
    ok = np.isfinite(p)
    out[ok] = multipletests(p[ok], method="fdr_bh")[1]
    return out


# ---------------------------------------------------------------- cohort loaders
def load_gse_patient():
    """GSE126078 patient samples only, with the patient identifier parsed from the title."""
    from common import load_gse126078
    g = load_gse126078()
    is_model = g.columns.str.contains("LuCaP|C42B|PC3|PacMet")
    g = g.loc[:, ~is_model]
    pat = pd.Series([c.split("_", 1)[1][:6] for c in g.columns], index=g.columns, name="patient")
    return g, pat


def load_beltran_annot():
    log, cs = load_beltran()
    m = cs.set_index("OTHER_SAMPLE_ID")
    idx = log.columns.intersection(m.index)
    ann = pd.DataFrame({
        "path": m.loc[idx, "CANCER_TYPE_DETAILED"],
        "patient": m.loc[idx, "PATIENT_ID"] if "PATIENT_ID" in m.columns else idx},
        index=idx)
    return log[idx], ann
