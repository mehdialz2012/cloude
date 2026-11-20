"""Shared loaders for the independent replication of the PSMA/FOLH1 manuscript.

Data were re-downloaded from primary sources by this replication (not taken from the
authors' deliverables):
  TCGA-PRAD   UCSC Xena HiSeqV2 (log2(norm_count+1)) + PRAD_clinicalMatrix + HM450
  SU2C/PCF    cBioPortal datahub prad_su2c_2019 (polyA FPKM) + clinical
  Beltran     cBioPortal datahub nepc_wcm_2016 (FPKM) + clinical
  GSE126078   GEO per-sample processed FPKM files
  DepMap      24Q4 figshare (CRISPRGeneEffect, OmicsExpressionProteinCodingGenesTPMLogp1, Model)
  HPA         rna_tissue_consensus.tsv
  CollecTRI   decoupler/omnipath
"""
import os
import numpy as np
import pandas as pd

DATA = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "data")
OUT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "out")
os.makedirs(OUT, exist_ok=True)

# Signatures exactly as specified in the manuscript Methods
AR20 = ['KLK3', 'KLK2', 'TMPRSS2', 'FKBP5', 'NKX3-1', 'PLPP1', 'PMEPA1', 'PART1', 'STEAP4',
        'ABCC4', 'ALDH1A3', 'GNMT', 'HERC3', 'INPP4B', 'MAP7', 'MPHOSPH9', 'SLC26A2', 'TPD52',
        'ZBTB10', 'AR']
NE12 = ['SYP', 'CHGA', 'NCAM1', 'ENO2', 'DLL3', 'ASCL1', 'NEUROD1', 'INSM1', 'POU3F2', 'EZH2',
        'AURKA', 'MYCN']
LIVER10 = ['ASGR1', 'LIPC', 'GYS2', 'PNPLA3', 'TFR2', 'ALB', 'APOB', 'HAMP', 'CYP3A4', 'SERPINA1']
STROMAL10 = ['COL1A1', 'COL1A2', 'COL3A1', 'DCN', 'LUM', 'FAP', 'ACTA2', 'TAGLN', 'VIM', 'SPARC']


def zrow(mat):
    """Gene-wise (row-wise) z-score across samples."""
    m = mat.mean(axis=1)
    s = mat.std(axis=1)
    return mat.sub(m, axis=0).div(s.replace(0, np.nan), axis=0)


def sig_score(logmat, genes):
    gs = [g for g in genes if g in logmat.index]
    return zrow(logmat.loc[gs]).mean(axis=0), len(gs)


def load_tcga():
    expr = pd.read_csv(f"{DATA}/tcga_prad_HiSeqV2.gz", sep="\t", index_col=0)
    clin = pd.read_csv(f"{DATA}/tcga_prad_clinicalMatrix", sep="\t", index_col=0, low_memory=False)
    return expr, clin


def load_su2c():
    expr = pd.read_csv(f"{DATA}/su2c_data_mrna_seq_fpkm_polya.txt", sep="\t", index_col=0, comment="#")
    expr = expr[~expr.index.duplicated(keep="first")]
    expr = expr.drop(columns=["Entrez_Gene_Id"], errors="ignore")
    expr = expr.select_dtypes(include=[np.number])
    log = np.log2(expr + 1)
    cs = pd.read_csv(f"{DATA}/su2c_data_clinical_sample.txt", sep="\t", comment="#", index_col=0)
    cp = pd.read_csv(f"{DATA}/su2c_data_clinical_patient.txt", sep="\t", comment="#", index_col=0)
    return log, cs, cp


def load_beltran():
    expr = pd.read_csv(f"{DATA}/nepc_data_mrna_seq_fpkm.txt", sep="\t", index_col=0, comment="#")
    expr = expr[~expr.index.duplicated(keep="first")]
    expr = expr.drop(columns=["Entrez_Gene_Id"], errors="ignore")
    expr = expr.select_dtypes(include=[np.number])
    log = np.log2(expr.clip(lower=0) + 1)   # manuscript: negative-floor artifact clipped at 0
    cs = pd.read_csv(f"{DATA}/nepc_data_clinical_sample.txt", sep="\t", comment="#", index_col=0)
    return log, cs


def load_gse126078():
    import glob
    files = sorted(glob.glob(f"{DATA}/gse126078/*_processed_data.txt.gz"))
    cols = {}
    for fp in files:
        name = os.path.basename(fp).replace("_processed_data.txt.gz", "")
        df = pd.read_csv(fp, sep="\t", usecols=["Gene_symbol", "FPKM"])
        cols[name] = df.dropna(subset=["Gene_symbol"]).groupby("Gene_symbol")["FPKM"].max()
    gse = pd.DataFrame(cols)
    return np.log2(gse + 1)
