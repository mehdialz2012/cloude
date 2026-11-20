"""Step 7 — data-integrity and statistical-independence checks on the cohorts as used.

C1  SU2C: samples vs patients (pseudo-replication in the 67-vs-67 DE and in the quartile split)
C2  GSE126078: what the 153 samples actually are, and patient-level structure
C3  the liver-score adjustment: does it do what it is claimed to do?
C4  GSEA permutation p-values (decoupler `times`) — the manuscript reports p<1e-3
C5  epithelial-adjusted DE for every nominated candidate (not only signature members)
C6  HPA version drift (number of consensus tissues)
C7  Beltran negative-floor clipping: sensitivity of the headline AUC
"""
import json
import glob
import os
import numpy as np
import pandas as pd
from scipy import stats as st
from statsmodels.stats.multitest import multipletests
from sklearn.metrics import roc_auc_score
import decoupler as dc

from common import (DATA, OUT, AR20, NE12, LIVER10, load_su2c, load_beltran,
                    load_gse126078, sig_score, zrow)

res = {}
su2c_log, su2c_cs, su2c_cp = load_su2c()
f = su2c_log.loc["FOLH1"]
low, high = f[f <= f.quantile(.25)].index, f[f >= f.quantile(.75)].index

# ---------------------------------------------------------------- C1
pid = su2c_cs["PATIENT_ID"].reindex(f.index)
print(f"[C1] SU2C RNA samples n={len(f)}; distinct patients n={pid.nunique()}; "
      f"patients with >1 sample: {int((pid.value_counts()>1).sum())}")
dups = pid.value_counts()
print(f"[C1] max samples per patient: {dups.max()}; "
      f"samples belonging to multi-sample patients: {int(dups[dups>1].sum())}")
print(f"[C1] PSMA-low group: {pid.loc[low].nunique()} patients / {len(low)} samples; "
      f"PSMA-high: {pid.loc[high].nunique()} patients / {len(high)} samples")
res["C1"] = dict(n_samples=int(len(f)), n_patients=int(pid.nunique()),
                 n_multi=int((dups > 1).sum()), max_per_patient=int(dups.max()),
                 low_patients=int(pid.loc[low].nunique()), high_patients=int(pid.loc[high].nunique()))
# patient-level (one sample per patient, highest-expressed first) re-test of a few key genes
first = pid.reset_index().drop_duplicates("PATIENT_ID")["Hugo_Symbol" if "Hugo_Symbol" in
                                                        pid.reset_index().columns else
                                                        pid.reset_index().columns[0]]
keep = pid.groupby(pid.values).apply(lambda s: s.index[0]).tolist()
fk = f.loc[keep]
lo2, hi2 = fk[fk <= fk.quantile(.25)].index, fk[fk >= fk.quantile(.75)].index
print(f"[C1] patient-level (1 sample/patient) split: n={len(fk)}, low {len(lo2)}, high {len(hi2)}")
for g in ["STEAP1", "PLAUR", "SLC16A3", "SLC7A5", "GPC2", "NCAM1"]:
    if g in su2c_log.index:
        U, p = st.mannwhitneyu(su2c_log.loc[g, lo2], su2c_log.loc[g, hi2])
        rbc = 2 * U / (len(lo2) * len(hi2)) - 1
        print(f"     {g:9s} patient-level rbc={rbc:+.3f} (p={p:.1e})")
        res[f"C1_{g}"] = dict(rbc=float(rbc), p=float(p))

# ---------------------------------------------------------------- C2
files = sorted(glob.glob(f"{DATA}/gse126078/*_processed_data.txt.gz"))
names = [os.path.basename(x).replace("_processed_data.txt.gz", "") for x in files]
cell = [n for n in names if any(k in n for k in ["C42B", "PC3", "PacMet"])]
lucap = [n for n in names if "LuCaP" in n]
pat = [n for n in names if n not in cell and n not in lucap]
print(f"\n[C2] GSE126078: {len(names)} files = {len(pat)} patient, {len(lucap)} LuCaP, "
      f"{len(cell)} 'cell line'")
print(f"[C2] the 12 'cell line' samples are siRNA/knockout perturbation experiments:")
for n in cell:
    print(f"       {n}")
# patient-level structure: sample names look like GSMxxxxxxx_<case>_<site>
stems = sorted(set("_".join(n.split("_")[1:3]) for n in pat))
cases = sorted(set(n.split("_")[1] for n in pat))
print(f"[C2] distinct leading identifiers among the 98 patient samples: {len(cases)}")
print(f"[C2] e.g. {cases[:12]}")
res["C2"] = dict(n_total=len(names), n_patient=len(pat), n_lucap=len(lucap), n_cell=len(cell),
                 cell_samples=cell, n_distinct_case_ids=len(cases), example_cases=cases[:20])

# ---------------------------------------------------------------- C3
liver_score = zrow(su2c_log.loc[[g for g in LIVER10 if g in su2c_log.index]]).mean(axis=0)
site = su2c_cs["TISSUE_SITE"].reindex(f.index)
is_liver = site.astype(str).str.contains("Liver", na=False)
print(f"\n[C3] liver-site samples: {int(is_liver.sum())}/{len(f)}")
print(f"[C3] liver score: liver biopsies mean={liver_score[is_liver].mean():+.2f}, "
      f"others mean={liver_score[~is_liver].mean():+.2f} "
      f"(AUC={roc_auc_score(is_liver.astype(int), liver_score):.3f})")
print(f"[C3] FOLH1 vs liver score WITHIN non-liver samples: "
      f"rho={st.spearmanr(f[~is_liver], liver_score[~is_liver])[0]:+.3f}")
print(f"[C3] FOLH1: liver biopsies median={f[is_liver].median():.2f} vs non-liver "
      f"{f[~is_liver].median():.2f} (MWU p={st.mannwhitneyu(f[is_liver], f[~is_liver]).pvalue:.1e})")
# after residualization, is site still predictive of the signature?
sig = pd.read_csv(f"{OUT}/s03_signature100.csv")
up50, dn50 = sig["up50"].dropna().tolist(), sig["dn50"].dropna().tolist()
z = zrow(su2c_log)
score = z.loc[[g for g in up50 if g in z.index]].mean(axis=0) - \
    z.loc[[g for g in dn50 if g in z.index]].mean(axis=0)
print(f"[C3] PSMA-low signature score: liver biopsies mean={score[is_liver].mean():+.2f} vs "
      f"others {score[~is_liver].mean():+.2f} (AUC for liver site = "
      f"{roc_auc_score(is_liver.astype(int), score):.3f})")
res["C3"] = dict(n_liver=int(is_liver.sum()),
                 liverscore_auc=float(roc_auc_score(is_liver.astype(int), liver_score)),
                 rho_within_nonliver=float(st.spearmanr(f[~is_liver], liver_score[~is_liver])[0]),
                 folh1_liver_median=float(f[is_liver].median()),
                 folh1_nonliver_median=float(f[~is_liver].median()),
                 sig_liver_auc=float(roc_auc_score(is_liver.astype(int), score)))

# ---------------------------------------------------------------- C4
de_adj = pd.read_csv(f"{OUT}/s03_de_adj.csv", index_col=0)
hm = dc.op.hallmark(organism="human")
rank = de_adj["rbc"].sort_values(ascending=False)
print("\n[C4] Hallmark GSEA permutation p-values at increasing permutation counts:")
for times in [1000, 10000]:
    s, p = dc.mt.gsea(data=rank.to_frame("rbc").T, net=hm, times=times, seed=1, verbose=False)
    g = s.T.copy(); g.columns = ["NES"]; g["p"] = p.T.values
    g = g.sort_values("NES", ascending=False)
    g["padj"] = multipletests(g["p"].clip(lower=1 / times), method="fdr_bh")[1]
    sel = g.loc[[i for i in ["E2F_TARGETS", "G2M_CHECKPOINT",
                             "EPITHELIAL_MESENCHYMAL_TRANSITION", "ANDROGEN_RESPONSE",
                             "OXIDATIVE_PHOSPHORYLATION"] if i in g.index]]
    print(f"  times={times}:")
    print(sel.round(4).to_string())
    res[f"C4_times{times}"] = sel.round(5).to_dict("index")

# ---------------------------------------------------------------- C5
EPI = ["EPCAM", "CDH1", "KRT8", "KRT18", "CLDN3", "CLDN4"]
comp_epi = sig_score(su2c_log, EPI)[0]
c = comp_epi.loc[su2c_log.columns].values
X = np.vstack([np.ones_like(c), c]).T
expressed = su2c_log.index[su2c_log.median(axis=1) >= 1.0]
Y = su2c_log.loc[expressed].values
beta, *_ = np.linalg.lstsq(X, Y.T, rcond=None)
Xe = pd.DataFrame(Y - (X @ beta).T, index=expressed, columns=su2c_log.columns)
rows = []
for g in ["PAQR8", "ATP2B4", "SLC12A4", "SLC16A3", "PLAUR", "SLC7A5", "ADAM19", "MRC2",
          "GPC2", "MUC1", "ACTN1", "VCAN", "MXRA5", "STEAP1", "STEAP2", "TACSTD2", "CD276"]:
    if g not in Xe.index:
        continue
    U, p = st.mannwhitneyu(Xe.loc[g, low], Xe.loc[g, high])
    rows.append(dict(gene=g, rbc_liver_adj=float(de_adj["rbc"].get(g, np.nan)),
                     rbc_epi_adj=2 * U / (len(low) * len(high)) - 1, p_epi_adj=p))
t5 = pd.DataFrame(rows).set_index("gene")
t5["padj_epi_adj"] = multipletests(t5["p_epi_adj"], method="fdr_bh")[1]
t5["shrinkage_pct"] = 100 * (1 - t5["rbc_epi_adj"] / t5["rbc_liver_adj"])
print("\n[C5] candidate effect sizes after adjusting for tumour (epithelial) content:")
print(t5.round(4).to_string())
res["C5"] = t5.round(4).to_dict("index")
t5.to_csv(f"{OUT}/s07_candidates_epi_adjusted.csv")

# ---------------------------------------------------------------- C6
hpa = pd.read_csv(f"{DATA}/rna_tissue_consensus.tsv", sep="\t")
hpa.columns = ["gene_id", "gene", "tissue", "ntpm"]
print(f"\n[C6] HPA consensus file: {hpa['tissue'].nunique()} tissues, "
      f"{hpa['gene'].nunique()} genes   [manuscript says 51 tissues]")
res["C6"] = dict(n_tissues=int(hpa["tissue"].nunique()), n_genes=int(hpa["gene"].nunique()))

# ---------------------------------------------------------------- C7
raw = pd.read_csv(f"{DATA}/nepc_data_mrna_seq_fpkm.txt", sep="\t", index_col=0, comment="#")
raw = raw[~raw.index.duplicated()].drop(columns=["Entrez_Gene_Id"], errors="ignore")
raw = raw.select_dtypes(include=[np.number])
nepc_cs = pd.read_csv(f"{DATA}/nepc_data_clinical_sample.txt", sep="\t", comment="#", index_col=0)
nmap = nepc_cs.set_index("OTHER_SAMPLE_ID")


def score_mat(mat):
    zz = zrow(mat)
    u = [g for g in up50 if g in mat.index and mat.loc[g].std() > 0]
    d = [g for g in dn50 if g in mat.index and mat.loc[g].std() > 0]
    return zz.loc[u].mean(axis=0) - zz.loc[d].mean(axis=0), len(u), len(d)


variants = {
    "clip at 0, log2(x+1)  [authors]": np.log2(raw.clip(lower=0) + 1),
    "shift by global min, log2(x+1)": np.log2(raw - np.nanmin(raw.values) + 1),
    "rank-based (per-sample ranks)": raw.rank(axis=0),
    "no transform (raw FPKM-like)": raw,
}
print("\n[C7] Beltran headline AUC under alternative handling of the negative floor:")
for nm, mat in variants.items():
    s, nu, nd = score_mat(mat)
    idx = s.index.intersection(nmap.index)
    yv = (nmap.loc[idx, "CANCER_TYPE_DETAILED"] == "Prostate Neuroendocrine Carcinoma").astype(int)
    a = roc_auc_score(yv, s[idx])
    ne, _ = sig_score(mat, NE12)
    a_ne = roc_auc_score(yv, ne[idx])
    print(f"   {nm:34s} signature AUC={a:.3f} (genes used {nu}+{nd})   NE-score AUC={a_ne:.3f}")
    res[f"C7::{nm}"] = dict(auc=float(a), genes=f"{nu}+{nd}", ne_auc=float(a_ne))

json.dump(res, open(f"{OUT}/s07_results.json", "w"), indent=1)
print("\nsaved s07 outputs")
