"""Phase 2, step 2 — independent validation in GSE126078 and Beltran.

The question is not whether a signature anti-correlates with FOLH1 again. It is whether
the same gene-level associations, the same programmes and the same candidate targets
reappear in independent metastatic cohorts once composition, site and lineage are
handled the same way.

GSE126078 carries the published Labrecque molecular phenotype calls, so the key clinical
question can be asked directly: does a PSMA-low program exist inside AR-positive,
non-neuroendocrine disease — the compartment that is actually underserved?
"""
import json
import numpy as np
import pandas as pd
from scipy import stats as st
import decoupler as dc

from common import DATA, OUT, AR20, NE12, load_beltran, sig_score, zrow, load_gse126078
from p2_common import composition_scores, ols_beta, design, fdr

res = {}
disc = pd.read_csv(f"{OUT}/p2_discovery_models.csv", index_col=0)

# ---------------------------------------------------------------- GSE126078
gse_all = load_gse126078()
meta = pd.read_csv(f"{DATA}/gse126078_meta.csv")
meta["key"] = meta["gsm"] + "_" + meta["title"]
name_map = {}
for c in gse_all.columns:
    gsm = c.split("_")[0]
    hit = meta[meta["gsm"] == gsm]
    if len(hit):
        name_map[c] = hit.iloc[0]
m = meta.set_index("gsm")
cols = [c for c in gse_all.columns if m.loc[c.split("_")[0], "patient"] is not np.nan
        and isinstance(m.loc[c.split("_")[0], "patient"], str)]
gse = gse_all[cols]
info = pd.DataFrame({
    "patient": [m.loc[c.split("_")[0], "patient"] for c in cols],
    "site": [m.loc[c.split("_")[0], "tumor site"] for c in cols],
    "phenotype": [m.loc[c.split("_")[0], "molecular phenotype"] for c in cols]},
    index=cols)
print(f"GSE126078 patient samples: {gse.shape[1]} from {info['patient'].nunique()} patients")
print(info["phenotype"].value_counts().to_string())
print(info["site"].value_counts().head(8).to_string())
res["gse_n_samples"] = int(gse.shape[1])
res["gse_n_patients"] = int(info["patient"].nunique())

fg = gse.loc["FOLH1"]
print("\n[1] FOLH1 by Labrecque molecular phenotype (median log2 FPKM+1):")
for ph, sub in info.groupby("phenotype"):
    v = fg[sub.index]
    print(f"    {ph:14s} n={len(v):3d}  median={v.median():5.2f}  "
          f"{int((v <= fg.quantile(.25)).sum())} of {len(v)} in the cohort's bottom quartile")
    res[f"gse_folh1_{ph}"] = float(v.median())

expr_g = gse.index[gse.median(axis=1) >= 1.0]
comp_g = composition_scores(gse)
comp_g["FOLH1_z"] = (fg - fg.mean()) / fg.std()
cov = ["immune", "stromal", "hepatocyte", "osteoblast"]
Yg = gse.loc[expr_g].values
Mg1 = design(comp_g, ["FOLH1_z"])
Mg2 = design(comp_g, ["FOLH1_z"] + cov, dummies={"site": info["site"]})
Mg3 = design(comp_g, ["FOLH1_z"] + cov + ["NE"], dummies={"site": info["site"]})
val = pd.DataFrame(index=expr_g)
for nm, D in [("M1", Mg1), ("M2", Mg2), ("M3", Mg3)]:
    b, se, t, p = ols_beta(Yg, D.values, coef=1)
    val[f"beta_{nm}"] = b
    val[f"p_{nm}"] = p
    val[f"fdr_{nm}"] = fdr(p)
    print(f"[2] GSE126078 {nm}: {int((val[f'fdr_{nm}']<.05).sum())} genes at FDR<0.05")
    res[f"gse_n_sig_{nm}"] = int((val[f"fdr_{nm}"] < .05).sum())

# patient-level (one observation per patient) — guards against pseudo-replication
pl = gse.loc[expr_g].T.groupby(info["patient"].values).mean().T
compl_ = composition_scores(pl)
fp = pl.loc["FOLH1"]
compl_["FOLH1_z"] = (fp - fp.mean()) / fp.std()
Dp = design(compl_, ["FOLH1_z"] + cov + ["NE"])
b, se, t, p = ols_beta(pl.values, Dp.values, coef=1)
val["beta_patient_M3"] = pd.Series(b, index=pl.index).reindex(val.index)
val["p_patient_M3"] = pd.Series(p, index=pl.index).reindex(val.index)
val["fdr_patient_M3"] = fdr(val["p_patient_M3"])
print(f"[2] GSE126078 patient-level (n={pl.shape[1]}) M3: "
      f"{int((val['fdr_patient_M3']<.05).sum())} genes at FDR<0.05")
res["gse_n_patients_model"] = int(pl.shape[1])
res["gse_n_sig_patient_M3"] = int((val["fdr_patient_M3"] < .05).sum())

# ---------------------------------------------------------------- Beltran
bel, bcs = load_beltran()
bm = bcs.set_index("OTHER_SAMPLE_ID")
bidx = bel.columns.intersection(bm.index)
bel = bel[bidx]
bpath = bm.loc[bidx, "CANCER_TYPE_DETAILED"]
fb = bel.loc["FOLH1"]
expr_b = bel.index[bel.median(axis=1) >= 1.0]
comp_b = composition_scores(bel)
comp_b["FOLH1_z"] = (fb - fb.mean()) / fb.std()
Db1 = design(comp_b, ["FOLH1_z"])
Db3 = design(comp_b, ["FOLH1_z"] + cov + ["NE"])
belt = pd.DataFrame(index=expr_b)
for nm, D in [("M1", Db1), ("M3", Db3)]:
    b, se, t, p = ols_beta(bel.loc[expr_b].values, D.values, coef=1)
    belt[f"beta_{nm}"] = b
    belt[f"p_{nm}"] = p
    belt[f"fdr_{nm}"] = fdr(p)
print(f"\n[3] Beltran n={bel.shape[1]} ({int((bpath=='Prostate Neuroendocrine Carcinoma').sum())} NEPC): "
      f"M1 {int((belt['fdr_M1']<.05).sum())} genes, M3 {int((belt['fdr_M3']<.05).sum())} at FDR<0.05")
res["beltran_n_sig_M1"] = int((belt["fdr_M1"] < .05).sum())
res["beltran_n_sig_M3"] = int((belt["fdr_M3"] < .05).sum())

# ---------------------------------------------------------------- cross-cohort
common = disc.index.intersection(val.index).intersection(belt.index)
X = pd.DataFrame({
    "su2c_M1": disc.loc[common, "beta_M1"], "su2c_M3": disc.loc[common, "beta_M3"],
    "gse_M1": val.loc[common, "beta_M1"], "gse_M3": val.loc[common, "beta_M3"],
    "gse_patient_M3": val.loc[common, "beta_patient_M3"],
    "beltran_M1": belt.loc[common, "beta_M1"], "beltran_M3": belt.loc[common, "beta_M3"]})
print(f"\n[4] genes measured in all three cohorts: {len(common)}")
print("[4] cross-cohort concordance of FOLH1 effect sizes (Spearman):")
cc = X.corr(method="spearman")
print(cc.round(3).to_string())
res["crosscohort"] = cc.round(4).to_dict()

sig_su2c = disc.loc[common, "fdr_M3"] < .05
rep_gse = (np.sign(X["su2c_M3"]) == np.sign(X["gse_M3"])) & (val.loc[common, "p_M3"] < .05)
rep_pat = (np.sign(X["su2c_M3"]) == np.sign(X["gse_patient_M3"])) & (val.loc[common, "p_patient_M3"] < .05)
rep_bel = (np.sign(X["su2c_M3"]) == np.sign(X["beltran_M3"])) & (belt.loc[common, "p_M3"] < .05)
core = sig_su2c & rep_gse
core_strict = sig_su2c & rep_gse & rep_bel
print(f"[4] SU2C-robust genes (M3 FDR<0.05) present in all cohorts: {int(sig_su2c.sum())}")
print(f"[4]   ... reproduced in GSE126078 (same sign, p<0.05): {int(core.sum())} "
      f"({100*core.sum()/max(sig_su2c.sum(),1):.0f}%)")
print(f"[4]   ... also at patient level: {int((sig_su2c & rep_pat).sum())}")
print(f"[4]   ... and in Beltran as well: {int(core_strict.sum())}")
res["n_su2c_robust"] = int(sig_su2c.sum())
res["n_rep_gse"] = int(core.sum())
res["n_rep_patient"] = int((sig_su2c & rep_pat).sum())
res["n_rep_all3"] = int(core_strict.sum())

X["su2c_fdr_M3"] = disc.loc[common, "fdr_M3"]
X["gse_p_M3"] = val.loc[common, "p_M3"]
X["gse_patient_p_M3"] = val.loc[common, "p_patient_M3"]
X["beltran_p_M3"] = belt.loc[common, "p_M3"]
X["reproducible"] = core
X["reproducible_all3"] = core_strict
X["direction"] = np.where(X["su2c_M3"] < 0, "up in PSMA-low", "down in PSMA-low")
X.to_csv(f"{OUT}/p2_crosscohort_genes.csv")

# ---------------------------------------------------------------- programmes
print("\n[5] Hallmark programme concordance across cohorts (fully adjusted ranking):")
hm = dc.op.hallmark(organism="human")
nes = {}
for nm, series in [("SU2C", -disc["beta_M3"]), ("GSE126078", -val["beta_M3"]),
                   ("Beltran", -belt["beta_M3"])]:
    r = series.dropna().sort_values(ascending=False)
    s, p = dc.mt.gsea(data=r.to_frame("x").T, net=hm, times=10000, seed=1, verbose=False)
    nes[nm] = s.T.iloc[:, 0]
N = pd.DataFrame(nes)
print(N.loc[N["SU2C"].abs().sort_values(ascending=False).index].head(10).round(2).to_string())
print(f"    SU2C vs GSE126078 NES rho = {st.spearmanr(N['SU2C'], N['GSE126078'])[0]:.3f}")
print(f"    SU2C vs Beltran    NES rho = {st.spearmanr(N['SU2C'], N['Beltran'])[0]:.3f}")
res["nes_rho_gse"] = float(st.spearmanr(N["SU2C"], N["GSE126078"])[0])
res["nes_rho_beltran"] = float(st.spearmanr(N["SU2C"], N["Beltran"])[0])
N.to_csv(f"{OUT}/p2_programme_concordance.csv")

# ---------------------------------------------------------------- AR+/NE- subgroup
print("\n[6] does a PSMA-low program exist inside AR-positive, non-NE disease?")
sub = info[info["phenotype"] == "ARpos_NEneg"].index
print(f"    ARpos_NEneg samples: {len(sub)} from {info.loc[sub,'patient'].nunique()} patients; "
      f"FOLH1 median={gse.loc['FOLH1', sub].median():.2f}, "
      f"range {gse.loc['FOLH1', sub].min():.2f}-{gse.loc['FOLH1', sub].max():.2f}")
gs_sub = gse[sub]
cs_sub = composition_scores(gs_sub)
fs = gs_sub.loc["FOLH1"]
cs_sub["FOLH1_z"] = (fs - fs.mean()) / fs.std()
Ds = design(cs_sub, ["FOLH1_z"] + cov, dummies={"site": info.loc[sub, "site"]})
b, se, t, p = ols_beta(gs_sub.loc[expr_g].values, Ds.values, coef=1)
subres = pd.DataFrame({"beta": b, "p": p}, index=expr_g)
subres["fdr"] = fdr(subres["p"])
print(f"    genes associated with FOLH1 within ARpos_NEneg (FDR<0.05): "
      f"{int((subres['fdr']<.05).sum())}")
ov = disc.index.intersection(subres.index)
print(f"    concordance of effect sizes with the SU2C discovery (M3): "
      f"rho={st.spearmanr(disc.loc[ov,'beta_M3'], subres.loc[ov,'beta'])[0]:.3f}")
res["arpos_neneg_n"] = int(len(sub))
res["arpos_neneg_n_sig"] = int((subres["fdr"] < .05).sum())
res["arpos_neneg_concordance"] = float(st.spearmanr(disc.loc[ov, "beta_M3"],
                                                    subres.loc[ov, "beta"])[0])
subres.to_csv(f"{OUT}/p2_arpos_neneg_subgroup.csv")

val.to_csv(f"{OUT}/p2_gse_models.csv")
belt.to_csv(f"{OUT}/p2_beltran_models.csv")
info.to_csv(f"{OUT}/p2_gse_info.csv")
json.dump(res, open(f"{OUT}/p2_02_results.json", "w"), indent=1)
print("\nsaved p2 step 2")
