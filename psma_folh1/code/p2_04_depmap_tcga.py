"""Phase 2, steps 4 and 5.

Step 4 — DepMap, explicitly exploratory. The dependency scan is repeated with the
cross-cohort-validated programme rather than the discovery-only signature, and is
reported with the lineage structure quantified rather than asserted away.

Step 5 — TCGA primary tumours as a secondary comparison: is the PSMA-low programme a
feature of advanced disease, or is it already present in treatment-naive primary cancer?
"""
import json
import numpy as np
import pandas as pd
from scipy import stats as st
from scipy.stats import rankdata, norm
from statsmodels.stats.multitest import multipletests
import decoupler as dc

from common import DATA, OUT, AR20, NE12, load_tcga, load_su2c, sig_score
from p2_common import composition_scores, ols_beta, design, fdr

res = {}

# ---------------------------------------------------------------- validated programme
X = pd.read_csv(f"{OUT}/p2_crosscohort_genes.csv", index_col=0)
rep = X[X["reproducible"]]
up = rep[rep["su2c_M3"] < 0].sort_values("su2c_M3").index.tolist()      # higher in PSMA-low
dn = rep[rep["su2c_M3"] > 0].sort_values("su2c_M3", ascending=False).index.tolist()
print(f"[4] cross-cohort validated programme: {len(up)} genes up in PSMA-low, {len(dn)} down")
res["prog_up"] = len(up)
res["prog_dn"] = len(dn)
pd.DataFrame({"gene": up + dn,
              "direction": ["up_in_PSMA_low"] * len(up) + ["down_in_PSMA_low"] * len(dn)}
             ).to_csv(f"{OUT}/p2_validated_programme.csv", index=False)

# ---------------------------------------------------------------- step 4 DepMap
model = pd.read_csv(f"{DATA}/depmap_Model.csv")
crispr = pd.read_csv(f"{DATA}/depmap_CRISPRGeneEffect.csv", index_col=0)
crispr.columns = [c.split(" ")[0] for c in crispr.columns]
crispr = crispr.loc[:, ~crispr.columns.duplicated()].astype("float32")
dexpr = pd.read_csv(f"{DATA}/depmap_ExpressionTPMLogp1.csv", index_col=0)
dexpr.columns = [c.split(" ")[0] for c in dexpr.columns]
dexpr = dexpr.loc[:, ~dexpr.columns.duplicated()].astype("float32")
ug = [g for g in up if g in dexpr.columns]
dg = [g for g in dn if g in dexpr.columns]
z = dexpr[ug + dg].apply(lambda x: (x - x.mean()) / x.std())
score = z[ug].mean(axis=1) - z[dg].mean(axis=1)
lin = model.set_index("ModelID")["OncotreeLineage"].reindex(score.index)
cm = crispr.index.intersection(score.index)
sc = score.loc[cm]
C = crispr.loc[cm]
print(f"[4] scored {len(cm)} DepMap models with the validated programme "
      f"({len(ug)}/{len(up)} up and {len(dg)}/{len(dn)} down genes measured)")

ls = pd.DataFrame({"score": score, "lineage": lin}).dropna()
eta = 1 - (ls.groupby("lineage")["score"].transform(lambda x: x - x.mean()).var()
           / ls["score"].var())
print(f"[4] variance of the programme score explained by lineage alone: {eta*100:.1f}%")
agg = ls.groupby("lineage")["score"].agg(["count", "mean"]).sort_values("mean")
print(f"[4] lowest lineages: {agg.head(3)['mean'].round(2).to_dict()}; "
      f"highest: {agg.tail(3)['mean'].round(2).to_dict()}")
res["depmap_lineage_var"] = float(eta)

scr = rankdata(sc.values); scr = (scr - scr.mean()) / scr.std()
rows = []
V = C.values
for j, g in enumerate(C.columns):
    v = V[:, j]
    m = ~np.isnan(v)
    if m.sum() < 500:
        continue
    rows.append((g, np.corrcoef(scr[m], rankdata(v[m]))[0, 1], int(m.sum())))
dep = pd.DataFrame(rows, columns=["gene", "rho", "n"]).set_index("gene")
dep["p"] = 2 * norm.sf(np.abs(dep["rho"] * np.sqrt(dep["n"] - 1)))
dep["fdr"] = multipletests(dep["p"], method="fdr_bh")[1]

# within-lineage
lin_c = lin.reindex(C.index)
keep = lin_c.value_counts()
keep = keep[keep >= 15].index
mask = lin_c.isin(keep).values
sc_l = sc[mask]
lin_l = lin_c[mask]
sc_dm = sc_l.groupby(lin_l.values).transform(lambda x: x - x.mean()).values
Cl = C.loc[sc_l.index]
codes = pd.Categorical(lin_l).codes
rows = []
for j, g in enumerate(Cl.columns):
    v = Cl.values[:, j]
    m = ~np.isnan(v)
    if m.sum() < 400:
        continue
    vv, lc = v[m], codes[m]
    s_ = np.bincount(lc, weights=vv); c_ = np.bincount(lc)
    rows.append((g, np.corrcoef(rankdata(sc_dm[m]), rankdata(vv - (s_ / np.maximum(c_, 1))[lc]))[0, 1]))
depwl = pd.DataFrame(rows, columns=["gene", "rho_within_lineage"]).set_index("gene")
dep = dep.join(depwl)
frac = (C < -1).mean(axis=0)
ce = set(frac[frac >= .9].index)
dep["common_essential"] = dep.index.isin(ce)
nn = dep[~dep["common_essential"]]
print("\n[4] top dependencies of programme-high (PSMA-low-like) cells, non-essential genes:")
top = nn.sort_values("rho").head(10)
for g, r in top.iterrows():
    rk_wl = int((nn["rho_within_lineage"] < r["rho_within_lineage"]).sum() + 1)
    print(f"      {g:10s} pan rho={r['rho']:+.3f} (rank {int((nn['rho']<r['rho']).sum()+1)})  "
          f"within-lineage {r['rho_within_lineage']:+.3f} (rank {rk_wl})")
for g in ["TXNRD1", "SEPSECS", "EEFSEC", "PSTK", "GPX4", "GAPDH", "FOXA1"]:
    if g in dep.index:
        res[f"dep_{g}"] = dict(pan=float(dep.loc[g, "rho"]),
                               wl=float(dep.loc[g, "rho_within_lineage"]),
                               rank_pan=int((nn["rho"] < dep.loc[g, "rho"]).sum() + 1),
                               rank_wl=int((nn["rho_within_lineage"]
                                            < dep.loc[g, "rho_within_lineage"]).sum() + 1))
        print(f"      [axis] {g:8s} pan {dep.loc[g,'rho']:+.3f} (rank {res[f'dep_{g}']['rank_pan']}) "
              f"| within-lineage {dep.loc[g,'rho_within_lineage']:+.3f} "
              f"(rank {res[f'dep_{g}']['rank_wl']})")
fo = dexpr.loc[C.index, "FOLH1"]
print("\n[4] tested directly against FOLH1 mRNA (the construct itself):")
for g in ["TXNRD1", "SEPSECS", "EEFSEC", "PSTK", "GPX4"]:
    if g in C.columns:
        v = C[g]
        m = v.notna() & fo.notna()
        r = st.spearmanr(fo[m], v[m])
        print(f"      {g:8s} rho(FOLH1, gene effect) = {r[0]:+.3f} (p={r[1]:.3f})  "
              f"[the claim predicts a positive value]")
        res[f"folh1_direct_{g}"] = float(r[0])
dep.to_csv(f"{OUT}/p2_depmap_validated_programme.csv")

# ---------------------------------------------------------------- step 5 TCGA
tcga, tclin = load_tcga()
tum = [c for c in tcga.columns if c.endswith("-01")]
T = tcga[tum]
ft = T.loc["FOLH1"]
su2c, cs, cp = load_su2c()
fs = su2c.loc["FOLH1"]
print(f"\n[5] TCGA primary tumours: n={T.shape[1]}; FOLH1 IQR "
      f"{ft.quantile(.25):.2f}-{ft.quantile(.75):.2f} (SD {ft.std():.2f})")
print(f"[5] SU2C mCRPC:            n={su2c.shape[1]}; FOLH1 IQR "
      f"{fs.quantile(.25):.2f}-{fs.quantile(.75):.2f} (SD {fs.std():.2f})")
res["tcga_folh1_sd"] = float(ft.std())
res["su2c_folh1_sd"] = float(fs.std())
expr_t = T.index[T.median(axis=1) >= 1.0]
ct = composition_scores(T)
ct["FOLH1_z"] = (ft - ft.mean()) / ft.std()
cov = ["immune", "stromal", "hepatocyte", "osteoblast"]
bt, _, _, pt = ols_beta(T.loc[expr_t].values, design(ct, ["FOLH1_z"] + cov + ["NE"]).values, coef=1)
tt = pd.DataFrame({"beta_M3": bt, "p_M3": pt}, index=expr_t)
tt["fdr_M3"] = fdr(tt["p_M3"])
disc = pd.read_csv(f"{OUT}/p2_discovery_models.csv", index_col=0)
ov = disc.index.intersection(tt.index)
r_prog = st.spearmanr(disc.loc[ov, "beta_M3"], tt.loc[ov, "beta_M3"])
print(f"[5] TCGA genes associated with FOLH1 (fully adjusted, FDR<0.05): "
      f"{int((tt['fdr_M3']<.05).sum())}")
print(f"[5] concordance of effect sizes with the mCRPC discovery: rho={r_prog[0]:.3f} "
      f"(p={r_prog[1]:.1e})")
core = X[X["reproducible"]].index.intersection(tt.index)
same = (np.sign(disc.loc[core, "beta_M3"]) == np.sign(tt.loc[core, "beta_M3"]))
print(f"[5] of the {len(core)} validated-programme genes measured in TCGA, "
      f"{int(same.sum())} ({100*same.mean():.0f}%) have the same direction in primary disease, "
      f"{int((same & (tt.loc[core,'p_M3']<.05)).sum())} also nominally significant")
res["tcga_n_sig"] = int((tt["fdr_M3"] < .05).sum())
res["tcga_concordance"] = float(r_prog[0])
res["tcga_same_direction_frac"] = float(same.mean())
hm = dc.op.hallmark(organism="human")
s_, p_ = dc.mt.gsea(data=(-tt["beta_M3"]).sort_values(ascending=False).to_frame("x").T,
                    net=hm, times=10000, seed=1, verbose=False)
N = pd.read_csv(f"{OUT}/p2_programme_concordance.csv", index_col=0)
N["TCGA_primary"] = s_.T.iloc[:, 0]
print("[5] key programmes, metastatic vs primary (NES):")
print(N.loc[[i for i in ["E2F_TARGETS", "G2M_CHECKPOINT", "ANDROGEN_RESPONSE",
                         "OXIDATIVE_PHOSPHORYLATION", "EPITHELIAL_MESENCHYMAL_TRANSITION"]
             if i in N.index]].round(2).to_string())
N.to_csv(f"{OUT}/p2_programme_concordance.csv")
tt.to_csv(f"{OUT}/p2_tcga_models.csv")
json.dump(res, open(f"{OUT}/p2_04_results.json", "w"), indent=1)
print("\nsaved p2 steps 4-5")
