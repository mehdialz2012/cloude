"""Phase 2, step 1 — discovery in SU2C mCRPC (n=266).

FOLH1 is treated as a continuous variable. Three nested models are fitted for every
expressed gene, so that each association can be labelled by what explains it:

  M1  expression ~ FOLH1                                   (raw association)
  M2  M1 + immune + stromal + hepatocyte + osteoblast      (tumour composition)
          + biopsy site
  M3  M2 + neuroendocrine score                            (lineage)

Genes significant in M3 are tumour-intrinsic and not explained by composition, biopsy
site or neuroendocrine differentiation. Cut-off sensitivity and the quartile contrast
are reported alongside, not instead.
"""
import json
import numpy as np
import pandas as pd
from scipy import stats as st
import decoupler as dc

from common import OUT, AR20, NE12, load_su2c, zrow, sig_score
from p2_common import composition_scores, ols_beta, design, fdr

res = {}
su2c, cs, cp = load_su2c()
f = su2c.loc["FOLH1"]
expressed = su2c.index[su2c.median(axis=1) >= 1.0]
X = su2c.loc[expressed]
print(f"SU2C: {su2c.shape[1]} samples, {len(expressed)} expressed genes")

comp = composition_scores(su2c)
comp["FOLH1"] = f
comp["FOLH1_z"] = (f - f.mean()) / f.std()
site = cs["TISSUE_SITE"].reindex(f.index).fillna("Unknown")
ne_feat = cs["NEUROENDOCRINE_FEATURES"].reindex(f.index)
comp["site"] = site
comp.to_csv(f"{OUT}/p2_su2c_covariates.csv")

print("\n[1] correlation of FOLH1 with each covariate (Spearman):")
for c in ["immune", "stromal", "hepatocyte", "osteoblast", "epithelial", "AR", "NE"]:
    r, p = st.spearmanr(comp["FOLH1"], comp[c])
    print(f"    {c:12s} rho={r:+.3f}  p={p:.1e}")
    res[f"rho_folh1_{c}"] = float(r)
print(f"    biopsy site: Kruskal p="
      f"{st.kruskal(*[g['FOLH1'].values for _, g in comp.groupby('site') if len(g) >= 5]).pvalue:.1e}")

# ---------------------------------------------------------------- nested models
cov2 = ["immune", "stromal", "hepatocyte", "osteoblast"]
M1 = design(comp, ["FOLH1_z"])
M2 = design(comp, ["FOLH1_z"] + cov2, dummies={"site": site})
M3 = design(comp, ["FOLH1_z"] + cov2 + ["NE"], dummies={"site": site})
Y = X[comp.index].values

tab = pd.DataFrame(index=expressed)
for name, D in [("M1", M1), ("M2", M2), ("M3", M3)]:
    b, se, t, p = ols_beta(Y, D.values, coef=1)
    tab[f"beta_{name}"] = b
    tab[f"p_{name}"] = p
    tab[f"fdr_{name}"] = fdr(p)
    n_sig = int((tab[f"fdr_{name}"] < .05).sum())
    print(f"[2] {name}: {n_sig} genes at FDR<0.05 "
          f"({int(((tab[f'fdr_{name}']<.05)&(b<0)).sum())} higher in PSMA-low, "
          f"{int(((tab[f'fdr_{name}']<.05)&(b>0)).sum())} lower)")
    res[f"n_sig_{name}"] = n_sig

sig3 = tab["fdr_M3"] < .05
sig1 = tab["fdr_M1"] < .05
sig2 = tab["fdr_M2"] < .05
tab["class"] = np.where(sig3, "robust",
                 np.where(sig2, "NE-dependent",
                   np.where(sig1, "composition/site-dependent", "ns")))
print("\n[3] what explains each association:")
print(tab["class"].value_counts().to_string())
res["classes"] = {k: int(v) for k, v in tab["class"].value_counts().items()}
print(f"[3] median |beta| shrinkage M1->M2: "
      f"{100*(1-(tab.loc[sig1,'beta_M2'].abs()/tab.loc[sig1,'beta_M1'].abs()).median()):.1f}%")
print(f"[3] median |beta| shrinkage M2->M3: "
      f"{100*(1-(tab.loc[sig2,'beta_M3'].abs()/tab.loc[sig2,'beta_M2'].abs()).median()):.1f}%")
res["shrink_M1_M2"] = float(100 * (1 - (tab.loc[sig1, "beta_M2"].abs()
                                        / tab.loc[sig1, "beta_M1"].abs()).median()))
res["shrink_M2_M3"] = float(100 * (1 - (tab.loc[sig2, "beta_M3"].abs()
                                        / tab.loc[sig2, "beta_M2"].abs()).median()))

# ---------------------------------------------------------------- cut-off sensitivity
print("\n[4] cut-off sensitivity (rank-biserial, PSMA-low vs PSMA-high):")


def de_cut(frac):
    lo = f[f <= f.quantile(frac)].index
    hi = f[f >= f.quantile(1 - frac)].index
    out = []
    for g in expressed:
        a, b_ = su2c.loc[g, lo].values, su2c.loc[g, hi].values
        U, p = st.mannwhitneyu(a, b_)
        out.append(2 * U / (len(a) * len(b_)) - 1)
    return pd.Series(out, index=expressed), len(lo), len(hi)


cuts = {}
for frac, nm in [(0.10, "bottom/top decile"), (0.25, "quartile (as published)"),
                 (1 / 3, "tercile"), (0.50, "median split")]:
    s, nlo, nhi = de_cut(frac)
    cuts[nm] = s
    print(f"    {nm:26s} n={nlo}/{nhi}  rho with continuous beta_M1 = "
          f"{st.spearmanr(s, -tab['beta_M1'])[0]:.3f}")
    res[f"cut_{nm}"] = float(st.spearmanr(s, -tab["beta_M1"])[0])
cut_df = pd.DataFrame(cuts)
cut_df["beta_M1"] = tab["beta_M1"]
cut_df["beta_M3"] = tab["beta_M3"]
cut_df.to_csv(f"{OUT}/p2_cutoff_sensitivity.csv")
corr = cut_df[list(cuts)].corr(method="spearman")
print("    pairwise concordance between cut-offs:")
print(corr.round(3).to_string())

# ---------------------------------------------------------------- programs
print("\n[5] Hallmark programmes (10,000 permutations):")
hm = dc.op.hallmark(organism="human")
gsea_out = {}
for name in ["M1", "M2", "M3"]:
    rank = (-tab[f"beta_{name}"]).sort_values(ascending=False)   # positive = up in PSMA-low
    s, p = dc.mt.gsea(data=rank.to_frame("x").T, net=hm, times=10000, seed=1, verbose=False)
    g = s.T.copy()
    g.columns = ["NES"]
    g["p"] = p.T.values
    g["fdr"] = fdr(g["p"].clip(lower=1e-4))
    gsea_out[name] = g
    if name == "M3":
        print("    top up in PSMA-low (M3, fully adjusted):")
        print(g.sort_values("NES", ascending=False).head(8).round(3).to_string())
        print("    top down:")
        print(g.sort_values("NES").head(6).round(3).to_string())
gs = pd.concat({k: v["NES"] for k, v in gsea_out.items()}, axis=1)
gs.columns = [f"NES_{c}" for c in gs.columns]
gs["fdr_M3"] = gsea_out["M3"]["fdr"]
gs = gs.sort_values("NES_M3", ascending=False)
gs.to_csv(f"{OUT}/p2_gsea_nested.csv")
print(f"    NES concordance M1 vs M3: rho={st.spearmanr(gs['NES_M1'], gs['NES_M3'])[0]:.3f}")
res["gsea_M1_M3_concordance"] = float(st.spearmanr(gs["NES_M1"], gs["NES_M3"])[0])

tab["median_expr"] = su2c.loc[expressed].median(axis=1)
tab.to_csv(f"{OUT}/p2_discovery_models.csv")
json.dump(res, open(f"{OUT}/p2_01_results.json", "w"), indent=1)
print("\nsaved p2 step 1")
