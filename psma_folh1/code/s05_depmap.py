"""Step 5 — DepMap dependency claims (manuscript Results section 5).

PART A (replication)
  A1  signature score across DepMap 24Q4 models; prostate-line ordering
  A2  genome-wide Spearman(score, Chronos gene effect); top hits; FDR count
  A3  common-essential exclusion and the "rank 1,2,4,15 of 17,293" claim
  A4  epithelial-restricted re-analysis; per-line fitness adjustment; GPX4 partials
PART B (independent validation)
  B1  what IS the score across lineages? (is "PSMA-low-like" a lineage axis?)
  B2  WITHIN-LINEAGE (lineage-demeaned) correlations — the decisive confounding test
  B3  haematopoietic-line removal
  B4  direct test of the actual construct: FOLH1 expression vs TXNRD1 dependency
  B5  prostate-lineage-only analysis, with benign lines shown
  B6  effect-size yardstick: where does rho=-0.23 sit among real selective dependencies?
"""
import json
import numpy as np
import pandas as pd
from scipy import stats as st
from scipy.stats import rankdata, norm
from statsmodels.stats.multitest import multipletests

from common import DATA, OUT

res = {}
AXIS = ["TXNRD1", "SEPSECS", "EEFSEC", "PSTK", "GPX4", "RHOA", "ANGPT2", "NCKAP1L",
        "GAPDH", "MTF1", "KEAP1", "NFE2L2", "FOXA1", "CCND1", "CDK4"]

# ------------------------------------------------------------------ load
model = pd.read_csv(f"{DATA}/depmap_Model.csv")
crispr = pd.read_csv(f"{DATA}/depmap_CRISPRGeneEffect.csv", index_col=0)
crispr.columns = [c.split(" ")[0] for c in crispr.columns]
crispr = crispr.loc[:, ~crispr.columns.duplicated()].astype("float32")
dexpr = pd.read_csv(f"{DATA}/depmap_ExpressionTPMLogp1.csv", index_col=0)
dexpr.columns = [c.split(" ")[0] for c in dexpr.columns]
dexpr = dexpr.loc[:, ~dexpr.columns.duplicated()].astype("float32")
print(f"CRISPR {crispr.shape}; expression {dexpr.shape}")
both = crispr.index.intersection(dexpr.index)
print(f"models with both: {len(both)}   [manuscript 1,103]")
res["n_models"] = int(len(both))

sig = pd.read_csv(f"{OUT}/s03_signature100.csv")
up50, dn50 = sig["up50"].dropna().tolist(), sig["dn50"].dropna().tolist()
ug = [g for g in up50 if g in dexpr.columns]
dg = [g for g in dn50 if g in dexpr.columns]
print(f"signature genes present in DepMap: up {len(ug)}/50, down {len(dg)}/50")
z = dexpr[ug + dg].apply(lambda x: (x - x.mean()) / x.std())
score = z[ug].mean(axis=1) - z[dg].mean(axis=1)

lin = model.set_index("ModelID")["OncotreeLineage"].reindex(score.index)
dis = model.set_index("ModelID")["OncotreePrimaryDisease"].reindex(score.index)
names = model.set_index("ModelID")["StrippedCellLineName"]

# ------------------------------------------------------------------ A1 prostate lines
pr = model[model["OncotreeLineage"] == "Prostate"]["ModelID"]
pr = [m for m in pr if m in score.index]
print(f"\n[A1] prostate models with expression: {len(pr)}")
pt = pd.DataFrame({"line": names.reindex(pr), "score": score.reindex(pr),
                   "FOLH1": dexpr.loc[pr, "FOLH1"],
                   "disease": dis.reindex(pr)}).sort_values("score", ascending=False)
print(pt.round(2).to_string(index=False))
res["A1_prostate"] = pt.reset_index().to_dict("records")

# ------------------------------------------------------------------ A2 genome-wide
cm = crispr.index.intersection(score.index)
sc = score.loc[cm]
C = crispr.loc[cm]
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
dep["z"] = dep["rho"] * np.sqrt(dep["n"] - 1)
dep["p"] = 2 * norm.sf(np.abs(dep["z"]))
dep["padj"] = multipletests(dep["p"], method="fdr_bh")[1]
print(f"\n[A2] pan-lineage n={len(cm)}; genes tested {len(dep)}; "
      f"FDR<0.05 {int((dep.padj<.05).sum())}")
print("[A2] top 12 negative (dependency of PSMA-low-like cells):")
print(dep.sort_values("rho").head(12)[["rho", "n", "padj"]].round(4).to_string())
print("[A2] top 8 positive:")
print(dep.sort_values("rho", ascending=False).head(8)[["rho", "n", "padj"]].round(4).to_string())
res["A2"] = dict(n_genes=int(len(dep)), n_fdr05=int((dep.padj < .05).sum()),
                 axis={g: float(dep.loc[g, "rho"]) for g in AXIS if g in dep.index})
dep.to_csv(f"{OUT}/s05_dep_genomewide.csv")

# ------------------------------------------------------------------ A3 common essentials
frac = (C < -1).mean(axis=0)
ce = set(frac[frac >= .9].index)
dep["common_essential"] = dep.index.isin(ce)
noness = dep[~dep["common_essential"]]
print(f"\n[A3] common essentials (Chronos<-1 in >=90% lines): {len(ce)}; "
      f"non-essential universe {len(noness)}   [manuscript 17,293]")
for g in ["SEPSECS", "TXNRD1", "EEFSEC", "PSTK", "GPX4", "GAPDH"]:
    if g in dep.index:
        rk = int((noness["rho"] < dep.loc[g, "rho"]).sum() + 1)
        print(f"   {g:8s} rho={dep.loc[g,'rho']:+.3f} rank {rk}/{len(noness)} "
              f"(common-essential: {g in ce})")
        res[f"A3_rank_{g}"] = rk
res["A3_n_noness"] = int(len(noness))

# ------------------------------------------------------------------ A4 epithelial / fitness
epi_markers = [g for g in ["EPCAM", "CDH1", "KRT8", "KRT18"] if g in dexpr.columns]
epc = dexpr.loc[score.index, epi_markers].mean(axis=1)
print(f"\n[A4] score vs epithelial marker score: rho={st.spearmanr(score, epc)[0]:.3f}")
epi_models = epc[epc > epc.median()].index.intersection(C.index)
sce = score.loc[epi_models]
scer = rankdata(sce.values); scer = (scer - scer.mean()) / scer.std()
mean_ge = C.mean(axis=1)
mr = rankdata(mean_ge.values); mr = (mr - mr.mean()) / mr.std()
sc_res = scr - np.polyval(np.polyfit(mr, scr, 1), mr)
print(f"[A4] per-line mean gene effect vs score: rho={st.spearmanr(mean_ge, sc)[0]:.3f}")
tbl = []
for g in AXIS:
    if g not in C.columns:
        continue
    v = C[g].values; m = ~np.isnan(v)
    r_pan = np.corrcoef(scr[m], rankdata(v[m]))[0, 1]
    r_fit = np.corrcoef(sc_res[m], rankdata(v[m]))[0, 1]
    ve = C.loc[epi_models, g].values; me = ~np.isnan(ve)
    r_epi = np.corrcoef(scer[me], rankdata(ve[me]))[0, 1]
    tbl.append(dict(gene=g, rho_pan=r_pan, rho_epithelial=r_epi, rho_fitness_adj=r_fit))
t4 = pd.DataFrame(tbl).set_index("gene")
print(f"[A4] epithelial-restricted n={len(epi_models)}   [manuscript 621]")
print(t4.round(3).to_string())
res["A4"] = dict(n_epi=int(len(epi_models)), rho_score_epi=float(st.spearmanr(score, epc)[0]),
                 rho_meanGE=float(st.spearmanr(mean_ge, sc)[0]), table=t4.to_dict("index"))


def partial_sp(x, y, cov):
    xr, yr, cr = rankdata(x), rankdata(y), rankdata(cov)
    rx = xr - np.polyval(np.polyfit(cr, xr, 1), cr)
    ry = yr - np.polyval(np.polyfit(cr, yr, 1), cr)
    return np.corrcoef(rx, ry)[0, 1]


print("\n[A4] partial correlations vs GPX4:")
for g in ["TXNRD1", "SEPSECS", "EEFSEC", "PSTK"]:
    v1, v2 = C[g].values, C["GPX4"].values
    m = ~np.isnan(v1) & ~np.isnan(v2)
    print(f"   {g}|GPX4={partial_sp(sc.values[m], v1[m], v2[m]):+.3f}   "
          f"GPX4|{g}={partial_sp(sc.values[m], v2[m], v1[m]):+.3f}")

# ================================================================ B1 lineage structure
print("\n==================== B1  WHAT IS THE SCORE ACROSS LINEAGES? ====================")
ls = pd.DataFrame({"score": score, "lineage": lin}).dropna()
agg = ls.groupby("lineage")["score"].agg(["count", "mean"]).sort_values("mean", ascending=False)
print(agg[agg["count"] >= 10].round(2).to_string())
kw = st.kruskal(*[g["score"].values for _, g in ls.groupby("lineage") if len(g) >= 5])
print(f"[B1] Kruskal-Wallis across lineages: H={kw.statistic:.0f}, p={kw.pvalue:.2e}")
eta = 1 - (ls.groupby("lineage")["score"].transform(lambda x: x - x.mean()).var() / ls["score"].var())
print(f"[B1] variance of the score explained by lineage alone = {eta*100:.1f}%")
res["B1"] = dict(kruskal_p=float(kw.pvalue), var_explained_by_lineage=float(eta),
                 lineage_means=agg[agg["count"] >= 10]["mean"].round(3).to_dict(),
                 prostate_rank=int(list(agg.index).index("Prostate")) + 1 if "Prostate" in agg.index else None,
                 n_lineages=int(ls["lineage"].nunique()))

# ================================================================ B2 within-lineage
print("\n==================== B2  WITHIN-LINEAGE (lineage-demeaned) ====================")
lin_c = lin.reindex(C.index)
keep_lin = lin_c.value_counts()
keep_lin = keep_lin[keep_lin >= 15].index
mask_l = lin_c.isin(keep_lin).values
sc_l = sc[mask_l]
lin_l = lin_c[mask_l]
# demean the score within lineage, and demean each gene's effect within lineage
sc_dm = sc_l.groupby(lin_l.values).transform(lambda x: x - x.mean())
rows = []
Cl = C.loc[sc_l.index]
for g in AXIS:
    if g not in Cl.columns:
        continue
    v = Cl[g]
    m = v.notna().values
    gdm = v[m].groupby(lin_l[m].values).transform(lambda x: x - x.mean())
    r = st.spearmanr(sc_dm[m], gdm)
    rows.append(dict(gene=g, rho_pan=float(dep.loc[g, "rho"]) if g in dep.index else np.nan,
                     rho_within_lineage=float(r[0]), p=float(r[1]), n=int(m.sum())))
t_dm = pd.DataFrame(rows).set_index("gene")
print(f"[B2] lineages with >=15 screened lines: {len(keep_lin)} covering {mask_l.sum()} models")
print(t_dm.round(4).to_string())
res["B2"] = dict(n_lineages=int(len(keep_lin)), n_models=int(mask_l.sum()),
                 table=t_dm.to_dict("index"))

# genome-wide within-lineage ranking for TXNRD1
print("[B2] genome-wide within-lineage recomputation (this takes a minute)...")
Vl = Cl.values
sc_dm_v = sc_dm.values
lin_codes = pd.Categorical(lin_l).codes
rows = []
for j, g in enumerate(Cl.columns):
    v = Vl[:, j]
    m = ~np.isnan(v)
    if m.sum() < 400:
        continue
    vv = v[m]
    lc = lin_codes[m]
    # demean within lineage
    s = np.bincount(lc, weights=vv); cnt = np.bincount(lc)
    vv_dm = vv - (s / np.maximum(cnt, 1))[lc]
    a = rankdata(sc_dm_v[m]); b = rankdata(vv_dm)
    rows.append((g, np.corrcoef(a, b)[0, 1], int(m.sum())))
dep_wl = pd.DataFrame(rows, columns=["gene", "rho_wl", "n"]).set_index("gene")
dep_wl["ce"] = dep_wl.index.isin(ce)
nn = dep_wl[~dep_wl["ce"]]
for g in ["TXNRD1", "SEPSECS", "EEFSEC", "PSTK", "GPX4", "GAPDH"]:
    if g in dep_wl.index:
        rk = int((nn["rho_wl"] < dep_wl.loc[g, "rho_wl"]).sum() + 1)
        print(f"   {g:8s} within-lineage rho={dep_wl.loc[g,'rho_wl']:+.3f}  rank {rk}/{len(nn)}  "
              f"(pan-lineage rank was {res.get('A3_rank_'+g,'?')})")
        res[f"B2_rank_{g}"] = rk
        res[f"B2_rho_{g}"] = float(dep_wl.loc[g, "rho_wl"])
dep_wl.to_csv(f"{OUT}/s05_dep_within_lineage.csv")

# ================================================================ B3 remove blood lineages
print("\n==================== B3  WITHOUT HAEMATOPOIETIC/LYMPHOID LINES ====================")
blood = lin_c.astype(str).str.contains("Lymphoid|Myeloid", na=False)
keep = (~blood).values
scb = sc[keep]
scbr = rankdata(scb.values); scbr = (scbr - scbr.mean()) / scbr.std()
Cb = C.loc[scb.index]
print(f"[B3] removed {int(blood.sum())} blood-lineage lines, n={len(scb)}")
rows = []
for g in AXIS:
    if g not in Cb.columns:
        continue
    v = Cb[g].values; m = ~np.isnan(v)
    rows.append(dict(gene=g, rho_pan=float(dep.loc[g, "rho"]),
                     rho_no_blood=float(np.corrcoef(scbr[m], rankdata(v[m]))[0, 1])))
t_b = pd.DataFrame(rows).set_index("gene")
print(t_b.round(3).to_string())
res["B3"] = dict(n_removed=int(blood.sum()), table=t_b.to_dict("index"))

# ================================================================ B4 the actual construct
print("\n==================== B4  THE ACTUAL CONSTRUCT: FOLH1 EXPRESSION ====================")
fo = dexpr.loc[C.index, "FOLH1"]
print(f"[B4] FOLH1 expression in DepMap models: median={fo.median():.2f}, "
      f"frac >1 log2TPM = {(fo>1).mean():.2f}, max={fo.max():.2f}")
rows = []
for g in ["TXNRD1", "SEPSECS", "EEFSEC", "PSTK", "GPX4", "FOXA1"]:
    if g not in C.columns:
        continue
    v = C[g]
    m = v.notna() & fo.notna()
    r_all = st.spearmanr(fo[m], v[m])
    expressed = m & (fo > 1)
    r_exp = st.spearmanr(fo[expressed], v[expressed]) if expressed.sum() > 30 else (np.nan, np.nan)
    rows.append(dict(gene=g, rho_FOLH1_vs_dependency=float(r_all[0]), p=float(r_all[1]),
                     n=int(m.sum()),
                     rho_in_FOLH1_expressing=float(r_exp[0]), n_expr=int(expressed.sum())))
t_f = pd.DataFrame(rows).set_index("gene")
print("[B4] direct correlation of FOLH1 mRNA with dependency (sign convention: negative rho")
print("     = more dependent when FOLH1 is HIGH; the manuscript's claim predicts POSITIVE rho)")
print(t_f.round(4).to_string())
res["B4"] = dict(folh1_median=float(fo.median()), frac_expressing=float((fo > 1).mean()),
                 table=t_f.to_dict("index"))

# ================================================================ B5 prostate only
print("\n==================== B5  PROSTATE LINEAGE ONLY ====================")
prc = [m for m in pr if m in C.index]
pp = pd.DataFrame({"line": names.reindex(prc), "score": score.reindex(prc),
                   "FOLH1": dexpr.loc[prc, "FOLH1"], "disease": dis.reindex(prc)})
for g in ["TXNRD1", "SEPSECS", "GPX4", "FOXA1"]:
    if g in C.columns:
        pp[g] = C.loc[prc, g].values
pp = pp.sort_values("score", ascending=False)
print(pp.round(2).to_string(index=False))
if len(prc) > 5:
    for g in ["TXNRD1", "SEPSECS", "GPX4"]:
        if g in pp.columns:
            r = st.spearmanr(pp["score"], pp[g], nan_policy="omit")
            r2 = st.spearmanr(pp["FOLH1"], pp[g], nan_policy="omit")
            print(f"[B5] prostate-only (n={pp[g].notna().sum()}): score vs {g} rho={r[0]:+.3f} "
                  f"(p={r[1]:.2f}) | FOLH1 vs {g} rho={r2[0]:+.3f} (p={r2[1]:.2f})")
            res[f"B5_{g}"] = dict(rho_score=float(r[0]), p_score=float(r[1]),
                                  rho_folh1=float(r2[0]), p_folh1=float(r2[1]))
res["B5_table"] = pp.reset_index().to_dict("records")

# ================================================================ B6 effect-size yardstick
print("\n==================== B6  EFFECT-SIZE YARDSTICK ====================")
print("    For reference: how strongly does a gene's OWN expression predict its OWN")
print("    dependency across the same 1,103 models? That is what a real selective")
print("    dependency looks like.")
bench = {}
for g in ["ERBB2", "EGFR", "BRAF", "KRAS", "SOX10", "MITF", "CTNNB1", "ESR1", "AR", "MDM2",
          "CCND1", "FOXA1", "TXNRD1", "GPX4", "SEPSECS"]:
    if g in C.columns and g in dexpr.columns:
        v = C[g]; e = dexpr.loc[C.index, g]
        m = v.notna() & e.notna()
        bench[g] = float(st.spearmanr(e[m], v[m])[0])
bs = pd.Series(bench).sort_values()
print(bs.round(3).to_string())
print(f"[B6] the manuscript's headline TXNRD1 signal is rho={dep.loc['TXNRD1','rho']:+.3f} "
      f"(R2={dep.loc['TXNRD1','rho']**2*100:.1f}% of variance)")
res["B6"] = dict(own_expression_vs_own_dependency=bench,
                 txnrd1_r2_pct=float(dep.loc["TXNRD1", "rho"] ** 2 * 100))

json.dump(res, open(f"{OUT}/s05_results.json", "w"), indent=1)
print("\nsaved s05 outputs")
