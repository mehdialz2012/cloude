"""Step 3 — the PSMA-low state: definition, signature, external validation, survival.

PART A (replication)
  A1  Q1 vs Q3 split in SU2C; liver-site enrichment; liver score correlation
  A2  liver-residualized DE; 1,002-gene set; raw-vs-adjusted concordance
  A3  100-gene signature; Beltran AUC; GSE126078 rho; GSEA; robustness controls
  A4  survival (log-rank, univariable + multivariable Cox)
PART B (independent validation)
  B1  TUMOUR COMPOSITION: epithelial / stromal / immune content scores vs FOLH1;
      how much of the signature is non-tumour cell content
  B2  IS THE SIGNATURE MORE THAN LINEAGE? head-to-head AUC/rho vs plain AR and NE scores
  B3  CIRCULARITY NULL: signatures derived the same way from random anchor genes
  B4  SITE/NE STRATIFICATION: does the state survive excluding liver biopsies and
      histologically NE tumours?
  B5  CONTINUOUS vs QUARTILE: is there a PSMA-low subpopulation at all?
"""
import json
import numpy as np
import pandas as pd
from scipy import stats as st
from scipy.stats import rankdata
from statsmodels.stats.multitest import multipletests
from sklearn.metrics import roc_auc_score
import statsmodels.api as sm
import decoupler as dc

from common import (DATA, OUT, AR20, NE12, LIVER10, STROMAL10,
                    load_su2c, load_beltran, load_gse126078, sig_score, zrow)

res = {}
su2c_log, su2c_cs, su2c_cp = load_su2c()
f = su2c_log.loc["FOLH1"]

# ---------------------------------------------------------------- A1
q25, q75 = f.quantile(.25), f.quantile(.75)
low, high = f[f <= q25].index, f[f >= q75].index
print(f"[A1] PSMA-low n={len(low)} (FOLH1<={q25:.2f}) vs PSMA-high n={len(high)} "
      f"(>={q75:.2f})   [manuscript 67/67]")
res["A1"] = dict(n_low=len(low), n_high=len(high), q25=float(q25), q75=float(q75))

site = su2c_cs["TISSUE_SITE"].reindex(f.index)
fq = pd.qcut(f, 4, labels=["Q1", "Q2", "Q3", "Q4"])
tab = pd.crosstab(fq, site, normalize="index").round(3)
print("\n[A1] tissue site by FOLH1 quartile (row %):")
print(tab.to_string())
liver_col = [c for c in tab.columns if "iver" in c]
if liver_col:
    print(f"[A1] liver fraction Q1={tab.loc['Q1', liver_col[0]]:.2f} vs "
          f"Q4={tab.loc['Q4', liver_col[0]]:.2f}   [manuscript 33% vs 1%]")
    res["A1_liver_q1"] = float(tab.loc["Q1", liver_col[0]])
    res["A1_liver_q4"] = float(tab.loc["Q4", liver_col[0]])

liver_score = zrow(su2c_log.loc[[g for g in LIVER10 if g in su2c_log.index]]).mean(axis=0)
r_lv, p_lv = st.spearmanr(f, liver_score)
print(f"[A1] FOLH1 vs liver score: rho={r_lv:.3f} (p={p_lv:.1e})   [manuscript -0.537]")
res["A1_rho_folh1_liver"] = float(r_lv)

# ---------------------------------------------------------------- A2
expressed = su2c_log.index[su2c_log.median(axis=1) >= 1.0]
print(f"\n[A2] expressed genes (median log2FPKM+1 >= 1): {len(expressed)}")


def residualize(mat, cov):
    """Row-wise OLS residualization on a single covariate (authors' np.polyfit approach)."""
    c = cov.loc[mat.columns].values
    X = np.vstack([np.ones_like(c), c]).T
    Y = mat.values
    beta, *_ = np.linalg.lstsq(X, Y.T, rcond=None)
    return pd.DataFrame(Y - (X @ beta).T, index=mat.index, columns=mat.columns)


def de_test(mat, lo, hi, genes):
    out = []
    for g in genes:
        a, b = mat.loc[g, lo].values, mat.loc[g, hi].values
        U, p = st.mannwhitneyu(a, b)
        out.append((g, 2 * U / (len(a) * len(b)) - 1, p))   # rbc>0 == higher in PSMA-low
    r = pd.DataFrame(out, columns=["gene", "rbc", "p"]).set_index("gene")
    r["padj"] = multipletests(r["p"], method="fdr_bh")[1]
    return r


Xr = residualize(su2c_log, liver_score)
de_raw = de_test(su2c_log, low, high, expressed)
de_adj = de_test(Xr, low, high, expressed)
r_c = st.spearmanr(de_raw["rbc"], de_adj["rbc"])[0]
sig_adj = de_adj[(de_adj["padj"] < .05) & (de_adj["rbc"].abs() > .4)]
print(f"[A2] raw vs liver-adjusted DE concordance rho={r_c:.3f}   [manuscript 0.955]")
print(f"[A2] adjusted DE genes: {len(sig_adj)} (up {int((sig_adj.rbc>.4).sum())}, "
      f"down {int((sig_adj.rbc<-.4).sum())})   [manuscript 1002 = 168 up + 834 down]")
res["A2"] = dict(n_expressed=int(len(expressed)), concordance=float(r_c), n_sig=int(len(sig_adj)),
                 n_up=int((sig_adj.rbc > .4).sum()), n_down=int((sig_adj.rbc < -.4).sum()))
de_adj.to_csv(f"{OUT}/s03_de_adj.csv")

up50 = sig_adj[sig_adj.rbc > 0].sort_values("rbc", ascending=False).head(50).index.tolist()
dn50 = sig_adj[sig_adj.rbc < 0].sort_values("rbc").head(50).index.tolist()
pd.DataFrame({"up50": up50, "dn50": dn50}).to_csv(f"{OUT}/s03_signature100.csv", index=False)
print(f"[A2] up50 head: {up50[:8]}")
print(f"[A2] dn50 head: {dn50[:8]}")
res["A2_up50"] = up50
res["A2_dn50"] = dn50


def psma_low_score(logmat, up, dn):
    z = zrow(logmat)
    ug = [g for g in up if g in logmat.index]
    dg = [g for g in dn if g in logmat.index]
    return z.loc[ug].mean(axis=0) - z.loc[dg].mean(axis=0), len(ug), len(dg)


# ---------------------------------------------------------------- A3 external validation
nepc_log, nepc_cs = load_beltran()
nepc_map = nepc_cs.set_index("OTHER_SAMPLE_ID")
sc_b, nu, nd = psma_low_score(nepc_log, up50, dn50)
# how many signature genes are destroyed by the negative-floor clipping?
const_genes = [g for g in up50 + dn50 if g in nepc_log.index and nepc_log.loc[g].std() == 0]
print(f"\n[A3] Beltran: signature genes matched up={nu}/50 dn={nd}/50; "
      f"constant-after-clipping (silently dropped) = {len(const_genes)}: {const_genes}")
res["A3_beltran_constant_genes"] = const_genes

bidx = sc_b.index.intersection(nepc_map.index)
y_ne = (nepc_map.loc[bidx, "CANCER_TYPE_DETAILED"] == "Prostate Neuroendocrine Carcinoma").astype(int)
auc = roc_auc_score(y_ne, sc_b[bidx])
p_auc = st.mannwhitneyu(sc_b[bidx][y_ne == 1], sc_b[bidx][y_ne == 0]).pvalue
sc_b_nf, _, _ = psma_low_score(nepc_log, up50, [g for g in dn50 if g != "FOLH1"])
auc_nf = roc_auc_score(y_ne, sc_b_nf[bidx])
print(f"[A3] Beltran AUC (NEPC vs CRPC-adeno) = {auc:.3f} (p={p_auc:.1e}); "
      f"FOLH1-free = {auc_nf:.3f}   [manuscript 0.82 / 0.814]")
res["A3_beltran_auc"] = float(auc)
res["A3_beltran_auc_nofolh1"] = float(auc_nf)
res["A3_beltran_p"] = float(p_auc)

gse_log = load_gse126078()
is_lucap = gse_log.columns.str.contains("LuCaP")
is_cell = gse_log.columns.str.contains("C42B|PC3|PacMet")
is_pat = ~is_lucap & ~is_cell
sc_g, _, _ = psma_low_score(gse_log, up50, dn50)
f_g = gse_log.loc["FOLH1"]
r_pool = st.spearmanr(sc_g, f_g)
r_pat = st.spearmanr(sc_g[is_pat], f_g[is_pat])
r_lu = st.spearmanr(sc_g[is_lucap], f_g[is_lucap])
sc_g_nf, _, _ = psma_low_score(gse_log, up50, [g for g in dn50 if g != "FOLH1"])
print(f"\n[A3] GSE126078 n={gse_log.shape[1]} (patient {is_pat.sum()}, LuCaP {is_lucap.sum()}, "
      f"cell lines {is_cell.sum()})   [manuscript 153 = 98/43/12]")
print(f"[A3] signature vs FOLH1: pooled rho={r_pool[0]:.3f} (p={r_pool[1]:.1e}); "
      f"patient-only {r_pat[0]:.3f}; LuCaP {r_lu[0]:.3f}; FOLH1-free pooled "
      f"{st.spearmanr(sc_g_nf, f_g)[0]:.3f}   [manuscript -0.757 / -0.824 / -0.753]")
res["A3_gse"] = dict(n=int(gse_log.shape[1]), n_pat=int(is_pat.sum()), n_lucap=int(is_lucap.sum()),
                     n_cell=int(is_cell.sum()), rho_pooled=float(r_pool[0]),
                     rho_patient=float(r_pat[0]), rho_lucap=float(r_lu[0]),
                     rho_nofolh1=float(st.spearmanr(sc_g_nf, f_g)[0]))
ne_g, _ = sig_score(gse_log, NE12)
ar_g, _ = sig_score(gse_log, AR20)
print(f"[A3] patient mCRPC: signature vs NE rho={st.spearmanr(sc_g[is_pat], ne_g[is_pat])[0]:.3f}, "
      f"vs AR rho={st.spearmanr(sc_g[is_pat], ar_g[is_pat])[0]:.3f}   [manuscript +0.470 / -0.879]")
res["A3_gse_ne"] = float(st.spearmanr(sc_g[is_pat], ne_g[is_pat])[0])
res["A3_gse_ar"] = float(st.spearmanr(sc_g[is_pat], ar_g[is_pat])[0])

# GSEA
hm = dc.op.hallmark(organism="human")
rank = de_adj["rbc"].sort_values(ascending=False)
sc_gsea, pv_gsea = dc.mt.gsea(data=rank.to_frame("rbc").T, net=hm, verbose=False)
g = sc_gsea.T.sort_values(by="rbc", ascending=False)
g.columns = ["NES"]
g["p"] = pv_gsea.T.loc[g.index].values
print("\n[A3] Hallmark GSEA top up / down:")
print(g.head(5).round(3).to_string())
print(g.tail(5).round(3).to_string())
res["A3_gsea_top"] = {k: float(v) for k, v in g["NES"].head(6).items()}
res["A3_gsea_bot"] = {k: float(v) for k, v in g["NES"].tail(6).items()}
g.to_csv(f"{OUT}/s03_gsea.csv")

# robustness
t1, t2 = f.quantile(1 / 3), f.quantile(2 / 3)
de_terc = de_test(Xr, f[f <= t1].index, f[f >= t2].index, de_adj.index)
print(f"\n[A3] tercile-vs-quartile DE concordance rho="
      f"{st.spearmanr(de_terc['rbc'], de_adj['rbc'])[0]:.3f}   [manuscript 0.966]")
sc_s, _, _ = psma_low_score(su2c_log, up50, dn50)
r_cont = st.spearmanr(sc_s, f)
print(f"[A3] signature vs continuous FOLH1 (SU2C) rho={r_cont[0]:.3f}   [manuscript -0.761]")
rng = np.random.default_rng(42)
nulls = []
allg = np.array(su2c_log.index)
for _ in range(200):
    s, _, _ = psma_low_score(su2c_log, rng.choice(allg, 50, replace=False).tolist(),
                             rng.choice(allg, 50, replace=False).tolist())
    nulls.append(st.spearmanr(s, f)[0])
print(f"[A3] shuffled-signature null |rho| 95th pct = {np.percentile(np.abs(nulls),95):.3f}  "
      f"[manuscript 0.301]")
res["A3_robust"] = dict(tercile_conc=float(st.spearmanr(de_terc['rbc'], de_adj['rbc'])[0]),
                        rho_continuous=float(r_cont[0]),
                        null95=float(np.percentile(np.abs(nulls), 95)))

# ---------------------------------------------------------------- A4 survival
from lifelines import CoxPHFitter
from lifelines.statistics import multivariate_logrank_test
surv = pd.DataFrame({"FOLH1": f}).join(su2c_cs[["PATIENT_ID"]], how="left")
surv = surv.join(su2c_cp[["OS_STATUS", "OS_MONTHS", "AGE_AT_DIAGNOSIS", "PSA"]], on="PATIENT_ID")
surv = surv.dropna(subset=["OS_STATUS", "OS_MONTHS"])
surv["event"] = (surv["OS_STATUS"] == "1:DECEASED").astype(int)
surv["OS_MONTHS"] = pd.to_numeric(surv["OS_MONTHS"], errors="coerce")
surv = surv.dropna(subset=["OS_MONTHS"])
surv["q"] = pd.qcut(surv["FOLH1"], 4, labels=["Q1", "Q2", "Q3", "Q4"])
lr = multivariate_logrank_test(surv["OS_MONTHS"], surv["q"], surv["event"])
print(f"\n[A4] survival set n={len(surv)}, events={int(surv.event.sum())}; "
      f"log-rank across quartiles p={lr.p_value:.3f}   [manuscript n=81, p=0.89]")
surv["FOLH1_z"] = (surv.FOLH1 - surv.FOLH1.mean()) / surv.FOLH1.std()
cu = CoxPHFitter().fit(surv[["OS_MONTHS", "event", "FOLH1_z"]], "OS_MONTHS", "event")
hr = float(np.exp(cu.params_["FOLH1_z"])); ci = np.exp(cu.confidence_intervals_.loc["FOLH1_z"].values)
print(f"[A4] univariable Cox HR per SD = {hr:.2f} ({ci[0]:.2f}-{ci[1]:.2f}), "
      f"p={cu.summary.loc['FOLH1_z','p']:.2f}   [manuscript 1.03 (0.77-1.38)]")
s2 = surv.dropna(subset=["AGE_AT_DIAGNOSIS"]).copy()
s2["PSA"] = pd.to_numeric(s2["PSA"], errors="coerce")
n_imp = int(s2["PSA"].isna().sum())
s2["PSA_log"] = np.log10(s2["PSA"].clip(lower=.1))
s2["PSA_log"] = s2["PSA_log"].fillna(s2["PSA_log"].median())
cm = CoxPHFitter().fit(s2[["OS_MONTHS", "event", "FOLH1_z", "AGE_AT_DIAGNOSIS", "PSA_log"]],
                       "OS_MONTHS", "event")
hrm = float(np.exp(cm.params_["FOLH1_z"])); cim = np.exp(cm.confidence_intervals_.loc["FOLH1_z"].values)
print(f"[A4] multivariable Cox n={len(s2)} (PSA imputed {n_imp}={100*n_imp/len(s2):.0f}%): "
      f"HR={hrm:.2f} ({cim[0]:.2f}-{cim[1]:.2f})   [manuscript n=69, 1.06 (0.76-1.47), 13%]")
# the signature itself as a survival predictor (not reported by the authors)
surv["sig_z"] = ((sc_s - sc_s.mean()) / sc_s.std()).reindex(surv.index)
cs = CoxPHFitter().fit(surv[["OS_MONTHS", "event", "sig_z"]].dropna(), "OS_MONTHS", "event")
print(f"[A4] univariable Cox on the PSMA-low SIGNATURE: HR per SD = "
      f"{np.exp(cs.params_['sig_z']):.2f} "
      f"({np.exp(cs.confidence_intervals_.loc['sig_z'].values)[0]:.2f}-"
      f"{np.exp(cs.confidence_intervals_.loc['sig_z'].values)[1]:.2f}), "
      f"p={cs.summary.loc['sig_z','p']:.2f}   [not reported in manuscript]")
res["A4"] = dict(n=len(surv), events=int(surv.event.sum()), logrank_p=float(lr.p_value),
                 hr_uni=hr, ci_uni=[float(ci[0]), float(ci[1])], n_mv=len(s2),
                 psa_imputed_frac=float(n_imp / len(s2)), hr_mv=hrm,
                 ci_mv=[float(cim[0]), float(cim[1])],
                 hr_signature=float(np.exp(cs.params_["sig_z"])),
                 p_signature=float(cs.summary.loc["sig_z", "p"]))

# ================================================================ B1 tumour composition
print("\n==================== B1  TUMOUR-COMPOSITION CONFOUND ====================")
EPI = ["EPCAM", "CDH1", "KRT8", "KRT18", "KRT18", "CLDN3", "CLDN4"]
IMM = ["PTPRC", "CD3E", "CD2", "CD53", "LCP1", "CSF1R", "AIF1", "ITGAM", "CD68", "TYROBP"]
comp = pd.DataFrame({
    "epithelial": sig_score(su2c_log, EPI)[0],
    "stromal": sig_score(su2c_log, STROMAL10)[0],
    "immune": sig_score(su2c_log, IMM)[0],
    "liver": liver_score,
    "AR": sig_score(su2c_log, AR20)[0],
    "NE": sig_score(su2c_log, NE12)[0],
    "FOLH1": f,
    "sig": sc_s})
for c in ["epithelial", "stromal", "immune", "liver"]:
    r_f = st.spearmanr(comp["FOLH1"], comp[c])
    r_s = st.spearmanr(comp["sig"], comp[c])
    print(f"[B1] {c:11s}: vs FOLH1 rho={r_f[0]:+.3f} (p={r_f[1]:.1e}) | "
          f"vs PSMA-low score rho={r_s[0]:+.3f} (p={r_s[1]:.1e})")
    res[f"B1_{c}"] = dict(rho_folh1=float(r_f[0]), p_folh1=float(r_f[1]),
                          rho_sig=float(r_s[0]), p_sig=float(r_s[1]))
print(f"[B1] epithelial content: PSMA-low mean={comp.loc[low,'epithelial'].mean():+.2f} vs "
      f"PSMA-high mean={comp.loc[high,'epithelial'].mean():+.2f} "
      f"(MWU p={st.mannwhitneyu(comp.loc[low,'epithelial'], comp.loc[high,'epithelial']).pvalue:.1e})")
res["B1_epi_low_mean"] = float(comp.loc[low, "epithelial"].mean())
res["B1_epi_high_mean"] = float(comp.loc[high, "epithelial"].mean())

# how many signature "up" genes are immune/stromal-associated?
up_imm = {g: float(st.spearmanr(su2c_log.loc[g], comp["immune"])[0]) for g in up50}
up_str = {g: float(st.spearmanr(su2c_log.loc[g], comp["stromal"])[0]) for g in up50}
n_imm = sum(1 for v in up_imm.values() if v > .5)
n_str = sum(1 for v in up_str.values() if v > .5)
print(f"[B1] of the 50 'up in PSMA-low' genes: {n_imm} correlate rho>0.5 with the immune score, "
      f"{n_str} with the stromal score")
res["B1_up50_immune_gt05"] = n_imm
res["B1_up50_stromal_gt05"] = n_str
res["B1_up50_immune_rho"] = up_imm
res["B1_up50_stromal_rho"] = up_str

# DE adjusted for epithelial content instead of liver
Xe = residualize(su2c_log, comp["epithelial"])
de_epi = de_test(Xe, low, high, pd.Index(up50 + dn50).intersection(expressed))
kept = de_epi[(de_epi.padj < .05) & (de_epi.rbc.abs() > .4)]
print(f"[B1] of the 100 signature genes, {len(kept)} remain FDR<0.05 & |rbc|>0.4 after "
      f"residualizing on EPITHELIAL content (up: "
      f"{int((kept.rbc>0).sum())}, down: {int((kept.rbc<0).sum())})")
res["B1_sig_genes_surviving_epi_adj"] = dict(n=int(len(kept)), up=int((kept.rbc > 0).sum()),
                                             down=int((kept.rbc < 0).sum()))
de_epi.to_csv(f"{OUT}/s03_de_epi_adjusted.csv")
comp.to_csv(f"{OUT}/s03_composition.csv")

# candidate surface targets vs composition
CAND = ["SLC16A3", "SLC7A5", "PLAUR", "GPC2", "ATP2B4", "SLC12A4", "ADAM19", "MRC2",
        "MUC1", "NCAM1", "SSTR2", "CEACAM5", "DLL3", "STEAP1", "STEAP2", "CD276", "TACSTD2"]
rows = []
for gg in CAND:
    if gg not in su2c_log.index:
        continue
    rows.append(dict(gene=gg, rbc_liver_adj=float(de_adj["rbc"].get(gg, np.nan)),
                     rho_immune=float(st.spearmanr(su2c_log.loc[gg], comp["immune"])[0]),
                     rho_stromal=float(st.spearmanr(su2c_log.loc[gg], comp["stromal"])[0]),
                     rho_epithelial=float(st.spearmanr(su2c_log.loc[gg], comp["epithelial"])[0]),
                     rbc_epi_adj=float(de_epi["rbc"].get(gg, np.nan))))
cand_comp = pd.DataFrame(rows).set_index("gene")
print("\n[B1] candidate targets vs cell-composition scores:")
print(cand_comp.round(3).to_string())
cand_comp.to_csv(f"{OUT}/s03_candidates_composition.csv")
res["B1_candidates"] = cand_comp.to_dict("index")

# ================================================================ B2 lineage head-to-head
print("\n==================== B2  SIGNATURE vs PLAIN LINEAGE SCORES ====================")
ne_b, _ = sig_score(nepc_log, NE12)
ar_b, _ = sig_score(nepc_log, AR20)
for nm, v in [("PSMA-low signature (100 genes)", sc_b), ("NE score (12 genes)", ne_b),
              ("AR score (20 genes, inverted)", -ar_b), ("FOLH1 alone (inverted)", -nepc_log.loc["FOLH1"])]:
    a = roc_auc_score(y_ne, v[bidx])
    print(f"[B2] Beltran NEPC-vs-adeno AUC | {nm:32s} = {a:.3f}")
    res[f"B2_beltran_auc::{nm}"] = float(a)
for nm, v in [("PSMA-low signature", sc_g), ("NE score", ne_g), ("AR score (inverted)", -ar_g)]:
    r = st.spearmanr(v[is_pat], f_g[is_pat])
    print(f"[B2] GSE126078 patient-only rho vs FOLH1 | {nm:22s} = {r[0]:+.3f} (p={r[1]:.1e})")
    res[f"B2_gse_rho::{nm}"] = float(r[0])
# partial correlations in the patient subset
dpat = pd.DataFrame({"sig": sc_g[is_pat], "folh1": f_g[is_pat], "AR": ar_g[is_pat], "NE": ne_g[is_pat]})
rs = sm.OLS(rankdata(dpat.sig), sm.add_constant(np.c_[rankdata(dpat.AR), rankdata(dpat.NE)])).fit().resid
rf = sm.OLS(rankdata(dpat.folh1), sm.add_constant(np.c_[rankdata(dpat.AR), rankdata(dpat.NE)])).fit().resid
print(f"[B2] GSE126078 patient-only: signature vs FOLH1 partial rho controlling AR+NE = "
      f"{st.pearsonr(rs, rf)[0]:+.3f}")
res["B2_gse_partial_rho"] = float(st.pearsonr(rs, rf)[0])

# ================================================================ B3 circularity null
print("\n==================== B3  'EXTERNAL VALIDATION' CIRCULARITY NULL ====================")
print("    A signature derived from the quartiles of ANY gene will anti-correlate with that")
print("    gene in another cohort. Repeat the whole pipeline for random anchor genes.")
rng = np.random.default_rng(7)
pool = [g for g in expressed if g in gse_log.index and g != "FOLH1"
        and su2c_log.loc[g].std() > 0.5]
anchors = rng.choice(pool, 15, replace=False)
null_rows = []
for ag in anchors:
    fa = su2c_log.loc[ag]
    lo, hi = fa[fa <= fa.quantile(.25)].index, fa[fa >= fa.quantile(.75)].index
    d = de_test(su2c_log, lo, hi, expressed)
    s = d[(d.padj < .05) & (d.rbc.abs() > .4)]
    if (s.rbc > 0).sum() < 20 or (s.rbc < 0).sum() < 20:
        continue
    u = s[s.rbc > 0].sort_values("rbc", ascending=False).head(50).index.tolist()
    dd = s[s.rbc < 0].sort_values("rbc").head(50).index.tolist()
    sg, _, _ = psma_low_score(gse_log, u, dd)
    r_ext = st.spearmanr(sg[is_pat], gse_log.loc[ag][is_pat])[0]
    sb, _, _ = psma_low_score(nepc_log, u, dd)
    a_ne = roc_auc_score(y_ne, sb[bidx])
    null_rows.append(dict(anchor=ag, n_sig=len(s), rho_gse_patient=r_ext, beltran_auc=a_ne))
    print(f"   anchor {ag:10s}: {len(s):5d} DE genes | GSE126078 patient rho={r_ext:+.3f} | "
          f"Beltran NEPC AUC={a_ne:.3f}")
nullc = pd.DataFrame(null_rows)
print(f"[B3] random-anchor signatures: median |rho| vs own anchor in GSE126078 = "
      f"{nullc['rho_gse_patient'].abs().median():.3f} (FOLH1-anchored: {abs(r_pat[0]):.3f})")
print(f"[B3] random-anchor signatures: Beltran NEPC AUC median = "
      f"{nullc['beltran_auc'].median():.3f}, max = {nullc['beltran_auc'].max():.3f} "
      f"(FOLH1-anchored: {auc:.3f})")
nullc.to_csv(f"{OUT}/s03_circularity_null.csv", index=False)
res["B3"] = dict(median_abs_rho=float(nullc["rho_gse_patient"].abs().median()),
                 median_auc=float(nullc["beltran_auc"].median()),
                 max_auc=float(nullc["beltran_auc"].max()),
                 table=nullc.to_dict("records"))

# ================================================================ B4 site / NE stratification
print("\n==================== B4  STRATIFIED RE-DERIVATION ====================")
ne_feat = su2c_cs["NEUROENDOCRINE_FEATURES"].reindex(f.index)
keep = f.index[(~site.astype(str).str.contains("Liver", na=False)) & (ne_feat != "Yes")]
fk = f.loc[keep]
lok, hik = fk[fk <= fk.quantile(.25)].index, fk[fk >= fk.quantile(.75)].index
print(f"[B4] excluding liver biopsies and NE-feature tumours: n={len(keep)} "
      f"(low {len(lok)}, high {len(hik)})")
de_k = de_test(su2c_log.loc[:, keep], lok, hik, expressed)
sig_k = de_k[(de_k.padj < .05) & (de_k.rbc.abs() > .4)]
ov = len(set(sig_k.index) & set(sig_adj.index))
print(f"[B4] DE genes in the restricted cohort: {len(sig_k)} "
      f"(up {int((sig_k.rbc>.4).sum())}, down {int((sig_k.rbc<-.4).sum())}); "
      f"overlap with the 1,002-gene set = {ov}")
print(f"[B4] rbc concordance (restricted vs full, all expressed genes) rho="
      f"{st.spearmanr(de_k['rbc'], de_adj['rbc'])[0]:.3f}")
# does the published signature still separate NEPC after removing NE genes from it?
sig_ne_members = [g for g in up50 + dn50 if g in NE12 + AR20]
print(f"[B4] signature genes that are themselves AR/NE panel members: {sig_ne_members}")
res["B4"] = dict(n_restricted=int(len(keep)), n_sig=int(len(sig_k)),
                 up=int((sig_k.rbc > .4).sum()), down=int((sig_k.rbc < -.4).sum()),
                 overlap_with_1002=int(ov),
                 rbc_concordance=float(st.spearmanr(de_k["rbc"], de_adj["rbc"])[0]),
                 sig_genes_in_AR_NE_panels=sig_ne_members)
de_k.to_csv(f"{OUT}/s03_de_restricted.csv")

# ================================================================ B5 continuous structure
print("\n==================== B5  IS THERE A PSMA-LOW SUBPOPULATION? ====================")
print(f"[B5] SU2C FOLH1 log2(FPKM+1): min={f.min():.2f} Q1={q25:.2f} median={f.median():.2f} "
      f"Q3={q75:.2f} max={f.max():.2f}")
print(f"[B5] Hartigan-style check — bimodality coefficient = "
      f"{(st.skew(f)**2 + 1)/st.kurtosis(f, fisher=False):.3f} (>0.555 suggests bimodality)")
dip = st.shapiro(f)
print(f"[B5] Shapiro-Wilk normality p={dip.pvalue:.1e}")
print(f"[B5] samples within +/-0.5 SD of the median: {int(((f-f.median()).abs()<.5*f.std()).sum())}"
      f"/{len(f)}; Q1/Q3 boundary gap = {q75-q25:.2f} log2 units "
      f"({2**(q75-q25):.1f}-fold in FPKM)")
res["B5"] = dict(min=float(f.min()), q1=float(q25), median=float(f.median()), q3=float(f.quantile(.75)),
                 max=float(f.max()),
                 bimodality_coef=float((st.skew(f)**2 + 1) / st.kurtosis(f, fisher=False)),
                 shapiro_p=float(dip.pvalue), fold_q1_q3=float(2**(q75 - q25)))

json.dump(res, open(f"{OUT}/s03_results.json", "w"), indent=1)
print("\nsaved s03 outputs")
