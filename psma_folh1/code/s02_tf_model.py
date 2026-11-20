"""Step 2 — the transcription-factor model of FOLH1 (manuscript Results section 2).

PART A (replication, same recipe as the authors' notebook)
  A1  CollecTRI/decoupleR-ULM TF activities in TCGA (all 550) and SU2C
  A2  TF-FOLH1 Spearman associations, cross-cohort concordance, sign counts, Fisher z
  A3  Model A (TF + 36 locus methylation probes) nested 5x5 CV in TCGA; TF-only comparison
  A4  Model B (TF-only) TCGA -> SU2C in z-space, bootstrap CI, AR-only baseline
PART B (independent validation tests the authors did not run)
  B1  TARGET LEAKAGE: FOLH1 is a member of the CollecTRI regulons of ELF5/ELF3/AR/VDR.
      Recompute everything with FOLH1 removed from every regulon.
  B2  NULL MODEL: size-matched random gene sets as pseudo-regulons -> same pipeline.
  B3  PLAIN-TRANSCRIPTOME BASELINE: elastic net on 766 random genes' expression.
  B4  LINEAGE-ONLY BASELINE: AR score, NE score, AR+NE.
  B5  CLINICAL UTILITY: does the model identify the PSMA-low quartile in SU2C?
"""
import json
import numpy as np
import pandas as pd
from scipy import stats as st
from statsmodels.stats.multitest import multipletests
from sklearn.linear_model import ElasticNetCV, LinearRegression
from sklearn.model_selection import KFold
from sklearn.metrics import r2_score, roc_auc_score
from sklearn.preprocessing import StandardScaler
import decoupler as dc

from common import (DATA, OUT, AR20, NE12, load_tcga, load_su2c, sig_score)

RNG = np.random.default_rng(1)
res = {}

# ------------------------------------------------------------------ data
tcga_log, _ = load_tcga()
su2c_log, su2c_cs, _ = load_su2c()
net = pd.read_csv(f"{DATA}/collectri.csv")
print(f"CollecTRI: {net.shape[0]} edges, {net['source'].nunique()} TFs")
folh1_tfs = net.loc[net["target"] == "FOLH1", "source"].tolist()
print(f"TFs whose CollecTRI regulon CONTAINS FOLH1: {folh1_tfs}")
for tf in folh1_tfs:
    n = (net["source"] == tf).sum()
    print(f"   {tf}: regulon size {n} -> FOLH1 is 1/{n} = {100/n:.1f}% of the set")
res["folh1_in_regulons"] = {tf: int((net['source'] == tf).sum()) for tf in folh1_tfs}

y_tcga = tcga_log.loc["FOLH1"]
y_su2c = su2c_log.loc["FOLH1"]


def tf_act(logmat, network):
    est, _ = dc.mt.ulm(data=logmat.T, net=network, verbose=False)
    return est


def assoc(acts, y):
    common = acts.index.intersection(y.index)
    out = [(tf, *st.spearmanr(acts.loc[common, tf], y[common])) for tf in acts.columns]
    r = pd.DataFrame(out, columns=["TF", "rho", "p"]).set_index("TF")
    r["padj"] = multipletests(r["p"], method="fdr_bh")[1]
    return r


def fisher_z(r1, n1, r2, n2):
    z = (np.arctanh(r1) - np.arctanh(r2)) / np.sqrt(1 / (n1 - 3) + 1 / (n2 - 3))
    return z, 2 * st.norm.sf(abs(z))


def zdf(df):
    return pd.DataFrame(StandardScaler().fit_transform(df), index=df.index, columns=df.columns)


def modelB(Xtr, ytr, Xte, yte, seed=1):
    """Elastic net trained on standardized TCGA, tested on standardized SU2C."""
    m = ElasticNetCV(l1_ratio=[0.5, 0.7, 0.9, 1.0], alphas=50, cv=5,
                     random_state=seed, max_iter=20000, n_jobs=4)
    m.fit(Xtr, ytr)
    pred = m.predict(Xte)
    return m, pred, r2_score(yte, pred), st.spearmanr(yte, pred)[0]


# ============================================================ A1/A2
tf_tcga = tf_act(tcga_log, net)
tf_su2c = tf_act(su2c_log, net)
print(f"\nTF activities: TCGA {tf_tcga.shape}, SU2C {tf_su2c.shape}")
a_t, a_s = assoc(tf_tcga, y_tcga), assoc(tf_su2c, y_su2c)
rep = a_t[["rho", "padj"]].join(a_s[["rho", "padj"]], lsuffix="_TCGA", rsuffix="_SU2C").dropna()
r_rep, p_rep = st.spearmanr(rep["rho_TCGA"], rep["rho_SU2C"])
same = np.sign(rep["rho_TCGA"]) == np.sign(rep["rho_SU2C"])
sig_both = (rep["padj_TCGA"] < .05) & (rep["padj_SU2C"] < .05)
print(f"\n[A2] n TFs compared = {len(rep)}")
print(f"[A2] cross-cohort concordance rho = {r_rep:.3f} (p={p_rep:.1e})   "
      f"[manuscript 0.502]")
print(f"[A2] significant & same sign in both = {int((sig_both & same).sum())}/{len(rep)}  "
      f"[manuscript 367/766]")
print(f"[A2] same sign regardless of significance = {int(same.sum())}/{len(rep)}  "
      f"[manuscript 594/766]")
print(f"[A2] significant in both, OPPOSITE sign = {int((sig_both & ~same).sum())}  [manuscript 39]")
res["A2"] = dict(n_tf=int(len(rep)), concordance_rho=float(r_rep),
                 sig_same=int((sig_both & same).sum()), same_sign=int(same.sum()),
                 sig_opposite=int((sig_both & ~same).sum()))
for tf in ["AR", "FOXA1", "MYCN", "EZH2", "HOXB13", "ELF5", "ELF3", "VDR"]:
    if tf in rep.index:
        z, p = fisher_z(rep.loc[tf, "rho_TCGA"], len(y_tcga), rep.loc[tf, "rho_SU2C"], len(y_su2c))
        print(f"   {tf:7s} TCGA rho={rep.loc[tf,'rho_TCGA']:+.3f}  SU2C rho={rep.loc[tf,'rho_SU2C']:+.3f}"
              f"  Fisher z={z:+.2f} p={p:.1e}")
        res[f"A2_{tf}"] = dict(rho_tcga=float(rep.loc[tf, "rho_TCGA"]),
                               rho_su2c=float(rep.loc[tf, "rho_SU2C"]), z=float(z), p=float(p))
rep.to_csv(f"{OUT}/s02_tf_assoc.csv")

# ============================================================ A3 Model A (TF + methylation)
meth = pd.read_csv(f"{DATA}/folh1_locus_meth.csv", index_col=0)   # probes x samples, built by s00
meth_tum = [c for c in meth.columns if c.endswith("-01")]
meth_p = meth[meth_tum]
meth_p.columns = [c[:12] for c in meth_tum]
meth_p = meth_p.loc[:, ~meth_p.columns.duplicated()]

tf_p = tf_tcga.copy()
tf_p.index = [c[:12] for c in tf_p.index]
tf_p = tf_p[~tf_p.index.duplicated(keep="first")]
yp = y_tcga.copy()
yp.index = [c[:12] for c in y_tcga.index]
yp = yp[~yp.index.duplicated(keep="first")]
commonA = tf_p.index.intersection(meth_p.columns).intersection(yp.index)
print(f"\n[A3] paired methylation/expression patients: {len(commonA)}  [manuscript 497]")

XA = pd.concat([tf_p.loc[commonA], meth_p[commonA].T.add_prefix("meth_")], axis=1)
XA = XA.loc[:, XA.notna().sum() > 0.9 * len(XA)].fillna(XA.median())
yA = yp.loc[commonA]
XA_z = zdf(XA)
yA_z = (yA - yA.mean()) / yA.std()
kf = KFold(n_splits=5, shuffle=True, random_state=1)
n_tf = tf_tcga.shape[1]


def nested_cv(X, y):
    pred = pd.Series(index=X.index, dtype=float)
    for tr, te in kf.split(X):
        m = ElasticNetCV(l1_ratio=[0.5, 0.7, 0.9, 1.0], alphas=50, cv=5,
                         random_state=1, max_iter=20000, n_jobs=4)
        m.fit(X.iloc[tr], y.iloc[tr])
        pred.iloc[te] = m.predict(X.iloc[te])
    return r2_score(y, pred), pred


r2A, _ = nested_cv(XA_z, yA_z)
r2A_tf, _ = nested_cv(XA_z.iloc[:, :n_tf], yA_z)
meth_cols = [c for c in XA_z.columns if c.startswith("meth_")]
r2A_meth, _ = nested_cv(XA_z[meth_cols], yA_z)
print(f"[A3] Model A  (TF+methylation) nested-CV R2 = {r2A:.3f}   [manuscript 0.717]")
print(f"[A3] Model A  (TF only)        nested-CV R2 = {r2A_tf:.3f}   [manuscript 0.712]")
print(f"[A3] methylation-only          nested-CV R2 = {r2A_meth:.3f}  [not reported]")
res["A3"] = dict(n=len(commonA), r2_tf_meth=float(r2A), r2_tf=float(r2A_tf),
                 r2_meth_only=float(r2A_meth))

# ============================================================ A4 Model B external
common_tfs = tf_tcga.columns.intersection(tf_su2c.columns)
Xtr, Xte = zdf(tf_tcga[common_tfs]), zdf(tf_su2c[common_tfs])
ytr = (y_tcga - y_tcga.mean()) / y_tcga.std()
yte = (y_su2c - y_su2c.mean()) / y_su2c.std()
mB, predB, r2B, rhoB = modelB(Xtr, ytr.loc[Xtr.index], Xte, yte.loc[Xte.index])
boots = [r2_score(yte.values[i], predB[i]) for i in
         (RNG.integers(0, len(yte), len(yte)) for _ in range(1000))]
print(f"\n[A4] Model B external R2 = {r2B:.3f}  95% CI [{np.percentile(boots,2.5):.3f},"
      f" {np.percentile(boots,97.5):.3f}]  rho={rhoB:.3f}   [manuscript R2 0.662 "
      f"CI 0.559-0.736 rho 0.802]")
# raw-space transfer (manuscript: R2 < 0)
mB_raw, pred_raw, r2B_raw, _ = modelB(tf_tcga[common_tfs], y_tcga, tf_su2c[common_tfs], y_su2c)
print(f"[A4] Model B raw-space external R2 = {r2B_raw:.2f}   [manuscript 'R2<0' / naive -12]")
base = LinearRegression().fit(zdf(pd.DataFrame({"AR": sig_score(tcga_log, AR20)[0]})), ytr)
r2_ar = r2_score(yte, base.predict(zdf(pd.DataFrame({"AR": sig_score(su2c_log, AR20)[0]}))))
print(f"[A4] AR-score-only baseline external R2 = {r2_ar:.3f}   [manuscript 0.382]")
sel = pd.Series(mB.coef_, index=common_tfs)
sel = sel[sel != 0].sort_values()
print(f"[A4] non-zero TFs: {len(sel)}; top +: {sel.tail(5).round(3).to_dict()}")
print(f"[A4]                      top -: {sel.head(5).round(3).to_dict()}")
res["A4"] = dict(r2_ext=float(r2B), ci=[float(np.percentile(boots, 2.5)), float(np.percentile(boots, 97.5))],
                 rho=float(rhoB), r2_raw=float(r2B_raw), r2_ar_only=float(r2_ar),
                 n_nonzero=int(len(sel)),
                 top_pos={k: float(v) for k, v in sel.tail(5).items()},
                 top_neg={k: float(v) for k, v in sel.head(5).items()})
sel.to_csv(f"{OUT}/s02_modelB_weights.csv")

# ============================================================ B1 leakage: drop FOLH1 from regulons
net_nl = net[net["target"] != "FOLH1"].copy()
tf_tcga_nl, tf_su2c_nl = tf_act(tcga_log, net_nl), tf_act(su2c_log, net_nl)
ct_nl = tf_tcga_nl.columns.intersection(tf_su2c_nl.columns)
mB2, predB2, r2B2, rhoB2 = modelB(zdf(tf_tcga_nl[ct_nl]), ytr.loc[tf_tcga_nl.index],
                                  zdf(tf_su2c_nl[ct_nl]), yte.loc[tf_su2c_nl.index])
XA2 = pd.concat([tf_tcga_nl.copy().set_axis([c[:12] for c in tf_tcga_nl.index])
                 .pipe(lambda d: d[~d.index.duplicated(keep="first")]).loc[commonA],
                 meth_p[commonA].T.add_prefix("meth_")], axis=1)
XA2 = XA2.loc[:, XA2.notna().sum() > 0.9 * len(XA2)].fillna(XA2.median())
r2A2, _ = nested_cv(zdf(XA2), yA_z)
a_t_nl = assoc(tf_tcga_nl, y_tcga)
a_s_nl = assoc(tf_su2c_nl, y_su2c)
print("\n==================== B1  TARGET-LEAKAGE TEST ====================")
print(f"[B1] Model A nested-CV R2:  with FOLH1 in regulons {r2A:.3f}  ->  without {r2A2:.3f}")
print(f"[B1] Model B external  R2:  with FOLH1 in regulons {r2B:.3f}  ->  without {r2B2:.3f}")
print(f"[B1] Model B external rho:  {rhoB:.3f} -> {rhoB2:.3f}")
sel2 = pd.Series(mB2.coef_, index=ct_nl)
sel2 = sel2[sel2 != 0].sort_values()
print(f"[B1] ELF5 weight: {sel.get('ELF5', 0.0):+.3f} -> {sel2.get('ELF5', 0.0):+.3f}")
for tf in folh1_tfs:
    if tf in a_t.index and tf in a_t_nl.index:
        print(f"[B1] {tf:5s} rho(FOLH1) TCGA {a_t.loc[tf,'rho']:+.3f} -> {a_t_nl.loc[tf,'rho']:+.3f}"
              f" | SU2C {a_s.loc[tf,'rho']:+.3f} -> {a_s_nl.loc[tf,'rho']:+.3f}")
res["B1"] = dict(r2A_noleak=float(r2A2), r2B_noleak=float(r2B2), rhoB_noleak=float(rhoB2),
                 elf5_w=float(sel.get("ELF5", 0.0)), elf5_w_noleak=float(sel2.get("ELF5", 0.0)),
                 per_tf={tf: dict(tcga=float(a_t.loc[tf, "rho"]), tcga_nl=float(a_t_nl.loc[tf, "rho"]),
                                  su2c=float(a_s.loc[tf, "rho"]), su2c_nl=float(a_s_nl.loc[tf, "rho"]))
                         for tf in folh1_tfs if tf in a_t.index})

# ============================================================ B2 random-regulon null
sizes = net_nl.groupby("source").size()
sizes = sizes[sizes >= 5]
genes_pool = np.array([g for g in tcga_log.index.intersection(su2c_log.index) if g != "FOLH1"])
N_NULL = 15
null_r2, null_rho = [], []
print(f"\n==================== B2  RANDOM-REGULON NULL (n={N_NULL}) ====================")
for it in range(N_NULL):
    rows = []
    for i, (tf, k) in enumerate(sizes.items()):
        tg = RNG.choice(genes_pool, size=int(min(k, len(genes_pool))), replace=False)
        rows.append(pd.DataFrame({"source": f"R{i}", "target": tg, "weight": 1.0}))
    rnet = pd.concat(rows, ignore_index=True)
    at, asu = tf_act(tcga_log, rnet), tf_act(su2c_log, rnet)
    cc = at.columns.intersection(asu.columns)
    _, _, r2n, rhon = modelB(zdf(at[cc]), ytr.loc[at.index], zdf(asu[cc]), yte.loc[asu.index])
    null_r2.append(r2n)
    null_rho.append(rhon)
    print(f"   null {it+1:2d}: external R2={r2n:.3f}  rho={rhon:.3f}")
print(f"[B2] random-regulon null external R2: mean {np.mean(null_r2):.3f}, "
      f"median {np.median(null_r2):.3f}, max {np.max(null_r2):.3f}  | real TF model {r2B:.3f}")
print(f"[B2] random-regulon null rho:         mean {np.mean(null_rho):.3f}, "
      f"max {np.max(null_rho):.3f}  | real {rhoB:.3f}")
res["B2"] = dict(n=N_NULL, r2_mean=float(np.mean(null_r2)), r2_median=float(np.median(null_r2)),
                 r2_max=float(np.max(null_r2)), r2_min=float(np.min(null_r2)),
                 rho_mean=float(np.mean(null_rho)), rho_max=float(np.max(null_rho)),
                 r2_all=[float(x) for x in null_r2])

# ============================================================ B3 random-gene expression baseline
print("\n==================== B3  RANDOM-GENE EXPRESSION BASELINE ====================")
b3 = []
for it in range(10):
    gs = RNG.choice(genes_pool, size=len(common_tfs), replace=False)
    Xt = zdf(tcga_log.loc[gs].T)
    Xs = zdf(su2c_log.loc[gs].T)
    _, _, r2r, rhor = modelB(Xt, ytr.loc[Xt.index], Xs, yte.loc[Xs.index])
    b3.append((r2r, rhor))
    print(f"   {len(common_tfs)} random genes, rep {it+1}: external R2={r2r:.3f}  rho={rhor:.3f}")
print(f"[B3] mean external R2 = {np.mean([x[0] for x in b3]):.3f}  | real TF model {r2B:.3f}")
res["B3"] = dict(r2_mean=float(np.mean([x[0] for x in b3])), r2_all=[float(x[0]) for x in b3],
                 rho_mean=float(np.mean([x[1] for x in b3])))

# ============================================================ B4 lineage-only baselines
print("\n==================== B4  LINEAGE-ONLY BASELINES ====================")
ar_t, _ = sig_score(tcga_log, AR20); ne_t, _ = sig_score(tcga_log, NE12)
ar_s, _ = sig_score(su2c_log, AR20); ne_s, _ = sig_score(su2c_log, NE12)
sets = {"AR score": ["AR"], "NE score": ["NE"], "AR + NE": ["AR", "NE"]}
Ft = zdf(pd.DataFrame({"AR": ar_t, "NE": ne_t}))
Fs = zdf(pd.DataFrame({"AR": ar_s, "NE": ne_s}))
for nm, cols in sets.items():
    m = LinearRegression().fit(Ft[cols], ytr.loc[Ft.index])
    p = m.predict(Fs[cols])
    print(f"   {nm:10s}: external R2={r2_score(yte.loc[Fs.index], p):.3f}  "
          f"rho={st.spearmanr(yte.loc[Fs.index], p)[0]:.3f}")
    res[f"B4_{nm}"] = dict(r2=float(r2_score(yte.loc[Fs.index], p)),
                           rho=float(st.spearmanr(yte.loc[Fs.index], p)[0]))
# AR + NE + TF model prediction: does the TF model add anything over lineage?
dfm = pd.DataFrame({"obs": yte.loc[Xte.index], "pred": predB,
                    "AR": ar_s.loc[Xte.index], "NE": ne_s.loc[Xte.index]})
pr = st.pearsonr(dfm["obs"], dfm["pred"])[0]
# partial correlation of prediction with observation, controlling AR and NE
import statsmodels.api as sm
rx = sm.OLS(dfm["pred"], sm.add_constant(dfm[["AR", "NE"]])).fit().resid
ry = sm.OLS(dfm["obs"], sm.add_constant(dfm[["AR", "NE"]])).fit().resid
print(f"[B4] Model B pred vs observed: r={pr:.3f}; partial r controlling AR+NE score = "
      f"{st.pearsonr(rx, ry)[0]:.3f}")
res["B4_partial_r_controlling_lineage"] = float(st.pearsonr(rx, ry)[0])

# ============================================================ B5 clinical utility
print("\n==================== B5  CLINICAL-UTILITY FRAMING ====================")
obs = yte.loc[Xte.index]
q1 = obs.quantile(0.25)
lab = (obs <= q1).astype(int)
auc = roc_auc_score(lab, -predB)
print(f"[B5] AUC for identifying the bottom-FOLH1 quartile from the TF model = {auc:.3f}")
# PPV at the model's own bottom quartile
thr = np.quantile(predB, 0.25)
pick = predB <= thr
print(f"[B5] if you call the model's bottom quartile 'PSMA-low': PPV={lab[pick].mean():.2f}, "
      f"sensitivity={lab[pick].sum()/lab.sum():.2f} (n picked={pick.sum()}, n true low={lab.sum()})")
err = (obs - predB)
print(f"[B5] residual SD in z-units = {err.std():.3f}  (i.e. +/-{1.96*err.std():.2f} SD "
      f"95% prediction interval on a standardized scale)")
res["B5"] = dict(auc_bottom_quartile=float(auc), ppv=float(lab[pick].mean()),
                 sens=float(lab[pick].sum() / lab.sum()), resid_sd_z=float(err.std()))

json.dump(res, open(f"{OUT}/s02_results.json", "w"), indent=1)
print("\nsaved s02 outputs")
