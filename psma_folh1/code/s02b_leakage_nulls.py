"""Step 2B — the validation tests for the TF model (part B of s02, run separately).

B0  AR-only external baseline + ELF5 weight (completes the replication)
B1  TARGET LEAKAGE: FOLH1 is a member of the CollecTRI regulons of ELF5/ELF3/AR/VDR.
    Recompute the whole model with FOLH1 removed from every regulon.
B2  NULL MODEL: size-matched random gene sets used as pseudo-regulons.
B3  PLAIN-TRANSCRIPTOME BASELINE: elastic net on N random genes' expression.
B4  LINEAGE-ONLY BASELINE: AR score, NE score, AR+NE.
B5  CLINICAL UTILITY: can the model identify the PSMA-low quartile in SU2C?
B6  RAW-SPACE TRANSFER (manuscript: R2<0) — single fit with a bounded iteration count.
"""
import json
import numpy as np
import pandas as pd
from scipy import stats as st
from statsmodels.stats.multitest import multipletests
import statsmodels.api as sm
from sklearn.linear_model import ElasticNetCV, LinearRegression
from sklearn.metrics import r2_score, roc_auc_score
from sklearn.preprocessing import StandardScaler
import decoupler as dc

from common import DATA, OUT, AR20, NE12, load_tcga, load_su2c, sig_score

RNG = np.random.default_rng(1)
res = {}
tcga_log, _ = load_tcga()
su2c_log, _, _ = load_su2c()
net = pd.read_csv(f"{DATA}/collectri.csv")
y_t, y_s = tcga_log.loc["FOLH1"], su2c_log.loc["FOLH1"]
yt = (y_t - y_t.mean()) / y_t.std()
ys = (y_s - y_s.mean()) / y_s.std()


def tf_act(mat, network):
    return dc.mt.ulm(data=mat.T, net=network, verbose=False)[0]


def zdf(d):
    return pd.DataFrame(StandardScaler().fit_transform(d), index=d.index, columns=d.columns)


def fit_transfer(Xtr, ytr, Xte, yte, max_iter=20000):
    m = ElasticNetCV(l1_ratio=[0.5, 0.7, 0.9, 1.0], alphas=50, cv=5, random_state=1,
                     max_iter=max_iter, n_jobs=4).fit(Xtr, ytr)
    p = m.predict(Xte)
    return m, p, r2_score(yte, p), st.spearmanr(yte, p)[0]


def assoc(acts, y):
    c = acts.index.intersection(y.index)
    r = pd.DataFrame([(tf, *st.spearmanr(acts.loc[c, tf], y[c])) for tf in acts.columns],
                     columns=["TF", "rho", "p"]).set_index("TF")
    r["padj"] = multipletests(r["p"], method="fdr_bh")[1]
    return r


# ---------------------------------------------------------------- reference model
tf_t, tf_s = tf_act(tcga_log, net), tf_act(su2c_log, net)
ct = tf_t.columns.intersection(tf_s.columns)
Xt, Xs = zdf(tf_t[ct]), zdf(tf_s[ct])
mB, predB, r2B, rhoB = fit_transfer(Xt, yt.loc[Xt.index], Xs, ys.loc[Xs.index])
w = pd.Series(mB.coef_, index=ct)
print(f"[ref] Model B external R2={r2B:.3f} rho={rhoB:.3f}; non-zero TFs={int((w!=0).sum())}")
print(f"[B0] ELF5 weight = {w.get('ELF5', 0):+.3f}   [manuscript +0.655 dominant positive]")
print(f"[B0] most negative weights: {w.nsmallest(4).round(3).to_dict()}   "
      f"[manuscript ZNF24, YBX1, SOX4]")
ar_t, _ = sig_score(tcga_log, AR20)
ar_s, _ = sig_score(su2c_log, AR20)
ne_t, _ = sig_score(tcga_log, NE12)
ne_s, _ = sig_score(su2c_log, NE12)
b = LinearRegression().fit(zdf(ar_t.to_frame("AR")), yt.loc[ar_t.index])
r2_ar = r2_score(ys.loc[ar_s.index], b.predict(zdf(ar_s.to_frame("AR"))))
print(f"[B0] AR-score-only external R2 = {r2_ar:.3f}   [manuscript 0.382]")
res["B0"] = dict(r2_ext=float(r2B), rho=float(rhoB), elf5_weight=float(w.get("ELF5", 0)),
                 n_nonzero=int((w != 0).sum()), r2_ar_only=float(r2_ar),
                 top_pos={k: float(v) for k, v in w.nlargest(6).items()},
                 top_neg={k: float(v) for k, v in w.nsmallest(6).items()})

# ================================================================ B1 leakage
print("\n==================== B1  TARGET-LEAKAGE TEST ====================")
net_nl = net[net["target"] != "FOLH1"]
tf_t2, tf_s2 = tf_act(tcga_log, net_nl), tf_act(su2c_log, net_nl)
ct2 = tf_t2.columns.intersection(tf_s2.columns)
mB2, predB2, r2B2, rhoB2 = fit_transfer(zdf(tf_t2[ct2]), yt.loc[tf_t2.index],
                                        zdf(tf_s2[ct2]), ys.loc[tf_s2.index])
w2 = pd.Series(mB2.coef_, index=ct2)
a1t, a1s = assoc(tf_t, y_t), assoc(tf_s, y_s)
a2t, a2s = assoc(tf_t2, y_t), assoc(tf_s2, y_s)
print(f"[B1] external R2  {r2B:.3f} -> {r2B2:.3f}   (rho {rhoB:.3f} -> {rhoB2:.3f})")
print(f"[B1] ELF5 weight  {w.get('ELF5',0):+.3f} -> {w2.get('ELF5',0):+.3f}")
for tf in ["ELF5", "ELF3", "AR", "VDR"]:
    print(f"[B1] {tf:5s} TF-FOLH1 rho: TCGA {a1t.loc[tf,'rho']:+.3f} -> {a2t.loc[tf,'rho']:+.3f}"
          f" | SU2C {a1s.loc[tf,'rho']:+.3f} -> {a2s.loc[tf,'rho']:+.3f}")
# how much of ELF5's activity is FOLH1 itself?
elf5_targets = net.loc[net["source"] == "ELF5", "target"].tolist()
present = [g for g in elf5_targets if g in su2c_log.index]
print(f"[B1] ELF5 regulon = {len(elf5_targets)} genes ({len(present)} measured in SU2C): {present}")
print(f"[B1] correlation between ELF5 activity and FOLH1 alone (SU2C) = "
      f"{a1s.loc['ELF5','rho']:+.3f}; with FOLH1 removed from the regulon = "
      f"{a2s.loc['ELF5','rho']:+.3f}")
res["B1"] = dict(r2_noleak=float(r2B2), rho_noleak=float(rhoB2),
                 elf5_w_noleak=float(w2.get("ELF5", 0)),
                 elf5_regulon=present,
                 per_tf={tf: dict(tcga=float(a1t.loc[tf, "rho"]), tcga_nl=float(a2t.loc[tf, "rho"]),
                                  su2c=float(a1s.loc[tf, "rho"]), su2c_nl=float(a2s.loc[tf, "rho"]))
                         for tf in ["ELF5", "ELF3", "AR", "VDR"]})

# ================================================================ B2 random-regulon null
print("\n==================== B2  RANDOM-REGULON NULL ====================")
sizes = net_nl.groupby("source").size()
sizes = sizes[sizes >= 5]
pool = np.array([g for g in tcga_log.index.intersection(su2c_log.index) if g != "FOLH1"])
N = 12
nr2, nrho = [], []
for i in range(N):
    rows = [pd.DataFrame({"source": f"R{j}", "target": RNG.choice(pool, int(k), replace=False),
                          "weight": 1.0}) for j, (tf, k) in enumerate(sizes.items())]
    rnet = pd.concat(rows, ignore_index=True)
    at, asu = tf_act(tcga_log, rnet), tf_act(su2c_log, rnet)
    cc = at.columns.intersection(asu.columns)
    _, _, r2n, rhon = fit_transfer(zdf(at[cc]), yt.loc[at.index], zdf(asu[cc]), ys.loc[asu.index])
    nr2.append(r2n); nrho.append(rhon)
    print(f"   null {i+1:2d}: {len(cc)} pseudo-regulons -> external R2={r2n:.3f}  rho={rhon:.3f}")
print(f"[B2] random-regulon null: R2 mean {np.mean(nr2):.3f} (range {min(nr2):.3f}-{max(nr2):.3f})"
      f" | real CollecTRI model {r2B:.3f}")
res["B2"] = dict(n=N, r2_mean=float(np.mean(nr2)), r2_min=float(min(nr2)), r2_max=float(max(nr2)),
                 rho_mean=float(np.mean(nrho)), r2_all=[float(x) for x in nr2])

# ================================================================ B3 random genes
print("\n==================== B3  RANDOM-GENE EXPRESSION BASELINE ====================")
b3 = []
for i in range(8):
    gs = RNG.choice(pool, len(ct), replace=False)
    A, B = zdf(tcga_log.loc[gs].T), zdf(su2c_log.loc[gs].T)
    _, _, r2r, rhor = fit_transfer(A, yt.loc[A.index], B, ys.loc[B.index])
    b3.append((r2r, rhor))
    print(f"   {len(ct)} random genes, rep {i+1}: external R2={r2r:.3f}  rho={rhor:.3f}")
print(f"[B3] mean external R2 = {np.mean([x[0] for x in b3]):.3f} | real TF model {r2B:.3f}")
res["B3"] = dict(r2_mean=float(np.mean([x[0] for x in b3])),
                 r2_all=[float(x[0]) for x in b3],
                 rho_mean=float(np.mean([x[1] for x in b3])))

# ================================================================ B4 lineage baselines
print("\n==================== B4  LINEAGE-ONLY BASELINES ====================")
Ft = zdf(pd.DataFrame({"AR": ar_t, "NE": ne_t}))
Fs = zdf(pd.DataFrame({"AR": ar_s, "NE": ne_s}))
for nm, cols in [("AR score", ["AR"]), ("NE score", ["NE"]), ("AR + NE", ["AR", "NE"])]:
    m = LinearRegression().fit(Ft[cols], yt.loc[Ft.index])
    p = m.predict(Fs[cols])
    print(f"   {nm:9s}: external R2={r2_score(ys.loc[Fs.index], p):+.3f}  "
          f"rho={st.spearmanr(ys.loc[Fs.index], p)[0]:+.3f}")
    res[f"B4_{nm}"] = dict(r2=float(r2_score(ys.loc[Fs.index], p)),
                           rho=float(st.spearmanr(ys.loc[Fs.index], p)[0]))
d = pd.DataFrame({"obs": ys.loc[Xs.index], "pred": predB,
                  "AR": ar_s.loc[Xs.index], "NE": ne_s.loc[Xs.index]})
rx = sm.OLS(d["pred"], sm.add_constant(d[["AR", "NE"]])).fit().resid
ry = sm.OLS(d["obs"], sm.add_constant(d[["AR", "NE"]])).fit().resid
print(f"[B4] predicted vs observed r={st.pearsonr(d.obs,d.pred)[0]:.3f}; "
      f"partial r controlling AR+NE scores = {st.pearsonr(rx,ry)[0]:.3f}")
res["B4_partial_r"] = float(st.pearsonr(rx, ry)[0])

# ================================================================ B5 clinical utility
print("\n==================== B5  CLINICAL-UTILITY FRAMING ====================")
obs = ys.loc[Xs.index]
lab = (obs <= obs.quantile(.25)).astype(int)
auc = roc_auc_score(lab, -predB)
thr = np.quantile(predB, .25)
pick = predB <= thr
err = obs - predB
print(f"[B5] AUC for identifying the bottom-FOLH1 quartile = {auc:.3f}")
print(f"[B5] calling the model's own bottom quartile 'PSMA-low': PPV={lab[pick].mean():.2f}, "
      f"sensitivity={lab[pick].sum()/lab.sum():.2f}")
print(f"[B5] residual SD = {err.std():.3f} SD units -> 95% prediction interval +/-"
      f"{1.96*err.std():.2f} SD of FOLH1 (= {2**(1.96*err.std()*y_s.std()):.1f}-fold in FPKM)")
res["B5"] = dict(auc=float(auc), ppv=float(lab[pick].mean()),
                 sens=float(lab[pick].sum() / lab.sum()), resid_sd=float(err.std()),
                 fold_pi=float(2 ** (1.96 * err.std() * y_s.std())))

# ================================================================ B6 raw-space
print("\n==================== B6  RAW-SPACE TRANSFER ====================")
try:
    _, _, r2raw, _ = fit_transfer(tf_t[ct], y_t.loc[tf_t.index], tf_s[ct], y_s.loc[tf_s.index],
                                  max_iter=1500)
    print(f"[B6] raw-space (unstandardized) external R2 = {r2raw:.2f}   "
          f"[manuscript 'R2<0'; naive transfer -12]")
    res["B6_r2_raw"] = float(r2raw)
except Exception as e:
    print(f"[B6] failed: {type(e).__name__}")

json.dump(res, open(f"{OUT}/s02b_results.json", "w"), indent=1)
print("\nsaved s02b outputs")
