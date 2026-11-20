"""B4 (lineage-only baselines) and B5 (clinical utility) — light-weight rerun."""
import json
import numpy as np, pandas as pd
from scipy import stats as st
import statsmodels.api as sm
from sklearn.linear_model import ElasticNetCV, LinearRegression
from sklearn.metrics import r2_score, roc_auc_score
from sklearn.preprocessing import StandardScaler
import decoupler as dc
from common import DATA, OUT, AR20, NE12, load_tcga, load_su2c, sig_score

res = {}
tcga_log, _ = load_tcga(); su2c_log, _, _ = load_su2c()
net = pd.read_csv(f"{DATA}/collectri.csv")
y_t, y_s = tcga_log.loc["FOLH1"], su2c_log.loc["FOLH1"]
yt = (y_t - y_t.mean())/y_t.std(); ys = (y_s - y_s.mean())/y_s.std()
zdf = lambda d: pd.DataFrame(StandardScaler().fit_transform(d), index=d.index, columns=d.columns)
ar_t,_ = sig_score(tcga_log, AR20); ne_t,_ = sig_score(tcga_log, NE12)
ar_s,_ = sig_score(su2c_log, AR20); ne_s,_ = sig_score(su2c_log, NE12)
Ft = zdf(pd.DataFrame({"AR": ar_t, "NE": ne_t})); Fs = zdf(pd.DataFrame({"AR": ar_s, "NE": ne_s}))
print("==== B4 lineage-only baselines (external TCGA -> SU2C, z-space) ====")
for nm, cols in [("AR score only", ["AR"]), ("NE score only", ["NE"]), ("AR + NE", ["AR","NE"])]:
    m = LinearRegression().fit(Ft[cols], yt.loc[Ft.index]); p = m.predict(Fs[cols])
    r2 = r2_score(ys.loc[Fs.index], p); rho = st.spearmanr(ys.loc[Fs.index], p)[0]
    print(f"   {nm:14s}: R2={r2:+.3f}  rho={rho:+.3f}")
    res[f"B4_{nm}"] = dict(r2=float(r2), rho=float(rho))
tf_t = dc.mt.ulm(data=tcga_log.T, net=net, verbose=False)[0]
tf_s = dc.mt.ulm(data=su2c_log.T, net=net, verbose=False)[0]
ct = tf_t.columns.intersection(tf_s.columns)
Xt, Xs = zdf(tf_t[ct]), zdf(tf_s[ct])
m = ElasticNetCV(l1_ratio=[0.5,0.7,0.9,1.0], alphas=50, cv=5, random_state=1,
                 max_iter=20000, n_jobs=1).fit(Xt, yt.loc[Xt.index])
pred = m.predict(Xs)
obs = ys.loc[Xs.index]
d = pd.DataFrame({"obs": obs, "pred": pred, "AR": ar_s.loc[Xs.index], "NE": ne_s.loc[Xs.index]})
rx = sm.OLS(d["pred"], sm.add_constant(d[["AR","NE"]])).fit().resid
ry = sm.OLS(d["obs"], sm.add_constant(d[["AR","NE"]])).fit().resid
print(f"   Model B pred vs obs r={st.pearsonr(d.obs,d.pred)[0]:.3f}; "
      f"partial r controlling AR+NE = {st.pearsonr(rx,ry)[0]:.3f}")
res["B4_partial_r"] = float(st.pearsonr(rx, ry)[0])
res["B4_r"] = float(st.pearsonr(d.obs, d.pred)[0])
print("\n==== B5 clinical utility ====")
lab = (obs <= obs.quantile(.25)).astype(int)
auc = roc_auc_score(lab, -pred); thr = np.quantile(pred, .25); pick = pred <= thr
err = obs - pred
print(f"   AUC for the bottom-FOLH1 quartile = {auc:.3f}")
print(f"   model's own bottom quartile: PPV={lab[pick].mean():.2f}, sens={lab[pick].sum()/lab.sum():.2f}")
print(f"   residual SD = {err.std():.3f} SD units; 95% PI = +/-{1.96*err.std():.2f} SD "
      f"= {2**(2*1.96*err.std()*y_s.std()):.0f}-fold span in FPKM")
res["B5"] = dict(auc=float(auc), ppv=float(lab[pick].mean()),
                 sens=float(lab[pick].sum()/lab.sum()), resid_sd=float(err.std()),
                 fold_span=float(2**(2*1.96*err.std()*y_s.std())))
json.dump(res, open(f"{OUT}/s02c_results.json","w"), indent=1)
print("\nsaved s02c")
