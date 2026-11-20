"""Phase 2, step 2b — Beltran used for its intended purpose: neuroendocrine focus.

In a 49-sample cohort that is 31% NEPC, FOLH1 and the NE score are collinear (r=-0.46)
and conditioning on NE acts as a suppressor: median |t| for the FOLH1 term rises from
1.38 (composition-adjusted) to 2.68 (NE-adjusted). The NE-adjusted fit is therefore
reported as unstable, and Beltran is used instead to ask the NE-specific questions.
"""
import json
import numpy as np, pandas as pd
from scipy import stats as st
from common import OUT, load_beltran, sig_score, AR20, NE12
from p2_common import composition_scores, design, ols_beta, fdr

res = {}
bel, bcs = load_beltran()
bm = bcs.set_index("OTHER_SAMPLE_ID")
idx = bel.columns.intersection(bm.index)
bel = bel[idx]
path = bm.loc[idx, "CANCER_TYPE_DETAILED"]
is_nepc = (path == "Prostate Neuroendocrine Carcinoma")
f = bel.loc["FOLH1"]
print(f"Beltran: {len(idx)} samples, {int(is_nepc.sum())} NEPC, {int((~is_nepc).sum())} CRPC-adeno")
print(f"FOLH1 median: NEPC {f[is_nepc].median():.2f} vs adeno {f[~is_nepc].median():.2f} "
      f"(MWU p={st.mannwhitneyu(f[is_nepc], f[~is_nepc]).pvalue:.1e})")
res["folh1_nepc"] = float(f[is_nepc].median()); res["folh1_adeno"] = float(f[~is_nepc].median())
res["folh1_nepc_p"] = float(st.mannwhitneyu(f[is_nepc], f[~is_nepc]).pvalue)
print(f"NEPC samples in the cohort's bottom FOLH1 quartile: "
      f"{int((f[is_nepc] <= f.quantile(.25)).sum())}/{int(is_nepc.sum())}; "
      f"adeno: {int((f[~is_nepc] <= f.quantile(.25)).sum())}/{int((~is_nepc).sum())}")

expr = bel.index[bel.median(axis=1) >= 1.0]
disc = pd.read_csv(f"{OUT}/p2_discovery_models.csv", index_col=0)
cov = ["immune", "stromal", "hepatocyte", "osteoblast"]

# (a) composition-adjusted whole cohort
comp = composition_scores(bel); comp["FOLH1_z"] = (f - f.mean()) / f.std()
b2, _, _, p2 = ols_beta(bel.loc[expr].values, design(comp, ["FOLH1_z"] + cov).values, coef=1)
out = pd.DataFrame({"beta_M2": b2, "p_M2": p2}, index=expr)
out["fdr_M2"] = fdr(out["p_M2"])

# (b) within CRPC-adenocarcinoma only — removes the NE axis by design, not by covariate
ad = bel.loc[:, ~is_nepc]
expr_a = ad.index[ad.median(axis=1) >= 1.0]
ca = composition_scores(ad); fa = ad.loc["FOLH1"]
ca["FOLH1_z"] = (fa - fa.mean()) / fa.std()
ba, _, _, pa = ols_beta(ad.loc[expr_a].values, design(ca, ["FOLH1_z"] + cov).values, coef=1)
adeno = pd.DataFrame({"beta": ba, "p": pa}, index=expr_a)
adeno["fdr"] = fdr(adeno["p"])
print(f"\n[a] whole cohort, composition-adjusted: {int((out['fdr_M2']<.05).sum())} genes FDR<0.05")
print(f"[b] within CRPC-adenocarcinoma only (n={ad.shape[1]}): "
      f"{int((adeno['fdr']<.05).sum())} genes FDR<0.05")
ov = disc.index.intersection(out.index)
ova = disc.index.intersection(adeno.index)
print(f"[c] concordance with SU2C discovery (M3 betas): whole cohort "
      f"rho={st.spearmanr(disc.loc[ov,'beta_M3'], out.loc[ov,'beta_M2'])[0]:.3f}; "
      f"adeno-only rho={st.spearmanr(disc.loc[ova,'beta_M3'], adeno.loc[ova,'beta'])[0]:.3f}")
res["n_sig_M2"] = int((out["fdr_M2"] < .05).sum())
res["n_sig_adeno"] = int((adeno["fdr"] < .05).sum())
res["n_adeno"] = int(ad.shape[1])
res["conc_whole"] = float(st.spearmanr(disc.loc[ov, "beta_M3"], out.loc[ov, "beta_M2"])[0])
res["conc_adeno"] = float(st.spearmanr(disc.loc[ova, "beta_M3"], adeno.loc[ova, "beta"])[0])
out.to_csv(f"{OUT}/p2_beltran_M2.csv"); adeno.to_csv(f"{OUT}/p2_beltran_adeno.csv")
json.dump(res, open(f"{OUT}/p2_02b_results.json", "w"), indent=1)
print("\nsaved p2 step 2b")
