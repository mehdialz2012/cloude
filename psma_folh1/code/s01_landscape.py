"""Step 1 — cohort composition + FOLH1 landscape (manuscript Results section 1).

Replicates: FOLH1-AR / FOLH1-NE correlations in TCGA, SU2C, Beltran; NE-feature
comparison; ARSI exposure; Gleason.
Adds (independent checks):
  (a) TCGA sample composition: are the n=550 "primary tumors" actually all tumors?
  (b) all landscape correlations recomputed on tumour-only TCGA
  (c) Beltran clipping sensitivity
"""
import json
import numpy as np
import pandas as pd
from scipy import stats as st
from common import (DATA, OUT, AR20, NE12, load_tcga, load_su2c, load_beltran, sig_score)

res = {}

# ---------------------------------------------------------------- TCGA
tcga_log, tcga_clin = load_tcga()
print(f"TCGA HiSeqV2 matrix: {tcga_log.shape}  (genes x samples)")
suffix = pd.Series([c.split("-")[-1] for c in tcga_log.columns]).value_counts()
print("TCGA sample-type suffixes:\n", suffix.to_string())
res["tcga_n_total"] = int(tcga_log.shape[1])
res["tcga_suffix_counts"] = {k: int(v) for k, v in suffix.items()}

tum = [c for c in tcga_log.columns if c.endswith("-01")]
nrm = [c for c in tcga_log.columns if c.endswith("-11")]
res["tcga_n_tumor"] = len(tum)
res["tcga_n_normal"] = len(nrm)
print(f"\n==> TCGA tumours (-01): {len(tum)}   normals (-11): {len(nrm)}")

# ---------------------------------------------------------------- SU2C / Beltran
su2c_log, su2c_cs, su2c_cp = load_su2c()
nepc_log, nepc_cs = load_beltran()
print(f"SU2C polyA FPKM: {su2c_log.shape}; Beltran FPKM: {nepc_log.shape}")
res["su2c_n"] = int(su2c_log.shape[1])
res["beltran_n"] = int(nepc_log.shape[1])

scores = {}
for name, mat in [("TCGA_all550", tcga_log), ("TCGA_tumoronly", tcga_log[tum]),
                  ("SU2C", su2c_log), ("Beltran", nepc_log)]:
    ar, n_ar = sig_score(mat, AR20)
    ne, n_ne = sig_score(mat, NE12)
    scores[name] = pd.DataFrame({"FOLH1": mat.loc["FOLH1"], "AR_score": ar, "NE_score": ne})
    print(f"{name}: n={mat.shape[1]}, AR genes={n_ar}, NE genes={n_ne}, "
          f"FOLH1 median={mat.loc['FOLH1'].median():.2f}")

rows = []
for name, df in scores.items():
    r_ar, p_ar = st.spearmanr(df["FOLH1"], df["AR_score"])
    r_ne, p_ne = st.spearmanr(df["FOLH1"], df["NE_score"])
    rows.append(dict(cohort=name, n=len(df), rho_FOLH1_AR=round(r_ar, 3), p_AR=p_ar,
                     rho_FOLH1_NE=round(r_ne, 3), p_NE=p_ne))
land = pd.DataFrame(rows)
print("\n=== FOLH1 landscape correlations ===")
print(land.to_string(index=False))
res["landscape"] = land.to_dict("records")

# Beltran by pathology
nepc_map = nepc_cs.set_index("OTHER_SAMPLE_ID")
common = scores["Beltran"].index.intersection(nepc_map.index)
val = scores["Beltran"].loc[common].copy()
val["path"] = nepc_map.loc[common, "CANCER_TYPE_DETAILED"]
print(f"\nBeltran matched samples: {len(common)}")
print(val.groupby("path")[["FOLH1", "AR_score", "NE_score"]].agg(["count", "mean"]).round(2).to_string())
res["beltran_matched"] = int(len(common))
res["beltran_path_counts"] = {k: int(v) for k, v in val["path"].value_counts().items()}

# SU2C NE features / ARSI
su2c = scores["SU2C"].join(
    su2c_cs[["NEUROENDOCRINE_FEATURES", "ABI_ENZA_EXPOSURE_STATUS", "TISSUE_SITE"]], how="left")
for col in ["NEUROENDOCRINE_FEATURES", "ABI_ENZA_EXPOSURE_STATUS"]:
    sub = su2c[[col, "FOLH1"]].dropna()
    groups = [(k, g["FOLH1"].values) for k, g in sub.groupby(col) if len(g) >= 5]
    H, p = st.kruskal(*[v for _, v in groups])
    print(f"\n{col}: Kruskal p={p:.2e} | " +
          ", ".join(f"{k}(n={len(v)}) med={np.median(v):.2f}" for k, v in groups))
    res[f"su2c_{col}_kruskal_p"] = float(p)
    res[f"su2c_{col}_medians"] = {str(k): float(np.median(v)) for k, v in groups}

# Gleason in TCGA tumours
gl = [c for c in tcga_clin.columns if "gleason" in c.lower()]
print("\nTCGA Gleason columns:", gl)
f_t = scores["TCGA_tumoronly"]["FOLH1"]
f_t.index = [c[:12] for c in f_t.index]
for gc in gl:
    sub = pd.DataFrame({"g": tcga_clin[gc]}).join(f_t.rename("FOLH1")).dropna()
    if sub["g"].nunique() < 2 or sub["g"].nunique() > 12:
        continue
    groups = [gg["FOLH1"].values for _, gg in sub.groupby("g")]
    H, p = st.kruskal(*groups)
    meds = {str(k): round(float(np.median(gg["FOLH1"])), 2) for k, gg in sub.groupby("g")}
    print(f"  {gc}: Kruskal p={p:.3e}  medians={meds}")
    res[f"tcga_gleason_{gc}_p"] = float(p)
    res[f"tcga_gleason_{gc}_medians"] = meds

# Beltran clipping sensitivity (negative-floor artifact)
raw = pd.read_csv(f"{DATA}/nepc_data_mrna_seq_fpkm.txt", sep="\t", index_col=0, comment="#")
raw = raw[~raw.index.duplicated(keep="first")].drop(columns=["Entrez_Gene_Id"], errors="ignore")
raw = raw.select_dtypes(include=[np.number])
v = raw.values.flatten()
print(f"\nBeltran raw matrix: frac<0 = {np.nanmean(v < 0):.4f}, min={np.nanmin(v):.2f}")
res["beltran_frac_negative"] = float(np.nanmean(v < 0))

for k, df in scores.items():
    df.to_csv(f"{OUT}/s01_scores_{k}.csv")
val.to_csv(f"{OUT}/s01_beltran_annotated.csv")
su2c.to_csv(f"{OUT}/s01_su2c_annotated.csv")
json.dump(res, open(f"{OUT}/s01_results.json", "w"), indent=1)
print("\nsaved s01 outputs")
