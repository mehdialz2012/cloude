"""Step 4 — target-switching map (manuscript Results section 4).

PART A (replication)
  A1  surfaceome union (SURFY + CSPA high-confidence + UniProt cell-membrane w/ TM or GPI)
  A2  surface genes up in PSMA-low (FDR<0.05, rbc>0.3)
  A3  known-target panel on raw and liver-adjusted matrices with BH-FDR
  A4  HPA vital-organ safety filter and the composite ranking
PART B (independent validation)
  B1  cell-type origin of each candidate (HPA single-cell consensus)
  B2  tumour-to-normal contrast — the quantity that actually matters for a ligand
  B3  absolute antigen abundance in PSMA-low tumours vs FOLH1 in PSMA-high
  B4  "down-regulated" vs "absent": what is retained in PSMA-low disease
"""
import json
import numpy as np
import pandas as pd
import pyreadr
from scipy import stats as st
from statsmodels.stats.multitest import multipletests

from common import DATA, OUT, LIVER10, load_su2c, zrow

res = {}
su2c_log, su2c_cs, _ = load_su2c()
f = su2c_log.loc["FOLH1"]
low, high = f[f <= f.quantile(.25)].index, f[f >= f.quantile(.75)].index
liver_score = zrow(su2c_log.loc[[g for g in LIVER10 if g in su2c_log.index]]).mean(axis=0)
de_adj = pd.read_csv(f"{OUT}/s03_de_adj.csv", index_col=0)

# ---------------------------------------------------------------- A1 surfaceome
surfy = pyreadr.read_r(f"{DATA}/surfy_gene.rda")["surfy_gene"]
cspa = pyreadr.read_r(f"{DATA}/cspa_gene.rda")["cspa_gene"]
sgenes = set(surfy["GENE"].dropna()) | set(surfy.get("UniProt gene", pd.Series(dtype=str)).dropna())
cspa1 = cspa[cspa["CSPA category"] == "1 - high confidence"] if "CSPA category" in cspa else cspa
cgenes = set(cspa1["GENE"].dropna()) | set(cspa1.get("ENTREZ gene symbol", pd.Series(dtype=str)).dropna())
up = pd.read_csv(f"{DATA}/uniprot_cellmembrane.tsv", sep="\t")


def pm_topo(r):
    tm = pd.notna(r.get("Transmembrane")) and "Helical" in str(r.get("Transmembrane"))
    gpi = pd.notna(r.get("Lipidation")) and "GPI-anchor" in str(r.get("Lipidation"))
    return tm or gpi


ugenes = set(g.split()[0] for g in up.loc[up.apply(pm_topo, axis=1), "Gene Names"].dropna())
surfaceome = sgenes | cgenes | ugenes
print(f"[A1] SURFY {len(sgenes)} | CSPA-hi {len(cgenes)} | UniProt-PM+TM/GPI {len(ugenes)} "
      f"| union {len(surfaceome)}   [manuscript 3,701]")
res["A1_surfaceome_n"] = len(surfaceome)

# ---------------------------------------------------------------- A2
surf_de = de_adj[de_adj.index.isin(surfaceome)]
surf_up = surf_de[(surf_de.padj < .05) & (surf_de.rbc > .3)].sort_values("rbc", ascending=False)
surf_dn = surf_de[(surf_de.padj < .05) & (surf_de.rbc < -.3)]
print(f"[A2] surfaceome genes tested {len(surf_de)}; UP in PSMA-low (FDR<0.05, rbc>0.3) "
      f"{len(surf_up)}   [manuscript 78]; DOWN {len(surf_dn)}")
print(surf_up.head(20).round(3).to_string())
res["A2"] = dict(n_tested=int(len(surf_de)), n_up=int(len(surf_up)), n_down=int(len(surf_dn)),
                 top20={k: float(v) for k, v in surf_up["rbc"].head(20).items()})
surf_up.to_csv(f"{OUT}/s04_surface_up.csv")

# ---------------------------------------------------------------- A3 known panel
panel = ['FOLH1', 'STEAP1', 'STEAP2', 'CD276', 'TACSTD2', 'DLL3', 'SSTR2', 'CEACAM5', 'PSCA',
         'GRPR', 'NCAM1', 'SEZ6', 'ENPP3', 'CD46', 'SLC44A4', 'ANO7', 'OR51E2', 'TRPM8', 'KLK3',
         'KLK2', 'ACPP', 'CLDN3', 'CLDN4', 'EPCAM', 'MUC1', 'ERBB2', 'EGFR', 'MET', 'AXL', 'INSM1']
lv = su2c_log.loc[[g for g in LIVER10 if g in su2c_log.index]].mean().values  # authors' unscaled mean
rows = []
for g in panel:
    if g not in su2c_log.index:
        continue
    med = float(su2c_log.loc[g].median())
    a, b = su2c_log.loc[g, low].values, su2c_log.loc[g, high].values
    U, p_raw = st.mannwhitneyu(a, b)
    rbc_raw = 2 * U / (len(a) * len(b)) - 1
    v = su2c_log.loc[g].values
    resid = pd.Series(v - np.polyval(np.polyfit(lv, v, 1), lv), index=su2c_log.columns)
    U2, p_adj = st.mannwhitneyu(resid[low], resid[high])
    rows.append((g, med, med >= 1, rbc_raw, p_raw, 2 * U2 / (len(a) * len(b)) - 1, p_adj))
kp = pd.DataFrame(rows, columns=["gene", "median_log2FPKM", "passes_filter", "rbc_raw", "p_raw",
                                 "rbc_liver_adj", "p_liver_adj"]).set_index("gene")
kp["padj_raw"] = multipletests(kp["p_raw"], method="fdr_bh")[1]
kp["padj_liver_adj"] = multipletests(kp["p_liver_adj"], method="fdr_bh")[1]
print("\n[A3] known-target panel:")
print(kp.round(4).to_string())
res["A3"] = kp.round(5).to_dict("index")
kp.to_csv(f"{OUT}/s04_known_panel.csv")

# ---------------------------------------------------------------- A4 HPA safety
hpa = pd.read_csv(f"{DATA}/rna_tissue_consensus.tsv", sep="\t")
hpa.columns = ["gene_id", "gene", "tissue", "ntpm"]
hpa["tl"] = hpa["tissue"].str.lower()
VITAL = ['heart muscle', 'liver', 'kidney', 'lung', 'bone marrow', 'salivary gland', 'pancreas',
         'small intestine', 'colon', 'cerebral cortex', 'cerebellum', 'basal ganglia',
         'hippocampal formation', 'midbrain', 'spinal cord', 'hypothalamus', 'amygdala',
         'retina', 'spleen', 'stomach 1', 'stomach']
vital = hpa[hpa["tl"].isin(VITAL)].groupby("gene")["ntpm"].max()
allnorm = hpa.groupby("gene")["ntpm"].max()
sal = hpa[hpa["tl"] == "salivary gland"].groupby("gene")["ntpm"].max()
kid = hpa[hpa["tl"] == "kidney"].groupby("gene")["ntpm"].max()
pro = hpa[hpa["tl"] == "prostate"].groupby("gene")["ntpm"].max()

CAND = list(surf_up.index[:40]) + ['MUC1', 'NCAM1', 'CEACAM5', 'SSTR2', 'GPC2', 'PLAUR',
                                   'SLC16A3', 'SLC7A5', 'FOLH1', 'STEAP1', 'STEAP2', 'TACSTD2',
                                   'CD276', 'DLL3', 'KLK3']
CAND = list(dict.fromkeys(CAND))
cand = pd.DataFrame(index=CAND)
cand["rbc"] = de_adj["rbc"].reindex(CAND)
cand["padj"] = de_adj["padj"].reindex(CAND)
cand["max_vital_nTPM"] = vital.reindex(CAND)
cand["max_any_normal_nTPM"] = allnorm.reindex(CAND)
cand["salivary"] = sal.reindex(CAND)
cand["kidney"] = kid.reindex(CAND)
cand["prostate_nTPM"] = pro.reindex(CAND)
cand["tumor_low_median"] = su2c_log.loc[[g for g in CAND if g in su2c_log.index], low]\
    .median(axis=1).reindex(CAND)
cand["tumor_high_median"] = su2c_log.loc[[g for g in CAND if g in su2c_log.index], high]\
    .median(axis=1).reindex(CAND)
c2 = cand.dropna(subset=["rbc"]).copy()
c2["safety"] = 1 / (1 + np.log1p(c2["max_vital_nTPM"]))
c2["composite"] = c2["rbc"] * c2["safety"] * np.log1p(c2["tumor_low_median"])
top = c2.sort_values("composite", ascending=False)
print("\n[A4] composite-ranked candidates (top 15):")
print(top.head(15)[["rbc", "padj", "max_vital_nTPM", "salivary", "kidney",
                    "tumor_low_median", "composite"]].round(3).to_string())
res["A4_top15"] = top.head(15).round(4).to_dict("index")
top.to_csv(f"{OUT}/s04_target_map.csv")

# ================================================================ B1 cell-type origin
print("\n==================== B1  CELL-TYPE ORIGIN (HPA single cell) ====================")
sc = pd.read_csv(f"{DATA}/rna_single_cell_type.tsv", sep="\t")
sc.columns = ["gene_id", "gene", "cell_type", "ncpm"]
IMMUNE_KEY = ['macrophage', 'monocyte', 'granulocyte', 'neutrophil', 't-cell', 'b-cell', 'nk-',
              'dendritic', 'kupffer', 'microglial', 'langerhans', 'plasma cell', 'erythroid',
              'hofbauer', 'mast']
STROMA_KEY = ['fibroblast', 'endothelial', 'smooth muscle', 'peritubular', 'mesenchy',
              'pericyte', 'adipocyte']
FOCUS = ['FOLH1', 'PLAUR', 'SLC16A3', 'SLC7A5', 'GPC2', 'ATP2B4', 'SLC12A4', 'ADAM19', 'MRC2',
         'STEAP1', 'STEAP2', 'TACSTD2', 'CD276', 'SSTR2', 'NCAM1', 'DLL3', 'MUC1']
rows = []
for g in FOCUS:
    s = sc[sc["gene"] == g]
    if not len(s):
        continue
    s = s.sort_values("ncpm", ascending=False)
    tot = s["ncpm"].sum()
    imm = s[s["cell_type"].str.lower().str.contains("|".join(IMMUNE_KEY))]["ncpm"].sum()
    stm = s[s["cell_type"].str.lower().str.contains("|".join(STROMA_KEY))]["ncpm"].sum()
    rows.append(dict(gene=g, top_cell_types="; ".join(s["cell_type"].head(3)),
                     frac_immune=imm / tot if tot else np.nan,
                     frac_stromal=stm / tot if tot else np.nan))
ct = pd.DataFrame(rows).set_index("gene")
pd.set_option("display.max_colwidth", 70)
print(ct.round(3).to_string())
res["B1_celltype"] = ct.to_dict("index")
ct.to_csv(f"{OUT}/s04_celltype_origin.csv")

# ================================================================ B2 tumour:normal contrast
print("\n==================== B2  TUMOUR-TO-NORMAL CONTRAST ====================")
print("    For a radioligand/ADC target the relevant quantity is tumour signal minus the")
print("    worst normal organ, not a rank-biserial correlation between tumour groups.")
b2 = cand.loc[[g for g in FOCUS if g in cand.index]].copy()
b2["log2_max_vital"] = np.log2(b2["max_vital_nTPM"] + 1)
b2["contrast_low_vs_vital"] = b2["tumor_low_median"] - b2["log2_max_vital"]
b2["contrast_high_vs_vital"] = b2["tumor_high_median"] - b2["log2_max_vital"]
print(b2[["rbc", "tumor_low_median", "tumor_high_median", "max_vital_nTPM",
          "log2_max_vital", "contrast_low_vs_vital", "contrast_high_vs_vital"]].round(2).to_string())
res["B2"] = b2.round(3).to_dict("index")

# ================================================================ B3/B4 retention
print("\n==================== B3/B4  WHAT IS ACTUALLY RETAINED IN PSMA-LOW? ====================")
ref = float(su2c_log.loc["FOLH1", high].median())
print(f"    reference: FOLH1 median in PSMA-high tumours = {ref:.2f} log2(FPKM+1)")
b4 = su2c_log.loc[[g for g in FOCUS if g in su2c_log.index]]
tab = pd.DataFrame({
    "median_PSMA_low": b4[low].median(axis=1),
    "median_PSMA_high": b4[high].median(axis=1),
    "pct_of_FOLH1_high_level": 100 * (2 ** b4[low].median(axis=1) - 1) / (2 ** ref - 1),
    "frac_low_samples_above_2log2": (b4[low] >= 2).mean(axis=1),
    "rbc_liver_adj": de_adj["rbc"].reindex(b4.index)})
print(tab.sort_values("median_PSMA_low", ascending=False).round(3).to_string())
res["B4"] = tab.round(4).to_dict("index")
tab.to_csv(f"{OUT}/s04_retention.csv")

json.dump(res, open(f"{OUT}/s04_results.json", "w"), indent=1)
print("\nsaved s04 outputs")
