"""Phase 2, step 3 — surfaceome: what is retained, lost or enriched in PSMA-low disease.

Statistical enrichment and therapeutic suitability are kept apart. Every surface gene is
scored on six pre-specified axes and then tiered; nothing is ranked by effect size alone.

  A  direction and cross-cohort reproducibility (SU2C M3, GSE126078 patient-level, Beltran adeno)
  B  antigen abundance in PSMA-low tumours (median expression and the fraction of tumours above
     a positivity threshold), in two cohorts
  C  tumour-to-normal contrast against the vital organs that limit radioligand dosing (HPA)
  D  tumour-cell specificity (HPA single-cell: fraction of expression in immune/stromal cells)
  E  protein-level evidence (added in step 3b from HPA IHC)
  F  availability of a clinical modality (curated)
"""
import json
import numpy as np
import pandas as pd
import pyreadr
from scipy import stats as st

from common import DATA, OUT, load_su2c, load_gse126078
from p2_common import fdr

res = {}

# ---------------------------------------------------------------- surfaceome
surfy = pyreadr.read_r(f"{DATA}/surfy_gene.rda")["surfy_gene"]
cspa = pyreadr.read_r(f"{DATA}/cspa_gene.rda")["cspa_gene"]
sg = set(surfy["GENE"].dropna()) | set(surfy.get("UniProt gene", pd.Series(dtype=str)).dropna())
c1 = cspa[cspa["CSPA category"] == "1 - high confidence"] if "CSPA category" in cspa else cspa
cg = set(c1["GENE"].dropna()) | set(c1.get("ENTREZ gene symbol", pd.Series(dtype=str)).dropna())
up = pd.read_csv(f"{DATA}/uniprot_cellmembrane.tsv", sep="\t")


def pm(r):
    tm = pd.notna(r.get("Transmembrane")) and "Helical" in str(r.get("Transmembrane"))
    gpi = pd.notna(r.get("Lipidation")) and "GPI-anchor" in str(r.get("Lipidation"))
    return tm or gpi


ug = set(g.split()[0] for g in up.loc[up.apply(pm, axis=1), "Gene Names"].dropna())
SURF = sg | cg | ug
print(f"surfaceome: {len(SURF)} genes "
      f"(SURFY {len(sg)}, CSPA high-confidence {len(cg)}, UniProt curated cell membrane "
      f"with TM/GPI {len(ug)})")
SRC = pd.DataFrame(index=sorted(SURF))
SRC["surfy"] = SRC.index.isin(sg)
SRC["cspa_hi"] = SRC.index.isin(cg)
SRC["uniprot_cellmembrane"] = SRC.index.isin(ug)
SRC.to_csv(f"{OUT}/p2_surfaceome_sources.csv")

# ---------------------------------------------------------------- A reproducibility
disc = pd.read_csv(f"{OUT}/p2_discovery_models.csv", index_col=0)
gse_m = pd.read_csv(f"{OUT}/p2_gse_models.csv", index_col=0)
bel_a = pd.read_csv(f"{OUT}/p2_beltran_adeno.csv", index_col=0)
idx = [g for g in disc.index if g in SURF]
T = pd.DataFrame(index=idx)
for c in ["surfy", "cspa_hi", "uniprot_cellmembrane"]:
    T[c] = SRC[c].reindex(idx)
T["membrane_confident"] = T["uniprot_cellmembrane"] | (T["cspa_hi"] & T["surfy"])
T["su2c_beta"] = disc.loc[idx, "beta_M3"]
T["su2c_fdr"] = disc.loc[idx, "fdr_M3"]
T["gse_beta"] = gse_m["beta_patient_M3"].reindex(idx)
T["gse_p"] = gse_m["p_patient_M3"].reindex(idx)
T["beltran_beta"] = bel_a["beta"].reindex(idx)
T["beltran_p"] = bel_a["p"].reindex(idx)
# beta is the coefficient on FOLH1 (z): negative = higher when FOLH1 is low
T["enriched_su2c"] = (T["su2c_fdr"] < .05) & (T["su2c_beta"] < 0)
T["lost_su2c"] = (T["su2c_fdr"] < .05) & (T["su2c_beta"] > 0)
T["rep_gse"] = (np.sign(T["su2c_beta"]) == np.sign(T["gse_beta"])) & (T["gse_p"] < .05)
T["rep_beltran"] = (np.sign(T["su2c_beta"]) == np.sign(T["beltran_beta"])) & (T["beltran_p"] < .05)
T["n_cohorts_reproduced"] = T[["rep_gse", "rep_beltran"]].sum(axis=1)

# ---------------------------------------------------------------- B abundance
su2c, cs, cp = load_su2c()
f = su2c.loc["FOLH1"]
low = f[f <= f.quantile(.25)].index
high = f[f >= f.quantile(.75)].index
POS = 2.0   # log2(FPKM+1) >= 2 ~ 3 FPKM: a pragmatic antigen-positivity threshold
T["su2c_low_median"] = su2c.loc[idx, low].median(axis=1)
T["su2c_high_median"] = su2c.loc[idx, high].median(axis=1)
T["su2c_low_pos_frac"] = (su2c.loc[idx, low] >= POS).mean(axis=1)

gse = load_gse126078()
info = pd.read_csv(f"{OUT}/p2_gse_info.csv", index_col=0)
gse = gse[info.index]
fg = gse.loc["FOLH1"]
glow = fg[fg <= fg.quantile(.25)].index
gidx = [g for g in idx if g in gse.index]
T["gse_low_median"] = gse.loc[gidx, glow].median(axis=1).reindex(idx)
T["gse_low_pos_frac"] = (gse.loc[gidx, glow] >= POS).mean(axis=1).reindex(idx)

# ---------------------------------------------------------------- C tumour-to-normal
hpa = pd.read_csv(f"{DATA}/rna_tissue_consensus.tsv", sep="\t")
hpa.columns = ["gid", "gene", "tissue", "ntpm"]
hpa["tl"] = hpa["tissue"].str.lower()
VITAL = ['heart muscle', 'liver', 'kidney', 'lung', 'bone marrow', 'salivary gland', 'pancreas',
         'small intestine', 'colon', 'cerebral cortex', 'cerebellum', 'basal ganglia',
         'hippocampal formation', 'midbrain', 'spinal cord', 'hypothalamus', 'amygdala',
         'retina', 'spleen', 'stomach 1', 'stomach', 'esophagus', 'duodenum', 'rectum',
         'thyroid gland', 'parathyroid gland', 'adrenal gland', 'urinary bladder']
vital = hpa[hpa["tl"].isin(VITAL)].groupby("gene")["ntpm"].max()
worst = hpa[hpa["tl"].isin(VITAL)].sort_values("ntpm").groupby("gene").last()["tissue"]
sal = hpa[hpa["tl"] == "salivary gland"].groupby("gene")["ntpm"].max()
kid = hpa[hpa["tl"] == "kidney"].groupby("gene")["ntpm"].max()
T["max_vital_nTPM"] = vital.reindex(idx)
T["worst_organ"] = worst.reindex(idx)
T["salivary_nTPM"] = sal.reindex(idx)
T["kidney_nTPM"] = kid.reindex(idx)
T["contrast_low"] = T["su2c_low_median"] - np.log2(T["max_vital_nTPM"] + 1)

# ---------------------------------------------------------------- D tumour-cell specificity
sc = pd.read_csv(f"{DATA}/rna_single_cell_type.tsv", sep="\t")
sc.columns = ["gid", "gene", "cell_type", "ncpm"]
IMM = ['macrophage', 'monocyte', 'granulocyte', 'neutrophil', 't-cell', 'b-cell', 'nk-',
       'dendritic', 'kupffer', 'microglial', 'langerhans', 'plasma cell', 'erythroid',
       'hofbauer', 'mast']
STR = ['fibroblast', 'endothelial', 'smooth muscle', 'peritubular', 'mesenchy', 'pericyte',
       'adipocyte']
sc["ct"] = sc["cell_type"].str.lower()
sc["is_imm"] = sc["ct"].str.contains("|".join(IMM))
sc["is_str"] = sc["ct"].str.contains("|".join(STR))
tot = sc.groupby("gene")["ncpm"].sum()
imm = sc[sc["is_imm"]].groupby("gene")["ncpm"].sum()
stm = sc[sc["is_str"]].groupby("gene")["ncpm"].sum()
T["frac_immune_sc"] = (imm.reindex(idx) / tot.reindex(idx))
T["frac_stromal_sc"] = (stm.reindex(idx) / tot.reindex(idx))
T["frac_nontumour_sc"] = T["frac_immune_sc"].fillna(0) + T["frac_stromal_sc"].fillna(0)

# ---------------------------------------------------------------- classification
def klass(r):
    if r["enriched_su2c"]:
        return "enriched"
    if r["lost_su2c"]:
        return "lost"
    return "unchanged"


T["state"] = T.apply(klass, axis=1)
T["retained"] = (T["su2c_low_pos_frac"] >= 0.5)
print("\n[A] surface genes by state in PSMA-low tumours (SU2C, fully adjusted):")
print(T["state"].value_counts().to_string())
print(f"[B] of those, expressed in >=50% of PSMA-low tumours: "
      f"{T.groupby('state')['retained'].sum().to_dict()}")
res["state_counts"] = {k: int(v) for k, v in T["state"].value_counts().items()}
res["retained_by_state"] = {k: int(v) for k, v in T.groupby("state")["retained"].sum().items()}

# ---------------------------------------------------------------- suitability screen
T["pass_abundance"] = T["su2c_low_pos_frac"] >= 0.5
T["pass_contrast"] = T["contrast_low"] > 0
T["pass_specificity"] = T["frac_nontumour_sc"] < 0.25
T["pass_reproducible"] = (T["n_cohorts_reproduced"] >= 1) | (T["state"] == "unchanged")
T["n_pass"] = T[["pass_abundance", "pass_contrast", "pass_specificity"]].sum(axis=1)
print(f"\n[C] surface genes passing all three suitability filters "
      f"(abundance, tumour-to-normal contrast, tumour-cell specificity): "
      f"{int((T['n_pass']==3).sum())}")
print(T[T["n_pass"] == 3].sort_values("su2c_low_median", ascending=False)
      .head(25)[["state", "su2c_low_median", "su2c_low_pos_frac", "contrast_low",
                 "frac_nontumour_sc", "max_vital_nTPM", "worst_organ"]].round(3).to_string())
res["n_pass_all3"] = int((T["n_pass"] == 3).sum())

PANEL = ['FOLH1', 'STEAP1', 'STEAP2', 'CD276', 'TACSTD2', 'DLL3', 'SSTR2', 'CEACAM5', 'PSCA',
         'GRPR', 'NCAM1', 'SEZ6', 'ENPP3', 'CD46', 'SLC44A4', 'ANO7', 'TRPM8', 'CLDN3',
         'CLDN4', 'EPCAM', 'MUC1', 'ERBB2', 'CCKBR', 'KLK3',
         'SLC16A3', 'SLC7A5', 'PLAUR', 'GPC2', 'ATP2B4', 'SLC12A4', 'ADAM19', 'MRC2',
         'PAQR8', 'SLC6A6', 'DCHS1', 'RELT', 'MXRA5', 'ACTN1', 'SLC20A1', 'UNC5B', 'VCAN']
panel = T.reindex([g for g in PANEL if g in T.index])
print("\n[D] the established and newly nominated targets, side by side:")
print(panel[["state", "su2c_beta", "su2c_fdr", "n_cohorts_reproduced", "su2c_low_median",
             "su2c_low_pos_frac", "contrast_low", "frac_nontumour_sc", "n_pass"]]
      .round(3).to_string())
T.to_csv(f"{OUT}/p2_surfaceome_full.csv")
panel.to_csv(f"{OUT}/p2_surfaceome_panel.csv")

# shortlist for protein-level follow-up
short = sorted(set(panel.index) |
               set(T[T["n_pass"] == 3].sort_values("su2c_low_median", ascending=False)
                   .head(20).index))
pd.Series(short, name="gene").to_csv(f"{OUT}/p2_protein_shortlist.csv", index=False)
print(f"\n[E] shortlist for protein-level evidence: {len(short)} genes")
res["shortlist"] = short
json.dump(res, open(f"{OUT}/p2_03_results.json", "w"), indent=1)
print("\nsaved p2 step 3")
