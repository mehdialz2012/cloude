"""Phase 2, step 3b — protein-level evidence and the final suitability tiering.

Adds, for the shortlist:
  - HPA immunohistochemistry in prostate cancer (staining levels across patients)
  - HPA immunohistochemistry in the organs that limit radioligand dosing
  - HPA subcellular localisation (is the protein actually at the plasma membrane?)
  - antibody reliability
and converts the tumour-to-normal contrast into the only scale that is interpretable
across platforms: PSMA's own contrast in PSMA-high disease.
"""
import json
import os
import re
import time

import numpy as np
import pandas as pd
import requests

from common import DATA, OUT

CACHE = f"{DATA}/hpa_xml"
os.makedirs(CACHE, exist_ok=True)
LEVELS = {"not detected": 0, "low": 1, "medium": 2, "high": 3}
DOSE_LIMITING = ["salivary gland", "kidney", "liver", "bone marrow", "lacrimal gland",
                 "duodenum", "small intestine", "colon", "pancreas", "spleen", "lung",
                 "heart muscle", "cerebral cortex"]

short = pd.read_csv(f"{OUT}/p2_protein_shortlist.csv")["gene"].tolist()
T = pd.read_csv(f"{OUT}/p2_surfaceome_full.csv", index_col=0)

# gene -> Ensembl from the full HPA table
hpa_full = pd.read_csv(f"{DATA}/proteinatlas.tsv", sep="\t", usecols=["Gene", "Ensembl",
                                                                     "Subcellular location",
                                                                     "Subcellular main location",
                                                                     "Reliability (IH)",
                                                                     "Protein class"])
g2e = hpa_full.set_index("Gene")["Ensembl"].to_dict()


def fetch(gene):
    ens = g2e.get(gene)
    if not ens:
        return None
    p = f"{CACHE}/{ens}.xml"
    if os.path.exists(p) and os.path.getsize(p) > 1000:
        return open(p, encoding="utf8", errors="ignore").read()
    for attempt in range(4):
        try:
            r = requests.get(f"https://www.proteinatlas.org/{ens}.xml", timeout=180)
            if r.status_code == 200 and len(r.content) > 1000:
                open(p, "wb").write(r.content)
                return r.text
        except Exception:
            pass
        time.sleep(5)
    return None


def parse(xml):
    """Return (normal tissue IHC max level per tissue, prostate-cancer staining counts)."""
    normal, cancer = {}, []
    for block in re.finditer(r'<tissueExpression[^>]*assayType="tissue"[^>]*>(.*?)</tissueExpression>',
                             xml, re.S):
        for d in re.finditer(r"<data>(.*?)</data>", block.group(1), re.S):
            b = d.group(1)
            tis = re.findall(r"<tissue[^>]*>([^<]*)</tissue>", b)
            lv = re.findall(r'<level type="expression"[^>]*>([^<]*)</level>', b)
            for t in tis:
                vals = [LEVELS.get(x.lower(), np.nan) for x in lv]
                vals = [v for v in vals if not np.isnan(v)]
                if vals:
                    normal[t.lower()] = max(normal.get(t.lower(), 0), max(vals))
    for block in re.finditer(r'<tissueExpression[^>]*assayType="cancer"[^>]*>(.*?)</tissueExpression>',
                             xml, re.S):
        for d in re.finditer(r"<data>(.*?)</data>", block.group(1), re.S):
            b = d.group(1)
            tis = re.findall(r"<tissue[^>]*>([^<]*)</tissue>", b)
            if not any("prostate" in t.lower() for t in tis):
                continue
            for pt in re.finditer(r"<patient>(.*?)</patient>", b, re.S):
                lv = re.search(r'<level type="staining"[^>]*>([^<]*)</level>', pt.group(1))
                if lv:
                    cancer.append(lv.group(1).lower())
    return normal, cancer


rows = []
for i, g in enumerate(short):
    xml = fetch(g)
    if xml is None:
        rows.append(dict(gene=g))
        print(f"  {g}: no HPA entry")
        continue
    normal, cancer = parse(xml)
    dl = {t: normal.get(t, np.nan) for t in DOSE_LIMITING}
    dl_max = np.nanmax([v for v in dl.values() if not np.isnan(v)]) if any(
        not np.isnan(v) for v in dl.values()) else np.nan
    worst = max([(v, t) for t, v in dl.items() if not np.isnan(v)], default=(np.nan, None))[1]
    n = len(cancer)
    frac_pos = np.mean([c != "not detected" for c in cancer]) if n else np.nan
    frac_strong = np.mean([c in ("medium", "high") for c in cancer]) if n else np.nan
    hr = hpa_full[hpa_full["Gene"] == g]
    loc = hr["Subcellular main location"].iloc[0] if len(hr) else None
    rel = hr["Reliability (IH)"].iloc[0] if len(hr) else None
    rows.append(dict(gene=g, ihc_prostate_cancer_n=n, ihc_pca_frac_positive=frac_pos,
                     ihc_pca_frac_medium_high=frac_strong,
                     ihc_dose_limiting_max=dl_max, ihc_worst_organ=worst,
                     ihc_salivary=normal.get("salivary gland", np.nan),
                     ihc_kidney=normal.get("kidney", np.nan),
                     ihc_liver=normal.get("liver", np.nan),
                     subcellular_main=loc, ihc_reliability=rel,
                     membrane_annotated=bool(loc and "membrane" in str(loc).lower())))
    print(f"  {g:10s} PCa IHC n={n:3d} pos={frac_pos if n else float('nan'):.2f} "
          f"med/high={frac_strong if n else float('nan'):.2f} | "
          f"dose-limiting max={dl_max} ({worst}) | {loc}")
P = pd.DataFrame(rows).set_index("gene")
P.to_csv(f"{OUT}/p2_protein_evidence.csv")

# ---------------------------------------------------------------- final table
M = T.loc[[g for g in short if g in T.index]].join(P)
bench = float(T.loc["FOLH1", "su2c_high_median"] - np.log2(T.loc["FOLH1", "max_vital_nTPM"] + 1))
M["contrast_vs_PSMA_benchmark"] = M["contrast_low"] - bench
print(f"\n[benchmark] PSMA's own tumour-to-normal contrast in PSMA-high disease = {bench:+.2f} log2")

# protein-level expression evidence; membrane evidence is a separate criterion
M["pass_protein"] = M["ihc_pca_frac_positive"] >= 0.5
# selectivity is only interpretable against the one target that is known to work:
# a candidate is credible if its tumour-to-normal contrast is within 2 log2 units
# (4-fold) of PSMA's own contrast in PSMA-high disease.
M["pass_selectivity"] = M["contrast_vs_PSMA_benchmark"] > -2.0
M["pass_membrane"] = M["membrane_confident"].fillna(False)
CRIT = ["pass_abundance", "pass_selectivity", "pass_specificity", "pass_protein"]
M["n_criteria"] = M[CRIT].sum(axis=1)

MODALITY = {
 "CD276": "B7-H3; ADC (ifinatamab deruxtecan) and CAR-T in trials; 177Lu/131I-omburtamab tested",
 "STEAP1": "ADC (vandortuzumab vedotin), T-cell engager (xaluritamig) in mCRPC trials",
 "STEAP2": "preclinical ADC/CAR-T only",
 "TACSTD2": "TROP2; approved ADC (sacituzumab govitecan) in other tumours",
 "SLC7A5": "LAT1; PET tracers; JPH203 phase 1",
 "SLC16A3": "MCT4; no clinical agent (AZD3965 targets MCT1)",
 "PLAUR": "uPAR; AE105 PET in phase 2; ADC/CAR-T preclinical",
 "GPC2": "ADC and CAR-T preclinical in neuroblastoma",
 "DLL3": "approved T-cell engager (tarlatamab) in SCLC; NEPC trials",
 "SSTR2": "approved RLT (177Lu-DOTATATE)",
 "CEACAM5": "ADC in trials; approved imaging agent",
 "NCAM1": "CD56; ADC tested in SCLC (lorvotuzumab)",
 "MUC1": "multiple ADC/CAR-T programmes",
 "ERBB2": "approved ADCs",
 "EPCAM": "catumaxomab (withdrawn); pan-epithelial",
 "CD46": "ADC (FOR46) in mCRPC phase 1",
 "SLC44A4": "ADC (ASG-5ME) tested in mCRPC",
 "PSCA": "CAR-T in mCRPC phase 1",
 "FOLH1": "the reference target: 177Lu-PSMA-617 (VISION, TheraP, PSMAfore)",
}
M["modality"] = pd.Series(MODALITY).reindex(M.index)


# HPA immunohistochemistry is antibody-dependent (STEAP1 and PSCA, both validated
# clinical targets, score 1/9 and 0/12 positive here), so protein data is carried as
# supporting evidence with its reliability flag rather than as a gate.
CORE = ["pass_abundance", "pass_selectivity", "pass_specificity"]
M["n_core"] = M[CORE].sum(axis=1)


def tier(r):
    if not r["pass_membrane"]:
        return "X - excluded (not confidently plasma membrane)"
    if r["n_core"] == 3:
        return "A - meets all suitability criteria"
    if r["pass_abundance"] and r["pass_selectivity"]:
        return "B - credible, specificity uncertain"
    if r["pass_abundance"] and r["pass_specificity"]:
        return "C - retained but poor tumour-to-normal contrast"
    if r["state"] == "enriched":
        return "D - enriched but not therapeutically suitable"
    return "E - not supported for PSMA-low disease"


M["tier"] = M.apply(tier, axis=1)
cols = ["tier", "state", "n_cohorts_reproduced", "su2c_low_median", "su2c_low_pos_frac",
        "contrast_vs_PSMA_benchmark", "frac_nontumour_sc", "ihc_pca_frac_positive",
        "ihc_pca_frac_medium_high", "ihc_reliability", "n_core"]
M = M.sort_values(["tier", "su2c_low_median"], ascending=[True, False])
print("\n[final] surface-target prioritisation for PSMA-low mCRPC")
print("        (beta < 0 = higher when FOLH1 is low; contrast is relative to PSMA's own)")
print(M[cols].round(3).to_string())
print("\n[tiers]")
print(M["tier"].value_counts().to_string())
M.to_csv(f"{OUT}/p2_target_prioritisation.csv")
json.dump({"benchmark_contrast": bench,
           "tiers": {k: int(v) for k, v in M["tier"].value_counts().items()},
           "criteria": {"abundance": ">=50% of PSMA-low tumours at log2(FPKM+1)>=2",
                        "selectivity": "tumour-to-normal contrast within 2 log2 of PSMA's own",
                        "specificity": "<25% of single-cell expression in immune/stromal cells",
                        "protein": "HPA IHC detected in >=50% of prostate cancers and plasma-membrane localisation"}},
          open(f"{OUT}/p2_03b_results.json", "w"), indent=1)
print("\nsaved p2 step 3b")
