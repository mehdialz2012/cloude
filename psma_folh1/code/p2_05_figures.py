"""Phase 2 — figures for the restructured manuscript."""
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
matplotlib.rcParams["font.family"] = ["DejaVu Sans"]
matplotlib.rcParams["svg.fonttype"] = "none"
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from scipy import stats as st

from common import OUT, load_su2c, load_gse126078

BLUE, ORANGE, AQUA, YELLOW = "#2a78d6", "#eb6834", "#1baf7a", "#eda100"
RED, VIOLET = "#e34948", "#6250d6"
INK, INK2, MUTED, GRID, SURFACE = "#0b0b0b", "#52514e", "#8a8984", "#e4e3df", "#fcfcfb"


def style(ax):
    ax.set_facecolor(SURFACE)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(GRID)
    ax.tick_params(colors=INK2, labelsize=8.5, length=3)
    return ax


def save(fig, name):
    for ext in ("png", "svg"):
        fig.savefig(f"{OUT}/{name}.{ext}", dpi=210 if ext == "png" else None,
                    bbox_inches="tight", facecolor=SURFACE)
    print("saved", name)


r1 = json.load(open(f"{OUT}/p2_01_results.json"))
r2 = json.load(open(f"{OUT}/p2_02_results.json"))
r4 = json.load(open(f"{OUT}/p2_04_results.json"))
disc = pd.read_csv(f"{OUT}/p2_discovery_models.csv", index_col=0)
gsea = pd.read_csv(f"{OUT}/p2_gsea_nested.csv", index_col=0)
cov = pd.read_csv(f"{OUT}/p2_su2c_covariates.csv", index_col=0)
cross = pd.read_csv(f"{OUT}/p2_crosscohort_genes.csv", index_col=0)
N = pd.read_csv(f"{OUT}/p2_programme_concordance.csv", index_col=0)
info = pd.read_csv(f"{OUT}/p2_gse_info.csv", index_col=0)
prio = pd.read_csv(f"{OUT}/p2_target_prioritisation.csv", index_col=0)
surf = pd.read_csv(f"{OUT}/p2_surfaceome_full.csv", index_col=0)
dep = pd.read_csv(f"{OUT}/p2_depmap_validated_programme.csv", index_col=0)

# ============================================================ FIGURE 1
fig, ax = plt.subplots(2, 2, figsize=(13.2, 9.0), facecolor=SURFACE)
for a in ax.ravel():
    style(a)

a = ax[0, 0]
f = cov["FOLH1"]
a.hist(f, bins=38, color=BLUE, alpha=.85, edgecolor=SURFACE, linewidth=.6)
for q, lab in [(.25, "Q1"), (.5, "median"), (.75, "Q3")]:
    a.axvline(f.quantile(q), color=INK2, lw=1, ls=(0, (4, 3)))
    a.text(f.quantile(q), a.get_ylim()[1] * .96, f" {lab}", fontsize=8, color=INK2, va="top")
a.set_xlabel("FOLH1, log2(FPKM+1) — SU2C mCRPC (n=266)", fontsize=9.5, color=INK2)
a.set_ylabel("tumours", fontsize=9.5, color=INK2)
a.set_title("a   FOLH1 is continuous, with a low tail", loc="left",
            fontsize=11, color=INK, fontweight="bold")
a.text(0, -.23, f"SD 2.71 · Q1–Q3 span {2**(f.quantile(.75)-f.quantile(.25)):.0f}-fold · "
                "quartile, tercile, median and decile cut-offs all agree\nwith the continuous "
                "model (ρ 0.91–0.94), so the cut-off is not what drives the result",
       transform=a.transAxes, ha="left", va="top", fontsize=8.3, color=MUTED)

a = ax[0, 1]
items = [("biopsy site (liver etc.)", abs(r1["rho_folh1_hepatocyte"])),
         ("osteoblast content", abs(r1["rho_folh1_osteoblast"])),
         ("immune content", abs(r1["rho_folh1_immune"])),
         ("stromal content", abs(r1["rho_folh1_stromal"])),
         ("neuroendocrine score", abs(r1["rho_folh1_NE"])),
         ("AR activity score", abs(r1["rho_folh1_AR"]))]
items.sort(key=lambda t: t[1])
y = np.arange(len(items))
cols = [ORANGE if "AR" in t[0] or "neuroendo" in t[0] else MUTED for t in items]
a.barh(y, [t[1] for t in items], height=.6, color=cols, edgecolor=SURFACE, linewidth=1.4)
for yy, (lab, v) in zip(y, items):
    a.text(v + .008, yy, f"{v:.2f}", va="center", fontsize=8.5, color=INK2)
a.set_yticks(y)
a.set_yticklabels([t[0] for t in items], fontsize=9, color=INK)
a.set_xlim(0, .66)
a.set_xlabel("|Spearman ρ| with FOLH1", fontsize=9.5, color=INK2)
a.set_title("b   What FOLH1 tracks in a bulk biopsy", loc="left",
            fontsize=11, color=INK, fontweight="bold")

a = ax[1, 0]
steps = ["raw\nassociation", "+ composition\n+ biopsy site", "+ neuroendocrine\ndifferentiation"]
vals = [r1["n_sig_M1"], r1["n_sig_M2"], r1["n_sig_M3"]]
bars = a.bar(range(3), vals, color=[MUTED, MUTED, BLUE], width=.6,
             edgecolor=SURFACE, linewidth=1.6)
for i, v in enumerate(vals):
    a.text(i, v + 90, f"{v:,}", ha="center", fontsize=10, color=INK, fontweight="bold")
for i in range(2):
    lost = vals[i] - vals[i + 1]
    a.annotate(f"−{lost:,}", xy=(i + .5, (vals[i] + vals[i + 1]) / 2),
               ha="center", fontsize=9, color=RED)
a.set_xticks(range(3))
a.set_xticklabels(steps, fontsize=8.8, color=INK)
a.set_ylabel("genes associated with FOLH1 (FDR<0.05)", fontsize=9.5, color=INK2)
a.set_ylim(0, 7400)
a.set_title("c   Three quarters of the association is confounding", loc="left",
            fontsize=11, color=INK, fontweight="bold")

a = ax[1, 1]
sel = ["E2F_TARGETS", "G2M_CHECKPOINT", "ANDROGEN_RESPONSE", "OXIDATIVE_PHOSPHORYLATION",
       "EPITHELIAL_MESENCHYMAL_TRANSITION", "INFLAMMATORY_RESPONSE", "HYPOXIA"]
sel = [s for s in sel if s in gsea.index]
w = .38
xs = np.arange(len(sel))
a.barh(xs + w / 2, gsea.loc[sel, "NES_M1"], w, color=MUTED, edgecolor=SURFACE, linewidth=1.3,
       label="raw association")
a.barh(xs - w / 2, gsea.loc[sel, "NES_M3"], w, color=BLUE, edgecolor=SURFACE, linewidth=1.3,
       label="fully adjusted")
a.axvline(0, color=INK2, lw=1)
a.set_yticks(xs)
a.set_yticklabels([s.replace("_", " ").title().replace("E2F", "E2F").replace("G2M", "G2M")
                   for s in sel], fontsize=8.6, color=INK)
a.set_xlabel("GSEA normalised enrichment score (positive = up in PSMA-low)",
             fontsize=9.3, color=INK2)
a.legend(fontsize=8.4, frameon=False, labelcolor=INK2, loc="upper left")
a.set_title("d   EMT, inflammation and hypoxia are composition, not tumour", loc="left",
            fontsize=11, color=INK, fontweight="bold")
fig.suptitle("Figure 1 · Defining the PSMA-low phenotype in metastatic disease",
             x=.012, y=1.0, ha="left", fontsize=13, color=INK, fontweight="bold")
plt.tight_layout(rect=[0, 0, 1, .965])
save(fig, "p2_fig1_discovery")

# ============================================================ FIGURE 2
fig, ax = plt.subplots(2, 2, figsize=(13.2, 9.0), facecolor=SURFACE)
for a in ax.ravel():
    style(a)

a = ax[0, 0]
gse = load_gse126078()[info.index]
fg = gse.loc["FOLH1"]
order = ["ARpos_NEpos", "ARpos_NEneg", "ARlow_NEneg", "ARneg_NEpos", "ARneg_NEneg"]
data = [fg[info[info["phenotype"] == p].index].values for p in order]
bp = a.boxplot(data, vert=True, widths=.55, patch_artist=True, showfliers=False)
for p, box in zip(order, bp["boxes"]):
    box.set(facecolor=BLUE if p.startswith("ARpos") else ORANGE, alpha=.75,
            edgecolor=INK2, linewidth=.8)
for med in bp["medians"]:
    med.set(color=INK, linewidth=1.4)
for i, d in enumerate(data):
    a.scatter(np.random.normal(i + 1, .07, len(d)), d, s=12, color=INK2, alpha=.5, zorder=3)
a.axhline(fg.quantile(.25), color=RED, lw=1.2, ls=(0, (4, 3)))
a.text(5.45, fg.quantile(.25), " cohort Q1", fontsize=8, color=RED, va="center")
a.set_xticklabels([f"{p}\n(n={len(d)})" for p, d in zip(order, data)], fontsize=8.2, color=INK)
a.set_ylabel("FOLH1, log2(FPKM+1)", fontsize=9.5, color=INK2)
a.set_title("a   PSMA-low is overwhelmingly AR-low/negative disease", loc="left",
            fontsize=11, color=INK, fontweight="bold")
a.text(0, -.19, "Only 6 of 60 AR-positive, non-NE tumours fall in the cohort's bottom quartile:\n"
                "the PSMA-low compartment is, in the main, AR-low and neuroendocrine disease.",
       transform=a.transAxes, ha="left", va="top", fontsize=8.3, color=MUTED)

a = ax[0, 1]
sub = cross.dropna(subset=["su2c_M3", "gse_patient_M3"])
a.scatter(sub["su2c_M3"], sub["gse_patient_M3"], s=5, color=MUTED, alpha=.25, lw=0)
rep = sub[sub["reproducible"]]
a.scatter(rep["su2c_M3"], rep["gse_patient_M3"], s=9, color=BLUE, alpha=.8, lw=0,
          label=f"reproduced (n={len(rep)})")
a.axhline(0, color=INK2, lw=.8)
a.axvline(0, color=INK2, lw=.8)
rho = st.spearmanr(sub["su2c_M3"], sub["gse_patient_M3"])[0]
a.set_xlim(-0.75, 1.35)
a.set_ylim(-0.75, 1.35)
a.set_xlabel("SU2C effect (β per SD of FOLH1, fully adjusted)", fontsize=9.3, color=INK2)
a.set_ylabel("GSE126078 patient-level effect", fontsize=9.3, color=INK2)
a.set_title(f"b   Independent validation: ρ = {rho:.2f}", loc="left",
            fontsize=11, color=INK, fontweight="bold")
a.legend(fontsize=8.4, frameon=False, labelcolor=INK2, loc="upper left")

a = ax[1, 0]
labels = ["raw\n(SU2C)", "adjusted\n(SU2C)", "reproduced in\nGSE126078",
          "also at\npatient level", "also in\nBeltran"]
vals = [r1["n_sig_M1"], r2["n_su2c_robust"], r2["n_rep_gse"], r2["n_rep_patient"],
        r2["n_rep_all3"]]
a.bar(range(5), vals, color=[MUTED, MUTED, BLUE, BLUE, AQUA], width=.62,
      edgecolor=SURFACE, linewidth=1.6)
for i, v in enumerate(vals):
    a.text(i, v + 90, f"{v:,}", ha="center", fontsize=9.5, color=INK, fontweight="bold")
a.set_xticks(range(5))
a.set_xticklabels(labels, fontsize=8.2, color=INK)
a.set_ylabel("genes", fontsize=9.5, color=INK2)
a.set_ylim(0, 7400)
a.set_title("c   From 6,438 associations to a 251-gene validated programme", loc="left",
            fontsize=11, color=INK, fontweight="bold")
a.text(.98, .55, "(the adjusted count is restricted to\ngenes measured in all three cohorts)",
       transform=a.transAxes, ha="right", fontsize=8.2, color=MUTED)

a = ax[1, 1]
keep = ["E2F_TARGETS", "G2M_CHECKPOINT", "ANDROGEN_RESPONSE", "OXIDATIVE_PHOSPHORYLATION",
        "PROTEIN_SECRETION", "EPITHELIAL_MESENCHYMAL_TRANSITION"]
keep = [k for k in keep if k in N.index]
cohorts = ["SU2C", "GSE126078", "Beltran", "TCGA_primary"]
M = N.loc[keep, cohorts]
im = a.imshow(M.values, cmap="RdBu_r", vmin=-3.2, vmax=3.2, aspect="auto")
a.set_xticks(range(len(cohorts)))
a.set_xticklabels(["SU2C\n(discovery)", "GSE126078\n(validation)", "Beltran\n(NE focus)",
                   "TCGA\n(primary)"], fontsize=8.2, color=INK)
a.set_yticks(range(len(keep)))
a.set_yticklabels([k.replace("_", " ").title() for k in keep], fontsize=8.4, color=INK)
for i in range(len(keep)):
    for j in range(len(cohorts)):
        v = M.values[i, j]
        a.text(j, i, f"{v:+.1f}", ha="center", va="center", fontsize=8.2,
               color="white" if abs(v) > 1.9 else INK)
plt.colorbar(im, ax=a, label="NES", shrink=.8)
a.set_title("d   Only proliferation and AR loss travel across cohorts", loc="left",
            fontsize=11, color=INK, fontweight="bold")
fig.suptitle("Figure 2 · Independent validation in metastatic cohorts",
             x=.012, y=1.0, ha="left", fontsize=13, color=INK, fontweight="bold")
plt.tight_layout(rect=[0, 0, 1, .965])
save(fig, "p2_fig2_validation")

# ============================================================ FIGURE 3
fig = plt.figure(figsize=(13.4, 9.6), facecolor=SURFACE)
gs = fig.add_gridspec(2, 2, height_ratios=[1, 1.15], hspace=.42, wspace=.26)

a = style(fig.add_subplot(gs[0, 0]))
counts = surf["state"].value_counts()
ret = surf.groupby("state")["retained"].sum()
order_s = ["enriched", "unchanged", "lost"]
x = np.arange(3)
a.bar(x, [counts.get(s, 0) for s in order_s], width=.6, color=MUTED,
      edgecolor=SURFACE, linewidth=1.6, label="all surface genes")
a.bar(x, [ret.get(s, 0) for s in order_s], width=.6, color=BLUE,
      edgecolor=SURFACE, linewidth=1.6, label="expressed in ≥50% of PSMA-low tumours")
for i, s in enumerate(order_s):
    a.text(i, counts.get(s, 0) + 22, f"{counts.get(s,0)}", ha="center", fontsize=9, color=INK)
    a.text(i, max(ret.get(s, 0) / 2, 18), f"{int(ret.get(s,0))}", ha="center", fontsize=8.6,
           color="white" if ret.get(s, 0) > 60 else INK, fontweight="bold")
a.set_xticks(x)
a.set_xticklabels(["enriched\nin PSMA-low", "unchanged", "lost\nin PSMA-low"],
                  fontsize=8.8, color=INK)
a.set_ylabel("cell-surface genes", fontsize=9.5, color=INK2)
a.legend(fontsize=8.2, frameon=False, labelcolor=INK2)
a.set_title("a   Most surface antigens are retained, not lost", loc="left",
            fontsize=11, color=INK, fontweight="bold")

a = style(fig.add_subplot(gs[0, 1]))
P = prio.dropna(subset=["contrast_vs_PSMA_benchmark", "su2c_low_median"])
P = P[P["pass_membrane"]]
colmap = {"A": AQUA, "B": BLUE, "C": MUTED, "D": ORANGE, "E": GRID}
cols = [colmap.get(str(t)[0], MUTED) for t in P["tier"]]
a.scatter(P["contrast_vs_PSMA_benchmark"], P["su2c_low_median"], s=52, c=cols,
          edgecolor=SURFACE, linewidth=1.2, zorder=3)
a.axvline(-2, color=INK2, lw=1, ls=(0, (4, 3)))
a.axhline(2, color=INK2, lw=1, ls=(0, (4, 3)))
a.text(-2.15, 1.15, "selectivity threshold\n(within 4-fold of PSMA)", fontsize=7.6,
       color=INK2, ha="right")
for g in ["CD276", "ABCB6", "STEAP1", "STEAP2", "TACSTD2", "SLC7A5", "SLC16A3", "PLAUR",
          "GPC2", "FOLH1", "MUC1", "EPCAM"]:
    if g in P.index:
        a.annotate(g, (P.loc[g, "contrast_vs_PSMA_benchmark"], P.loc[g, "su2c_low_median"]),
                   xytext=(5, 4), textcoords="offset points", fontsize=8.2, color=INK)
a.set_xlabel("tumour-to-normal contrast, relative to PSMA's own (log2)", fontsize=9.3, color=INK2)
a.set_ylabel("abundance in PSMA-low tumours\nmedian log2(FPKM+1)", fontsize=9.3, color=INK2)
a.set_title("b   Abundance and selectivity, side by side", loc="left",
            fontsize=11, color=INK, fontweight="bold")
a.legend(handles=[Patch(color=AQUA, label="A · meets all criteria"),
                  Patch(color=MUTED, label="C · retained, poor contrast"),
                  Patch(color=ORANGE, label="D · enriched, unsuitable"),
                  Patch(color=GRID, label="E · not supported")],
         fontsize=7.8, frameon=False, labelcolor=INK2, loc="lower left")

a = style(fig.add_subplot(gs[1, :]))
panel = ["CD276", "STEAP1", "STEAP2", "TACSTD2", "SLC44A4", "CD46", "EPCAM", "ERBB2",
         "SLC7A5", "SLC16A3", "ATP2B4", "GPC2", "MUC1", "PLAUR", "PSCA", "FOLH1"]
panel = [g for g in panel if g in prio.index]
crit = ["pass_abundance", "pass_selectivity", "pass_specificity", "pass_protein",
        "pass_membrane"]
names = ["abundance\nin PSMA-low", "selectivity\nvs normal organs", "tumour-cell\nspecificity",
         "protein detected\n(HPA IHC)", "plasma-membrane\nevidence"]
Mx = prio.loc[panel, crit].astype(float).values.T
a.imshow(Mx, cmap=matplotlib.colors.ListedColormap([GRID, AQUA]), aspect="auto", vmin=0, vmax=1)
a.set_xticks(range(len(panel)))
a.set_xticklabels(panel, fontsize=8.6, color=INK, rotation=42, ha="right")
a.set_yticks(range(len(crit)))
a.set_yticklabels(names, fontsize=8.4, color=INK)
for i in range(len(crit)):
    for j in range(len(panel)):
        a.text(j, i, "✓" if Mx[i, j] == 1 else "✗", ha="center", va="center",
               fontsize=9.5, color="white" if Mx[i, j] == 1 else MUTED)
for j, g in enumerate(panel):
    a.text(j, -0.78, str(prio.loc[g, "tier"])[0], ha="center", fontsize=9,
           color=INK, fontweight="bold")
a.text(-0.9, -0.78, "tier", ha="right", fontsize=8.4, color=INK2)
a.set_title("c   Pre-specified suitability criteria, established targets and new nominations side by side",
            loc="left", fontsize=11, color=INK, fontweight="bold", pad=22)
fig.suptitle("Figure 3 · Which surface targets could actually cover PSMA-low disease",
             x=.012, y=1.0, ha="left", fontsize=13, color=INK, fontweight="bold")
save(fig, "p2_fig3_targets")

# ============================================================ FIGURE 4
fig, ax = plt.subplots(1, 3, figsize=(15.0, 4.6), facecolor=SURFACE)
for a in ax:
    style(a)

a = ax[0]
genes = ["SEPSECS", "TXNRD1", "EEFSEC", "PSTK", "GPX4"]
old = [1, 2, 4, 15, 13]
new = [r4[f"dep_{g}"]["rank_pan"] for g in genes]
x = np.arange(len(genes))
a.scatter(x, old, s=70, color=ORANGE, edgecolor=SURFACE, linewidth=1.3, zorder=4,
          label="discovery-only signature (as published)")
a.scatter(x, new, s=70, color=BLUE, edgecolor=SURFACE, linewidth=1.3, zorder=4,
          label="cross-cohort validated programme")
for xi, (o, n) in enumerate(zip(old, new)):
    a.plot([xi, xi], [o, n], color=GRID, lw=1.6, zorder=2)
a.set_yscale("log")
a.set_ylim(.7, 2e4)
a.invert_yaxis()
a.set_xticks(x)
a.set_xticklabels(genes, fontsize=8.8, color=INK)
a.set_ylabel("dependency rank among 17,293\nnon-essential genes", fontsize=9.3, color=INK2)
a.legend(fontsize=8, frameon=False, labelcolor=INK2, loc="lower left")
a.set_title("a   The selenocysteine axis does not\n      survive validation", loc="left",
            fontsize=10.5, color=INK, fontweight="bold")

a = ax[1]
nn = dep[~dep["common_essential"]].sort_values("rho").head(10)
y = np.arange(len(nn))[::-1]
a.barh(y, -nn["rho"], height=.62, color=BLUE, edgecolor=SURFACE, linewidth=1.4)
a.set_yticks(y)
a.set_yticklabels(nn.index, fontsize=8.6, color=INK)
for yy, v in zip(y, -nn["rho"]):
    a.text(v + .004, yy, f"{v:.2f}", va="center", fontsize=8, color=INK2)
a.set_xlabel("dependency enrichment in programme-high cells (−ρ)", fontsize=9, color=INK2)
a.set_title("b   What the validated programme\n      actually nominates", loc="left",
            fontsize=10.5, color=INK, fontweight="bold")
a.text(0, -.17, "Purine and one-carbon synthesis, and p53 status — i.e. proliferation,\n"
                "not a selective vulnerability of PSMA-low disease.",
       transform=a.transAxes, ha="left", va="top", fontsize=8, color=MUTED)

a = ax[2]
tt = pd.read_csv(f"{OUT}/p2_tcga_models.csv", index_col=0)
ov = disc.index.intersection(tt.index)
a.scatter(disc.loc[ov, "beta_M3"], tt.loc[ov, "beta_M3"], s=4, color=MUTED, alpha=.2, lw=0)
prog = pd.read_csv(f"{OUT}/p2_validated_programme.csv")["gene"]
pg = [g for g in prog if g in ov]
a.scatter(disc.loc[pg, "beta_M3"], tt.loc[pg, "beta_M3"], s=11, color=BLUE, alpha=.85, lw=0,
          label=f"validated programme (n={len(pg)})")
a.axhline(0, color=INK2, lw=.8)
a.axvline(0, color=INK2, lw=.8)
a.set_xlim(-0.8, 1.4)
a.set_ylim(-0.6, 1.0)
a.set_xlabel("effect in mCRPC (SU2C)", fontsize=9.3, color=INK2)
a.set_ylabel("effect in primary disease (TCGA)", fontsize=9.3, color=INK2)
a.legend(fontsize=8, frameon=False, labelcolor=INK2, loc="upper left")
a.set_title(f"c   The programme is already present\n      in primary tumours (ρ = "
            f"{r4['tcga_concordance']:.2f})", loc="left",
            fontsize=10.5, color=INK, fontweight="bold")
a.text(0, -.17, f"{100*r4['tcga_same_direction_frac']:.0f}% of validated-programme genes keep their\n"
                "direction in treatment-naive primary disease.",
       transform=a.transAxes, ha="left", va="top", fontsize=8, color=MUTED)
fig.suptitle("Figure 4 · Exploratory dependencies, and primary versus metastatic disease",
             x=.008, y=1.03, ha="left", fontsize=13, color=INK, fontweight="bold")
plt.tight_layout(rect=[0, .06, 1, .95])
save(fig, "p2_fig4_depmap_tcga")
print("all figures done")
