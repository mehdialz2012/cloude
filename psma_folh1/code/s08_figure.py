"""Step 8 — one summary figure for the four validation tests that change the conclusions."""
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
matplotlib.rcParams["font.family"] = ["DejaVu Sans"]
matplotlib.rcParams["svg.fonttype"] = "none"
import matplotlib.pyplot as plt

from common import OUT

# validated categorical slots 1-3 (light surface) + diverging pair + neutrals
BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"
RED = "#e34948"
INK, INK2, MUTED = "#0b0b0b", "#52514e", "#8a8984"
GRID = "#e4e3df"
SURFACE = "#fcfcfb"

s03 = json.load(open(f"{OUT}/s03_results.json"))
s04 = json.load(open(f"{OUT}/s04_results.json"))
s05 = json.load(open(f"{OUT}/s05_results.json"))
nullc = pd.read_csv(f"{OUT}/s03_circularity_null.csv")

fig, axes = plt.subplots(2, 2, figsize=(14.0, 11.2), facecolor=SURFACE)
for ax in axes.ravel():
    ax.set_facecolor(SURFACE)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(GRID)
    ax.tick_params(colors=INK2, labelsize=9, length=3)

# ---------------------------------------------------------------- A
ax = axes[0, 0]
items = [("Neuroendocrine score\n(12 genes, standard panel)", s03["B2_beltran_auc::NE score (12 genes)"], AQUA),
         ("PSMA-low signature\n(100 genes, this paper)", s03["B2_beltran_auc::PSMA-low signature (100 genes)"], BLUE),
         ("AR activity score\n(20 genes, inverted)", s03["B2_beltran_auc::AR score (20 genes, inverted)"], MUTED),
         ("FOLH1 expression alone\n(inverted)", s03["B2_beltran_auc::FOLH1 alone (inverted)"], MUTED)]
items.sort(key=lambda t: t[1])
ypos = np.arange(len(items))
for i, (lab, v, col) in enumerate(items):
    ax.barh(i, v - 0.5, left=0.5, height=0.52, color=col,
            edgecolor=SURFACE, linewidth=2)
    ax.text(v + 0.006, i, f"{v:.3f}", va="center", ha="left", fontsize=10,
            color=INK, fontweight="bold")
ax.set_yticks(ypos)
ax.set_yticklabels([t[0] for t in items], fontsize=9, color=INK)
ax.set_xlim(0.5, 1.06)
ax.set_xlabel("AUC, NEPC vs CRPC-adenocarcinoma (Beltran, n=49)", fontsize=9.5, color=INK2)
ax.axvline(0.5, color=GRID, lw=1)
ax.set_title("A   The external validation is a lineage read-out",
             loc="left", fontsize=11.5, color=INK, fontweight="bold")
ax.text(0.0, -0.26, "A standard 12-gene NE panel classifies the same samples better than the "
                    "100-gene signature.\nRobust to how this cohort's negative "
                    "expression floor is handled.",
        transform=ax.transAxes, ha="left", va="top", fontsize=8.6, color=MUTED)

# ---------------------------------------------------------------- B
ax = axes[0, 1]
vals = nullc["rho_gse_patient"].abs().values
rng = np.random.default_rng(3)
ax.scatter(vals, 0.04 * rng.normal(size=len(vals)), s=52, color=MUTED, alpha=.75,
           edgecolor=SURFACE, linewidth=1.5, zorder=3,
           label=f"signatures built the same way\nfrom {len(vals)} random anchor genes")
obs = abs(s03["A3_gse"]["rho_patient"])
ax.scatter([obs], [0], s=190, marker="D", color=ORANGE, edgecolor=SURFACE, linewidth=2,
           zorder=5, label="FOLH1-anchored signature\n(the published result)")
ax.annotate(f"FOLH1\nρ = −{obs:.3f}", xy=(obs, 0), xytext=(obs - 0.015, 0.16),
            ha="center", fontsize=9.5, color=INK, fontweight="bold",
            arrowprops=dict(arrowstyle="-", lw=0.8, color=INK2))
med = np.median(vals)
ax.axvline(med, color=MUTED, lw=1.2, ls=(0, (4, 3)))
ax.text(med, -0.19, f"median of the null\n|ρ| = {med:.3f}", ha="center", va="top",
        fontsize=8.5, color=MUTED)
ax.set_ylim(-0.26, 0.26)
ax.set_yticks([])
ax.set_xlim(0.3, 0.95)
ax.spines["left"].set_visible(False)
ax.set_xlabel("|Spearman ρ| between the signature and its own anchor gene\n"
              "in GSE126078 patient samples", fontsize=9.5, color=INK2)
ax.set_title("B   Anti-correlating with FOLH1 in a second cohort is not validation",
             loc="left", fontsize=11.5, color=INK, fontweight="bold")
ax.legend(fontsize=8.5, frameon=False, loc="upper left", labelcolor=INK2, handletextpad=.6)
ax.text(0.0, -0.26, "Any gene's quartile-derived signature anti-correlates with that gene in another\n"
                    "cohort. The published value sits inside that null.",
        transform=ax.transAxes, ha="left", va="top", fontsize=8.6, color=MUTED)

# ---------------------------------------------------------------- C
ax = axes[1, 0]
b2 = pd.DataFrame(s04["B2"]).T
order = ["CD276", "GPC2", "STEAP1", "STEAP2", "SLC12A4", "SLC16A3", "ATP2B4", "ADAM19",
         "SLC7A5", "MRC2", "TACSTD2", "NCAM1", "PLAUR", "MUC1"]
order = [g for g in order if g in b2.index]
labmap = {"CD276": "CD276 (B7-H3)", "TACSTD2": "TACSTD2 (TROP2)", "SLC16A3": "SLC16A3 (MCT4)",
          "SLC7A5": "SLC7A5 (LAT1)", "PLAUR": "PLAUR (uPAR)", "MRC2": "MRC2 (uPARAP)",
          "NCAM1": "NCAM1 (CD56)"}
v = list(b2.loc[order, "contrast_low_vs_vital"].astype(float).values)
ref = float(b2.loc["FOLH1", "contrast_high_vs_vital"])
order = order + ["FOLH1 in PSMA-high\n(reference target)"]
v = np.array(v + [ref])
y = np.arange(len(order))[::-1]
cols = [BLUE if x > 0 else RED for x in v]
ax.barh(y, v, height=0.6, color=cols, edgecolor=SURFACE, linewidth=1.6)
for yy, xx in zip(y, v):
    ax.text(xx + (0.18 if xx > 0 else -0.18), yy, f"{xx:+.1f}", va="center",
            ha="left" if xx > 0 else "right", fontsize=8.5, color=INK2)
ax.set_yticks(y)
ax.set_yticklabels([labmap.get(g, g) for g in order], fontsize=9, color=INK)
ax.axvline(0, color=INK2, lw=1)
ax.set_xlim(-8.6, 4.2)
ax.set_xlabel("tumour-to-normal contrast, log2 units\n"
              "(median tumour expression − highest vital-organ expression, HPA)",
              fontsize=9.5, color=INK2)
ax.set_title("C   Every nominated new target is out-expressed by a vital organ",
             loc="left", fontsize=11.5, color=INK, fontweight="bold")
ax.text(0.0, -0.30, "Only CD276 (B7-H3) keeps positive contrast in PSMA-low tumours. "
                    "The three best-expressed\nsurface targets there are TROP2, B7-H3 and STEAP1 "
                    "- the luminal antigens the paper rules out.",
        transform=ax.transAxes, ha="left", va="top", fontsize=8.6, color=MUTED)

# ---------------------------------------------------------------- D
ax = axes[1, 1]
genes = ["TXNRD1", "SEPSECS", "EEFSEC", "PSTK", "GPX4", "GAPDH"]
pan = [-s05["A2"]["axis"][g] for g in genes]
wl = [-s05[f"B2_rho_{g}"] for g in genes]
x = np.arange(len(genes))
w = 0.36
ax.bar(x - w / 2, pan, w, color=BLUE, edgecolor=SURFACE, linewidth=1.6,
       label="pan-lineage (as published)")
ax.bar(x + w / 2, wl, w, color=ORANGE, edgecolor=SURFACE, linewidth=1.6,
       label="within-lineage (lineage-demeaned)")
rk = [s05[f"B2_rank_{g}"] for g in genes]
for xi, (a, b, r) in enumerate(zip(pan, wl, rk)):
    ax.text(xi - w / 2, a + .004, f"{a:.2f}", ha="center", fontsize=8, color=INK2)
    ax.text(xi + w / 2, b + .004, f"{b:.2f}  (rank {r:,})", ha="center", fontsize=7.8,
            color=ORANGE, rotation=90, va="bottom")
# sign convention: the axis is enrichment in PSMA-LOW-like cells. FOLH1 mRNA is the
# inverse of "PSMA-low", so its enrichment value is rho(FOLH1, gene effect) itself.
fol = s05["B4"]["table"]["TXNRD1"]["rho_FOLH1_vs_dependency"]
ax.axhline(0, color=INK2, lw=1)
ax.axhline(fol, color=MUTED, lw=1.4, ls=(0, (5, 3)))
ax.text(-0.45, fol - .008, f"measured against FOLH1 mRNA itself: {fol:+.2f}",
        ha="left", va="top", fontsize=8.4, color=MUTED)
ax.set_xticks(x)
ax.set_xticklabels(genes, fontsize=9, color=INK)
ax.set_ylabel("dependency enrichment in\nPSMA-low-like cells (−ρ)", fontsize=9.5, color=INK2)

ax.set_ylim(-0.155, 0.36)
ax.set_title("D   The selenocysteine 'axis' collapses to TXNRD1 alone",
             loc="left", fontsize=11.5, color=INK, fontweight="bold")
ax.legend(fontsize=8.5, frameon=False, loc="upper right", labelcolor=INK2)
ax.text(0.0, -0.30, "Lineage alone explains 29.8% of the score across DepMap; prostate is the "
                    "lowest of 23 lineages.\nTested directly against FOLH1 mRNA the correlation is "
                    "near zero and of the wrong sign (dashed line).",
        transform=ax.transAxes, fontsize=8.6, color=MUTED, va="top", ha="left")

fig.suptitle("Independent validation: the four tests that change the conclusions",
             x=0.012, y=0.985, ha="left", fontsize=14, color=INK, fontweight="bold")
fig.text(0.012, 0.957, "Every published number replicated exactly (claim_by_claim_verdicts.csv). "
                       "These are new tests, not failures to reproduce.",
         ha="left", fontsize=9.5, color=MUTED)
plt.tight_layout(rect=[0, 0.03, 1, 0.935])
fig.subplots_adjust(hspace=0.62, wspace=0.30)
fig.text(0.012, 0.008, "Replication data: TCGA-PRAD (Xena) · SU2C/PCF prad_su2c_2019 · Beltran "
                       "nepc_wcm_2016 · GSE126078 · DepMap 24Q4 · HPA v24 · CollecTRI. "
                       "Re-downloaded from primary sources and re-analysed with independent code.",
         fontsize=7.6, color=MUTED)
for ext in ("png", "svg"):
    fig.savefig(f"{OUT}/fig_validation_summary.{ext}", dpi=220 if ext == "png" else None,
                bbox_inches="tight", facecolor=SURFACE)
print("saved fig_validation_summary.png/.svg")
