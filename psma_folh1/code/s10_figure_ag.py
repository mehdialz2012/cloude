"""Step 10 — AlphaGenome validation figure (2 panels)."""
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
matplotlib.rcParams["font.family"] = ["DejaVu Sans"]
matplotlib.rcParams["svg.fonttype"] = "none"
import matplotlib.pyplot as plt
from common import OUT

BLUE, ORANGE, AQUA, RED = "#2a78d6", "#eb6834", "#1baf7a", "#e34948"
INK, INK2, MUTED, GRID, SURFACE = "#0b0b0b", "#52514e", "#8a8984", "#e4e3df", "#fcfcfb"

al = pd.read_csv(f"{OUT}/s09_allele_specificity.csv")
enh = pd.read_csv(f"{OUT}/s09_ism_enhancer.csv").sort_values("pos")
ctl = pd.read_csv(f"{OUT}/s09_ism_control.csv").sort_values("pos")

fig, axes = plt.subplots(1, 2, figsize=(13.4, 5.0), facecolor=SURFACE,
                         gridspec_kw={"width_ratios": [1, 1.35]})
for ax in axes:
    ax.set_facecolor(SURFACE)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    for sp in ("left", "bottom"):
        ax.spines[sp].set_color(GRID)
    ax.tick_params(colors=INK2, labelsize=9, length=3)

# ---------------------------------------------------------------- A
ax = axes[0]
P = 49_208_484
sub = al[al["pos"].between(P - 2, P + 2)].copy()
sub["offset"] = sub["pos"] - P
sub = sub.sort_values(["offset", "alt"])
labels, vals, cols = [], [], []
for _, r in sub.iterrows():
    at_cpg = r["pos"] == P
    labels.append(f"{r['ref']}>{r['alt']}" + ("" if at_cpg else f"  ({r['offset']:+d})"))
    vals.append(r["expr_log2fc"])
    cols.append(ORANGE if at_cpg else MUTED)
y = np.arange(len(vals))[::-1]
ax.barh(y, vals, height=0.62, color=cols, edgecolor=SURFACE, linewidth=1.4)
for yy, v, c in zip(y, vals, cols):
    ax.text(v + (0.008 if v >= 0 else -0.008), yy, f"{v:+.3f}", va="center",
            ha="left" if v >= 0 else "right", fontsize=8,
            color=INK if c == ORANGE else INK2,
            fontweight="bold" if c == ORANGE else "normal")
ax.set_yticks(y)
ax.set_yticklabels(labels, fontsize=8.6, color=INK)
ax.axvline(0, color=INK2, lw=1)
ax.set_xlim(-0.145, 0.37)
ax.set_xlabel("predicted effect on FOLH1 expression (log2FC)", fontsize=9.5, color=INK2)
ax.set_title("A   The cg25346915 result is allele-specific, not CpG disruption",
             loc="left", fontsize=11.5, color=INK, fontweight="bold")
ax.text(0.0, -0.24, "Orange = the interrogated cytosine (chr11:49,208,484); grey = the four "
                    "neighbouring bases. All three substitutions\nat that cytosine abolish the "
                    "CpG, yet C>T gives +0.288, C>G +0.086 and C>A +0.006. AlphaGenome reads DNA\n"
                    "sequence and has no methylation input channel.",
        transform=ax.transAxes, ha="left", va="top", fontsize=8.6, color=MUTED)

# ---------------------------------------------------------------- B
ax = axes[1]
band = ctl["expr_log2fc"].abs().max()
ax.axhspan(-band, band, color=AQUA, alpha=0.16, lw=0,
           label=f"matched peak-free control region\n(86 positions, max |effect| {band:.3f})")
ax.axhline(0, color=INK2, lw=1)
x = enh["pos"].values
v = enh["expr_log2fc"].values
ax.vlines(x, 0, v, color=GRID, lw=1)
ax.scatter(x, v, s=26, c=[ORANGE if abs(t) > 0.05 else BLUE for t in v],
           edgecolor=SURFACE, linewidth=0.8, zorder=4)
core = enh[enh["expr_log2fc"].abs() > 0.09]
ax.annotate(f"sensitive core\n{core['pos'].min():,}–{core['pos'].max():,}",
            xy=(49_178_177, -0.2325), xytext=(49_178_240, -0.185),
            fontsize=8.8, color=INK, ha="left",
            arrowprops=dict(arrowstyle="-|>", mutation_scale=8, lw=0.7, color=INK2))
ax.set_xlabel("chr11 position (hg38) — 4-bp tiling scan across the nominated gene-body element",
              fontsize=9.5, color=INK2)
ax.set_ylabel("predicted effect on FOLH1\nexpression (log2FC)", fontsize=9.5, color=INK2)
ax.set_title("B   The element itself survives a background control the paper lacks",
             loc="left", fontsize=11.5, color=INK, fontweight="bold")
ax.legend(fontsize=8.4, frameon=False, loc="upper left", labelcolor=INK2)
ax.ticklabel_format(axis="x", style="plain", useOffset=False)
ax.set_xticks(np.arange(49_178_000, 49_178_400, 100))
ax.set_xticklabels([f"{t:,}" for t in np.arange(49_178_000, 49_178_400, 100)], fontsize=8.5)
ax.text(0.0, -0.24, "Effect dispersion is 7.1x that of the control region (Levene p = 2.3e-8); "
                    "13 of 86 element positions\nexceed |0.05| versus 0 of 86 controls.",
        transform=ax.transAxes, ha="left", va="top", fontsize=8.6, color=MUTED)

fig.suptitle("AlphaGenome replication: what holds and what does not",
             x=0.012, y=0.995, ha="left", fontsize=13, color=INK, fontweight="bold")
plt.tight_layout(rect=[0, 0.09, 1, 0.945])
for ext in ("png", "svg"):
    fig.savefig(f"{OUT}/fig_alphagenome_validation.{ext}",
                dpi=220 if ext == "png" else None, bbox_inches="tight", facecolor=SURFACE)
print("saved fig_alphagenome_validation.png/.svg")
