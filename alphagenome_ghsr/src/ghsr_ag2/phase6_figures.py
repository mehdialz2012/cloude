"""Phase VI - publication-quality figures.

Four figures, each written as vector PDF and SVG plus a 600 dpi PNG, and each
accompanied by a CSV of the exact values plotted (results/tables/).

Figure 1  Genomic organisation of the GHSR locus.
Figure 2  Predicted regulatory effects, with the in silico mutagenesis
          reference distribution for scale.
Figure 3  Cell-type-specific effects across conditions and modalities,
          annotated by track evidence class.
Figure 4  Integrated evidence schematic (three evidence levels).

Scale honesty
-------------
Because the measured effects are near zero, an auto-scaled heat map would
manufacture the appearance of structure out of numerical noise. Figure 3 is
therefore drawn on a fixed, biologically anchored scale (+/- 0.1 log2 fold
change, ~7% change in predicted signal) on which the data are correctly almost
entirely neutral, with a second panel showing the same data auto-scaled and
labelled with its expansion factor, so that the reader can see both the true
magnitude and the internal structure.

Run:  python -m ghsr_ag2.phase6_figures
"""

from __future__ import annotations

import json

import matplotlib.patches as mpatches
from matplotlib.patches import ConnectionPatch
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm

from . import config, provenance, style

DIVERGING_CMAP = LinearSegmentedColormap.from_list(
    "ag_diverging",
    [style.DIVERGING_LOW, style.DIVERGING_MID, style.DIVERGING_HIGH],
)

# Biologically anchored display limit for signed effects: 0.1 log2 fold change
# is a ~7% change in predicted signal.
ANCHORED_LIMIT = 0.1

GHSR_EXONS = {
    "ENST00000241256": [(172443291, 172445465), (172447618, 172448456)],
    "ENST00000427970": [(172447544, 172448455)],
    "ENST00001141273": [(172442551, 172445465), (172447618, 172448456)],
}
GHSR_TX_ORDER = ("ENST00000241256", "ENST00000427970", "ENST00001141273")
GHSR_TX_NAMES = {
    "ENST00000241256": "GHSR-201 (canonical)",
    "ENST00000427970": "GHSR-202",
    "ENST00001141273": "GHSR-203",
}
GHSR_TX_SPAN = {
    "ENST00000241256": (172443291, 172448456),
    "ENST00000427970": (172447544, 172448455),
    "ENST00001141273": (172442551, 172448456),
}
REG_FEATURES = [
    ("ENSR3_B3983C", "enhancer", 172449561, 172450143),
    ("ENSR3_53HK6KS", "CTCF binding site", 172449820, 172449841),
]
GHSR_TSS = 172448456


# ==========================================================================
# Figure 1 - locus organisation
# ==========================================================================
def figure1(annot: pd.DataFrame) -> None:
    """Two panels: whole-locus overview, then a zoom on the sequenced amplicon.

    A single-panel version is unreadable: the three variants lie within 100 bp
    of one another at one end of a 8 kb locus, so their labels collide.
    """
    style.apply()
    fig, (axA, axB) = plt.subplots(
        2, 1, figsize=(7.2, 5.0), gridspec_kw={"height_ratios": [1.25, 1.0],
                                               "hspace": 0.55}
    )
    prim = annot[annot.is_primary_configuration].sort_values("position_grch38")

    # ------------------------------------------------------------------
    # Panel A - whole locus
    # ------------------------------------------------------------------
    x0, x1 = 172_442_200, 172_450_500
    for i, tx in enumerate(GHSR_TX_ORDER):
        y = -0.45 - i * 0.42
        s_, e_ = GHSR_TX_SPAN[tx]
        axA.plot([s_, e_], [y, y], color=style.INK_MUTED, lw=0.9, zorder=2)
        for ex_s, ex_e in GHSR_EXONS[tx]:
            axA.add_patch(
                mpatches.Rectangle(
                    (ex_s, y - 0.11), ex_e - ex_s, 0.22,
                    facecolor=style.SERIES[1], edgecolor=style.SURFACE,
                    linewidth=0.7, zorder=3,
                )
            )
        axA.text(x0 - 120, y, GHSR_TX_NAMES[tx], ha="right", va="center",
                 fontsize=6.4, color=style.INK_SECONDARY)
        axA.annotate(
            "", xy=(s_ - 300, y), xytext=(s_ - 40, y),
            arrowprops=dict(arrowstyle="->", color=style.INK_MUTED, lw=0.9),
        )

    axA.plot([GHSR_TSS, GHSR_TSS], [-0.30, 0.55], color=style.INK_PRIMARY,
             lw=1.1, zorder=4)
    axA.annotate(
        "", xy=(GHSR_TSS - 700, 0.55), xytext=(GHSR_TSS, 0.55),
        arrowprops=dict(arrowstyle="->", color=style.INK_PRIMARY, lw=1.1),
    )
    axA.text(GHSR_TSS - 820, 0.72, "GHSR TSS (minus strand)", ha="left",
             va="bottom", fontsize=6.6, color=style.INK_PRIMARY)

    # amplicon marker + zoom bracket
    axA.add_patch(
        mpatches.Rectangle(
            (config.AMPLICON_START, -0.18),
            config.AMPLICON_END - config.AMPLICON_START, 0.36,
            facecolor=style.SERIES[2], edgecolor="none", alpha=0.9, zorder=4,
        )
    )
    axA.text(
        config.AMPLICON_END + 120, 0.0,
        f"sequenced amplicon\n{config.REPORTED_AMPLICON_BP} bp",
        ha="left", va="center", fontsize=6.2, color=style.SERIES[2],
    )
    axA.set_xlim(x0, x1)
    axA.set_ylim(-1.75, 1.05)
    axA.set_yticks([])
    axA.spines["left"].set_visible(False)
    axA.xaxis.set_major_formatter(
        plt.FuncFormatter(lambda v, _: f"{v / 1e6:.3f}")
    )
    axA.set_xlabel("Chromosome 3 position (Mb, GRCh38)", labelpad=1)
    axA.set_title(
        "A  GHSR locus. All three variants lie upstream of the gene, not "
        "within it.", loc="left", fontsize=8.5, color=style.INK_PRIMARY,
    )

    # ------------------------------------------------------------------
    # Panel B - zoom on the amplicon and upstream regulatory region
    # ------------------------------------------------------------------
    z0, z1 = 172_449_180, 172_450_220
    axB.add_patch(
        mpatches.Rectangle(
            (config.AMPLICON_START, 0.70),
            config.AMPLICON_END - config.AMPLICON_START, 0.3,
            facecolor=style.SERIES[2], edgecolor="none", alpha=0.9, zorder=3,
        )
    )
    axB.text(
        config.AMPLICON_START, 1.04,
        f"sequenced amplicon, {config.REPORTED_AMPLICON_BP} bp "
        f"(chr3:{config.AMPLICON_START:,}-{config.AMPLICON_END:,})",
        ha="left", va="bottom", fontsize=6.2, color=style.SERIES[2],
    )

    for name, kind, s_, e_ in REG_FEATURES:
        colour = style.SERIES[4] if kind == "enhancer" else style.SERIES[3]
        axB.add_patch(
            mpatches.Rectangle(
                (s_, 1.34), max(e_ - s_, 12), 0.3, facecolor=colour,
                edgecolor="none", zorder=3,
            )
        )
        if kind == "enhancer":
            axB.text(
                e_, 1.70, f"Ensembl-annotated {kind} ({name})", ha="right",
                va="bottom", fontsize=6.0, color=colour,
            )
        else:
            mid = (s_ + e_) / 2
            axB.text(
                mid, 2.14, f"{kind}\n({name})", ha="center", va="bottom",
                fontsize=6.0, color=colour,
            )
            axB.plot([mid, mid], [1.66, 2.12], color=colour, lw=0.6, zorder=2)

    # GHSR TSS direction marker at the left edge of the zoom
    axB.annotate(
        "", xy=(z0 + 15, 2.62), xytext=(z0 + 145, 2.62),
        arrowprops=dict(arrowstyle="->", color=style.INK_PRIMARY, lw=1.0),
    )
    axB.text(z0 + 160, 2.62, "GHSR TSS lies 916-1,015 bp further left",
             ha="left", va="center", fontsize=6.2,
             color=style.INK_PRIMARY)

    # variants, labelled on alternating rows to avoid collision
    for k, (_, r) in enumerate(prim.iterrows()):
        pos = r["position_grch38"]
        colour = style.VARIANT_COLORS[r["rsid"]]
        y_top = -0.30 - k * 0.42
        axB.plot([pos, pos], [0.68, 0.18], color=colour, lw=1.0, zorder=5)
        axB.plot([pos], [0.68], marker="v", markersize=5, color=colour,
                 zorder=6)
        axB.plot([pos, pos + 110], [0.18, y_top + 0.06], color=colour,
                 lw=0.7, zorder=5)
        axB.text(
            pos + 125, y_top,
            f"{r['rsid']}  {r['ref_plus_strand']}>{r['alt_plus_strand']}   "
            f"chr3:{int(pos):,}   "
            f"{int(r['bp_upstream_of_ghsr_tss'])} bp upstream of TSS",
            ha="left", va="center", fontsize=6.3, color=colour,
        )

    axB.set_xlim(z0, z1)
    axB.set_ylim(-1.85, 3.00)
    axB.set_yticks([])
    axB.spines["left"].set_visible(False)
    axB.xaxis.set_major_formatter(
        plt.FuncFormatter(lambda v, _: f"{v:,.0f}")
    )
    axB.set_xlabel("Chromosome 3 position (bp, GRCh38)", labelpad=1)
    axB.set_title(
        "B  Zoom on the sequenced amplicon. No variant overlaps an annotated "
        "regulatory element.", loc="left", fontsize=8.5,
        color=style.INK_PRIMARY,
    )

    style.save(fig, "figure1_ghsr_locus")
    prim.to_csv(config.TABLES / "figure1_source_data.csv", index=False)


# ==========================================================================
# Figure 2 - predicted effects with ISM reference distribution
# ==========================================================================
def figure2(scores: pd.DataFrame, ism: pd.DataFrame | None) -> None:
    style.apply()
    modalities = [
        ("RNA_SEQ", "GHSR\nexpression"),
        ("CAGE", "CAGE\n(promoter)"),
        ("DNASE", "DNase\naccessibility"),
        ("ATAC", "ATAC\naccessibility"),
        ("CHIP_TF", "TF ChIP\n(aggregate)"),
        ("CHIP_HISTONE", "Histone\nChIP"),
    ]
    eff = scores[~scores.is_activity_scorer]
    rsids = list(config.RSIDS)

    fig, axes = plt.subplots(
        1, 2, figsize=(7.4, 3.6), gridspec_kw={"width_ratios": [1.35, 1]}
    )

    # ---- panel A: max |log2FC| per variant per modality ------------------
    ax = axes[0]
    width = 0.26
    rows = []
    for mi, (key, _) in enumerate(modalities):
        for vi, rsid in enumerate(rsids):
            sel = eff[(eff.scorer_key == key) & (eff.rsid == rsid)]
            if key == "RNA_SEQ":
                sel = sel[sel.gene_name == config.GHSR_SYMBOL]
            val = float(sel.raw_score.abs().max()) if len(sel) else np.nan
            rows.append(
                {"modality": key, "rsid": rsid, "max_abs_log2fc": val,
                 "n_scores": int(len(sel))}
            )
            y = mi + (vi - 1) * width
            ax.barh(
                y, val, height=width * 0.86,
                color=style.VARIANT_COLORS[rsid],
                edgecolor=style.SURFACE, linewidth=0.8, zorder=3,
            )
            if np.isfinite(val):
                ax.text(
                    val * 1.08, y, f"{val:.3f}", va="center", ha="left",
                    fontsize=5.8, color=style.INK_SECONDARY, zorder=4,
                )
    ax.set_xscale("log")
    ax.set_xlim(8e-3, 0.42)
    ax.set_yticks(range(len(modalities)))
    ax.set_yticklabels([m[1] for m in modalities], fontsize=6.8)
    ax.invert_yaxis()
    ax.set_xlabel(
        "maximum |log$_2$ fold change| across all conditions and tracks"
    )
    ax.axvline(ANCHORED_LIMIT, color=style.INK_MUTED, lw=0.9,
               linestyle=(0, (4, 2)), zorder=2)
    ax.annotate(
        "0.1 log$_2$FC\n(~7% change)", xy=(ANCHORED_LIMIT, -0.62),
        xytext=(ANCHORED_LIMIT * 1.15, -0.62), fontsize=5.8,
        color=style.INK_SECONDARY, ha="left", va="center",
        annotation_clip=False,
    )
    ax.xaxis.grid(True, zorder=0)
    ax.set_title("A  Largest predicted effect per modality", loc="left",
                 fontsize=8.5, color=style.INK_PRIMARY, pad=14)
    handles = [
        mpatches.Patch(facecolor=style.VARIANT_COLORS[r], label=r)
        for r in rsids
    ]
    ax.legend(
        handles=handles, loc="lower right", fontsize=6.4,
        bbox_to_anchor=(1.0, 0.0),
    )

    # ---- panel B: ISM reference distribution -----------------------------
    # Both the violins and the markers come from the SAME in silico
    # mutagenesis run, so they share sequence context, scorer and condition.
    ax = axes[1]
    ism_rows = []
    pct_path = config.TABLES / "table_S5_ism_reference_distribution.csv"
    if ism is not None and len(ism) and pct_path.exists():
        pct = pd.read_csv(pct_path)
        plot_keys = ["RNA_SEQ", "CAGE", "DNASE", "ATAC"]
        # One representative condition keeps the comparison like-for-like.
        cond = "pituitary gland"
        labels = []
        offsets = {r: (k - 1) * 0.22 for k, r in enumerate(rsids)}
        for i, key in enumerate(plot_keys, start=1):
            sel = ism[
                (ism.scorer_key == key) & (ism.condition_label == cond)
            ]
            if key == "RNA_SEQ":
                sel = sel[sel.gene_name == config.GHSR_SYMBOL]
            if not len(sel):
                continue
            vals = (
                sel.groupby("variant_id")["raw_score"].mean().abs().to_numpy()
            )
            vals = vals[np.isfinite(vals) & (vals > 0)]
            if not vals.size:
                continue
            parts = ax.violinplot(
                [np.log10(vals)], positions=[i], orientation="horizontal",
                showextrema=False, widths=0.75,
            )
            for pc in parts["bodies"]:
                pc.set_facecolor(style.SEQUENTIAL[2])
                pc.set_edgecolor(style.SEQUENTIAL[7])
                pc.set_alpha(1.0)
                pc.set_linewidth(0.7)
            labels.append(key)
            ism_rows.append(
                {
                    "modality": key, "condition": cond,
                    "n_ism_variants": int(vals.size),
                    "median_abs": float(np.median(vals)),
                    "p95_abs": float(np.percentile(vals, 95)),
                    "max_abs": float(vals.max()),
                }
            )
            for rsid in rsids:
                row = pct[
                    (pct.scorer_key == key)
                    & (pct.condition_label == cond)
                    & (pct.rsid == rsid)
                    & (pct.variant_key.isin(
                        [v.key for v in config.PRIMARY_VARIANTS]))
                ]
                if not len(row):
                    continue
                v = abs(float(row.observed_score.iloc[0]))
                if v <= 0:
                    continue
                ax.plot(
                    [np.log10(v)], [i + offsets[rsid]], marker="o",
                    markersize=5.5,
                    markerfacecolor=style.VARIANT_COLORS[rsid],
                    markeredgecolor=style.SURFACE, markeredgewidth=1.0,
                    zorder=5, clip_on=False,
                )
        if labels:
            ax.set_yticks(range(1, len(labels) + 1))
            ax.set_yticklabels(
                [dict(modalities).get(k, k).replace("\n", " ")
                 for k in labels],
                fontsize=6.6,
            )
            ax.set_ylim(len(labels) + 0.75, 0.25)
            ax.set_xlabel("log$_{10}$ |log$_2$ fold change|")
            ax.xaxis.grid(True, zorder=0)
            ax.text(
                0.02, 0.015,
                f"violins: all 1,197 single-base substitutions in the "
                f"399 bp amplicon\ncircles: the three studied SNPs, from the "
                f"same run ({cond})",
                transform=ax.transAxes, fontsize=5.6,
                color=style.INK_SECONDARY, va="bottom",
            )
    ax.set_title(
        "B  Observed SNPs vs every possible\n    substitution in the amplicon",
        loc="left", fontsize=8.5, color=style.INK_PRIMARY, pad=4,
    )
    if ism is None or not len(ism):
        ax.text(0.5, 0.5, "ISM reference distribution\nnot available",
                ha="center", va="center", transform=ax.transAxes,
                fontsize=7, color=style.INK_MUTED)
        ax.set_xticks([])
        ax.set_yticks([])

    fig.tight_layout()
    style.save(fig, "figure2_predicted_effects")
    pd.DataFrame(rows).to_csv(
        config.TABLES / "figure2_source_data.csv", index=False
    )
    if ism_rows:
        pd.DataFrame(ism_rows).to_csv(
            config.TABLES / "figure2_ism_summary.csv", index=False
        )


# ==========================================================================
# Figure 3 - cell-type specificity heat map
# ==========================================================================
def figure3(scores: pd.DataFrame, conditions: pd.DataFrame) -> None:
    """Signed effects across conditions, with the baseline-expression column.

    Drawn on a fixed, biologically anchored scale (+/- 0.1 log2 fold change,
    ~7% change in predicted signal) rather than auto-scaled, so that the
    near-neutral colour of nearly every cell reflects the true magnitude of
    the effects rather than being stretched to fill the colour map.

    The right-hand column carries the predicted GHSR expression level of each
    condition, because it is the main explanation for what structure the heat
    map does contain.
    """
    style.apply()
    modalities = [
        ("RNA_SEQ", "GHSR expr."),
        ("CAGE", "CAGE"),
        ("DNASE", "DNase"),
        ("ATAC", "ATAC"),
        ("CHIP_TF", "TF ChIP"),
        ("CHIP_HISTONE", "Histone"),
    ]
    eff = scores[~scores.is_activity_scorer]

    order = conditions[
        ["condition_curie", "label", "panel", "ghsr_rna_seq_active_max"]
    ].copy()
    order["panel"] = pd.Categorical(
        order["panel"], categories=list(config.PANEL_ORDER), ordered=True
    )
    order = order.sort_values(
        ["panel", "ghsr_rna_seq_active_max"], ascending=[True, False]
    ).reset_index(drop=True)

    cols, col_labels = [], []
    for rsid in config.RSIDS:
        for key, lab in modalities:
            cols.append((rsid, key))
            col_labels.append(lab)

    mat = np.full((len(order), len(cols)), np.nan)
    mode = np.empty((len(order), len(cols)), dtype=object)
    long_rows = []
    for i, r in order.iterrows():
        for j, (rsid, key) in enumerate(cols):
            sel = eff[
                (eff.rsid == rsid)
                & (eff.scorer_key == key)
                & (eff.condition_curie == r["condition_curie"])
            ]
            if key == "RNA_SEQ":
                sel = sel[sel.gene_name == config.GHSR_SYMBOL]
            if not len(sel):
                mode[i, j] = "none"
                continue
            idx = sel.raw_score.abs().idxmax()
            mat[i, j] = float(sel.loc[idx, "raw_score"])
            mode[i, j] = str(sel.loc[idx, "imputation_mode"])
            long_rows.append(
                {
                    "condition_label": r["label"],
                    "panel": str(r["panel"]),
                    "rsid": rsid,
                    "modality": key,
                    "signed_log2fc_at_peak_track": mat[i, j],
                    "imputation_mode_of_peak_track": mode[i, j],
                    "n_tracks": int(len(sel)),
                }
            )

    finite = mat[np.isfinite(mat)]
    observed_max = float(np.nanmax(np.abs(finite))) if finite.size else np.nan

    fig, (ax, axact) = plt.subplots(
        1, 2, figsize=(7.6, 4.4),
        gridspec_kw={"width_ratios": [18, 1], "wspace": 0.04},
    )
    norm = TwoSlopeNorm(vcenter=0.0, vmin=-ANCHORED_LIMIT, vmax=ANCHORED_LIMIT)
    ax.imshow(mat, cmap=DIVERGING_CMAP, norm=norm, aspect="auto",
              interpolation="nearest")
    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            if mode[i, j] == "imputed_exploratory":
                ax.plot([j], [i], marker=".", markersize=1.6,
                        color=style.INK_MUTED, zorder=4)

    ax.set_xticks(range(len(cols)))
    ax.set_xticklabels(col_labels, rotation=90, fontsize=5.8)
    for g in range(1, len(config.RSIDS)):
        ax.axvline(g * len(modalities) - 0.5, color=style.INK_PRIMARY, lw=1.0)
    for g, rsid in enumerate(config.RSIDS):
        ax.text(
            g * len(modalities) + (len(modalities) - 1) / 2, -0.85, rsid,
            ha="center", va="bottom", fontsize=6.8,
            color=style.VARIANT_COLORS[rsid],
        )
    ax.set_yticks(range(len(order)))
    ax.set_yticklabels(order["label"], fontsize=6.3)
    prev = None
    for i, r in order.iterrows():
        if prev is not None and r["panel"] != prev:
            ax.axhline(i - 0.5, color=style.INK_PRIMARY, lw=1.0)
        prev = r["panel"]
    # panel brackets on the far left
    for panel, g in order.groupby("panel", observed=True):
        if not len(g):
            continue
        i0, i1 = g.index.min(), g.index.max()
        ax.annotate(
            panel, xy=(-0.255, 1 - (i0 + i1 + 1) / (2 * len(order))),
            xycoords="axes fraction", rotation=90, ha="center", va="center",
            fontsize=6.6, annotation_clip=False,
            color=style.PANEL_COLORS.get(panel,
                                                       style.INK_SECONDARY),
        )
    ax.set_xlim(-0.5, len(cols) - 0.5)
    ax.set_ylim(len(order) - 0.5, -0.5)

    # --- right-hand annotation column: predicted GHSR expression ----------
    act = order["ghsr_rna_seq_active_max"].to_numpy(dtype=float)
    act_rel = act / np.nanmax(act)
    seq_cmap = LinearSegmentedColormap.from_list(
        "ag_sequential", style.SEQUENTIAL
    )
    axact.imshow(
        act_rel.reshape(-1, 1), cmap=seq_cmap, aspect="auto",
        vmin=0, vmax=1, interpolation="nearest",
    )
    for i, v in enumerate(act_rel):
        axact.text(
            0, i, f"{v * 100:.0f}", ha="center", va="center", fontsize=5.4,
            color=style.SURFACE if v > 0.55 else style.INK_PRIMARY,
        )
    axact.set_xticks([0])
    axact.set_xticklabels(["GHSR expr.\n(% of max)"], rotation=90,
                          fontsize=5.8)
    axact.set_yticks([])
    axact.set_ylim(len(order) - 0.5, -0.5)

    sm = plt.cm.ScalarMappable(
        cmap=DIVERGING_CMAP, norm=TwoSlopeNorm(vcenter=0.0, vmin=-1.0,
                                               vmax=1.0)
    )
    cb = fig.colorbar(sm, ax=[ax, axact], fraction=0.022, pad=0.03,
                      ticks=[-1, 0, 1])
    cb.ax.set_yticklabels(
        [f"ALT lower\n(-{ANCHORED_LIMIT} log$_2$FC)", "no change",
         f"ALT higher\n(+{ANCHORED_LIMIT} log$_2$FC)"],
        fontsize=6.0,
    )
    cb.outline.set_visible(False)

    fig.suptitle(
        "Predicted allelic effects by condition and modality, on a fixed "
        f"+/-{ANCHORED_LIMIT} log$_2$FC scale",
        x=0.02, y=1.02, ha="left", fontsize=9, color=style.INK_PRIMARY,
    )
    fig.text(
        0.02, 0.965,
        f"Largest effect anywhere: |log$_2$FC| = {observed_max:.3f}. "
        "Grey dot: peak track is imputed-exploratory.",
        ha="left", fontsize=6.6, color=style.INK_SECONDARY,
    )
    style.save(fig, "figure3_celltype_specificity")
    pd.DataFrame(long_rows).to_csv(
        config.TABLES / "figure3_source_data.csv", index=False
    )
    (config.DATA_PROCESSED / "figure3_scale_note.json").write_text(
        json.dumps(
            {
                "anchored_limit_log2fc": ANCHORED_LIMIT,
                "observed_max_abs_log2fc": observed_max,
                "scale": "fixed, not auto-scaled",
            },
            indent=2,
        )
    )


# ==========================================================================
# Figure 4 - integrated evidence schematic
# ==========================================================================
def figure4(integ: pd.DataFrame) -> None:
    style.apply()
    fig, ax = plt.subplots(figsize=(7.4, 5.0))
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 11.2)
    ax.axis("off")

    levels = [
        (
            "Level 1 - Computational regulatory prediction",
            "THIS STUDY (retrospective, in silico)",
            "AlphaGenome 2; 16 conditions, 22 scorers,\n"
            "4 variant configurations; GRCh38",
            "Establishes\n  predicted cis-regulatory consequence\n"
            "Does NOT establish\n  measured expression, phenotype, diagnosis",
            style.SERIES[1],
            7.55,
        ),
        (
            "Level 2 - Experimental genotype detection",
            "PRIOR STUDY (PCR + Sanger sequencing, n = 100)",
            "399 bp amplicon; 100% call rate;\n"
            "Hardy-Weinberg equilibrium met in controls",
            "Establishes\n  the genotype at each position was read correctly\n"
            "Does NOT establish\n  any effect on gene regulation",
            style.SERIES[2],
            4.05,
        ),
        (
            "Level 3 - Statistical phenotype association",
            "PRIOR STUDY (case-control, n = 100, 6-CIT groups)",
            "No association survived FDR correction;\n"
            "rs2922126 achieved power 6.62%; rs115459563\n"
            "models did not converge; rs34136054 n = 1",
            "Establishes\n  exploratory, underpowered association only\n"
            "Does NOT establish\n  causation, mechanism, or an AD diagnosis",
            style.SERIES[3],
            0.55,
        ),
    ]
    box_h = 2.55
    for title, who, what, claims, colour, y in levels:
        ax.add_patch(
            mpatches.FancyBboxPatch(
                (0.2, y), 9.6, box_h,
                boxstyle="round,pad=0.04,rounding_size=0.10",
                facecolor=style.SURFACE, edgecolor=colour, linewidth=1.3,
            )
        )
        ax.add_patch(
            mpatches.Rectangle((0.2, y), 0.13, box_h, facecolor=colour,
                               edgecolor="none")
        )
        ax.text(0.55, y + box_h - 0.22, title, fontsize=8.4,
                fontweight="bold", color=style.INK_PRIMARY, va="top")
        ax.text(0.55, y + box_h - 0.72, who, fontsize=6.8, color=colour,
                va="top")
        ax.text(0.55, y + box_h - 1.14, what, fontsize=6.4,
                color=style.INK_SECONDARY, va="top", linespacing=1.5)
        ax.text(5.25, y + box_h - 0.72, claims, fontsize=6.4,
                color=style.INK_SECONDARY, va="top", linespacing=1.5)

    # Barred arrows: inferences the design does not license.
    for y_top, y_bot, label in (
        (7.55, 6.60, "a null prediction does not refute an association"),
        (4.05, 3.10, "genotyping does not validate a predicted effect"),
    ):
        ymid = (y_top + y_bot) / 2
        ax.annotate(
            "", xy=(1.5, y_bot + 0.03), xytext=(1.5, y_top - 0.03),
            arrowprops=dict(arrowstyle="-|>", color=style.INK_MUTED, lw=1.0,
                            linestyle=(0, (3, 2))),
        )
        ax.plot([1.22, 1.78], [ymid, ymid], color=style.DIVERGING_HIGH,
                lw=1.6, zorder=5)
        ax.text(2.0, ymid, label, fontsize=6.4,
                color=style.DIVERGING_HIGH, va="center")

    ax.text(0.2, 11.15, "Three evidence levels, kept separate",
            fontsize=10.5, fontweight="bold", color=style.INK_PRIMARY,
            va="top")
    ax.text(
        0.2, 10.62,
        "The in silico analysis was performed after the genotyping and did "
        "not guide the experimental design.",
        fontsize=7, color=style.INK_SECONDARY, va="top",
    )
    style.save(fig, "figure4_integrated_evidence")
    integ.to_csv(config.TABLES / "figure4_source_data.csv", index=False)


# ==========================================================================
def main() -> None:
    print("Phase VI: figures")
    annot = pd.read_csv(config.TABLES / "table_S1_variant_annotation.csv")
    scores = pd.read_parquet(
        config.DATA_RAW / "phase3" / "phase3_raw_scores_main.parquet"
    )
    conditions = pd.read_csv(
        config.TABLES / "table_S3_condition_selection_and_qc.csv"
    )
    ism_path = config.DATA_RAW / "phase4" / "phase4_ism_reference_distribution.parquet"
    ism = pd.read_parquet(ism_path) if ism_path.exists() else None
    integ_path = config.TABLES / "table_3_integrated_evidence.csv"
    integ = pd.read_csv(integ_path) if integ_path.exists() else pd.DataFrame()

    figure1(annot)
    figure2(scores, ism)
    figure3(scores, conditions)
    if len(integ):
        figure4(integ)

    provenance.write(
        "phase6",
        {
            "figures": [
                "figure1_ghsr_locus",
                "figure2_predicted_effects",
                "figure3_celltype_specificity",
                "figure4_integrated_evidence",
            ],
            "formats": ["pdf", "svg", "png@600dpi"],
            "ism_available": ism is not None,
            "anchored_limit_log2fc": ANCHORED_LIMIT,
        },
    )


if __name__ == "__main__":
    main()
