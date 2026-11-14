"""Shared figure styling.

Colours come from a validated categorical palette; the subset used here was
checked with the data-visualisation palette validator on the light chart
surface (`#fcfcfb`):

    lightness band PASS, chroma floor PASS,
    worst adjacent CVD dE 9.1 (protan) PASS,
    worst adjacent normal-vision dE 22.9 PASS,
    contrast WARN for aqua/yellow -> mitigated by direct value labels and by
    shipping every figure's underlying numbers as a CSV table.

Figures are rendered for print (light surface only); each figure is written as
vector PDF and SVG plus a 600 dpi PNG.
"""

from __future__ import annotations

import matplotlib as mpl
import matplotlib.pyplot as plt

from . import config

# --- categorical slots (fixed order, never cycled) ------------------------
SERIES = {
    1: "#2a78d6",  # blue
    2: "#eb6834",  # orange
    3: "#1baf7a",  # aqua
    4: "#eda100",  # yellow
}

# --- diverging pair for signed effects ------------------------------------
DIVERGING_LOW = "#2a78d6"   # blue  = ALT decreases signal
DIVERGING_MID = "#f0efec"   # neutral gray midpoint
DIVERGING_HIGH = "#e34948"  # red   = ALT increases signal

# --- sequential blue ramp -------------------------------------------------
SEQUENTIAL = [
    "#cde2fb", "#b7d3f6", "#9ec5f4", "#86b6ef", "#6da7ec",
    "#5598e7", "#3987e5", "#2a78d6", "#256abf", "#1c5cab",
]

# --- chrome and ink -------------------------------------------------------
SURFACE = "#fcfcfb"
INK_PRIMARY = "#0b0b0b"
INK_SECONDARY = "#52514e"
INK_MUTED = "#898781"
GRIDLINE = "#e1e0d9"
BASELINE = "#c3c2b7"

# --- panel colours (one categorical slot per biological panel) ------------
PANEL_COLORS = {
    "Brain / neuronal": SERIES[1],
    "Pituitary": SERIES[2],
    "Bone": SERIES[3],
    "Context tissue": INK_MUTED,
}

# --- variant colours (one slot per rsID, assigned by entity not rank) -----
VARIANT_COLORS = {
    "rs2922126": SERIES[1],
    "rs115459563": SERIES[2],
    "rs34136054": SERIES[3],
}

# --- evidence-class hatching (the accessibility channel) ------------------
IMPUTATION_HATCH = {
    "observed": "",
    "imputed_supported": "///",
    "imputed_exploratory": "xxx",
}


def apply() -> None:
    mpl.rcParams.update(
        {
            "figure.facecolor": SURFACE,
            "axes.facecolor": SURFACE,
            "savefig.facecolor": SURFACE,
            "font.family": "sans-serif",
            "font.sans-serif": ["DejaVu Sans"],
            "font.size": 8,
            "axes.titlesize": 9,
            "axes.labelsize": 8,
            "axes.labelcolor": INK_SECONDARY,
            "axes.edgecolor": BASELINE,
            "axes.linewidth": 0.8,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.grid": False,
            "grid.color": GRIDLINE,
            "grid.linewidth": 0.6,
            "xtick.color": INK_MUTED,
            "ytick.color": INK_MUTED,
            "xtick.labelcolor": INK_SECONDARY,
            "ytick.labelcolor": INK_SECONDARY,
            "xtick.major.width": 0.8,
            "ytick.major.width": 0.8,
            "legend.frameon": False,
            "legend.fontsize": 7,
            "lines.linewidth": 1.6,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
        }
    )


def save(fig: plt.Figure, name: str) -> list[str]:
    """Write a figure as PDF, SVG and 600 dpi PNG; return the paths."""
    paths = []
    for ext, kw in (
        ("pdf", {}),
        ("svg", {}),
        ("png", {"dpi": 600}),
    ):
        p = config.FIGURES / f"{name}.{ext}"
        fig.savefig(p, bbox_inches="tight", **kw)
        paths.append(str(p))
    plt.close(fig)
    print(f"[figure] {name}: " + ", ".join(f".{p.rsplit('.', 1)[1]}" for p in paths))
    return paths
