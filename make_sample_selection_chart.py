from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "manuscript" / "AJSE_submission"

plt.rcParams.update(
    {
        "font.family": "serif",
        "font.serif": ["STIX Two Text", "STIXGeneral", "Times New Roman", "DejaVu Serif"],
        "mathtext.fontset": "stix",
        "font.size": 9,
    }
)

INK = "#222222"
TEAL = "#3F7480"
LIGHT = "#F4F4F2"
PALE_TEAL = "#E7EFF0"
PALE_WARM = "#F4ECE8"
MUTED = "#666666"


def box(ax, xy, wh, title, lines, face=LIGHT, edge=INK, lw=0.9, title_color=INK):
    x, y = xy
    w, h = wh
    patch = FancyBboxPatch(
        (x, y), w, h,
        boxstyle="round,pad=0.012,rounding_size=0.012",
        linewidth=lw, edgecolor=edge, facecolor=face,
    )
    ax.add_patch(patch)
    ax.text(x + 0.018, y + h - 0.035, title, ha="left", va="top",
            fontsize=9.4, fontweight="bold", color=title_color)
    ax.text(x + 0.018, y + h - 0.090, lines, ha="left", va="top",
            fontsize=8.2, color=INK, linespacing=1.32)
    return patch


def arrow(ax, start, end, color=MUTED, lw=1.0):
    ax.add_patch(FancyArrowPatch(start, end, arrowstyle="-|>", mutation_scale=9,
                                 linewidth=lw, color=color, shrinkA=2, shrinkB=2,
                                 connectionstyle="arc3,rad=0"))


fig, ax = plt.subplots(figsize=(12.0, 6.6))
ax.set_xlim(0, 1)
ax.set_ylim(0, 1)
ax.axis("off")

# Primary ACI sample and the three analysis-specific branches.
box(ax, (0.035, 0.735), (0.20, 0.16), "ACI 133 source survey", "242 RC buildings\n234 unique coordinates", face="#FFFFFF")

box(ax, (0.315, 0.765), (0.235, 0.17), "Primary model assessment", "Baseline LFINF + baseline design\nExact GEM height + two mainshocks\n$n=109$ buildings", face=PALE_TEAL, edge=TEAL, lw=1.35, title_color=TEAL)
box(ax, (0.315, 0.515), (0.235, 0.17), "Event-scope sensitivity", "Primary rules + common coverage\nTwo versus four events\n$n=60$ buildings", face="#FFFFFF")
box(ax, (0.315, 0.265), (0.235, 0.19), "Taxonomy sensitivity", "3 systems × 3 designs × 2 height policies\n18 scenarios; nearest-height intersection\n$n=91$ buildings", face="#FFFFFF")

arrow(ax, (0.235, 0.815), (0.315, 0.850), color=TEAL, lw=1.25)
ax.plot([0.432, 0.432], [0.765, 0.715], color=MUTED, lw=0.9)
arrow(ax, (0.432, 0.715), (0.432, 0.685))
arrow(ax, (0.235, 0.795), (0.315, 0.365))

box(ax, (0.635, 0.785), (0.30, 0.13), "Primary outputs", "Loss and fragility performance; residual analysis", face=PALE_TEAL, edge=TEAL, lw=1.1, title_color=TEAL)
box(ax, (0.635, 0.555), (0.30, 0.13), "Coverage-controlled comparison", "Same buildings isolate event scope from grid coverage", face="#FFFFFF")
box(ax, (0.635, 0.315), (0.30, 0.15), "Fixed-sample comparison", "Exact-height common intersection: $n=3$ (insufficient)\nMatched exact-height checks: $n=46$ and $n=105$", face="#FFFFFF")

arrow(ax, (0.550, 0.850), (0.635, 0.850), color=TEAL, lw=1.25)
arrow(ax, (0.550, 0.600), (0.635, 0.620))
arrow(ax, (0.550, 0.360), (0.635, 0.390))

# Clarify the subset used only for binned calibration graphics.
arrow(ax, (0.550, 0.790), (0.635, 0.735))
ax.text(0.645, 0.735, "IM-bin diagnostic: 85 of 109 use SA(0.6 s);\nfive groups of 17 buildings",
        ha="left", va="center", fontsize=8.1, color=MUTED)

# Separate replication dataset: this is not a continuation of the ACI filtering path.
ax.plot([0.035, 0.935], [0.205, 0.205], color="#B8B8B8", lw=0.8)
ax.text(0.035, 0.180, "Independent aggregate replication", ha="left", va="top",
        fontsize=9.2, fontweight="bold", color=INK)
box(ax, (0.035, 0.025), (0.20, 0.105), "Gaziantep source", "1,308 records", face=PALE_WARM)
box(ax, (0.315, 0.025), (0.235, 0.105), "Eligible RC records", "739 RC; 676 complete and post-1974", face="#FFFFFF")
box(ax, (0.635, 0.025), (0.30, 0.105), "Replication analysis", "Same baseline exact-height rules: $n=647$", face=PALE_WARM)
arrow(ax, (0.235, 0.077), (0.315, 0.077))
arrow(ax, (0.550, 0.077), (0.635, 0.077))

fig.savefig(OUT / "Fig_sample_selection.png", dpi=450, bbox_inches="tight", facecolor="white")
fig.savefig(OUT / "Fig_sample_selection.pdf", bbox_inches="tight", facecolor="white")
plt.close(fig)
