"""Generate title-free publication figures for AJSE; captions remain in the manuscript."""
from pathlib import Path
import shutil

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import pandas as pd

ROOT = Path(__file__).resolve().parent
DATA = ROOT / "outputs/aci133_exploratory"
DEST = ROOT / "outputs/manuscript_figures"
BASE = "lfinf__baseline__exact"
MAPPING = "code_C_as_DS4"

plt.rcParams.update({"font.family": "sans-serif", "font.size": 9,
                     "pdf.fonttype": 42, "ps.fonttype": 42})


def save(fig, stem):
    DEST.mkdir(parents=True, exist_ok=True)
    fig.savefig(DEST / f"{stem}.pdf", bbox_inches="tight")
    fig.savefig(DEST / f"{stem}.png", dpi=600, bbox_inches="tight")
    number = stem.split("_", 1)[0]
    shutil.copy2(DEST / f"{stem}.png", ROOT / "manuscript" / "AJSE_submission" / f"{number}.png")
    plt.close(fig)


def main():
    buildings = pd.read_csv(ROOT / "data/processed/aci133_buildings.csv")
    counts = buildings.structural_damage_code.value_counts().reindex(["N", "L", "M", "S", "C"], fill_value=0)
    fig, ax = plt.subplots(figsize=(3.35, 2.5))
    counts.plot.bar(ax=ax, color="#3c6e71", rot=0)
    ax.set(xlabel="ACI structural-damage code", ylabel="Buildings")
    fig.tight_layout(); save(fig, "Fig2_damage_distribution")

    calibration = pd.read_csv(DATA / "fragility_calibration.csv")
    c = calibration[(calibration.scenario == BASE) &
                    (calibration.event_scope == "four_survey_preceding_events") &
                    (calibration.damage_mapping == MAPPING) &
                    (calibration.method == "fragility_max")]
    fig, ax = plt.subplots(figsize=(3.45, 3.15))
    ax.plot([0, 1], [0, 1], color="#8B949E", lw=1.1, ls=(0, (4, 3)), zorder=1)
    ax.scatter(c.predicted_rate, c.observed_rate, s=38, color="#28536B",
               edgecolor="white", linewidth=0.6, zorder=3, clip_on=False)
    offsets = {"DS1+": (-31, 6), "DS2+": (6, 7), "DS3+": (7, 6), "DS4+": (7, 8)}
    for row in c.itertuples(index=False):
        ax.annotate(row.threshold, (row.predicted_rate, row.observed_rate),
                    xytext=offsets.get(row.threshold, (5, 5)), textcoords="offset points",
                    fontsize=8, annotation_clip=False)
    ax.set(xlim=(-0.035, 1.035), ylim=(-0.035, 1.035),
           xlabel="Predicted exceedance probability", ylabel="Observed frequency")
    ax.set_xticks([0, .2, .4, .6, .8, 1.0])
    ax.set_yticks([0, .2, .4, .6, .8, 1.0])
    ax.grid(color="#D9E1E8", linewidth=0.55, alpha=0.8)
    ax.set_axisbelow(True)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout(); save(fig, "Fig3_fragility_calibration")

    bins = pd.read_csv(DATA / "im_bin_calibration.csv")
    b = bins[(bins.scenario == BASE) &
             (bins.event_scope == "four_survey_preceding_events") &
             (bins.damage_mapping == MAPPING) & bins.imt.eq("SA(0.6)")]
    fig, ax = plt.subplots(figsize=(4.2, 3.0))
    ax.plot(b.mean_im_g, b.predicted_loss_mean, color="#35607D", lw=1.8,
            marker="o", ms=5.2, label="GEM prediction")
    ax.plot(b.mean_im_g, b.observed_loss_proxy_mean, color="#A65A6A", lw=1.8,
            ls=(0, (5, 2.5)), marker="D", ms=4.8, label="Damage-derived proxy")
    ax.set(xlabel="Mean SA(0.6 s) within bin (g)", ylabel="Mean loss ratio")
    ax.grid(axis="y", color="#D9E1E8", linewidth=0.55)
    ax.set_axisbelow(True)
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(frameon=False, fontsize=8, loc="best")
    fig.tight_layout(); save(fig, "Fig4_im_bin_calibration")

    metrics = pd.read_csv(DATA / "metrics.csv")
    m = metrics[(metrics.method == "vulnerability_max") &
                (metrics.event_scope == "four_survey_preceding_events") &
                (metrics.damage_mapping == MAPPING)].sort_values("bias_pred_obs").copy()
    system_names = {"lfinf": "Infilled frame", "lfm": "Moment frame", "lwal": "Wall system"}
    system_colors = {"lfinf": "#35607D", "lfm": "#B07A3E", "lwal": "#4F7D68"}
    system_markers = {"exact": "o", "nearest": "s"}
    parsed = m.scenario.str.split("__", expand=True)
    m["system"] = parsed[0]
    m["design"] = parsed[1]
    m["height"] = parsed[2]
    m["display"] = [f"{system_names[s]}  |  {d.title()}  |  {h.title()}"
                    for s, d, h in zip(m.system, m.design, m.height)]

    fig, ax = plt.subplots(figsize=(7.0, 5.1))
    y = range(len(m))
    ax.hlines(y, 0, m.bias_pred_obs, color="#D7DEE5", lw=1.4, zorder=1)
    for yi, row in enumerate(m.itertuples(index=False)):
        ax.scatter(row.bias_pred_obs, yi, s=48, marker=system_markers[row.height],
                   color=system_colors[row.system], edgecolor="white", linewidth=0.6, zorder=3)
        ax.text(row.bias_pred_obs + 0.025, yi, f"n={int(row.n)}", va="center",
                fontsize=7.2, color="#37474F")
    ax.axvline(1, color="#4F5963", lw=1.1, ls=(0, (4, 3)))
    ax.set_yticks(list(y), m.display)
    for tick, scenario in zip(ax.get_yticklabels(), m.scenario):
        if scenario == BASE:
            tick.set_fontweight("bold")
    ax.set_xlim(0, 1.48)
    ax.set(xlabel="Predicted / observed damage-derived loss")
    ax.grid(axis="x", color="#E1E6EA", linewidth=0.55)
    ax.set_axisbelow(True)
    ax.spines[["top", "right", "left"]].set_visible(False)
    ax.tick_params(axis="y", length=0, labelsize=8)
    legend = [Line2D([0], [0], marker="o", color="none", markerfacecolor=system_colors[key],
                     markeredgecolor="white", markersize=6.5, label=name)
              for key, name in system_names.items()]
    ax.legend(handles=legend, frameon=False, ncol=3, fontsize=8,
              loc="lower right", title="Structural-system assumption", title_fontsize=8)
    fig.tight_layout(); save(fig, "Fig5_taxonomy_sensitivity")
    print(f"Wrote publication figures to {DEST}")


if __name__ == "__main__":
    main()
