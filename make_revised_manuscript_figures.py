"""Generate publication figures for the revised AJSE manuscript."""
from pathlib import Path
import json
import shutil

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import transforms
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
RAW = ROOT / "data/raw"
DATA = ROOT / "outputs/manuscript_enhanced"
DEST = ROOT / "outputs/manuscript_figures"
SUBMISSION = ROOT / "manuscript/AJSE_submission"
TURKEY_ADM1 = RAW / "boundaries/geoBoundaries-TUR-ADM1_simplified.geojson"
SYRIA_ADM0 = RAW / "boundaries/geoBoundaries-SYR-ADM0_simplified.geojson"

# Analytical figures use the same restrained palette as the study-area map.
# Teal consistently denotes GEM output and terracotta denotes observations.
PRIMARY = "#3F7480"
OBSERVED = "#B85F47"
PRIMARY_FILL = "#688F97"
OBSERVED_FILL = "#C97A61"
NEUTRAL_FILL = "#969696"
HIGHLIGHT = OBSERVED
INK = "#263746"
INK_LIGHT = "#7A8387"
WARM_INK = OBSERVED
COLORS = {"N": "#D7DEE2", "L": "#8DAF9F", "M": "#D1B363",
          "S": "#C66D52", "C": "#713247"}
SYSTEM_COLORS = {"lfinf": PRIMARY_FILL, "lfm": OBSERVED_FILL, "lwal": NEUTRAL_FILL}

plt.rcParams.update({"font.family": "serif", "font.serif": ["STIXGeneral"],
                     "mathtext.fontset": "stix", "font.size": 8.5,
                     "axes.labelcolor": "#263746", "text.color": "#263746",
                     "hatch.linewidth": .65,
                     "pdf.fonttype": 42, "ps.fonttype": 42})


def save(fig, stem):
    DEST.mkdir(parents=True, exist_ok=True)
    fig.savefig(DEST / f"{stem}.pdf", bbox_inches="tight")
    fig.savefig(DEST / f"{stem}.png", dpi=600, bbox_inches="tight")
    number = stem.split("_", 1)[0]
    shutil.copy2(DEST / f"{stem}.png", SUBMISSION / f"{number}.png")
    plt.close(fig)


def journal_axes(ax, keep_grid=False):
    """Apply the adopted AJSE-style frame to analytical panels."""
    if not keep_grid:
        ax.grid(False)
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_color("black")
        spine.set_linewidth(.75)
    ax.tick_params(direction="in", width=.7, length=3)


def iter_lines(geometry):
    kind, coords = geometry.get("type"), geometry.get("coordinates", [])
    if kind == "LineString":
        yield coords
    elif kind == "MultiLineString":
        yield from coords
    elif kind == "Polygon":
        yield from coords
    elif kind == "MultiPolygon":
        for polygon in coords:
            yield from polygon
    elif kind == "GeometryCollection":
        for item in geometry.get("geometries", []):
            yield from iter_lines(item)


def plot_study_map():
    b = pd.read_csv(ROOT / "data/processed/aci133_buildings.csv")
    # Give the spatial information a full-width panel.  The damage-count bars
    # are deliberately exported as a separate, single-column figure below.
    fig, ax = plt.subplots(figsize=(7.2, 4.65))
    # Distinguish water from neighbouring land explicitly.  The axes
    # background represents the Mediterranean, while the Syria ADM0 polygon
    # is drawn as contextual land before the Türkiye province layer.
    ax.set_facecolor("#EAF3F6")
    if SYRIA_ADM0.exists():
        syria = json.loads(SYRIA_ADM0.read_text(encoding="utf-8"))
        for feature in syria.get("features", []):
            for line in iter_lines(feature["geometry"]):
                arr = np.asarray(line)
                if arr.ndim == 2 and arr.shape[1] >= 2:
                    ax.fill(arr[:, 0], arr[:, 1], facecolor="#F7F6F2",
                            edgecolor="#AEB8BD", lw=.42, zorder=-1)
    # Use one consistent high-resolution source for both the land mask and
    # province boundaries. The pinned file represents Turkey's 81 ADM1 units
    # (2021), avoiding mismatched outlines from a coarse country basemap.
    if TURKEY_ADM1.exists():
        provinces = json.loads(TURKEY_ADM1.read_text(encoding="utf-8"))
        for feature in provinces.get("features", []):
            for line in iter_lines(feature["geometry"]):
                arr = np.asarray(line)
                if arr.ndim == 2 and arr.shape[1] >= 2:
                    ax.fill(arr[:, 0], arr[:, 1], facecolor="#F4F1EA",
                            edgecolor="#BCC5CA", lw=.34, alpha=.98, zorder=0)

    for key, color in [("event1", "#A84642"), ("event2", "#684E70")]:
        rupture = json.loads((RAW / f"{key}_rupture.json").read_text(encoding="utf-8"))
        for feature in rupture.get("features", []):
            for line in iter_lines(feature["geometry"]):
                arr = np.asarray(line)
                if arr.ndim == 2 and arr.shape[1] >= 2:
                    ax.plot(arr[:, 0], arr[:, 1], color=color, lw=1.35, alpha=.9, zorder=2)

    for key, marker, color, label in [
        ("event1", "*", "#A84642", "6 Feb Mw 7.8"),
        ("event2", "*", "#684E70", "6 Feb Mw 7.5"),
        ("event3", "^", "#C4863E", "6 Feb Mw 6.7"),
        ("event4", "^", "#3F7480", "20 Feb Mw 6.3")]:
        event = json.loads((RAW / f"{key}.json").read_text(encoding="utf-8"))
        lon, lat = event["geometry"]["coordinates"][:2]
        ax.scatter(lon, lat, marker=marker, s=85 if marker == "*" else 38,
                   color=color, ec="white", lw=.6, zorder=5, label=label)

    info = json.loads((RAW / "event4_info.json").read_text(encoding="utf-8"))
    mi = info["output"]["map_information"]
    xmin, ymin = float(mi["min"]["longitude"]), float(mi["min"]["latitude"])
    xmax, ymax = float(mi["max"]["longitude"]), float(mi["max"]["latitude"])
    ax.add_patch(Rectangle((xmin, ymin), xmax-xmin, ymax-ymin, fill=False,
                           ec="#3F7480", lw=1, ls=(0, (4, 3)), zorder=1,
                           label="20 Feb grid extent"))

    # Plot unique sites directly. Records sharing an exact coordinate are
    # separated by fixed display-space offsets so that all 13 records in the
    # five coincident-coordinate groups remain visible without random jitter.
    duplicated = b.duplicated(["longitude", "latitude"], keep=False)
    for code in ["N", "L", "M", "S", "C"]:
        z = b[(b.structural_damage_code == code) & ~duplicated]
        ax.scatter(z.longitude, z.latitude, s=20, color=COLORS[code],
                   ec="#263746", lw=.25, alpha=.9, zorder=4)
    for _, group in b[duplicated].groupby(["longitude", "latitude"], sort=True):
        group = group.sort_values("building_id")
        angles = np.linspace(0, 2*np.pi, len(group), endpoint=False)
        for row, angle in zip(group.itertuples(), angles):
            offset = transforms.ScaledTranslation(
                5.5*np.cos(angle)/72, 5.5*np.sin(angle)/72, fig.dpi_scale_trans)
            ax.scatter([row.longitude], [row.latitude], s=20,
                       color=COLORS[row.structural_damage_code], ec="#263746",
                       lw=.25, alpha=.95, zorder=4, transform=ax.transData + offset)
    ax.set(xlim=(34.9, 39.3), ylim=(35.7, 38.7), xlabel="Longitude (°E)", ylabel="Latitude (°N)")
    ax.grid(color="#E2E7EA", lw=.45)
    ax.set_aspect(1/np.cos(np.deg2rad(37.2)))
    event_handles, event_labels = ax.get_legend_handles_labels()
    event_legend = ax.legend(event_handles, event_labels, loc="upper left",
                             frameon=True, facecolor="white", edgecolor="#C7CED3",
                             framealpha=.92, fontsize=7.1, ncol=2,
                             columnspacing=.9, handletextpad=.4, borderpad=.45)
    ax.add_artist(event_legend)
    damage_handles = [Line2D([0], [0], marker="o", linestyle="none", markersize=5.2,
                             markerfacecolor=COLORS[code], markeredgecolor="#263746",
                             markeredgewidth=.35, label=code)
                      for code in ["N", "L", "M", "S", "C"]]
    ax.legend(handles=damage_handles, title="ACI 133 damage code", loc="lower right",
              frameon=True, facecolor="white", edgecolor="#C7CED3", framealpha=.92,
              fontsize=7.2, title_fontsize=7.3, ncol=5, columnspacing=.65,
              handletextpad=.25, borderpad=.45)
    fig.tight_layout()
    save(fig, "Fig2_study_area")


def plot_damage_counts():
    b = pd.read_csv(ROOT / "data/processed/aci133_buildings.csv")
    counts = b.structural_damage_code.value_counts().reindex(["N", "L", "M", "S", "C"], fill_value=0)
    fig, axb = plt.subplots(figsize=(3.45, 2.75))
    # Category identity is already encoded by the N--C x-axis labels; a
    # single neutral fill avoids redundant colour coding and prints cleanly.
    axb.bar(counts.index, counts.values, color=PRIMARY_FILL,
            edgecolor=PRIMARY, linewidth=.8)
    for i, value in enumerate(counts.values):
        axb.text(i, value + 3, str(value), ha="center", fontsize=8)
    axb.set(xlabel="ACI 133 structural-damage code", ylabel="Buildings", ylim=(0, 165))
    journal_axes(axb)
    fig.tight_layout()
    save(fig, "Fig2b_damage_counts")


def plot_taxonomy():
    d = pd.read_csv(DATA / "Table_common_sample_taxonomy.csv")
    d = d[d.comparison == "nearest_all_taxonomies"].copy()
    parts = d.scenario.str.split("__", expand=True)
    d["system"], d["design"] = parts[0], parts[1]
    order = ["lower", "baseline", "upper"]
    systems = ["lfinf", "lfm", "lwal"]
    names = {"lfinf": "Infilled frame", "lfm": "Moment frame", "lwal": "Wall system"}
    fig, ax = plt.subplots(figsize=(6.5, 3.2))
    x = np.arange(len(order)); width = .23
    edge_system = {"lfinf": PRIMARY, "lfm": OBSERVED, "lwal": "#666666"}
    for j, system in enumerate(systems):
        z = d[d.system == system].set_index("design").reindex(order)
        ax.bar(x + (j-1)*width, z.bias_pred_obs, width=width, color=SYSTEM_COLORS[system],
               label=names[system], edgecolor=edge_system[system], linewidth=.8)
        for xi, value in zip(x + (j-1)*width, z.bias_pred_obs):
            ax.text(xi, value + .018, f"{value:.2f}", ha="center", fontsize=7)
    ax.axhline(1, color="#4F5963", lw=1, ls=(0, (4, 3)))
    ax.set_xticks(x, ["Lower", "Baseline", "Upper"])
    ax.set(xlabel="Design/enforcement interpretation", ylabel="Predicted / observed loss proxy",
           ylim=(0, 1.08))
    ax.text(.99, .80, "Common sample: n=91\nNearest-height sensitivity",
            transform=ax.transAxes, ha="right", va="top", fontsize=7.3,
            bbox=dict(boxstyle="round,pad=.25", fc="white", ec="#BFC8CE", lw=.6))
    journal_axes(ax)
    ax.legend(frameon=False, ncol=3, fontsize=7.5, loc="upper center")
    fig.tight_layout(); save(fig, "Fig3_taxonomy_common_sample")


def plot_calibration():
    c = pd.read_csv(DATA / "Table_fragility_calibration_primary.csv")
    fig, ax = plt.subplots(figsize=(3.7, 3.35))
    ax.plot([0, 1], [0, 1], color="#8B949E", lw=1, ls=(0, (4, 3)), zorder=1)
    xerr = np.vstack([c.predicted_rate-c.pred_ci_low, c.pred_ci_high-c.predicted_rate])
    yerr = np.vstack([c.observed_rate-c.obs_ci_low, c.obs_ci_high-c.observed_rate])
    ax.errorbar(c.predicted_rate, c.observed_rate, xerr=xerr, yerr=yerr, fmt="o",
                color=PRIMARY, ecolor=PRIMARY, markerfacecolor="white",
                markeredgewidth=1.1, capsize=2.5, ms=5.5, lw=1, zorder=3)
    offsets = {"DS1+": (-34, 6), "DS2+": (6, 5), "DS3+": (6, 5), "DS4+": (8, 7)}
    for row in c.itertuples():
        ax.annotate(row.threshold, (row.predicted_rate, row.observed_rate),
                    xytext=offsets[row.threshold], textcoords="offset points", fontsize=7.5)
    ax.set(xlim=(-.02, 1.02), ylim=(-.02, 1.02),
           xlabel="Mean predicted exceedance probability",
           ylabel="Observed exceedance frequency")
    ax.text(.04, .94, "Two mainshocks; n=109", transform=ax.transAxes, fontsize=7.2)
    journal_axes(ax)
    fig.tight_layout(); save(fig, "Fig4_fragility_calibration")


def plot_im_bins():
    b = pd.read_csv(DATA / "Table_im_bin_primary.csv")
    fig, axes = plt.subplots(1, 2, figsize=(7.0, 3.15), sharex=True)
    specs = [
        (axes[0], "predicted_loss", "observed_loss", "Mean loss ratio", "(a) Loss-ratio calibration"),
        (axes[1], "predicted_ds3plus", "observed_ds3plus", "DS3+ probability / frequency", "(b) DS3+ calibration")]
    aliases = {"predicted_loss": "pred_loss", "observed_loss": "obs_loss",
               "predicted_ds3plus": "pred_ds3", "observed_ds3plus": "obs_ds3"}
    x = np.arange(1, len(b) + 1)
    tick_labels = [f"G{i}\n{value:.2f} g" for i, value in enumerate(b.mean_im_g, start=1)]
    for ax, pred, obs, ylabel, panel_title in specs:
        series = [(pred, "GEM prediction", PRIMARY, "o", -0.09),
                  (obs, "Observed", OBSERVED, "D", 0.09)]
        for col, label, color, marker, x_offset in series:
            low, high = b[f"{aliases[col]}_low"], b[f"{aliases[col]}_high"]
            err = np.vstack([b[col]-low, high-b[col]])
            ax.errorbar(x + x_offset, b[col], yerr=err, fmt=marker,
                        color=color, ecolor=color, elinewidth=1.15,
                        markerfacecolor="white" if col == pred else color,
                        markeredgewidth=1.15, ms=5, capsize=2.5,
                        linestyle="none", label=label, zorder=3)
        ax.set_xticks(x, tick_labels)
        ax.set(xlim=(0.55, 5.45), xlabel="SA(0.6 s) group and group mean", ylabel=ylabel)
        upper = max(b[f"{aliases[pred]}_high"].max(),
                    b[f"{aliases[obs]}_high"].max())
        # A small lower margin keeps near-zero G1 markers visibly separate
        # from the bottom axis without changing the displayed estimates.
        ax.set_ylim(bottom=-0.035, top=upper * 1.12)
        ax.set_title(panel_title, loc="left", fontsize=8.2, fontweight="bold", pad=5)
        journal_axes(ax)

    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", bbox_to_anchor=(.5, 1.0),
               frameon=False, fontsize=7.5, ncol=2, columnspacing=1.6)
    fig.tight_layout(w_pad=1.7, rect=(0, 0, 1, .91)); save(fig, "Fig5_im_bin_calibration")


def plot_beta():
    d = pd.read_csv(DATA / "Table_beta_sensitivity.csv")
    fig, axes = plt.subplots(1, 2, figsize=(6.8, 2.8))
    for ax, col, label in [(axes[0], "brier", "Multiclass Brier score"),
                           (axes[1], "ranked_probability_score", "Ranked probability score")]:
        ax.plot(d.beta_scale, d[col], color=PRIMARY, lw=1.7)
        best = d.loc[d[col].idxmin()]
        ax.scatter(best.beta_scale, best[col], color=OBSERVED, s=28, zorder=3)
        ax.axvline(1, color="#6B747C", lw=1, ls=(0, (4, 3)))
        # Keep the label box clear of the response curve in both panels.  The
        # boundary minimum in panel (b) needs a higher leftward callout.
        label_offset = (-20, 25) if best.beta_scale < 2.8 else (-48, 60)
        ax.annotate(f"minimum at {best.beta_scale:.2g}×",
                    (best.beta_scale, best[col]), xytext=label_offset,
                    textcoords="offset points", fontsize=7, ha="center", va="bottom",
                    bbox=dict(boxstyle="round,pad=.18", fc="white", ec="#C7CED3", lw=.5),
                    arrowprops=dict(arrowstyle="-", color="#6B747C", lw=.65,
                                    shrinkA=2, shrinkB=4))
        ax.set(xlabel="Multiplier on published total beta", ylabel=label)
        journal_axes(ax)
    fig.tight_layout(w_pad=1.6); save(fig, "Fig6_beta_sensitivity")


def plot_replication():
    a = pd.read_csv(DATA / "Table_fragility_calibration_primary.csv")
    g = pd.read_csv(DATA / "Table_gaziantep_fragility_calibration.csv")
    x = np.arange(4); width = .19
    fig, ax = plt.subplots(figsize=(6.4, 3.0))
    ax.bar(x-1.5*width, a.predicted_rate, width, color=PRIMARY_FILL,
           edgecolor=PRIMARY, linewidth=.9, label="ACI: GEM")
    ax.bar(x-.5*width, a.observed_rate, width, color=OBSERVED_FILL,
           edgecolor=OBSERVED, linewidth=.9, label="ACI: observed")
    ax.bar(x+.5*width, g.predicted_rate, width, color="white",
           edgecolor=PRIMARY, linewidth=.9, hatch="////", label="Gaziantep: GEM")
    ax.bar(x+1.5*width, g.observed_rate, width, color="white",
           edgecolor=OBSERVED, linewidth=.9, hatch="////", label="Gaziantep: observed")
    ax.set_xticks(x, a.threshold)
    ax.set(xlabel="Damage-state exceedance", ylabel="Probability or observed frequency", ylim=(0, 1.0))
    journal_axes(ax)
    ax.legend(frameon=False, ncol=2, fontsize=7.2)
    fig.tight_layout(); save(fig, "Fig7_replication")


def main():
    plot_study_map(); plot_damage_counts(); plot_taxonomy(); plot_calibration()
    plot_im_bins(); plot_beta(); plot_replication()
    print(f"Wrote revised publication figures to {DEST}")


if __name__ == "__main__":
    main()
