"""Analyses added after an internal reviewer-style audit.

All ACI outputs are reproducible from distributed inputs.  No Gaziantep
building-level output is written because that source is CC BY-NC-ND 4.0.
"""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from aci133_analysis import DLR, fragility_exceedance, load_fragilities
from fragility_tools import score_states, state_probabilities


ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "outputs" / "aci133_exploratory" / "building_results.csv"
BUILDINGS = ROOT / "data" / "processed" / "aci133_buildings.csv"
OUT = ROOT / "outputs" / "manuscript_enhanced"
SUBMISSION = ROOT / "manuscript" / "AJSE_submission"
BASE = "lfinf__baseline__exact"
SCOPE = "two_main_events"
SEED = 20230929


def write(frame: pd.DataFrame, name: str) -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    frame.to_csv(OUT / name, index=False, encoding="utf-8-sig")


def primary(results: pd.DataFrame) -> pd.DataFrame:
    return results[(results.scenario == BASE) &
                   (results.event_scope == SCOPE) &
                   (results.status == "eligible")].copy()


def score_reference_and_rank_ci(results: pd.DataFrame) -> None:
    g = primary(results)
    y = g.observed_ds.to_numpy(int)
    p = g[[f"fragility_max_p{i}" for i in range(5)]].to_numpy(float)
    model = score_states(p, y)
    frequencies = np.bincount(y, minlength=5) / len(y)
    climatology = np.tile(frequencies, (len(y), 1))
    reference = score_states(climatology, y)
    rows = []
    for metric in ["brier", "ranked_probability_score", "negative_log_score"]:
        rows.append({
            "metric": metric,
            "model_score": model[metric],
            "climatology_score": reference[metric],
            "skill_score_1_minus_model_over_reference": 1 - model[metric] / reference[metric],
        })
    write(pd.DataFrame(rows), "Table_climatology_skill.csv")

    observed_loss = g.observed_ds.map(DLR).to_numpy(float)
    predicted_loss = g.predicted_lr_max.to_numpy(float)
    point = spearmanr(predicted_loss, observed_loss).statistic
    rng = np.random.default_rng(SEED)
    boot = []
    for _ in range(10000):
        idx = rng.integers(0, len(g), len(g))
        value = spearmanr(predicted_loss[idx], observed_loss[idx]).statistic
        if np.isfinite(value):
            boot.append(value)
    write(pd.DataFrame([{
        "n": len(g), "spearman_rho": point,
        "bootstrap_ci_low": np.quantile(boot, 0.025),
        "bootstrap_ci_high": np.quantile(boot, 0.975),
        "bootstrap_replicates": len(boot),
    }]), "Table_spearman_bootstrap.csv")


def median_beta_surface(results: pd.DataFrame) -> None:
    g = primary(results)
    y = g.observed_ds.to_numpy(int)
    fragilities = load_fragilities()
    theta_scales = np.round(np.arange(0.40, 1.601, 0.05), 2)
    beta_scales = np.round(np.arange(0.50, 5.001, 0.10), 2)
    rows = []
    for theta_scale in theta_scales:
        for beta_scale in beta_scales:
            probabilities = []
            for row in g.itertuples():
                base_taxonomy = row.taxonomy.rsplit("/", 1)[0]
                pars = fragilities[(base_taxonomy, row.imt)]
                scaled = {
                    "median": pars["median"] * theta_scale,
                    "beta": pars["beta"] * beta_scale,
                }
                exceed, _ = fragility_exceedance(scaled, [row.max_im_g])
                probabilities.append(state_probabilities(exceed)[0])
            p = np.asarray(probabilities)
            scores = score_states(p, y)
            rows.append({"theta_scale": theta_scale, "beta_scale": beta_scale, **scores})
    surface = pd.DataFrame(rows)
    write(surface, "Table_theta_beta_surface.csv")

    summary = []
    for metric in ["brier", "ranked_probability_score", "negative_log_score"]:
        best = surface.loc[surface[metric].idxmin()]
        summary.append({"metric": metric, "theta_scale": best.theta_scale,
                        "beta_scale": best.beta_scale, "minimum_score": best[metric],
                        "theta_at_boundary": best.theta_scale in (theta_scales.min(), theta_scales.max()),
                        "beta_at_boundary": best.beta_scale in (beta_scales.min(), beta_scales.max())})
    write(pd.DataFrame(summary), "Table_theta_beta_optima.csv")

    # A compact two-panel surface for Brier and RPS, the metrics used for skill.
    plt.rcParams.update({"font.family": "serif", "font.serif": ["STIXGeneral"],
                         "mathtext.fontset": "stix", "font.size": 8.5,
                         "pdf.fonttype": 42, "ps.fonttype": 42})
    fig, axes = plt.subplots(1, 2, figsize=(6.7, 3.0), constrained_layout=True)
    for ax, metric, title in zip(
            axes, ["brier", "ranked_probability_score"],
            ["Multiclass Brier score", "Ranked probability score"]):
        z = surface.pivot(index="beta_scale", columns="theta_scale", values=metric)
        im = ax.contourf(z.columns, z.index, z.to_numpy(), levels=18, cmap="viridis_r")
        best = surface.loc[surface[metric].idxmin()]
        ax.plot(1, 1, marker="x", color="white", mec="black", mew=0.7, ms=6,
                label="Published")
        ax.plot(best.theta_scale, best.beta_scale, marker="o", color="#B85F47",
                mec="white", mew=0.7, ms=5, label="Minimum")
        ax.set(xlabel=r"Multiplier on published median $\theta$",
               ylabel=r"Multiplier on published total $\beta$", title=title)
        ax.legend(frameon=False, fontsize=7, loc="upper right")
        fig.colorbar(im, ax=ax, pad=0.02)
    fig.savefig(SUBMISSION / "Fig_theta_beta_surface.pdf", bbox_inches="tight")
    fig.savefig(SUBMISSION / "Fig_theta_beta_surface.png", dpi=600, bbox_inches="tight")
    plt.close(fig)


def selection_and_g1(results: pd.DataFrame, buildings: pd.DataFrame) -> None:
    base = results[(results.scenario == BASE) & (results.event_scope == SCOPE)].copy()
    fields = ["building_id", "status", "observed_ds", "structural_damage_code",
              "storeys_reported", "age_years_reported"]
    base = base[fields].drop_duplicates("building_id")
    rows = []
    for name, z in [("full_source", base),
                    ("primary_exact_height", base[base.status == "eligible"]),
                    ("excluded_height_unmapped", base[base.status != "eligible"])]:
        storeys = z.storeys_reported.to_numpy(float)
        ds = z.observed_ds.to_numpy(int)
        rows.append({
            "cohort": name, "n": len(z),
            "storeys_median": np.median(storeys),
            "storeys_q1": np.quantile(storeys, 0.25),
            "storeys_q3": np.quantile(storeys, 0.75),
            "DS0_n": int((ds == 0).sum()), "DS1_n": int((ds == 1).sum()),
            "DS2_n": int((ds == 2).sum()), "DS3_n": int((ds == 3).sum()),
            "DS4_n": int((ds == 4).sum()), "DS3plus_percent": 100 * (ds >= 3).mean(),
        })
    write(pd.DataFrame(rows), "Table_exact_height_selection.csv")

    # Composition of the lowest SA(0.6 s) group, needed to rule out unit errors.
    g = primary(results)
    g = g[g.imt == "SA(0.6)"].copy()
    g["im_group"] = pd.qcut(g.max_im_g, 5, labels=False, duplicates="drop") + 1
    g1 = g[g.im_group == 1].copy()
    design = g1.taxonomy.str.extract(r"CR/LFINF/([^/]+)/H:")[0]
    height = g1.taxonomy.str.extract(r"/H:([^/]+)/")[0]
    detail = pd.DataFrame({
        "category": ["n", "mean_SA06_g", "mean_predicted_loss", "mean_observed_loss",
                     "mean_predicted_DS3plus", "observed_DS3plus_frequency"],
        "value": [len(g1), g1.max_im_g.mean(), g1.predicted_lr_max.mean(),
                  g1.observed_ds.map(DLR).mean(),
                  (g1.fragility_max_p3 + g1.fragility_max_p4).mean(),
                  (g1.observed_ds >= 3).mean()],
    })
    for key, count in design.value_counts().items():
        detail.loc[len(detail)] = [f"design_{key}_n", count]
    for key, count in height.value_counts().items():
        detail.loc[len(detail)] = [f"height_H{key}_n", count]
    write(detail, "Table_G1_composition.csv")

    years = buildings.construction_year_estimated
    write(pd.DataFrame([{
        "n": len(buildings), "minimum_estimated_construction_year": years.min(),
        "pre_1975_count": int((years < 1975).sum()),
    }]), "Table_ACI_construction_year_check.csv")


def main() -> None:
    results = pd.read_csv(RESULTS)
    buildings = pd.read_csv(BUILDINGS)
    score_reference_and_rank_ci(results)
    median_beta_surface(results)
    selection_and_g1(results, buildings)
    print("Reviewer-enhancement outputs written to", OUT)


if __name__ == "__main__":
    main()
