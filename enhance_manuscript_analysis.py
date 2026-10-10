"""Publication-focused analyses for the AJSE manuscript.

This module leaves both source workbooks unchanged.  It produces only aggregate
outputs for the CC BY-NC-ND Gaziantep dataset and never writes derived
building-level Gaziantep records.
"""
from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from scipy.stats import spearmanr

from aci133_analysis import (
    DESIGN_MAPS,
    DLR,
    IM_FIELDS,
    RAW,
    fragility_exceedance,
    load_fragilities,
    vulnerability_metrics,
)
from fragility_tools import score_states, state_probabilities
from pipeline import evaluate, read_grid, read_model, sample_grid


ROOT = Path(__file__).resolve().parent
ACI_OUT = ROOT / "outputs/aci133_exploratory"
OUT = ROOT / "outputs/manuscript_enhanced"
GAZIANTEP_SOURCE = RAW / "rapor_final_v2.xlsx"
BASE = "lfinf__baseline__exact"
PRIMARY_SCOPE = "two_main_events"
MAPPING = "code_C_as_DS4"
EVENT_KEYS = ["event1", "event2"]
EVENT_NAMES = {"event1": "6 Feb Mw 7.8", "event2": "6 Feb Mw 7.5"}


def write(frame: pd.DataFrame, name: str) -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    frame.to_csv(OUT / name, index=False, encoding="utf-8-sig")


def mapped_ds(values: pd.Series) -> np.ndarray:
    return values.to_numpy(int)


def fragility_scores(frame: pd.DataFrame, prefix: str = "fragility_max_p") -> dict:
    y = mapped_ds(frame.observed_ds)
    p = frame[[f"{prefix}{i}" for i in range(5)]].to_numpy(float)
    scores = score_states(p, y)
    scores["accuracy_argmax"] = float((p.argmax(axis=1) == y).mean())
    return scores


def primary_tables(results: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for scope in ["two_main_events", "four_survey_preceding_events"]:
        g = results[(results.scenario == BASE) & (results.event_scope == scope) &
                    (results.status == "eligible")].copy()
        y = mapped_ds(g.observed_ds)
        yloss = pd.Series(y).map(DLR).to_numpy(float)
        rows.append({"event_scope": scope, "output": "vulnerability",
                     **vulnerability_metrics(yloss, g.predicted_lr_max)})
        rows.append({"event_scope": scope, "output": "fragility", "n": len(g),
                     **fragility_scores(g)})
    out = pd.DataFrame(rows)
    write(out, "Table_primary_performance.csv")
    return out


def common_taxonomy(results: pd.DataFrame) -> pd.DataFrame:
    """Hold buildings fixed while changing taxonomy assumptions.

    The all-scenario exact-height intersection has only three buildings.  The
    headline factorial comparison therefore uses nearest-height mappings on the
    common 91-building intersection.  Exact-height comparisons are retained in
    matched blocks that vary one taxonomy dimension at a time.
    """
    e = results[(results.status == "eligible") &
                (results.event_scope == PRIMARY_SCOPE)].copy()
    rows = []

    nearest = e[e.height_policy == "nearest"]
    scenarios = sorted(nearest.scenario.unique())
    common = set.intersection(*[set(nearest[nearest.scenario == s].building_id) for s in scenarios])
    for scenario in scenarios:
        g = nearest[(nearest.scenario == scenario) & nearest.building_id.isin(common)].copy()
        yloss = g.observed_ds.map(DLR).to_numpy(float)
        rows.append({"comparison": "nearest_all_taxonomies", "block": "common_91",
                     "scenario": scenario, "common_n": len(common),
                     **vulnerability_metrics(yloss, g.predicted_lr_max),
                     **{f"fragility_{k}": v for k, v in fragility_scores(g).items()}})

    exact = e[e.height_policy == "exact"]
    for design in sorted(exact.design_scenario.unique()):
        z = exact[exact.design_scenario == design]
        ss = sorted(z.scenario.unique())
        common_ids = set.intersection(*[set(z[z.scenario == s].building_id) for s in ss])
        for scenario in ss:
            g = z[(z.scenario == scenario) & z.building_id.isin(common_ids)]
            rows.append({"comparison": "exact_system_effect", "block": design,
                         "scenario": scenario, "common_n": len(common_ids),
                         **vulnerability_metrics(g.observed_ds.map(DLR), g.predicted_lr_max)})
    for system in sorted(exact.system_scenario.unique()):
        z = exact[exact.system_scenario == system]
        ss = sorted(z.scenario.unique())
        common_ids = set.intersection(*[set(z[z.scenario == s].building_id) for s in ss])
        for scenario in ss:
            g = z[(z.scenario == scenario) & z.building_id.isin(common_ids)]
            rows.append({"comparison": "exact_design_effect", "block": system,
                         "scenario": scenario, "common_n": len(common_ids),
                         **vulnerability_metrics(g.observed_ds.map(DLR), g.predicted_lr_max)})
    out = pd.DataFrame(rows)
    write(out, "Table_common_sample_taxonomy.csv")
    return out


def cluster_bootstrap(results: pd.DataFrame, buildings: pd.DataFrame) -> pd.DataFrame:
    g = results[(results.scenario == BASE) & (results.event_scope == PRIMARY_SCOPE) &
                (results.status == "eligible")].merge(
                    buildings[["building_id", "latitude", "longitude"]], on="building_id", how="left")
    g["coordinate_cluster"] = g.latitude.round(6).astype(str) + "_" + g.longitude.round(6).astype(str)
    g["spatial_block"] = np.floor(g.latitude / 0.10).astype(int).astype(str) + "_" + \
        np.floor(g.longitude / 0.10).astype(int).astype(str)
    rng = np.random.default_rng(20230929)
    rows = []
    for scheme, col in [("building", None), ("coordinate_cluster", "coordinate_cluster"),
                        ("0.1_degree_spatial_block", "spatial_block")]:
        factors = []
        if col is None:
            groups = [np.array([i]) for i in range(len(g))]
        else:
            groups = [np.asarray(idx, int) for idx in g.groupby(col).indices.values()]
        for _ in range(5000):
            chosen = rng.integers(0, len(groups), len(groups))
            idx = np.concatenate([groups[i] for i in chosen])
            pred = g.predicted_lr_max.to_numpy(float)[idx].sum()
            obs = g.observed_ds.map(DLR).to_numpy(float)[idx].sum()
            if pred > 0:
                factors.append(obs / pred)
        point = g.observed_ds.map(DLR).sum() / g.predicted_lr_max.sum()
        lo, hi = np.quantile(factors, [0.025, 0.975])
        rows.append({"bootstrap_scheme": scheme, "n_buildings": len(g),
                     "n_sampling_units": len(groups), "observed_over_predicted": point,
                     "ci_low": lo, "ci_high": hi, "replicates": len(factors)})
    out = pd.DataFrame(rows)
    write(out, "Table_bootstrap_sensitivity.csv")
    return out


def beta_sensitivity(results: pd.DataFrame) -> pd.DataFrame:
    fragilities = load_fragilities()
    g = results[(results.scenario == BASE) & (results.event_scope == PRIMARY_SCOPE) &
                (results.status == "eligible")].copy()
    y = mapped_ds(g.observed_ds)
    rows = []
    for scale in np.round(np.arange(0.50, 3.001, 0.05), 2):
        probabilities = []
        for row in g.itertuples():
            base_taxonomy = row.taxonomy.rsplit("/", 1)[0]
            parameters = fragilities[(base_taxonomy, row.imt)]
            scaled = {"median": parameters["median"], "beta": parameters["beta"] * scale}
            exceed, _ = fragility_exceedance(scaled, [row.max_im_g])
            probabilities.append(state_probabilities(exceed)[0])
        p = np.asarray(probabilities)
        scores = score_states(p, y)
        scores["accuracy_argmax"] = float((p.argmax(axis=1) == y).mean())
        rows.append({"beta_scale": scale, "n": len(g), **scores})
    out = pd.DataFrame(rows)
    for metric in ["brier", "ranked_probability_score", "negative_log_score"]:
        out[f"best_{metric}"] = out[metric].eq(out[metric].min())
    write(out, "Table_beta_sensitivity.csv")
    return out


def primary_calibration_with_intervals(results: pd.DataFrame) -> pd.DataFrame:
    g = results[(results.scenario == BASE) & (results.event_scope == PRIMARY_SCOPE) &
                (results.status == "eligible")].copy()
    y = mapped_ds(g.observed_ds)
    p = g[[f"fragility_max_p{i}" for i in range(5)]].to_numpy(float)
    exceed = 1 - np.cumsum(p, axis=1)[:, :4]
    rng = np.random.default_rng(20230929)
    rows = []
    for ds in range(1, 5):
        obs = (y >= ds).astype(float)
        pred = exceed[:, ds-1]
        boot_obs, boot_pred = [], []
        for _ in range(5000):
            idx = rng.integers(0, len(g), len(g))
            boot_obs.append(obs[idx].mean())
            boot_pred.append(pred[idx].mean())
        rows.append({"threshold": f"DS{ds}+", "n": len(g),
                     "predicted_rate": pred.mean(),
                     "pred_ci_low": np.quantile(boot_pred, 0.025),
                     "pred_ci_high": np.quantile(boot_pred, 0.975),
                     "observed_rate": obs.mean(),
                     "obs_ci_low": np.quantile(boot_obs, 0.025),
                     "obs_ci_high": np.quantile(boot_obs, 0.975),
                     "binary_brier": np.mean((pred-obs)**2)})
    out = pd.DataFrame(rows)
    write(out, "Table_fragility_calibration_primary.csv")
    return out


def primary_im_bins(results: pd.DataFrame) -> pd.DataFrame:
    g = results[(results.scenario == BASE) & (results.event_scope == PRIMARY_SCOPE) &
                (results.status == "eligible") & (results.imt == "SA(0.6)")].copy()
    g["bin"], edges = pd.qcut(g.max_im_g, 5, duplicates="drop", retbins=True)
    g["im_group"] = g["bin"].cat.codes + 1
    g["bin_lower_g"] = g["im_group"].map(
        {i + 1: edges[i] for i in range(len(edges) - 1)})
    g["bin_upper_g"] = g["im_group"].map(
        {i + 1: edges[i + 1] for i in range(len(edges) - 1)})
    g["predicted_ds3plus"] = g.fragility_max_p3 + g.fragility_max_p4
    membership = g[["building_id", "im_group", "max_im_g", "bin_lower_g",
                    "bin_upper_g", "structural_damage_code", "observed_ds",
                    "predicted_lr_max", "predicted_ds3plus"]].copy()
    membership = membership.rename(columns={
        "max_im_g": "maximum_two_mainshock_sa06_g",
        "observed_ds": "observed_damage_state",
        "predicted_lr_max": "predicted_loss_ratio"})
    membership = membership.sort_values(
        ["im_group", "maximum_two_mainshock_sa06_g", "building_id"])
    write(membership, "Table_im_bin_membership.csv")
    rng = np.random.default_rng(20230929)
    rows = []
    for interval, z in g.groupby("bin", observed=True):
        obs_loss = z.observed_ds.map(DLR).to_numpy(float)
        pred_loss = z.predicted_lr_max.to_numpy(float)
        obs_ds3 = (z.observed_ds.to_numpy(int) >= 3).astype(float)
        pred_ds3 = z.fragility_max_p3.add(z.fragility_max_p4).to_numpy(float)
        boot = {"obs_loss": [], "pred_loss": [], "obs_ds3": [], "pred_ds3": []}
        for _ in range(3000):
            idx = rng.integers(0, len(z), len(z))
            for key, values in [("obs_loss", obs_loss), ("pred_loss", pred_loss),
                                ("obs_ds3", obs_ds3), ("pred_ds3", pred_ds3)]:
                boot[key].append(values[idx].mean())
        row = {"im_bin": str(interval), "n": len(z), "mean_im_g": z.max_im_g.mean(),
               "observed_loss": obs_loss.mean(), "predicted_loss": pred_loss.mean(),
               "observed_ds3plus": obs_ds3.mean(), "predicted_ds3plus": pred_ds3.mean()}
        for key in boot:
            row[f"{key}_low"], row[f"{key}_high"] = np.quantile(boot[key], [0.025, 0.975])
        rows.append(row)
    out = pd.DataFrame(rows)
    write(out, "Table_im_bin_primary.csv")
    return out


def event_inventory(ims: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for key in ["event1", "event3", "event2", "event4"]:
        event = json.loads((RAW / f"{key}.json").read_text(encoding="utf-8"))
        props = event["properties"]
        coords = event["geometry"]["coordinates"]
        z = ims[ims.event_key == key]
        coverage = z.groupby("imt").im_g.apply(lambda x: int(x.notna().sum())).to_dict()
        deriv = z.groupby("imt").im_derivation.first().to_dict()
        rows.append({"event_key": key, "usgs_id": event["id"], "magnitude": props["mag"],
                     "origin_utc": pd.to_datetime(props["time"], unit="ms", utc=True),
                     "longitude": coords[0], "latitude": coords[1],
                     "PGA_coverage": coverage.get("PGA", 0),
                     "SA03_coverage": coverage.get("SA(0.3)", 0),
                     "SA06_coverage": coverage.get("SA(0.6)", 0),
                     "SA10_coverage": coverage.get("SA(1.0)", 0),
                     "SA06_derivation": deriv.get("SA(0.6)")})
    out = pd.DataFrame(rows)
    write(out, "Table_event_inventory.csv")
    return out


def year_to_era(year: float) -> str | None:
    if not np.isfinite(year) or year < 1975:
        return None
    if year <= 1997:
        return "1975_1997"
    if year <= 2006:
        return "1998_2006"
    if year <= 2018:
        return "2007_2018"
    return "2019_plus"


def gaziantep_replication(aci_buildings: pd.DataFrame) -> pd.DataFrame:
    """Aggregate-only replication using the ND-licensed Gaziantep workbook."""
    source = pd.read_excel(GAZIANTEP_SOURCE, sheet_name="full_data")
    damage_map = {"Non-damaged": 0, "Slight Damage": 1, "Moderate Damage": 2,
                  "Heavy Damage": 3, "Urgently Demolish": 4, "Collapsed": 4}
    g = source[(source["Structural System"] == "Reinforced Concrete") &
               source["Damage"].isin(damage_map)].copy()
    g["observed_ds"] = g["Damage"].map(damage_map)
    g["era"] = pd.to_numeric(g["Construction Year"], errors="coerce").map(year_to_era)
    g["storeys"] = pd.to_numeric(g["Story Number"], errors="coerce")
    g["latitude"] = pd.to_numeric(g["Latitude"], errors="coerce")
    g["longitude"] = pd.to_numeric(g["Longitude"], errors="coerce")
    g = g[g.era.notna() & g.latitude.notna() & g.longitude.notna()].copy()

    model = read_model(RAW / "gem_vulnerability_structural.xml")
    fragilities = load_fragilities()
    grid_cache = {}
    for event in EVENT_KEYS:
        grid_cache[event] = read_grid(RAW / f"{event}_grid.xml")[:2]

    # Sample every required IM once per event.  Re-sampling the full grid for
    # each individual building is unnecessarily expensive.
    im_samples = {}
    for event in EVENT_KEYS:
        grid, units = grid_cache[event]
        for imt, field in IM_FIELDS.items():
            values = sample_grid(grid, field, g.latitude, g.longitude)
            if units[field] == "%g":
                values = values / 100.0
            im_samples[(event, imt)] = values

    records = []
    for position, row in enumerate(g.itertuples()):
        design = DESIGN_MAPS["baseline"][row.era]
        prefix = f"CR/LFINF/{design}/H:"
        allowed = sorted({int(k.split("/H:")[1].split("/")[0]) for k in model if k.startswith(prefix)})
        if not float(row.storeys).is_integer() or int(row.storeys) not in allowed:
            continue
        taxonomy = f"CR/LFINF/{design}/H:{int(row.storeys)}/RES"
        if taxonomy not in model:
            continue
        fun = model[taxonomy]
        imt = fun["imt"]
        values = [im_samples[(event, imt)][position] for event in EVENT_KEYS]
        if not np.isfinite(values).all():
            continue
        max_im = max(values)
        pred, _, outside = evaluate(fun, [max_im], outside="nan")
        if outside[0] or not np.isfinite(pred[0]):
            continue
        base = taxonomy.rsplit("/", 1)[0]
        exceed, _ = fragility_exceedance(fragilities[(base, imt)], [max_im])
        probs = state_probabilities(exceed)[0]
        records.append({"latitude": row.latitude, "longitude": row.longitude,
                        "observed_ds": row.observed_ds, "predicted_loss": pred[0],
                        **{f"p{i}": probs[i] for i in range(5)}})
    r = pd.DataFrame(records)

    lat0 = np.deg2rad(np.nanmean(np.r_[aci_buildings.latitude, r.latitude]))
    scale = np.array([111.32, 111.32 * np.cos(lat0)])
    tree = cKDTree(aci_buildings[["latitude", "longitude"]].to_numpy() * scale)
    distance, _ = tree.query(r[["latitude", "longitude"]].to_numpy() * scale)
    r["aci_nearest_km"] = distance

    rows = []
    for cohort, z in [("all_eligible", r), ("excluding_within_50m_of_ACI133", r[r.aci_nearest_km >= 0.05])]:
        y = z.observed_ds.to_numpy(int)
        yloss = pd.Series(y).map(DLR).to_numpy(float)
        p = z[[f"p{i}" for i in range(5)]].to_numpy(float)
        vm = vulnerability_metrics(yloss, z.predicted_loss)
        fs = score_states(p, y)
        fs["accuracy_argmax"] = float((p.argmax(axis=1) == y).mean())
        rows.append({"cohort": cohort, **vm, **{f"fragility_{k}": v for k, v in fs.items()},
                     "near_ACI_50m_excluded": int((r.aci_nearest_km < 0.05).sum()),
                     "source_records": len(source), "rc_evaluated_records": len(g),
                     "license": "CC BY-NC-ND 4.0; aggregate results only"})
    out = pd.DataFrame(rows)
    write(out, "Table_gaziantep_replication.csv")

    cal_rows = []
    y = r.observed_ds.to_numpy(int)
    p = r[[f"p{i}" for i in range(5)]].to_numpy(float)
    exceed = 1 - np.cumsum(p, axis=1)[:, :4]
    for ds in range(1, 5):
        cal_rows.append({"threshold": f"DS{ds}+", "n": len(r),
                         "predicted_rate": exceed[:, ds-1].mean(),
                         "observed_rate": (y >= ds).mean()})
    write(pd.DataFrame(cal_rows), "Table_gaziantep_fragility_calibration.csv")
    return out


def main() -> None:
    results = pd.read_csv(ACI_OUT / "building_results.csv")
    buildings = pd.read_csv(ROOT / "data/processed/aci133_buildings.csv")
    ims = pd.read_csv(ROOT / "data/processed/aci133_im_assignments.csv")
    primary_tables(results)
    common_taxonomy(results)
    cluster_bootstrap(results, buildings)
    beta_sensitivity(results)
    primary_calibration_with_intervals(results)
    primary_im_bins(results)
    event_inventory(ims)
    gaziantep_replication(buildings)
    manifest = {
        "primary_analysis": "two mainshocks; LFINF baseline design; exact height; code C as DS4",
        "taxonomy_common_sample": "nearest-height intersection across nine system-design scenarios",
        "beta_statement": "Total published beta scaled; the 0.30 component is not separately updated.",
        "gaziantep_distribution": "Only aggregate results written; no derived ND-licensed microdata exported.",
    }
    (OUT / "analysis_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(f"Enhanced manuscript analyses written to {OUT}")


if __name__ == "__main__":
    main()
