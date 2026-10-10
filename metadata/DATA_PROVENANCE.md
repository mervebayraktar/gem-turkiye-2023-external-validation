# Data provenance and selection rationale

## ACI133 building data

The local file `data/raw/aci133/ACI133v3.xlsx` was downloaded from Zenodo record `10.5281/zenodo.13386343`, version v1, published 28 August 2024 by Mohammad Amin Hariri-Ardebili. Direct file endpoint:

`https://zenodo.org/api/records/13386343/files/ACI133v3.xlsx/content`

The verified file checksum is `md5:c912ccae688adf77339571c621f31e2d`. The record is licensed CC BY 4.0.

This is not described as the untouched original ACI/ASCE survey workbook. The Zenodo record states that the building and damage data were primarily developed under the ACI Committee 133 initiative and presented by Pujol et al. (2024), while several missing values were supplemented by the Zenodo dataset developer using further investigation and engineering judgment. The workbook marks 84 supplemented cells affecting 42 buildings; these flags are retained.

The counts N=19, L=42, M=24, S=151 and C=6 come directly from the workbook column `Structural Damage (5-class)`. The analysis did not create these source codes. The conversion N/L/M/S/C to GEM DS0-DS4 and the damage-to-loss ratios are analysis mappings. The meaning of C is unresolved because the workbook also contains a separate `Critical Damage` field; C-to-DS4 and C-to-DS3 are therefore separate sensitivities.

## GEM vulnerability and fragility models

The country vulnerability XML files were downloaded from the public GEM repository:

`https://github.com/gem/global_vulnerability_model`

The analysis is pinned to commit `5974372ac3f4a99f25d0649eb030fbe596f23b36`, under `Europe/Turkey`. The structural file used in the analysis is:

`https://raw.githubusercontent.com/gem/global_vulnerability_model/5974372ac3f4a99f25d0649eb030fbe596f23b36/Europe/Turkey/vulnerability_structural.xml`

The repository license is CC BY-NC-SA 4.0.

Fragility NRML files are not distributed in the repository. Tabulated fragility parameters for PGA, SA(0.3), SA(0.6) and SA(1.0) were downloaded from the official OpenQuake vulnerability documentation under `world_v4/GEM/summaries`. All 186 base structural taxonomies used by the Turkey structural XML have parameter matches. The study proceeds with these openly published official tables and records the absence of an explicit release-pairing statement as a provenance limitation; no private clarification is required for the planned analysis.

## Earthquake selection

Four events were selected because they occurred before the field survey and are represented as event-specific intensity-measure blocks in the ACI133/Hariri-Ardebili dataset and associated study, rather than through an arbitrary magnitude cut-off:

| Analysis label | USGS event ID | USGS magnitude/time (UTC) | Role |
|---|---|---|---|
| event1 | us6000jllz | Mw 7.8, 2023-02-06 01:17:34 | first mainshock |
| event3 | us6000jlm1 | Mw 6.7, 2023-02-06 01:28:15 | strong early aftershock |
| event2 | us6000jlqa | Mw 7.5, 2023-02-06 10:24:48 | second mainshock |
| event4 | us6000jqcn | Mw 6.3, 2023-02-20 17:04:29 | later damaging event before survey |

USGS magnitudes are retained because the ground-motion grids are USGS ShakeMap products. AFAD reports different values for some events (notably Mw 7.7 and 7.6 for the mainshocks and approximately Mw 6.4 for the 20 February event); the manuscript must state the magnitude agency whenever values are quoted.

Two event scopes are evaluated: the two mainshocks and all four survey-preceding events. Maximum IM is the primary practical representation. It is not a physical sequential-damage model. Thousands of smaller aftershocks are not individually included because the source dataset does not supply a consistent event-specific IM block for each and intact-building fragility functions do not define damaged-state transitions.

## Vs30

`V_S30 (m/s)` is a direct column in `ACI133v3.xlsx`; the analysis stores it as `vs30_source_mps`. It was not calculated by the present workflow. The associated paper indicates that the velocity profile used in its ShakeMap/kriging comparison was provided by ShakeMap, so this field must not be described as a building-site field measurement unless the authors provide contrary metadata.

AFAD has published a nationwide GIS-based predicted-Vs30 dataset developed under UDAP-Ç-20-01, with 70 x 70 m rasters for best-fit, one-standard-deviation-low and quantile scenarios. This is suitable as an independent Vs30 sensitivity source and should be cited to the AFAD (2025) project report and Şahin et al. (2024), DOI `10.1785/0220230321`. AFAD cautions that the product is not a parcel-scale site investigation.

The AFAD Türkiye Earthquake Hazard Map is a different product: its hazard maps use the reference condition Vs30=760 m/s and do not represent local site-specific Vs30. It should not be substituted for the predicted-Vs30 raster.

Recommended use in this project:

1. retain the workbook/ShakeMap-related Vs30 as the source baseline;
2. sample AFAD best-fit, 1-sigma-low and quantile rasters at all building coordinates;
3. compare the AFAD and source values and repeat residual-correlation diagnostics;
4. do not directly replace Vs30 inside already site-adjusted ShakeMap IMs, which could double-apply site effects;
5. if AFAD Vs30 is to alter ground-motion estimates, rebuild the ground-motion field from TADAS records or a documented GMM/site-amplification workflow.
