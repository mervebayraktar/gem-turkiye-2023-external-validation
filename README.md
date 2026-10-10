# External assessment of GEM v2026 functions in Türkiye

This repository contains the reproducible analysis package for a building-level external assessment of the GEM v2026 structural fragility and vulnerability functions using damage observed after the 2023 Kahramanmaraş earthquake sequence.

## Scope

The workflow links the open ACI 133 reconnaissance dataset to versioned USGS ShakeMaps and published GEM functions. It evaluates vulnerability means, fragility probabilities, taxonomy sensitivity, event-coverage sensitivity, calibration by intensity group, and joint median-dispersion score surfaces. An independent Gaziantep dataset is used only for aggregate replication because its licence prohibits distribution of transformed microdata.

## Repository contents

- Repository root: acquisition, model parsing, analysis, sensitivity and figure scripts.
- `config/`: damage, loss and taxonomy mapping assumptions.
- `outputs/tables/`: non-restricted aggregate results used in the manuscript.
- `outputs/figures/`: selected generated figures.
- `supplementary/`: exact IM-group membership table referenced by the manuscript.
- `metadata/`: provenance, source versions, licences and checksums.

## Data access

Third-party source data and model files are not redistributed here.

- ACI 133: Zenodo record [10.5281/zenodo.13386343](https://doi.org/10.5281/zenodo.13386343), CC BY 4.0.
- GEM global vulnerability model: [gem/global_vulnerability_model](https://github.com/gem/global_vulnerability_model), commit `5974372ac3f4a99f25d0649eb030fbe596f23b36`, CC BY-NC-SA 4.0.
- GEM fragility parameter summaries: [OpenQuake vulnerability documentation](https://docs.openquake.org/vulnerability/fragility_models.html).
- USGS ShakeMaps: events `us6000jllz`, `us6000jlm1`, `us6000jlqa`, and `us6000jqcn`.
- Gaziantep replication data: [10.7910/DVN/OQCPWF](https://doi.org/10.7910/DVN/OQCPWF), CC BY-NC-ND 4.0. No source or transformed building-level Gaziantep records are included.

See `metadata/DATA_PROVENANCE.md` for version and interpretation details.

## Environment

The analysis was developed with Python 3.11. Create an isolated environment and install:

```bash
python -m venv .venv
python -m pip install -r requirements.txt
```

The scripts expect downloaded inputs under a project-local `data/raw/` directory. Paths are resolved relative to the repository; no user-specific absolute paths are required. Aggregate manuscript enhancement tables can be regenerated after the source inputs have been obtained:

```bash
python aci133_analysis.py
python enhance_manuscript_analysis.py
python reviewer_enhancements.py
```

## Interpretation

The reported factors are conditional on the surveyed samples, taxonomy mappings and maximum-IM event representation. They are not national correction factors. The published 0.30 modelling component cannot be updated separately from the available public parameter tables.

## Licence

Original code in this repository is released under the MIT License. Tables and documentation authored for this study are available under CC BY 4.0 unless a file states otherwise. Third-party inputs retain their original licences and are referenced rather than redistributed.

## Citation

Citation metadata are provided in `CITATION.cff`. A version-specific DOI will be added after the first archived release.

