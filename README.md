# Wesselman Woods Tree-Ring Climate Analysis

This repository contains the analysis workflow used for the Wesselman Woods thesis chapter on tree growth, climate response, deer relationships, and species-level chronologies.

The project is organized around one main runner script:

- `wwn_main_run_all_outputs.R`

That script runs three analyses in sequence:

1. `wwn_handoff_integrated_analysis.R`
2. `wwn_species_deer_pipeline.R`
3. `wwn_species_climate_analysis.R`

## What This Project Produces

The workflow generates:

- integrated chronology and climate outputs
- MODIS-era and long-term climate correlation plots
- ENSO response outputs
- deer vs climate and deer vs growth outputs
- species-level standard and residual chronologies
- species-by-species climate matrices and regression panels
- species marker-year figures and story figures
- climate-to-climate, species-to-species, and species-to-climate heatmaps

Main output folders:

- `wwn_handoff_outputs`
- `wwn_species_deer_outputs`
- `wwn_species_climate_outputs`

## Required Input Files

The main workflow expects these files in the repository root:

- `BulkClean.rwl`
- `dated8.rwl`
- `ACER1.rwl`
- `CACA.rwl`
- `CESP.rwl`
- `FRSP.rwl`
- `JUNI.rwl`
- `LIST.rwl`
- `LITU.rwl`
- `NYSY.rwl`
- `PRSP.rwl`
- `QUER.rwl`
- `SAAL.rwl`
- `ULSP.rwl`
- `chronology_climate.xlsx`
- `pdsi.csv`
- `PRISM_ppt_tmin_tmean_tmax_tdmean_vpdmin_vpdmax_stable_4km_189501_202401_37.9810_-87.5064.csv`
- `ee-chart.csv`
- `enso_annual_nino34.csv`
- `wwn_deer.xlsx`
- `documented_flood_drought_years.csv`

Optional but included in this workspace:

- `WWSAS_Thesis_ISU.Rproj`
- `wwn_treedata.xlsx`
- `wwn_treedata_sample_table.R`
- `wwn_enso_rwi_analysis.R`
- `wwn_enso_rwi_nonlinear_analysis.R`

## R Packages

The scripts rely on these packages:

- `ggplot2`
- `dplyr`
- `dplR`
- `readxl`
- `zoo`
- `SPEI`
- `patchwork`

Some helper scripts only require a subset of these packages.

## How To Run

Run the full workflow from the repository root:

```r
Rscript wwn_main_run_all_outputs.R
```

You can also run the scripts individually:

```r
Rscript wwn_handoff_integrated_analysis.R
Rscript wwn_species_deer_pipeline.R
Rscript wwn_species_climate_analysis.R
```

## Recommended GitHub Structure

For a thesis-linked repository, the cleanest approach is:

- upload the scripts
- upload the input data files needed for reproducibility
- optionally upload the output folders if you want readers to see the finished figures without rerunning anything

See `GITHUB_UPLOAD_CHECKLIST.md` for a file-by-file checklist.

## Notes

- The workflow was rerun successfully from `wwn_main_run_all_outputs.R` on March 24, 2026.
- The current short MODIS tables include both `ClimateComposite` and `ENSO_NINO34`.
- Best-lag species heatmaps now show the winning lag directly in each tile as `L0`, `L1`, or `L2`.
- `enso_annual_nino34.csv` is the active ENSO source used by the current workflow.
- `dated8.rwl` still triggers parser fallback warnings in `dplR`, but the full workflow completes successfully and writes outputs.

## Repository Purpose

This repository is intended as a reproducible record of the thesis analysis, not as a polished R package. It is designed so that another student or committee member can inspect the scripts, rerun the workflow, and recover the major figures and tables used in the thesis.
