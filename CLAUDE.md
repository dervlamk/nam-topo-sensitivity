# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Research codebase for studying the sensitivity of North American Monsoon (NAM) precipitation to topographic boundary conditions. Analyzes JAS (July–August–September) precipitation across FLOR topography sensitivity experiments and HighResMIP (HRMIP) historical and future simulations.

## Environment

This project runs on NASA's NCCS Discover supercomputer. Notebooks should be run from the `notebooks/` directory. The Python environment includes xarray, numpy, matplotlib, cartopy, metpy, scipy, regionmask, pyproj, shapely, cmocean, seaborn, cf_xarray, cftime, and netCDF4.

Cartopy requires a pre-existing data directory; all notebooks set:
```python
cartopy.config['data_dir'] = "/discover/nobackup/projects/jh_tutorials/JH_examples/JH_datafiles/Cartopy"
cartopy.config['pre_existing_data_dir'] = "/discover/nobackup/projects/jh_tutorials/JH_examples/JH_datafiles/Cartopy"
```

## Data Paths

All external data lives under:
```
dpath0 = '/discover/nobackup/projects/giss/baldwin_nip/dmkumar'
```

Key sub-paths:
- Obs precip: `{dpath0}/obs_data/prec/` (IMERG, TRMM, GPCP, GPCC)
- FLOR model output: `{dpath0}/FLOR/{run}/{case}/flor.{run}.precip.monthly.nc`
  - Runs: `ctrl`, `hicam` (HI_cam), `hitopo` (HI_gbl)
  - Cases: `pi` (pre-industrial), `2xco2`
- Topography files: `{dpath0}/topo_files/` (ETOPO5 obs; FLOR ctrl/cam/hitopo ZSURFs)
- HighResMIP data: `{dpath0}/hrmip/highresSST-present/{var}/{var}.Amon.highresSST-present.{cfg}.nc`

Pre-processed HRMIP climatology files are stored locally in `hrmip/` (excluded from git via `*.nc` in `.gitignore`).

## Custom Functions

Notebooks import shared utilities from a `py_functions/` directory. The import path varies by notebook:
- Most notebooks (one level deep in `notebooks/`): `../py_functions/`
- The main figure notebook `hrmip_historical_prec_biases_figs_2_s1.ipynb` uses `../../py_functions/`

Key modules:
- `map_plot_tools.py` — `quick_map()` and other cartopy-based map helpers
- `colorbar_funcs.py` — `get_settings(field, diff)` returns colormap, levels, norms
- `data_funcs.py` — `lonFlip()` for converting longitude conventions (0:360 ↔ -180:180)
- `stats_funcs.py` — `sigtest()` and `sigtest2n()` for Student's t-test significance masking

## Notebook Architecture

All analysis is in `notebooks/`. Figure outputs go to `figs/`.

| Notebook | Purpose |
|---|---|
| `hrmip_historical_prec_biases_figs_2_s1.ipynb` | Main paper figures: HRMIP historical JAS precip biases vs. IMERG, including transect profiles and MME maps (Figs. 2 & S1) |
| `nam_precip_flor_pi.ipynb` | FLOR pre-industrial JAS precip bias and topography sensitivity |
| `nam_precip_flor_2xCO2.ipynb` | FLOR 2xCO2 future precip change; ctrl vs. hitopo |
| `nam_precip_hrmip_future.ipynb` | HRMIP future (2030–2050, 2040–2050) precip change |
| `nam_synoptic.ipynb` | Synoptic-scale analysis |
| `nam_winds.ipynb` | Wind field analysis |
| `slp_flor_pi.ipynb` | Sea level pressure (FLOR PI) |
| `slp_llj_flor.ipynb` | SLP and low-level jet analysis |
| `qs_flor.ipynb` | Specific humidity (FLOR) |
| `h850_ws850_flor.ipynb` | 850 hPa heights and wind speed |

### Common Analysis Pattern

All notebooks follow this structure:
1. Import packages + custom functions
2. Define data paths and load/pre-process datasets (units conversion: `* 86400` for mm/s→mm/day, `* 24` for mm/hr→mm/day; `lonFlip()` for FLOR output)
3. Compute JAS seasonal means — use `jas_seasonal_mean()` (weighted by `days_in_month`) and `jas_yearly_mean()` (per-year means for significance testing)
4. Regrid obs to model grid via `xr.DataArray.interp()`
5. Significance testing via `sigtest2n()` or `sigtest()`
6. Produce matplotlib/cartopy figures, save to `figs/` as both `.pdf` and `.png`

### HRMIP Model Configuration Maps

The main HRMIP notebook uses two dicts, `lr_map` and `hr_map`, keyed by modeling group name (15 groups), mapping to lists of low-res and high-res configuration IDs respectively. Model data is stored in nested dicts `dat['lr'][cfg]` and `dat['hr'][cfg]`.

### NAM Sub-domains

Two analysis domains are defined via `regionmask`:
- **core**: rotated box centered at (-107°, 25°), 4.5°×10°, 35° tilt
- **north**: box centered at (-110.5°, 32.5°), 8°×6°, no tilt
