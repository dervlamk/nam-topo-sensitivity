# Data Manifest

Raw input data files required to replicate the analyses in this repository,
**grouped by experiment** (e.g. FLOR ctrl, FLOR hitopo, FLOR hicam) rather than
by variable type.

All files live under a single top-level directory:

```
dpath0 = /discover/nobackup/projects/giss/baldwin_nip/dmkumar
```

Paths below are relative to `dpath0`. Files **written** by the notebooks
(caches in `../highresmip_diagnostics/`, figure outputs in `figs/`) are not
raw inputs and are not listed here. A flat list of full paths is also provided
in `DATA_MANIFEST.csv` (the `group` column matches the section headings below).

## Observations

### Observed precipitation — `obs_data/prec/`

| File | Product | Used by |
|---|---|---|
| `imerg.gn.timeseries.2001-2018.nc` | IMERG | flor_pi, flor_2xCO2, HRMIP historical, HRMIP future |
| `pr_TRMM-L3_v7-7A_199801-201312.nc` | TRMM | flor_pi, HRMIP historical, HRMIP future |
| `gpcp.precip.1979-2018.monthly.nc` | GPCP | flor_pi, HRMIP historical |
| `gpcc.precip.mon.rate.0.25x0.25.v2020.nc` | GPCC | flor_pi, HRMIP historical |

*(`chirps-v3.0.monthly.nc` appears commented-out in `nam_precip_flor_pi.ipynb` and is not actually used.)*

### Observed reanalysis (MERRA-2) — `obs_data/merra2/`

| File | Variable | Used by |
|---|---|---|
| `merra2.U.1980-2022.monthly.nc` | U winds (lev 700, 925) | 700mb_psi, slp_gocllj |
| `merra2.V.1980-2022.monthly.nc` | V winds (lev 700, 925) | 700mb_psi, slp_gocllj |
| `merra2.H.1980-2022.monthly.nc` | Geopotential height (lev 500, 700) | 500mb_geopotential, 700mb_psi |
| `merra2.SLP.1980-2022.monthly.nc` | Sea level pressure | slp_gocllj |

### Observed topography — `topo_files/`

| File | Source | Used by |
|---|---|---|
| `obs.etopo5.zsurf.nc` | ETOPO5 observed | all FLOR notebooks, HRMIP, topography |

> **Grid registration — checked, not the OIPC bug.** `topography_obs_models_comparisons.ipynb`
> plots `obs.etopo5.zsurf.nc`'s topography against `ax.coastlines()` (Natural Earth vector
> data), and it initially looked offset — the same symptom documented for the OIPC isoscape
> in `nam-dD-lig/analysis-repo/DATA_MANIFEST.md` under "OIPC grid registration". Measured
> with `find_coastline_offset()` (`py_functions/data_funcs.py`), which sweeps a rigid lon/lat
> shift and scores agreement between the field's land mask and Natural Earth polygons,
> independently over four coastal windows:
>
> | window | agreement | shift lon | shift lat |
> |---|---|---|---|
> | Baja Pacific coast | 0.983 | +0.0000 | −0.0000 |
> | Gulf of California | 0.962 | +0.0000 | −0.0000 |
> | Gulf of Mexico | 0.990 | +0.1667 | −0.0000 |
> | Pacific NW | 0.991 | +0.0000 | −0.0000 |
>
> Three of four windows land on zero shift with 96–99% agreement — correctly-registered
> data. Gulf of Mexico's +0.1667° is a lone outlier (that coast's low-lying deltas/barrier
> islands are hard for ETOPO5 to resolve well); a real registration bug would show the same
> shift everywhere, so this isn't one. **`ETOPO05_X`/`ETOPO05_Y` are not offset — no
> coordinate correction applied.**
>
> The actual cause: this notebook's `ax.pcolormesh(field.lon, field.lat, field, ...)` calls
> didn't pass `shading`. With coordinate arrays matching the data array's shape exactly,
> matplotlib's legacy default drops the last row/column of data and treats the given lon/lat
> as cell edges instead of centers — shifting the mesh by about half a grid cell (~0.04° at
> ETOPO05's 5′ resolution), enough to look like a coastline mismatch on a regional map. Fixed
> 2026-08-24 by adding `shading='auto'` to every topo `pcolormesh()` call in the notebook.

## FLOR experiments — `FLOR/{run}/{case}/` and `topo_files/`

Each FLOR experiment is grouped together below with all of its variables and its
ZSURF topography boundary condition. Cases: `pi` (pre-industrial), `2xco2`.

### FLOR ctrl

| File | Variable | Used by |
|---|---|---|
| `FLOR/ctrl/pi/flor.ctrl.precip.monthly.nc` | Precip (PI) | flor_pi, flor_2xCO2 |
| `FLOR/ctrl/2xco2/flor.ctrl.2xco2.precip.monthly.nc` | Precip (2×CO₂) | flor_2xCO2 |
| `FLOR/ctrl/pi/flor.ctrl.h500.monthly.nc` | 500 hPa height | 500mb_geopotential |
| `FLOR/ctrl/pi/flor.ctrl.h700.monthly.nc` | 700 hPa height | 700mb_psi |
| `FLOR/ctrl/pi/flor.ctrl.u700.monthly.nc` | 700 hPa U wind | 700mb_psi |
| `FLOR/ctrl/pi/flor.ctrl.v700.monthly.nc` | 700 hPa V wind | 700mb_psi |
| `FLOR/ctrl/pi/flor.ctrl.u925.monthly.nc` | 925 hPa U wind | slp_gocllj |
| `FLOR/ctrl/pi/flor.ctrl.v925.monthly.nc` | 925 hPa V wind | slp_gocllj |
| `FLOR/ctrl/pi/flor.ctrl.slp.monthly.nc` | Sea level pressure | slp_gocllj |
| `topo_files/flor.ctrl.zsurf.nc` | ZSURF topography | all FLOR notebooks, topography |

### FLOR hitopo

| File | Variable | Used by |
|---|---|---|
| `FLOR/hitopo/pi/flor.hitopo.precip.monthly.nc` | Precip (PI) | flor_pi, flor_2xCO2 |
| `FLOR/hitopo/2xco2/flor.hitopo.2xco2.precip.monthly.nc` | Precip (2×CO₂) | flor_2xCO2 |
| `FLOR/hitopo/pi/flor.hitopo.h500.monthly.nc` | 500 hPa height | 500mb_geopotential |
| `FLOR/hitopo/pi/flor.hitopo.h700.monthly.nc` | 700 hPa height | 700mb_psi |
| `FLOR/hitopo/pi/flor.hitopo.u700.monthly.nc` | 700 hPa U wind | 700mb_psi |
| `FLOR/hitopo/pi/flor.hitopo.v700.monthly.nc` | 700 hPa V wind | 700mb_psi |
| `FLOR/hitopo/pi/flor.hitopo.u925.monthly.nc` | 925 hPa U wind | slp_gocllj |
| `FLOR/hitopo/pi/flor.hitopo.v925.monthly.nc` | 925 hPa V wind | slp_gocllj |
| `FLOR/hitopo/pi/flor.hitopo.slp.monthly.nc` | Sea level pressure | slp_gocllj |
| `topo_files/flor.hitopo.zsurf.nc` | ZSURF topography | all FLOR notebooks, topography |

### FLOR hicam

| File | Variable | Used by |
|---|---|---|
| `FLOR/hicam/pi/flor.hicam.precip.monthly.nc` | Precip (PI) | flor_pi |
| `FLOR/hicam/pi/flor.hicam.h500.monthly.nc` | 500 hPa height | 500mb_geopotential |
| `FLOR/hicam/pi/flor.hicam.h700.monthly.nc` | 700 hPa height | 700mb_psi |
| `FLOR/hicam/pi/flor.hicam.u700.monthly.nc` | 700 hPa U wind | 700mb_psi |
| `FLOR/hicam/pi/flor.hicam.v700.monthly.nc` | 700 hPa V wind | 700mb_psi |
| `FLOR/hicam/pi/flor.hicam.u925.monthly.nc` | 925 hPa U wind | slp_gocllj |
| `FLOR/hicam/pi/flor.hicam.v925.monthly.nc` | 925 hPa V wind | slp_gocllj |
| `FLOR/hicam/pi/flor.hicam.slp.monthly.nc` | Sea level pressure | slp_gocllj |
| `topo_files/flor.hicam.zsurf.nc` | ZSURF topography | flor_pi, 700mb_psi, slp_gocllj, topography |

> **Verify:** `500mb_geopotential_flor_pi.ipynb` references `flor.cam.zsurf.nc`,
> while every other notebook uses `flor.hicam.zsurf.nc`. These may be the same
> file under two names, or one may be a typo — confirm before publishing.

### FLOR common

| File | Variable | Used by |
|---|---|---|
| `FLOR/flor.land_mask.nc` | Land mask (shared by all runs) | topography_comparisons |

## HighResMIP precipitation — `hrmip/{case}/pr/`

File pattern: `hrmip/{case}/pr/pr.Amon.{case}.{cfg}.nc`. Only precipitation (`pr`)
is used, so these are grouped by HighResMIP experiment (`highresSST-present`,
`highresSST-future`).

### highresSST-present (historical)

Used by `HighResMIP_prec_historical_biases.ipynb` (all configs) and, for the
7-model subset, by `HighResMIP_prec_future_delta.ipynb`.

**Low-resolution configs:**
`CAM-MPAS-LR`, `CAMS-CSM1-0-LR`, `CMCC-CM2-HR4`, `CNRM-CM6-1`, `ECMWF-IFS-LR`,
`FGOALS-f3-L`, `GFDL-CM4`, `HadGEM3-GC31-LM`, `INM-CM5-0`, `IPSL-CM6A-LR`,
`MPI-ESM1-2-LR`, `MRI-ESM2-0`

**High-resolution configs:**
`CAM-MPAS-HR`, `CAMS-CSM1-0`, `CMCC-CM2-VHR4`, `CNRM-CM6-1-HR`, `EC-Earth3P`,
`EC-Earth3P-HR`, `ECMWF-IFS-HR`, `FGOALS-f3-H`, `GFDL-CM4C192`,
`HadGEM3-GC31-MM`, `HadGEM3-GC31-HM`, `HiRAM-SIT-LR`, `HiRAM-SIT-HR`,
`INM-CM5-H`, `IPSL-CM6A-ATM-HR`, `MPI-ESM1-2-HR`, `MPI-ESM1-2-XR`,
`MRI-AGCM3-2-H`, `MRI-AGCM3-2-S`, `NICAM16-7S`, `NICAM16-8S`

### highresSST-future

Used by `HighResMIP_prec_future_delta.ipynb` (7-model LR/HR paired subset):

| LR config | HR config |
|---|---|
| `CAM-MPAS-LR` | `CAM-MPAS-HR` |
| `CMCC-CM2-HR4` | `CMCC-CM2-VHR4` |
| `EC-Earth3P` | `EC-Earth3P-HR` |
| `FGOALS-f3-L` | `FGOALS-f3-H` |
| `HadGEM3-GC31-LM` | `HadGEM3-GC31-HM` |
| `HiRAM-SIT-LR` | `HiRAM-SIT-HR` |
| `MPI-ESM1-2-HR` | `MPI-ESM1-2-XR` |

> The future-delta notebook also reads the matching `highresSST-present` files
> for these same 7 model pairs (already listed above), so there are no
> additional present-day files beyond those listed.
