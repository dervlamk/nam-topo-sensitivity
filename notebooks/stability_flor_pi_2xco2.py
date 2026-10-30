#!/usr/bin/env python
"""Convective stability (buoyancy index b = h_10m - h*) in FLOR: PI and 2xCO2.

Script version of notebooks/stability_flor_pi_2xco2.ipynb. FLOR ctrl and hitopo, PI and 2xco2 (see `runs`/`cases` below).
Run from any directory: python notebooks/stability_flor_pi_2xco2.py   (figures go to ../figs relative to notebooks/)

# Convective stability (buoyancy index) in FLOR: PI and 2xCO$_2$

Buoyancy index of a near-surface parcel lifted (moist) adiabatically, following the methods of the paper:

$$b = h_{10\,\mathrm{m}} - h^{*}, \qquad h = c_p T + g z + L q, \qquad h_{10\,\mathrm{m}} = c_p T_{10\,\mathrm{m}} + g z_{10\,\mathrm{m}} + L q_{10\,\mathrm{m}}$$

$h^*$ is the saturation moist static energy of the environment ($q \to q^*(T,p)$, August-Roche-Magnus). Since $c_p(T_c-T) = (h_c-h^*)/(1+\gamma)$ with $\gamma$ positive and O(1), $b$ is roughly twice the parcel buoyancy; $b>0$ means upward acceleration. $b$ is evaluated on pressure levels above the lifted condensation level (LCL) and we look at the anomaly $\Delta b$ between a perturbed and a control simulation.

**Status of this iteration:** FLOR `ctrl` and `hitopo` (HIgbl), `pi` and `2xco2` cases (`runs` / `cases`). Setting `cases = ['pi']` or `runs = ['ctrl']` runs a subset; the 2xCO2 and hitopo steps are gated on `RUN_2XCO2` / `RUN_HITOPO`.

## Input files

All under `dpath0/FLOR/{run}/{case}/`; `pi` files are `flor.{run}.{var}.monthly.nc`, `2xco2` files are `flor.{run}.2xco2.{var}.monthly.nc`.
**File and variable names for everything marked "confirm" are assumptions; check them with `ls` / `ncdump -h` on Discover.**

| Input | Variable | Used for | Status |
|---|---|---|---|
| `temp` (full pressure levels) | $T(p)$ | $h^*$ | new, confirm |
| `z_full` (full pressure levels) | $z(p)$ | $g z$ term of $h^*$ | new, confirm |
| `t_ref` | $T_{10\,\mathrm{m}}$ (model reference level) | $h_{10\,\mathrm{m}}$, LCL | exists (used by `warming_rate.py`), path unrecorded |
| `q_ref` | $q_{10\,\mathrm{m}}$ | $h_{10\,\mathrm{m}}$, LCL dewpoint | new, confirm |
| `ps` | surface pressure | LCL pressure | new, confirm |
| `topo_files/flor.{run}.zsurf.nc` | `ZSURF` | $z_{10\,\mathrm{m}}$ = ZSURF + 10 m | exists |
| RAS convective mass flux (optional) | e.g. `mc_full` | convection intensity (last section) | new, confirm |

`sphum` (3-D specific humidity) is **not** needed: $h^*$ uses the saturation humidity computed from $T$ and $p$.

## Constants

| Constant | Value | Source |
|---|---|---|
| $c_p$ | 1004 J K$^{-1}$ kg$^{-1}$ | methods |
| $L$ | 2.5$\times$10$^6$ J kg$^{-1}$ | methods |
| $z_{10\,\mathrm{m}}$ | 10 m above the surface | methods (model reference height) |
| $g$ | 9.80665 m s$^{-2}$ | assumed (methods give no value) |
| Magnus $A, B, C$ | 6.1094 hPa, 17.625, 243.04 $^\circ$C | assumed standard August-Roche-Magnus values |
| $\epsilon = R_d/R_v$ | 0.622 | assumed |
| $\kappa = R_d/c_p$ | 0.2854 | Bolton (1980), for the LCL pressure |
| LCL temperature constants | 56 K, 800 | Bolton (1980) Eq. 15 |
"""

import os
import sys
import warnings
warnings.filterwarnings('ignore')
warnings.simplefilter('ignore')

import xarray as xr
xr.set_options(keep_attrs=True)
import numpy as np
np.seterr(divide='ignore', invalid='ignore')

import cartopy
cartopy.config['data_dir'] = "/discover/nobackup/projects/jh_tutorials/JH_examples/JH_datafiles/Cartopy"
cartopy.config['pre_existing_data_dir'] = "/discover/nobackup/projects/jh_tutorials/JH_examples/JH_datafiles/Cartopy"
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER

import matplotlib as mpl
mpl.use('Agg')  # non-interactive backend: figures are only saved to file
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.ticker as mticker

# settings

# run relative to notebooks/ so '../py_functions' and '../figs' resolve as in the notebooks
os.chdir(os.path.dirname(os.path.abspath(__file__)))

# add path to custom functions
module_path = os.path.abspath(os.path.join('..'))
if module_path not in sys.path:
    sys.path.append(module_path+"/py_functions")
from colorbar_funcs import *
from data_funcs import *
from stats_funcs import *
from domain_funcs import *

dpath0='/discover/nobackup/projects/giss/baldwin_nip/dmkumar' # top level data directory

# ## Constants and thermodynamic functions

### +++ CONSTANTS +++ ###

CP    = 1004.0      # J K-1 kg-1, isobaric specific heat of dry air  (methods)
LV    = 2.5e6       # J kg-1, latent heat of condensation            (methods)
Z_REF = 10.0        # m, model reference height of t_ref / q_ref     (methods)
G     = 9.80665     # m s-2, gravitational acceleration              (assumed)
EPS   = 0.622       # Rd/Rv                                          (assumed)
KAPPA = 0.2854      # Rd/cp, Bolton (1980); LCL pressure             (Bolton 1980)
# August-Roche-Magnus saturation vapor pressure: es[hPa] = MAG_A * exp(MAG_B*Tc / (Tc + MAG_C)), Tc in degC
MAG_A, MAG_B, MAG_C = 6.1094, 17.625, 243.04   # hPa, -, degC          (assumed)
T0 = 273.15         # K, 0 degC

for name in ['CP','LV','Z_REF','G','EPS','KAPPA','MAG_A','MAG_B','MAG_C']:
    print(f'{name:6s} = {globals()[name]}')

### +++ THERMODYNAMIC FUNCTIONS +++ ###
# Inputs: T in K, p in Pa, q in kg/kg, z in m. Work on floats, numpy arrays and xarray objects.

def esat(T):
    """August-Roche-Magnus saturation vapor pressure [Pa] over liquid water. T in K."""
    Tc = T - T0
    return 100.0 * MAG_A * np.exp(MAG_B * Tc / (Tc + MAG_C))


def qsat(T, p):
    """Saturation specific humidity [kg/kg] from Magnus es(T). T in K, p in Pa."""
    es = esat(T)
    return EPS * es / (p - (1.0 - EPS) * es)


def gamma(T, p):
    """gamma = (L/cp) * (dq*/dT)_p, from the analytic derivative of qsat. Not used in b, which is h10m - h* directly
    (b ~ (1+gamma) * cp * (Tc - T))."""
    Tc = T - T0
    es = esat(T)
    des_dT = es * MAG_B * MAG_C / (Tc + MAG_C)**2
    dq_des = EPS * p / (p - (1.0 - EPS) * es)**2
    return (LV / CP) * dq_des * des_dT


def mse(T, z, q):
    """Moist static energy h = cp*T + g*z + L*q  [J/kg]."""
    return CP * T + G * z + LV * q


def mse_sat(T, z, p):
    """Saturation moist static energy h* = cp*T + g*z + L*q*(T,p)  [J/kg]."""
    return mse(T, z, qsat(T, p))


def dewpoint(q, p):
    """Dewpoint [K] from specific humidity q [kg/kg] and pressure p [Pa] (inverting the Magnus formula)."""
    e = q * p / (EPS + (1.0 - EPS) * q)          # vapor pressure [Pa]
    x = np.log(e / (100.0 * MAG_A))
    return MAG_C * x / (MAG_B - x) + T0


def lcl_temperature(T, Td):
    """LCL temperature [K], Bolton (1980) Eq. 15. T and Td in K."""
    return 1.0 / (1.0 / (Td - 56.0) + np.log(T / Td) / 800.0) + 56.0


def lcl_pressure(t_ref, q_ref, ps):
    """LCL pressure [Pa] of near-surface air: dewpoint from (q_ref, ps), T_L from Bolton Eq. 15,
    then dry-adiabatic (Poisson) ascent p_L = ps * (T_L / t_ref)**(1/KAPPA)."""
    Td = np.minimum(dewpoint(q_ref, ps), t_ref)  # guard against supersaturated reference values
    TL = lcl_temperature(t_ref, Td)
    return ps * (TL / t_ref)**(1.0 / KAPPA)


def buoyancy_index(t_ref, q_ref, zsurf, temp, zg):
    """Buoyancy index b = h_10m - h* [J/kg] on the pressure levels of `temp`.

    Parameters
    ----------
    t_ref, q_ref : near-surface (reference-height) temperature [K] and specific humidity [kg/kg], dims (time, lat, lon)
    zsurf        : surface height [m], dims (lat, lon); z_10m = zsurf + Z_REF
    temp, zg     : environment temperature [K] and height [m], dims (time, plev, lat, lon), plev in Pa

    Returns b on all levels with valid data; levels below the LCL / ground are masked later (see `below_lcl`),
    because masking monthly values would bias the days-in-month weighted JAS means.
    """
    h10 = mse(t_ref, zsurf + Z_REF, q_ref)
    hstar = mse_sat(temp, zg, temp['plev'])
    return (h10 - hstar).transpose('time', 'plev', 'lat', 'lon')


def below_lcl(plev, plcl):
    """True where a pressure level lies below (at higher pressure than) the LCL."""
    return plev > plcl

### +++ SANITY CHECKS +++ ###

assert 0.0215 < qsat(300., 1e5) < 0.0230,            qsat(300., 1e5)           # ~0.0222 kg/kg
assert 1.0 < gamma(300., 1e5) < 6.0,                 gamma(300., 1e5)          # positive, O(1)
assert abs(dewpoint(qsat(290., 9e4), 9e4) - 290.) < 1e-3                       # dewpoint of saturated air = T
assert abs(lcl_temperature(300., 290.) - 287.7) < 0.3, lcl_temperature(300., 290.)  # LCL colder than Td
assert abs(lcl_temperature(290., 290.) - 290.) < 1e-6                          # saturated: LCL at surface
assert abs(lcl_pressure(300., qsat(300., 1e5), 1e5) - 1e5) < 1.0               # saturated: p_LCL = ps
assert 7e4 < lcl_pressure(300., 0.012, 1e5) < 9.5e4                            # dry-ish: LCL well above the surface

# saturated parcel at the same T and z as the environment has b = 0
_p  = xr.DataArray([1e5], dims='plev', coords={'plev':[1e5]})
_t  = xr.DataArray(np.full((1,1,1,1), 300.), dims=('time','plev','lat','lon'), coords={'plev':[1e5]})
_z  = xr.DataArray(np.full((1,1,1,1), 1000.), dims=('time','plev','lat','lon'), coords={'plev':[1e5]})
_tr = xr.DataArray(np.full((1,1,1), 300.), dims=('time','lat','lon'))
_qr = qsat(_tr, 1e5)
_zs = xr.DataArray(np.full((1,1), 1000.-Z_REF), dims=('lat','lon'))
assert float(np.abs(buoyancy_index(_tr, _qr, _zs, _t, _z)).max()) < 1e-6   # .max() because the result is a (1,1,1,1) array
print('All checks passed.')

# ## Data paths

### +++ DATA PATHS +++ ###

cases=['pi','2xco2']
RUN_2XCO2 = '2xco2' in cases   # gates every 2xCO2-dependent step below
RUN_HITOPO = True              # also run FLOR hitopo (HIgbl) and the HIgbl-CTRL comparison
runs=['ctrl','hitopo'] if RUN_HITOPO else ['ctrl']
RUN_LABEL = {'ctrl':'CTRL', 'hitopo':r'HI$_{gbl}$'}   # plot labels

vars_3d=['temp','z_full']  # full pressure-level fields, CONFIRM file/variable names
vars_2d=['t_ref','q_ref','ps']   # CONFIRM file/variable names (t_ref: see warming_rate.py)
PLEV_NAME = 'level'        # pressure-level dimension in the 3-D files (hPa: 1000,925,850,700,600,500,400,300,250,200,150,100,70,50,30,20,10)

files={ 'ctrl':   { 'pi':{}, '2xco2':{} },
        'hitopo': { 'pi':{}, '2xco2':{} }
      }

Q_REF_YEARS = {'pi':('0101','0200'), '2xco2':('0200','0299')}   # periods of the year-tagged q_ref files made by q_ref_from_rh.py

def flor_path(key, case, varn):
    """Monthly FLOR file. q_ref: the year-tagged file from q_ref_from_rh.py if it exists, else the untagged name."""
    tag = f'flor.{key}' if case=='pi' else f'flor.{key}.{case}'
    if varn=='q_ref':
        tagged = f'{dpath0}/FLOR/{key}/{case}/{tag}.q_ref.{"-".join(Q_REF_YEARS[case])}.monthly.nc'
        if os.path.exists(tagged):
            return tagged
    return f'{dpath0}/FLOR/{key}/{case}/{tag}.{varn}.monthly.nc'

for key in runs:
    for case in cases:
        for varn in vars_3d+vars_2d:
            files[key][case][varn] = flor_path(key, case, varn)
        print(f'{key} {case} q_ref file: {files[key][case]["q_ref"]}')

files['topo'] = {}
for run in runs:
    files['topo'][run] = f'{dpath0}/topo_files/flor.{run}.zsurf.nc'

# ## Organize data

## +++ ORGANIZE DATA +++ ##

# --- sub-setting bounds --- #
latmin, latmax = 5, 70
lonmin, lonmax = 200, 300
# years and figure directory can be overridden from the environment (used by sbatch_stability_flor_debug.sh)
yr_min, yr_max = os.environ.get('STAB_YR_MIN', '0101'), os.environ.get('STAB_YR_MAX', '0200') # PI years
yr_min_2x, yr_max_2x = '0200', '0299' # 2xCO2 years (CO2 held at 2x after doubling at ~yr 170)
FIG_DIR = os.environ.get('STAB_FIG_DIR', '../figs')
os.makedirs(FIG_DIR, exist_ok=True)
plev_min = 10000.  # Pa, drop levels above 100 hPa -> this is only to save compute requirements, not done by Pascale et al. originally

def load_var(path, varn, case, is_3d=False):
    """Open one FLOR variable, select the case's years and JAS months (the only months used downstream),
    rename to lon/lat/plev (plev in Pa), and subset to the analysis box."""
    ds = xr.open_dataset(path, chunks={'time': 12})[varn]
    yrs = (yr_min_2x, yr_max_2x) if case=='2xco2' else (yr_min, yr_max)
    ds = ds.sel(time=slice(*yrs))
    ds = ds.sel(time=ds['time.month'].isin([7,8,9]))
    ds = ds.rename({'grid_xt':'lon','grid_yt':'lat'})
    if is_3d:
        ds = ds.rename({PLEV_NAME:'plev'})
        if float(ds.plev.max()) < 2000.:       # levels given in hPa -> Pa
            ds = ds.assign_coords(plev=ds.plev*100.)
        ds = ds.sel(plev=slice(plev_min, None)) if float(ds.plev[0]) < float(ds.plev[-1]) else ds.sel(plev=slice(None, plev_min))
        ds = ds.transpose('time','plev','lat','lon')
    else:
        ds = ds.transpose('time','lat','lon')
    return ds.sel(lon=slice(lonmin,lonmax), lat=slice(latmin,latmax)).astype('float32')

dat={ 'ctrl':   { 'pi':{}, '2xco2':{} },
      'hitopo': { 'pi':{}, '2xco2':{} }
    }

print('Working on...')
for key in runs:
    print(f'{key}')
    for case in cases:
        for varn in vars_3d:
            dat[key][case][varn] = load_var(files[key][case][varn], varn, case, is_3d=True)
        for varn in vars_2d:
            dat[key][case][varn] = load_var(files[key][case][varn], varn, case)

topo = {}
for run in runs:
    ds = xr.open_dataset(files['topo'][run]).ZSURF.rename({'GRID_XT':'lon', 'GRID_YT':'lat'})
    topo[run] = ds.clip(min=0).sel(lat=slice(latmin,latmax), lon=slice(lonmin,lonmax)).astype('float32')  # m, ocean -> 0
    del ds

print('Done.')

# ## Buoyancy index and JAS means
#
# $b$ is computed **monthly** (the saturation humidity is nonlinear in $T$) and then averaged with the repo's days-in-month weighted JAS means. The levels below the LCL are masked *after* averaging, using the JAS-mean LCL pressure; masking monthly values would bias the weighted means. Because ground-level points lie below the LCL, the same mask removes below-ground levels. For $\Delta b$ the mask is the union of the two cases so both fields are defined on the same points.

### +++ BUOYANCY INDEX AND SEASONAL MEANS +++ ###

b_jas     = { key:{ 'pi':{}, '2xco2':{} } for key in runs }
b_ann     = { key:{ 'pi':{}, '2xco2':{} } for key in runs }
plcl_jas  = { key:{ 'pi':{}, '2xco2':{} } for key in runs }

print('Calculating b for...')
for key in runs:
    print(f'{key}')
    for case in cases:
        d = dat[key][case]
        b  = buoyancy_index(d['t_ref'], d['q_ref'], topo[key], d['temp'], d['z_full'])
        pl = lcl_pressure(d['t_ref'], d['q_ref'], d['ps'])
        b_jas[key][case]    = jas_seasonal_mean(b).load()
        b_ann[key][case]    = jas_yearly_mean(b).load().transpose('year','plev','lat','lon')
        plcl_jas[key][case] = jas_seasonal_mean(pl).load()
        del b, pl

# mask levels below the LCL: union over every run and case, so that all fields (and the differences between
# them: 2xCO2-PI, HIgbl-CTRL) are defined on common points
plev = b_jas['ctrl']['pi'].plev
mask = None
for key in runs:
    for case in cases:
        m = below_lcl(plev, plcl_jas[key][case])
        mask = m if mask is None else (mask | m)
for key in runs:
    for case in cases:
        b_jas[key][case] = b_jas[key][case].where(~mask)
        b_ann[key][case] = b_ann[key][case].where(~mask)
print('Done.')

# ## Significance testing: 2xCO$_2$ $-$ PI

## +++ COMPARING 2xCO2 TO PI +++ ##

if RUN_2XCO2:   # needs the 2xCO2 case
    future_diff      = {}
    future_diff_mask = {}
    future_ptvals    = {}

    print('Significance testing for:')
    for key in runs:
        print(f'{key}')
        diff_, diff_mask_, ptvals_ = sigtest2n(b_ann[key]['2xco2'], b_ann[key]['pi'],
                                               b_jas[key]['2xco2'], b_jas[key]['pi'])
        future_diff[key] = diff_
        future_diff_mask[key] = diff_mask_
        future_ptvals[key] = ptvals_

# --- HIgbl - CTRL (the LCL mask above is the union over ctrl AND hitopo, so both are on common points) --- #
if RUN_HITOPO:
    topo_diff, topo_diff_mask, topo_ptvals = {}, {}, {}
    for case in cases:
        topo_diff[case], topo_diff_mask[case], topo_ptvals[case] = sigtest2n(
            b_ann['hitopo'][case], b_ann['ctrl'][case], b_jas['hitopo'][case], b_jas['ctrl'][case])
print('Done.')

coords_map, regions = nam_regions()
domains = list(coords_map.keys())

# ## FIGS

# ### JAS buoyancy index $b$ in PI (CTRL) at selected levels

# --- Settings --- #
text_kw={'color':'k', 'weight':'bold', 'size':24, 'horizontalalignment':'center', 'verticalalignment':'bottom'}
text_kw3={'color':'k', 'weight':'bold', 'size':30, 'horizontalalignment':'center', 'verticalalignment':'bottom'}
plot_levels=[850, 700, 500]       # hPa
letters=['a','b','c','d','e','f']
tx=-109
ty=38.1
cmap=plt.get_cmap('RdBu_r')       # red = more unstable (b>0)
bvmin, bvmax = -10, 10            # kJ/kg, JAS-mean b
dvmin, dvmax = -2, 2              # kJ/kg, Delta b
trans=ccrs.PlateCarree()
proj=ccrs.PlateCarree()
map_bnds=[243,259,18,38]
step=1                            # stippling stride

def at_level(da, hpa):
    """Select the model level nearest to `hpa` and convert J/kg -> kJ/kg."""
    return da.sel(plev=hpa*100., method='nearest') / 1000.

def level_maps(fields, masks, titles, vmin, vmax, nlev, label, extend='both'):
    """Row of maps (one per entry in `fields`). `masks` (or None) are significance masks: stipple where NaN."""
    norm=mpl.colors.BoundaryNorm(np.linspace(vmin, vmax, nlev), cmap.N)
    lon=fields[0].lon; lat=fields[0].lat
    lon2d, lat2d = np.meshgrid(lon, lat)
    fig, axs = plt.subplots(nrows=1, ncols=len(fields), figsize=(7*len(fields),10), layout='constrained', subplot_kw={'projection':proj})
    for i, ax in enumerate(np.atleast_1d(axs)):
        cf=ax.pcolormesh(lon, lat, fields[i], cmap=cmap, norm=norm, transform=trans)
        if masks is not None:
            sig_mask = masks[i].isnull().values
            ax.scatter(lon2d[sig_mask][::step], lat2d[sig_mask][::step], s=30, c='k', marker='.', linewidths=0, transform=trans, zorder=100)
        ax.text(tx, ty, titles[i], **text_kw)
        ax.text(-117, 38.25, letters[i], **text_kw3)
        poly1=patches.Polygon(coords_map['south'], closed=True, ec='red', fc='none', lw=3, ls='--', transform=trans, zorder=100); ax.add_patch(poly1)
        poly2=patches.Polygon(coords_map['north'], closed=True, ec='firebrick', fc='none', lw=3, ls='--', transform=trans, zorder=100); ax.add_patch(poly2)
        ax.coastlines(color='k', lw=1)
        ax.set_extent(map_bnds, crs=trans)
        gl=ax.gridlines(crs=trans, lw=.5, colors='black', alpha=0.5, linestyle='--', zorder=10, draw_labels=True)
        gl.xlocator = mticker.MultipleLocator(4)
        gl.ylocator = mticker.MultipleLocator(3)
        gl.bottom_labels=True; gl.left_labels=(i==0); gl.top_labels=False; gl.right_labels=False
        gl.xformatter=LONGITUDE_FORMATTER
        gl.yformatter=LATITUDE_FORMATTER
        gl.xlabel_style={'color': 'black', 'weight': 'normal', 'size':16}
        gl.ylabel_style={'color': 'black', 'weight': 'normal', 'size':16}
    cax=fig.add_axes([1.01, .1, 0.02, 0.8])
    cbar=fig.colorbar(cf, ticks=np.linspace(vmin,vmax,11), orientation='vertical', extend=extend, cax=cax)
    cbar.set_label(label, labelpad=40, rotation=270, size=18, fontweight='normal', ha='center')
    cbar.ax.tick_params(labelsize=18)
    return fig

# --- PI b maps (one figure per run; the CTRL one keeps the original file name) --- #
for key in runs:
    fields=[at_level(b_jas[key]['pi'], l) for l in plot_levels]
    titles=[f'{RUN_LABEL[key]} {l} hPa' for l in plot_levels]
    fig=level_maps(fields, None, titles, bvmin, bvmax, 21, r'JAS $b = h_{10m} - h^{*}$ [kJ kg$^{-1}$]')
    sfx = '' if key=='ctrl' else f'_{key}'
    plt.savefig(FIG_DIR + f'/b_flor_pi_baseline_maps{sfx}.png', transparent=False, bbox_inches='tight')
    plt.savefig(FIG_DIR + f'/b_flor_pi_baseline_maps{sfx}.pdf', transparent=False, bbox_inches='tight')
    plt.close('all')

# --- HIgbl - CTRL in PI (stippling = NOT significant) --- #
if RUN_HITOPO:
    fields=[at_level(topo_diff['pi'], l) for l in plot_levels]
    masks =[at_level(topo_diff_mask['pi'], l) for l in plot_levels]
    titles=[rf'HI$_{{gbl}}-$CTRL {l} hPa' for l in plot_levels]
    fig=level_maps(fields, masks, titles, dvmin, dvmax, 21, r'PI HI$_{gbl}-$CTRL $\Delta b$ [kJ kg$^{-1}$]')
    plt.savefig(FIG_DIR + '/b_flor_pi_hitopo_minus_ctrl_maps.png', transparent=False, bbox_inches='tight')
    plt.savefig(FIG_DIR + '/b_flor_pi_hitopo_minus_ctrl_maps.pdf', transparent=False, bbox_inches='tight')
    plt.close('all')

# ### 2xCO$_2$$-$PI change in JAS buoyancy index $\Delta b$

if RUN_2XCO2:   # needs the 2xCO2 case
    # --- 2xCO2 - PI Delta b maps (CTRL); stippling = NOT significant (p>0.1) --- #
    for key in runs:
        fields=[at_level(future_diff[key], l) for l in plot_levels]
        masks =[at_level(future_diff_mask[key], l) for l in plot_levels]
        titles=[f'{RUN_LABEL[key]} {l} hPa' for l in plot_levels]
        fig=level_maps(fields, masks, titles, dvmin, dvmax, 21, r'2xCO$_2-$PI $\Delta b$ [kJ kg$^{-1}$]')
        sfx = '' if key=='ctrl' else f'_{key}'
        plt.savefig(FIG_DIR + f'/b_flor_2xco2_delta_maps{sfx}.png', transparent=False, bbox_inches='tight')
        plt.savefig(FIG_DIR + f'/b_flor_2xco2_delta_maps{sfx}.pdf', transparent=False, bbox_inches='tight')
        plt.close('all')

    # difference of differences, HIgbl - CTRL warming response (no significance test)
    if RUN_HITOPO:
        fields=[at_level(future_diff['hitopo'] - future_diff['ctrl'], l) for l in plot_levels]
        titles=[rf'HI$_{{gbl}}-$CTRL {l} hPa' for l in plot_levels]
        fig=level_maps(fields, None, titles, dvmin, dvmax, 21, r'$\Delta b$ response: HI$_{gbl}-$CTRL [kJ kg$^{-1}$]')
        plt.savefig(FIG_DIR + '/b_flor_2xco2_delta_hitopo_minus_ctrl_maps.png', transparent=False, bbox_inches='tight')
        plt.savefig(FIG_DIR + '/b_flor_2xco2_delta_hitopo_minus_ctrl_maps.pdf', transparent=False, bbox_inches='tight')
        plt.close('all')

# ### Regional-mean $\Delta b$ profiles (NAM south/core and north domains)

if RUN_2XCO2:   # needs the 2xCO2 case
    # --- area-weighted regional means of Delta b vs pressure --- #
    fig, axs = plt.subplots(nrows=1, ncols=len(domains), figsize=(6*len(domains),7), sharey=True, layout='constrained')
    for ax, dom in zip(np.atleast_1d(axs), domains):
        for key, color in zip(runs, ['k', 'tab:blue']):   # ctrl, hitopo
            prof = regional_weighted_mean(future_diff[key], regions).sel(region=dom) / 1000.
            ax.plot(prof, prof.plev/100., color=color, lw=3, label=key.upper() if key=='ctrl' else r'HI$_{gbl}$')
        ax.axvline(0, color='gray', lw=1)
        ax.set_title(dom, fontsize=20)
        ax.set_xlabel(r'2xCO$_2-$PI $\Delta b$ [kJ kg$^{-1}$]', fontsize=16)
        ax.tick_params(labelsize=14)
    axs[0].set_ylabel('Pressure [hPa]', fontsize=16)
    axs[0].invert_yaxis()
    axs[0].legend(fontsize=14)

    plt.savefig(FIG_DIR + '/b_flor_2xco2_delta_profiles.png', transparent=False, bbox_inches='tight')
    plt.savefig(FIG_DIR + '/b_flor_2xco2_delta_profiles.pdf', transparent=False, bbox_inches='tight')
    plt.close('all')

# ### Longitude$-$pressure transect of $\Delta b$

if RUN_2XCO2:   # needs the 2xCO2 case
    # --- transect: latitude-band mean, longitude vs pressure --- #
    tr_lat = (22, 32)       # degN band averaged
    tr_lon = (240, 265)     # degE extent shown

    def band_mean(da):
        return da.sel(lat=slice(*tr_lat), lon=slice(*tr_lon)).mean('lat')

    tr_diff={}; tr_mask={}
    for key in runs:
        tr_diff[key], tr_mask[key], _ = sigtest2n(band_mean(b_ann[key]['2xco2']), band_mean(b_ann[key]['pi']),
                                                     band_mean(b_jas[key]['2xco2']), band_mean(b_jas[key]['pi']))

    norm=mpl.colors.BoundaryNorm(np.linspace(dvmin, dvmax, 21), cmap.N)
    fig, axs = plt.subplots(nrows=1, ncols=len(runs), figsize=(9*len(runs),7), sharey=True, layout='constrained', squeeze=False)
    for ax, key in zip(axs[0], runs):
        d = tr_diff[key]/1000.
        cf = ax.pcolormesh(d.lon, d.plev/100., d, cmap=cmap, norm=norm)
        lon2d, p2d = np.meshgrid(d.lon, d.plev/100.)
        ns = tr_mask[key].isnull().values
        ax.scatter(lon2d[ns], p2d[ns], s=25, c='k', marker='.', linewidths=0, zorder=100)
        ax.set_title(RUN_LABEL[key], fontsize=20)
        ax.set_xlabel('Longitude [$^\\circ$E]', fontsize=16)
        ax.tick_params(labelsize=14)
    axs[0,0].set_ylabel('Pressure [hPa]', fontsize=16)
    axs[0,0].invert_yaxis()
    cbar=fig.colorbar(cf, ax=axs, extend='both', ticks=np.linspace(dvmin, dvmax, 11), shrink=0.9)
    cbar.set_label(r'2xCO$_2-$PI $\Delta b$ [kJ kg$^{-1}$]', size=16)

    plt.savefig(FIG_DIR + '/b_flor_2xco2_delta_transect.png', transparent=False, bbox_inches='tight')
    plt.savefig(FIG_DIR + '/b_flor_2xco2_delta_transect.pdf', transparent=False, bbox_inches='tight')
    plt.close('all')

# ### Latitude$-$pressure transect of $b$ (108$^\circ$W)
#
# Following Fig. 3 of Pascale et al.: black contours = PI climatological JAS-mean $b$ (dashed = negative; thick zero contour = level of free convection, LFC), black line = LCL (levels below it are masked), black fill = model terrain (JAS-mean surface pressure; the Sierra Madre Occidental near 28$^\circ$N). With only the PI case enabled the shading is the PI-mean $b$ itself; when `2xco2` is enabled the shading becomes the 2xCO$_2-$PI $\Delta b$ (stippling = significant, blue outline = land points with a significant negative $\Delta b$). Edit `tr_lon0`, `tr_lat_rng` and the levels at the top of the cell.

## +++ LATITUDE-PRESSURE TRANSECT +++ ##

RUN_2XCO2 = '2xco2' in cases
tr_lon0 = 252.                  # degE (108 W): longitude of the transect (nearest model grid column is used)
tr_lat_rng = (5, 40)            # degN range shown
tr_pmin, tr_pmax = 100, 1000    # hPa axis limits
stipple_significant = True      # True: stipple where the change IS significant (as in the paper); False: where it is NOT (as in the maps above)
b_clim_levels = [-40, -30, -20, -15, -10, -5, 5, 10, 15, 20, 30, 40]   # kJ/kg, black contours of the PI climatological b (0 = LFC is drawn separately)
dtr_levels = [-5, -4, -3, -2, -1, -0.5, -0.1, 0.1, 0.5, 1, 2, 3, 4, 5]  # kJ/kg, shading of Delta b (2xCO2 - PI)
dtr_cmap = plt.get_cmap('RdYlBu', len(dtr_levels) - 1)                  # red = more negative buoyancy (stronger inhibition)

ps_jas = {key: jas_seasonal_mean(dat[key]['pi']['ps']).load() for key in runs}   # JAS-mean surface pressure [Pa]: the terrain (SMO)

def col(da):
    """Column nearest `tr_lon0`, restricted to `tr_lat_rng`."""
    return da.sel(lon=tr_lon0, method='nearest').sel(lat=slice(*tr_lat_rng))

def lat_pres_transect(key):
    """Latitude-pressure transect (as Fig. 3a-c of Pascale et al.): shading = 2xCO2-PI change in b (PI b if 2xCO2 is
    not enabled), black contours = PI climatological b (dashed = negative, thick = zero = LFC), black line = LCL,
    black fill = terrain, blue outline = land points with a significant negative change in b."""
    b_pi = col(b_jas[key]['pi']) / 1000.                      # (plev, lat), kJ/kg; NaN below the LCL / ground
    lat, p_hpa = b_pi.lat.values, b_pi.plev.values / 100.
    lonw = 360. - float(b_pi.lon)
    plcl = col(plcl_jas[key]['pi']).values / 100.             # hPa
    ps_c = col(ps_jas[key]).values / 100.                     # hPa
    land = col(topo[key]) > 0                                 # land proxy: model surface height > 0 m

    fig, ax = plt.subplots(figsize=(7, 6), layout='constrained')
    if RUN_2XCO2:
        d = col(future_diff[key]) / 1000.
        sig = col(future_diff_mask[key]).notnull()            # True where significant (p <= 0.1)
        cf = ax.contourf(lat, p_hpa, d.values, levels=dtr_levels, cmap=dtr_cmap, extend='both')
        cbar_label = r'2xCO$_2-$PI $\Delta b$ [kJ kg$^{-1}$]'
        st = (sig.values if stipple_significant else ~sig.values) & d.notnull().values
        pp, ll = np.meshgrid(p_hpa, lat, indexing='ij')
        ax.scatter(ll[st], pp[st], s=8, c='k', marker='.', linewidths=0, zorder=5)
        neg_sig = ((d < 0) & sig & land).astype(float)        # significant negative change over land
        if float(neg_sig.max()) > 0:
            ax.contour(lat, p_hpa, neg_sig.values, levels=[0.5], colors='tab:blue', linewidths=2.5, zorder=6)
    else:
        cf = ax.contourf(lat, p_hpa, b_pi.values, levels=np.linspace(bvmin, bvmax, 21), cmap=cmap, extend='both')
        cbar_label = r'PI $b$ [kJ kg$^{-1}$]'

    cs = ax.contour(lat, p_hpa, b_pi.values, levels=b_clim_levels, colors='k', linewidths=1.2, zorder=4)   # negative: dashed
    ax.clabel(cs, fmt='%g', fontsize=9)
    lfc = ax.contour(lat, p_hpa, b_pi.values, levels=[0], colors='k', linewidths=3, zorder=4)
    ax.clabel(lfc, fmt={0: 'LFC'}, fontsize=11)
    ax.plot(lat, plcl, color='k', lw=1, zorder=4, label='LCL')
    ax.fill_between(lat, ps_c, tr_pmax, where=ps_c < tr_pmax, color='k', zorder=7)
    if np.nanmin(ps_c) < tr_pmax:   # label the highest terrain along the transect
        ax.text(lat[np.nanargmin(ps_c)], (np.nanmin(ps_c) + tr_pmax) / 2, 'SMO', color='w', weight='bold', ha='center', va='center', zorder=8)
    ax.set_ylim(tr_pmax, tr_pmin); ax.set_xlim(*tr_lat_rng)
    ax.set_xlabel(r'Latitude [$^\circ$N]', fontsize=14); ax.set_ylabel('Pressure [hPa]', fontsize=14)
    ax.set_title(f'{RUN_LABEL[key]} JAS, {lonw:.0f}$^\\circ$W', fontsize=16)
    ax.tick_params(labelsize=12); ax.legend(loc='upper left', fontsize=10)
    cbar = fig.colorbar(cf, ax=ax, orientation='horizontal', ticks=cf.levels, shrink=0.9)
    cbar.set_label(cbar_label, size=13); cbar.ax.tick_params(labelsize=10)
    return fig, lonw

for key in runs:
    fig, lonw = lat_pres_transect(key)
    stem = f'{FIG_DIR}/b_flor_' + ('2xco2_delta' if RUN_2XCO2 else 'pi') + f'_lat_pres_transect_{lonw:.0f}W'
    plt.savefig(stem + '.png', transparent=False, bbox_inches='tight')
    plt.savefig(stem + '.pdf', transparent=False, bbox_inches='tight')
plt.close('all')

# ## Convective mass flux (optional)
#
# Change in the diagnosed RAS cumulus mass flux. Disabled by default: the file and variable name are not recorded in the repo (confirm on Discover; GFDL `mc_full` is a common name for the net convective mass flux).

RUN_MASSFLUX = False
MF_NAME = 'mc_full'       # CONFIRM: variable (and file tag) of the diagnosed convective mass flux, on pressure levels

if RUN_MASSFLUX and RUN_2XCO2:
    mf_jas, mf_ann = {}, {}
    for key in runs:
        mf_jas[key], mf_ann[key] = {}, {}
        for case in cases:
            f = f'{dpath0}/FLOR/{key}/{case}/flor.{key}.{"" if case=="pi" else case+"."}{MF_NAME}.monthly.nc'
            ds = load_var(f, MF_NAME, case, is_3d=True)
            mf_jas[key][case] = jas_seasonal_mean(ds).load()
            mf_ann[key][case] = jas_yearly_mean(ds).load().transpose('year','plev','lat','lon')
    mf_diff, mf_mask = {}, {}
    for key in runs:
        mf_diff[key], mf_mask[key], _ = sigtest2n(mf_ann[key]['2xco2'], mf_ann[key]['pi'], mf_jas[key]['2xco2'], mf_jas[key]['pi'])
    # profiles over the NAM domains
    fig, axs = plt.subplots(nrows=1, ncols=len(domains), figsize=(6*len(domains),7), sharey=True, layout='constrained')
    for ax, dom in zip(np.atleast_1d(axs), domains):
        for key in runs:
            prof = regional_weighted_mean(mf_diff[key], regions).sel(region=dom)
            ax.plot(prof, prof.plev/100., lw=3, label=key)
        ax.axvline(0, color='gray', lw=1); ax.set_title(dom, fontsize=20)
        ax.set_xlabel(f'2xCO$_2-$PI $\\Delta$ {MF_NAME} [{mf_jas[runs[0]]["pi"].attrs.get("units","")}]', fontsize=16)
    axs[0].set_ylabel('Pressure [hPa]', fontsize=16); axs[0].invert_yaxis(); axs[0].legend(fontsize=14)
    plt.savefig(FIG_DIR + '/massflux_flor_2xco2_delta_profiles.png', transparent=False, bbox_inches='tight')
    plt.savefig(FIG_DIR + '/massflux_flor_2xco2_delta_profiles.pdf', transparent=False, bbox_inches='tight')
    plt.close('all')
