#!/usr/bin/env python
"""Compute FLOR reference-height specific humidity (q_ref) from rh_ref, t_ref and ps.

Reproduces the method of the old qs_flor.ipynb (repo history, commit b94369f):
    e  = RH * es(T)                      es: saturation vapor pressure (Bolton 1980, same as metpy)
    w  = Rd * e / (Rv * (p - e))         mixing ratio
    q  = w / (1 + w)                     specific humidity
with RH = rh_ref/100, T = t_ref [K], p = ps [Pa].

Input files (in dpath0/FLOR/{run}/{case}/):
    pi    : flor.{run}.{var}.monthly.nc
    2xco2 : flor.{run}.2xco2.{var}.monthly.nc          var in rh_ref, t_ref, ps
Output (same directory): flor.{run}[.{case}].q_ref.{START}-{END}.monthly.nc, variable q_ref [kg/kg].
The years are part of the file name so that runs for different periods do not collide. The older
untagged file flor.ctrl.q_ref.monthly.nc covers PI years 0201-0300 only (it was the last 100 years of the
run), so it is never read or overwritten by default.

Default years: pi 0101-0200, 2xco2 0200-0299 (the periods used in the analysis notebooks).

Usage:
    python q_ref_from_rh.py RUN CASE [--start YYYY --end YYYY | --all-years] [--out FILE] [--overwrite]
Examples:
    python q_ref_from_rh.py ctrl pi                        # PI years 0101-0200 -> flor.ctrl.q_ref.0101-0200.monthly.nc
    python q_ref_from_rh.py ctrl pi --start 0201 --end 0300   # the period of the existing untagged file
    python q_ref_from_rh.py ctrl 2xco2                     # 2xCO2 years 0200-0299
    python q_ref_from_rh.py ctrl pi --out /path/to/scratch/q_ref_test.nc   # write elsewhere (e.g. to test)
The script refuses to overwrite an existing file unless --overwrite is given.
"""
import argparse
import os
import sys

import numpy as np
import xarray as xr

# ---- user-defined settings ----
DPATH0 = '/discover/nobackup/projects/giss/baldwin_nip/dmkumar'  # top level data directory
VARNS = ['rh_ref', 't_ref', 'ps']
DEFAULT_YEARS = {'pi': ('0101', '0200'), '2xco2': ('0200', '0299')}  # same periods as the analysis notebooks
RD = 287.05   # J kg-1 K-1, gas constant of dry air (as in qs_flor.ipynb)
RV = 461.50   # J kg-1 K-1, gas constant of water vapor (as in qs_flor.ipynb)


def flor_file(run, case, varn, years=None):
    """Path of a FLOR monthly file, following the repo's naming convention.
    `years` = (start, end) inserts the period into the name (used for the q_ref output)."""
    tag = f'flor.{run}' if case == 'pi' else f'flor.{run}.{case}'
    yr = f'.{years[0]}-{years[1]}' if years else ''
    return f'{DPATH0}/FLOR/{run}/{case}/{tag}.{varn}{yr}.monthly.nc'


def saturation_vapor_pressure(T):
    """Saturation vapor pressure [Pa] over liquid water, Bolton (1980); matches metpy's default. T in K."""
    return 611.2 * np.exp(17.67 * (T - 273.15) / (T - 29.65))


def specific_humidity_from_rh(RH, T, p):
    """Specific humidity [kg/kg] from relative humidity (fraction), temperature [K] and pressure [Pa]."""
    e = RH * saturation_vapor_pressure(T)
    w = (e * RD) / (RV * (p - e))
    return w / (w + 1.0)


def load(run, case, varn, start, end):
    """Open one variable, optionally restricted to years start..end (e.g. '0101', '0200')."""
    path = flor_file(run, case, varn)
    if not os.path.exists(path):
        sys.exit(f'Missing input file: {path}')
    da = xr.open_dataset(path, use_cftime=True, chunks={'time': 120})[varn]
    if start is not None or end is not None:
        avail = f'{da.time.values[0]} to {da.time.values[-1]}'
        da = da.sel(time=slice(start, end))
        if da.sizes['time'] == 0:
            sys.exit(f'{path}: no time steps in {start}-{end}; the file covers {avail}.')
    return da


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('run', help="FLOR run, e.g. ctrl or hitopo")
    ap.add_argument('case', choices=['pi', '2xco2'])
    ap.add_argument('--start', default=None, help="first year, e.g. 0101 (default: pi 0101, 2xco2 0200)")
    ap.add_argument('--end', default=None, help="last year, e.g. 0200 (default: pi 0200, 2xco2 0299)")
    ap.add_argument('--all-years', action='store_true', help="use every time step in the files (output name has no years)")
    ap.add_argument('--out', default=None, help="output file (default: q_ref file next to the inputs, with the years in the name)")
    ap.add_argument('--overwrite', action='store_true', help="allow overwriting an existing output file")
    args = ap.parse_args()

    if args.all_years:
        start = end = None
    else:
        start = args.start or DEFAULT_YEARS[args.case][0]
        end = args.end or DEFAULT_YEARS[args.case][1]
    years = (start, end) if start is not None else None

    out = args.out or flor_file(args.run, args.case, 'q_ref', years=years)
    if os.path.exists(out) and not args.overwrite:
        sys.exit(f'{out} already exists; use --out for a different file or --overwrite to replace it.')

    dat = {v: load(args.run, args.case, v, start, end) for v in VARNS}

    # all three inputs must be on identical time/space axes
    try:
        dat['rh_ref'], dat['t_ref'], dat['ps'] = xr.align(dat['rh_ref'], dat['t_ref'], dat['ps'], join='exact')
    except ValueError as err:
        sys.exit(f'rh_ref, t_ref and ps do not share the same coordinates (check years/grids): {err}')
    if dat['t_ref'].sizes['time'] == 0:
        sys.exit('No time steps selected; check --start/--end against the file time axes.')

    # unit sanity checks on a small sample (rh in %, t in K, ps in Pa)
    rh_max = float(dat['rh_ref'].isel(time=0).max())
    t_mean = float(dat['t_ref'].isel(time=0).mean())
    ps_mean = float(dat['ps'].isel(time=0).mean())
    if rh_max <= 1.5:
        sys.exit(f'rh_ref max is {rh_max:.2f}: looks like a fraction, expected percent; refusing to divide by 100.')
    if not 200 < t_mean < 330:
        sys.exit(f't_ref mean is {t_mean:.1f}: expected Kelvin.')
    if not 5e4 < ps_mean < 1.1e5:
        sys.exit(f'ps mean is {ps_mean:.1f}: expected Pa.')

    q = specific_humidity_from_rh(dat['rh_ref'] / 100.0, dat['t_ref'], dat['ps'])
    q.name = 'q_ref'
    q.attrs = {'long_name': 'surface (2m) specific humidity', 'units': 'kg/kg',
               'method': 'calculated from relative humidity (rh_ref, t_ref, ps) with q_ref_from_rh.py'}

    print(f'{args.run} {args.case}: {q.sizes["time"]} time steps, '
          f'{q.time.values[0]} to {q.time.values[-1]} -> {out}')
    q.to_netcdf(out, mode='w')

    # report range of the result
    qmin, qmax = float(q.min()), float(q.max())
    print(f'q_ref range: {qmin:.5f} to {qmax:.5f} kg/kg')


if __name__ == '__main__':
    main()
