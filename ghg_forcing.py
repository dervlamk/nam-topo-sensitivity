"""
Radiative forcing from mean greenhouse gas concentrations over user-defined time intervals.

Reads the Meinshausen et al. (2017, 2020) CMIP6 concentration tables (whitespace-delimited .txt) and
converts the interval-MEAN concentrations to forcing (W m-2) relative to 1750 using the
expressions in Table 1 of NOAA's AGGI page (https://gml.noaa.gov/aggi/aggi.html), which follow
Etminan et al. (2016), Meinshausen et al. (2020), Forster et al. (2021) and Smith et al. (2021).
The forcing is computed from the mean concentrations, not as the mean of annual forcings.

    CO2:  dF = (a1' + c1*sqrt(N)) * ln(C/C0) * e1,   a1' = d1 + a1*(C-C0)^2 + b1*(C-C0)
    CH4:  dF = (a2*sqrt(M) + b2*sqrt(N) + d2) * (sqrt(M) - sqrt(M0)) * e2
    N2O:  dF = (a3*sqrt(C) + b3*sqrt(N) + c3*sqrt(M) + d3) * (sqrt(N) - sqrt(N0)) * e3

with C = CO2 [ppm], M = CH4 [ppb], N = N2O [ppb]. The halogenated gases are not included.

Usage
-----
    python ghg_forcing.py 1850:1879 1995:2014
    python ghg_forcing.py 1995:2014 --csv forcing.csv
    python ghg_forcing.py 1850:1869 --experiment 1pctco2
    python ghg_forcing.py 1995:2014 2030:2050 --experiment ssp585

Each interval is START:END in years, inclusive. As a module:

    from ghg_forcing import load_ghg, mean_concentrations, forcing_from_concentrations
    ghg = load_ghg()
    forcing_from_concentrations(mean_concentrations(ghg, 1995, 2014))
"""

import argparse
import os

import numpy as np
import pandas as pd

DEFAULT_DATA_DIR = '/Users/dervlamk/data/cmip6'
GASES = ['CO2', 'CH4', 'N2O']  # columns used; halogens are ignored
HISTORICAL_FILE = 'Meinshausen2017_GHG-concentrations_1850-2014_CMIP6_historical.txt'
ONEPCT_FILE = 'Meinshausen2017_GHG-concentrations_CMIP6_1pctCO2.txt'
SCENARIO_FILE = 'Meinshausen2020_CMIP6_{scenario}_GHG.txt'  # 2015-2500 global means
SCENARIOS = ['ssp245', 'ssp585']
EXPERIMENTS = ['historical', '1pctco2'] + SCENARIOS

# --- AGGI Table 1 constants (NOAA AGGI, Table 1) --- #
# 1750 (unperturbed) abundances. These are the AGGI baseline values, which the Table 1 coefficients
# were fitted with, and they are used on purpose. They differ from the 1750 values in
# Meinshausen2017_GHG-concentrations_1750_reference.txt (CO2 277.15 ppm, CH4 731.41 ppb,
# N2O 273.87 ppb), which this script does NOT read. For 1995-2014, switching to the Meinshausen
# 1750 values changes the CO2 forcing by about +0.02 W m-2 and the N2O forcing by about
# -0.014 W m-2 (total about +0.01 W m-2).
C0, M0, N0 = 278.3, 729.2, 270.1   # CO2 [ppm], CH4 [ppb], N2O [ppb]

CO2_COEF = dict(a=-2.48e-07, b=7.59e-04, c=-2.15e-03, d=5.2488, e=1.05)
CH4_COEF = dict(a=-8.96e-05, b=-1.25e-04, d=0.045194, e=0.86)
N2O_COEF = dict(a=-3.42e-04, b=2.546e-04, c=-2.44e-04, d=0.12173, e=1.07)


# ------------------------------------------------------------------ #
#                              data I/O                              #
# ------------------------------------------------------------------ #
def _read_table(path):
    df = pd.read_csv(path, comment='#', sep=r'\s+')
    df['year'] = df['year'].astype(int)
    return df.set_index('year')


def load_ghg(data_dir=DEFAULT_DATA_DIR, experiment='historical'):
    """Annual global-mean concentrations, indexed by year.

    experiment='historical': Table 6 (1850-2014).
    experiment='1pctco2': Table 9 CO2 (1850-2100); the table only provides CO2 and states that
        all other gases follow piControl, so CH4 and N2O are held at the 1850 values of
        Table 6 (the piControl/1850 concentrations).
    experiment='ssp245' / 'ssp585': the historical record (Table 6, 1850-2014) followed by the
        scenario concentrations (2015-2500), so intervals may span 1850-2500.
    """
    hist = _read_table(os.path.join(data_dir, HISTORICAL_FILE))[GASES]
    if experiment == 'historical':
        return hist
    if experiment == '1pctco2':
        co2 = _read_table(os.path.join(data_dir, ONEPCT_FILE))
        out = pd.DataFrame(index=co2.index)
        for col in hist.columns:
            out[col] = hist.loc[1850, col]
        out['CO2'] = co2['CO2']
        return out
    if experiment in SCENARIOS:
        scen = _read_table(os.path.join(data_dir, SCENARIO_FILE.format(scenario=experiment)))[GASES]
        if scen.index.min() != hist.index.max() + 1:
            raise ValueError(f'scenario starts in {scen.index.min()} but the historical record '
                             f'ends in {hist.index.max()}; they are not contiguous')
        return pd.concat([hist, scen])
    raise ValueError(f'experiment must be one of {EXPERIMENTS}, got {experiment!r}')


def mean_concentrations(ghg, start_year, end_year):
    """Mean concentration of every column over start_year..end_year (inclusive)."""
    if start_year > end_year:
        raise ValueError(f'interval start {start_year} is after end {end_year}')
    if start_year < ghg.index.min() or end_year > ghg.index.max():
        raise ValueError(f'interval {start_year}-{end_year} is outside the data '
                         f'range {ghg.index.min()}-{ghg.index.max()}')
    return ghg.loc[start_year:end_year].mean()


# ------------------------------------------------------------------ #
#                       forcing (AGGI Table 1)                       #
# ------------------------------------------------------------------ #
def forcing_co2(C, N):
    """CO2 forcing [W m-2]; C in ppm, N = N2O in ppb."""
    k = CO2_COEF
    alpha = k['d'] + k['a'] * (C - C0) ** 2 + k['b'] * (C - C0)
    return (alpha + k['c'] * np.sqrt(N)) * np.log(C / C0) * k['e']


def forcing_ch4(M, N):
    """CH4 forcing [W m-2]; M in ppb, N = N2O in ppb."""
    k = CH4_COEF
    return (k['a'] * np.sqrt(M) + k['b'] * np.sqrt(N) + k['d']) * (np.sqrt(M) - np.sqrt(M0)) * k['e']


def forcing_n2o(C, M, N):
    """N2O forcing [W m-2]; C = CO2 in ppm, M = CH4 in ppb, N in ppb."""
    k = N2O_COEF
    return (k['a'] * np.sqrt(C) + k['b'] * np.sqrt(N) + k['c'] * np.sqrt(M) + k['d']) \
        * (np.sqrt(N) - np.sqrt(N0)) * k['e']


def forcing_from_concentrations(conc):
    """Per-gas and total forcing [W m-2] from a Series (or dict) of (mean) concentrations
    with CO2 [ppm], CH4 [ppb] and N2O [ppb]."""
    C, M, N = conc['CO2'], conc['CH4'], conc['N2O']
    out = pd.Series({
        'CO2': forcing_co2(C, N),
        'CH4': forcing_ch4(M, N),
        'N2O': forcing_n2o(C, M, N),
    })
    out['Total'] = out.sum()
    return out


# ------------------------------------------------------------------ #
#                                 CLI                                #
# ------------------------------------------------------------------ #
def _parse_interval(text):
    try:
        a, b = text.split(':')
        return int(a), int(b)
    except ValueError:
        raise argparse.ArgumentTypeError(f"interval must look like START:END (e.g. 1995:2014), got {text!r}")


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('intervals', nargs='+', type=_parse_interval, metavar='START:END',
                   help='inclusive year intervals, e.g. 1850:1879 1995:2014')
    p.add_argument('--data-dir', default=DEFAULT_DATA_DIR, help='directory with the Meinshausen .txt files')
    p.add_argument('--experiment', choices=EXPERIMENTS, default='historical',
                   help='ssp245/ssp585 = historical (to 2014) followed by the scenario (2015-2500)')
    p.add_argument('--csv', help='optional path to write the results table')
    args = p.parse_args(argv)

    ghg = load_ghg(args.data_dir, args.experiment)
    conc_rows, force_rows = {}, {}
    for y0, y1 in args.intervals:
        label = f'{y0}-{y1}'
        conc = mean_concentrations(ghg, y0, y1)
        conc_rows[label] = conc
        force_rows[label] = forcing_from_concentrations(conc)

    conc_df = pd.DataFrame(conc_rows).T
    force_df = pd.DataFrame(force_rows).T
    pd.set_option('display.width', 200)
    print(f'\nExperiment: {args.experiment}')
    print('\nMean concentrations (CO2 ppm, CH4/N2O ppb):')
    print(conc_df[GASES].round(2).to_string())
    print('\nRadiative forcing relative to 1750 [W m-2]:')
    print(force_df.round(3).to_string())

    if args.csv:
        out = pd.concat({'conc': conc_df, 'forcing_Wm2': force_df}, axis=1)
        out.to_csv(args.csv, index_label='interval')
        print(f'\nWrote {args.csv}')


if __name__ == '__main__':
    main()
