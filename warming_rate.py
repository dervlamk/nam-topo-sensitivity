#!/usr/bin/env python
"""Average warming rate (degrees/year) of t_ref over given year intervals,
plus a line plot of the annual global-mean series with fitted trends.

Usage: python warming_rate.py file1.nc file2.nc
"""
import sys

import matplotlib.pyplot as plt
import numpy as np
import xarray as xr

# ---- user-defined settings -------------------------------------------------
LABELS = ["CTRL", "HIgbl"]  # legend label for each input file
COLORS = ["tab:peru", "tab:red"]         # line color for each input file
TITLE = "Response of Annual Global Mean Surface Temperature to CO2 Forcing" # plot title
YLABEL = "Global-mean T-2m (deg C)"         # y-axis label
OUTFILE = "figs/t_ref_warming_rate"      # saved as .png and .pdf
INTERVALS = [(100, 170), (200, 299)]     # model years, inclusive
# -----------------------------------------------------------------------------


def annual_global_mean(path):
    """Area-weighted (cos lat) global mean of t_ref, averaged to annual."""
    da = xr.open_dataset(path, use_cftime=True)["t_ref"]
    w = np.cos(np.deg2rad(da["grid_yt"]))
    gm = da.weighted(w).mean(("grid_yt", "grid_xt"))
    return gm.groupby("time.year").mean("time")


def fit_trend(annual, y0, y1):
    """Least-squares (slope, intercept) over years y0..y1 inclusive."""
    seg = annual.sel(year=slice(y0, y1))
    if seg.sizes["year"] < 2:
        raise ValueError(f"no data in {y0}-{y1} (years available: "
                         f"{int(annual.year[0])}-{int(annual.year[-1])})")
    return np.polyfit(seg["year"].values, seg.values, 1)


if __name__ == "__main__":
    paths = sys.argv[1:3]
    fig, ax = plt.subplots(figsize=(8, 4.5))

    for path, label, color in zip(paths, LABELS, COLORS):
        annual = annual_global_mean(path)
        ax.plot(annual.year, annual, color=color, lw=1, label=label)
        print(path)
        for y0, y1 in INTERVALS:
            slope, icpt = fit_trend(annual, y0, y1)
            yrs = np.array([y0, y1])
            ax.plot(yrs, slope * yrs + icpt, color=color, lw=2.5, ls="--")
            ax.annotate(f"{slope:+.4f} °C/yr", xy=(yrs.mean(), slope * yrs.mean() + icpt),
                        xytext=(0, 8), textcoords="offset points",
                        ha="center", color=color, fontsize=9, fontweight="bold")
            print(f"  years {y0}-{y1}: {slope:+.5f} deg/yr")

    ax.set_xlabel("Model year")
    ax.set_ylabel(YLABEL)
    ax.set_title(TITLE)
    ax.legend(frameon=False)
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(f"{OUTFILE}.{ext}", dpi=200, bbox_inches="tight")
