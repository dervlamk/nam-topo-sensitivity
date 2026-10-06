import os
import sys
import xarray as xr
import cf_xarray
import netCDF4 as nc
import numpy as np
import metpy.calc as mp
from datetime import datetime


#############################

def get_xy_coords(var):
    """
    Get lon and lat arrays without knowing coordinate names
    """
    if isinstance(var, xr.DataArray):
        x,y=var.metpy.coordinates('x','y')
        return(x,y)
    if isinstance(var, xr.Dataset):
        print('This is a dataset. Please use an xarray DataArray')

def get_season(season='ann'):
    """
    Index months to average over to derive an annual or seasonal mean
        - can only be applied to monthly climatologies
    """
    if season in ['ANNUAL', 'ANN', 'ann']:
        mons = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]
    if season in ['DJF', 'djf']:
        mons = [0, 1, 11]
    if season in ['JFM', 'jfm']:
        mons = [0, 1, 2]
    if season in ['MAM', 'mam']:
        mons = [2, 3, 4]
    if season in ['JJA', 'jja']:
        mons = [5, 6, 7]
    if season in ['JJAS', 'jjas']:
        mons = [5, 6, 7, 8]
    if season in ['JAS', 'jas']:
        mons = [6, 7, 8]
    if season in ['SON', 'son']:
        mons = [8, 9, 10]
    if season==None:
        pass
    return mons

def lonFlip(var):
    """
    Convert longitude values from the -180:180 to 0:360 convention or vice versa.

    ** Works for both global data due to auto-detection of longitude convention **
    More efficient than rolling: only relabels coordinates + sorts.

    Parameters
    ----------
    var : xr.DataArray or xr.Dataset
    """

    #=== Get var info
    try:
        lon_name = var.cf.axes["X"][0]
    except KeyError:
        # fallback: find coordinate with 'lon' in its name
        lon_name = [c for c in var.coords if 'lon' in c.lower()][0]
    # extract lon array
    lon=var[lon_name]

    #=== Detect current longitude convention and wrap values
    if lon.min() < 0:
        # -180:180 -> 0:360
        new_lon = lon % 360
        target_range = "0:360"
    else:
        # 0:360 -> -180:180
        new_lon = ((lon + 180) % 360) - 180
        target_range = "-180:180"

    #=== Assign and sort
    var = var.assign_coords({lon_name: new_lon}).sortby(lon_name)

    #=== Add history
    timestamp = datetime.now().strftime("%B %d, %Y, %r")
    hist_message = f"wrapped longitudes to {target_range} on {timestamp}"
    if isinstance(var, xr.DataArray):
        var.attrs["history"] = hist_message
    else: # Dataset
        var.attrs["history"] = var.attrs.get("history","") + "\n" + hist_message

    return var

def longitude_flip(var):
    """
    Convert longitude values from the -180:180 to 0:360 convention or vice versa.
        
    ** Only works for global data. Do not apply to data with a clipped longitude range **
    ** Is computationally expensive because it forces a particular layout of the data in memory **
    ** Probably not optimal to use for most scenarios **
        
    Parameters
    ----------
    var : Data Array
    """    
    # get var info
    x,_=get_xy_coords(var) # extract original longitude values
    lon_name=x.name        # store name of longitude coordinate
    nx=len(x)              # longitude resolution
    
    # determine longitude format and create an array of new lons in opposite convention
    if min(x)<0: 
        # if there are negative values, data is -180:180 and need to switch to 0:360
        new_lons=np.linspace((min(x)+180), (max(x)+180), nx)
    elif max(x)>180:
        # if the max value is >180, data is in 0:360 format and need to switch to -180:180
        new_lons=np.linspace((min(x)-180), (max(x)-180), nx)
        
    # shift the data by 180° of longitude
    nshift=nx//2
    var=var.roll({lon_name: nshift}, roll_coords=False)
            
    # update longitude coord with new values
    var=var.assign_coords({lon_name: new_lons})
    
    # add attributes documenting change
    timestamp=datetime.now().strftime("%B %d, %Y, %r")
    var.attrs['history']=f'flipped longitudes {timestamp}'
    var.attrs['original_lons']=x.values
    
    return(var)

def regrid_like(ref, var):
    """
    Regrid data to match a reference 
    ** Only works for global data **
    
    Parameters
    ----------
    ref : reference array
    var : variable array to regrid
    """ 
    x_ref,y_ref=get_xy_coords(ref) # lat lon coords from reference variable
    x_var,_=get_xy_coords(var) # lat lon coords from variable to be regridded
    
    # if longitudes are referenced differently, flip lons of variable
    if np.sign(min(x_ref)) != np.sign(min(x_var)):
        var=longitude_flip(var)
        x_var,y_var=get_xy_coords(var)
    
    # rename coordinates to match reference
    var=var.rename({x_var.name:x_ref.name, y_var.name:y_ref.name})
    
    # interpolate var data
    var_regridded=var.interp_like(ref, method='linear')
    return(var_regridded)

def _cell_edges(centers):
    """Cell edges for a 1D rectilinear coordinate, placed at midpoints between centers
    (outer edges extrapolated by half a grid step). Returns (lower, upper) arrays."""
    c = np.asarray(centers, dtype=float)
    mid = 0.5 * (c[1:] + c[:-1])
    edges = np.concatenate([[c[0] - (mid[0] - c[0])], mid, [c[-1] + (c[-1] - mid[-1])]])
    return np.minimum(edges[:-1], edges[1:]), np.maximum(edges[:-1], edges[1:])


def conservative_regrid(da, lat, lon, min_coverage=0.5, lat_name='lat', lon_name='lon'):
    """
    First-order conservative (area-weighted) regridding of a regular lat/lon field
    onto another regular lat/lon grid. Each target cell is the spherical-area-weighted
    mean of the source cells overlapping it, so unlike bilinear interpolation it
    represents the target cell's area mean when the source is finer than the target.

    Source cells that are NaN (e.g. ocean in land-only products) are excluded and the
    weights re-normalized over valid source area. Target cells for which less than
    `min_coverage` of the cell area is covered by valid source data are set to NaN.
    Both grids must be 1D, regular, and use the same longitude convention
    (e.g. 0:360); any extra dimensions (time, year, ...) are preserved.

    Parameters
    ----------
    da : xr.DataArray with dims including lat_name and lon_name
    lat, lon : 1D target lat/lon values (array-like or DataArray)
    min_coverage : minimum fraction (0-1) of a target cell covered by valid source data
    """
    lat_out = np.asarray(lat, dtype=float)
    lon_out = np.asarray(lon, dtype=float)

    def overlap(src, dst, spherical):
        s_lo, s_hi = _cell_edges(src)
        d_lo, d_hi = _cell_edges(dst)
        if spherical:  # latitude bands: area proportional to difference in sin(lat)
            f = lambda x: np.sin(np.deg2rad(np.clip(x, -90, 90)))
        else:
            f = lambda x: x
        lo = np.maximum(d_lo[:, None], s_lo[None, :])
        hi = np.minimum(d_hi[:, None], s_hi[None, :])
        return np.clip(f(hi) - f(lo), 0, None)  # (n_dst, n_src)

    w_lat = xr.DataArray(overlap(da[lat_name].values, lat_out, True), dims=['lat_out', 'lat_in'])
    w_lon = xr.DataArray(overlap(da[lon_name].values, lon_out, False), dims=['lon_out', 'lon_in'])
    # total area of each target cell, to measure how much of it the source covers
    t_lat = xr.DataArray(overlap(lat_out, lat_out, True).diagonal(), dims=['lat_out'])
    t_lon = xr.DataArray(overlap(lon_out, lon_out, False).diagonal(), dims=['lon_out'])

    src = da.rename({lat_name: 'lat_in', lon_name: 'lon_in'}).drop_vars(['lat_in', 'lon_in'], errors='ignore')
    valid = src.notnull()

    def wsum(x):  # separable weights: (lat_out x lat_in) @ field @ (lon_in x lon_out)
        return xr.apply_ufunc(
            lambda f, wy, wx: wy @ f @ wx.T, x, w_lat, w_lon,
            input_core_dims=[['lat_in', 'lon_in'], ['lat_out', 'lat_in'], ['lon_out', 'lon_in']],
            output_core_dims=[['lat_out', 'lon_out']],
            dask='parallelized', output_dtypes=[float],
            dask_gufunc_kwargs={'output_sizes': {'lat_out': lat_out.size, 'lon_out': lon_out.size}},
        )

    num = wsum(src.fillna(0))
    den = wsum(valid.astype(float))
    out = (num / den).where(den >= min_coverage * (t_lat * t_lon))
    out = out.rename({'lat_out': lat_name, 'lon_out': lon_name})
    out = out.assign_coords({lat_name: lat_out, lon_name: lon_out})
    out = out.transpose(*da.dims)
    out.attrs = da.attrs
    return out


def lat_weighted_mean(var, lat_name='lat', lon_name='lon'):
    lats = var[lat_name]
    # determine weight based on latitude value
    weights = np.cos(np.deg2rad(lats))
    weights.name = 'weights'
    # calculate area-weighted values
    weighted_var = var.weighted(weights)
    # calculate global mean of weighted data
    weighted_mean = weighted_var.mean(dim=[lat_name,lon_name], keep_attrs=True)
    return(weighted_mean)

def latitude_weighted_mean(var):
    """
    Calculate the mean of geospatial data taking into account unequal grid cell area
    """
    # get x and y coordinate data
    lons,lats = get_xy_coords(var)
    # determine weight based on latitude value
    weights = np.cos(np.deg2rad(lats))
    weights.name = 'weights'
    # calculate area-weighted values
    weighted_var = var.weighted(weights)
    # calculate global mean of weighted data
    weighted_mean = weighted_var.mean(dim=[lats.name,lons.name], keep_attrs=True)
    return(weighted_mean)

def season_mean(ds, calendar='standard'):
	seas_mean = ds.groupby('time.season').mean(dim='time') #sum(dim='time')
	return seas_mean

def annual_season_mean(ds): #, calendar='standard'):
	ds_seasonal = {}
	ann_seasonal_mean = {}
	for season in ['DJF','MAM','JJA','SON']:
		# extract data for season
		ds_seasonal[season] = ds.where(ds['time.season'] == season)
		# I'm not sure what this step is doing
		if season == 'DJF':
			ds_seasonal[season] = ds_seasonal[season].shift(time=1)
		# get timeseries of seasonal mean climatologies
		ann_seasonal_mean[season] = ds_seasonal[season].groupby('time.year').mean(dim='time') #sum(dim='time')
		# but, cut first year of DJF timeseries as there is no Dec data from year -1
		if season == 'DJF':
			year_min = ann_seasonal_mean[season].year.min()+1
			year_max = ann_seasonal_mean[season].year.max()
			ann_seasonal_mean[season] = ann_seasonal_mean[season].sel(year=slice(year_min, year_max)) #.isel(year=slice(1,len(ann_seas_mean['DJF'].year)))
	return ann_seasonal_mean


def jas_seasonal_mean(ds):
    """Climatological JAS mean weighted by days_in_month. Assumes no leap years."""
    jas = ds.sel(time=ds['time.month'].isin([7, 8, 9]))
    weights = jas['time'].dt.days_in_month
    return (jas * weights).sum(dim='time') / weights.sum(dim='time')


def jas_yearly_mean(ds):
    """Per-year JAS mean weighted by days_in_month. Assumes no leap years."""
    jas = ds.sel(time=ds['time.month'].isin([7, 8, 9]))
    weights = jas['time'].dt.days_in_month
    return ((jas * weights).groupby('time.year').sum(dim='time') /
            weights.groupby('time.year').sum(dim='time'))


def match_lat_lon_names(ds):
    """Rename non-standard lat/lon coordinate names to 'lat' and 'lon'."""
    for lat_name in ['y', 'latitude', 'nav_lat']:
        if lat_name in ds.coords and 'lat' not in ds.coords:
            ds = ds.rename({lat_name: 'lat'})
    for lon_name in ['x', 'longitude', 'nav_lon']:
        if lon_name in ds.coords and 'lon' not in ds.coords:
            ds = ds.rename({lon_name: 'lon'})
    return ds


def windSpd(u, v):
    """Compute wind speed magnitude from u and v components."""
    return np.sqrt(u**2 + v**2)


def find_coastline_offset(field, windows, land_thresh=0, lon_name='lon', lat_name='lat',
                           resolution='50m', max_shift_deg=1.0):
    """
    Empirically measure the rigid lon/lat shift that best aligns a gridded field's
    implied land/ocean mask with the cartopy Natural Earth coastline.

    Diagnoses the grid-registration bug documented for the OIPC isoscape in
    nam-dD-lig/analysis-repo/DATA_MANIFEST.md ("OIPC grid registration"): some gridded
    products carry a coordinate array offset from the true lon/lat of each cell by a
    fixed number of grid cells, so the field disagrees with vector coastlines even
    though the data values themselves are fine. Sweeping a shift and re-scoring against
    Natural Earth is how that offset was found and corrected there; this generalizes
    the same check to any land/ocean-separable field (e.g. ETOPO05 topography).

    Parameters
    ----------
    field : xr.DataArray on a regular lat/lon grid, land > land_thresh, ocean/missing
            values <= land_thresh (e.g. raw elevation with negative or NaN ocean cells)
    windows : dict of {name: (lon_min, lon_max, lat_min, lat_max)}, degrees matching
              field's lon convention. Use several disjoint coastal sub-regions: a shift
              that agrees across all of them is a global registration issue; a shift
              that differs by window is local coastline-resolution disagreement, not a
              registration bug.
    land_thresh : threshold separating land from ocean/missing in `field`
    resolution : Natural Earth polygon resolution ('110m', '50m', '10m')
    max_shift_deg : search radius in degrees around zero shift

    Returns
    -------
    dict of {window_name: (best_agreement_fraction, dlon, dlat)}
    """
    import cartopy.io.shapereader as shpreader
    from shapely.ops import unary_union
    try:
        from shapely import contains_xy
    except ImportError:
        from shapely.vectorized import contains as _contains
        def contains_xy(geom, x, y):
            return _contains(geom, x, y)

    shp = shpreader.natural_earth(resolution=resolution, category='physical', name='land')
    land_union = unary_union(list(shpreader.Reader(shp).geometries()))

    lon = field[lon_name].values
    lat = field[lat_name].values
    dlon = float(lon[1] - lon[0])
    dlat = float(lat[1] - lat[0])
    lon2d, lat2d = np.meshgrid(lon, lat)
    lon2d_signed = ((lon2d + 180) % 360) - 180  # Natural Earth uses -180:180

    field_land = np.asarray(field.values > land_thresh)

    shifts_lon = np.arange(-max_shift_deg, max_shift_deg + dlon / 2, dlon)
    shifts_lat = np.arange(-max_shift_deg, max_shift_deg + dlat / 2, dlat)

    results = {}
    for name, (lon_min, lon_max, lat_min, lat_max) in windows.items():
        sub = ((lon2d >= lon_min) & (lon2d <= lon_max) &
               (lat2d >= lat_min) & (lat2d <= lat_max))
        sub_lon = lon2d_signed[sub]
        sub_lat = lat2d[sub]
        sub_field_land = field_land[sub]

        best = (-1.0, 0.0, 0.0)
        for slon in shifts_lon:
            for slat in shifts_lat:
                ref_land = contains_xy(land_union, sub_lon + slon, sub_lat + slat)
                agree = np.mean(ref_land == sub_field_land)
                if agree > best[0]:
                    best = (agree, round(float(slon), 6), round(float(slat), 6))
        results[name] = best

    return results
