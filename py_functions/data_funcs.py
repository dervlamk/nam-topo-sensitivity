import os
import sys
import xarray as xr
#import cf_xarray as cf
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
