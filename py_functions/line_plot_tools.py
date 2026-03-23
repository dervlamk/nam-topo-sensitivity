import xarray as xr
import numpy as np
import metpy.calc as mp

import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from shapely.geometry.polygon import LinearRing

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.ticker as mticker
from matplotlib import cm
from matplotlib.colors import ListedColormap,LinearSegmentedColormap
import cmocean
import cmocean.cm as cmo
import colorcet as cc


#############################


def get_season(season='ann'):
    """
    This indexes which months to average over to derive an annual or seasonal mean
        - can only be applied to monthly climatologies
    """
    mons = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]
    if season in ['ANNUAL', 'ANN', 'ann']:
        mons = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]
    if season in ['DJF', 'djf']:
        mons = [0, 1, 11]
    if season in ['JJA', 'jja']:
        mons = [5, 6, 7]
    if season in ['JJAS', 'jjas']:
        mons = [5, 6, 7, 8]
    if season in ['JFM', 'jfm']:
        mons = [0, 1, 2]
    if season in ['JAS', 'jas']:
        mons = [6, 7, 8]
    if season==None:
        mons = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]
    return mons

####

def get_xy_coords(var):
    """
    This extracts the lat and lon coordinates without having to know the specific coordinate names
    """
    if isinstance(var, xr.DataArray):
        x = var.metpy.x
        y = var.metpy.y
        return(x,y)
    if isinstance(var, xr.Dataset):
        print('This is a dataset. Please use an xarray DataArray')

###

def weighted_global_mean_1d(var):
    """
    This weights var data to account for unequal grid cell areas that are a function of latitude,
    then calculates a global mean timeseries of zonal mean [time x lat] data
    """
    # find latitude variable
    lats=var.metpy.y
    # determine weight based on latitude value
    weights = np.cos(np.deg2rad(lats))
    weights.name = 'weights'
    # calculate area-weighted values
    weighted_var = var.weighted(weights)
    # calculate global mean of weighted data
    weighted_global_mean = weighted_var.mean(dim=[lats.name], keep_attrs=True)
    return(weighted_global_mean)

###

def weighted_global_mean_md(var):
    """
    This weights var data to account for unequal grid cell areas that are a function of latitude,
    then calculates a global mean timeseries for [time x lat x lon x n] data
    """
    # find latitude & longitude variable
    lons,lats=get_xy_coords(var)
    # determine weight based on latitude value
    weights = np.cos(np.deg2rad(lats))
    weights.name = 'weights'
    # calculate area-weighted values
    weighted_var = var.weighted(weights)
    # calculate global mean of weighted data
    weighted_global_mean = weighted_var.mean(dim=[lats.name,lons.name], keep_attrs=True)
    return(weighted_global_mean)

###

def plot_tseries(var):
    """
    Plot a time series
    """
    # get time variable (only works if time is first dimension)
    t_name = var.dims[0]
    time = var[t_name]
    
    # make plot
    fig = plt.figure(figsize=(9, 6))
    ax = plt.subplot(111)

    ax.plot(time, var, c='r', lw=2, linestyle='-', label=var.name) 
    
    # set axes information
    xlb = t_name
    if 'units' in var.attrs:
        ylb = var.attrs['units']
    ax.set(xlabel=xlb, ylabel=ylb)
    
    # legend
    ax.legend(loc=0, frameon=False, fontsize=10)

    ###
    plt.show()

###

def plot_zonal_mean(var, season='ann', col='k', xlim=None, ylim=None, mask=None, bounds=None, ax=None):
    """
    Plot zonal mean data
    """
    # use current axes if ax is not specified
    if ax is None:
        ax = plt.gca()

    ## get var info
    _, lats = get_xy_coords(var) # get lat coord without having to know the specific coordinate name
    mons = get_season(season=season) # index months to average over based on season var
    var_avg = var[mons].mean(dim='month', keep_attrs=True) # calculate seasonal or annual average
    """
    if mask!=None:
        # if a mask is provided, mask the data
        var_masked = var_avg.where(mask==1,np.nan) 
        # then calculate zonal mean
        var_zonal_mean = var_masked.mean(dim='lon', keep_attrs=True)
    else:
        # if not, just calculate zonal mean
    """
    var_zonal_mean = var_avg.mean(dim='lon', keep_attrs=True)
    
    ## draw plot
    ax.plot(lats, var_zonal_mean, c=col, lw=2, linestyle='-', label=var.source) 
    
    ## axes labels
    # if the var has a units attributes, use that for the y-axis label
    if 'units' in var.attrs:
        ylb = var.attrs['units']
    # otherwise, leave the y-label axis blank
    else:
        ylb = ' '
    xlb = 'LATITUDE (deg N)'
    ## set axes information
    ax.set(xlabel=xlb, ylabel=ylb)
    
    # legend
    legend_properties = {'weight':'bold',
                         'size':10}
    ax.legend(loc=0, frameon=False, labelcolor='linecolor', prop=legend_properties)



"""
    ## axes limits
    # if xlims are provided, use those
    if xlim!=None:
        xlim=xlim 
    # otherwise, use the min and max of the data's latitude coordinate
    else:
        xlim=[lats.min(),lats.max()] 
    # if ylims are provided, use those
    if ylim!=None:
        ylim=ylim 
    # otherwise, use the min and max of the dataset
    else:
        ylim=[var_zonal_mean.min(),var_zonal_mean.max()]
"""
