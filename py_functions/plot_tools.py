import sys
sys.dont_write_bytecode = True

import xarray as xr
import numpy as np
import metpy.calc as mp
from misc_functions import *

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

def plot_field(var, cmap=None, vmin=None, vmax=None, season='ann', mask=None, boundaries=None):    
    ## get var info
    get_xy_coords(var) # get lat & lon coords without having to know coordinate names
    lons = x
    lats = y
    mons = get_season(season=season) # index months to average over based on season var
    var_avg = var[mons].mean(dim='month', keep_attrs=True) # find seasonal or annual var mean
    
    ## initialize figure
    trans = ccrs.PlateCarree()
    proj = ccrs.PlateCarree()
    fig = plt.figure(figsize=(9, 6))
    ax = plt.subplot(111, projection=proj)
    
    ## map & boundary info
    # add optional masks
    if mask != None:
        if mask in ['land', 'Land']:
            ax.add_feature(cfeature.LAND, fc=(1, 1, 1), zorder=2) # masks continents
        if mask in ['ocean', 'Ocean']:
            ax.add_feature(cfeature.OCEAN, fc=(1, 1, 1), zorder=2) # masks oceans
    # map boundaries
    if boundaries != None:
        ax.set_extent(boundaries, crs=trans) # clips map extent according to designated boundaries
    else:
        ax.set_global() # make a global map
    # add coastlines
    ax.coastlines()
    # grid line specs
    gl = ax.gridlines(crs=trans, lw=.5, colors='black', alpha=1.0, linestyle='--', zorder=10, draw_labels=True)
    gl.top_labels=False
    gl.right_labels=False
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    gl.xlabel_style = {'color': 'black', 'weight': 'bold'}
    gl.ylabel_style = {'color': 'black', 'weight': 'bold'}
        
    ## colormamp info
    if ((cmap==None) and (var_avg.min()<0)):
        cmap=plt.colormaps['RdBu']
    elif ((cmap==None) and (var_avg.min()>=0)):
        cmap=plt.colormaps['viridis_r']
        
    ## plot contours
    if ((vmin!=None) and (vmax!=None)):
        vn=vmin
        vx=vmax
        n_lvls = 2*(vx-vn)+1
        levels = np.linspace(vn, vx, n_lvls)
        norm = mpl.colors.BoundaryNorm(levels, cmap.N)
        cf = mpl.cm.ScalarMappable(norm=norm, cmap=cmap)
        ax.pcolormesh(lons, lats, var_avg, cmap=cmap, vmin=vn, vmax=vx, transform=trans)
    else:
        if vmin==None:
            vn=var_avg.min()
        if vmax==None:
            vx=var_avg.max()
        cf = ax.pcolormesh(lons, lats, var_avg, cmap=cmap, vmin=vn, vmax=vx, transform=trans)
        
    ## colorbar info
    cbar = fig.colorbar(cf,
                        extend='both',
                        orientation='vertical',
                        shrink=0.525,
                        location='right',
                        pad=0.025,
                        ax=ax
                        )
    if 'units' in var.attrs:
        clb = var.attrs['units']
        cbar.set_label(clb, labelpad=-15, y=1.125, rotation=0)
    cbar.ax.tick_params(labelsize=10)
    for tick in cbar.ax.yaxis.get_major_ticks():
        tick.label2.set_fontweight('bold')
    
    ###
    plt.show()
