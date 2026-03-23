import sys
sys.dont_write_bytecode = True

import xarray as xr
xr.set_options(keep_attrs=True)
import numpy as np
import metpy.calc as mp

import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from shapely.geometry.polygon import LinearRing

import matplotlib as mpl
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.ticker as mticker
from matplotlib import cm
import matplotlib.colors as mcolors
from matplotlib.colors import ListedColormap,LinearSegmentedColormap
import cmocean

from colorbar_funcs import *
from data_funcs import *

#############################

def quick_map(var, season='ann', field=None, diff=False, cbar=True, mask=None, mask_color='w', boundaries=None, label=None, ax=None, save=False, ofile=''):
    """
    Quickly plot up climate field data based preset values for various vars
    - only works for monthly climatology data where the time dimension is labeled 'month'
    """
    ### +++ GET VAR INFO +++ ###
    lons,lats = get_xy_coords(var) # get lat & lon coords without having to know coordinate names

    ## need to figure out how to write an if statement to only do this if there is a time/month variable
    if season!=None:
        mons = get_season(season=season) # index months to average over based on season var
        var_avg = var[mons].mean(dim='month', keep_attrs=True) # find seasonal or annual var mean
    else:
        var_avg=var

    ### +++ INITIALIZE FIGURE +++ ###
    trans = ccrs.PlateCarree()
    # create a new figure if no axes are designated
    if ax==None:
        fig = plt.figure(figsize=(9, 6))
        proj = ccrs.PlateCarree()
        ax = plt.subplot(111, projection=proj)

    ### +++ MAP & BOUNDARY INFO +++ ###
    # add optional masks
    if mask != None:
        if mask in ['land', 'Land']:
            ax.add_feature(cfeature.LAND, fc=mask_color, zorder=2) # masks continents
        if mask in ['ocean', 'Ocean']:
            ax.add_feature(cfeature.OCEAN, fc=mask_color, zorder=2) # masks oceans
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

    ### +++ COLORMAP INFO +++ ###
    cmap, vmin, vmax, cf = get_settings(field=field, diff=diff)

    ### +++ PLOT CONTOURS +++ ###
    ax.pcolormesh(lons, lats, var_avg, cmap=cmap, vmin=vmin, vmax=vmax, transform=trans)

    ### +++ ADD FIG LABEL +++ ###
    if label==None:
        label=var_avg.name
    else:
        label=label
    t=ax.annotate(label,
                  xy=(0.015, 0.025), xycoords='axes fraction',
                  annotation_clip=False,
                  fontsize=12, fontweight='bold', ha='left', va='bottom', zorder=100)
    t.set_bbox(dict(facecolor='white', alpha=1, edgecolor='k'))

    ### +++ COLORBAR INFO +++ ###
    if cbar==True:
        cbar = plt.colorbar(cf,
                            extend='both',
                            orientation='vertical',
                            shrink=0.5,
                            location='right',
                            pad=0.01,
                            ax=ax
                            )
        if 'units' in var.attrs:
            clb = '['+var.attrs['units']+']'
            cbar.set_label(clb, labelpad=15, rotation=270, fontweight='bold', ha='center')
        cbar.ax.tick_params(labelsize=10)
        for tick in cbar.ax.yaxis.get_major_ticks():
            tick.label2.set_fontweight('bold')
        cbar.ax.minorticks_off()
    else:
        pass

    if save is True:
            plt.savefig('{}.pdf'.format(ofile),
                        dpi=None, facecolor='w', edgecolor='w',
                        orientation='portrait', papertype=None,
                        format=None, transparent=False,
                        bbox_inches='tight', pad_inches=0.1)


###

def custom_map(var, season='ann', vmin=None, vmax=None, scale=None, cmap=None, cbar=True, mask=None, mask_color='w', boundaries=None, label=None, ax=None, save=False, ofile=''):
    
    ### +++ GET VAR INFO +++ ###
    lons,lats = get_xy_coords(var) # get lat & lon coords without having to know coordinate names
    if season!=None:
        mons = get_season(season=season) # index months to average over based on season var
        var_avg = var[mons].mean(dim='month', keep_attrs=True) # find seasonal or annual var mean
    else:
        var_avg=var

    ### +++ INITIALIZE FIGURE +++ ###
    trans = ccrs.PlateCarree()
    # create a new figure if no axes are designated
    if ax==None:
        fig = plt.figure(figsize=(9, 6))
        proj = ccrs.PlateCarree()
        ax = plt.subplot(111, projection=proj)

    ### +++ MAP & BOUNDARY INFO +++ ###
    # add optional masks
    if mask != None:
        if mask in ['land', 'Land']:
            ax.add_feature(cfeature.LAND, fc=mask_color, zorder=2) # masks continents
        if mask in ['ocean', 'Ocean']:
            ax.add_feature(cfeature.OCEAN, fc=mask_color, zorder=2) # masks oceans
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

    ### +++ COLORMAP INFO +++ ###
    if ((cmap==None) and (var_avg.min()<0)):
        cmap=plt.colormaps['RdBu']
    elif ((cmap==None) and (var_avg.min()>=0)):
        cmap=plt.colormaps['viridis_r']

    ### +++ PLOT CONTOURS +++ ###
    if scale==None:
        scale=0
    if ((vmin!=None) and (vmax!=None)):
        vn=vmin*(10**scale)
        vx=vmax*(10**scale)
        n_lvls = 2*(vx-vn)+1
        if n_lvls<=7:
            n_lvls=(n_lvls*2)+1
        levels = np.linspace(vmin, vmax, int(n_lvls))
        norm = mpl.colors.BoundaryNorm(levels, cmap.N)
        cf = ax.pcolormesh(lons, lats, var_avg, cmap=cmap, norm=norm, transform=trans)
    else:
        if vmin==None:
            vmin=var_avg.min()
        if vmax==None:
            vmax=var_avg.max()
        cf = ax.pcolormesh(lons, lats, var_avg, cmap=cmap, norm=norm, transform=trans)
    
    ### +++ ADD FIG LABEL +++ ###
    if label==None:
        label=var_avg.name
    else:
        label=label
    t=ax.annotate(label,
                  xy=(0.015, 0.025), xycoords='axes fraction',
                  annotation_clip=False,
                  fontsize=12, fontweight='bold', ha='left', va='bottom', zorder=100)
    t.set_bbox(dict(facecolor='white', alpha=1, edgecolor='k'))

    ## colorbar info
    if cbar==True:
        cbar = plt.colorbar(cf,
                            extend='both',
                            orientation='vertical',
                            shrink=0.5,
                            location='right',
                            pad=0.01,
                            ax=ax
                            )
        if 'units' in var.attrs:
            clb = '['+var.attrs['units']+']'
            cbar.set_label(clb, labelpad=15, rotation=270, fontweight='bold', ha='center')  
        cbar.ax.tick_params(labelsize=10)
        for tick in cbar.ax.yaxis.get_major_ticks():
            tick.label2.set_fontweight('bold')
        cbar.ax.minorticks_off()
    else:
        pass

    if save is True:
            plt.savefig('{}.pdf'.format(ofile),
                        dpi=None, facecolor='w', edgecolor='w',
                        orientation='portrait', papertype=None,
                        format=None, transparent=False,
                        bbox_inches='tight', pad_inches=0.1)


###

def wind_vectors(u, v, season='ann', color='k', scalef=15, w=0.005, skip_n=2, key_length=5,
                 mask=None, mask_color=None, boundaries=None, label=None, grid=False, ax=None):
    """
    Make a map of wind vectors from U and V components
    !!!! Cannot use this function with a projection that has a central_longitude=180 !!!!
    """
    ### +++ GET VAR INFO +++ ###
    lons,lats = get_xy_coords(u) # get lat & lon coords without having to know coordinate names
    #z_name = u.metpy.vertical.name # get vertical coordinate
    if season!=None:
        mons=get_season(season=season) # index months to average over based on season var
        u_avg=u[mons].mean(dim=['month'], keep_attrs=True)
        v_avg=v[mons].mean(dim=['month'], keep_attrs=True)
    else:
        u_avg=u
        v_avg=v

    ### +++ INITIALIZE FIGURE +++ ###
    trans = ccrs.PlateCarree()
    # create a new figure if no axes are designated
    if ax==None:
        fig = plt.figure(figsize=(9, 6))
        proj = ccrs.PlateCarree()
        ax = plt.subplot(111, projection=proj)

    ### +++ MAP & BOUNDARY INFO +++ ###
    # add optional masks
    if mask != None:
        if mask in ['land', 'Land']:
            ax.add_feature(cfeature.LAND, fc=mask_color, zorder=2) # masks continents
        if mask in ['ocean', 'Ocean']:
            ax.add_feature(cfeature.OCEAN, fc=mask_color, zorder=2) # masks oceans
    # map boundaries
    if boundaries != None:
        ax.set_extent(boundaries, crs=trans) # clips map extent according to designated boundaries
    else:
        ax.set_global() # make a global map
    # add coastlines
    ax.coastlines()
    if grid==True:
        # grid line specs
        gl = ax.gridlines(crs=trans, lw=.5, colors='black', alpha=1.0, linestyle='--', zorder=10, draw_labels=True)
        gl.top_labels=False
        gl.right_labels=False
        gl.xformatter = LONGITUDE_FORMATTER
        gl.yformatter = LATITUDE_FORMATTER
        gl.xlabel_style = {'color': 'black', 'weight': 'bold'}
        gl.ylabel_style = {'color': 'black', 'weight': 'bold'}
    else:
        ax.gridlines(crs=trans, alpha=0, colors=None, draw_labels=False)

    # Draw vectors
    q1=ax.quiver(lons[::skip_n], lats[::skip_n], u_avg[::skip_n,::skip_n], v_avg[::skip_n,::skip_n],
                 color=color, width=w, scale=scalef, scale_units='inches', units='height', transform=trans)
    ax.quiverkey(q1, .925, 1.02, key_length, rf'{key_length} m/s', labelcolor=color, labelpos='W')

###

def plot_field_diff(var1, var2, season=None,
                    vmin_fld=None, vmax_fld=None, cmap_fld=None,
                    vmin_diff=None, vmax_diff=None, cmap_diff=None,
                    mask=None, boundaries=None, cl=None,
                    var1_name='', var2_name='',
                    save=False, ofile=None):
    
    ## initialize figure
    trans = ccrs.PlateCarree()
    if cl != None:
        proj = ccrs.PlateCarree(central_longitude=cl)
    else:
        proj = ccrs.PlateCarree()
    fig, ax = plt.subplots(nrows=1, ncols=3, figsize=(20,4), subplot_kw={'projection': proj})
    fig.subplots_adjust(bottom=0.05, top=0.95)

    ## get var info
    lons, lats = get_xy_coords(var1) # get lat & lon coords without having to know coordinate names

    if season != None:
        mons = get_season(season=season) # index months to average over based on season var
        # find seasonal or annual var mean
        var_avg1 = var1[mons].mean(dim='month', keep_attrs=True) 
        var_avg2 = var2[mons].mean(dim='month', keep_attrs=True) 
        plt.suptitle((var1.attrs['long_name']+' ('+season+')'), fontsize=16, fontweight='bold', ha='center', va='center')
    else:
        var_avg1 = var1
        var_avg2 = var2
        plt.suptitle((var1.attrs['long_name']), fontsize=16, fontweight='bold', ha='center', va='center')

    ## colormamp info
    if ((cmap_fld==None) and (var_avg1.min()<0)):
        cmap_fld=plt.colormaps['RdBu']
    if ((cmap_fld==None) and (var_avg1.min()>=0)):
        cmap_fld=plt.colormaps['viridis_r']
    if cmap_diff==None:
        cmap_diff=plt.colormaps['RdBu']
    
    ## plot contours for fields
    if ((vmin_fld!=None) and (vmax_fld!=None)):
        vn_fld=vmin_fld
        vx_fld=vmax_fld
        n_lvls_fld = 2*(vx_fld-vn_fld)+1
        levels_fld = np.linspace(vn_fld, vx_fld, n_lvls_fld)
        norm_fld = mpl.colors.BoundaryNorm(levels_fld, cmap_fld.N)
        cf_fld = ax[0].pcolormesh(lons, lats, var_avg1, cmap=cmap_fld, norm=norm_fld, transform=trans)
        t=ax[0].annotate(var1_name, xy=(0.015, 0.025), xycoords='axes fraction', annotation_clip=False, fontsize=10, fontweight='bold', ha='left', va='bottom')
        t.set_bbox(dict(facecolor='white', alpha=1, edgecolor='k'))
        ax[1].pcolormesh(lons, lats, var_avg2, cmap=cmap_fld, norm=norm_fld, transform=trans)
        t=ax[1].annotate(var2_name, xy=(0.015, 0.025), xycoords='axes fraction', annotation_clip=False, fontsize=10, fontweight='bold', ha='left', va='bottom')
        t.set_bbox(dict(facecolor='white', alpha=1, edgecolor='k'))
    else:
        if vmin_fld==None:
            vn_fld=var_avg1.min()
        if vmax_fld==None:
            vx_fld=var_avg1.max()
        cf = ax[0].pcolormesh(lons, lats, var_avg1, cmap=cmap_fld, vmin=vn_fld, vmax=vx_fld, transform=trans)
        t=ax[0].annotate(var1_name, xy=(0.015, 0.025), xycoords='axes fraction', annotation_clip=False, fontsize=10, fontweight='bold', ha='left', va='bottom')
        t.set_bbox(dict(facecolor='white', alpha=1, edgecolor='k'))
        ax[1].pcolormesh(lons, lats, var_avg2, cmap=cmap_fld, vmin=vn_fld, vmax=vx_fld, transform=trans)
        t=ax[1].annotate(var2_name, xy=(0.015, 0.025), xycoords='axes fraction', annotation_clip=False, fontsize=10, fontweight='bold', ha='left', va='bottom')
        t.set_bbox(dict(facecolor='white', alpha=1, edgecolor='k'))
        
    ## plot contours for difference
    var_diff=var_avg2-var_avg1
    if ((vmin_diff!=None) and (vmax_diff!=None)):
        vn_diff=vmin_diff
        vx_diff=vmax_diff
        n_lvls_diff = 2*(vx_diff-vn_diff)+1
        levels_diff = np.linspace(vn_diff, vx_diff, n_lvls_diff)
        norm_diff = mpl.colors.BoundaryNorm(levels_diff, cmap_diff.N)
        cf_diff = ax[2].pcolormesh(lons, lats, var_diff, cmap=cmap_diff, norm=norm_diff, transform=trans)
        t=ax[2].annotate(f'{var2_name}$-${var1_name}', xy=(0.015, 0.025), xycoords='axes fraction', annotation_clip=False, fontsize=10, fontweight='bold', ha='left', va='bottom')
        t.set_bbox(dict(facecolor='white', alpha=1, edgecolor='k'))
    else:
        if vmin_diff==None:
            vn_diff=var_diff.min()
        if vmax_diff==None:
            vx_diff=var_diff.max()
        cf = ax[2].pcolormesh(lons, lats, var_diff, cmap=cmap_diff, vmin=vn_diff, vmax=vx_diff, transform=trans)
        t=ax[2].annotate(f'{var2_name}$-${var1_name}', xy=(0.015, 0.025), xycoords='axes fraction', annotation_clip=False, fontsize=10, fontweight='bold', ha='left', va='bottom')
        t.set_bbox(dict(facecolor='white', alpha=1, edgecolor='k'))
    
    for i in [0,1,2]:
        ## map & boundary info
        # add optional masks
        if mask != None:
            if mask in ['land', 'Land']:
                ax[i].add_feature(cfeature.LAND, fc=(1, 1, 1), zorder=2) # masks continents
            if mask in ['ocean', 'Ocean']:
                ax[i].add_feature(cfeature.OCEAN, fc=(1, 1, 1), zorder=2) # masks oceans
        # map boundaries
        if boundaries != None:
            ax[i].set_extent(boundaries, crs=trans) # clips map extent according to designated boundaries
        else:
            ax[i].set_global() # make a global map
        # add coastlines
        ax[i].coastlines()
        # grid line specs
        gl = ax[i].gridlines(crs=trans,
                             lw=.5,
                             colors='black',
                             alpha=1.0,
                             linestyle='--',
                             zorder=10,
                             draw_labels=True)
        gl.top_labels=False
        gl.right_labels=False
        gl.xformatter = LONGITUDE_FORMATTER
        gl.yformatter = LATITUDE_FORMATTER
        gl.xlabel_style = {'color': 'black', 'weight': 'bold'}
        gl.ylabel_style = {'color': 'black', 'weight': 'bold'}
    
        ## colorbar info
        if 'units' in var_avg1.attrs:
            clb = var_avg1.attrs['units']
        if i in [0,1]:
            cbar = plt.colorbar(cf_fld,
                                extend='both',
                                orientation='horizontal',
                                shrink=0.9,
                                location='bottom',
                                pad=0.075,
                                ax=ax[i])
            cbar.set_label(clb, labelpad=5, y=1, rotation=0, fontweight='bold')
            cbar.ax.tick_params(labelsize=10)
            for tick in cbar.ax.xaxis.get_major_ticks():
                tick.label1.set_fontweight('bold')
        else:
            cbar = plt.colorbar(cf_diff,
                                extend='both',
                                orientation='horizontal',
                                shrink=0.9,
                                location='bottom',
                                pad=0.075,
                                ax=ax[i])
            cbar.set_label(clb, labelpad=5, y=1, rotation=0, fontweight='bold')
            cbar.ax.tick_params(labelsize=10)
            for tick in cbar.ax.xaxis.get_major_ticks():
                tick.label1.set_fontweight('bold')

    if save is True:
            plt.savefig('{}.pdf'.format(ofile),
                        dpi=None, facecolor='w', edgecolor='w',
                        orientation='portrait', papertype=None,
                        format=None, transparent=False,
                        bbox_inches='tight', pad_inches=0.1)
