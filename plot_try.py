# -*- coding: utf-8 -*-
"""
Spyder Editor

This is the script file to check difference.
"""

import numpy as np
import xarray as xr
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
from cartopy.util import add_cyclic_point
import matplotlib.gridspec as gridspec
import scipy
import os
figurePlot=False
# Data process
# original process
# ttttttt
# Function definition
# To get phase-averaged values
def read_merge_nc(filedir):
    filelist      = os.listdir(filedir)
    filenum       = True
    print('These files are read and merged:\n')
    for file_ in filelist:
        if file_.endswith('.nc') and int(file_[-8:-6])>30:
            print(file_)
            dstmp = xr.open_dataset(filedir+file_)
            if filenum:
                dsmge = dstmp
                filenum = False
            dsmge = xr.merge([dsmge, dstmp])
    return dsmge
def get_field(da, field_script, average_time):
    if len(average_time) == 0:
        da_processed = da[field_script].mean(dim=['time'])
    else:
        da_season = da[field_script].groupby('time.season').mean(dim='time')
        if average_time == 'MAM':
            da_processed = da_season.loc[da_season['season']=='MAM']
        elif average_time == 'JJA':
            da_processed = da_season.loc[da_season['season']=='JJA']
        elif average_time == 'SOM':
            da_processed = da_season.loc[da_season['season']=='SOM']
        elif average_time == 'DJF':
            da_processed = da_season.loc[da_season['season']=='DJF']
    if 'PREC' in field_script:
        return 24*60*60*da_processed
    else:
        return da_processed
# To get full lons without blanks
def da_full_lons(da):
    lons = da.lon
    lats = da.lat
    da_cycle, lon_cycle = add_cyclic_point(da,coord=lons)
    lons_cycle, lats_cycle   = np.meshgrid(lon_cycle, lats)
    return da_cycle, lons_cycle, lats_cycle
# To prepare uniform plot process
def da_plot(ax, projection, lons, lats, da, title_script, cmap_script, cbar_script):
    ax.set_extent([-180, 180, -90, 90], crs=projection)
    ax.coastlines()
    ax.set_global()
    ax.set_xticks(np.linspace(-180, 180, 5), crs=projection)
    ax.set_yticks(np.linspace(-90, 90, 5), crs=projection)
    lon_formatter = LongitudeFormatter(zero_direction_label=True)
    lat_formatter = LatitudeFormatter()
    ax.xaxis.set_major_formatter(lon_formatter)
    ax.yaxis.set_major_formatter(lat_formatter)
    pcm=ax.contourf(lons, lats, da,
                    transform=projection,
                    cmap=cmap_script, extend='both')
    ax.set_title(title_script)
    #pos = ax.get_position()
    #cax = fig.add_axes([pos.x0, 0.06, pos.x1-pos.x0, 0.02])
    plt.colorbar(pcm, orientation='vertical', extend='both')
    pcm.colorbar.set_label(cbar_script)
# Data pre-process
# Read in data
print('NC file processing...')
filebase                = '/Users/zhehao/Documents/CLM/f09_g16-veg_test/nc_output/f09_g16/'
f09_g16                 = read_merge_nc(filebase)
fileveg                 = '/Users/zhehao/Documents/CLM/f09_g16-veg_test/nc_output/f09_g16_veg/'
f09_g16_veg             = read_merge_nc(fileveg)
# T and Precc
time_required           = 'JJA'
print('Required field slicing (in %s)...'%(time_required))
f09_g16_t, lons, lats   = da_full_lons(get_field(f09_g16,'TS',time_required))
f09_g16_t_v, lons, lats = da_full_lons(get_field(f09_g16_veg,'TS',time_required))
f09_g16_p, lons, lats   = da_full_lons(get_field(f09_g16,'PRECC',time_required))
f09_g16_p_v, lons, lats = da_full_lons(get_field(f09_g16_veg,'PRECC',time_required))
# ttest
print('T-test validating...')
f09_g16_ts_da, lons, lats = da_full_lons(f09_g16['TS']);
f09_g16_ts_veg_da, lons, lats = da_full_lons(f09_g16_veg['TS']);
t_ts, p_ts = scipy.stats.ttest_ind(f09_g16_ts_da, f09_g16_ts_veg_da, equal_var=True, nan_policy='propagate')
f09_g16_pr_da, lons, lats = da_full_lons(f09_g16['PRECC']);
f09_g16_pr_veg_da, lons, lats = da_full_lons(f09_g16_veg['PRECC']);
t_pr, p_pr = scipy.stats.ttest_ind(f09_g16_pr_da, f09_g16_pr_veg_da, equal_var=True, nan_policy='propagate')
# Figure plot
# Plot preparation
print('Figure plotting...')
projection = ccrs.PlateCarree(central_longitude=0)
plt.close('all')
# Temperature plot: delta, original, veg_changed
fig = plt.figure(figsize=(10,12))
gs1 = gridspec.GridSpec(3, 1)
ax1 = fig.add_subplot(gs1[0:2], projection=projection)
da_plot(ax1, projection, lons, lats, f09_g16_t_v[0]-f09_g16_t[0], 'dTS (%s)'%(time_required), 'RdBu_r', 'Delta Temperature (K)')
dot = ax1.contourf(lons,lats, p_ts,[np.min(p_ts),0.05,np.max(p_ts)], zorder=1,hatches=['.', None],colors="none", transform=ccrs.PlateCarree())
gs2 = gridspec.GridSpec(3, 2)
ax2 = fig.add_subplot(gs2[4], projection=projection)
da_plot(ax2, projection, lons, lats, f09_g16_t[0], 'TS (%s)'%(time_required), 'RdBu_r', 'Temperature (K)')
ax3 = fig.add_subplot(gs2[5], projection=projection)
da_plot(ax3, projection, lons, lats, f09_g16_t_v[0], 'TS_veg (%s)'%(time_required), 'RdBu_r', 'Temperature (K)')
plt.tight_layout()
fig.savefig('./f09_g16_T_%s.jpeg'%(time_required), dpi=600, bbox_inches='tight')
# Precipitation plot: delta, original, veg_changed
fig = plt.figure(figsize=(10,12))
gs1 = gridspec.GridSpec(3, 1)
ax1 = fig.add_subplot(gs1[0:2], projection=projection)
da_plot(ax1, projection, lons, lats, f09_g16_p_v[0]-f09_g16_p[0], 'dPRECC (%s)'%(time_required), 'BrBG', 'Delta Precipitation (mm/d)')
dot = ax1.contourf(lons,lats, p_pr,[np.min(p_pr),0.05,np.max(p_pr)], zorder=1,hatches=['.', None],colors="none", transform=ccrs.PlateCarree())
gs2 = gridspec.GridSpec(3, 2)
ax2 = fig.add_subplot(gs2[4], projection=projection)
da_plot(ax2, projection, lons, lats, f09_g16_p[0], 'PRECC (%s)'%(time_required), 'BrBG', 'Precipitation (mm/d)')
ax3 = fig.add_subplot(gs2[5], projection=projection)
da_plot(ax3, projection, lons, lats, f09_g16_p_v[0], 'PRECC_veg (%s)'%(time_required), 'BrBG', 'Precipitation (mm/d)')
plt.tight_layout()
plt.show()
if figurePlot:
    fig.savefig('./f09_g16_PRECC_%s.jpeg'%(time_required), dpi=600, bbox_inches='tight')
#


