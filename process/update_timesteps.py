# -*- coding: utf-8 -*-
"""
Created on Mon May 13 16:27:48 2019

@author: engelen
"""

import xarray as xr
import numpy as np
import cftime
import sys, re, os
from glob import glob

def natural_sort(l): 
    convert = lambda text: int(text) if text.isdigit() else text.lower() 
    alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ] 
    return sorted(l, key = alphanum_key)

#%%Path management
res_folder = sys.argv[1]
#res_folder = r"g:\synthdelta\test_cftime_update"

nc_paths = natural_sort(glob(os.path.join(res_folder, "results_*.nc")))

#%%Processing

##We do not use xr.open_mfdataset, as as of now encoding attributes are lost when using open_mfdataset
ds_ls = [xr.open_dataset(path, use_cftime=True, 
                         chunks=dict(x = -1, y=-1, layer=10, time=-1)
                         ) for path in nc_paths]

start_year = cftime.num2date(0., ds_ls[0].time.encoding["units"], 
                                 ds_ls[0].time.encoding["calendar"])

timedeltas = [(ds.time.values - start_year) for ds in ds_ls]

#Update timedeltas with ends
ends = np.cumsum([timedeltas[i][-1] for i in range(len(timedeltas)-1)])
timedeltas[1:] = [(timedeltas[i+1] + ends[i]) for i in range(len(timedeltas)-1)]

units = "days since 0000-01-01"
calendar = '365_day'

#After 26000 years, 1 year off. For the sake of presentation, we round to decades.
years = [np.round([timedelta.days/365.25 for timedelta in timedelta_arr], -1) for timedelta_arr in timedeltas]
years = [cftime.num2date(np.array(years_arr)*365, units=units, calendar=calendar) for years_arr in years]

for i, year_arr in enumerate(years):
    ds_ls[i]["time"] = year_arr
    ds_ls[i]["time"].encoding["units"] = units
    ds_ls[i]["time"].encoding["calendar"] = calendar
    ds_ls[i].transpose("time", "layer", "y", "x")

for i, nc_path in enumerate(nc_paths):
#    nu_nc_path = nc_path.replace("results", "nu_results")
#    ds_ls[i].to_netcdf(nu_nc_path)
    #Can't save time as unlimited times. Luckily does not matter for Paraview.
    ds_ls[i].to_netcdf(nc_path)