# -*- coding: utf-8 -*-
"""
Created on Thu Nov 14 15:17:53 2019

@author: engelen
"""

import xarray as xr
from glob import glob
import os, sys
import cftime
import numpy as np
from copy import deepcopy

if len(sys.argv) > 1:
    res_fol  = sys.argv[1]
    test = False
else:
    #Local testing on my own windows laptop
#    res_fol = r"g:\synthdelta\results\test_output\synth_SD_i175_m24_7120041"
    res_fol = r"g:\synthdelta\results\test_output\synth_SD_i029_m24_7065039"
    test = True

res_paths = glob(os.path.join(res_fol, "results_[0-9][0-9][0-9].nc"))
res_paths.sort()

ds_ls = [xr.open_dataset(path) for path in res_paths]

times_peek = ds_ls[0].time.values
dt = times_peek[2].year - times_peek[1].year
assert(dt==1000.)

n_ts = [ds_ls[i].time.shape[0] for i in range(len(ds_ls))]

units = "days since 0000-01-01"
calendar = '365_day'

t_start = 0
times = []
for n in n_ts:
    times_sub = np.arange(1, n+1)*dt+t_start
    t_start = times_sub[-1]
    times.append(times_sub)

times = [cftime.num2date(
        np.array(time)*365, units=units, calendar=calendar
        ) for time in times]

#%%Close files linked to data
for ds in ds_ls:
    ds = ds.load()
    ds.close()

#%%Create backup for test before adapting the Dataset
if test:
    ds_ls_backup = deepcopy(ds_ls)
    
#%%Override timesteps
for i, time in enumerate(times):
    ds_ls[i]["time"] = time
    ds_ls[i]["time"].encoding["units"] = units
    ds_ls[i]["time"].encoding["calendar"] = calendar
    ds_ls[i]=ds_ls[i].transpose("time", "z", "y", "x")
    
#%%Test
if test:
    for i, ds in enumerate(ds_ls):
        assert(ds.time.equals(ds_ls_backup[i].time))

#%%Override timesteps
for i, ds in enumerate(ds_ls):
    ds.to_netcdf(res_paths[i])