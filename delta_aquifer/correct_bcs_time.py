# -*- coding: utf-8 -*-
"""
Created on Fri Nov 22 12:53:45 2019

@author: engelen
"""

import sys, os
import xarray as xr
import numpy as np
import cftime

if len(sys.argv) > 1:
    model_fol  = sys.argv[1]
else:
    #Local testing on my own windows laptop
    model_fol = r"g:\synthdelta\results\test_output\synth_SD_i200_m24_7248732"
    
#%%Path management
path = os.path.join(model_fol, "input", "data", "bcs.nc")

#%%Load
ds = xr.open_dataset(path)
ds.load()
ds.close()

units = "days since 0000-01-01"
calendar = '365_day'

#%%Correct timesteps
timedeltas = ds.time.values - ds.time[0].values
years = np.round([timedelta.days/365 for timedelta in timedeltas], -1)
years = cftime.num2date(np.array(years)*365, units=units, calendar=calendar)

ds = ds.assign_coords(time=years)
ds["time"].encoding["units"] = units
ds["time"].encoding["calendar"] = calendar

#%%Save
ds.to_netcdf(path)