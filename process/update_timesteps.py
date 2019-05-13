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
#res_folder = sys.argv[1]
res_folder = r"g:\synthdelta\test_cftime_update"

nc_paths = natural_sort(glob(os.path.join(res_folder, "*.nc")))

#%%Processing
ds_ls = [xr.open_dataset(path, use_cftime=True) for path in nc_paths]

start_year = cftime.num2date(0., ds_ls[0].time.encoding["units"], 
                                 ds_ls[0].time.encoding["calendar"])

timedeltas = [(ds.time.values - start_year) for ds in ds_ls]

#Update timedeltas with ends
ends = np.cumsum([timedeltas[i][-1] for i in range(len(timedeltas)-1)])
timedeltas[1:] = [(timedeltas[i+1] + ends[i]) for i in range(len(timedeltas)-1)]

years = [[np.round(timedelta.days/365.25, 0) for timedelta in timedelta_arr] for timedelta_arr in timedeltas]