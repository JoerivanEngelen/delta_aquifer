# -*- coding: utf-8 -*-
"""
Created on Tue May  7 11:24:25 2019

@author: engelen
"""

import imod
import xarray as xr
import sys, os
from collections import OrderedDict

modelfol = sys.argv[1]
mod_nr = int(sys.argv[2])
s_nr = int(sys.argv[3])
##For testing
#modelfol = r"g:\synthdelta\test_output\test_nspecies"
#mod_nr = 0
#s_nr = 12

#%%Path management
varnames = ["bdgfff", "bdgflf", "bdgfrf", "conc", "dcdt", "head"]
path_template = os.path.join(modelfol, "results" ,"{}", "*_p{:03d}.[iI][dD][fF]")

nc_path = os.path.join(modelfol, "results", "results_{:03d}_p{:03d}.nc".format(mod_nr, s_nr))

#%%Process
#import re
#pattern = r"(?P<name>[\w-]+?)_(?P<time>[\w-]+?)(_c(?P<species>[\w-]+))?_l(?P<layer>[\w-]+?)_p\d*"
#pattern = re.compile(pattern)

ds = OrderedDict()

for varname in varnames:
    path = path_template.format(varname, s_nr)
    
    if varname == "conc":
        pattern = r"{name}_{time}_c{species}_l{layer}_p\d*"
    else:
        pattern = r"{name}_{time}_l{layer}_p\d*"
    try:
        ds[varname] = imod.idf.open(path, use_cftime=True, pattern=pattern)
    except FileNotFoundError:
        print("Could not find files for: %s" % varname)

ds = xr.merge(ds.values())
ds.time.encoding["units"] = r"days since 1999-01-01"
ds.time.encoding["calendar"] = r"proleptic_gregorian"
ds.to_netcdf(nc_path, unlimited_dims=["time"])